"""G11: batched section verification is semantics-preserving.

- verify_section_edges_batched() must return records identical to the
  unbatched scalar reference path on a battery of section
  configurations (same verdicts, same errors, same trim decisions).
- Batching must actually batch: one batched call covering N section
  edges performs one classifier setup per face, not per sample.
- _point_segment_distance() must match a hand-computed float reference.

This test FAILS before the G11 implementation exists (ImportError on
the batched entry point) and PASSES after.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import brepkernel.intersection as ix_mod
from brepkernel.intersection import (
    _adaptive_edge_samples,
    _exact_curve_on_surface_distance,
    _p2,
    _repair_same_parameter,
    _run_section_engine,
    _section_tolerance_ceiling,
    _surface_d1,
    verify_section_edges_batched,  # noqa: F401  (fails before G11)
)
from brepkernel.step_ingest import index_shape

from OCP.BRep import BRep_Builder, BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCP import BRepClass as bcm
import OCP.BRepClass as _bcreal  # real module; bcm is OCP's lazy proxy
from OCP.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder,
                             BRepPrimAPI_MakeSphere)
from OCP.Geom import Geom_BSplineSurface, Geom_Plane
from OCP.TopAbs import TopAbs_FACE, TopAbs_IN, TopAbs_ON
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS, TopoDS_Shell
from OCP.gp import gp_Ax2, gp_Dir, gp_Pln, gp_Pnt, gp_Pnt2d

try:
    from OCP.TColgp import TColgp_Array2OfPnt
    from OCP.TColStd import TColStd_Array1OfReal, TColStd_Array1OfInteger
except ImportError:
    from OCP.collections import (
        Array2_gp_Pnt as TColgp_Array2OfPnt,
        Array1_double as TColStd_Array1OfReal,
        Array1_int as TColStd_Array1OfInteger,
    )


PASS = True


def check(name, cond, detail=""):
    global PASS
    print("[%s] %s %s" % ("PASS" if cond else "FAIL", name, detail))
    PASS = PASS and bool(cond)


# ---------------------------------------------------------------------------
# Battery construction
# ---------------------------------------------------------------------------

def _saddle_face():
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Shell
    n, p = 4, 3
    P = TColgp_Array2OfPnt(1, n, 1, n)
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            x = (i - 1) / (n - 1) * 2.0 - 1.0
            y = (j - 1) / (n - 1) * 2.0 - 1.0
            z = 0.3 * (x * x - y * y)
            P.SetValue(i, j, gp_Pnt(float(x), float(y), float(z)))
    K = TColStd_Array1OfReal(1, 2)
    M = TColStd_Array1OfInteger(1, 2)
    K.SetValue(1, 0.0)
    K.SetValue(2, 1.0)
    M.SetValue(1, 4)
    M.SetValue(2, 4)
    surf = Geom_BSplineSurface(P, K, K, M, M, p, p, False, False)
    face = BRepBuilderAPI_MakeFace(surf, 1e-7).Face()
    return face


def _tilted_plane_face():
    plane = Geom_Plane(gp_Pln(gp_Pnt(0, 0, 0.15), gp_Dir(0.2, 0, 1)))
    return BRepBuilderAPI_MakeFace(plane, -3, 3, -3, 3, 1e-7).Face()


def _shell(face):
    s = TopoDS_Shell()
    b = BRep_Builder()
    b.MakeShell(s)
    b.Add(s, face)
    return s


def _face_pairs_with_edges():
    """Yield (name, fa, fb, edges) with at least one section edge each."""
    configs = []

    # NURBS saddle vs tilted plane (analytic-ish freeform).
    ma = index_shape(_shell(_saddle_face()))
    mb = index_shape(_shell(_tilted_plane_face()))
    configs.append(("saddle_vs_plane", ma.faces[0], mb.faces[0]))

    # Crossing cylinders (analytic curved / curved).
    cyl_a = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(0, 0, -2), gp_Dir(0, 0, 1)), 1.0, 4.0).Shape()
    cyl_b = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(-2, 0, 0), gp_Dir(1, 0, 0)), 0.6, 4.0).Shape()
    ia, ib = index_shape(cyl_a), index_shape(cyl_b)
    configs.append(("crossing_cyls", ia.faces[0], ib.faces[0]))

    # Box vs box, overlapping (planar / planar with trim).
    box_a = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), gp_Pnt(2, 2, 2)).Shape()
    box_b = BRepPrimAPI_MakeBox(gp_Pnt(1, 1, 1), gp_Pnt(3, 3, 3)).Shape()
    ja, jb = index_shape(box_a), index_shape(box_b)
    configs.append(("boxes_overlap", ja.faces[0], jb.faces[2]))

    # Sphere vs box (curved / planar).
    sph = BRepPrimAPI_MakeSphere(gp_Pnt(1, 1, 1), 1.2).Shape()
    ka, kb = index_shape(sph), index_shape(box_a)
    configs.append(("sphere_vs_box", ka.faces[0], kb.faces[0]))

    out = []
    for name, fa, fb, in configs:
        edges, _verts = _run_section_engine(
            fa, fb, approximation=True, fuzzy=0.0,
            parallel=False, use_obb=True)
        if edges:
            out.append((name, fa, fb, edges))
            print("battery: %s -> %d section edge(s)" % (name, len(edges)))
    return out


# ---------------------------------------------------------------------------
# Scalar reference: the obviously-correct per-sample loop (mirrors the
# pre-G11 _verify_section_edge sample loop, including fresh classifier
# construction per sample).
# ---------------------------------------------------------------------------

def _scalar_reference(edge, fa, fb, edge_index, *, base_tol,
                      chord_tol, tangent_sin_tol, max_section_tol):
    et = float(BRep_Tool.Tolerance_s(edge))
    ft = max(float(BRep_Tool.Tolerance_s(fa.face)),
             float(BRep_Tool.Tolerance_s(fb.face)))
    tol_limit = _section_tolerance_ceiling(
        fa, fb, base_tol=base_tol, max_section_tol=max_section_tol)
    if et > tol_limit or ft > tol_limit:
        raise RuntimeError("battery edge exceeds tolerance ceiling")
    verify_tol = max(float(base_tol), 2.0 * et, 2.0 * ft)
    edge, _repaired = _repair_same_parameter(edge, verify_tol)
    et = float(BRep_Tool.Tolerance_s(edge))
    verify_tol = max(float(base_tol), 2.0 * et, 2.0 * ft)
    pc_a = BRep_Tool.CurveOnSurface_s(edge, fa.face, 0.0, 0.0)
    pc_b = BRep_Tool.CurveOnSurface_s(edge, fb.face, 0.0, 0.0)
    assert pc_a is not None and pc_b is not None
    c3 = BRepAdaptor_Curve(edge)
    local_chord = (max(2.0 * verify_tol, 1e-9)
                   if chord_tol is None else float(chord_tol))
    ts, xyz = _adaptive_edge_samples(edge, local_chord)
    sa = BRepAdaptor_Surface(fa.face)
    sb = BRepAdaptor_Surface(fb.face)

    uva, uvb, err_a, err_b, err_cross, trans = [], [], [], [], [], []
    trim_ok = True
    for t, p in zip(ts, xyz):
        qa = pc_a.Value(float(t))
        qb = pc_b.Value(float(t))
        ua, ub = _p2(qa), _p2(qb)
        pa, sua, sva = _surface_d1(sa, *ua)
        pb, sub, svb = _surface_d1(sb, *ub)
        err_a.append(float(np.linalg.norm(p - pa)))
        err_b.append(float(np.linalg.norm(p - pb)))
        err_cross.append(float(np.linalg.norm(pa - pb)))
        na = np.cross(sua, sva)
        nb = np.cross(sub, svb)
        den = float(np.linalg.norm(na) * np.linalg.norm(nb))
        sinang = (0.0 if den <= 1e-300 else
                  float(np.linalg.norm(np.cross(na, nb)) / den))
        trans.append(sinang)
        ca = _bcreal.BRepClass_FaceClassifier(
            fa.face, gp_Pnt2d(float(ua[0]), float(ua[1])), verify_tol, True)
        cb = _bcreal.BRepClass_FaceClassifier(
            fb.face, gp_Pnt2d(float(ub[0]), float(ub[1])), verify_tol, True)
        trim_ok = trim_ok and (
            ca.State() in (TopAbs_IN, TopAbs_ON)
            and cb.State() in (TopAbs_IN, TopAbs_ON))
        uva.append(ua)
        uvb.append(ub)
    return {
        "max_a": max(err_a, default=0.0),
        "max_b": max(err_b, default=0.0),
        "max_cross": max(err_cross, default=0.0),
        "trim_ok": trim_ok,
        "min_tr": min(trans, default=1.0),
        "max_tr": max(trans, default=1.0),
        "uv_a": np.vstack(uva) if uva else np.zeros((0, 2)),
        "uv_b": np.vstack(uvb) if uvb else np.zeros((0, 2)),
        "n": len(ts),
    }


def _close(a, b):
    # Batched NumPy reductions may differ from the scalar loop by 1 ULP;
    # that cannot change any geometric decision (thresholds are >= 1e-9).
    return a == b or math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-18)


def _records_equal(name, rec, ref):
    check(name + " max_surface_error_a",
          _close(rec.max_surface_error_a, ref["max_a"]),
          "%r vs %r" % (rec.max_surface_error_a, ref["max_a"]))
    check(name + " max_surface_error_b",
          _close(rec.max_surface_error_b, ref["max_b"]))
    check(name + " max_cross_surface_error",
          _close(rec.max_cross_surface_error, ref["max_cross"]))
    check(name + " trim_ok", rec.trim_ok == ref["trim_ok"])
    check(name + " min_transversality",
          _close(rec.min_transversality, ref["min_tr"]))
    check(name + " max_transversality",
          _close(rec.max_transversality, ref["max_tr"]))
    check(name + " uv_a",
          np.allclose(rec.uv_a, ref["uv_a"], rtol=0, atol=0))
    check(name + " uv_b",
          np.allclose(rec.uv_b, ref["uv_b"], rtol=0, atol=0))
    check(name + " sample count",
          rec.parameters.shape[0] == ref["n"])


class _CountingClassifier(bcm.BRepClass_FaceClassifier):
    constructions = 0

    def __init__(self, *args, **kwargs):
        _CountingClassifier.constructions += 1
        super().__init__(*args, **kwargs)


def main():
    kw = dict(base_tol=1e-7, chord_tol=None, tangent_sin_tol=1e-4,
              max_section_tol=None)
    battery = _face_pairs_with_edges()
    check("battery non-empty", len(battery) >= 3,
          "%d configs" % len(battery))
    if not battery:
        print("checks=1 failures=1")
        return 1

    # 1. Batched records identical to the scalar reference.
    for name, fa, fb, edges in battery:
        items = [(e, i) for i, e in enumerate(edges)]
        batched = verify_section_edges_batched(items, fa, fb, **kw)
        check(name + " record count", len(batched) == len(edges))
        for (e, i), rec in zip(items, batched):
            ref = _scalar_reference(e, fa, fb, i, **kw)
            _records_equal("%s edge %d" % (name, i), rec, ref)

    # 2. Batching actually batches: one call over N edges builds one
    # classifier per face, not per sample.
    name, fa, fb, edges = battery[0]
    items = [(e, i) for i, e in enumerate(edges)]
    _CountingClassifier.constructions = 0
    _real_classifier = _bcreal.BRepClass_FaceClassifier
    _bcreal.BRepClass_FaceClassifier = _CountingClassifier
    try:
        verify_section_edges_batched(items, fa, fb, **kw)
        n_batched = _CountingClassifier.constructions
    finally:
        _bcreal.BRepClass_FaceClassifier = _real_classifier
    check("batched: one classifier per face",
          n_batched == 2, "constructions=%d" % n_batched)

    _CountingClassifier.constructions = 0
    _bcreal.BRepClass_FaceClassifier = _CountingClassifier
    try:
        for e, i in items:
            _scalar_reference(e, fa, fb, i, **kw)
        n_scalar = _CountingClassifier.constructions
    finally:
        _bcreal.BRepClass_FaceClassifier = _real_classifier
    check("reference builds more classifiers", n_scalar > n_batched,
          "scalar=%d batched=%d" % (n_scalar, n_batched))

    # 3. _point_segment_distance matches hand-computed float arithmetic
    # to within a couple of ULPs (NumPy dot/norm may order the last
    # rounding differently; decisions compare against >= 1e-9).
    rng = np.random.default_rng(11)
    worst_rel = 0.0
    for _ in range(2000):
        p = rng.normal(size=3) * 10
        a = rng.normal(size=3) * 10
        b = a + rng.normal(size=3) * 3
        got = ix_mod._point_segment_distance(p, a, b)
        ab = b - a
        d = float(ab @ ab)
        if d <= 1e-300:
            want = float(np.linalg.norm(p - a))
        else:
            t = min(1.0, max(0.0, float(((p - a) @ ab) / d)))
            want = float(np.linalg.norm(p - (a + t * ab)))
        denom = max(abs(want), 1e-300)
        worst_rel = max(worst_rel, abs(got - want) / denom)
    check("point_segment_distance within 2 ULPs of float reference",
          worst_rel <= 2e-15, "worst rel diff=%r" % worst_rel)

    print("overall: %s" % ("PASS" if PASS else "FAIL"))
    return 0 if PASS else 1


if __name__ == "__main__":
    raise SystemExit(main())
