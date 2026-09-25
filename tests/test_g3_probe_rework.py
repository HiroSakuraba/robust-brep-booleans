"""G3 rework: per-interval completeness matching + raw tolerance ceiling.

t1 (partial branch omission): truncate one verified section edge to half
its range -- only PART of an intersection branch is dropped, not a whole
loop -- and require the probe to refuse with SectionCompletenessMismatch.
The pre-rework rule accepted this via its approximation_gap blanket
accept (one near sample was enough), so this test FAILS on pre-rework
code and PASSES after the rework.

t2 (RawIntersectionToleranceTooLoose): artificially inflate the raw
curve's OCCT-reported tolerance and require the probe to refuse with
RawIntersectionToleranceTooLoose instead of widening the match window.
The pre-rework code folded 2*curve_tol into tol_i, so this test FAILS
on pre-rework code and PASSES after the rework.
"""
import dataclasses
import sys

sys.path.insert(0, "src")

import brepkernel.intersection as ix_mod
from brepkernel.intersection import (
    IntersectionError,
    _raw_intersector_completeness_probe,
    _run_section_engine,
    _verify_section_edge,
)
from brepkernel.step_ingest import index_shape

from OCP.BRep import BRep_Builder
from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.BRepBuilderAPI import (
    BRepBuilderAPI_MakeEdge,
    BRepBuilderAPI_MakeFace,
    BRepBuilderAPI_NurbsConvert,
)
from OCP.Geom import Geom_BSplineSurface, Geom_Plane
from OCP.TopAbs import TopAbs_FACE
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS, TopoDS_Shell
from OCP.gp import gp_Dir, gp_Pln, gp_Pnt

try:
    from OCP.TColgp import TColgp_Array2OfPnt
    from OCP.TColStd import TColStd_Array1OfReal, TColStd_Array1OfInteger
except ImportError:
    from OCP.collections import (
        Array2_gp_Pnt as TColgp_Array2OfPnt,
        Array1_double as TColStd_Array1OfReal,
        Array1_int as TColStd_Array1OfInteger,
    )


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def _shell(face):
    s = TopoDS_Shell()
    b = BRep_Builder()
    b.MakeShell(s)
    b.Add(s, face)
    return s


def _saddle():
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
    return Geom_BSplineSurface(P, K, K, M, M, p, p, False, False)


def _tilted_nurbs_plane():
    plane = Geom_Plane(gp_Pln(gp_Pnt(0, 0, 0.15), gp_Dir(0.2, 0, 1)))
    pf = BRepBuilderAPI_MakeFace(plane, -3, 3, -3, 3, 1e-7).Face()
    conv = BRepBuilderAPI_NurbsConvert(_shell(pf), True)
    assert conv.IsDone()
    ex = TopExp_Explorer(conv.Shape(), TopAbs_FACE)
    return TopoDS.Face(ex.Current())


def _pair():
    """One-branch NURBS saddle vs tilted NURBS plane pair, verified."""
    a = index_shape(
        _shell(BRepBuilderAPI_MakeFace(_saddle(), 1e-7).Face()))
    b = index_shape(_shell(_tilted_nurbs_plane()))
    fa, fb = a.faces[0], b.faces[0]
    edges, _ = _run_section_engine(
        fa, fb, approximation=True, fuzzy=0.0,
        parallel=False, use_obb=True)
    verified = [
        _verify_section_edge(
            e, fa, fb, i, base_tol=1e-7, chord_tol=None,
            tangent_sin_tol=1e-4, max_section_tol=None)
        for i, e in enumerate(edges)
    ]
    return fa, fb, verified


def t1_partial_branch_omission_refused():
    """Dropping half of one branch must refuse (not a whole loop)."""
    fa, fb, verified = _pair()
    if len(verified) != 1:
        return check("t1 setup has exactly one section edge", False,
                     f"edges={len(verified)}")
    rec = verified[0]
    c = BRepAdaptor_Curve(rec.edge)
    t0, t1 = c.FirstParameter(), c.LastParameter()
    half = BRepBuilderAPI_MakeEdge(
        c.Curve().Curve(), t0, 0.5 * (t0 + t1)).Edge()
    # The probe matches raw samples against rec.edge geometry; swap in
    # the 3D-only half edge (no p-curve surgery needed for the probe).
    partial = dataclasses.replace(rec, edge=half)
    try:
        _raw_intersector_completeness_probe(
            fa, fb, [partial], base_tol=1e-7, fuzzy=0.0, parallel=False)
    except IntersectionError as e:
        return check("t1 partial branch omission refused",
                     getattr(e, "kind", "") == "SectionCompletenessMismatch",
                     f"kind={getattr(e, 'kind', '?')}")
    return check("t1 partial branch omission refused", False,
                 "probe accepted a half-dropped branch")


def t2_raw_tolerance_too_loose_refused():
    """An artificially inflated raw tolerance must refuse loudly."""
    if not hasattr(ix_mod, "_raw_curve_tolerance"):
        return check("t2 inflated raw tolerance refused", False,
                     "probe has no raw-tolerance hook: inflated "
                     "tolerances are silently absorbed into the window")
    fa, fb, verified = _pair()
    if len(verified) != 1:
        return check("t2 setup has exactly one section edge", False,
                     f"edges={len(verified)}")
    real = ix_mod._raw_curve_tolerance
    ix_mod._raw_curve_tolerance = lambda ic: 1.0  # noqa: E731
    try:
        _raw_intersector_completeness_probe(
            fa, fb, verified, base_tol=1e-7, fuzzy=0.0, parallel=False)
    except IntersectionError as e:
        ok = getattr(e, "kind", "") == "RawIntersectionToleranceTooLoose"
        detail = f"kind={getattr(e, 'kind', '?')} message={str(e)[:100]}"
    else:
        ok, detail = False, "probe accepted an inflated raw tolerance"
    finally:
        ix_mod._raw_curve_tolerance = real
    return check("t2 inflated raw tolerance refused", ok, detail)


def main():
    ok = True
    ok &= t1_partial_branch_omission_refused()
    ok &= t2_raw_tolerance_too_loose_refused()
    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
