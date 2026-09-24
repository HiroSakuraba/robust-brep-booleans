"""Freeform/NURBS acceleration tests (24 Sept 2026).

Pins the Tier B/C foundations:
- NumPy rational B-spline evaluation/derivatives agree with OCCT;
- positive-weight knot-span control hulls are conservative;
- local Gauss-Newton projection recovers known normal offsets;
- curvature planning refines high-curvature spans without being used as a
  certificate;
- span-pair broad phase rejects far geometry and culls partial overlap;
- STEP/B-rep indexing preserves Solid->Shell->Face provenance.
"""
import sys
import numpy as np

sys.path.insert(0, "src")

from brepkernel.freeform import (FreeformFaceAccel, NurbsPatchIndex,
                                 candidate_patch_pairs,
                                 nurbs_from_occt_face)
from brepkernel.step_ingest import index_shape, candidate_face_pairs

try:
    # OCP <= 7.x bindings
    from OCP.TColgp import TColgp_Array2OfPnt
    from OCP.TColStd import (TColStd_Array1OfReal,
                             TColStd_Array1OfInteger,
                             TColStd_Array2OfReal)
except ImportError:
    # OCP 8.x moved generated collection types into OCP.collections.
    from OCP.collections import (
        Array2_gp_Pnt as TColgp_Array2OfPnt,
        Array1_double as TColStd_Array1OfReal,
        Array1_int as TColStd_Array1OfInteger,
        Array2_double as TColStd_Array2OfReal,
    )
from OCP.gp import gp_Pnt, gp_Vec
from OCP.Geom import Geom_BSplineSurface
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.BRep import BRep_Builder
from OCP.TopoDS import TopoDS_Shell


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def make_surface(shift=0.0):
    n, deg = 7, 3
    P = TColgp_Array2OfPnt(1, n, 1, n)
    W = TColStd_Array2OfReal(1, n, 1, n)
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            x = (i - 1) / (n - 1) * 4.0 + shift
            y = (j - 1) / (n - 1) * 3.0
            z = 0.25 * np.sin(x * 1.2) * np.cos(y * 1.1) + 0.05 * x * y
            P.SetValue(i, j, gp_Pnt(float(x), float(y), float(z)))
            W.SetValue(i, j, float(1.0 + 0.08 * np.sin(i + j)))
    ks = [0.0, 0.25, 0.5, 0.75, 1.0]
    ms = [4, 1, 1, 1, 4]
    U = TColStd_Array1OfReal(1, len(ks))
    V = TColStd_Array1OfReal(1, len(ks))
    MU = TColStd_Array1OfInteger(1, len(ms))
    MV = TColStd_Array1OfInteger(1, len(ms))
    for k, x in enumerate(ks, 1):
        U.SetValue(k, x); V.SetValue(k, x)
    for k, x in enumerate(ms, 1):
        MU.SetValue(k, x); MV.SetValue(k, x)
    return Geom_BSplineSurface(P, W, U, V, MU, MV,
                               deg, deg, False, False)


def make_face(shift=0.0):
    bs = make_surface(shift)
    return bs, BRepBuilderAPI_MakeFace(bs, 1e-7).Face()


def make_shell(shift=0.0):
    _, face = make_face(shift)
    sh = TopoDS_Shell()
    b = BRep_Builder()
    b.MakeShell(sh)
    b.Add(sh, face)
    return sh


def t1_eval_matches_occt():
    bs, face = make_face()
    ns = nurbs_from_occt_face(face)
    rng = np.random.default_rng(20260924)
    emax = dmax = 0.0
    for _ in range(200):
        u, v = map(float, rng.random(2))
        S, Su, Sv = ns.eval_d1(u, v)
        P, Uv, Vv = gp_Pnt(), gp_Vec(), gp_Vec()
        bs.D1(u, v, P, Uv, Vv)
        Q = np.array([P.X(), P.Y(), P.Z()])
        Qu = np.array([Uv.X(), Uv.Y(), Uv.Z()])
        Qv = np.array([Vv.X(), Vv.Y(), Vv.Z()])
        emax = max(emax, float(np.linalg.norm(S - Q)))
        dmax = max(dmax, float(np.linalg.norm(Su - Qu)),
                   float(np.linalg.norm(Sv - Qv)))
    return check("NURBS evaluator == OCCT", emax < 1e-11 and dmax < 1e-10,
                 f"value={emax:.3e} d1={dmax:.3e}")


def t2_span_hulls_are_conservative():
    _, face = make_face()
    ns = nurbs_from_occt_face(face)
    idx = NurbsPatchIndex(ns)
    ok = check("16 knot-span records", len(idx.records) == 16,
               f"n={len(idx.records)}")
    rng = np.random.default_rng(7)
    worst = -1e300
    for r in idx.records:
        for _ in range(40):
            u = float(rng.uniform(r.u0, r.u1))
            v = float(rng.uniform(r.v0, r.v1))
            p = ns.value(u, v)
            worst = max(worst, float(np.max(np.maximum(r.lo - p,
                                                       p - r.hi))))
    ok &= check("control-hull AABB contains samples", worst <= 1e-12,
                f"worst_outside={worst:.3e}")
    return ok


def t3_projection_recovers_offset():
    _, face = make_face()
    a = FreeformFaceAccel.from_occt_face(face)
    ok = True
    for u, v in ((0.13, 0.37), (0.48, 0.62), (0.84, 0.22)):
        S, Su, Sv = a.nurbs.eval_d1(u, v)
        n = np.cross(Su, Sv)
        n /= np.linalg.norm(n)
        want = 0.003
        r = a.project_point(S + n * want)
        ok &= (abs(r.distance - want) < 1e-8 and r.trimmed_in is True
               and np.linalg.norm(np.array(r.uv) - [u, v]) < 1e-7)
    return check("local NURBS projector", ok)


def t4_curvature_plan_is_finite():
    _, face = make_face()
    a = FreeformFaceAccel.from_occt_face(face)
    plan = a.refinement_plan(1e-3, max_edge=0.5)
    ok = (len(plan) == len(a.index.records)
          and all(np.isfinite(r["priority"]) for r in plan)
          and all(0.0 <= r["recommended_edge"] <= 0.5 for r in plan))
    return check("curvature refinement plan", ok,
                 f"patches={len(plan)} top_priority={plan[0]['priority']:.3g}")


def t5_patch_pair_culling():
    _, fa = make_face(0.0)
    _, ff = make_face(10.0)
    _, fn = make_face(1.5)
    ia = NurbsPatchIndex(nurbs_from_occt_face(fa))
    iff = NurbsPatchIndex(nurbs_from_occt_face(ff))
    inn = NurbsPatchIndex(nurbs_from_occt_face(fn))
    far = candidate_patch_pairs(ia, iff)
    near = candidate_patch_pairs(ia, inn)
    naive = len(ia.records) * len(inn.records)
    ok = len(far) == 0 and 0 < len(near) < naive
    return check("knot-span broad phase", ok,
                 f"far={len(far)} near={len(near)} naive={naive}")


def t6_step_hierarchy_and_face_broadphase():
    m = index_shape(BRepPrimAPI_MakeBox(1, 2, 3).Shape())
    ok = (len(m.solids) == 1 and len(m.shells) == 1 and len(m.faces) == 6
          and all(f.solid_id == 0 and f.shell_id == 0 for f in m.faces))
    ok &= check("STEP/BRep hierarchy preserved", ok,
                f"solids={len(m.solids)} shells={len(m.shells)} "
                f"faces={len(m.faces)}")

    a = index_shape(make_shell(0.0))
    b = index_shape(make_shell(10.0))
    c = index_shape(make_shell(1.5))
    far = candidate_face_pairs(a, b)
    near = candidate_face_pairs(a, c)
    ok2 = (len(a.faces) == 1 and a.faces[0].freeform is not None
           and len(far) == 0 and len(near) == 1
           and near[0]["patch_pairs"] is not None
           and len(near[0]["patch_pairs"]) > 0)
    ok &= check("face broad phase + NURBS refinement", ok2,
                f"far={len(far)} near={len(near)}")
    return ok


def main():
    ok = True
    ok &= t1_eval_matches_occt()
    ok &= t2_span_hulls_are_conservative()
    ok &= t3_projection_recovers_offset()
    ok &= t4_curvature_plan_is_finite()
    ok &= t5_patch_pair_culling()
    ok &= t6_step_hierarchy_and_face_broadphase()
    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
