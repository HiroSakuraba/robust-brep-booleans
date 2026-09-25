"""Verified freeform intersection tests (24 Sept 2026).

These tests pin the next Tier B/C invariant: an accepted section edge must
have a coherent 3D curve plus p-curves on both trimmed input faces, and
point/near-tangent contacts must not silently become "disjoint".
"""

import sys
import numpy as np

sys.path.insert(0, "src")

from brepkernel.intersection import (
    IntersectionError, intersect_models, section_face_pair,
)
from brepkernel.step_ingest import index_shape

from OCP.BRep import BRep_Builder
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCP.BRepPrimAPI import BRepPrimAPI_MakeSphere
from OCP.Geom import Geom_BSplineSurface, Geom_Plane
try:
    # OCP <= 7.x bindings
    from OCP.TColgp import TColgp_Array2OfPnt
    from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
except ImportError:
    # OCP 8.x generated collections
    from OCP.collections import (
        Array2_gp_Pnt as TColgp_Array2OfPnt,
        Array1_int as TColStd_Array1OfInteger,
        Array1_double as TColStd_Array1OfReal,
    )
from OCP.TopoDS import TopoDS_Shell
from OCP.gp import gp_Dir, gp_Pln, gp_Pnt


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def _shell(face):
    s = TopoDS_Shell()
    b = BRep_Builder()
    b.MakeShell(s)
    b.Add(s, face)
    return s


def _bspline():
    # One smooth cubic patch: x/y span [-1,1], z is a saddle-like graph.
    n = 4
    p = 3
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


def _vertical_plane(x=0.0):
    surf = Geom_Plane(gp_Pln(gp_Pnt(float(x), 0, 0), gp_Dir(1, 0, 0)))
    return BRepBuilderAPI_MakeFace(surf, -2, 2, -2, 2, 1e-7).Face()


def _horizontal_plane(z=0.0):
    surf = Geom_Plane(gp_Pln(gp_Pnt(0, 0, float(z)), gp_Dir(0, 0, 1)))
    return BRepBuilderAPI_MakeFace(surf, -2, 2, -2, 2, 1e-7).Face()


def t1_transverse_curve_has_verified_pcurves():
    bs = _bspline()
    fa = BRepBuilderAPI_MakeFace(bs, 1e-7).Face()
    fb = _vertical_plane()
    a = index_shape(_shell(fa))
    b = index_shape(_shell(fb))
    r = intersect_models(a, b, base_tol=1e-7, chord_tol=1e-5)
    ok = check("t1 one candidate/section call",
               r.candidate_pairs == 1 and r.section_calls == 1,
               f"candidates={r.candidate_pairs} calls={r.section_calls}")
    ok &= check("t1 section curve verified",
                r.verified_edges >= 1 and r.pairs[0].status == "curve",
                f"edges={r.verified_edges} status={r.pairs[0].status}")
    for e in r.pairs[0].edges:
        ok &= check("t1 pcurves/trim coherent",
                    e.trim_ok and len(e.uv_a) == len(e.parameters)
                    and len(e.uv_b) == len(e.parameters))
        ok &= check("t1 3d/surface error bounded",
                    max(e.max_surface_error_a,
                        e.max_surface_error_b,
                        e.max_cross_surface_error) <= e.verify_tolerance,
                    f"tol={e.verify_tolerance:.3e} "
                    f"errA={e.max_surface_error_a:.3e} "
                    f"errB={e.max_surface_error_b:.3e}")
        ok &= check("t1 transverse not tangent",
                    e.min_transversality > 0.1,
                    f"sin(theta)_min={e.min_transversality:.6g}")
    return ok


def t2_trimmed_face_limits_section_pcurve():
    bs = _bspline()
    # The x=0 plane corresponds to approximately u=.5. Trimming v to
    # [.2,.8] must shorten the section edge; accepting the infinite/
    # untrimmed surface curve would be wrong.
    fa = BRepBuilderAPI_MakeFace(bs, 0.2, 0.8, 0.2, 0.8, 1e-7).Face()
    fb = _vertical_plane()
    a = index_shape(_shell(fa))
    b = index_shape(_shell(fb))
    r = intersect_models(a, b, base_tol=1e-7, chord_tol=1e-5)
    ok = r.verified_edges >= 1
    for e in r.pairs[0].edges:
        ok &= bool(np.all(e.uv_a[:, 0] >= 0.2 - e.verify_tolerance)
                   and np.all(e.uv_a[:, 0] <= 0.8 + e.verify_tolerance)
                   and np.all(e.uv_a[:, 1] >= 0.2 - e.verify_tolerance)
                   and np.all(e.uv_a[:, 1] <= 0.8 + e.verify_tolerance)
                   and e.trim_ok)
    return check("t2 section respects trim domain", ok,
                 f"edges={r.verified_edges}")


def t3_far_faces_cost_zero_section_calls():
    bs = _bspline()
    fa = BRepBuilderAPI_MakeFace(bs, 1e-7).Face()
    fb = _vertical_plane(10.0)
    a = index_shape(_shell(fa))
    b = index_shape(_shell(fb))
    r = intersect_models(a, b)
    return check("t3 broadphase skips far faces",
                 r.candidate_pairs == 0 and r.section_calls == 0,
                 f"candidates={r.candidate_pairs} calls={r.section_calls}")


def t4_tangent_contact_not_disjoint():
    # Sphere/plane tangency is zero-dimensional: Section returns a vertex,
    # not a curve. The kernel must preserve that as topology-sensitive.
    sphere = BRepPrimAPI_MakeSphere(1.0).Shape()
    plane = _horizontal_plane(1.0)
    a = index_shape(sphere)
    b = index_shape(_shell(plane))
    r = intersect_models(a, b, broadphase_pad=1e-8,
                         contact_tol=1e-6)
    statuses = [x.status for x in r.pairs]
    return check("t4 tangent contact survives",
                 "point_contact" in statuses,
                 f"statuses={statuses}")


def t5_near_tangent_gap_is_ambiguous_not_disjoint():
    # Above the sphere by 5e-6: no section curve, but deliberately inside
    # the caller's contact band. This must escalate rather than be culled.
    sphere = BRepPrimAPI_MakeSphere(1.0).Shape()
    plane = _horizontal_plane(1.0 + 5e-6)
    a = index_shape(sphere)
    b = index_shape(_shell(plane))
    r = intersect_models(a, b, broadphase_pad=1e-5,
                         contact_tol=1e-5)
    statuses = [x.status for x in r.pairs]
    return check("t5 near contact is ambiguous",
                 "ambiguous_contact" in statuses
                 and r.ambiguous_contacts >= 1,
                 f"statuses={statuses}")



def t6_loose_section_tolerance_refuses():
    """The verifier may not widen itself to an arbitrarily loose OCCT edge.

    An unrealistically tight caller ceiling is used to pin the refusal path:
    the ordinary section succeeds, but the same geometry must refuse when the
    allowed section/face tolerance is below OCCT's B-rep tolerance.
    """
    bs = _bspline()
    fa = BRepBuilderAPI_MakeFace(bs, 1e-7).Face()
    fb = _vertical_plane()
    a = index_shape(_shell(fa))
    b = index_shape(_shell(fb))
    try:
        intersect_models(
            a, b, base_tol=1e-7, chord_tol=1e-5,
            max_section_tol=1e-12)
    except IntersectionError as e:
        return check(
            "t6 loose section tolerance refuses",
            getattr(e, "kind", "") == "SectionToleranceTooLoose",
            f"kind={getattr(e, 'kind', '?')} message={e}")
    return check("t6 loose section tolerance refuses", False, "no refusal")

def main():
    ok = True
    ok &= t1_transverse_curve_has_verified_pcurves()
    ok &= t2_trimmed_face_limits_section_pcurve()
    ok &= t3_far_faces_cost_zero_section_calls()
    ok &= t4_tangent_contact_not_disjoint()
    ok &= t5_near_tangent_gap_is_ambiguous_not_disjoint()
    ok &= t6_loose_section_tolerance_refuses()
    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
