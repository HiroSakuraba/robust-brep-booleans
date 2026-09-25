"""Local freeform face-splitting regressions (24 Sept 2026).

Pins the next Tier B/C step: verified section p-curves may locally partition
affected faces, while unaffected or topology-ambiguous contacts do not trigger
speculative splitting.
"""

import sys
import numpy as np

sys.path.insert(0, "src")

from brepkernel.intersection import intersect_models
from brepkernel.split import split_models
from brepkernel.step_ingest import index_shape

from OCP.BRep import BRep_Builder
from OCP.BRepBuilderAPI import (
    BRepBuilderAPI_MakeFace,
    BRepBuilderAPI_NurbsConvert,
)
from OCP.BRepPrimAPI import (
    BRepPrimAPI_MakeBox,
    BRepPrimAPI_MakeSphere,
    BRepPrimAPI_MakeTorus,
)
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
    n = 4
    P = TColgp_Array2OfPnt(1, n, 1, n)
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            x = (i - 1) / (n - 1) * 2.0 - 1.0
            y = (j - 1) / (n - 1) * 2.0 - 1.0
            z = 0.3 * (x * x - y * y)
            P.SetValue(i, j, gp_Pnt(float(x), float(y), float(z)))
    K = TColStd_Array1OfReal(1, 2)
    M = TColStd_Array1OfInteger(1, 2)
    K.SetValue(1, 0.0); K.SetValue(2, 1.0)
    M.SetValue(1, 4); M.SetValue(2, 4)
    return Geom_BSplineSurface(P, K, K, M, M, 3, 3, False, False)


def _vertical_plane(x=0.0):
    s = Geom_Plane(gp_Pln(gp_Pnt(float(x), 0, 0), gp_Dir(1, 0, 0)))
    return BRepBuilderAPI_MakeFace(s, -2, 2, -2, 2, 1e-7).Face()


def _horizontal_plane(z=0.0):
    s = Geom_Plane(gp_Pln(gp_Pnt(0, 0, float(z)), gp_Dir(0, 0, 1)))
    return BRepBuilderAPI_MakeFace(s, -2, 2, -2, 2, 1e-7).Face()


def t1_only_affected_trimmed_face_splits():
    bs = _bspline()
    fa = BRepBuilderAPI_MakeFace(bs, 0.2, 0.8, 0.2, 0.8, 1e-7).Face()
    fb = _vertical_plane()
    a = index_shape(_shell(fa))
    b = index_shape(_shell(fb))
    ix = intersect_models(a, b, base_tol=1e-7, chord_tol=1e-5)
    sp = split_models(a, b, ix, base_tol=1e-7)

    ra = sp.faces_a[0]
    rb = sp.faces_b[0]
    ok = check("t1 affected NURBS face split",
               ra.status == "split" and len(ra.pieces) == 2,
               f"status={ra.status} pieces={len(ra.pieces)}")
    ok &= check("t1 parent area conserved",
                abs(ra.area_error)
                <= max(2e-6 * max(ra.area_before, 1.0), 1e-6),
                f"before={ra.area_before:.12g} "
                f"after={ra.area_after:.12g} err={ra.area_error:.3e}")
    ok &= check("t1 split provenance retained",
                all(p.parent_face_id == a.faces[0].face_id
                    for p in ra.pieces))
    # The isolated section edge ends inside the much larger plane face; it is
    # not a complete splitting contour there. Local splitter correctly leaves
    # that face as one piece. In a closed solid other adjacent section edges
    # later complete the contour.
    ok &= check("t1 incomplete plane contour not invented",
                len(rb.pieces) == 1,
                f"status={rb.status} pieces={len(rb.pieces)}")
    ok &= check("t1 no unresolved transverse contact",
                not sp.unresolved_contacts,
                f"unresolved={sp.unresolved_contacts}")
    return ok


def t2_far_faces_are_bit_identical_passthrough():
    bs = _bspline()
    fa = BRepBuilderAPI_MakeFace(bs, 1e-7).Face()
    fb = _vertical_plane(10.0)
    a = index_shape(_shell(fa))
    b = index_shape(_shell(fb))
    ix = intersect_models(a, b)
    sp = split_models(a, b, ix)
    ok = check("t2 zero split calls", sp.split_calls == 0,
               f"calls={sp.split_calls}")
    ok &= check("t2 A unchanged",
                sp.faces_a[0].pieces[0].unchanged
                and sp.faces_a[0].pieces[0].face.IsSame(fa))
    ok &= check("t2 B unchanged",
                sp.faces_b[0].pieces[0].unchanged
                and sp.faces_b[0].pieces[0].face.IsSame(fb))
    return ok


def t3_point_tangency_blocks_speculative_split():
    sphere = BRepPrimAPI_MakeSphere(1.0).Shape()
    plane = _horizontal_plane(1.0)
    a = index_shape(sphere)
    b = index_shape(_shell(plane))
    ix = intersect_models(a, b, broadphase_pad=1e-8,
                          contact_tol=1e-6)
    sp = split_models(a, b, ix)
    statuses = [x[2] for x in sp.unresolved_contacts]
    return check("t3 tangent remains unresolved",
                 "point_contact" in statuses and sp.split_calls == 0,
                 f"unresolved={sp.unresolved_contacts} "
                 f"split_calls={sp.split_calls}")


def t4_near_tangent_gap_blocks_speculative_split():
    sphere = BRepPrimAPI_MakeSphere(1.0).Shape()
    plane = _horizontal_plane(1.0 + 5e-6)
    a = index_shape(sphere)
    b = index_shape(_shell(plane))
    ix = intersect_models(a, b, broadphase_pad=1e-5,
                          contact_tol=1e-5)
    sp = split_models(a, b, ix)
    statuses = [x[2] for x in sp.unresolved_contacts]
    return check("t4 near gap remains unresolved",
                 "ambiguous_contact" in statuses and sp.split_calls == 0,
                 f"unresolved={sp.unresolved_contacts}")



def t5_torus_seam_routes_existing_boundary_operand_specifically():
    """A seam loop is not a new cut on the seam-side face.

    The torus/plane pair has two verified loops. One lies on the torus's
    existing periodic seam; the other is interior. Default splitting must use
    only the interior loop on the torus while using both loops on the cutter.
    """
    tor0 = BRepPrimAPI_MakeTorus(3.0, 1.0).Shape()
    conv = BRepBuilderAPI_NurbsConvert(tor0, True)
    assert conv.IsDone()
    torus = conv.Shape()
    cutter = BRepPrimAPI_MakeBox(
        gp_Pnt(-5.0, -5.0, -2.0),
        gp_Pnt(5.0, 5.0, 0.0)).Shape()

    a = index_shape(torus)
    b = index_shape(cutter)
    ix = intersect_models(a, b, base_tol=1e-7)
    curve_pairs = [p for p in ix.pairs if p.edges]
    if len(curve_pairs) != 1 or len(curve_pairs[0].edges) != 2:
        return check(
            "t5 setup has two torus section loops",
            False,
            f"pairs={[(p.face_a,p.face_b,len(p.edges)) for p in ix.pairs]}")

    risks = [e.risk_flags for e in curve_pairs[0].edges]
    seam_edges = sum("seam_on_a" in r for r in risks)
    ok = check(
        "t5 setup has exactly one torus seam loop",
        seam_edges == 1,
        f"risks={risks}")

    sp = split_models(a, b, ix, base_tol=1e-7)
    ra = [r for r in sp.faces_a if r.source_edges]
    rb = [r for r in sp.faces_b if r.source_edges]

    ok &= check(
        "t5 seam-side torus uses only interior loop",
        not sp.unresolved_contacts
        and len(ra) == 1
        and ra[0].source_edges == 1
        and len(ra[0].pieces) >= 2,
        f"A={[(r.status,len(r.pieces),r.source_edges,r.area_error) for r in ra]} "
        f"unresolved={sp.unresolved_contacts}")
    ok &= check(
        "t5 opposite cutter uses both loops",
        len(rb) == 1
        and rb[0].source_edges == 2
        and len(rb[0].pieces) >= 3,
        f"B={[(r.status,len(r.pieces),r.source_edges,r.area_error) for r in rb]}")
    return ok

def main():
    ok = True
    ok &= t1_only_affected_trimmed_face_splits()
    ok &= t2_far_faces_are_bit_identical_passthrough()
    ok &= t3_point_tangency_blocks_speculative_split()
    ok &= t4_near_tangent_gap_blocks_speculative_split()
    ok &= t5_torus_seam_routes_existing_boundary_operand_specifically()
    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
