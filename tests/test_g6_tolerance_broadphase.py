"""G6 regression: the broad phase is tolerance-aware.

A face with a large OCCT tolerance near another solid must not be silently
treated as disjoint. Box A gets one face raised to a 1e-3 tolerance; box B
sits 5e-4 away from that face - inside the tolerance band, so the contact
is ambiguous. The operation must REFUSE with a typed refusal (I1); it must
never accept the pair as disjoint.
"""
import sys

sys.path.insert(0, "src")

from brepkernel import boolean_brep, BRepAmbiguousResult
from brepkernel.step_ingest import (
    _shape_bbox,
    candidate_face_pairs,
    index_shape,
)

from OCP.Bnd import Bnd_Box
from OCP.BRep import BRep_Builder, BRep_Tool
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.TopAbs import TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Pnt

FACE_TOL = 1e-3
GAP = 5e-4


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def _face_bbox(face):
    b = Bnd_Box()
    BRepBndLib.AddOptimal_s(face, b, False, False)
    p0, p1 = b.CornerMin(), b.CornerMax()
    return (p0.X(), p0.Y(), p0.Z()), (p1.X(), p1.Y(), p1.Z())


def make_high_tol_pair():
    """Box A with one 1e-3-tolerance face; box B GAP away from that face."""
    a = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 1.0, 1.0, 1.0).Shape()
    target = None
    ex = TopExp_Explorer(a, TopAbs_FACE)
    while ex.More():
        face = TopoDS.Face(ex.Current())
        lo, hi = _face_bbox(face)
        if lo[0] >= 0.999 and hi[0] <= 1.001:
            target = face
            break
        ex.Next()
    assert target is not None, "no +X face found on box A"
    builder = BRep_Builder()
    builder.UpdateFace(target, FACE_TOL)
    ee = TopExp_Explorer(target, TopAbs_EDGE)
    while ee.More():
        builder.UpdateEdge(TopoDS.Edge(ee.Current()), FACE_TOL)
        ev = TopExp_Explorer(ee.Current(), TopAbs_VERTEX)
        while ev.More():
            builder.UpdateVertex(TopoDS.Vertex(ev.Current()), FACE_TOL)
            ev.Next()
        ee.Next()
    got = float(BRep_Tool.Tolerance_s(target))
    assert got >= FACE_TOL, f"face tolerance did not stick: {got}"
    b = BRepPrimAPI_MakeBox(gp_Pnt(1.0 + GAP, 0, 0), 1.0, 1.0, 1.0).Shape()
    return a, b


def t1_pair_is_candidate():
    """The high-tolerance face pair must survive the broad phase."""
    a, b = make_high_tol_pair()
    ma, mb = index_shape(a), index_shape(b)
    cands = candidate_face_pairs(ma, mb, pad=0.0)
    face_tols = [float(BRep_Tool.Tolerance_s(f.face)) for f in ma.faces]
    ok = check("g6 pair is a candidate even with zero pad",
               len(cands) >= 1,
               f"candidates={len(cands)}")
    ok &= check("g6 face tolerance visible to the kernel",
                max(face_tols) >= FACE_TOL,
                f"max_face_tol={max(face_tols):.3g}")
    return ok


def t2_union_refuses_typed():
    """Union must refuse with a typed refusal, never accept as disjoint.

    Refusal-kind choice, documented per the G6 task: the plan lists
    ambiguous_contact or NearCoincidentFaces, but OCCT's section engine
    merges the 5e-4 gap inside the face's 1e-3 tolerance and emits section
    edges, so the pair never reaches the near-contact classifier. The
    section verifier then refuses with SectionToleranceTooLoose because
    the face tolerance (1e-3) exceeds the acceptance ceiling
    (max(128*base_tol, 1e-8*scale) = 1.28e-5). That is the most accurate
    existing kind: the section evidence genuinely cannot be certified on
    a face this loose. UnresolvedContact would describe a later stage the
    pair never reaches.
    """
    a, b = make_high_tol_pair()
    try:
        out, report = boolean_brep(a, b, "union")
    except BRepAmbiguousResult as exc:
        refusal = exc.report.get("refusal", {})
        pad = exc.report.get("broadphase_pad")
        ok = check("g6 union refuses with typed kind",
                   refusal.get("kind") == "SectionToleranceTooLoose"
                   and refusal.get("stage") == "intersection",
                   f"refusal={refusal}")
        ok &= check("g6 pad recorded in report and covers the tolerance",
                    pad is not None and pad >= 2.0 * FACE_TOL,
                    f"broadphase_pad={pad}")
        return ok
    detail = (f"accepted={report.get('accepted')} "
              f"candidates={report['stages']['intersection'].get('candidate_face_pairs')}")
    return check("g6 union must not accept an ambiguous contact", False,
                 f"accepted a wrong/disjoint result: {detail}")


def main():
    ok = True
    ok &= t1_pair_is_candidate()
    ok &= t2_union_refuses_typed()
    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
