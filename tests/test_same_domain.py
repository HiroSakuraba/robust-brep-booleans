"""Strict same-domain equivalence regressions."""
import sys

sys.path.insert(0, "src")

from brepkernel.same_domain import same_domain_shapes

from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCP.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeSphere
from OCP.TopAbs import TopAbs_FACE
from OCP.TopExp import TopExp_Explorer
from OCP.gp import gp_Pnt


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def face_count(shape):
    n = 0
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        n += 1
        ex.Next()
    return n


def t1_independent_boxes_match():
    a = BRepPrimAPI_MakeBox(1.0, 2.0, 3.0).Shape()
    b = BRepPrimAPI_MakeBox(1.0, 2.0, 3.0).Shape()
    assert not a.IsSame(b)
    r = same_domain_shapes(a, b)
    return check(
        "sd1 independent identical boxes",
        r.equivalent and len(r.matches) == 6,
        f"equivalent={r.equivalent} reason={r.reason} matches={len(r.matches)}")


def t2_shifted_box_rejected():
    a = BRepPrimAPI_MakeBox(
        gp_Pnt(0, 0, 0), gp_Pnt(1, 2, 3)).Shape()
    b = BRepPrimAPI_MakeBox(
        gp_Pnt(0.01, 0, 0), gp_Pnt(1.01, 2, 3)).Shape()
    r = same_domain_shapes(a, b)
    return check(
        "sd2 shifted box rejected without canonicalization",
        not r.equivalent
        and not r.canonicalized
        and r.reason == "model bounding boxes differ",
        f"reason={r.reason} canonicalized={r.canonicalized}")


def t3_reversed_same_tshape_rejected():
    a = BRepPrimAPI_MakeBox(1.0, 2.0, 3.0).Shape()
    b = a.Reversed()
    # This pins the reason pipeline must use IsEqual rather than IsSame.
    same_ignores_orientation = a.IsSame(b)
    equal_respects_orientation = a.IsEqual(b)
    r = same_domain_shapes(a, b)
    return check(
        "sd3 reversed same TShape rejected without canonicalization",
        same_ignores_orientation and not equal_respects_orientation
        and not r.equivalent
        and not r.canonicalized
        and r.reason == "global material orientation differs",
        f"IsSame={same_ignores_orientation} IsEqual={equal_respects_orientation} "
        f"reason={r.reason} canonicalized={r.canonicalized}")


def t4_independent_nurbs_spheres_match():
    a0 = BRepPrimAPI_MakeSphere(1.0).Shape()
    b0 = BRepPrimAPI_MakeSphere(1.0).Shape()
    a = BRepBuilderAPI_NurbsConvert(a0, True).Shape()
    b = BRepBuilderAPI_NurbsConvert(b0, True).Shape()
    assert not a.IsSame(b)
    r = same_domain_shapes(a, b)
    return check(
        "sd4 independent NURBS spheres",
        r.equivalent and len(r.matches) == 1,
        f"equivalent={r.equivalent} reason={r.reason} matches={len(r.matches)}")



def t5_different_face_decomposition_matches_after_canonicalization():
    a = BRepPrimAPI_MakeBox(2.0, 1.0, 1.0).Shape()
    left = BRepPrimAPI_MakeBox(
        gp_Pnt(0, 0, 0), gp_Pnt(1, 1, 1)).Shape()
    right = BRepPrimAPI_MakeBox(
        gp_Pnt(1, 0, 0), gp_Pnt(2, 1, 1)).Shape()
    f = BRepAlgoAPI_Fuse(left, right)
    f.Build()
    assert f.IsDone()
    b = f.Shape()

    ca, cb = face_count(a), face_count(b)
    r = same_domain_shapes(a, b)
    canon = r.canonical_b
    cb_after = face_count(b)
    return check(
        "sd5 different decomposition canonicalized on copy",
        ca != cb
        and r.equivalent
        and r.canonicalized
        and canon is not None
        and canon.faces_before == cb
        and canon.faces_after == ca
        and cb_after == cb,
        f"faces={ca}/{cb} original_after={cb_after} "
        f"equivalent={r.equivalent} canonicalized={r.canonicalized} "
        f"canonB={canon}")


def t6_same_bbox_volume_different_material_rejected():
    # A: cubes on the SW and NE corners of a 2x2x1 bounding box.
    a0 = BRepPrimAPI_MakeBox(
        gp_Pnt(0, 0, 0), gp_Pnt(1, 1, 1)).Shape()
    a1 = BRepPrimAPI_MakeBox(
        gp_Pnt(1, 1, 0), gp_Pnt(2, 2, 1)).Shape()
    fa = BRepAlgoAPI_Fuse(a0, a1)
    fa.Build()
    assert fa.IsDone()
    a = fa.Shape()

    # B: same total volume and same global bbox, but cubes occupy the other
    # two corners. Bbox/volume agreement must not be mistaken for equivalence.
    b0 = BRepPrimAPI_MakeBox(
        gp_Pnt(0, 1, 0), gp_Pnt(1, 2, 1)).Shape()
    b1 = BRepPrimAPI_MakeBox(
        gp_Pnt(1, 0, 0), gp_Pnt(2, 1, 1)).Shape()
    fb = BRepAlgoAPI_Fuse(b0, b1)
    fb.Build()
    assert fb.IsDone()
    b = fb.Shape()

    r = same_domain_shapes(a, b)
    return check(
        "sd6 same bbox+volume different material rejected",
        not r.equivalent,
        f"reason={r.reason} canonicalized={r.canonicalized} "
        f"candidate_counts={r.candidate_counts}")

def main():
    ok = True
    ok &= t1_independent_boxes_match()
    ok &= t2_shifted_box_rejected()
    ok &= t3_reversed_same_tshape_rejected()
    ok &= t4_independent_nurbs_spheres_match()
    ok &= t5_different_face_decomposition_matches_after_canonicalization()
    ok &= t6_same_bbox_volume_different_material_rejected()
    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
