"""Strict same-domain equivalence regressions."""
import sys

sys.path.insert(0, "src")

from brepkernel.same_domain import same_domain_shapes

from OCP.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeSphere
from OCP.gp import gp_Pnt


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


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
        "sd2 shifted box rejected",
        not r.equivalent,
        f"reason={r.reason}")


def t3_reversed_same_tshape_rejected():
    a = BRepPrimAPI_MakeBox(1.0, 2.0, 3.0).Shape()
    b = a.Reversed()
    # This pins the reason pipeline must use IsEqual rather than IsSame.
    same_ignores_orientation = a.IsSame(b)
    equal_respects_orientation = a.IsEqual(b)
    r = same_domain_shapes(a, b)
    return check(
        "sd3 reversed same TShape rejected",
        same_ignores_orientation and not equal_respects_orientation
        and not r.equivalent,
        f"IsSame={same_ignores_orientation} IsEqual={equal_respects_orientation} "
        f"reason={r.reason}")


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


def main():
    ok = True
    ok &= t1_independent_boxes_match()
    ok &= t2_shifted_box_rejected()
    ok &= t3_reversed_same_tshape_rejected()
    ok &= t4_independent_nurbs_spheres_match()
    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
