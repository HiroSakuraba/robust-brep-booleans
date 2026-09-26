"""C3: a union that only touches along an edge must refuse in every pose.

Main 9e45a72 refused it axis-aligned (NonManifoldResult) but accepted it
after generic rotations as one shell whose boundary touches itself.
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from brepkernel import boolean_brep, BRepAmbiguousResult  # noqa: E402
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform  # noqa: E402
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox  # noqa: E402
from OCP.gp import gp_Ax1, gp_Dir, gp_Pnt, gp_Trsf, gp_Vec  # noqa: E402


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def rigid(axis, deg, tx=0.0, ty=0.0, tz=0.0):
    t = gp_Trsf()
    n = math.sqrt(sum(c * c for c in axis))
    t.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(*(c / n for c in axis))), math.radians(deg))
    t2 = gp_Trsf()
    t2.SetTranslation(gp_Vec(tx, ty, tz))
    t.Multiply(t2)
    return t


POSES = [("identity", gp_Trsf()), ("rot30", rigid((1, 2, 3), 30)), ("rotz90", rigid((0, 0, 1), 90)),
         ("r37", rigid((1, 2, 3), 37, 10, -7, 3)), ("r113", rigid((-2, 1, 0.5), 113, -50, 21, -13))]


def main():
    ok = True
    a0 = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), gp_Pnt(1, 1, 1)).Shape()
    b0 = BRepPrimAPI_MakeBox(gp_Pnt(1, 1, 0), gp_Pnt(2, 2, 1)).Shape()
    for name, T in POSES:
        a = BRepBuilderAPI_Transform(a0, T, True).Shape()
        b = BRepBuilderAPI_Transform(b0, T, True).Shape()
        try:
            boolean_brep(a, b, "union")
            ok &= check(f"{name} edge-touch union refuses", False, "accepted")
        except BRepAmbiguousResult as exc:
            kind = exc.report.get("refusal", {}).get("kind")
            ok &= check(f"{name} edge-touch union refuses", kind == "NonManifoldResult", kind)
        try:
            out, rep = boolean_brep(a, b, "union", allow_nonmanifold=True)
            ok &= check(f"{name} allow_nonmanifold returns it with a note", True,
                        str(rep["stages"]["assembly"].get("notes", ""))[:80])
        except BRepAmbiguousResult as exc:
            ok &= check(f"{name} allow_nonmanifold returns it with a note", False,
                        exc.report.get("refusal", {}).get("kind"))
    print("ALL PASS" if ok else "FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
