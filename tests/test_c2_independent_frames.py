"""C2: flush features built in their own frame on a rotated part.

The part is rotated by an OCCT transform; the tool is built directly in
the rotated frame from independently computed math (numpy Rodrigues).
Both describe the same mathematical frame, but at many angles the face
normals differ in the last bits. Main 9e45a72 refused all four features
at 25, 55, 90 and 135 degrees; they must be accepted with exact volumes.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from brepkernel import boolean_brep, BRepAmbiguousResult  # noqa: E402
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform  # noqa: E402
from OCP.BRepGProp import BRepGProp  # noqa: E402
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder  # noqa: E402
from OCP.GProp import GProp_GProps  # noqa: E402
from OCP.gp import gp_Ax1, gp_Ax2, gp_Dir, gp_Pnt, gp_Trsf  # noqa: E402


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def vol(s):
    g = GProp_GProps()
    BRepGProp.VolumePropertiesGK_s(s, g, 1e-10, True, True, False, False, False)
    return g.Mass()


def run(angle_deg):
    ang = math.radians(angle_deg)
    axis = np.array([1.0, 2.0, 3.0]) / math.sqrt(14.0)
    T = gp_Trsf()
    T.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(*axis)), ang)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    R = np.eye(3) + math.sin(ang) * K + (1 - math.cos(ang)) * K @ K

    def frame(origin):
        o, z, x = R @ np.array(origin, float), R @ np.array([0, 0, 1.0]), R @ np.array([1.0, 0, 0])
        return gp_Ax2(gp_Pnt(*o), gp_Dir(*z), gp_Dir(*x))

    part = BRepBuilderAPI_Transform(
        BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), gp_Pnt(2, 2, 2)).Shape(), T, True).Shape()
    cases = [
        ("boss", BRepPrimAPI_MakeCylinder(frame((1, 1, 2)), 0.5, 1).Shape(), "union", 8 + math.pi * .25),
        ("blind hole", BRepPrimAPI_MakeCylinder(frame((1, 1, 1)), 0.5, 1).Shape(), "difference", 8 - math.pi * .25),
        ("flush pocket", BRepPrimAPI_MakeBox(frame((.5, .5, 1)), 1, 1, 1).Shape(), "difference", 7.0),
        ("stacked block", BRepPrimAPI_MakeBox(frame((.5, .5, 2)), 1, 1, 1).Shape(), "union", 9.0),
    ]
    ok = True
    for name, tool, op, exact in cases:
        try:
            out, _ = boolean_brep(part, tool, op)
            v = vol(out)
            ok &= check(f"{angle_deg} deg {name}", abs(v - exact) <= 1e-7 * exact,
                        f"vol={v:.10g} exact={exact:.10g}")
        except BRepAmbiguousResult as exc:
            r = exc.report.get("refusal", {})
            ok &= check(f"{angle_deg} deg {name}", False, f"refused {r.get('kind')}")
    return ok


def main():
    ok = True
    for a in (25, 55, 90, 135):
        ok &= run(a)
    print("ALL PASS" if ok else "FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
