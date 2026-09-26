"""C4/C5: volume routing and full-precision reported volumes.

- A two-box union moved about 2000 units and rotated took 242 s on main
  9e45a72 (Gauss-Kronrod stalling on rotated planar faces). It must now
  finish quickly.
- Analytic solids: routed volume equals the closed form to 1e-12
  relative, near and far from the origin.
- The reported result volume matches a full-precision integral.
"""
import math
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from brepkernel import boolean_brep  # noqa: E402
from brepkernel.assembly import _shape_volume, clear_volume_cache  # noqa: E402
from OCP.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert, BRepBuilderAPI_Transform  # noqa: E402
from OCP.BRepGProp import BRepGProp  # noqa: E402
from OCP.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCone,  # noqa: E402
                             BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeSphere,
                             BRepPrimAPI_MakeTorus)
from OCP.GProp import GProp_GProps  # noqa: E402
from OCP.gp import gp_Ax1, gp_Dir, gp_Pnt, gp_Trsf, gp_Vec  # noqa: E402


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def far(s):
    t = gp_Trsf()
    t.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0.3, -0.8, 0.5)), math.radians(251))
    t2 = gp_Trsf()
    t2.SetTranslation(gp_Vec(1000, 2000, -1500))
    t.Multiply(t2)
    return BRepBuilderAPI_Transform(s, t, True).Shape()


def main():
    ok = True
    a = far(BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), gp_Pnt(1, 1, 1)).Shape())
    b = far(BRepPrimAPI_MakeBox(gp_Pnt(1, .25, .25), gp_Pnt(2, .75, .75)).Shape())
    t0 = time.perf_counter()
    out, rep = boolean_brep(a, b, "union")
    dt = time.perf_counter() - t0
    ok &= check("far rotated shared-face union is fast", dt < 20.0, f"{dt:.2f}s (was 242 s)")
    ok &= check("far rotated shared-face union volume", abs(rep["stages"]["assembly"]["volume"] - 1.25) < 1e-12,
                f"{rep['stages']['assembly']['volume']!r}")
    solids = {
        "sphere": (BRepPrimAPI_MakeSphere(1.0).Shape(), 4 / 3 * math.pi),
        "cylinder": (BRepPrimAPI_MakeCylinder(1.0, 2.0).Shape(), 2 * math.pi),
        "cone": (BRepPrimAPI_MakeCone(1.0, 0.5, 2.0).Shape(), math.pi * 2 / 3 * 1.75),
        "torus": (BRepPrimAPI_MakeTorus(3.0, 1.0).Shape(), 6 * math.pi ** 2),
    }
    for name, (s, exact) in solids.items():
        for pose, shp in (("origin", s), ("far", far(s))):
            clear_volume_cache()
            v = _shape_volume(shp)
            ok &= check(f"{name} {pose} routed volume", abs(v - exact) <= 1e-12 * exact,
                        f"rel err {abs(v - exact) / exact:.1e}")
    na = BRepBuilderAPI_NurbsConvert(BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), 1).Shape(), True).Shape()
    nb = BRepBuilderAPI_NurbsConvert(BRepPrimAPI_MakeSphere(gp_Pnt(1, 0, 0), 1).Shape(), True).Shape()
    out, rep = boolean_brep(na, nb, "union")
    g = GProp_GProps()
    BRepGProp.VolumePropertiesGK_s(out, g, 1e-10, True, True, False, False, False)
    r = rep["stages"]["assembly"]["volume"]
    ok &= check("reported NURBS union volume is full precision",
                abs(r - g.Mass()) <= 1e-12 * g.Mass(), f"rel diff {abs(r - g.Mass()) / g.Mass():.1e}")
    print("ALL PASS" if ok else "FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
