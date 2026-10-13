"""T4: degenerate (pole) edges are not a self-touching boundary.

Box union a sphere centred on its corner: the remaining 7/8 sphere face
crosses the south pole twice, so the result carries two distinct
degenerated edges on one vertex. The C3 self-touch check counted them as a
coincident edge pair and refused NonManifoldResult. Points cannot make a
boundary touch itself along a curve; the result is an ordinary manifold
solid. C3's real edge-touch cases must keep refusing (test_c3_self_touch).
"""
import math

import pytest

from brepkernel import boolean_brep
from brepkernel.assembly import _shape_volume
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeSphere
from OCP.gp import gp_Ax1, gp_Dir, gp_Pnt, gp_Trsf

R = 1.2
SPHERE = 4.0 / 3.0 * math.pi * R ** 3
EXPECT = {"union": 8 + 7 * SPHERE / 8, "difference": 8 - SPHERE / 8,
          "intersection": SPHERE / 8}


def rotated(s, ang):
    if not ang:
        return s
    t = gp_Trsf()
    t.SetRotation(gp_Ax1(gp_Pnt(0.3, -0.2, 0.1), gp_Dir(1, 2, 3)), ang)
    return BRepBuilderAPI_Transform(s, t, True).Shape()


@pytest.mark.parametrize("op", ["union", "difference", "intersection"])
@pytest.mark.parametrize("ang", [0.0, 0.7])
def test_corner_sphere(op, ang):
    box = rotated(BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), gp_Pnt(2, 2, 2)).Shape(), ang)
    sph = rotated(BRepPrimAPI_MakeSphere(gp_Pnt(2, 2, 2), R).Shape(), ang)
    out, _rep = boolean_brep(box, sph, op)
    v = _shape_volume(out, 1e-12)
    assert abs(v - EXPECT[op]) <= 1e-8 * EXPECT[op], (v, EXPECT[op])
