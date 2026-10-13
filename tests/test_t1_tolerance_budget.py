"""T1: the result may not widen tolerances beyond what the pipeline certified.

OCCT's sewing runs SameParameter, which can widen a tolerance to hide a
curve that is off its surface while the shape still passes BRepCheck.
The verification stage records the result's max tolerance against the
certified budget and refuses ResultToleranceGrowth above the limit.
"""
import pytest

import brepkernel.pipeline as pipeline
from brepkernel import BRepAmbiguousResult, boolean_brep
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeSphere
from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

BOX = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), gp_Pnt(2, 2, 2)).Shape()
CASES = [
    ("box_minus_sphere", BOX, BRepPrimAPI_MakeSphere(gp_Pnt(2, 2, 2), 1.2).Shape(), "difference"),
    ("through_hole", BOX, BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(1, 1, -1), gp_Dir(0, 0, 1)), 0.5, 4).Shape(), "difference"),
    ("cyl_T_union", BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -2), gp_Dir(0, 0, 1)), 1, 4).Shape(),
     BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-2, 0, 0), gp_Dir(1, 0, 0)), 0.6, 4).Shape(), "union"),
]


@pytest.mark.parametrize("name,a,b,op", CASES, ids=[c[0] for c in CASES])
def test_budget_is_recorded_and_respected(name, a, b, op):
    _out, rep = boolean_brep(a, b, op)
    v = rep["stages"]["verification"]
    assert v["certified_tolerance"] > 0
    assert v["result_max_tolerance"] > 0
    assert v["tolerance_ratio"] <= pipeline.RESULT_TOLERANCE_GROWTH_LIMIT
    # measured worst case is 1.05; anything near the limit deserves a look
    assert v["tolerance_ratio"] < 1.5, v


def test_growth_beyond_the_limit_refuses(monkeypatch):
    monkeypatch.setattr(pipeline, "RESULT_TOLERANCE_GROWTH_LIMIT", 1e-3)
    a, b, op = CASES[0][1:]
    with pytest.raises(BRepAmbiguousResult) as ei:
        boolean_brep(a, b, op, fast_paths=False)
    assert ei.value.report["refusal"]["kind"] == "ResultToleranceGrowth"
    assert ei.value.report["refusal"]["stage"] == "verification"
