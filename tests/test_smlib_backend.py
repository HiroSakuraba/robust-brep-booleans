"""Optional second engine: NVIDIA SMLib through usd-brep.

Tests marked needs_smlib skip when usd-brep is not installed (the fast CI
matrix); the usd-brep CI job runs them. The others check the behaviour
without SMLib, so they run everywhere.
"""
import pytest

import brepkernel._occt_compat  # noqa: F401
import brepkernel.usd.smlib_backend as sb
from brepkernel import BRepAmbiguousResult, boolean_brep
from brepkernel.evidence import validate_evidence
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCP.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCP.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder,
                             BRepPrimAPI_MakeSphere, BRepPrimAPI_MakeTorus)
from OCP.TopAbs import TopAbs_EDGE
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

needs_smlib = pytest.mark.skipif(not sb.available(), reason="usd-brep not installed")


def box(x0, y0, z0, x1, y1, z1):
    return BRepPrimAPI_MakeBox(gp_Pnt(x0, y0, z0), gp_Pnt(x1, y1, z1)).Shape()


def sph(c, r):
    return BRepPrimAPI_MakeSphere(gp_Pnt(*c), r).Shape()


def cyl(p, d, r, h):
    return BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(*p), gp_Dir(*d)), r, h).Shape()


def fillet(s, r):
    f = BRepFilletAPI_MakeFillet(s)
    ex = TopExp_Explorer(s, TopAbs_EDGE)
    while ex.More():
        f.Add(r, TopoDS.Edge(ex.Current()))
        ex.Next()
    return f.Shape()


B = box(0, 0, 0, 2, 2, 2)
S = sph((2, 2, 2), 1.2)

AGREE_CASES = [
    ("box_minus_sphere", B, S, "difference"),
    ("through_hole", B, cyl((1, 1, -1), (0, 0, 1), 0.5, 4), "difference"),
    ("shared_face_union", B, box(2, 0, 0, 4, 2, 2), "union"),
    ("cyl_T_union", cyl((0, 0, -2), (0, 0, 1), 1, 4), cyl((-2, 0, 0), (1, 0, 0), 0.6, 4), "union"),
    ("fillet_minus_cyl", fillet(B, 0.3), cyl((1, 1, -1), (0, 0, 1), 0.4, 4), "difference"),
    ("torus_minus_box", BRepPrimAPI_MakeTorus(2.0, 0.5).Shape(), box(0, 0, -1, 3, 3, 1),
     "difference"),
    ("split_into_two", box(0, 0, 0, 4, 2, 2), box(1.5, -1, -1, 2.5, 3, 3), "difference"),
    ("disjoint_intersection", B, sph((5, 5, 5), 1), "intersection"),
]


@needs_smlib
@pytest.mark.parametrize("name,a,b,op", AGREE_CASES, ids=[c[0] for c in AGREE_CASES])
def test_cross_check_agrees(name, a, b, op):
    out, _rep = boolean_brep(a, b, op)
    cc = sb.cross_check(a, b, op, out)
    assert cc["status"] == "agree", cc


@needs_smlib
def test_wrong_result_is_caught():
    wrong, _ = boolean_brep(B, S, "intersection")
    cc = sb.cross_check(B, S, "difference", wrong)
    assert cc["status"] == "disagree"
    assert cc["comparison"]["points"]["disagreements"] > 100


@needs_smlib
def test_to_smlib_keeps_cavities():
    U = sb._load()
    hollow = BRepAlgoAPI_Cut(box(0, 0, 0, 4, 4, 4), sph((2, 2, 2), 1)).Shape()
    brep = sb.to_smlib(hollow)
    assert brep.is_manifold_solid()
    assert U.material_census(brep) == {"solid_count": 1, "void_count": 1}


@needs_smlib
def test_smlib_boolean_round_trip_volume():
    import math
    out = sb.smlib_boolean(B, S, "difference")
    exact = 8 - (4 / 3 * math.pi * 1.2 ** 3) / 8
    assert abs(sb._volume(out) - exact) < 1e-6 * exact


@needs_smlib
def test_pipeline_records_stage_and_certification():
    out, rep = boolean_brep(B, S, "difference", independent_engine="smlib")
    ie = rep["stages"]["independent_engine"]
    assert ie["status"] == "agree" and ie["engine"] == "smlib"
    cert = rep["evidence"]["certification"]["independent_engine"]
    assert cert["status"] == "agree"
    assert rep["evidence"]["operation"]["params"]["independent_engine"] == "smlib"
    assert validate_evidence(rep["evidence"]) == []


@needs_smlib
def test_pipeline_refuses_on_disagreement(monkeypatch):
    real = sb.smlib_boolean
    monkeypatch.setattr(sb, "smlib_boolean",
                        lambda a, b, op, **kw: real(a, b, "intersection", **kw))
    with pytest.raises(BRepAmbiguousResult) as ei:
        boolean_brep(B, S, "difference", independent_engine="smlib")
    rep = ei.value.report
    assert rep["refusal"]["kind"] == "IndependentEngineDisagreement"
    assert rep["evidence"]["certification"]["independent_engine"]["status"] == "disagree"


# --- behaviour without SMLib (runs everywhere) ------------------------------

def test_default_runs_do_not_mention_the_option():
    _out, rep = boolean_brep(B, S, "difference")
    assert "independent_engine" not in rep["stages"]
    assert "independent_engine" not in rep["evidence"]["operation"]["params"]
    assert "independent_engine" not in rep["evidence"]["certification"]


def test_unavailable_is_recorded_not_refused(monkeypatch):
    monkeypatch.setattr(sb, "available", lambda: False)
    _out, rep = boolean_brep(B, S, "difference", independent_engine="smlib")
    assert rep["stages"]["independent_engine"]["status"] == "unavailable"


def test_unavailable_refuses_when_strict(monkeypatch):
    monkeypatch.setattr(sb, "available", lambda: False)
    with pytest.raises(BRepAmbiguousResult) as ei:
        boolean_brep(B, S, "difference", independent_engine="smlib",
                     independent_engine_strict=True)
    assert ei.value.report["refusal"]["kind"] == "IndependentEngineInconclusive"


def test_unknown_engine_is_a_usage_error():
    with pytest.raises(ValueError, match="independent_engine"):
        boolean_brep(B, S, "difference", independent_engine="parasolid")
