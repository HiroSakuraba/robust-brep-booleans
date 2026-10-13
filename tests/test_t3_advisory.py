"""T3: opt-in, uncertified second opinions on refusals.

advisory_on_refusal=True attaches what OCCT's own Boolean and (if
installed) SMLib return, with a triage label. It never changes the
verdict and is absent unless asked for.
"""
import pytest

import brepkernel.usd.smlib_backend as sb
from brepkernel import BRepAmbiguousResult, boolean_brep
from brepkernel.advisory import second_opinions, triage
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

# Equal-radius crossing cylinders (Steinmetz solid, volume 16/3): a
# genuine tangency that brepkernel refuses by design.
CA = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -2), gp_Dir(0, 0, 1)), 1, 4).Shape()
CB = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-2, 0, 0), gp_Dir(1, 0, 0)), 1, 4).Shape()


def _refusal(**kw):
    with pytest.raises(BRepAmbiguousResult) as ei:
        boolean_brep(CA, CB, "intersection", **kw)
    return ei.value


def test_off_by_default():
    exc = _refusal()
    assert "advisory" not in exc.report
    assert not hasattr(exc, "advisory_shapes")


def test_occt_opinion_without_smlib(monkeypatch):
    monkeypatch.setattr(sb, "available", lambda: False)
    exc = _refusal(advisory_on_refusal=True)
    adv = exc.report["advisory"]
    assert adv["certified"] is False
    assert adv["engines"]["smlib"]["status"] == "unavailable"
    assert abs(adv["engines"]["occt"]["volume"] - 16 / 3) < 1e-6
    assert adv["triage"] == "unclear"
    assert set(exc.advisory_shapes) == {"occt"}
    assert exc.report["evidence"]["operation"]["params"]["advisory_on_refusal"] is True


@pytest.mark.skipif(not sb.available(), reason="usd-brep not installed")
def test_two_engines_agree_on_steinmetz():
    exc = _refusal(advisory_on_refusal=True)
    adv = exc.report["advisory"]
    assert adv["engines"]["smlib"]["status"] == "answered"
    assert adv["triage"] == "suspected_over_refusal"
    assert exc.report["refusal"]["kind"]          # verdict unchanged


def test_accept_path_untouched():
    a = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), gp_Pnt(2, 2, 2)).Shape()
    b = BRepPrimAPI_MakeBox(gp_Pnt(1, 1, 1), gp_Pnt(3, 3, 3)).Shape()
    _out, rep = boolean_brep(a, b, "union", advisory_on_refusal=True)
    assert "advisory" not in rep


def test_triage_thresholds():
    assert triage(None) == "unclear"
    assert triage(1e-9) == "suspected_over_refusal"
    assert triage(1e-5) == "unclear"
    assert triage(0.07) == "justified"


def test_bad_op_is_a_usage_error():
    with pytest.raises(ValueError):
        second_opinions(CA, CB, "xor")
