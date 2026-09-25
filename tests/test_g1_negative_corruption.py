"""G1 negative test: a corrupted patch keep/discard decision must be caught.

Corruption: monkeypatch brepkernel.assembly._decision_rule so the first
patch decision (operand A's only face, classified outside B) flips its keep
flag from True to False in a disjoint-sphere union. The corrupted pipeline
then returns sphere B alone: closed, manifold, B-rep valid, and inside the
pre-G1 volume bounds (result volume == max(va, vb)), so the pre-G1 checks
accept it silently. The G1 guards must catch it:

  - the independent Boolean arbiter reports kernel_errors > 0;
  - boolean_brep(..., crosscheck_ops=True) raises a typed refusal with
    kind=OperationIdentityFailed, because
    vol(AuB) + vol(AnB) != vol(A) + vol(B) for the corrupted union.
"""
import math
import sys

sys.path.insert(0, "src")
sys.path.insert(0, "tests")

import _arbiter
import brepkernel.assembly as _asm
from brepkernel import BRepAmbiguousResult, boolean_brep

from OCP.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert
from OCP.BRepPrimAPI import BRepPrimAPI_MakeSphere
from OCP.gp import gp_Pnt


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def make_pair():
    a0 = BRepPrimAPI_MakeSphere(gp_Pnt(-2, 0, 0), 1.0).Shape()
    b0 = BRepPrimAPI_MakeSphere(gp_Pnt(2, 0, 0), 1.0).Shape()
    a = BRepBuilderAPI_NurbsConvert(a0, True).Shape()
    b = BRepBuilderAPI_NurbsConvert(b0, True).Shape()
    return a, b


class FlipFirstKeep:
    """Context manager flipping the first _decision_rule keep flag."""

    def __init__(self):
        self.real = _asm._decision_rule
        self.calls = []

    def __enter__(self):
        real = self.real
        calls = self.calls

        def flipped(operation, operand, classification):
            keep, rev = real(operation, operand, classification)
            calls.append((operation, operand, classification, keep, rev))
            if len(calls) == 1:
                return (not keep), rev
            return keep, rev

        _asm._decision_rule = flipped
        return self

    def __exit__(self, *exc):
        _asm._decision_rule = self.real
        return False


def t1_pre_g1_checks_pass_silently():
    """Documents the gap: pre-G1 checks accept the corrupted result."""
    a, b = make_pair()
    with FlipFirstKeep() as flip:
        out, report = boolean_brep(a, b, "union")
    ver = report["stages"]["verification"]
    true_union = 8.0 * math.pi / 3.0
    ok = check(
        "g1 corrupted union accepted by pre-G1 checks",
        report["accepted"]
        and ver["brep_valid"]
        and ver["closed"]
        and ver["volume_bounds_ok"]
        and abs(ver["result_volume"] - true_union) > 1.0,
        f"accepted={report['accepted']} "
        f"volume={ver['result_volume']:.6f} true={true_union:.6f} "
        f"flipped={flip.calls[0] if flip.calls else None}")
    return ok


def t2_arbiter_catches_corruption():
    """The independent Boolean arbiter alone catches the corruption
    (crosscheck disabled)."""
    a, b = make_pair()
    with FlipFirstKeep():
        out, report = boolean_brep(a, b, "union")
    assert report["accepted"], "corrupted run must be accepted, not refused"
    res, dt = _arbiter.raw_audit(a, b, out, "union", seed=7)
    return check(
        "g1 arbiter catches corrupted keep decision",
        res["kernel_errors"] > 0,
        f"kernel_errors={res['kernel_errors']} checked={res['checked']} "
        f"audit_ms={dt * 1000.0:.0f}")


def t3_crosscheck_catches_corruption():
    """crosscheck_ops=True turns the corruption into a typed refusal."""
    a, b = make_pair()
    try:
        with FlipFirstKeep():
            boolean_brep(a, b, "union", crosscheck_ops=True)
    except BRepAmbiguousResult as e:
        refusal = e.report.get("refusal", {})
        cross = e.report.get("stages", {}).get("crosscheck", {})
        return check(
            "g1 crosscheck_ops catches corrupted keep decision",
            refusal.get("kind") == "OperationIdentityFailed"
            and refusal.get("stage") == "crosscheck",
            f"refusal={refusal} crosscheck={cross}")
    except TypeError as e:
        return check(
            "g1 crosscheck_ops catches corrupted keep decision",
            False,
            f"crosscheck_ops not implemented: {e}")
    return check(
        "g1 crosscheck_ops catches corrupted keep decision",
        False,
        "corrupted union accepted even with crosscheck_ops=True")


def main():
    ok = True
    ok &= t1_pre_g1_checks_pass_silently()
    ok &= t2_arbiter_catches_corruption()
    ok &= t3_crosscheck_catches_corruption()
    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
