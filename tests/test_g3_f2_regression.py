"""G3 regression: F2 completeness-probe false alarm.

NURBS cone INTERSECT NURBS sphere. Before the G3 tolerance rework this case
refused with SectionCompletenessMismatch (raw=3, trimmed=2,
max_distance=9.41089e-05, tol=4e-05) even though the unmatched raw component
is the same branch as a verified section edge: an approximation gap, not a
missing branch. After G3 the operation must be accepted, agree with the OCCT
boolean oracle within 1e-6 relative, pass the independent membership audit
with zero errors, and record the accepted gap in the stage report.

This test must FAIL on the pre-G3 code (refusal) and PASS after the fix.
"""
import os
import sys

sys.path.insert(0, "src")
sys.path.insert(0, "tests")

import _arbiter

from brepkernel import BRepAmbiguousResult, boolean_brep

from OCP.BRep import BRep_Builder
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common
from OCP.BRepGProp import BRepGProp
from OCP.BRepTools import BRepTools
from OCP.GProp import GProp_GProps
from OCP.TopoDS import TopoDS_Shape

DATA = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "data", "review_20260925")


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def load(name):
    s = TopoDS_Shape()
    if not BRepTools.Read_s(
            s, os.path.join(DATA, name), BRep_Builder()):
        raise RuntimeError(f"could not read {name}")
    return s


def volume(shape):
    p = GProp_GProps()
    BRepGProp.VolumePropertiesGK_s(shape, p, 1e-10, True, True, False,
                                   False, False)
    return float(p.Mass())


def main():
    ok = True
    a = load("completeness_false_alarm_A.brep")
    b = load("completeness_false_alarm_B.brep")

    try:
        out, report = boolean_brep(a, b, "intersection")
        accepted = True
        refusal = {}
    except BRepAmbiguousResult as exc:
        accepted = False
        refusal = exc.report.get("refusal", {})
        out, report = None, exc.report

    ok &= check(
        "f2 cone/sphere intersection accepted",
        accepted,
        f"refusal kind={refusal.get('kind')} "
        f"message={refusal.get('message', '')[:120]}"
        if not accepted else "")
    if not accepted:
        print("\nSOME FAILURES (operation refused; remaining checks skipped)")
        return 1

    oracle = BRepAlgoAPI_Common(a, b)
    oracle.Build()
    ov = volume(oracle.Shape())
    kv = volume(out)
    rel = abs(kv - ov) / max(abs(ov), 1e-300)
    ok &= check(
        "f2 volume agrees with OCCT oracle within 1e-6 relative",
        oracle.IsDone() and ov > 0.0 and rel <= 1e-6,
        f"kernel={kv:.9f} occt={ov:.9f} rel_err={rel:.3g}")

    ok &= _arbiter.check_accepted(
        "f2", check, a, b, out, "intersection", n=300)[0]

    ix = report.get("stages", {}).get("intersection", {})
    comps = ix.get("completeness_components", [])
    # The component that the pre-G3 probe refused: its worst sample sits at
    # 9.41089e-05, above the old 4e-05 limit. It must now be recorded with
    # both coverage numbers and an accepting verdict.
    f2comp = [c for c in comps
              if c["max_sample_distance"] > c["match_tolerance"]]
    ok &= check(
        "f2 probe records per-component coverage numbers",
        len(f2comp) == 1 and len(comps) >= 1,
        f"components={len(comps)} above_old_limit={len(f2comp)}")
    for c in f2comp:
        ok &= check(
            "f2 former false alarm accepted by the coverage rule",
            c["verdict"] in ("matched", "approximation_gap")
            and c["coverage_within_tol"] >= 0.95
            and c["max_sample_distance"] <= c["separation_band"]
            and c["nearest_distance"] < c["separation_band"],
            f"verdict={c['verdict']} "
            f"coverage={c['coverage_within_tol']:.4f} "
            f"nearest={c['nearest_distance']:.6g} "
            f"max={c['max_sample_distance']:.6g} "
            f"tol_i={c['match_tolerance']:.6g} "
            f"band={c['separation_band']:.6g}")

    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
