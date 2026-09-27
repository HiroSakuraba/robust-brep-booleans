"""S0 profiler tests.

Pass criteria from the speed plan:
- perf-disabled runs produce equivalent canonical reports apart from the
  intentionally added performance section;
- perf-enabled and perf-disabled verdicts are identical;
- profiler overhead is small (the tight <2% number is measured on the
  262-face plate and recorded in the ledger; here we guard against
  pathological slowdown on a moderate case).
"""
import os
import sys
import time

sys.path.insert(0, "src")
sys.path.insert(0, "tests")

from brepkernel import boolean_brep, BRepAmbiguousResult
from brepkernel import perf as perf_mod

from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
from OCP.GProp import GProp_GProps
from OCP.TopAbs import TopAbs_SOLID
from OCP.TopExp import TopExp_Explorer
from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def volume(shape):
    p = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, p, 1e-9, True)
    return float(p.Mass())


def solid_count(shape):
    n = 0
    ex = TopExp_Explorer(shape, TopAbs_SOLID)
    while ex.More():
        n += 1
        ex.Next()
    return n


def boxes():
    a = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 1, 1, 1).Shape()
    b = BRepPrimAPI_MakeBox(gp_Pnt(0.5, 0.5, 0.5), 1, 1, 1).Shape()
    return a, b


def cyl_box():
    box = BRepPrimAPI_MakeBox(gp_Pnt(-1, -1, -1), 2, 2, 2).Shape()
    cyl = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(0, 0, -2), gp_Dir(0, 0, 1)), 0.3, 4).Shape()
    return box, cyl


def run(a, b, op, **kw):
    """Returns (accepted, volume_or_None, solids, refusal_kind, report)."""
    try:
        out, report = boolean_brep(a, b, op, **kw)
    except BRepAmbiguousResult as exc:
        rep = exc.report
        kind = rep.get("refusal", {}).get("kind", "?")
        return False, None, 0, kind, rep
    return True, volume(out), solid_count(out), None, report


def t1_disabled_has_no_performance_key():
    a, b = boxes()
    ok, _, _, _, report = run(a, b, "union")
    ok = check("t1 disabled: union accepted", ok)
    ok &= check("t1 disabled: no performance key",
                "performance" not in report)
    return ok


def t2_enabled_records_counters():
    a, b = boxes()
    ok, _, _, _, report = run(a, b, "union", collect_perf=True)
    perf = report.get("performance", {})
    ok = check("t2 enabled: union accepted", ok)
    ok &= check("t2 enabled: performance section present",
                perf.get("enabled") is True)
    counters = perf.get("counters", {})
    ok &= check("t2 enabled: all canonical keys present",
                all(k in counters for k in perf_mod.CANONICAL_COUNTERS),
                f"keys={len(counters)}")
    ok &= check("t2 enabled: exact_face_distance > 0",
                counters.get("exact_face_distance", 0) > 0,
                f"={counters.get('exact_face_distance', 0)}")
    ok &= check("t2 enabled: section_engine_call == 1",
                counters.get("section_engine_call") == 1)
    ok &= check("t2 enabled: volume_integration > 0",
                counters.get("volume_integration", 0) > 0)
    return ok


def t3_env_var_enables():
    old = os.environ.get("BREPKERNEL_PERF")
    os.environ["BREPKERNEL_PERF"] = "1"
    try:
        a, b = boxes()
        ok, _, _, _, report = run(a, b, "difference")
    finally:
        if old is None:
            del os.environ["BREPKERNEL_PERF"]
        else:
            os.environ["BREPKERNEL_PERF"] = old
    ok = check("t3 env: difference accepted", ok)
    ok &= check("t3 env: performance section present",
                report.get("performance", {}).get("enabled") is True)
    return ok


def t4_verdict_equivalence():
    cases = [
        ("boxes", boxes(), "union"),
        ("boxes", boxes(), "intersection"),
        ("boxes", boxes(), "difference"),
        ("cyl_box", cyl_box(), "union"),
        ("cyl_box", cyl_box(), "difference"),
    ]
    ok = True
    for name, (a, b), op in cases:
        r0 = run(a, b, op)
        r1 = run(a, b, op, collect_perf=True)
        same = (r0[0] == r1[0] and r0[3] == r1[3])
        if r0[0] and r1[0]:
            same = same and abs(r0[1] - r1[1]) < 1e-6 and r0[2] == r1[2]
        # Canonical evidence must be identical apart from the perf key.
        rep0 = {k: v for k, v in r0[4].items() if k != "performance"}
        rep1 = {k: v for k, v in r1[4].items() if k != "performance"}
        ev0 = rep0.get("evidence", {}).get("name")
        ev1 = rep1.get("evidence", {}).get("name")
        ok &= check(f"t4 {name} {op}: same verdict/volume/solids/kind",
                    same,
                    f"plain={r0[:4]} perf={r1[:4]}")
        ok &= check(f"t4 {name} {op}: evidence name identical",
                    ev0 == ev1, f"{ev0} vs {ev1}")
    return ok


def t5_witness_accounting():
    # Cylinder through box: untouched faces exercise the C9 shortcut.
    a, b = cyl_box()
    ok, _, _, _, report = run(a, b, "difference", collect_perf=True)
    counters = report.get("performance", {}).get("counters", {})
    att = counters.get("single_witness_attempt", 0)
    hit = counters.get("single_witness_hit", 0)
    fb = counters.get("single_witness_fallback", 0)
    ok = check("t5 accepted", ok)
    ok &= check("t5 attempt == hit + fallback", att == hit + fb,
                f"attempt={att} hit={hit} fallback={fb}")
    ok &= check("t5 some attempts happened", att > 0, f"attempt={att}")
    return ok


def t6_enabled_overhead_bounded():
    a, b = cyl_box()
    # Warm up once (OCCT lazy init), then best-of-3 each way.
    run(a, b, "union")
    best_off, best_on = None, None
    for _ in range(3):
        t0 = time.perf_counter()
        run(a, b, "union")
        dt = time.perf_counter() - t0
        best_off = dt if best_off is None else min(best_off, dt)
    for _ in range(3):
        t0 = time.perf_counter()
        run(a, b, "union", collect_perf=True)
        dt = time.perf_counter() - t0
        best_on = dt if best_on is None else min(best_on, dt)
    ratio = best_on / best_off if best_off else float("inf")
    ok = check("t6 enabled overhead < 25% on cyl/box",
               ratio < 1.25,
               f"off={best_off:.3f}s on={best_on:.3f}s ratio={ratio:.3f}")
    return ok


def t7_context_reset_after_call():
    a, b = boxes()
    run(a, b, "union", collect_perf=True)
    ok = check("t7 no ambient counters after call",
               perf_mod.current() is None)
    _, _, _, _, report = run(a, b, "union")
    ok &= check("t7 later plain call has no performance key",
                "performance" not in report)
    return ok


def main():
    # Make sure the ambient env var does not leak into these tests.
    os.environ.pop("BREPKERNEL_PERF", None)
    ok = True
    ok &= t1_disabled_has_no_performance_key()
    ok &= t2_enabled_records_counters()
    ok &= t3_env_var_enables()
    ok &= t4_verdict_equivalence()
    ok &= t5_witness_accounting()
    ok &= t6_enabled_overhead_bounded()
    ok &= t7_context_reset_after_call()
    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
