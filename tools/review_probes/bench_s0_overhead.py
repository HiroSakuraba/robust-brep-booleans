#!/usr/bin/env python3
"""S0 overhead measurement: 262-face plate minus slot, perf off vs on.

Best-of-3 each way, mirroring tools/review_probes/bench_plate256.py.
Also prints the counter snapshot from the perf-enabled run.
"""
import os
import sys
import time

# Worktree src first: the venv's editable brepkernel points at a stale tree.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

PLATE = os.path.join(
    os.path.dirname(__file__), "..", "..", "goals",
    "robust-b-rep-booleans-prototype", "hidden_files", "plate256.brep")


def load_shapes():
    from OCP.BRepTools import BRepTools
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Shape
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt
    plate = TopoDS_Shape()
    ok = BRepTools.Read_s(plate, PLATE, BRep_Builder())
    assert ok, f"cannot read {PLATE}"
    slot = BRepPrimAPI_MakeBox(gp_Pnt(3.0, 3.5, -0.25), 2.0, 1.0, 1.0).Shape()
    return plate, slot


def main():
    from brepkernel.pipeline import boolean_brep
    plate, slot = load_shapes()
    # Warm up.
    boolean_brep(plate, slot, "difference")

    def run(collect):
        t0 = time.perf_counter()
        out, report = boolean_brep(plate, slot, "difference",
                                   collect_perf=collect)
        dt = time.perf_counter() - t0
        return dt, report

    off, on = [], []
    on_reports = []
    # Interleaved to control for machine-load drift.
    for _ in range(3):
        dt, _ = run(False)
        off.append(dt)
        dt, rep = run(True)
        on.append(dt)
        on_reports.append(rep)
    best_off, best_on = min(off), min(on)
    print(f"perf off best-of-3: {best_off:.2f}s  samples="
          f"{[f'{x:.2f}' for x in sorted(off)]}")
    print(f"perf on  best-of-3: {best_on:.2f}s  samples="
          f"{[f'{x:.2f}' for x in sorted(on)]}")
    print(f"enabled overhead: {(best_on - best_off) / best_off * 100:.2f}%")
    print("counters (fastest perf-on run):")
    counters = on_reports[on.index(best_on)]["performance"]["counters"]
    for k, v in counters.items():
        if v:
            print(f"  {k}: {v}")


if __name__ == "__main__":
    raise SystemExit(main())
