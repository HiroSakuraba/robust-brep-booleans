"""Performance budgets for representative boolean_brep() cases.

Runs each pinned case best-of-N and records seconds. Intended output is
docs/perf_budget.json with budget = 2x measured (best-of-3), committed by
G13 on the reference machine. PR0 lands the tool only; do NOT commit
machine-specific budgets yet.

Cases time boolean_brep() itself; shapes are built once and reused across
reps (documented assumption: the pipeline does not mutate its inputs).

Usage:
  PYTHONPATH=src python tools/review_probes/perf_budget.py --reps 3
  PYTHONPATH=src python tools/review_probes/perf_budget.py --cases boxes_union,box_minus_tilted_cyl --out /tmp/b.json
  PYTHONPATH=src python tools/review_probes/perf_budget.py --include-slow
"""
from __future__ import annotations

import argparse
import json
import math
import os
import platform
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from brepkernel import boolean_brep  # noqa: E402

from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder  # noqa: E402
from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt  # noqa: E402


def _box(x0, y0, z0, x1, y1, z1):
    return BRepPrimAPI_MakeBox(gp_Pnt(x0, y0, z0), gp_Pnt(x1, y1, z1)).Shape()


def _cyl(p, d, r, h):
    return BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(*p), gp_Dir(*d)), r, h).Shape()


def _perforated_plate(n):
    plate = _box(0, 0, 0, 10, 10, 1)
    holes = []
    for i in range(n):
        for j in range(n):
            holes.append(_cyl((0.5 + i, 0.5 + j, -0.5), (0, 0, 1), 0.3, 2.0))
    return plate, holes


def case_boxes_union():
    return (_box(0, 0, 0, 1, 1, 1), _box(0.5, 0.5, 0.5, 1.5, 1.5, 1.5), "union")


def case_box_minus_tilted_cyl():
    ang = math.radians(15)
    d = (math.sin(ang), 0.0, math.cos(ang))
    return (_box(0, 0, 0, 2, 2, 2), _cyl((1, 1, -1), d, 0.5, 4), "difference")


def case_plate_64_holes():
    plate, holes = _perforated_plate(8)
    return (plate, holes, "difference")


def case_plate_256_holes():
    plate, holes = _perforated_plate(16)
    return (plate, holes, "difference")


CASES = {
    "boxes_union": (case_boxes_union, False),
    "box_minus_tilted_cyl": (case_box_minus_tilted_cyl, False),
    "plate_64_holes": (case_plate_64_holes, False),
    "plate_256_holes": (case_plate_256_holes, True),
}


def run_case(name, build, reps):
    a, b, op = build()
    tools = b if isinstance(b, list) else None
    best = None
    for _ in range(reps):
        if tools is not None:
            # sequential multi-tool difference, timed as one boolean workload
            from brepkernel import boolean_brep as _bb
            t0 = time.perf_counter()
            shape = a
            for tool in tools:
                shape, _ = _bb(shape, tool, op)
        else:
            t0 = time.perf_counter()
            boolean_brep(a, b, op)
        dt = time.perf_counter() - t0
        best = dt if best is None else min(best, dt)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--cases", default=None,
                    help="comma-separated subset of: " + ",".join(CASES))
    ap.add_argument("--include-slow", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    names = args.cases.split(",") if args.cases else list(CASES)
    results = {}
    for name in names:
        build, slow = CASES[name]
        if slow and not args.include_slow:
            print("skip (slow): %s" % name)
            continue
        secs = run_case(name, build, args.reps)
        results[name] = {"seconds": round(secs, 3), "budget": round(2 * secs, 3)}
        print("%-22s best-of-%d: %.3fs  budget: %.3fs" %
              (name, args.reps, secs, 2 * secs), flush=True)

    import OCP
    payload = {
        "machine": platform.node(),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "ocp": getattr(OCP, "__version__", "unknown"),
        "reps": args.reps,
        "cases": results,
    }
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(payload, fh, indent=1)
        print("wrote", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
