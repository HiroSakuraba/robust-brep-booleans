#!/usr/bin/env python3
"""Performance budget enforcement for the B-rep Boolean kernel (G13).

Pinned cases are timed best-of-3 and compared against docs/perf_budget.json,
where each budget is 2x the measured time on the reference runner.  The
nightly job fails if any pinned case exceeds its budget.

Usage:
    python tools/review_probes/perf_budget.py --record   # write docs/perf_budget.json
    python tools/review_probes/perf_budget.py --check    # fail if over budget
    python tools/review_probes/perf_budget.py --check --cases quick  # subset
"""

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

BUDGET_PATH = os.path.join(os.path.dirname(__file__), "..", "..",
                           "docs", "perf_budget.json")

# Pinned cases: (name, builder_key, op).  Builders live in _perf_cases so
# this tool stays importable without OCCT until a case actually runs.
CASES = [
    ("box_union", "two_boxes", "union"),
    ("box_difference", "two_boxes", "difference"),
    ("box_intersection", "two_boxes", "intersection"),
    ("cylinder_through_box", "cyl_box", "union"),
    ("hole_plate_4", "hole_plate_4", "difference"),
    ("rotated_box_union", "rotated_boxes", "union"),
]


def build_case(key):
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt
    import math
    if key == "two_boxes":
        a = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 1, 1, 1).Shape()
        b = BRepPrimAPI_MakeBox(gp_Pnt(0.5, 0.5, 0.5), 1, 1, 1).Shape()
        return a, b
    if key == "cyl_box":
        box = BRepPrimAPI_MakeBox(gp_Pnt(-1, -1, -1), 2, 2, 2).Shape()
        cyl = BRepPrimAPI_MakeCylinder(
            gp_Ax2(gp_Pnt(0, 0, -2), gp_Dir(0, 0, 1)), 0.3, 4).Shape()
        return box, cyl
    if key == "hole_plate_4":
        from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
        plate = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 4, 4, 0.5).Shape()
        drilled = plate
        for ix in (1, 3):
            for iy in (1, 3):
                h = BRepPrimAPI_MakeCylinder(
                    gp_Ax2(gp_Pnt(ix, iy, -0.1), gp_Dir(0, 0, 1)),
                    0.25, 0.7).Shape()
                cut = BRepAlgoAPI_Cut(drilled, h)
                cut.Build()
                drilled = cut.Shape()
        tool = BRepPrimAPI_MakeBox(gp_Pnt(1.5, 1.5, -0.5), 1, 1, 1.5).Shape()
        return drilled, tool
    if key == "rotated_boxes":
        from OCP.gp import gp_Trsf, gp_Ax1
        from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
        a = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 1, 1, 1).Shape()
        t = gp_Trsf()
        t.SetRotation(gp_Ax1(gp_Pnt(0.5, 0.5, 0.5), gp_Dir(1, 1, 1)),
                      math.pi / 7)
        b = BRepBuilderAPI_Transform(
            BRepPrimAPI_MakeBox(gp_Pnt(0.4, 0.4, 0.4), 1, 1, 1).Shape(),
            t, True).Shape()
        return a, b
    raise KeyError(key)


def time_case(name, key, op, repeats=3):
    from brepkernel.pipeline import boolean_brep
    best = None
    for _ in range(repeats):
        a, b = build_case(key)
        t0 = time.perf_counter()
        boolean_brep(a, b, op)
        dt = time.perf_counter() - t0
        best = dt if best is None else min(best, dt)
    return best


def load_budgets():
    path = os.path.abspath(BUDGET_PATH)
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)


def cmd_record(args):
    budgets = {"cases": {}, "note": (
        "Best-of-3 seconds on the recording runner; budget = 2x measured. "
        "Nightly fails if a pinned case exceeds its budget.")}
    for name, key, op in CASES:
        if args.cases != "all" and name not in args.cases.split(","):
            continue
        print(f"timing {name} ...", flush=True)
        best = time_case(name, key, op)
        budgets["cases"][name] = {
            "op": op, "measured_s": round(best, 3),
            "budget_s": round(2.0 * best, 3),
        }
        print(f"  {name}: best-of-3 {best:.3f}s -> budget {2*best:.3f}s")
    path = os.path.abspath(BUDGET_PATH)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(budgets, f, indent=2)
        f.write("\n")
    print(f"wrote {path}")


def cmd_check(args):
    budgets = load_budgets().get("cases", {})
    if not budgets:
        print("no budgets recorded; run --record first")
        return 2
    failures = []
    for name, key, op in CASES:
        if args.cases != "all" and name not in args.cases.split(","):
            continue
        if name not in budgets:
            print(f"SKIP {name}: no recorded budget")
            continue
        budget = budgets[name]["budget_s"]
        print(f"checking {name} (budget {budget:.3f}s) ...", flush=True)
        best = time_case(name, key, op)
        status = "OK" if best <= budget else "OVER BUDGET"
        print(f"  {name}: best-of-3 {best:.3f}s [{status}]")
        if best > budget:
            failures.append((name, best, budget))
    if failures:
        print("BUDGET FAILURES:")
        for name, best, budget in failures:
            print(f"  {name}: {best:.3f}s > {budget:.3f}s")
        return 1
    print("all pinned cases within budget")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--record", action="store_true")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--cases", default="all",
                    help="'all' or comma-separated case names")
    args = ap.parse_args()
    if args.record:
        cmd_record(args)
    elif args.check:
        sys.exit(cmd_check(args))
    else:
        ap.print_help()
        sys.exit(2)


if __name__ == "__main__":
    main()
