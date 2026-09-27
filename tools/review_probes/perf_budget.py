#!/usr/bin/env python3
"""Performance budget enforcement for the B-rep Boolean kernel (G13, S9).

Pinned cases are timed best-of-3 and compared against docs/perf_budget.json.
Each entry carries measured_s and budget_multiplier (default 2.0); the
nightly job fails if any pinned case exceeds measured * multiplier.
With --perf, expensive-call counters are recorded beside the time and
max_calls ceilings are enforced too, so optimizations cannot silently
disappear (S9: work avoided, not just time saved).

Case classes (S9):
  everyday    boxes, box/cylinder, rotated parts (general regressions)
  large       64/256-hole plates minus slot (C6/C9/regional classification)
  nurbs       converted-sphere section (section verification/completeness)
  separated   1000+ faces per operand (S3 whole-operation fast path)
  session     10 tools against one prepared base (S6 prepared-model reuse)
  many        16/64 cutters via boolean_brep_many (S7 batch scaling)

Usage:
    python tools/review_probes/perf_budget.py --record [--perf]
    python tools/review_probes/perf_budget.py --check [--perf]
    python tools/review_probes/perf_budget.py --check --cases quick
"""
import argparse
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

BUDGET_PATH = os.path.join(os.path.dirname(__file__), "..", "..",
                           "docs", "perf_budget.json")

DEFAULT_MULTIPLIER = 2.0   # pre-S9 entries
S9_MULTIPLIER = 1.5        # S9 budget format per the speed plan

# (name, class, runner, builder_key, op)
CASES = [
    # everyday (G13 originals)
    ("box_union", "everyday", "brep", "two_boxes", "union"),
    ("box_difference", "everyday", "brep", "two_boxes", "difference"),
    ("box_intersection", "everyday", "brep", "two_boxes", "intersection"),
    ("cylinder_through_box", "everyday", "brep", "cyl_box", "union"),
    ("hole_plate_4", "everyday", "brep", "hole_plate_4", "difference"),
    ("rotated_box_union", "everyday", "brep", "rotated_boxes", "union"),
    # large untouched assembly (S9)
    ("plate64_minus_slot", "large", "brep", "plate64", "difference"),
    ("plate256_minus_slot", "large", "brep", "plate256", "difference"),
    # NURBS section (S9)
    ("nurbs_sphere_cap", "nurbs", "brep", "nurbs_sphere_cap", "difference"),
    # separated large models (S9)
    ("separated_1020_union", "separated", "brep", "separated_1020", "union"),
    # repeated session (S9)
    ("session_10_warm", "session", "session", "session_10", "difference"),
    # many tools (S9)
    ("many16", "many", "many", "many16", "difference"),
    ("many64", "many", "many", "many64", "difference"),
]

QUICK = {"box_union", "box_difference", "cylinder_through_box",
         "nurbs_sphere_cap", "many16"}


def _plate_cache(name):
    return os.path.abspath(os.path.join(
        os.path.dirname(__file__), "..", "..", "goals",
        "robust-b-rep-booleans-prototype", "hidden_files", name))


def _drill_plate(n_holes, cache_name):
    """Drilled plate, cached to a BREP file (build once, reuse)."""
    from OCP.TopoDS import TopoDS_Shape
    from OCP.BRepTools import BRepTools
    from OCP.BRep import BRep_Builder
    cache = _plate_cache(cache_name)
    if os.path.exists(cache):
        shp = TopoDS_Shape()
        BRepTools.Read_s(shp, cache, BRep_Builder())
        return shp
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
    from OCP.BRep import BRep_Builder as _BB
    from OCP.TopoDS import TopoDS_Compound
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
    n = int(math.sqrt(n_holes))
    plate = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 8, 8, 0.5).Shape()
    comp = TopoDS_Compound()
    b = _BB()
    b.MakeCompound(comp)
    step = 8.0 / n
    for i in range(n):
        for j in range(n):
            x = step / 2 + i * step
            y = step / 2 + j * step
            cyl = BRepPrimAPI_MakeCylinder(
                gp_Ax2(gp_Pnt(x, y, -0.1), gp_Dir(0, 0, 1)),
                0.15, 0.7).Shape()
            b.Add(comp, cyl)
    cut = BRepAlgoAPI_Cut(plate, comp)
    cut.Build()
    if not cut.IsDone():
        raise RuntimeError("OCCT plate drilling failed")
    drilled = cut.Shape()
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    BRepTools.Write_s(drilled, cache)
    return drilled


def build_case(key):
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt
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
    if key == "plate64":
        drilled = _drill_plate(64, "plate64.brep")
        slot = BRepPrimAPI_MakeBox(gp_Pnt(3.0, 3.5, -0.25),
                                  2.0, 1.0, 1.0).Shape()
        return drilled, slot
    if key == "plate256":
        drilled = _drill_plate(256, "plate256.brep")
        slot = BRepPrimAPI_MakeBox(gp_Pnt(3.0, 3.5, -0.25),
                                  2.0, 1.0, 1.0).Shape()
        return drilled, slot
    if key == "nurbs_sphere_cap":
        # Same geometry as tests/test_nurbs_boolean_end_to_end.py.
        from OCP.BRepPrimAPI import BRepPrimAPI_MakeSphere
        from OCP.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert
        sph = BRepBuilderAPI_NurbsConvert(
            BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), 1.0).Shape(),
            True).Shape()
        cutter = BRepPrimAPI_MakeBox(gp_Pnt(-2, -2, -2),
                                    gp_Pnt(2, 2, 0.9999)).Shape()
        return sph, cutter
    if key == "separated_1020":
        sys.path.insert(0, os.path.dirname(__file__))
        from bench_s3_separated import build_compound
        return build_compound(170, 0.0), build_compound(170, 50.0)
    raise KeyError(key)


def _hole_tool(i, j, n):
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt
    x = 0.5 + i * (8.0 / n)
    y = 0.5 + j * (8.0 / n)
    return BRepPrimAPI_MakeBox(gp_Pnt(x - 0.15, y - 0.15, -0.25),
                              0.3, 0.3, 1.0).Shape()


def _merge_counters(dst, src):
    for k, v in (src or {}).items():
        dst[k] = dst.get(k, 0) + v


def run_brep(key, op, collect_perf):
    from brepkernel.pipeline import boolean_brep
    a, b = build_case(key)
    t0 = time.perf_counter()
    # fast_path_shadow=False: the budget pins the pipeline's algorithmic
    # performance, not the development-time S3 shadow re-run (which costs
    # ~274s on the separated case and would drown the fast path it guards).
    _, report = boolean_brep(a, b, op, collect_perf=collect_perf,
                             fast_path_shadow=False)
    dt = time.perf_counter() - t0
    if not report.get("accepted"):
        raise RuntimeError(f"budget case {key} refused; refusing cases "
                           "cannot pin a timing budget")
    counters = (report.get("performance", {}).get("counters")
                if collect_perf else None)
    return dt, counters


def run_session(key, op, collect_perf):
    # 10 small tools against one prepared base: pins S6 reuse.
    from brepkernel.session import BooleanSession
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt
    plate = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 8, 8, 0.5).Shape()
    tools = [_hole_tool(i % 5, i // 5, 5) for i in range(10)]
    sess = BooleanSession(plate)
    counters = {}
    t0 = time.perf_counter()
    try:
        for t in tools:
            _, rep = sess.boolean(t, op, collect_perf=collect_perf)
            if collect_perf:
                _merge_counters(
                    counters, rep.get("performance", {}).get("counters"))
    finally:
        sess.close()
    return time.perf_counter() - t0, (counters or None)


def run_many(key, op, collect_perf):
    # N cutters in one optimized call: pins S7 batch scaling.
    from brepkernel.session import boolean_brep_many
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt
    n = {"many16": 4, "many64": 8}[key]
    plate = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 8, 8, 0.5).Shape()
    tools = [_hole_tool(i, j, n) for i in range(n) for j in range(n)]
    t0 = time.perf_counter()
    _, summary = boolean_brep_many(plate, tools, op, optimized=True,
                                   collect_perf=collect_perf)
    dt = time.perf_counter() - t0
    counters = {}
    if collect_perf:
        for step in summary.get("steps", []):
            _merge_counters(
                counters, step.get("performance", {}).get("counters"))
    return dt, (counters or None)


RUNNERS = {"brep": run_brep, "session": run_session, "many": run_many}


def time_case(name, runner, key, op, repeats=3, collect_perf=False):
    best = None
    last_counters = None
    for _ in range(repeats):
        dt, counters = RUNNERS[runner](key, op, collect_perf)
        if collect_perf:
            last_counters = counters
        best = dt if best is None else min(best, dt)
    return best, last_counters


def budget_of(entry):
    # S9 format: measured_s * budget_multiplier.  Pre-S9 entries carry
    # an absolute budget_s.
    if "budget_s" in entry:
        return entry["budget_s"]
    return entry["measured_s"] * entry.get("budget_multiplier",
                                           DEFAULT_MULTIPLIER)


def load_budgets():
    path = os.path.abspath(BUDGET_PATH)
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)


def cmd_record(args):
    budgets = {"cases": {}, "note": (
        "Best-of-3 seconds on the recording runner; budget = measured_s * "
        "budget_multiplier (1.5 for S9 cases, 2.0 for older entries). "
        "Nightly fails if a pinned case exceeds its budget. "
        "max_calls ceilings (with --perf) pin expensive-call counts.")}
    for name, cls, runner, key, op in CASES:
        if args.cases == "quick":
            if name not in QUICK:
                continue
        elif args.cases != "all" and name not in args.cases.split(","):
            continue
        print(f"timing {name} [{cls}] ...", flush=True)
        best, counters = time_case(name, runner, key, op,
                                   collect_perf=args.perf)
        mult = S9_MULTIPLIER if cls != "everyday" else DEFAULT_MULTIPLIER
        entry = {"class": cls, "op": op,
                 "measured_s": round(best, 3),
                 "budget_multiplier": mult,
                 "budget_s": round(mult * best, 3)}
        if args.perf and counters:
            entry["counters"] = counters
            entry["max_calls"] = {k: int(math.ceil(v * 1.5)) + 1
                                  for k, v in sorted(counters.items())
                                  if v > 0}
        budgets["cases"][name] = entry
        print(f"  {name}: best-of-3 {best:.3f}s -> budget "
              f"{mult * best:.3f}s", flush=True)
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
    for name, cls, runner, key, op in CASES:
        if args.cases == "quick":
            if name not in QUICK:
                continue
        elif args.cases != "all" and name not in args.cases.split(","):
            continue
        if name not in budgets:
            print(f"SKIP {name}: no recorded budget")
            continue
        budget = budget_of(budgets[name])
        print(f"checking {name} (budget {budget:.3f}s) ...", flush=True)
        best, counters = time_case(name, runner, key, op,
                                   collect_perf=args.perf)
        status = "OK" if best <= budget else "OVER BUDGET"
        print(f"  {name}: best-of-3 {best:.3f}s [{status}]")
        if best > budget:
            failures.append((name, best, budget))
        max_calls = budgets[name].get("max_calls") or {}
        if args.perf and counters and max_calls:
            over = [(k, counters.get(k, 0), lim)
                    for k, lim in max_calls.items()
                    if counters.get(k, 0) > lim]
            for k, got, lim in over:
                print(f"  {name}: counter {k} {got} > max_calls {lim} "
                      f"[OVER BUDGET]")
                failures.append((f"{name}:{k}", got, lim))
    if failures:
        print("BUDGET FAILURES:")
        for name, best, budget in failures:
            print(f"  {name}: {best:.3f} > {budget:.3f}")
        return 1
    print("all pinned cases within budget")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--record", action="store_true")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--cases", default="all",
                    help="'all', 'quick', or comma-separated case names")
    ap.add_argument("--perf", action="store_true",
                    help="also record/check S0 expensive-call counters "
                         "beside time")
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
