#!/usr/bin/env python3
"""S8 decision-rule profile: stage breakdown of the 262-face plate benchmark.

Plan S8: parallelize only if independent face-pair narrow-phase work
(intersection stage) still exceeds ~35% of wall time.
"""
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from bench_plate256 import build_drilled, CACHE, count_faces  # noqa
from OCP.TopoDS import TopoDS_Shape
from OCP.BRepTools import BRepTools
from OCP.BRep import BRep_Builder


def main():
    if os.path.exists(CACHE):
        shp = TopoDS_Shape()
        BRepTools.Read_s(shp, CACHE, BRep_Builder())
        drilled = shp
    else:
        drilled = build_drilled()
    print("plate faces:", count_faces(drilled), flush=True)

    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt
    slot = BRepPrimAPI_MakeBox(gp_Pnt(3.0, 3.5, -0.25), 2.0, 1.0, 1.0).Shape()

    from brepkernel.pipeline import boolean_brep
    t0 = time.perf_counter()
    _shape, report = boolean_brep(drilled, slot, "difference")
    wall = time.perf_counter() - t0
    print(f"wall: {wall:.2f}s accepted={report['accepted']}", flush=True)

    tm = report.get("timings_ms", {})
    keys = sorted(tm)
    print("timing keys:", keys, flush=True)
    total = 0.0
    for k in keys:
        v = tm[k]
        if isinstance(v, (int, float)):
            total += v
            print(f"  {k}: {v/1000:.3f}s ({100*v/1000/wall:.1f}% of wall)",
                  flush=True)
    print(f"  sum: {total/1000:.3f}s", flush=True)

    stages = report.get("stages", {})
    inter = stages.get("intersection", {})
    print("intersection stage keys:", sorted(inter.keys()), flush=True)
    for k in ("n_pairs", "n_candidate_pairs", "n_section_calls",
              "section_calls", "pairs_considered"):
        if k in inter:
            print(f"  intersection.{k} = {inter[k]}", flush=True)


if __name__ == "__main__":
    main()
