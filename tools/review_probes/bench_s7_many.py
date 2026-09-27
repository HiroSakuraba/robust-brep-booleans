#!/usr/bin/env python3
"""S7 benchmark: N-hole pattern, Phase 1 sequential reference.

Builds a clean 8x8x0.5 box plate and an NxN grid of through-hole box
tools, then runs boolean_brep_many (Phase 1: plain sequential loop).
Reports total time, per-op stage medians, and the final face count.
"""
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
sys.path.insert(0, os.path.dirname(__file__))

from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.gp import gp_Pnt
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_FACE

from brepkernel.session import boolean_brep_many


def count_faces(shape):
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    n = 0
    while ex.More():
        n += 1
        ex.Next()
    return n


def plain_plate():
    return BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 8.0, 8.0, 0.5).Shape()


def hole_tool(i, j, n):
    step = 8.0 / n
    x = step / 2 + i * step
    y = step / 2 + j * step
    return BRepPrimAPI_MakeBox(gp_Pnt(x - 0.15, y - 0.15, -0.25),
                              0.3, 0.3, 1.0).Shape()


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8  # n x n holes
    base = plain_plate()
    tools = [hole_tool(i, j, n) for i in range(n) for j in range(n)]
    print(f"{n*n} hole tools on plain plate", flush=True)

    t0 = time.perf_counter()
    out, summary = boolean_brep_many(base, tools, "difference")
    total = time.perf_counter() - t0

    assert summary["optimized"] is False
    assert all(r["accepted"] for r in summary["steps"])
    stages = {}
    for r in summary["steps"]:
        for k, v in r["timings_ms"].items():
            stages.setdefault(k, []).append(v)
    print(f"sequential reference: {total:.1f}s for {n*n} holes "
          f"({total/(n*n):.2f}s/hole)", flush=True)
    for k in sorted(stages):
        vals = sorted(stages[k])
        med = vals[len(vals) // 2]
        print(f"  per-op median {k}: {med:.1f}ms", flush=True)
    print(f"final faces: {count_faces(out)}", flush=True)


if __name__ == "__main__":
    main()
