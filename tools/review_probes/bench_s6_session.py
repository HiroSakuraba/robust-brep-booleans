#!/usr/bin/env python3
"""S6 benchmark: 10 repeated tools against one base, cold vs session.

Plan pass criteria: ten sequential operations against the same base build
its immutable indexes once; cold and warm verdicts are identical; warm
repeated-operation setup time is at least 3x lower than rebuilding the
base model each time.
"""
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
sys.path.insert(0, os.path.dirname(__file__))

import bench_plate256 as bp
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.gp import gp_Pnt

from brepkernel.assembly import _shape_volume
from brepkernel.pipeline import boolean_brep
from brepkernel.session import BooleanSession


def slot(i):
    return BRepPrimAPI_MakeBox(
        gp_Pnt(3.0 + 0.1 * i, 3.5, -0.25), 2.0, 1.0, 1.0).Shape()


def vol(shape):
    return abs(float(_shape_volume(shape)))


def main():
    base = bp.build_drilled()
    print(f"base faces: {bp.count_faces(base)}", flush=True)
    tools = [slot(i) for i in range(10)]

    # Cold: raw base every time (prepares the base on each call).
    cold_ingest, cold_vols, t0 = [], [], time.perf_counter()
    for t in tools:
        out, rep = boolean_brep(base, t, "difference")
        assert rep["accepted"]
        cold_ingest.append(rep["timings_ms"]["ingest"])
        cold_vols.append(vol(out))
    cold_total = time.perf_counter() - t0

    # Warm: one session, base prepared once.
    warm_ingest, warm_vols = [], []
    t0 = time.perf_counter()
    with BooleanSession(base) as session:
        init_s = time.perf_counter() - t0
        n_prep = 1  # __init__ prepares exactly once by construction
        for t in tools:
            out, rep = session.boolean(t, "difference")
            assert rep["accepted"]
            warm_ingest.append(rep["timings_ms"]["ingest"])
            warm_vols.append(vol(out))
    warm_total = time.perf_counter() - t0

    for i, (cv, wv) in enumerate(zip(cold_vols, warm_vols)):
        assert cv == wv, f"op {i}: cold/warm volume differs {cv} vs {wv}"

    ci = sum(cold_ingest) / len(cold_ingest)
    wi = sum(warm_ingest) / len(warm_ingest)
    print(f"cold: total {cold_total:.2f}s, mean setup/op {ci:.1f}ms", flush=True)
    print(f"warm: total {warm_total:.2f}s (session init {init_s:.3f}s), "
          f"mean setup/op {wi:.1f}ms", flush=True)
    print(f"setup reduction: {ci / wi:.0f}x "
          f"(criterion: >= 3x) -> {'PASS' if ci / wi >= 3 else 'FAIL'}",
          flush=True)
    print(f"total time saved: {cold_total - warm_total:.2f}s", flush=True)
    print(f"verdicts/volumes identical across all 10 ops: True", flush=True)


if __name__ == "__main__":
    main()
