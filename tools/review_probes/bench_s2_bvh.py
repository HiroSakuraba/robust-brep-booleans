"""S2 crossover benchmark: vector scan vs BVH query time across box counts.

For each N in SIZES, builds one HybridBoxIndex forced onto the vector path
and one forced onto the BVH path over the same boxes, then reports the
median per-query time over Q random point/radius queries (best of 3 rounds).
Also sweeps BVH_LEAF_SIZE at N=4096 and records one-time build cost.

Usage: PYTHONPATH=src:tests <venv>/bin/python tools/review_probes/bench_s2_bvh.py
"""
import sys
import time

sys.path.insert(0, "src")
sys.path.insert(0, "tests")

import numpy as np

from brepkernel import spatial
from brepkernel.spatial import HybridBoxIndex, scan_box_indices

SIZES = [32, 64, 128, 256, 512, 1024, 4096]
LEAF_SWEEP = [8, 16, 32]
Q = 2000
ROUNDS = 3


def random_boxes(n, rng):
    cen = rng.uniform(0.0, 10.0, size=(n, 3))
    half = rng.uniform(0.01, 0.5, size=(n, 3))
    return np.concatenate([cen - half, cen + half], axis=1)


def random_queries(rng):
    pts = rng.uniform(-1.0, 11.0, size=(Q, 3))
    # Production-like radii: cap = 20*tol with tol in [1e-7, 1e-4], i.e.
    # essentially point queries with a tiny epsilon; plus a large-radius tail.
    tiny = 10.0 ** rng.uniform(-7.0, -2.7, size=(4 * Q // 5,))
    large = rng.uniform(2e-3, 0.1, size=(Q - 4 * Q // 5,))
    radii = np.concatenate([tiny, large])
    rng.shuffle(radii)
    return pts, radii


def median_query_time(idx, pts, radii):
    best = float("inf")
    for _ in range(ROUNDS):
        t0 = time.perf_counter()
        for p, r in zip(pts, radii):
            idx.query_point(p, r)
        dt = (time.perf_counter() - t0) / Q
        best = min(best, dt)
    return best


def main():
    rng = np.random.default_rng(20260926)
    print(f"{'N':>6} {'vector_us':>10} {'bvh_us':>10} {'speedup':>8}  build_ms(bvh)")
    results = {}
    for n in SIZES:
        boxes = random_boxes(n, rng)
        pts, radii = random_queries(rng)
        # correctness spot check on this workload
        spatial.BVH_THRESHOLD = 10 ** 9
        iv = HybridBoxIndex(boxes)
        spatial.BVH_THRESHOLD = 0
        ib = HybridBoxIndex(boxes)
        assert not iv.uses_bvh and ib.uses_bvh
        for p, r in zip(pts[:200], radii[:200]):
            a = iv.query_point(p, r)
            b = ib.query_point(p, r)
            assert np.array_equal(a, b), f"mismatch at N={n}"
        tv = median_query_time(iv, pts, radii)
        t0 = time.perf_counter()
        spatial.BVH_THRESHOLD = 0
        ib2 = HybridBoxIndex(boxes)
        build_ms = (time.perf_counter() - t0) * 1e3
        tb = median_query_time(ib2, pts, radii)
        results[n] = (tv, tb)
        print(f"{n:>6} {tv * 1e6:>10.2f} {tb * 1e6:>10.2f} {tv / tb:>8.2f}x  {build_ms:>8.2f}")

    print("\nleaf-size sweep at N=4096 (BVH query us):")
    boxes = random_boxes(4096, rng)
    pts, radii = random_queries(rng)
    for leaf in LEAF_SWEEP:
        spatial.BVH_LEAF_SIZE = leaf
        spatial.BVH_THRESHOLD = 0
        ib = HybridBoxIndex(boxes)
        tb = median_query_time(ib, pts, radii)
        print(f"  leaf={leaf:<3} {tb * 1e6:.2f} us")
    spatial.BVH_LEAF_SIZE = 16

    tv, tb = results[4096]
    print(f"\n4096-box speedup: {tv / tb:.2f}x (gate needs >= 3x or BVH stays disabled)")


if __name__ == "__main__":
    main()
