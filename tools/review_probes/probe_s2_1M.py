"""S2 gate probe: 1M randomized queries, BVH set vs vector-scan set.

Pass criterion from the plan: the BVH candidate set must contain every
candidate the vector scan returns (missing = fail; extra allowed).  The
implementation targets exact set equality, which is asserted here.

Usage: PYTHONPATH=src:tests <venv>/bin/python tools/review_probes/probe_s2_1M.py
"""
import sys
import time

sys.path.insert(0, "src")
sys.path.insert(0, "tests")

import numpy as np

from brepkernel import spatial
from brepkernel.spatial import HybridBoxIndex, scan_box_indices

SIZES = [256, 512, 1024, 2048, 4096]
PER_SIZE = 200_000  # 5 * 200k = 1M total
CHUNK = 2000


def main():
    rng = np.random.default_rng(777)
    total = 0
    missing = 0
    extra = 0
    t0 = time.perf_counter()
    for n in SIZES:
        cen = rng.uniform(0.0, 10.0, size=(n, 3))
        half = rng.uniform(0.0, 0.5, size=(n, 3))
        boxes = np.concatenate([cen - half, cen + half], axis=1)
        boxes[: n // 10, 3:] = boxes[: n // 10, :3]  # degenerate zero-volume
        spatial.BVH_THRESHOLD = 0
        idx = HybridBoxIndex(boxes)
        assert idx.uses_bvh, f"expected BVH at N={n}"
        corners = boxes[rng.integers(0, n, size=(CHUNK,)), :3]  # on-boundary pts
        for _ in range(PER_SIZE // CHUNK):
            pts = rng.uniform(-1.0, 11.0, size=(CHUNK, 3))
            pts[: CHUNK // 10] = corners[: CHUNK // 10]
            radii = 10.0 ** rng.uniform(-9.0, -1.0, size=(CHUNK,))
            radii[: CHUNK // 20] = 0.0  # exact-zero radius
            for p, r in zip(pts, radii):
                a = scan_box_indices(boxes, p, r)
                b = idx.query_point(p, float(r))
                total += 1
                if len(a) != len(b) or not np.array_equal(a, b):
                    sa, sb = set(a.tolist()), set(b.tolist())
                    missing += len(sa - sb)
                    extra += len(sb - sa)
    dt = time.perf_counter() - t0
    print(f"queries={total} missing={missing} extra={extra} time={dt:.1f}s")
    print("PASS" if missing == 0 else "FAIL")


if __name__ == "__main__":
    main()
