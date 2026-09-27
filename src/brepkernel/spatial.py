"""Hybrid AABB proximity index (S2 of the speed plan).

One implementation for face-box and edge-box point/radius queries, replacing
the bare O(N) scan at large N.

Below ``BVH_THRESHOLD`` boxes the index *is* the vectorized NumPy scan: a
tree cannot beat it there.  At or above the threshold it builds a
median-split bounding-volume hierarchy over the boxes (no third-party
dependency).  Queries traverse the tree and apply the identical per-box test
at the leaves, so the returned candidate set is exactly the vector scan's
set -- a fortiori the conservative superset the plan requires.  ``BRepExtrema``
stays the exact test downstream.

The threshold is empirical; see ``tools/review_probes/bench_s2_bvh.py`` for
the crossover measurement that set it.

Pruning safety (why the tree can never miss a candidate the scan finds):
node bounds are exact componentwise min/max of the child boxes, and
correctly-rounded subtraction is monotone, so whenever a box passes the
per-box test ``q >= lo - r, q <= hi + r``, every ancestor node's expanded
union box passes the same test and is never pruned.
"""

from __future__ import annotations

import numpy as np

from .perf import count as _perf_count

# Build a BVH at/above this many boxes; below it the vector scan wins.
# Calibrated by tools/review_probes/bench_s2_bvh.py on the reference runner:
# with production-like radii the BVH first beats the scan at N ~= 32 and the
# lead grows monotonically (1.31x at 256, 1.72x at 512, 2.28x at 1024,
# 5.12x at 4096, clearing the plan's 3x gate at 4096).  The threshold sits
# at 256 for margin across runners; absolute scan times below that are
# already <= ~28 us per query.
BVH_THRESHOLD = 256

# Max boxes stored in one BVH leaf.  Leaf-size sweep at N=4096 on the
# reference runner: 8 -> 47.6 us, 16 -> 44.9 us, 32 -> 43.3 us per query.
BVH_LEAF_SIZE = 32


def scan_box_indices(boxes: np.ndarray, p, radius: float) -> np.ndarray:
    """Indices of boxes within `radius` of point p (conservative, O(N) scan)."""
    q = np.asarray(p, dtype=np.float64)
    return np.nonzero(np.all((q >= boxes[:, :3] - radius)
                             & (q <= boxes[:, 3:] + radius), axis=1))[0]


class HybridBoxIndex:
    """Point/radius proximity index over an (n, 6) box array.

    The index is immutable after construction and safe to share across
    queries.  Boxes are copied, made C-contiguous, and marked read-only
    on the way in (so later mutation of the caller's array, or of the
    stored boxes, cannot desynchronize the BVH node bounds built from
    them); only indices are returned.
    """

    def __init__(self, boxes: np.ndarray):
        b = np.array(boxes, dtype=np.float64, copy=True, order="C"
                     ).reshape(-1, 6)
        b.setflags(write=False)
        self._boxes = b
        # BVH node storage; all None when the vector scan serves queries.
        self._node_lo: np.ndarray | None = None
        self._node_hi: np.ndarray | None = None
        self._node_left: np.ndarray | None = None
        self._node_right: np.ndarray | None = None
        self._node_start: np.ndarray | None = None
        self._node_count: np.ndarray | None = None
        self._order: np.ndarray | None = None
        if len(b) >= BVH_THRESHOLD:
            self._build_bvh()

    def __len__(self) -> int:
        return len(self._boxes)

    @property
    def uses_bvh(self) -> bool:
        """True when queries traverse the BVH rather than the flat scan."""
        return self._node_lo is not None

    @property
    def boxes(self) -> np.ndarray:
        return self._boxes

    # -- construction ----------------------------------------------------
    def _build_bvh(self) -> None:
        boxes = self._boxes
        order = np.arange(len(boxes))
        los: list = []
        his: list = []
        lefts: list[int] = []
        rights: list[int] = []
        starts: list[int] = []
        counts: list[int] = []

        def rec(l: int, r: int) -> int:
            idx = order[l:r]
            blo = boxes[idx, :3].min(axis=0)
            bhi = boxes[idx, 3:].max(axis=0)
            node = len(los)
            los.append(blo)
            his.append(bhi)
            lefts.append(-1)
            rights.append(-1)
            starts.append(l)
            counts.append(r - l)
            if r - l > BVH_LEAF_SIZE:
                axis = int(np.argmax(bhi - blo))
                cen = (boxes[idx, :3] + boxes[idx, 3:]) * 0.5
                m = (l + r) // 2
                part = np.argpartition(cen[:, axis], m - l)
                order[l:r] = idx[part]
                lefts[node] = rec(l, m)
                rights[node] = rec(m, r)
            return node

        rec(0, len(boxes))
        self._node_lo = np.stack(los)
        self._node_hi = np.stack(his)
        self._node_left = np.asarray(lefts, dtype=np.int64)
        self._node_right = np.asarray(rights, dtype=np.int64)
        self._node_start = np.asarray(starts, dtype=np.int64)
        self._node_count = np.asarray(counts, dtype=np.int64)
        self._order = order

    # -- queries ---------------------------------------------------------
    def query_point(self, p, radius: float) -> np.ndarray:
        """Indices of boxes within `radius` of point p.

        The returned set is identical to the vector scan's set (sorted
        ascending, like ``np.nonzero``).
        """
        q = np.asarray(p, dtype=np.float64).reshape(3)
        if self._node_lo is None:
            _perf_count("vector_scan_query")
            return scan_box_indices(self._boxes, q, radius)
        _perf_count("bvh_query")
        nlo = self._node_lo
        nhi = self._node_hi
        nleft = self._node_left
        nright = self._node_right
        nstart = self._node_start
        ncount = self._node_count
        order = self._order
        boxes = self._boxes
        # Scalar comparisons in the hot loop: each node visit is a few
        # hundred ns this way versus ~2-3 us for vectorized np.any on
        # throwaway temporaries.  (q0, q1, q2, r are plain floats.)
        q0, q1, q2 = float(q[0]), float(q[1]), float(q[2])
        r = float(radius)
        hits: list[np.ndarray] = []
        stack = [0]
        while stack:
            nd = stack.pop()
            if (q0 < nlo[nd, 0] - r or q0 > nhi[nd, 0] + r
                    or q1 < nlo[nd, 1] - r or q1 > nhi[nd, 1] + r
                    or q2 < nlo[nd, 2] - r or q2 > nhi[nd, 2] + r):
                continue
            left = int(nleft[nd])
            if left < 0:
                s = int(nstart[nd])
                leaf = order[s:s + int(ncount[nd])]
                b = boxes[leaf]
                keep = np.all((q >= b[:, :3] - r)
                              & (q <= b[:, 3:] + r), axis=1)
                if bool(np.any(keep)):
                    hits.append(leaf[keep])
            else:
                stack.append(left)
                stack.append(int(nright[nd]))
        if not hits:
            return np.empty(0, dtype=np.intp)
        return np.sort(np.concatenate(hits)).astype(np.intp)
