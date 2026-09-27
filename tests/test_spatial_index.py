"""S2 HybridBoxIndex tests.

Pass criteria from the speed plan:
- for one million randomized point/radius queries over synthetic boxes, the
  BVH candidate set contains every candidate the vector scan returns
  (missing candidates are not allowed);
- exact downstream distance verdicts match the existing implementation on
  the full suite (the suite re-run; t4 checks verdict equivalence here);
- at 4096 boxes the median proximity query is at least 3x faster than the
  full vector scan on the same runner (measured by
  tools/review_probes/bench_s2_bvh.py; t2 pins the calibrated threshold).

The full 1M-query gate runs as tools/review_probes/probe_s2_1M.py; t1 below
is the committed, faster version of the same check.
"""
import sys

sys.path.insert(0, "src")
sys.path.insert(0, "tests")

import numpy as np

from brepkernel import boolean_brep, BRepAmbiguousResult
from brepkernel import perf as perf_mod
from brepkernel.prepared import prepare_brep
from brepkernel.spatial import (
    BVH_THRESHOLD,
    HybridBoxIndex,
    scan_box_indices,
)

from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
from OCP.GProp import GProp_GProps
from OCP.TopAbs import TopAbs_SOLID
from OCP.TopExp import TopExp_Explorer
from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def volume(shape):
    p = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, p, 1e-9, True)
    return float(p.Mass())


def solid_count(shape):
    n = 0
    ex = TopExp_Explorer(shape, TopAbs_SOLID)
    while ex.More():
        n += 1
        ex.Next()
    return n


def boxes():
    a = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 1, 1, 1).Shape()
    b = BRepPrimAPI_MakeBox(gp_Pnt(0.5, 0.5, 0.5), 1, 1, 1).Shape()
    return a, b


def cyl_box():
    box = BRepPrimAPI_MakeBox(gp_Pnt(-1, -1, -1), 2, 2, 2).Shape()
    cyl = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(0, 0, -2), gp_Dir(0, 0, 1)), 0.3, 4).Shape()
    return box, cyl


def run(a, b, op, **kw):
    """Returns (accepted, volume_or_None, solids, refusal_kind, report)."""
    try:
        out, report = boolean_brep(a, b, op, **kw)
    except BRepAmbiguousResult as exc:
        rep = exc.report
        kind = rep.get("refusal", {}).get("kind", "?")
        return False, None, 0, kind, rep
    return True, volume(out), solid_count(out), None, report


def random_boxes(n, rng, degenerate=False):
    cen = rng.uniform(0.0, 10.0, size=(n, 3))
    half = rng.uniform(0.0, 0.5, size=(n, 3))
    boxes = np.concatenate([cen - half, cen + half], axis=1)
    if degenerate:
        boxes[: n // 10, 3:] = boxes[: n // 10, :3]  # zero-volume boxes
        boxes[n // 10: n // 5] = boxes[0]  # duplicated boxes
    return boxes


def t1_randomized_correctness():
    """BVH query set == vector scan set (exact; missing is forbidden)."""
    from brepkernel import spatial as spatial_mod
    rng = np.random.default_rng(20260926)
    ok = True
    n_checked = 0
    old_threshold = spatial_mod.BVH_THRESHOLD
    try:
        for n in (64, 255, 256, 1024, 4096):
            boxes = random_boxes(n, rng, degenerate=True)
            spatial_mod.BVH_THRESHOLD = 0  # force the BVH path
            idx = HybridBoxIndex(boxes)
            ok &= check(f"t1 N={n}: BVH active", idx.uses_bvh)
            corners = boxes[rng.integers(0, n, size=(500,)), :3]
            for i in range(4000):
                if i % 10 == 0:
                    p = corners[i % 500]  # exactly on box boundaries
                else:
                    p = rng.uniform(-1.0, 11.0, size=(3,))
                r = float(10.0 ** rng.uniform(-9.0, -1.0))
                if i % 50 == 0:
                    r = 0.0
                a = scan_box_indices(boxes, p, r)
                b = idx.query_point(p, r)
                same = np.array_equal(a, b)
                n_checked += 1
                if not same:
                    sa, sb = set(a.tolist()), set(b.tolist())
                    ok &= check(f"t1 N={n} query {i}: sets equal", False,
                                f"missing={sorted(sa - sb)[:5]} "
                                f"extra={sorted(sb - sa)[:5]}")
                    break
            else:
                ok &= check(f"t1 N={n}: 4000 queries identical", True)
    finally:
        spatial_mod.BVH_THRESHOLD = old_threshold
    ok &= check("t1 total queries checked", n_checked == 5 * 4000,
                f"={n_checked}")
    return ok


def t2_threshold_behavior():
    """Below the calibrated threshold the flat scan serves; above, the BVH."""
    rng = np.random.default_rng(7)
    ok = True
    small = HybridBoxIndex(random_boxes(BVH_THRESHOLD - 1, rng))
    big = HybridBoxIndex(random_boxes(BVH_THRESHOLD, rng))
    ok &= check("t2 below threshold: vector path",
                not small.uses_bvh, f"N={BVH_THRESHOLD - 1}")
    ok &= check("t2 at threshold: BVH path",
                big.uses_bvh, f"N={BVH_THRESHOLD}")
    # The path taken is observable through the perf counters.
    counters = perf_mod.PerfCounters()
    with perf_mod.scoped(counters):
        p = rng.uniform(0.0, 10.0, size=(3,))
        small.query_point(p, 1e-3)
        big.query_point(p, 1e-3)
    c = counters.counts
    ok &= check("t2 small query counted as vector scan",
                c.get("vector_scan_query", 0) == 1, f"={c}")
    ok &= check("t2 large query counted as BVH",
                c.get("bvh_query", 0) == 1, f"={c}")
    return ok


def t3_degenerate_inputs():
    """Empty, singleton, huge-offset, and on-boundary queries."""
    ok = True
    empty = HybridBoxIndex(np.empty((0, 6)))
    ok &= check("t3 empty index: no BVH, empty result",
                not empty.uses_bvh
                and empty.query_point([0, 0, 0], 1.0).shape == (0,))
    one = HybridBoxIndex(np.array([[0.0, 0.0, 0.0, 1.0, 1.0, 1.0]]))
    ok &= check("t3 single box hit",
                one.query_point([0.5, 0.5, 0.5], 0.0).tolist() == [0])
    ok &= check("t3 single box miss",
                one.query_point([5.0, 5.0, 5.0], 0.1).shape == (0,))
    # Huge offsets: the pruning proof is fp-exact, check it at 1e8 scale.
    rng = np.random.default_rng(11)
    big = random_boxes(300, rng) + 1e8
    from brepkernel import spatial as spatial_mod
    old = spatial_mod.BVH_THRESHOLD
    try:
        spatial_mod.BVH_THRESHOLD = 0
        idx = HybridBoxIndex(big)
        for _ in range(2000):
            p = rng.uniform(1e8 - 1.0, 1e8 + 11.0, size=(3,))
            r = float(10.0 ** rng.uniform(-7.0, 1.0))
            ok &= np.array_equal(scan_box_indices(big, p, r),
                                 idx.query_point(p, r))
            if not ok:
                break
    finally:
        spatial_mod.BVH_THRESHOLD = old
    ok &= check("t3 1e8-offset: 2000 queries identical", ok)
    # Infinite radius returns everything, like the scan.
    boxes = random_boxes(300, rng)
    idx = HybridBoxIndex(boxes)
    got = idx.query_point([5.0, 5.0, 5.0], float("inf"))
    ok &= check("t3 inf radius returns all",
                np.array_equal(got, np.arange(300)))
    return ok


def t4_integration_verdict_equivalence():
    """Prepared booleans with S2 live: verdicts identical, zero rebuilds."""
    ok = True
    cases = [("boxes", boxes()), ("cyl_box", cyl_box())]
    for cname, (sa, sb) in cases:
        for op in ("union", "intersection", "difference"):
            r_raw = run(sa, sb, op)
            pa, pb = prepare_brep(sa), prepare_brep(sb)
            try:
                out, rep = boolean_brep(pa, pb, op, collect_perf=True)
                ok_i, vol_i = True, volume(out)
            except BRepAmbiguousResult as exc:
                ok_i, vol_i, rep = False, None, exc.report
            same = (r_raw[0] == ok_i
                    and (not r_raw[0] or r_raw[1] == vol_i))
            ok &= check(f"t4 {cname}/{op}: S2 verdict == raw", same,
                        f"raw={r_raw[:3]} pre={(ok_i, vol_i)}")
            c = rep["performance"]["counters"]
            ok &= check(f"t4 {cname}/{op}: spatial index queried",
                        c.get("bvh_query", 0) + c.get("vector_scan_query", 0) >= 1,
                        f"={c.get('bvh_query', 0)}/{c.get('vector_scan_query', 0)}")
            ok &= check(f"t4 {cname}/{op}: zero face-box rebuilds",
                        c.get("face_box_build", -1) == 0)
            # The prepared deduped-edge path is taken (S1); the only
            # edge_box_builds left are the result-side single-solid
            # classifiers, which S1 deliberately leaves unprepared
            # (result geometry does not exist at prepare time).
            ok &= check(f"t4 {cname}/{op}: prepared edge path taken",
                        c.get("prepared_edge_index_hit", 0) >= 1,
                        f"={c.get('prepared_edge_index_hit', 0)}")
    return ok


def t5_mismatched_index_falls_back():
    """A stale index never changes a distance; the old paths take over."""
    from brepkernel import assembly as asm
    from brepkernel.step_ingest import index_shape
    ok = True
    sa, _ = boxes()
    model = index_shape(sa)
    pts = np.array([[0.25, 0.25, 0.25], [5.0, 5.0, 5.0]])
    ref = asm._point_boundary_distances(pts, model, cap=0.5)
    # Wrong-length index: falls back to the per-call build, same distances.
    bad = HybridBoxIndex(np.zeros((3, 6)))
    got = asm._point_boundary_distances(pts, model, cap=0.5, face_index=bad)
    ok &= check("t5 bad face_index falls back, distances unchanged",
                np.allclose(ref, got), f"ref={ref} got={got}")
    # Matching index over the prepared boxes: identical distances.
    pa = prepare_brep(sa)
    got2 = asm._point_boundary_distances(
        pts, model, cap=0.5, face_index=pa.face_index)
    ok &= check("t5 prepared face_index: distances identical",
                np.allclose(ref, got2), f"ref={ref} got={got2}")
    # The prepared boxes path (S1) still works alongside.
    got3 = asm._point_boundary_distances(
        pts, model, cap=0.5, face_boxes=pa.face_boxes)
    ok &= check("t5 S1 face_boxes path unchanged",
                np.allclose(ref, got3), f"ref={ref} got={got3}")
    return ok


TESTS = [t1_randomized_correctness,
         t2_threshold_behavior,
         t3_degenerate_inputs,
         t4_integration_verdict_equivalence,
         t5_mismatched_index_falls_back]


def main():
    results = []
    for t in TESTS:
        try:
            results.append((t.__name__, bool(t())))
        except Exception as exc:  # noqa: BLE001 - report, don't abort
            print(f"[FAIL] {t.__name__} raised {type(exc).__name__}: {exc}")
            results.append((t.__name__, False))
    failed = [n for n, r in results if not r]
    print(f"\n{len(results) - len(failed)}/{len(results)} test groups passed")
    if failed:
        print("FAILED:", ", ".join(failed))
        raise SystemExit(1)


if __name__ == "__main__":
    main()
