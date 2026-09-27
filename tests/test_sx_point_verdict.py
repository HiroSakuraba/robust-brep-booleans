#!/usr/bin/env python3
"""SX: point-verdict bbox short-circuit tests.

When a witness point is farther than `tol` outside the other model's
conservative bbox, the verdict is provably "outside" (the solid lies
inside its bbox; OCCT's ON band is exactly `tol`) and both classifiers
are skipped.  The short-circuit is conservative: it fires only when the
geometric proof holds; every other point takes the full dual path.
"""
import sys

sys.path.insert(0, "src")
sys.path.insert(0, "tests")

import numpy as np
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.gp import gp_Pnt

from brepkernel import perf as _perf_mod
from brepkernel.assembly import (
    _agreed_point_verdict,
    _MultiRayClassifier,
    _classify_point_in_model,
    AssemblyError,
)
from brepkernel.prepared import ensure_prepared
from brepkernel.query import QueryContext


def _box(x, y, z, dx, dy, dz):
    return BRepPrimAPI_MakeBox(gp_Pnt(x, y, z), dx, dy, dz).Shape()


def _model_and_ray(x, y, z, dx, dy, dz, tol=1e-7):
    model = ensure_prepared(_box(x, y, z, dx, dy, dz), base_tol=tol).model
    ray = _MultiRayClassifier([sr.solid for sr in model.solids], tol)
    ctx = QueryContext(base_tol=tol)
    return model, ray, ctx


def _dual_reference(point, model, tol, ray):
    """The verdict the full dual path would return (no short-circuit)."""
    occt = _classify_point_in_model(np.asarray(point, float), model, tol)
    if occt not in ("inside", "outside"):
        return occt
    independent = ray.classify(np.asarray(point, float))
    if independent == occt:
        return occt
    return "DISAGREE"


def test_bbox_shortcircuit_triggers_and_skips_classifiers():
    """Far-outside point -> "outside" with neither classifier evaluated."""
    model, ray, ctx = _model_and_ray(0, 0, 0, 1, 1, 1)
    pt = np.array([50.0, 50.0, 50.0])
    with _perf_mod.scoped(_perf_mod.PerfCounters()) as pc:
        v = _agreed_point_verdict(pt, model, 1e-7, ray, ctx)
    assert v == "outside", f"verdict: {v}"
    snap = pc.snapshot()
    assert snap.get("point_verdict_bbox_shortcircuit", 0) == 1, snap
    assert snap.get("solid_classifier_eval", 0) == 0, snap
    assert snap.get("ray_intersector_perform", 0) == 0, snap
    print("[PASS] test_bbox_shortcircuit_triggers_and_skips_classifiers")


def test_shortcircuit_matches_dual_path():
    """On a battery of points, short-circuit == full dual verdict."""
    model, ray, ctx = _model_and_ray(0, 0, 0, 4, 4, 4)
    pts = [
        [100.0, 0.0, 0.0], [-100.0, -100.0, -100.0], [2.0, 2.0, 100.0],
        [2.0, 2.0, 2.0], [0.5, 0.5, 0.5], [3.999, 2.0, 2.0],
        [4.0 + 5e-8, 2.0, 2.0],   # within tol of the bbox: dual path
        [10.0, 10.0, 10.0], [4.5, 4.5, 4.5],
    ]
    for p in pts:
        ctx2 = QueryContext(base_tol=1e-7)
        v = _agreed_point_verdict(np.array(p), model, 1e-7, ray, ctx2)
        ref = _dual_reference(p, model, 1e-7, ray)
        assert v == ref, f"point {p}: short-circuit path gave {v}, dual gave {ref}"
    print("[PASS] test_shortcircuit_matches_dual_path")


def test_within_tol_of_bbox_stays_on_dual_path():
    """Point within tol of the bbox must NOT take the short-circuit."""
    model, ray, ctx = _model_and_ray(0, 0, 0, 1, 1, 1)
    pt = np.array([1.0 + 5e-8, 0.5, 0.5])  # 5e-8 outside the face, tol=1e-7
    with _perf_mod.scoped(_perf_mod.PerfCounters()) as pc:
        v = _agreed_point_verdict(pt, model, 1e-7, ray, ctx)
    snap = pc.snapshot()
    assert snap.get("point_verdict_bbox_shortcircuit", 0) == 0, snap
    assert v in ("inside", "outside", "boundary", "unknown"), v
    print(f"[PASS] test_within_tol_of_bbox_stays_on_dual_path (verdict={v})")


def test_inside_point_stays_on_dual_path():
    """A genuinely inside point takes the full dual path."""
    model, ray, ctx = _model_and_ray(0, 0, 0, 4, 4, 4)
    pt = np.array([2.0, 2.0, 2.0])
    with _perf_mod.scoped(_perf_mod.PerfCounters()) as pc:
        v = _agreed_point_verdict(pt, model, 1e-7, ray, ctx)
    snap = pc.snapshot()
    assert snap.get("point_verdict_bbox_shortcircuit", 0) == 0, snap
    assert v == "inside", f"verdict: {v}"
    print("[PASS] test_inside_point_stays_on_dual_path")


def test_shell_only_model_still_raises():
    """A model with no solids must raise, not silently say outside."""
    model, ray, ctx = _model_and_ray(0, 0, 0, 1, 1, 1)
    model.solids = []
    pt = np.array([50.0, 50.0, 50.0])
    try:
        _agreed_point_verdict(pt, model, 1e-7, ray, ctx)
    except AssemblyError as e:
        assert "closed OCCT solids" in str(e), str(e)
        print("[PASS] test_shell_only_model_still_raises")
        return
    raise AssertionError("shell-only model did not raise")


def test_no_ctx_uncached_path():
    """ctx=None (S3 fast-path probe) still short-circuits correctly."""
    model, ray, _ = _model_and_ray(0, 0, 0, 1, 1, 1)
    pt = np.array([50.0, 50.0, 50.0])
    with _perf_mod.scoped(_perf_mod.PerfCounters()) as pc:
        v = _agreed_point_verdict(pt, model, 1e-7, ray, None)
    assert v == "outside", f"verdict: {v}"
    assert pc.snapshot().get("point_verdict_bbox_shortcircuit", 0) == 1
    print("[PASS] test_no_ctx_uncached_path")


def test_memo_stores_shortcircuit_verdict():
    """A repeated far point hits the ctx memo (no recompute at all)."""
    model, ray, ctx = _model_and_ray(0, 0, 0, 1, 1, 1)
    pt = np.array([50.0, 50.0, 50.0])
    with _perf_mod.scoped(_perf_mod.PerfCounters()) as pc:
        v1 = _agreed_point_verdict(pt, model, 1e-7, ray, ctx)
        v2 = _agreed_point_verdict(pt, model, 1e-7, ray, ctx)
    assert (v1, v2) == ("outside", "outside")
    snap = pc.snapshot()
    assert snap.get("point_verdict_bbox_shortcircuit", 0) == 1, snap
    assert snap.get("ctx_point_verdict_hit", 0) == 1, snap
    print("[PASS] test_memo_stores_shortcircuit_verdict")


def main():
    test_bbox_shortcircuit_triggers_and_skips_classifiers()
    test_shortcircuit_matches_dual_path()
    test_within_tol_of_bbox_stays_on_dual_path()
    test_inside_point_stays_on_dual_path()
    test_shell_only_model_still_raises()
    test_no_ctx_uncached_path()
    test_memo_stores_shortcircuit_verdict()
    print("ALL PASS")


if __name__ == "__main__":
    main()
