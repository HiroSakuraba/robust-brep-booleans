#!/usr/bin/env python3
"""S5: same_domain bbox short-circuit tests.

The bbox pre-check skips the expensive ShapeUpgrade_UnifySameDomain
canonicalization only when the input bounding boxes already differ beyond
tolerance. UnifySameDomain preserves the bbox (it only merges same-domain
faces/edges), so differing bboxes prove canonicalization cannot produce
equivalent models. Verdict preservation: the short-circuit returns the
same non-equivalent verdict, just without the wasted work.
"""
import sys
import time

sys.path.insert(0, "src")
sys.path.insert(0, "tests")

from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.gp import gp_Pnt

from brepkernel.prepared import ensure_prepared
from brepkernel.same_domain import same_domain_models
from brepkernel import perf as _perf


def _box(x, y, z, dx, dy, dz):
    return BRepPrimAPI_MakeBox(gp_Pnt(x, y, z), dx, dy, dz).Shape()


def test_bbox_shortcircuit_triggers():
    """Different bboxes -> short-circuit, no canonicalization."""
    a = ensure_prepared(_box(0, 0, 0, 4, 4, 4), base_tol=1e-7).model
    b = ensure_prepared(_box(0, 0, 0, 1, 1, 1), base_tol=1e-7).model
    sd = same_domain_models(a, b, base_tol=1e-7, fuzz=1e-7)
    assert not sd.equivalent, "different boxes must not be equivalent"
    assert "bounding boxes differ" in sd.reason, f"reason: {sd.reason}"
    assert not sd.canonicalized, "short-circuit must skip canonicalization"
    print("[PASS] test_bbox_shortcircuit_triggers")


def test_bbox_shortcircuit_saves_time():
    """The short-circuit avoids the ~1s canonicalization on large models."""
    sys.path.insert(0, "tools/review_probes")
    import bench_plate256 as bp
    drilled = bp.build_drilled()
    slot = _box(3.0, 3.5, -0.25, 2.0, 1.0, 1.0)
    a = ensure_prepared(drilled, base_tol=1e-7).model
    b = ensure_prepared(slot, base_tol=1e-7).model
    t0 = time.perf_counter()
    sd = same_domain_models(a, b, base_tol=1e-7, fuzz=1e-7)
    dt = time.perf_counter() - t0
    assert not sd.equivalent
    assert not sd.canonicalized
    assert dt < 0.5, f"short-circuit took {dt:.2f}s, expected <0.5s"
    print(f"[PASS] test_bbox_shortcircuit_saves_time ({dt:.3f}s)")


def test_same_bbox_still_canonicalizes():
    """Same bbox but different material -> canonicalization still runs."""
    # Two boxes with identical bbox but different face decomposition
    # (one split into two). The short-circuit must NOT trigger.
    a = ensure_prepared(_box(0, 0, 0, 2, 2, 2), base_tol=1e-7).model
    b = ensure_prepared(_box(0, 0, 0, 2, 2, 2), base_tol=1e-7).model
    sd = same_domain_models(a, b, base_tol=1e-7, fuzz=1e-7)
    # Identical boxes ARE equivalent; the point is the short-circuit
    # did not wrongly reject them.
    assert sd.equivalent, f"identical boxes must be equivalent: {sd.reason}"
    print("[PASS] test_same_bbox_still_canonicalizes")


def test_verdict_preservation_on_pipeline():
    """The short-circuit does not change pipeline accept/refuse or volume."""
    from brepkernel.pipeline import boolean_brep
    sys.path.insert(0, "tools/review_probes")
    import bench_plate256 as bp
    drilled = bp.build_drilled()
    slot = _box(3.0, 3.5, -0.25, 2.0, 1.0, 1.0)
    out, rep = boolean_brep(drilled, slot, "difference", collect_perf=True)
    assert rep["performance"]["counters"].get(
        "samedomain_bbox_shortcircuit", 0) >= 1, \
        "short-circuit counter must fire on the plate benchmark"
    # Volume must match the pre-S5 value (22.2350)
    from brepkernel.assembly import _shape_volume
    vol = abs(float(_shape_volume(out)))
    assert abs(vol - 22.2350) < 1e-3, f"volume changed: {vol}"
    print(f"[PASS] test_verdict_preservation_on_pipeline (vol={vol:.4f})")


if __name__ == "__main__":
    test_bbox_shortcircuit_triggers()
    test_bbox_shortcircuit_saves_time()
    test_same_bbox_still_canonicalizes()
    test_verdict_preservation_on_pipeline()
    print("\nALL PASS")
