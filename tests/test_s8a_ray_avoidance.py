"""S8A work-avoidance: signed line cast, per-solid AABB pruning, classifier reuse.

Behavioral contract: every S8A optimization is provably equivalent to the
old path --
  * _cast_line_pair == _cast_bidirectional on every (point, direction)
    (the full suite additionally runs in _LINE_CAST_SHADOW mode, which
    raises LineCastShadowMismatch on any divergence);
  * AABB pruning only skips solids the segment/point provably cannot
    touch (conservative boxes);
  * reused BRepClass3d_SolidClassifier instances are reset by Perform.
"""
import numpy as np
import pytest

from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.gp import gp_Pnt

from brepkernel import assembly as A
from brepkernel.assembly import (
    _MultiRayClassifier,
    _conservative_boxes,
    _occt_point_verdict,
    _segment_hits_aabb,
)
from brepkernel import perf as _perf_mod


def _box(x0, y0, z0, dx, dy, dz):
    return BRepPrimAPI_MakeBox(
        gp_Pnt(x0, y0, z0), gp_Pnt(x0 + dx, y0 + dy, z0 + dz)).Shape()


@pytest.fixture()
def two_solids():
    # Two disjoint boxes 10 units apart: the far solid is the pruning target.
    return [_box(0, 0, 0, 2, 2, 2), _box(12, 0, 0, 2, 2, 2)]


def test_line_cast_matches_two_ray(two_solids):
    """_cast_line_pair agrees with _cast_bidirectional everywhere probed,
    including near-tangent and near-edge lines."""
    tol = 1e-7
    ray = _MultiRayClassifier(two_solids, tol)
    pts = [
        np.array([1.0, 1.0, 1.0]),      # inside solid 0
        np.array([13.0, 1.0, 1.0]),    # inside solid 1
        np.array([7.0, 1.0, 1.0]),     # between them
        np.array([1.0, 1.0, 5.0]),     # above solid 0
        np.array([2.0 + 1e-4, 1.0, 1.0]),  # just outside solid 0's +x face
        np.array([1.0, 2.0 - 1e-4, 1.0]),  # just inside near an edge
        np.array([-5.0, -5.0, -5.0]),  # far away
    ]
    dirs = [np.array(d, float) for d in A._RAY_DIRECTIONS]
    for p in pts:
        for d in dirs:
            old = ray._cast_bidirectional(p, d)
            new = ray._cast_line_pair(p, d)
            assert old == new, (p, d, old, new)


def test_line_cast_shadow_flag_raises_nothing(two_solids):
    """classify() with the shadow flag on runs clean on tricky points."""
    old = A._LINE_CAST_SHADOW
    A._LINE_CAST_SHADOW = True
    try:
        ray = _MultiRayClassifier(two_solids, 1e-7)
        for p in (np.array([1.0, 1.0, 1.0]),
                  np.array([7.0, 1.0, 1.0]),
                  np.array([13.0, 1.0, 1.0])):
            assert ray.classify(p) in ("inside", "outside", "unknown")
    finally:
        A._LINE_CAST_SHADOW = old


def test_ray_aabb_prune_fires_and_stays_correct(two_solids):
    """Far solids are skipped (counter fires) and verdicts stay right."""
    tol = 1e-7
    ray = _MultiRayClassifier(two_solids, tol)
    with _perf_mod.scoped(_perf_mod.PerfCounters()) as pc:
        v = ray.classify(np.array([1.0, 1.0, 1.0]))
        assert v == "inside"
        v = ray.classify(np.array([7.0, 5.0, 1.0]))
        assert v == "outside"
    snap = pc.snapshot()
    assert snap.get("ray_aabb_prune", 0) > 0, snap


def test_segment_hits_aabb():
    box = np.array([0.0, 0.0, 0.0, 2.0, 2.0, 2.0])
    assert _segment_hits_aabb(np.array([-5.0, 1.0, 1.0]),
                              np.array([1.0, 0.0, 0.0]), box, 1e100)
    assert not _segment_hits_aabb(np.array([-5.0, 1.0, 1.0]),
                                  np.array([0.0, 1.0, 0.0]), box, 1e100)
    assert not _segment_hits_aabb(np.array([-5.0, 1.0, 1.0]),
                                  np.array([1.0, 0.0, 0.0]), box, 1.0)
    # Point inside the box always "hits".
    assert _segment_hits_aabb(np.array([1.0, 1.0, 1.0]),
                              np.array([0.3, -0.2, 0.9]), box, 1e100)


def test_occt_box_prune_matches_unpruned(two_solids):
    """Per-solid box pruning never changes the OCCT verdict."""
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    tol = 1e-7
    boxes = _conservative_boxes(two_solids)
    classifiers = [BRepClass3d_SolidClassifier(s) for s in two_solids]
    pts = [
        np.array([1.0, 1.0, 1.0]),
        np.array([13.0, 1.0, 1.0]),
        np.array([7.0, 1.0, 1.0]),
        np.array([2.0 + 5e-8, 1.0, 1.0]),  # within tol of solid 0 -> boundary
        np.array([-3.0, -3.0, -3.0]),
        np.array([12.0, 1.0, 3.0]),
    ]
    for p in pts:
        plain = _occt_point_verdict(p, two_solids, tol)
        pruned = _occt_point_verdict(p, two_solids, tol,
                                     _solid_boxes=boxes)
        reused = _occt_point_verdict(p, two_solids, tol,
                                     _solid_boxes=boxes,
                                     _classifiers=classifiers)
        assert pruned == plain, (p, plain, pruned)
        assert reused == plain, (p, plain, reused)


def test_occt_prune_counter_fires(two_solids):
    boxes = _conservative_boxes(two_solids)
    with _perf_mod.scoped(_perf_mod.PerfCounters()) as pc:
        # Far from both solids: both pruned.
        v = _occt_point_verdict(np.array([-50.0, 0.0, 0.0]), two_solids,
                                1e-7, _solid_boxes=boxes)
    assert v == "outside"
    assert pc.snapshot().get("occt_box_prune", 0) >= 2


def test_classifier_reuse_across_tolerances(two_solids):
    """One reused classifier serves several Perform tolerances faithfully."""
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    classifiers = [BRepClass3d_SolidClassifier(s) for s in two_solids]
    p = np.array([1.0, 1.0, 1.0])
    for tol in (1e-7, 1e-5, 1e-9):
        v = _occt_point_verdict(p, two_solids, tol,
                                _classifiers=classifiers)
        assert v == "inside"
    p2 = np.array([7.0, 1.0, 1.0])
    for tol in (1e-7, 1e-5, 1e-9):
        v = _occt_point_verdict(p2, two_solids, tol,
                                _classifiers=classifiers)
        assert v == "outside"


def test_s8a_end_to_end_volumes():
    """Boolean outcomes are unchanged by the S8A stack (analytic volumes)."""
    from brepkernel.pipeline import boolean_brep
    a = _box(0, 0, 0, 4, 4, 4)
    b = _box(2, 2, 2, 4, 4, 4)

    def vol(shape):
        from OCP.GProp import GProp_GProps
        from OCP.BRepGProp import BRepGProp
        pr = GProp_GProps()
        BRepGProp.VolumeProperties_s(shape, pr)
        return pr.Mass()

    _, rep = boolean_brep(a, b, "union")
    assert rep["accepted"]
    assert vol(_) == pytest.approx(64 + 64 - 8, rel=1e-6)
    _, rep = boolean_brep(a, b, "intersection")
    assert rep["accepted"]
    assert vol(_) == pytest.approx(8, rel=1e-6)
    _, rep = boolean_brep(a, b, "difference")
    assert rep["accepted"]
    assert vol(_) == pytest.approx(64 - 8, rel=1e-6)
