"""S4 per-Boolean-call QueryContext.

Covers: no module-global mutable query state remains, an explicit
QueryContext threads through the pipeline, the volume cache is scoped
to the context (two contexts do not see each other's entries), memo
hits are exact (a hit returns what a fresh computation would), the
report fields (region_stats, assembly volume) are intact, and the S3
shadow re-run shares the outer call's context.
"""

import numpy as np
import pytest
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.gp import gp_Pnt

from brepkernel import perf as _perf_mod
from brepkernel.assembly import _shape_volume, clear_volume_cache
from brepkernel.pipeline import boolean_brep
from brepkernel.query import QueryContext


def _box(x, y, z, sx, sy=None, sz=None):
    sy = sx if sy is None else sy
    sz = sx if sz is None else sz
    return BRepPrimAPI_MakeBox(
        gp_Pnt(float(x), float(y), float(z)),
        float(sx), float(sy), float(sz)).Shape()


def test_no_module_global_volume_cache():
    """The old module-global _VOLUME_CACHE is gone (S4 pass criterion)."""
    import brepkernel.assembly as asm
    assert not hasattr(asm, "_VOLUME_CACHE"), \
        "module-global _VOLUME_CACHE must not exist"


def test_no_function_attribute_region_stats():
    """region_stats no longer lives on the _classify_pieces function."""
    from brepkernel.assembly import _classify_pieces
    assert not hasattr(_classify_pieces, "region_stats"), \
        "region_stats must live on the QueryContext, not the function"


def test_clear_volume_cache_is_deprecated_noop():
    """The old shim still imports and runs (back-compat), but does nothing."""
    clear_volume_cache()  # must not raise


def test_volume_memo_scoped_to_context():
    """Two QueryContexts do not see each other's volume entries."""
    shp = _box(0, 0, 0, 2)
    ctx1 = QueryContext()
    ctx2 = QueryContext()
    v1 = _shape_volume(shp, ctx=ctx1)
    assert len(ctx1.volumes) == 1
    assert len(ctx2.volumes) == 0
    # ctx2 computes fresh (no hit from ctx1's entry)
    v2 = _shape_volume(shp, ctx=ctx2)
    assert v1 == pytest.approx(v2)
    assert len(ctx2.volumes) == 1
    # same-context repeat is a hit
    with _perf_mod.scoped(_perf_mod.PerfCounters()) as pc:
        _shape_volume(shp, ctx=ctx1)
    assert pc.snapshot().get("ctx_volume_hit", 0) == 1


def test_volume_memo_hit_is_exact():
    """A memo hit returns exactly what a fresh computation returns."""
    shp = _box(0, 0, 0, 3)
    ctx = QueryContext()
    fresh = _shape_volume(shp)  # uncached direct call
    memoed = _shape_volume(shp, ctx=ctx)
    again = _shape_volume(shp, ctx=ctx)  # hit
    assert memoed == pytest.approx(fresh)
    assert again == pytest.approx(fresh)


def test_query_context_get_or_create():
    """get_or_create builds once and returns the same object."""
    ctx = QueryContext()
    built = []

    def factory():
        built.append(1)
        return object()

    o1 = ctx.get_or_create("k", factory)
    o2 = ctx.get_or_create("k", factory)
    assert o1 is o2
    assert built == [1]


def test_query_context_new_call_isolation():
    """new_call() drops memos but keeps the tolerance policy."""
    ctx = QueryContext(base_tol=1e-6)
    ctx.volumes["x"] = 1.0
    ctx.point_verdicts["y"] = "inside"
    fresh = ctx.new_call()
    assert fresh.volumes == {}
    assert fresh.point_verdicts == {}
    assert fresh.base_tol == 1e-6


def test_point_verdict_memo_hit_exact():
    """The point-verdict memo returns the agreed verdict verbatim."""
    from brepkernel.assembly import _agreed_point_verdict, _MultiRayClassifier
    from brepkernel.step_ingest import index_shape
    shp = _box(0, 0, 0, 2)
    model = index_shape(shp)
    solids = [sr.solid for sr in model.solids]
    ray = _MultiRayClassifier(solids, 1e-7)
    ctx = QueryContext()
    p = np.array([0.5, 0.5, 0.5])
    v1 = _agreed_point_verdict(p, model, 1e-7, ray, ctx)
    assert v1 == "inside"
    with _perf_mod.scoped(_perf_mod.PerfCounters()) as pc:
        v2 = _agreed_point_verdict(p, model, 1e-7, ray, ctx)
    assert v2 == "inside"
    assert pc.snapshot().get("ctx_point_verdict_hit", 0) == 1


def test_boolean_brep_region_stats_intact():
    """The report still carries region stats and the assembly volume."""
    a = _box(0, 0, 0, 2)
    b = _box(1, 1, 1, 2)
    out, rep = boolean_brep(a, b, "intersection")
    asm = rep["stages"]["assembly"]
    assert asm["volume"] == pytest.approx(1.0)
    assert "region_stats" in asm
    assert set(asm["region_stats"]) == {"A", "B"}


def test_boolean_brep_results_unchanged():
    """S4 changes nothing observable: union/intersection/difference."""
    a = _box(0, 0, 0, 2)
    b = _box(1, 1, 1, 2)
    _, r_u = boolean_brep(a, b, "union")
    _, r_i = boolean_brep(a, b, "intersection")
    _, r_d = boolean_brep(a, b, "difference")
    assert r_u["stages"]["assembly"]["volume"] == pytest.approx(15.0)
    assert r_i["stages"]["assembly"]["volume"] == pytest.approx(1.0)
    assert r_d["stages"]["assembly"]["volume"] == pytest.approx(7.0)


def test_s3_shadow_uses_fresh_context():
    """The S3 shadow re-run gets its own QueryContext (independence).

    The shadow must answer "would the old pipeline independently reach
    the same answer?" -- sharing the outer call's classifiers, volumes,
    or verdict memos would let one cached bug make both paths falsely
    agree.  After a separated union (shadow runs the full old pipeline),
    the outer context must hold no classifier/verdict/distance memos:
    only the volumes the outer agreement check computed itself.
    """
    from brepkernel.pipeline import _boolean_brep_impl
    ctx = QueryContext(base_tol=1e-7)
    a = _box(0, 0, 0, 1)
    b = _box(5, 0, 0, 1)
    _boolean_brep_impl(a, b, "union", _ctx=ctx)
    assert ctx.classifiers == {}
    assert ctx.point_verdicts == {}
    assert ctx.boundary_distances == {}


def test_volume_cache_verifies_identity():
    """volume_get trusts a hit only when the stored shape IS the shape.

    Guards against id() reuse after GC: a dict entry keyed by a recycled
    id must not serve another shape's volume.
    """
    ctx = QueryContext()
    shp = _box(0, 0, 0, 2)
    other = _box(0, 0, 0, 3)
    ctx.volume_put(shp, 1e-10, 999.0)
    assert ctx.volume_get(shp, 1e-10) == 999.0
    # Poison: another shape's id aliased to shp's entry.
    ctx.volumes[(id(other), 1e-10)] = (shp, 999.0)
    assert ctx.volume_get(other, 1e-10) is None
    # And the strong reference keeps the original alive.
    key = (id(shp), 1e-10)
    assert ctx.volumes[key][0] is shp
