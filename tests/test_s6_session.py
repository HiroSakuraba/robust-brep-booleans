#!/usr/bin/env python3
"""S6: BooleanSession tests.

A session prepares its base once and reuses the immutable indexes across
repeated operations; every operation still gets a fresh QueryContext
(S4) and the cold/warm verdicts must be identical (R1).
"""
import sys

sys.path.insert(0, "src")
sys.path.insert(0, "tests")

import pytest
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.gp import gp_Pnt

from brepkernel.assembly import _shape_volume
from brepkernel.pipeline import BRepAmbiguousResult, boolean_brep
from brepkernel.prepared import prepare_brep
from brepkernel.session import (
    BRepSessionClosed,
    BooleanSession,
    boolean_brep_many,
)


def _box(x, y, z, s=1.0):
    return BRepPrimAPI_MakeBox(gp_Pnt(x, y, z), s, s, s).Shape()


def _plate():
    sys.path.insert(0, "tools/review_probes")
    import bench_plate256 as bp
    return bp.build_drilled()


def _slot(dx=0.0):
    return BRepPrimAPI_MakeBox(
        gp_Pnt(3.0 + dx, 3.5, -0.25), 2.0, 1.0, 1.0).Shape()


def _volume(shape):
    return abs(float(_shape_volume(shape)))


def _ev_name(report):
    recs = report.get("evidence", {}).get("records", [])
    return recs[0].get("name") if recs else None


def test_cold_warm_verdict_identity():
    # R1: the session must not change accept/refuse, volume, or evidence.
    base, tool = _plate(), _slot()
    out_cold, rep_cold = boolean_brep(base, tool, "difference")
    with BooleanSession(base) as session:
        out_warm, rep_warm = session.boolean(tool, "difference")
    assert rep_cold["accepted"] == rep_warm["accepted"] is True
    assert _volume(out_cold) == _volume(out_warm)
    assert _ev_name(rep_cold) == _ev_name(rep_warm)


def test_base_prepared_once_across_ten_ops(monkeypatch):
    # Plan pass criterion: ten sequential ops build the base's immutable
    # indexes once.
    import brepkernel.prepared as prep_mod

    base = _plate()
    calls = {"base": 0}
    real_ensure = prep_mod.ensure_prepared

    def counting(shape_or_model, *, base_tol=1e-7):
        if shape_or_model is base:
            calls["base"] += 1
        return real_ensure(shape_or_model, base_tol=base_tol)

    monkeypatch.setattr(prep_mod, "ensure_prepared", counting)
    with BooleanSession(base) as session:
        for i in range(10):
            out, rep = session.boolean(_slot(dx=0.05 * i), "difference")
            assert rep["accepted"]
        assert session.op_count == 10
    assert calls["base"] == 1, f"base prepared {calls['base']}x, want 1"


def test_warm_ingest_skips_base_rebuild():
    # The warm op must not redo the ~300 ms base preparation.
    # (timings_ms is in milliseconds.)
    base, tool = _plate(), _slot()
    _, rep_cold = boolean_brep(base, tool, "difference")
    with BooleanSession(base) as session:
        _, rep_warm = session.boolean(tool, "difference")
    cold_ingest = rep_cold["timings_ms"]["ingest"]
    warm_ingest = rep_warm["timings_ms"]["ingest"]
    assert cold_ingest > 100, f"cold ingest suspiciously fast: {cold_ingest}"
    assert warm_ingest < 50, f"warm ingest rebuilt the base: {warm_ingest}"


def test_each_op_gets_fresh_query_context(monkeypatch):
    # S4 invariant: one Boolean call, one QueryContext. Count
    # instantiations with the shadow/fast paths disabled so each op
    # creates exactly one.
    import brepkernel.query as query_mod

    real_qc = query_mod.QueryContext
    count = {"n": 0}

    class CountingQC(real_qc):
        def __init__(self, *a, **k):
            count["n"] += 1
            super().__init__(*a, **k)

    monkeypatch.setattr(query_mod, "QueryContext", CountingQC)
    base = _box(0, 0, 0)
    with BooleanSession(base) as session:
        kw = {"fast_paths": False, "fast_path_shadow": False}
        session.boolean(_box(0.5, 0.5, 0.5), "union", **kw)
        session.boolean(_box(0.25, 0.25, 0.25), "difference", **kw)
    assert count["n"] == 2, f"QueryContext created {count['n']}x for 2 ops"


def test_context_manager_and_close():
    base = _box(0, 0, 0)
    session = BooleanSession(base)
    assert not session.closed
    out, rep = session.boolean(_box(0.5, 0.5, 0.5), "union")
    assert rep["accepted"] and session.op_count == 1
    session.close()
    assert session.closed
    with pytest.raises(BRepSessionClosed):
        session.boolean(_box(0.5, 0.5, 0.5), "union")

    with BooleanSession(base) as s2:
        assert not s2.closed
    assert s2.closed
    with pytest.raises(BRepSessionClosed):
        s2.boolean(_box(0.5, 0.5, 0.5), "union")


def test_base_tol_mismatch_rejected():
    with BooleanSession(_box(0, 0, 0)) as session:
        with pytest.raises(ValueError):
            session.boolean(_box(0.5, 0.5, 0.5), "union", base_tol=1e-6)


def test_session_survives_refusal(monkeypatch):
    # A refused op must not poison the session.
    import brepkernel.pipeline as pl

    real_impl = pl._boolean_brep_impl
    calls = {"n": 0}

    def fail_once(shapeA, shapeB, op, **kw):
        calls["n"] += 1
        if calls["n"] == 1:
            raise BRepAmbiguousResult(
                "simulated refusal", {"stages": {}, "accepted": False}, None)
        return real_impl(shapeA, shapeB, op, **kw)

    monkeypatch.setattr(pl, "_boolean_brep_impl", fail_once)
    with BooleanSession(_box(0, 0, 0)) as session:
        with pytest.raises(BRepAmbiguousResult):
            session.boolean(_box(0.5, 0.5, 0.5), "union")
        assert session.op_count == 0
        out, rep = session.boolean(_box(0.5, 0.5, 0.5), "union")
        assert rep["accepted"]
        assert session.op_count == 1


def test_prepared_tool_accepted():
    base = _box(0, 0, 0)
    with BooleanSession(base) as session:
        out_raw, _ = session.boolean(_box(0.5, 0.5, 0.5), "union")
        out_pre, rep_pre = session.boolean(
            prepare_brep(_box(0.5, 0.5, 0.5)), "union")
    assert rep_pre["accepted"]
    assert _volume(out_raw) == _volume(out_pre)


def test_session_pipeline_defaults_and_override():
    # Defaults apply to every op; per-call kwargs override them.
    base, far = _box(0, 0, 0), _box(5, 0, 0)
    with BooleanSession(base, collect_perf=True) as session:
        _, rep = session.boolean(far, "union")
        assert "performance" in rep
    with BooleanSession(base, fast_paths=False) as session:
        _, rep_off = session.boolean(far, "union")
        assert "s3_fast_path" not in rep_off["stages"]
        _, rep_on = session.boolean(far, "union", fast_paths=True)
        assert rep_on["stages"]["s3_fast_path"]["path"] == "A_disjoint"


def test_boolean_brep_many_matches_sequential():
    base = _box(0, 0, 0)
    tools = [_box(0.5, 0.5, 0.5), _box(-0.5, 0.2, 0.2), _box(0.2, -0.5, 0.3)]
    out, combined = boolean_brep_many(base, tools, "union")
    assert combined["optimized"] is False
    assert len(combined["steps"]) == 3

    ref = base
    for tool, step_rep in zip(tools, combined["steps"]):
        ref, ref_rep = boolean_brep(ref, tool, "union")
        assert step_rep["accepted"] == ref_rep["accepted"] is True
    assert _volume(out) == _volume(ref)
