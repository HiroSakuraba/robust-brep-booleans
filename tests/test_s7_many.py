"""S7 Phase 2 tests: optimized boolean_brep_many vs the sequential reference.

R1: zero accept/refuse changes, no real volume change (1e-9 band admits
only float summation-order noise), same face/solid counts.
"""
import math
import random
import sys

import pytest

from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID

from brepkernel.many import boolean_brep_many
from brepkernel.session import boolean_brep_many as reexported_many
from brepkernel.pipeline import BRepAmbiguousResult
from brepkernel.assembly import _shape_volume


def plate():
    return BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 8.0, 8.0, 0.5).Shape()


def count(shape, what):
    ex = TopExp_Explorer(shape, what)
    n = 0
    while ex.More():
        n += 1
        ex.Next()
    return n


def vol(shape):
    return abs(float(_shape_volume(shape)))


def run_many(*args, **kwargs):
    """Run boolean_brep_many, catching refusals as (None, exc)."""
    try:
        return boolean_brep_many(*args, **kwargs), None
    except BRepAmbiguousResult as e:
        return None, e


def must_accept(res, exc):
    """Unpack a run_many result, failing the test on refusal."""
    assert exc is None, f"unexpected refusal: {exc}"
    return res  # (out, summary)


def refusal_key(exc):
    # Stage + underlying cause type identify the refusal kind (R1: zero
    # refusal-kind changes).
    return (str(exc).split(":")[0], type(exc.cause).__name__)


def assert_same_result(ref_res, ref_exc, opt_res, opt_exc):
    if ref_exc is not None or opt_exc is not None:
        # A refusal on both sides with the same kind is equivalence.
        assert ref_exc is not None and opt_exc is not None, \
            f"refusal mismatch: ref={ref_exc} opt={opt_exc}"
        assert refusal_key(ref_exc) == refusal_key(opt_exc), \
            f"refusal-kind mismatch: {refusal_key(ref_exc)} vs {refusal_key(opt_exc)}"
        return
    (ref_out, ref_sum), (opt_out, opt_sum) = ref_res, opt_res
    ra = all(r["accepted"] for r in ref_sum["steps"])
    oa = all(r["accepted"] for r in opt_sum["steps"])
    assert ra == oa, "acceptance mismatch"
    assert opt_sum["optimized"] in (True, False)
    if ra:
        rv, ov = vol(ref_out), vol(opt_out)
        assert abs(rv - ov) <= 1e-9 * max(1.0, rv), \
            f"volume changed: {rv} vs {ov}"
        assert count(ref_out, TopAbs_FACE) == count(opt_out, TopAbs_FACE)
        assert count(ref_out, TopAbs_SOLID) == count(opt_out, TopAbs_SOLID)


def random_tool(rng):
    """Seeded random tool: box or cylinder, on/near/off the plate."""
    x = rng.uniform(-2.0, 10.0)
    y = rng.uniform(-2.0, 10.0)
    z = rng.uniform(-0.6, 0.4)
    s = rng.uniform(0.2, 0.6)
    if rng.random() < 0.5:
        return BRepPrimAPI_MakeBox(gp_Pnt(x, y, z), s, s,
                                   rng.uniform(0.4, 1.2)).Shape()
    return BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(x, y, z - 0.2), gp_Dir(0, 0, 1)),
        s / 2, rng.uniform(0.4, 1.4)).Shape()


def test_reexport_from_session():
    assert reexported_many is boolean_brep_many


def test_zero_tools_returns_base():
    b = plate()
    out, summary = boolean_brep_many(b, [], "difference", optimized=True)
    assert summary["optimized"] is True
    assert summary["steps"] == []
    assert vol(out) == vol(b)


def test_single_tool_matches_reference():
    b = plate()
    t = BRepPrimAPI_MakeBox(gp_Pnt(2, 2, -0.25), 0.4, 0.4, 1.0).Shape()
    r_res, r_exc = run_many(b, [t], "difference")
    o_res, o_exc = run_many(b, [t], "difference", optimized=True)
    assert_same_result(r_res, r_exc, o_res, o_exc)


def test_all_tools_skipped_difference():
    b = plate()
    far = [BRepPrimAPI_MakeBox(gp_Pnt(20 + i, 20, 20), 1, 1, 1).Shape()
           for i in range(3)]
    o_res, o_exc = run_many(b, far, "difference", optimized=True)
    o_out, o_sum = must_accept(o_res, o_exc)
    assert o_sum["n_skipped"] == 3
    assert o_sum["steps"] == []
    assert vol(o_out) == vol(b)
    assert count(o_out, TopAbs_FACE) == count(b, TopAbs_FACE)


def test_some_tools_skipped_difference():
    b = plate()
    hit = BRepPrimAPI_MakeBox(gp_Pnt(2, 2, -0.25), 0.4, 0.4, 1.0).Shape()
    far = BRepPrimAPI_MakeBox(gp_Pnt(30, 30, 30), 1, 1, 1).Shape()
    r_res, r_exc = run_many(b, [hit, far], "difference")
    o_res, o_exc = run_many(b, [hit, far], "difference",
                                     optimized=True)
    o_out, o_sum = must_accept(o_res, o_exc)
    assert o_sum["n_skipped"] == 1
    assert_same_result(r_res, r_exc, o_res, o_exc)


def test_overlapping_tools_fold_difference():
    # Two overlapping hole tools: one interaction group -> one certified
    # tool-union, then a single difference.
    b = plate()
    t1 = BRepPrimAPI_MakeBox(gp_Pnt(2, 2, -0.25), 0.5, 0.5, 1.0).Shape()
    t2 = BRepPrimAPI_MakeBox(gp_Pnt(2.2, 2.2, -0.25), 0.5, 0.5, 1.0).Shape()
    r_res, r_exc = run_many(b, [t1, t2], "difference")
    o_res, o_exc = run_many(b, [t1, t2], "difference",
                                     optimized=True)
    o_out, o_sum = must_accept(o_res, o_exc)
    assert o_sum["n_groups"] == 1
    assert o_sum["n_combine_ops"] == 1
    assert_same_result(r_res, r_exc, o_res, o_exc)


def test_union_with_far_tool():
    b = plate()
    near = BRepPrimAPI_MakeBox(gp_Pnt(2, 2, -0.25), 0.4, 0.4, 1.0).Shape()
    far = BRepPrimAPI_MakeBox(gp_Pnt(20, 20, 20), 1, 1, 1).Shape()
    r_res, r_exc = run_many(b, [near, far], "union")
    o_res, o_exc = run_many(b, [near, far], "union",
                                     optimized=True)
    o_out, o_sum = must_accept(o_res, o_exc)
    assert o_sum["optimized"] is True
    assert o_sum["n_skipped"] == 0  # union never skips
    assert_same_result(r_res, r_exc, o_res, o_exc)
    assert count(o_out, TopAbs_SOLID) == 2  # disjoint union kept


def test_intersection_falls_back_to_reference():
    b = plate()
    t1 = BRepPrimAPI_MakeBox(gp_Pnt(1, 1, -0.25), 2.0, 2.0, 1.0).Shape()
    t2 = BRepPrimAPI_MakeBox(gp_Pnt(2, 2, -0.25), 2.0, 2.0, 1.0).Shape()
    r_res, r_exc = run_many(b, [t1, t2], "intersection")
    o_res, o_exc = run_many(b, [t1, t2], "intersection",
                                     optimized=True)
    o_out, o_sum = must_accept(o_res, o_exc)
    assert o_sum["optimized"] is False  # documented fallback
    assert_same_result(r_res, r_exc, o_res, o_exc)


def test_face_attribution_four_holes():
    b = plate()
    tools = [BRepPrimAPI_MakeBox(gp_Pnt(1 + 2 * i - 0.15, 1 + 2 * j - 0.15,
                                        -0.25), 0.3, 0.3, 1.0).Shape()
             for i in range(2) for j in range(2)]
    o_res, o_exc = run_many(b, tools, "difference",
                                     optimized=True)
    o_out, o_sum = must_accept(o_res, o_exc)
    attr = o_sum["face_attribution"]
    assert len(attr) == count(o_out, TopAbs_FACE)
    # 6 base faces unattributed, 4 wall faces per hole attributed.
    assert attr.count(-1) == 6
    for k in range(4):
        assert attr.count(k) == 4, f"tool {k}: {attr.count(k)} faces"
    # ancestry entries line up with the input tools, none skipped.
    anc = o_sum["tool_ancestry"]
    assert len(anc) == 4
    assert all(not e["skipped"] for e in anc)
    assert sorted(e["tool_index"] for e in anc) == [0, 1, 2, 3]


def test_shadow_reference_serves_reference():
    b = plate()
    tools = [BRepPrimAPI_MakeBox(gp_Pnt(1 + 2 * i - 0.15, 1.5 - 0.15,
                                        -0.25), 0.3, 0.3, 1.0).Shape()
             for i in range(3)]
    o_res, o_exc = run_many(b, tools, "difference",
                                     optimized=True, shadow_reference=True)
    o_out, o_sum = must_accept(o_res, o_exc)
    r_res, r_exc = run_many(b, tools, "difference")
    r_out, r_sum = must_accept(r_res, r_exc)
    assert o_sum["shadow_reference"] is True
    assert o_sum["shadow_matched"] is True
    assert vol(o_out) == vol(r_out)


SEEDS = list(range(20))


@pytest.mark.parametrize("seed", SEEDS)
def test_randomized_equivalence(seed):
    """Plan gate: 20 randomized tool sets match the sequential reference."""
    rng = random.Random(1000 + seed)
    op = rng.choice(["difference", "difference", "union", "intersection"])
    n = rng.randint(2, 5)
    tools = [random_tool(rng) for _ in range(n)]
    # Sometimes force two tools to overlap each other (fold path).
    if n >= 3 and rng.random() < 0.5:
        t = tools[0]
        tools[1] = BRepPrimAPI_MakeBox(gp_Pnt(3.9, 3.9, -0.2),
                                       0.5, 0.5, 1.0).Shape()
    b = plate()
    r_res, r_exc = run_many(b, tools, op)
    o_res, o_exc = run_many(b, tools, op, optimized=True)
    if o_exc is None:
        _, o_sum = o_res
        if op == "intersection":
            assert o_sum["optimized"] is False
        else:
            assert o_sum["optimized"] is True
    assert_same_result(r_res, r_exc, o_res, o_exc)


def test_tangent_tool_contact_conservative():
    # Box sitting exactly on the plate top: boxes touch -> never skipped,
    # both paths run the same pipeline op.
    b = plate()
    t = BRepPrimAPI_MakeBox(gp_Pnt(2, 2, 0.5), 1.0, 1.0, 1.0).Shape()
    r_res, r_exc = run_many(b, [t], "difference")
    o_res, o_exc = run_many(b, [t], "difference", optimized=True)
    o_out, o_sum = must_accept(o_res, o_exc)
    assert o_sum["n_skipped"] == 0
    assert_same_result(r_res, r_exc, o_res, o_exc)
