"""S3 whole-operation fast paths.

Covers: separated battery matches the old pipeline on all 3 ops,
near-touching/touching never fires, nested box/cylinder/sphere/NURBS
batteries give the same verdicts and volumes, multi-component
containment, mixed configurations fall back, shadow mode agrees, and
the flags-off path is unchanged.
"""

import pytest
from OCP.BRep import BRep_Builder
from OCP.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepPrimAPI import (
    BRepPrimAPI_MakeBox,
    BRepPrimAPI_MakeCylinder,
    BRepPrimAPI_MakeSphere,
)
from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt
from OCP.TopAbs import TopAbs_FACE, TopAbs_SHELL, TopAbs_SOLID
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS_Compound

from brepkernel.assembly import _shape_volume
from brepkernel.pipeline import BRepAmbiguousResult, boolean_brep


def _box(x, y, z, sx, sy=None, sz=None):
    sy = sx if sy is None else sy
    sz = sx if sz is None else sz
    return BRepPrimAPI_MakeBox(
        gp_Pnt(float(x), float(y), float(z)),
        float(sx), float(sy), float(sz)).Shape()


def _cylinder(r, h, x=0.0, y=0.0, z=0.0):
    return BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(float(x), float(y), float(z)), gp_Dir(0, 0, 1)),
        float(r), float(h)).Shape()


def _sphere(r, x=0.0, y=0.0, z=0.0):
    return BRepPrimAPI_MakeSphere(
        gp_Pnt(float(x), float(y), float(z)), float(r)).Shape()


def _to_nurbs(shape):
    conv = BRepBuilderAPI_NurbsConvert(shape, True)
    assert conv.IsDone()
    out = conv.Shape()
    assert BRepCheck_Analyzer(out, True).IsValid()
    return out


def _compound(shapes):
    builder = BRep_Builder()
    comp = TopoDS_Compound()
    builder.MakeCompound(comp)
    for s in shapes:
        builder.Add(comp, s)
    return comp


def _counts(shape):
    out = []
    for kind in (TopAbs_SOLID, TopAbs_SHELL, TopAbs_FACE):
        ex = TopExp_Explorer(shape, kind)
        n = 0
        while ex.More():
            n += 1
            ex.Next()
        out.append(n)
    return tuple(out)


def _volume(shape):
    if shape.IsNull() or _counts(shape)[0] == 0:
        return 0.0
    return abs(float(_shape_volume(shape)))


def _run(a, b, op, **kw):
    kw.setdefault("crosscheck_ops", False)
    return boolean_brep(a, b, op, **kw)


def _assert_same_result(out_new, rep_new, out_old, rep_old):
    assert rep_new["accepted"] and rep_old["accepted"]
    assert _counts(out_new) == _counts(out_old)
    vn, vo = _volume(out_new), _volume(out_old)
    assert abs(vn - vo) <= 1e-9 * max(1.0, abs(vo), abs(vn))


def _s3_stage(rep):
    return rep["stages"].get("s3_fast_path", {})


# --------------------------------------------------------------------------
# Fast path A: separated components

SEPARATED_CONFIGS = [
    ("boxes", _box(0, 0, 0, 1.0), _box(5, 0, 0, 1.0)),
    ("box_cylinder", _box(0, 0, 0, 1.0), _cylinder(0.5, 2.0, x=4.0)),
    ("offset_diag", _box(0, 0, 0, 2.0), _box(10, -7, 3, 1.5)),
    ("multi_a", _compound([_box(0, 0, 0, 1.0), _box(3, 0, 0, 1.0)]),
     _box(8, 8, 8, 2.0)),
]


@pytest.mark.parametrize("name,a,b", SEPARATED_CONFIGS)
@pytest.mark.parametrize("op", ["intersection", "difference", "union"])
def test_s3a_separated_matches_old_path(name, a, b, op):
    out_new, rep_new = _run(a, b, op, fast_paths=True,
                            fast_path_shadow=False)
    out_old, rep_old = _run(a, b, op, fast_paths=False)
    _assert_same_result(out_new, rep_new, out_old, rep_old)
    st = _s3_stage(rep_new)
    assert st.get("path") == "A_disjoint"
    assert "s3_internal_error" not in st


@pytest.mark.parametrize("op", ["intersection", "difference", "union"])
def test_s3a_expected_volumes(op):
    a, b = _box(0, 0, 0, 2.0), _box(10, 0, 0, 1.0)
    out, rep = _run(a, b, op, fast_paths=True, fast_path_shadow=False)
    v = _volume(out)
    expected = {"intersection": 0.0, "difference": 8.0, "union": 9.0}[op]
    assert abs(v - expected) < 1e-9
    assert _s3_stage(rep)["resolution"] == {
        "intersection": "empty", "difference": "A",
        "union": "union"}[op]


def test_s3a_touching_boxes_do_not_fire():
    a, b = _box(0, 0, 0, 1.0), _box(1, 0, 0, 1.0)  # share a face
    for op in ("intersection", "difference", "union"):
        out_new, rep_new = _run(a, b, op, fast_paths=True,
                                fast_path_shadow=False)
        st = _s3_stage(rep_new)
        assert st.get("path") is None, st
        assert st.get("reason") == "candidates_exist"
        out_old, rep_old = _run(a, b, op, fast_paths=False)
        _assert_same_result(out_new, rep_new, out_old, rep_old)


def test_s3a_gap_outside_pads_fires():
    # 1e-4 gap: far outside the default ~5e-7 broad-phase pads, so no
    # candidates exist and the fast path correctly fires (the old
    # pipeline would also find nothing to section).
    a, b = _box(0, 0, 0, 1.0), _box(1 + 1e-4, 0, 0, 1.0)
    for op in ("intersection", "difference", "union"):
        out_new, rep_new = _run(a, b, op, fast_paths=True,
                                fast_path_shadow=False)
        assert _s3_stage(rep_new).get("path") == "A_disjoint"
        out_old, rep_old = _run(a, b, op, fast_paths=False)
        _assert_same_result(out_new, rep_new, out_old, rep_old)


def test_s3a_near_touching_gap_does_not_fire():
    # 1e-4 gap with an explicit 1e-3 broad-phase pad: the gap is far
    # outside the coincidence band (so the old pipeline proceeds) but
    # inside the pads -> candidates exist -> the fast path must not fire.
    a, b = _box(0, 0, 0, 1.0), _box(1 + 1e-4, 0, 0, 1.0)
    kw = dict(broadphase_pad=1e-3)
    for op in ("intersection", "difference", "union"):
        out_new, rep_new = _run(a, b, op, fast_paths=True,
                                fast_path_shadow=False, **kw)
        st = _s3_stage(rep_new)
        assert st.get("path") is None, st
        assert st.get("reason") == "candidates_exist"
        out_old, rep_old = _run(a, b, op, fast_paths=False, **kw)
        _assert_same_result(out_new, rep_new, out_old, rep_old)


# --------------------------------------------------------------------------
# Fast path B: containment with zero broad-phase candidates

def test_s3b_nested_box_all_ops():
    a, b = _box(0, 0, 0, 4.0), _box(1, 1, 1, 1.0)
    for op, res, vol in (("intersection", "B", 1.0),
                         ("union", "A", 64.0)):
        out, rep = _run(a, b, op, fast_paths=True, fast_path_shadow=False)
        st = _s3_stage(rep)
        assert st.get("path") == "B_containment", st
        assert st.get("relation") == "B_in_A"
        assert st.get("resolution") == res
        assert abs(_volume(out) - vol) < 1e-9
        assert "s3_internal_error" not in st
    # difference with B_in_A needs real face splitting -> falls back.
    out, rep = _run(a, b, "difference", fast_paths=True,
                    fast_path_shadow=False)
    assert _s3_stage(rep).get("path") is None
    assert _s3_stage(rep).get("reason") == "unsupported_difference_cavity"
    out_old, rep_old = _run(a, b, "difference", fast_paths=False)
    _assert_same_result(out, rep, out_old, rep_old)
    assert abs(_volume(out) - 63.0) < 1e-9


def test_s3b_a_in_b():
    a, b = _box(1, 1, 1, 1.0), _box(0, 0, 0, 4.0)  # A inside B
    for op, res, vol in (("intersection", "A", 1.0),
                         ("union", "B", 64.0),
                         ("difference", "empty", 0.0)):
        out, rep = _run(a, b, op, fast_paths=True, fast_path_shadow=False)
        st = _s3_stage(rep)
        assert st.get("path") == "B_containment", st
        assert st.get("relation") == "A_in_B"
        assert st.get("resolution") == res
        assert abs(_volume(out) - vol) < 1e-9


@pytest.mark.parametrize("inner", [
    _box(1, 1, 1, 1.0),
    _cylinder(0.4, 1.0, x=2.0, y=2.0, z=1.5),
    _sphere(0.5, x=2.0, y=2.0, z=2.0),
    _to_nurbs(_box(1, 1, 1, 1.0)),
])
@pytest.mark.parametrize("op", ["intersection", "union"])
def test_s3b_nested_geometry_battery(inner, op):
    outer = _box(0, 0, 0, 4.0)
    out_new, rep_new = _run(outer, inner, op, fast_paths=True,
                            fast_path_shadow=False)
    out_old, rep_old = _run(outer, inner, op, fast_paths=False)
    _assert_same_result(out_new, rep_new, out_old, rep_old)
    st = _s3_stage(rep_new)
    assert st.get("path") == "B_containment", st
    assert st.get("relation") == "B_in_A"
    expected_vol = 64.0 if op == "union" else _volume(inner)
    assert abs(_volume(out_new) - expected_vol) < 1e-6


def test_s3b_multicomponent_containment():
    outer = _box(0, 0, 0, 6.0)
    inner = _compound([_box(1, 1, 1, 1.0), _box(4, 4, 4, 1.0)])
    for op, vol, solids in (("intersection", 2.0, 2),
                            ("union", 216.0, 1)):
        out, rep = _run(outer, inner, op, fast_paths=True,
                        fast_path_shadow=False)
        st = _s3_stage(rep)
        assert st.get("path") == "B_containment", st
        assert st.get("relation") == "B_in_A"
        assert abs(_volume(out) - vol) < 1e-9
        assert _counts(out)[0] == solids


def _ell_solid():
    # Single L-shaped solid: its box covers the notch, but no face box
    # reaches a small box parked in the notch -> B/disjoint exercises.
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse

    fuse = BRepAlgoAPI_Fuse(_box(0, 0, 0, 4, 1, 1), _box(0, 0, 0, 1, 4, 1))
    fuse.Build()
    assert fuse.IsDone()
    return fuse.Shape()


def test_s3b_disjoint_overlapping_boxes():
    # Small box parked in the L notch: component boxes overlap, zero
    # face candidates, solids disjoint.
    ell = _ell_solid()
    parked = _box(2, 2, 0, 1.0, 1.0, 1.0)
    for op, res, vol in (("intersection", "empty", 0.0),
                         ("difference", "A", None),
                         ("union", "union", None)):
        out_new, rep_new = _run(ell, parked, op, fast_paths=True,
                                fast_path_shadow=False)
        out_old, rep_old = _run(ell, parked, op, fast_paths=False)
        _assert_same_result(out_new, rep_new, out_old, rep_old)
        st = _s3_stage(rep_new)
        assert st.get("path") == "B_containment", st
        assert st.get("relation") == "disjoint"
        assert st.get("resolution") == res
        if vol is not None:
            assert abs(_volume(out_new) - vol) < 1e-9


def test_s3_mixed_configuration_falls_back():
    # One component of A overlaps B, another is far away -> candidates.
    a = _compound([_box(0, 0, 0, 1.0), _box(20, 0, 0, 1.0)])
    b = _box(0.5, 0, 0, 1.0)
    for op in ("intersection", "difference", "union"):
        out_new, rep_new = _run(a, b, op, fast_paths=True,
                                fast_path_shadow=False)
        st = _s3_stage(rep_new)
        assert st.get("path") is None, st
        assert st.get("reason") == "candidates_exist"
        out_old, rep_old = _run(a, b, op, fast_paths=False)
        _assert_same_result(out_new, rep_new, out_old, rep_old)


# --------------------------------------------------------------------------
# Shadow mode and flags

def test_s3_shadow_agrees_on_battery():
    configs = [
        (_box(0, 0, 0, 1.0), _box(5, 0, 0, 1.0)),       # A_disjoint
        (_box(0, 0, 0, 4.0), _box(1, 1, 1, 1.0)),       # B_containment
        (_box(0, 0, 0, 1.0), _box(0.5, 0, 0, 1.0)),     # overlap -> no fire
    ]
    for a, b in configs:
        for op in ("intersection", "difference", "union"):
            out, rep = _run(a, b, op, fast_paths=True,
                            fast_path_shadow=True)
            st = _s3_stage(rep)
            assert rep["accepted"]
            assert "s3_internal_error" not in st
            if st.get("path") is not None:
                assert st.get("shadow") is True
                assert st.get("shadow_agreement") is True
            out_old, rep_old = _run(a, b, op, fast_paths=False)
            _assert_same_result(out, rep, out_old, rep_old)


def test_s3_shadow_old_refusal_propagates(monkeypatch):
    # R1: if the old pipeline refuses while the fast path fired, the
    # shadow must refuse too (the fast result is discarded).  Simulated
    # by making the shadow invocation raise.
    import brepkernel.pipeline as pl

    real_impl = pl._boolean_brep_impl

    def fake_impl(shapeA, shapeB, op, **kw):
        if kw.get("fast_paths") is False and kw.get("fast_path_shadow") is False:
            raise BRepAmbiguousResult(
                "simulated old-pipeline refusal",
                {"stages": {}, "accepted": False}, None)
        return real_impl(shapeA, shapeB, op, **kw)

    monkeypatch.setattr(pl, "_boolean_brep_impl", fake_impl)
    a, b = _box(0, 0, 0, 1.0), _box(5, 0, 0, 1.0)
    with pytest.raises(BRepAmbiguousResult) as ei:
        _run(a, b, "union", fast_paths=True, fast_path_shadow=True)
    st = ei.value.report["stages"].get("s3_fast_path", {})
    assert st.get("path") == "A_disjoint"
    assert st.get("shadow_agreement") is False


def test_s3_disabled_leaves_no_stage():
    a, b = _box(0, 0, 0, 1.0), _box(5, 0, 0, 1.0)
    out, rep = _run(a, b, "union", fast_paths=False)
    assert "s3_fast_path" not in rep["stages"]
    assert rep["accepted"]


def test_s3_overlapping_unchanged_with_flags_on():
    # Ordinary intersecting boxes: no fast path, same result as before.
    a, b = _box(0, 0, 0, 2.0), _box(1, 1, 1, 2.0)
    for op in ("intersection", "difference", "union"):
        out_new, rep_new = _run(a, b, op, fast_paths=True,
                                fast_path_shadow=True)
        st = _s3_stage(rep_new)
        assert st.get("path") is None, st
        out_old, rep_old = _run(a, b, op, fast_paths=False)
        _assert_same_result(out_new, rep_new, out_old, rep_old)
