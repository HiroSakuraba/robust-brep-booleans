"""G12a: call-scoped assembly indexes.

(a) _EdgeSet / _AssemblyIndexes return exactly the same candidate sets
    as the linear scans they replace (pairwise predicate equality over
    a battery of real shapes).
(b) Indexes are call-scoped: no module-level or default-argument mutable
    state leaks between calls; repeated calls are identical.

Semantics-preserving only: this gate must not change any verdict.
"""
import inspect
import os
import sys
from types import SimpleNamespace

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import brepkernel.assembly as asm  # noqa: E402
from brepkernel.assembly import (  # noqa: E402
    _AssemblyIndexes,
    _EdgeSet,
    _edge_matches_ref_edge,
    _edge_matches_section,
    _edge_on_face_boundary,
    _shape_has_edge,
    _unique_edges,
)
from brepkernel import boolean_brep, BRepAmbiguousResult  # noqa: E402

from OCP.BRep import BRep_Tool  # noqa: E402
from OCP.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder,  # noqa: E402
                             BRepPrimAPI_MakeSphere)
from OCP.TopExp import TopExp_Explorer  # noqa: E402
from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID  # noqa: E402
from OCP.TopoDS import TopoDS  # noqa: E402
from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt  # noqa: E402
from OCP.BRepGProp import BRepGProp  # noqa: E402
from OCP.GProp import GProp_GProps  # noqa: E402


def _box(x0, y0, z0, x1, y1, z1):
    return BRepPrimAPI_MakeBox(gp_Pnt(x0, y0, z0), gp_Pnt(x1, y1, z1)).Shape()


def _cyl(px, py, pz, dx, dy, dz, r, h):
    return BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(px, py, pz), gp_Dir(dx, dy, dz)), r, h).Shape()


def _sph(cx, cy, cz, r):
    return BRepPrimAPI_MakeSphere(gp_Pnt(cx, cy, cz), r).Shape()


def _faces(shape):
    out = []
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        out.append(TopoDS.Face(ex.Current()))
        ex.Next()
    return out


def _has_solid(shape):
    return TopExp_Explorer(shape, TopAbs_SOLID).More()


def _volume(shape):
    if not _has_solid(shape):
        return 0.0
    g = GProp_GProps()
    BRepGProp.VolumePropertiesGK_s(shape, g, 1e-10, True, True, False, False, False)
    return float(g.Mass())


def _verdict(a, b, op):
    try:
        out, rep = boolean_brep(a, b, op)
        st = rep.get("stages", {}).get("assembly", {})
        lin = st.get("edge_lineage", {}).get("records", [])
        return ("accept", _volume(out), lin)
    except BRepAmbiguousResult as exc:
        r = (exc.report or {}).get("refusal", {}) or {}
        return ("refuse", r.get("kind"), r.get("stage"))
    except Exception as exc:  # noqa: BLE001
        return ("crash", type(exc).__name__, str(exc)[:100])


def _battery_shapes():
    shapes = [
        _box(0, 0, 0, 1, 1, 1),
        _box(0.5, 0.25, -0.5, 1.5, 1.25, 0.75),
        _cyl(0, 0, -1, 0, 0, 1, 0.5, 2.0),
        _sph(0, 0, 0, 0.8),
    ]
    ops = [
        (_box(0, 0, 0, 2, 2, 2), _box(1, 1, 1, 3, 3, 3), "union"),
        (_box(0, 0, 0, 2, 2, 2), _cyl(1, 1, -1, 0, 0, 1, 0.5, 4), "difference"),
        (_box(0, 0, 0, 2, 2, 2), _box(0, 0, 0, 1, 1, 1), "intersection"),
    ]
    for a, b, op in ops:
        try:
            out, _ = boolean_brep(a, b, op)
            shapes.append(out)
        except Exception:  # noqa: BLE001 - battery only needs successes
            pass
    return shapes


_sec_counter = [0]


def _sec_like(edge):
    _sec_counter[0] += 1
    return SimpleNamespace(
        edge=edge,
        verify_tolerance=float(BRep_Tool.Tolerance_s(edge)),
        edge_tolerance=float(BRep_Tool.Tolerance_s(edge)),
        face_a=-1, face_b=-1, edge_index=_sec_counter[0])


def test_edgeset_matches_shape_has_edge():
    shapes = _battery_shapes()
    n = 0
    # Pairwise property: _EdgeSet membership == _shape_has_edge, for
    # every (shape, subshape, edge) combination in the battery.
    for s1 in shapes:
        for s2 in shapes:
            es2 = _EdgeSet(s2)
            for sub in [s2] + _faces(s2):
                ess = es2 if sub is s2 else _EdgeSet(sub)
                for e in _unique_edges(s1):
                    n += 1
                    assert (e in ess) == _shape_has_edge(sub, e), \
                        "EdgeSet mismatch"
    print("edgeset pairs checked:", n)
    assert n > 200


def test_indexed_matchers_equal_linear():
    shapes = _battery_shapes()
    tol = 1e-7
    n = 0
    for rs in shapes:
        refs = []
        for other in shapes:
            for e in _unique_edges(other):
                refs.append(_sec_like(e))
        tools = [{"edge": s.edge, "operand": "T",
                  "face_id": j, "edge_id": j}
                 for j, s in enumerate(refs)]
        idx = _AssemblyIndexes(rs, [], refs, tools, {}, tol)
        edges = _unique_edges(rs)
        assert idx.n_result_edges == len(edges)
        for i, e in enumerate(edges):
            for j, sec in enumerate(refs):
                n += 1
                assert idx.edge_matches_section(i, j, tol) == \
                    _edge_matches_section(e, sec, tol), \
                    "section matcher mismatch at (%d, %d)" % (i, j)
            for j, sec in enumerate(refs):
                n += 1
                assert idx.edge_matches_tool(i, j, tol) == \
                    _edge_matches_ref_edge(
                        e, sec.edge, tol, tol, tol), \
                    "tool matcher mismatch at (%d, %d)" % (i, j)
    print("matcher pairs checked:", n)
    assert n > 500


def test_indexed_face_boundary_equal_linear():
    shapes = _battery_shapes()
    tol = 1e-7
    n = 0
    for rs in shapes:
        faces = _faces(rs)
        if not faces:
            continue
        idx = _AssemblyIndexes(rs, [], [], [], {}, tol)
        for fi, f in enumerate(faces):
            idx.register_parent_face(("T", fi), f)
        for i, e in enumerate(_unique_edges(rs)):
            for fi, f in enumerate(faces):
                n += 1
                assert idx.edge_on_face_boundary(
                    i, ("T", fi), tol) == _edge_on_face_boundary(e, f, tol), \
                    "face-boundary mismatch at (%d, %d)" % (i, fi)
    print("face-boundary pairs checked:", n)
    assert n > 100


def test_no_module_level_state_leaks():
    before_names = set(vars(asm))
    cases = [
        (_box(0, 0, 0, 2, 2, 2), _box(1, 1, 1, 3, 3, 3), "union"),
        (_box(0, 0, 0, 2, 2, 2), _cyl(1, 1, -1, 0, 0, 1, 0.5, 4), "difference"),
        (_box(0, 0, 0, 1, 1, 1), _box(1, 0, 0, 2, 1, 1), "union"),
    ]
    first = None
    for a, b, op in cases:
        v = _verdict(a, b, op)
        if first is None:
            first = v
    # Same op twice: identical verdict, volume, lineage sizes.
    r1 = _verdict(*cases[0])
    r2 = _verdict(*cases[0])
    assert r1 == r2, "repeat call diverged: %r vs %r" % (r1, r2)
    after_names = set(vars(asm))
    assert after_names == before_names, \
        "new module attrs: %s" % sorted(after_names - before_names)
    for v in vars(asm).values():
        assert not isinstance(v, (_EdgeSet, _AssemblyIndexes)), \
            "index object leaked to module level"
    # No mutable default-argument caches on the new helpers.
    for fn in (_EdgeSet.__init__, _AssemblyIndexes.__init__,
               _AssemblyIndexes.edge_matches_section,
               _AssemblyIndexes.edge_matches_tool,
               _AssemblyIndexes.edge_on_face_boundary):
        for p in inspect.signature(fn).parameters.values():
            d = p.default
            assert not isinstance(d, (list, dict, set, bytearray)), \
                "mutable default on %s.%s" % (fn.__qualname__, p.name)


def main():
    test_edgeset_matches_shape_has_edge()
    test_indexed_matchers_equal_linear()
    test_indexed_face_boundary_equal_linear()
    test_no_module_level_state_leaks()
    print("ALL G12a TESTS PASS")


if __name__ == "__main__":
    main()
