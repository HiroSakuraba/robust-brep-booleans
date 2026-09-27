"""S1 PreparedBRep tests.

Pass criteria from the speed plan:
- raw-shape and prepared-shape calls are verdict-equivalent on every
  test and fuzz battery;
- a second Boolean using the same prepared base performs zero
  face-box rebuilds and zero edge-index rebuilds for that base;
- preparing the model does not tessellate or mutate the input shape.
"""
import sys

sys.path.insert(0, "src")
sys.path.insert(0, "tests")

import numpy as np

from brepkernel import boolean_brep, BRepAmbiguousResult
from brepkernel import perf as perf_mod
from brepkernel.prepared import prepare_brep, prepare_model, ensure_prepared

from OCP.BRep import BRep_Tool
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
from OCP.GProp import GProp_GProps
from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID
from OCP.TopExp import TopExp_Explorer
from OCP.TopLoc import TopLoc_Location
from OCP.TopoDS import TopoDS
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


def ev_name(report):
    ev = report.get("evidence", {})
    recs = ev.get("records", [])
    return recs[0].get("name") if recs else None


def t1_raw_prepared_equivalence():
    ok = True
    cases = [("boxes", boxes()), ("cyl_box", cyl_box())]
    for cname, (sa, sb) in cases:
        for op in ("union", "intersection", "difference"):
            r_raw = run(sa, sb, op)
            pa, pb = prepare_brep(sa), prepare_brep(sb)
            r_pre = run(pa, pb, op)
            r_mix = run(sa, pb, op)
            same = (r_raw[:4] == r_pre[:4] == r_mix[:4])
            ok &= check(f"t1 {cname}/{op}: raw==prepared==mixed",
                        same, f"raw={r_raw[:4]} pre={r_pre[:4]}")
            if r_raw[0]:
                ok &= check(f"t1 {cname}/{op}: volume identical",
                            r_raw[1] == r_pre[1] == r_mix[1],
                            f"vol={r_raw[1]!r}")
                ok &= check(f"t1 {cname}/{op}: evidence name identical",
                            ev_name(r_raw[4]) == ev_name(r_pre[4]),
                            f"{ev_name(r_raw[4])}")
    return ok


def t2_zero_rebuilds_on_reuse():
    sa, sb = boxes()
    ok, _, _, _, raw_rep = run(sa, sb, "union", collect_perf=True)
    raw_c = raw_rep["performance"]["counters"]
    ok = check("t2 raw run accepted", ok)
    # The raw call prepares internally (inside the perf scope), so its
    # one face-box build is counted; nothing may rebuild after that.
    ok &= check("t2 raw run built face boxes once",
                raw_c.get("face_box_build", 0) >= 1,
                f"={raw_c.get('face_box_build', 0)}")

    pa, pb = prepare_brep(sa), prepare_brep(sb)
    for i in (1, 2):
        ok_i, _, _, _, rep = run(pa, pb, "union", collect_perf=True)
        c = rep["performance"]["counters"]
        ok &= check(f"t2 prepared run {i} accepted", ok_i)
        ok &= check(f"t2 prepared run {i}: zero face-box rebuilds",
                    c.get("face_box_build", -1) == 0,
                    f"={c.get('face_box_build', -1)}")
        ok &= check(f"t2 prepared run {i}: face-box hits recorded",
                    c.get("prepared_face_box_hit", 0) >= 1,
                    f"={c.get('prepared_face_box_hit', 0)}")
        ok &= check(f"t2 prepared run {i}: edge-index hits recorded",
                    c.get("prepared_edge_index_hit", 0) >= 1,
                    f"={c.get('prepared_edge_index_hit', 0)}")
    return ok


def t3_ensure_prepared_semantics():
    from brepkernel.step_ingest import BRepModel, index_shape
    sa, _ = boxes()
    pa = prepare_brep(sa)
    ok = check("t3 PreparedBRep passes through identical",
               ensure_prepared(pa) is pa)
    m = index_shape(sa)
    pm = ensure_prepared(m)
    ok &= check("t3 BRepModel wrapped without re-indexing",
                pm.model is m)
    ok &= check("t3 raw shape prepared fresh",
                ensure_prepared(sa).model is not m)
    ok &= check("t3 wrapped model verdict-equivalent",
                run(pm, prepare_brep(sa), "union")[0])
    return ok


def _triangulated_faces(shape):
    out = []
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        loc = TopLoc_Location()
        if BRep_Tool.Triangulation_s(TopoDS.Face(ex.Current()), loc) is not None:
            out.append(True)
        ex.Next()
    return out


def _face_count(shape):
    n = 0
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        n += 1
        ex.Next()
    return n


def t4_prepare_does_not_tessellate_or_mutate():
    sa, sb = cyl_box()
    n_a, n_b = _face_count(sa), _face_count(sb)
    before_a = _triangulated_faces(sa)
    before_b = _triangulated_faces(sb)
    vol_a, vol_b = volume(sa), volume(sb)
    pa, pb = prepare_brep(sa), prepare_brep(sb)
    after_a = _triangulated_faces(sa)
    after_b = _triangulated_faces(sb)
    ok = check("t4 no face triangulated by prepare (A)",
               before_a == after_a == [],
               f"faces={n_a}")
    ok &= check("t4 no face triangulated by prepare (B)",
                before_b == after_b == [],
                f"faces={n_b}")
    ok &= check("t4 input volumes unchanged",
                volume(sa) == vol_a and volume(sb) == vol_b)
    ok &= check("t4 prepared model face counts sane",
                pa.n_faces == n_a and pb.n_faces == n_b,
                f"{pa.n_faces}/{pb.n_faces}")
    return ok

def t5_face_adjacency_properties():
    sa, _ = boxes()
    pa = prepare_brep(sa)
    adj = pa.face_adjacency
    n = pa.n_faces
    ok = check("t5 one entry per face", len(adj) == n, f"n={n}")
    ok &= check("t5 box faces each adjacent to 4",
                all(len(a) == 4 for a in adj),
                f"{sorted(len(a) for a in adj)}")
    sym = all(i in adj[j] for i in range(n) for j in adj[i])
    ok &= check("t5 adjacency symmetric", sym)
    # Independent recomputation: naive all-pairs IsSame over face edges.
    from OCP.TopAbs import TopAbs_EDGE
    face_edges = []
    for fr in pa.model.faces:
        edges = []
        ex = TopExp_Explorer(fr.face, TopAbs_EDGE)
        while ex.More():
            edges.append(TopoDS.Edge(ex.Current()))
            ex.Next()
        face_edges.append(edges)
    expect = [set() for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if any(ea.IsSame(eb) for ea in face_edges[i]
                   for eb in face_edges[j]):
                expect[i].add(j)
                expect[j].add(i)
    match = all(set(adj[i]) == expect[i] for i in range(n))
    ok &= check("t5 matches naive IsSame recomputation", match)
    return ok


def t6_consumers_use_prepared_without_rebuild():
    from brepkernel.assembly import (_MultiRayClassifier,
                                     _point_boundary_distances)
    sa, _ = boxes()
    pa = prepare_brep(sa)
    model = pa.model
    pts = np.array([[0.5, 0.5, 2.0]])
    with perf_mod.scoped(perf_mod.PerfCounters()):
        _point_boundary_distances(pts, model, cap=1.0,
                                  face_boxes=pa.face_boxes)
        snap = perf_mod.current().snapshot()
    ok = check("t6 prepared face boxes: no rebuild",
               snap.get("face_box_build", -1) == 0)
    ok &= check("t6 prepared face boxes: hit recorded",
                snap.get("prepared_face_box_hit", 0) == 1)
    with perf_mod.scoped(perf_mod.PerfCounters()):
        _point_boundary_distances(pts, model, cap=1.0)
        snap2 = perf_mod.current().snapshot()
    ok &= check("t6 fallback still builds without prepared",
                snap2.get("face_box_build", 0) >= 1)
    with perf_mod.scoped(perf_mod.PerfCounters()):
        _MultiRayClassifier([sr.solid for sr in model.solids], 1e-7,
                            edges=pa.edges, edge_boxes=pa.edge_boxes)
        snap3 = perf_mod.current().snapshot()
    ok &= check("t6 prepared edges: no rebuild",
                snap3.get("edge_box_build", -1) == 0)
    ok &= check("t6 prepared edges: hit recorded",
                snap3.get("prepared_edge_index_hit", 0) == 1)
    with perf_mod.scoped(perf_mod.PerfCounters()):
        _MultiRayClassifier([sr.solid for sr in model.solids], 1e-7)
        snap4 = perf_mod.current().snapshot()
    ok &= check("t6 classifier fallback still builds",
                snap4.get("edge_box_build", 0) >= 1)
    return ok


def t7_prepared_boxes_bit_identical():
    # The prepared boxes must equal what the per-call path builds.
    from brepkernel.assembly import _conservative_boxes
    sa, _ = boxes()
    pa = prepare_brep(sa)
    faces = [fr.face for fr in pa.model.faces]
    rebuilt = _conservative_boxes(faces)
    ok = check("t7 face boxes bit-identical to per-call build",
               np.array_equal(pa.face_boxes, rebuilt))
    ok &= check("t7 pads identical to face_broadphase_pads",
                np.array_equal(pa.face_pads(4e-7),
                               __import__("brepkernel.step_ingest",
                                          fromlist=["face_broadphase_pads"]
                                          ).face_broadphase_pads(
                                   pa.model, 4e-7)))
    return ok


TESTS = [t1_raw_prepared_equivalence, t2_zero_rebuilds_on_reuse,
         t3_ensure_prepared_semantics, t4_prepare_does_not_tessellate_or_mutate,
         t5_face_adjacency_properties, t6_consumers_use_prepared_without_rebuild,
         t7_prepared_boxes_bit_identical]


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
