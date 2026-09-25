"""End-to-end B-rep assembly regressions for the Tier B/C path.

The tests exercise:
1. transverse sphere/sphere union, intersection and A-B;
2. disjoint union as two material solids;
3. contained subtraction as one solid with one cavity shell;
4. coincident/tangent ambiguity refusing rather than guessing.

The OCCT boolean engine is used only as a regression oracle for volume; the
assembled result itself comes from:
index -> verified intersection -> local split -> patch classify -> sew/build.
"""
import math
import sys

sys.path.insert(0, "src")

from brepkernel.assembly import (
    assemble_boolean, AssemblyError, _build_nested_solids,
    _classify_pieces, _shell_records,
)
from brepkernel.intersection import intersect_models
from brepkernel.split import (
    FaceSplitResult, ModelSplitResult, SplitFacePiece, split_models,
)
from brepkernel.step_ingest import index_shape

from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeSphere
from OCP.GProp import GProp_GProps
from OCP.TopAbs import TopAbs_SHELL
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Pnt


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def vol(shape):
    p = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, p)
    return float(p.Mass())


def pipeline(sa, sb, op, **kw):
    a = index_shape(sa)
    b = index_shape(sb)
    ix = intersect_models(a, b, broadphase_pad=kw.get("broadphase_pad", 0.0),
                          base_tol=1e-7, chord_tol=1e-5,
                          contact_tol=kw.get("contact_tol", None))
    sp = split_models(a, b, ix, base_tol=1e-7)
    return a, b, ix, sp, assemble_boolean(a, b, sp, op, base_tol=1e-7)


def oracle(sa, sb, op):
    if op == "union":
        x = BRepAlgoAPI_Fuse(sa, sb)
    elif op == "intersection":
        x = BRepAlgoAPI_Common(sa, sb)
    else:
        x = BRepAlgoAPI_Cut(sa, sb)
    x.Build()
    assert x.IsDone()
    return x.Shape()


def t1_overlap_spheres_all_ops():
    a = BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), 1.0).Shape()
    b = BRepPrimAPI_MakeSphere(gp_Pnt(1, 0, 0), 1.0).Shape()
    expected = {
        "intersection": 5.0 * math.pi / 12.0,
        "union": 9.0 * math.pi / 4.0,
        "difference": 11.0 * math.pi / 12.0,
    }
    ok = True
    for op in ("union", "intersection", "difference"):
        ma, mb, ix, sp, r = pipeline(a, b, op)
        ov = vol(oracle(a, b, op))
        err = abs(r.volume - expected[op])
        oerr = abs(r.volume - ov)
        ok &= check(
            f"t1 {op} valid closed result",
            r.free_edges == 0 and r.multiple_edges == 0
            and BRepCheck_Analyzer(r.shape, True).IsValid(),
            f"faces={r.selected_faces} shells={len(r.shells)} "
            f"solids={len(r.solids)}")
        ok &= check(
            f"t1 {op} analytic volume",
            err < 2e-6,
            f"got={r.volume:.12g} expected={expected[op]:.12g} "
            f"err={err:.3e}")
        ok &= check(
            f"t1 {op} agrees OCCT oracle",
            oerr < 2e-6,
            f"assembled={r.volume:.12g} oracle={ov:.12g} "
            f"err={oerr:.3e}")
        ok &= check(
            f"t1 {op} provenance decisions",
            len(r.decisions) >= 4
            and r.selected_faces == sum(int(d.keep) for d in r.decisions)
            and any(d.operand == "A" and d.keep for d in r.decisions)
            and any(d.operand == "B" and d.keep for d in r.decisions)
            and all(d.classification in ("inside", "outside")
                    for d in r.decisions)
            and (op != "difference"
                 or any(d.operand == "B" and d.keep
                        and d.reverse_for_difference for d in r.decisions)),
            f"decisions={[(d.operand,d.classification,d.keep) for d in r.decisions]}")
        section_lineage = [
            e for e in r.edge_lineage
            if e.provenance_kind == "boolean_section"]
        ok &= check(
            f"t1 {op} verified edge lineage",
            bool(section_lineage)
            and all(e.verified_pcurves and e.intersection_refs
                    for e in section_lineage)
            and any(set(e.operands) == {"A", "B"}
                    for e in section_lineage),
            f"lineage={[(e.result_edge_index,e.provenance_kind,e.operands,e.intersection_refs) for e in r.edge_lineage]}")
        payload_keys = {
            (p.face_a, p.face_b, p.section_edge_index)
            for p in r.section_payloads}
        lineage_refs = {
            ref for e in section_lineage for ref in e.intersection_refs}
        ok &= check(
            f"t1 {op} full section payload persistence",
            bool(r.section_payloads)
            and lineage_refs.issubset(payload_keys)
            and all(
                len(p.parameters) == len(p.xyz)
                == len(p.uv_a) == len(p.uv_b)
                and p.xyz.ndim == 2 and p.xyz.shape[1] == 3
                and p.uv_a.ndim == 2 and p.uv_a.shape[1] == 2
                and p.uv_b.ndim == 2 and p.uv_b.shape[1] == 2
                and len(p.parameters) >= 2
                for p in r.section_payloads)
            and any(p.result_edge_indices for p in r.section_payloads),
            f"payloads={[(p.face_a,p.face_b,p.section_edge_index,len(p.parameters),p.result_edge_indices) for p in r.section_payloads]}")
    return ok


def t2_disjoint_union_two_solids():
    a = BRepPrimAPI_MakeSphere(gp_Pnt(-2, 0, 0), 1.0).Shape()
    b = BRepPrimAPI_MakeSphere(gp_Pnt(2, 0, 0), 1.0).Shape()
    ma, mb, ix, sp, r = pipeline(a, b, "union")
    want = 8.0 * math.pi / 3.0
    ok = check("t2 broadphase skipped disjoint pair",
               ix.candidate_pairs == 0 and ix.section_calls == 0,
               f"candidates={ix.candidate_pairs} calls={ix.section_calls}")
    ok &= check("t2 disjoint union has two solids",
                len(r.solids) == 2 and len(r.shells) == 2,
                f"shells={len(r.shells)} solids={len(r.solids)}")
    ok &= check("t2 disjoint volume",
                abs(r.volume - want) < 2e-6,
                f"got={r.volume:.12g} expected={want:.12g}")
    ok &= check("t2 no invented section lineage",
                all(not e.intersection_refs and not e.verified_pcurves
                    for e in r.edge_lineage),
                f"lineage={[(e.result_edge_index,e.provenance_kind) for e in r.edge_lineage]}")
    ok &= check("t2 no invented section payloads",
                len(r.section_payloads) == 0,
                f"payloads={len(r.section_payloads)}")
    return ok


def t3_contained_difference_builds_cavity():
    outer = BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), 2.0).Shape()
    inner = BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), 0.5).Shape()
    ma, mb, ix, sp, r = pipeline(outer, inner, "difference")
    want = 4.0 * math.pi / 3.0 * (8.0 - 0.125)
    ok = check("t3 no intersection needed for containment",
               ix.verified_edges == 0,
               f"edges={ix.verified_edges}")
    ok &= check("t3 one material solid with cavity",
                len(r.solids) == 1
                and len(r.shells) == 2
                and len(r.solids[0].cavity_shells) == 1,
                f"shells={[(x.shell_index,x.depth,x.parent_shell) for x in r.shells]} "
                f"cavities={r.solids[0].cavity_shells if r.solids else None}")
    ok &= check("t3 cavity volume",
                abs(r.volume - want) < 3e-6,
                f"got={r.volume:.12g} expected={want:.12g}")
    ov = vol(oracle(outer, inner, "difference"))
    ok &= check("t3 cavity agrees OCCT oracle",
                abs(r.volume - ov) < 3e-6,
                f"assembled={r.volume:.12g} oracle={ov:.12g}")
    return ok


def t4_disjoint_intersection_is_empty():
    a = BRepPrimAPI_MakeSphere(gp_Pnt(-2, 0, 0), 1.0).Shape()
    b = BRepPrimAPI_MakeSphere(gp_Pnt(2, 0, 0), 1.0).Shape()
    ma, mb, ix, sp, r = pipeline(a, b, "intersection")
    return check("t4 disjoint intersection empty",
                 r.is_empty and r.volume == 0.0 and len(r.solids) == 0,
                 f"selected={r.selected_faces} volume={r.volume}")


def t5_exact_tangency_refuses_global_assembly():
    a = BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), 1.0).Shape()
    b = BRepPrimAPI_MakeSphere(gp_Pnt(2, 0, 0), 1.0).Shape()
    ma = index_shape(a)
    mb = index_shape(b)
    ix = intersect_models(ma, mb, broadphase_pad=1e-8,
                          base_tol=1e-7, contact_tol=1e-6)
    sp = split_models(ma, mb, ix, base_tol=1e-7)
    try:
        assemble_boolean(ma, mb, sp, "union", base_tol=1e-7)
    except AssemblyError as e:
        return check("t5 tangent assembly refuses",
                     getattr(e, "kind", "") == "UnresolvedContact",
                     f"kind={getattr(e, 'kind', '?')}")
    return check("t5 tangent assembly refuses", False, "no refusal")



def t6_unsplit_straddling_patch_refuses_multiwitness():
    """An intentionally unsplit face crossing another solid must not be
    classified from one lucky interior point.

    This simulates the dangerous downstream state caused by a missed section:
    the full sphere face spans both inside and outside of a box occupying the
    x >= 0.2 half-space. Multi-witness classification must refuse.
    """
    sphere = BRepPrimAPI_MakeSphere(
        gp_Pnt(0, 0, 0), 1.0).Shape()
    cutter = BRepPrimAPI_MakeBox(
        gp_Pnt(0.2, -2.0, -2.0),
        gp_Pnt(2.0, 2.0, 2.0)).Shape()
    ma = index_shape(sphere)
    mb = index_shape(cutter)

    face = ma.faces[0].face
    piece = SplitFacePiece(
        parent_face_id=ma.faces[0].face_id,
        piece_index=0,
        face=face,
        area=0.0,
        uv_witness=None,
        unchanged=True)
    fr = FaceSplitResult(
        parent_face_id=ma.faces[0].face_id,
        status="unchanged",
        pieces=[piece],
        source_edges=0,
        area_before=0.0,
        area_after=0.0,
        area_error=0.0)
    sp = ModelSplitResult(
        faces_a=[fr], faces_b=[],
        split_calls=0,
        affected_faces_a=0,
        affected_faces_b=0,
        unresolved_contacts=[])

    try:
        _classify_pieces(ma, mb, sp, "union", 1e-7)
    except AssemblyError as e:
        return check(
            "t6 straddling patch refuses multi-witness",
            getattr(e, "kind", "") == "PatchClassificationInconsistent",
            f"kind={getattr(e, 'kind', '?')} message={e}")
    return check(
        "t6 straddling patch refuses multi-witness",
        False, "classification unexpectedly accepted")


def _first_shell(shape):
    ex = TopExp_Explorer(shape, TopAbs_SHELL)
    assert ex.More()
    return TopoDS.Shell(ex.Current())


def t7_two_cavity_shell_nesting_is_consistent():
    """Multiple disjoint cavities must share the same outer parent.

    This directly exercises the shell-containment hierarchy rather than
    relying on one center point or a single-cavity special case.
    """
    outer = BRepPrimAPI_MakeSphere(
        gp_Pnt(0, 0, 0), 3.0).Shape()
    c1 = BRepPrimAPI_MakeSphere(
        gp_Pnt(-1.0, 0, 0), 0.5).Shape()
    c2 = BRepPrimAPI_MakeSphere(
        gp_Pnt(1.0, 0, 0), 0.5).Shape()

    records = _shell_records(
        [_first_shell(outer), _first_shell(c1), _first_shell(c2)],
        1e-7)
    solids = _build_nested_solids(records)

    depth0 = [r for r in records if r.depth == 0]
    depth1 = [r for r in records if r.depth == 1]
    ok = check(
        "t7 two cavities share one outer shell",
        len(depth0) == 1 and len(depth1) == 2
        and all(r.parent_shell == depth0[0].shell_index for r in depth1),
        f"records={[(r.shell_index,r.depth,r.parent_shell) for r in records]}")
    ok &= check(
        "t7 one solid with two cavities",
        len(solids) == 1 and len(solids[0].cavity_shells) == 2,
        f"solids={[(s.outer_shell,s.cavity_shells) for s in solids]}")
    want = 4.0 * math.pi / 3.0 * (3.0**3 - 2.0 * 0.5**3)
    ok &= check(
        "t7 two-cavity volume",
        abs(solids[0].volume - want) < 5e-6,
        f"got={solids[0].volume:.12g} expected={want:.12g}")
    return ok

def main():
    ok = True
    ok &= t1_overlap_spheres_all_ops()
    ok &= t2_disjoint_union_two_solids()
    ok &= t3_contained_difference_builds_cavity()
    ok &= t4_disjoint_intersection_is_empty()
    ok &= t5_exact_tangency_refuses_global_assembly()
    ok &= t6_unsplit_straddling_patch_refuses_multiwitness()
    ok &= t7_two_cavity_shell_nesting_is_consistent()
    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
