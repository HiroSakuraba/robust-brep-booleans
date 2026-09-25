"""Adversarial Tier B/C NURBS corpus.

New coverage beyond the baseline sphere/rounded-box cases:
1. periodic NURBS cylinder intersected across its seam-bearing side;
2. spherical-cap slivers down to 1e-5 feature height;
3. an oblique cut through a STEP-round-tripped, all-NURBS filleted solid;
4. a STEP-round-tripped periodic cylindrical face already trimmed by an
   earlier Boolean, then cut again on the opposite side;
5. a single NURBS torus/plane face pair with two disconnected intersection
   loops, pinning missed-branch detection.

The contract is intentionally asymmetric: straightforward transverse cases
must succeed and match an independent OCCT oracle; the smallest sliver may
instead produce a typed refusal, but may never silently return an inaccurate
valid-looking solid.
"""
import math
import os
import sys
import tempfile

sys.path.insert(0, "src")

from brepkernel import BRepAmbiguousResult, boolean_brep
from brepkernel.intersection import intersect_models
from brepkernel.step_ingest import index_shape

from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut
from OCP.BRepBuilderAPI import (
    BRepBuilderAPI_NurbsConvert,
    BRepBuilderAPI_Transform,
)
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import (
    BRepPrimAPI_MakeBox,
    BRepPrimAPI_MakeCylinder,
    BRepPrimAPI_MakeSphere,
    BRepPrimAPI_MakeTorus,
)
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import (
    STEPControl_AsIs,
    STEPControl_Reader,
    STEPControl_Writer,
)
from OCP.TopAbs import TopAbs_EDGE
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Ax1, gp_Dir, gp_Pnt, gp_Trsf


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def volume(shape):
    g = GProp_GProps()
    err = BRepGProp.VolumePropertiesGK_s(
        shape, g, 1e-10, True, True, False, False, False)
    assert float(err) >= 0.0
    return float(g.Mass())


def to_nurbs(shape):
    c = BRepBuilderAPI_NurbsConvert(shape, True)
    assert c.IsDone()
    out = c.Shape()
    assert BRepCheck_Analyzer(out, True).IsValid()
    return out


def cut_oracle(a, b):
    x = BRepAlgoAPI_Cut(a, b)
    x.Build()
    assert x.IsDone()
    out = x.Shape()
    assert BRepCheck_Analyzer(out, True).IsValid()
    return out


def common_oracle(a, b):
    x = BRepAlgoAPI_Common(a, b)
    x.Build()
    assert x.IsDone()
    out = x.Shape()
    assert BRepCheck_Analyzer(out, True).IsValid()
    return out


def all_edges(shape):
    out = []
    ex = TopExp_Explorer(shape, TopAbs_EDGE)
    while ex.More():
        e = TopoDS.Edge(ex.Current())
        if not any(e.IsSame(x) for x in out):
            out.append(e)
        ex.Next()
    return out


def step_roundtrip(shape):
    fd, path = tempfile.mkstemp(suffix=".step")
    os.close(fd)
    try:
        w = STEPControl_Writer()
        assert w.Transfer(shape, STEPControl_AsIs) == IFSelect_RetDone
        assert w.Write(path) == IFSelect_RetDone
        r = STEPControl_Reader()
        assert r.ReadFile(path) == IFSelect_RetDone
        assert r.TransferRoots() > 0
        out = r.OneShape()
        assert not out.IsNull()
        assert BRepCheck_Analyzer(out, True).IsValid()
        return out
    finally:
        try:
            os.remove(path)
        except OSError:
            pass


def rounded_nurbs_step():
    box = BRepPrimAPI_MakeBox(2.0, 1.5, 1.0).Shape()
    fillet = BRepFilletAPI_MakeFillet(box)
    edges = all_edges(box)
    assert len(edges) == 12
    for e in edges:
        fillet.Add(0.12, e)
    fillet.Build()
    assert fillet.IsDone()
    rounded = fillet.Shape()
    assert BRepCheck_Analyzer(rounded, True).IsValid()
    return step_roundtrip(to_nurbs(rounded))


def assert_oracle_close(name, out, oracle, report,
                        rel=4e-6, abs_tol=5e-10):
    rv = volume(out)
    ov = volume(oracle)
    err = abs(rv - ov)
    lim = max(abs_tol, rel * abs(ov))
    return check(
        name,
        report["accepted"]
        and BRepCheck_Analyzer(out, True).IsValid()
        and err <= lim,
        f"assembled={rv:.12g} oracle={ov:.12g} "
        f"err={err:.3e} lim={lim:.3e}")


def t1_periodic_nurbs_cylinder_transverse_cut():
    cyl = to_nurbs(BRepPrimAPI_MakeCylinder(1.0, 2.0).Shape())
    model = index_shape(cyl)

    periodic = []
    for fr in model.faces:
        if "BSpline" not in fr.surface_type:
            continue
        bs = BRepAdaptor_Surface(fr.face).BSpline()
        if bs.IsUPeriodic() or bs.IsVPeriodic():
            periodic.append(fr.face_id)

    ok = check(
        "a1 periodic NURBS support preserved",
        bool(periodic) and len(model.nurbs_faces) == len(model.faces),
        f"periodic_faces={periodic} "
        f"accels={len(model.nurbs_faces)}/{len(model.faces)}")

    cutter = BRepPrimAPI_MakeBox(
        gp_Pnt(0.25, -1.5, -0.5),
        gp_Pnt(1.5, 1.5, 2.5)).Shape()
    out, report = boolean_brep(cyl, cutter, "difference")
    oracle = cut_oracle(cyl, cutter)

    ix = report["stages"]["intersection"]
    asm = report["stages"]["assembly"]
    ok &= check(
        "a1 periodic cut exercises verified sections",
        ix["candidate_face_pairs"] >= 1
        and ix["verified_edges"] >= 1
        and ix["ambiguous_contacts"] == 0,
        f"intersection={ix}")
    ok &= check(
        "a1 periodic result closed/manifold",
        asm["free_edges"] == 0
        and asm["multiple_edges"] == 0
        and asm["edge_lineage"]["unattributed_edges"] == 0,
        f"assembly={{'free':{asm['free_edges']},"
        f"'multiple':{asm['multiple_edges']},"
        f"'unattributed':{asm['edge_lineage']['unattributed_edges']}}}")
    ok &= assert_oracle_close("a1 periodic cut oracle", out, oracle, report)
    return ok


def _sliver_case(height):
    sphere = to_nurbs(
        BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), 1.0).Shape())
    z = 1.0 - float(height)
    cutter = BRepPrimAPI_MakeBox(
        gp_Pnt(-2.0, -2.0, -2.0),
        gp_Pnt(2.0, 2.0, z)).Shape()
    oracle = cut_oracle(sphere, cutter)
    ov = volume(oracle)

    try:
        out, report = boolean_brep(
            sphere, cutter, "difference",
            base_tol=1e-7,
            contact_tol=4e-7)
    except BRepAmbiguousResult as exc:
        return False, None, ov, exc.report

    rv = volume(out)
    err = abs(rv - ov)
    lim = max(5e-11, 8e-4 * abs(ov))
    good = (
        report["accepted"]
        and BRepCheck_Analyzer(out, True).IsValid()
        and err <= lim)
    return good, (rv, err, lim), ov, report


def t2_sliver_scale_sweep():
    ok = True
    for h in (1e-3, 1e-4):
        good, info, ov, report = _sliver_case(h)
        ok &= check(
            f"a2 cap h={h:g} accepted accurately",
            good,
            (f"result={info} oracle={ov:.12g} "
             f"refusal={report.get('refusal') if report else None}"))

    # At 1e-5 we permit a deliberate refusal because the feature is only
    # 100x base_tol high. Acceptance still has to agree with the oracle.
    good, info, ov, report = _sliver_case(1e-5)
    if good:
        ok &= check(
            "a2 cap h=1e-5 accepted accurately",
            True,
            f"result={info} oracle={ov:.12g}")
    else:
        refusal = report.get("refusal", {}) if report else {}
        safe = {
            "UnresolvedContact",
            "SectionToleranceTooLoose",
            "SectionGeometryMismatch",
            "PatchClassificationInconsistent",
            "InsufficientPatchWitnesses",
            "InsufficientShellWitnesses",
            "OpenAssembly",
            "SewingInvalid",
            "SolidInvalid",
            "FinalBRepInvalid",
        }
        ok &= check(
            "a2 cap h=1e-5 typed refusal",
            refusal.get("kind") in safe,
            f"oracle={ov:.12g} refusal={refusal}")
    return ok


def t3_oblique_imported_blend_cut():
    source = rounded_nurbs_step()
    model = index_shape(source)
    ok = check(
        "a3 imported blend corpus is all freeform",
        len(model.faces) >= 20
        and len(model.nurbs_faces) == len(model.faces),
        f"faces={len(model.faces)} nurbs={len(model.nurbs_faces)}")

    base = BRepPrimAPI_MakeBox(
        gp_Pnt(0.92, -0.6, -0.5),
        gp_Pnt(2.8, 2.1, 1.5)).Shape()
    tr = gp_Trsf()
    tr.SetRotation(
        gp_Ax1(gp_Pnt(1.0, 0.75, 0.5), gp_Dir(0, 0, 1)),
        math.radians(11.0))
    cutter = BRepBuilderAPI_Transform(base, tr, True, False).Shape()
    assert BRepCheck_Analyzer(cutter, True).IsValid()

    out, report = boolean_brep(
        source, cutter, "difference",
        base_tol=1e-7)
    oracle = cut_oracle(source, cutter)

    ix = report["stages"]["intersection"]
    sp = report["stages"]["split"]
    asm = report["stages"]["assembly"]
    ok &= check(
        "a3 oblique cut has broad multi-face workset",
        ix["candidate_face_pairs"] >= 6
        and ix["verified_edges"] >= 6
        and ix["ambiguous_contacts"] == 0
        and sp["affected_faces_A"] >= 6,
        f"intersection={ix} split={sp}")
    ok &= check(
        "a3 oblique result provenance complete",
        asm["edge_lineage"]["boolean_section_edges"] >= 6
        and asm["edge_lineage"]["unattributed_edges"] == 0
        and asm["free_edges"] == 0
        and asm["multiple_edges"] == 0,
        f"lineage={asm['edge_lineage']}")
    ok &= assert_oracle_close(
        "a3 oblique blend oracle", out, oracle, report,
        rel=5e-6, abs_tol=1e-9)
    return ok



def t4_pretrimmed_periodic_step_roundtrip():
    """Periodic support + inherited trim/seam state from a prior Boolean.

    Build an analytic cylinder, carve a local side notch, convert the already
    trimmed result to NURBS, round-trip through STEP, then perform a second
    transverse cut on the opposite side. This exercises parameter-frame
    recovery after prior topology work rather than only a pristine periodic
    primitive.
    """
    analytic = BRepPrimAPI_MakeCylinder(1.0, 2.0).Shape()
    notch = BRepPrimAPI_MakeBox(
        gp_Pnt(0.55, -0.28, 0.35),
        gp_Pnt(1.45, 0.28, 1.65)).Shape()
    first = BRepAlgoAPI_Cut(analytic, notch)
    first.Build()
    assert first.IsDone()
    trimmed = first.Shape()
    assert BRepCheck_Analyzer(trimmed, True).IsValid()

    source = step_roundtrip(to_nurbs(trimmed))
    model = index_shape(source)
    periodic = []
    for fr in model.faces:
        if "BSpline" not in fr.surface_type:
            continue
        bs = BRepAdaptor_Surface(fr.face).BSpline()
        if bs.IsUPeriodic() or bs.IsVPeriodic():
            periodic.append(fr.face_id)

    ok = check(
        "a4 pretrimmed STEP retains periodic freeform support",
        bool(periodic)
        and len(model.nurbs_faces) == len(model.faces)
        and len(model.faces) >= 7,
        f"faces={len(model.faces)} periodic={periodic} "
        f"accels={len(model.nurbs_faces)}")

    second = BRepPrimAPI_MakeBox(
        gp_Pnt(-1.45, -1.2, 0.25),
        gp_Pnt(-0.20, 1.2, 1.75)).Shape()
    out, report = boolean_brep(source, second, "difference")
    oracle = cut_oracle(source, second)

    ix = report["stages"]["intersection"]
    asm = report["stages"]["assembly"]
    ok &= check(
        "a4 inherited trim produces verified sections",
        ix["candidate_face_pairs"] >= 3
        and ix["verified_edges"] >= 3
        and ix["ambiguous_contacts"] == 0,
        f"intersection={ix}")
    ok &= check(
        "a4 inherited trim provenance complete",
        asm["free_edges"] == 0
        and asm["multiple_edges"] == 0
        and asm["edge_lineage"]["unattributed_edges"] == 0
        and asm["edge_lineage"]["boolean_section_edges"] >= 3,
        f"assembly={{'free':{asm['free_edges']},"
        f"'multiple':{asm['multiple_edges']},"
        f"'section_edges':{asm['edge_lineage']['boolean_section_edges']},"
        f"'unattributed':{asm['edge_lineage']['unattributed_edges']}}}")
    ok &= assert_oracle_close(
        "a4 pretrimmed periodic STEP oracle",
        out, oracle, report,
        rel=5e-6, abs_tol=1e-9)
    return ok


def t5_two_loop_torus_periodic_seam_boolean():
    """Two disconnected loops, one periodic seam, must form a real Boolean.

    At z=0 a torus with major radius 3 and minor radius 1 meets the cutter's
    top plane in two circles. One loop is already the torus's periodic seam,
    so it is reused as existing topology on the torus side and remains a split
    tool on the cutter side. The second loop splits both operands normally.
    """
    torus = to_nurbs(BRepPrimAPI_MakeTorus(3.0, 1.0).Shape())
    cutter = BRepPrimAPI_MakeBox(
        gp_Pnt(-5.0, -5.0, -2.0),
        gp_Pnt(5.0, 5.0, 0.0)).Shape()

    ma = index_shape(torus)
    mb = index_shape(cutter)
    ixr = intersect_models(ma, mb, base_tol=1e-7)
    curve_pairs = [p for p in ixr.pairs if p.edges]

    ok = check(
        "a5 one face pair keeps two disconnected loops",
        ixr.candidate_pairs == 1
        and ixr.section_calls == 1
        and ixr.verified_edges == 2
        and len(curve_pairs) == 1
        and len(curve_pairs[0].edges) == 2
        and ixr.completeness_probes == 1
        and ixr.raw_curve_count >= 2
        and ixr.raw_trimmed_components >= 2
        and ixr.raw_unmatched_components == 0
        and ixr.ambiguous_contacts == 0,
        f"candidate_pairs={ixr.candidate_pairs} "
        f"verified={ixr.verified_edges} raw={ixr.raw_curve_count} "
        f"trimmed={ixr.raw_trimmed_components} "
        f"unmatched={ixr.raw_unmatched_components} "
        f"max_d={ixr.completeness_max_distance:.3e}")

    risks = [e.risk_flags for e in curve_pairs[0].edges]
    ok &= check(
        "a5 exactly one loop reuses torus seam",
        sum("seam_on_a" in r for r in risks) == 1,
        f"risks={risks}")

    out, report = boolean_brep(torus, cutter, "difference")
    oracle = cut_oracle(torus, cutter)
    sp = report["stages"]["split"]
    asm = report["stages"]["assembly"]
    ver = report["stages"]["verification"]

    ok &= check(
        "a5 seam-aware split completes",
        not sp["unresolved_contacts"]
        and sp["affected_faces_A"] == 1
        and sp["affected_faces_B"] == 1
        and sp["reused_seam_edges_A"] == 1
        and sp["reused_seam_edges_B"] == 0
        and sp["shared_seam_refusals"] == 0,
        f"split={sp}")
    ok &= check(
        "a5 seam-aware result closed and auditable",
        asm["free_edges"] == 0
        and asm["multiple_edges"] == 0
        and asm["edge_lineage"]["boolean_section_edges"] == 2
        and asm["edge_lineage"]["unattributed_edges"] == 0
        and ver["brep_valid"]
        and ver["closed"]
        and ver["manifold_edges"]
        and ver["complete_edge_lineage"]
        and ver["volume_bounds_ok"],
        f"assembly={asm['edge_lineage']} verification={ver}")
    ok &= assert_oracle_close(
        "a5 periodic-seam torus oracle",
        out, oracle, report,
        rel=5e-6, abs_tol=1e-9)
    return ok


def t6_reversed_operand_torus_intersection():
    """Prove seam routing through assembly when the torus is operand B."""
    torus = to_nurbs(BRepPrimAPI_MakeTorus(3.0, 1.0).Shape())
    cutter = BRepPrimAPI_MakeBox(
        gp_Pnt(-5.0, -5.0, -2.0),
        gp_Pnt(5.0, 5.0, 0.0)).Shape()

    out, report = boolean_brep(cutter, torus, "intersection")
    oracle = common_oracle(cutter, torus)

    ix = report["stages"]["intersection"]
    sp = report["stages"]["split"]
    asm = report["stages"]["assembly"]
    ver = report["stages"]["verification"]

    ok = check(
        "a6 reversed operand seam intersection workset",
        ix["candidate_face_pairs"] == 1
        and ix["verified_edges"] == 2
        and ix["raw_trimmed_components"] >= 2
        and ix["raw_unmatched_components"] == 0,
        f"intersection={ix}")
    ok &= check(
        "a6 reversed operand split/assembly complete",
        not sp["unresolved_contacts"]
        and sp["affected_faces_A"] == 1
        and sp["affected_faces_B"] == 1
        and sp["reused_seam_edges_A"] == 0
        and sp["reused_seam_edges_B"] == 1
        and sp["shared_seam_refusals"] == 0
        and asm["free_edges"] == 0
        and asm["multiple_edges"] == 0
        and asm["edge_lineage"]["boolean_section_edges"] == 2
        and asm["edge_lineage"]["unattributed_edges"] == 0
        and ver["brep_valid"]
        and ver["closed"]
        and ver["manifold_edges"]
        and ver["complete_edge_lineage"]
        and ver["volume_bounds_ok"],
        f"split={sp} lineage={asm['edge_lineage']} verification={ver}")
    ok &= assert_oracle_close(
        "a6 reversed operand torus oracle",
        out, oracle, report,
        rel=5e-6, abs_tol=1e-9)
    return ok

def main():
    ok = True
    ok &= t1_periodic_nurbs_cylinder_transverse_cut()
    ok &= t2_sliver_scale_sweep()
    ok &= t3_oblique_imported_blend_cut()
    ok &= t4_pretrimmed_periodic_step_roundtrip()
    ok &= t5_two_loop_torus_periodic_seam_boolean()
    ok &= t6_reversed_operand_torus_intersection()
    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
