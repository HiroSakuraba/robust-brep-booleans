"""True NURBS end-to-end Boolean regression.

Analytic spheres are converted by OCCT into B-spline/NURBS B-reps before they
enter this kernel.  The result must still traverse the same verified pipeline
and reproduce the analytic/OCCT volume without falling back to a mesh-defined
answer.
"""
import math
import sys

sys.path.insert(0, "src")

from brepkernel import BRepAmbiguousResult, boolean_brep
from brepkernel.assembly import assemble_boolean
from brepkernel.intersection import intersect_models
from brepkernel.split import split_models
from brepkernel.step_ingest import index_shape

from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCP.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeSphere
from OCP.GProp import GProp_GProps
from OCP.gp import gp_Pnt


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def volume(shape):
    p = GProp_GProps()
    err = BRepGProp.VolumePropertiesGK_s(
        shape, p, 1e-10, True, True, False, False, False)
    assert float(err) >= 0.0
    return float(p.Mass())


def to_nurbs(shape):
    c = BRepBuilderAPI_NurbsConvert(shape, True)
    assert c.IsDone()
    return c.Shape()



def t2_tiny_nurbs_cap_accepts_accurately_or_refuses():
    """A 1e-4-high NURBS spherical cap must never become a silent wrong solid.

    Acceptance requires a valid B-rep and close agreement with an independent
    OCCT cut oracle. Refusal is also correct for this deliberately difficult
    near-tangent/tiny-feature regime, but it must be a typed Tier B/C refusal.
    """
    sphere = to_nurbs(
        BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), 1.0).Shape())
    cutter = BRepPrimAPI_MakeBox(
        gp_Pnt(-2.0, -2.0, -2.0),
        gp_Pnt(2.0, 2.0, 0.9999)).Shape()

    oracle = BRepAlgoAPI_Cut(sphere, cutter)
    oracle.Build()
    assert oracle.IsDone()
    ov = volume(oracle.Shape())
    assert ov > 0.0

    try:
        out, report = boolean_brep(
            sphere, cutter, "difference",
            base_tol=1e-7, chord_tol=2e-7,
            contact_tol=4e-7)
    except BRepAmbiguousResult as exc:
        refusal = exc.report.get("refusal", {})
        safe_kinds = {
            "UnresolvedContact",
            "SectionToleranceTooLoose",
            "PatchClassificationInconsistent",
            "InsufficientPatchWitnesses",
            "InsufficientShellWitnesses",
            "OpenAssembly",
            "SewingInvalid",
            "SolidInvalid",
        }
        return check(
            "n2 tiny cap typed refusal",
            refusal.get("kind") in safe_kinds
            and not exc.report.get("accepted", False),
            f"refusal={refusal}")

    rv = volume(out)
    abs_err = abs(rv - ov)
    rel_err = abs_err / max(abs(ov), 1e-30)
    return check(
        "n2 tiny cap accepted accurately",
        report["accepted"]
        and BRepCheck_Analyzer(out, True).IsValid()
        and abs_err <= max(5e-10, 5e-4 * abs(ov)),
        f"assembled={rv:.12g} oracle={ov:.12g} "
        f"abs_err={abs_err:.3e} rel_err={rel_err:.3e}")

def main():
    analytic_a = BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), 1.0).Shape()
    analytic_b = BRepPrimAPI_MakeSphere(gp_Pnt(1, 0, 0), 1.0).Shape()
    a = to_nurbs(analytic_a)
    b = to_nurbs(analytic_b)

    ma = index_shape(a)
    mb = index_shape(b)

    ok = check(
        "n1 inputs are B-spline surfaces",
        ma.faces and mb.faces
        and all("BSpline" in f.surface_type for f in ma.faces)
        and all("BSpline" in f.surface_type for f in mb.faces),
        f"A={[f.surface_type for f in ma.faces]} "
        f"B={[f.surface_type for f in mb.faces]}")

    ok &= check(
        "n1 periodic NURBS accelerators built",
        len(ma.nurbs_faces) == len(ma.faces)
        and len(mb.nurbs_faces) == len(mb.faces)
        and all(len(fr.freeform.index.records) > 1
                for fr in ma.faces + mb.faces),
        f"A={len(ma.nurbs_faces)}/{len(ma.faces)} "
        f"B={len(mb.nurbs_faces)}/{len(mb.faces)} "
        f"patchesA={[len(fr.freeform.index.records) if fr.freeform else 0 for fr in ma.faces]} "
        f"patchesB={[len(fr.freeform.index.records) if fr.freeform else 0 for fr in mb.faces]}")

    ix = intersect_models(
        ma, mb, base_tol=1e-7, chord_tol=1e-5,
        tangent_sin_tol=1e-4)
    sp = split_models(ma, mb, ix, base_tol=1e-7)
    r = assemble_boolean(ma, mb, sp, "union", base_tol=1e-7)

    expected = 9.0 * math.pi / 4.0
    oracle = BRepAlgoAPI_Fuse(a, b)
    oracle.Build()
    assert oracle.IsDone()
    ov = volume(oracle.Shape())

    ok &= check(
        "n1 verified section exists",
        ix.verified_edges >= 1 and not ix.has_ambiguous_contact,
        f"edges={ix.verified_edges} ambiguous={ix.ambiguous_contacts}")
    ok &= check(
        "n1 global assembly valid",
        r.free_edges == 0 and r.multiple_edges == 0
        and BRepCheck_Analyzer(r.shape, True).IsValid(),
        f"selected={r.selected_faces} shells={len(r.shells)} "
        f"solids={len(r.solids)}")
    ok &= check(
        "n1 analytic volume",
        abs(r.volume - expected) < 3e-6,
        f"assembled={r.volume:.12g} expected={expected:.12g}")
    ok &= check(
        "n1 OCCT NURBS oracle volume",
        abs(r.volume - ov) < 3e-6,
        f"assembled={r.volume:.12g} oracle={ov:.12g}")

    # The result must keep explicit input-face lineage through patch selection.
    kept_a = [d for d in r.decisions if d.operand == "A" and d.keep]
    kept_b = [d for d in r.decisions if d.operand == "B" and d.keep]
    ok &= check(
        "n1 bilateral provenance",
        bool(kept_a) and bool(kept_b)
        and all(d.sewed_face is not None for d in kept_a + kept_b),
        f"keptA={len(kept_a)} keptB={len(kept_b)}")

    ok &= t2_tiny_nurbs_cap_accepts_accurately_or_refuses()
    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
