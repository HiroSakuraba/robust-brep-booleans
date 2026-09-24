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

from brepkernel.assembly import assemble_boolean, AssemblyError
from brepkernel.intersection import intersect_models
from brepkernel.split import split_models
from brepkernel.step_ingest import index_shape

from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import BRepPrimAPI_MakeSphere
from OCP.GProp import GProp_GProps
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
            len(r.decisions) == 4
            and all(d.classification in ("inside", "outside")
                    for d in r.decisions),
            f"decisions={[(d.operand,d.classification,d.keep) for d in r.decisions]}")
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


def main():
    ok = True
    ok &= t1_overlap_spheres_all_ops()
    ok &= t2_disjoint_union_two_solids()
    ok &= t3_contained_difference_builds_cavity()
    ok &= t4_disjoint_intersection_is_empty()
    ok &= t5_exact_tangency_refuses_global_assembly()
    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
