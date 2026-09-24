"""Public Tier B/C pipeline regressions.

These tests exercise the one-call API rather than manually chaining the
internal stages.
"""
import math
import sys

sys.path.insert(0, "src")

from brepkernel import boolean_brep, BRepAmbiguousResult

from OCP.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeSphere
from OCP.GProp import GProp_GProps
from OCP.TopAbs import TopAbs_SOLID
from OCP.TopExp import TopExp_Explorer
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


def solid_count(shape):
    n = 0
    ex = TopExp_Explorer(shape, TopAbs_SOLID)
    while ex.More():
        n += 1
        ex.Next()
    return n


def t1_one_call_true_nurbs_union():
    a0 = BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), 1.0).Shape()
    b0 = BRepPrimAPI_MakeSphere(gp_Pnt(1, 0, 0), 1.0).Shape()
    a = BRepBuilderAPI_NurbsConvert(a0, True).Shape()
    b = BRepBuilderAPI_NurbsConvert(b0, True).Shape()

    out, report = boolean_brep(a, b, "union")
    ing = report["stages"]["ingest"]
    ix = report["stages"]["intersection"]
    asm = report["stages"]["assembly"]
    ok = check("p1 accepted", report["accepted"])
    ok &= check(
        "p1 periodic accelerators visible in report",
        ing["A"]["freeform_accels"] == 1
        and ing["B"]["freeform_accels"] == 1
        and ing["A"]["local_patches"] >= 8
        and ing["B"]["local_patches"] >= 8,
        f"ingest={ing}")
    ok &= check(
        "p1 conservative exact workset",
        ix["candidate_face_pairs"] == 1 and ix["section_calls"] == 1
        and ix["verified_edges"] >= 1
        and ix["ambiguous_contacts"] == 0,
        f"intersection={ix}")
    want = 9.0 * math.pi / 4.0
    ok &= check(
        "p1 result volume",
        abs(volume(out) - want) < 3e-6
        and abs(asm["volume"] - want) < 3e-6,
        f"shape={volume(out):.12g} report={asm['volume']:.12g}")
    ok &= check(
        "p1 final verification",
        report["stages"]["verification"]["brep_valid"]
        and report["stages"]["verification"]["closed"]
        and report["stages"]["verification"]["manifold_edges"])
    return ok


def t2_exact_identity_fast_path():
    a = BRepPrimAPI_MakeSphere(1.0).Shape()
    u, ru = boolean_brep(a, a, "union")
    d, rd = boolean_brep(a, a, "difference")
    ok = check(
        "p2 union exact identity",
        u.IsSame(a) and ru["accepted"]
        and ru["stages"]["identity"]["resolution"] == "A"
        and "intersection" not in ru["stages"])
    ok &= check(
        "p2 difference exact empty",
        rd["accepted"]
        and rd["stages"]["identity"]["resolution"] == "empty"
        and solid_count(d) == 0)
    return ok


def t3_tangent_contact_structured_refusal():
    a = BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), 1.0).Shape()
    b = BRepPrimAPI_MakeSphere(gp_Pnt(2, 0, 0), 1.0).Shape()
    try:
        boolean_brep(a, b, "union")
    except BRepAmbiguousResult as e:
        refusal = e.report.get("refusal", {})
        split = e.report["stages"].get("split", {})
        return check(
            "p3 tangent refuses with report",
            not e.report["accepted"]
            and refusal.get("stage") == "assembly"
            and refusal.get("kind") == "UnresolvedContact"
            and split.get("unresolved_contacts"),
            f"refusal={refusal} split={split}")
    return check("p3 tangent refuses with report", False, "no refusal")



def t4_independent_same_domain_fast_path():
    a = BRepPrimAPI_MakeBox(1.0, 2.0, 3.0).Shape()
    b = BRepPrimAPI_MakeBox(1.0, 2.0, 3.0).Shape()
    assert not a.IsSame(b)
    u, ru = boolean_brep(a, b, "union")
    d, rd = boolean_brep(a, b, "difference")
    sd_u = ru["stages"].get("same_domain", {})
    sd_d = rd["stages"].get("same_domain", {})
    ok = check(
        "p4 independent same-domain union",
        ru["accepted"] and sd_u.get("equivalent")
        and sd_u.get("matched_faces") == 6
        and sd_u.get("resolution") == "A"
        and "intersection" not in ru["stages"]
        and solid_count(u) == 1,
        f"same_domain={sd_u}")
    ok &= check(
        "p4 independent same-domain difference",
        rd["accepted"] and sd_d.get("equivalent")
        and sd_d.get("resolution") == "empty"
        and "intersection" not in rd["stages"]
        and solid_count(d) == 0,
        f"same_domain={sd_d}")
    return ok

def main():
    ok = True
    ok &= t1_one_call_true_nurbs_union()
    ok &= t2_exact_identity_fast_path()
    ok &= t3_tangent_contact_structured_refusal()
    ok &= t4_independent_same_domain_fast_path()
    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
