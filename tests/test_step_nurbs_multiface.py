"""Imported multi-face NURBS STEP Boolean regression.

Build a rounded solid, convert its surfaces to NURBS, round-trip it through a
real STEP file, then run the Tier B/C Boolean path against a transverse cutter.
The goal is to exercise import topology, multiple trimmed freeform faces,
intersection, local splitting, classification, assembly and provenance in one
case that is materially harder than two one-face spheres.
"""
import os
import sys
import tempfile

sys.path.insert(0, "src")

from brepkernel import boolean_brep
from brepkernel.step_ingest import index_shape

from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCP.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_AsIs, STEPControl_Reader, STEPControl_Writer
from OCP.TopAbs import TopAbs_EDGE
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Pnt


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def volume(shape):
    g = GProp_GProps()
    err = BRepGProp.VolumePropertiesGK_s(
        shape, g, 1e-9, True, True, False, False, False)
    assert float(err) >= 0.0
    return float(g.Mass())


def all_edges(shape):
    out = []
    ee = TopExp_Explorer(shape, TopAbs_EDGE)
    while ee.More():
        e = TopoDS.Edge(ee.Current())
        if not any(e.IsSame(x) for x in out):
            out.append(e)
        ee.Next()
    return out


def rounded_nurbs_shape():
    box = BRepPrimAPI_MakeBox(2.0, 1.5, 1.0).Shape()
    fillet = BRepFilletAPI_MakeFillet(box)
    es = all_edges(box)
    assert len(es) == 12
    for e in es:
        fillet.Add(0.12, e)
    fillet.Build()
    assert fillet.IsDone()
    rounded = fillet.Shape()
    assert BRepCheck_Analyzer(rounded, True).IsValid()

    conv = BRepBuilderAPI_NurbsConvert(rounded, True)
    assert conv.IsDone()
    nurbs = conv.Shape()
    assert BRepCheck_Analyzer(nurbs, True).IsValid()
    return nurbs


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


def main():
    source = step_roundtrip(rounded_nurbs_shape())
    model = index_shape(source)

    ok = check(
        "step1 imported multi-face solid",
        len(model.solids) == 1 and len(model.faces) >= 8,
        f"solids={len(model.solids)} faces={len(model.faces)}")
    ok &= check(
        "step1 imported freeform coverage",
        len(model.nurbs_faces) >= 4,
        f"nurbs={len(model.nurbs_faces)}/{len(model.faces)} "
        f"types={[f.surface_type for f in model.faces]}")

    cutter = BRepPrimAPI_MakeBox(
        gp_Pnt(1.0, -0.4, -0.4),
        gp_Pnt(2.5, 1.9, 1.4)).Shape()

    out, report = boolean_brep(source, cutter, "difference")
    ix = report["stages"]["intersection"]
    sp = report["stages"]["split"]
    asm = report["stages"]["assembly"]

    ok &= check(
        "step2 multi-face workset",
        ix["candidate_face_pairs"] >= 4
        and ix["section_calls"] >= 4
        and ix["verified_edges"] >= 4
        and ix["ambiguous_contacts"] == 0,
        f"intersection={ix}")
    ok &= check(
        "step2 local splitting",
        sp["affected_faces_A"] >= 4
        and sp["split_calls"] >= 4
        and not sp["unresolved_contacts"],
        f"split={sp}")
    ok &= check(
        "step2 assembled valid solid",
        report["accepted"]
        and asm["solids"] == 1
        and asm["free_edges"] == 0
        and asm["multiple_edges"] == 0
        and BRepCheck_Analyzer(out, True).IsValid(),
        f"assembly={{'faces':{asm['selected_faces']},'shells':{asm['shells']},"
        f"'solids':{asm['solids']},'free':{asm['free_edges']},"
        f"'multiple':{asm['multiple_edges']}}}")

    oracle = BRepAlgoAPI_Cut(source, cutter)
    oracle.Build()
    assert oracle.IsDone()
    ov = volume(oracle.Shape())
    rv = volume(out)
    ok &= check(
        "step2 volume agrees OCCT oracle",
        abs(rv - ov) <= max(3e-6, 3e-6 * abs(ov)),
        f"assembled={rv:.12g} oracle={ov:.12g} err={abs(rv-ov):.3e}")

    lin = asm["edge_lineage"]
    payloads = asm["section_payloads"]
    ok &= check(
        "step2 imported result provenance",
        lin["boolean_section_edges"] >= 4
        and lin["unattributed_edges"] == 0
        and len(payloads) >= 4
        and all(p["samples"] >= 2 for p in payloads)
        and any(p["result_edges"] for p in payloads),
        f"lineage={lin} payloads={payloads}")

    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
