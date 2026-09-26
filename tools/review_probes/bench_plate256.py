#!/usr/bin/env python3
"""G12 benchmark: 262-face plate (256 holes) minus slot, best-of-3.

Rebuilds the holed plate with OCCT (cached to a BREP file), then times
brepkernel boolean_brep(drilled, slot, "difference").
"""
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                           "..", "..", "src"))

CACHE = os.path.join(os.path.dirname(__file__), "..", "..",
                     "goals", "robust-b-rep-booleans-prototype",
                     "hidden_files", "plate256.brep")
CACHE = os.path.abspath(CACHE)


def build_drilled():
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Compound
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut

    plate = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 8, 8, 0.5).Shape()
    comp = TopoDS_Compound()
    b = BRep_Builder()
    b.MakeCompound(comp)
    for i in range(16):
        for j in range(16):
            x = 0.25 + i * 0.5
            y = 0.25 + j * 0.5
            cyl = BRepPrimAPI_MakeCylinder(
                gp_Ax2(gp_Pnt(x, y, -0.1), gp_Dir(0, 0, 1)),
                0.15, 0.7).Shape()
            b.Add(comp, cyl)
    cut = BRepAlgoAPI_Cut(plate, comp)
    cut.Build()
    if not cut.IsDone():
        raise RuntimeError("OCCT plate drilling failed")
    return cut.Shape()


def count_faces(shape):
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopAbs import TopAbs_FACE
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    n = 0
    while ex.More():
        n += 1
        ex.Next()
    return n


def main():
    from OCP.TopoDS import TopoDS_Shape
    from OCP.BRepTools import BRepTools
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Copy
    from OCP.BRep import BRep_Builder
    if os.path.exists(CACHE):
        shp = TopoDS_Shape()
        BRepTools.Read_s(shp, CACHE, BRep_Builder())
        drilled = shp
        print(f"loaded cached plate from {CACHE}", flush=True)
    else:
        print("drilling plate with OCCT ...", flush=True)
        t0 = time.perf_counter()
        drilled = build_drilled()
        print(f"drilled in {time.perf_counter() - t0:.1f}s", flush=True)
        os.makedirs(os.path.dirname(CACHE), exist_ok=True)
        from OCP.BRepTools import BRepTools as _BT
        _BT.Write_s(drilled, CACHE)
        print(f"cached to {CACHE}", flush=True)
    print("plate faces:", count_faces(drilled), flush=True)

    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt
    slot = BRepPrimAPI_MakeBox(gp_Pnt(3.0, 3.5, -0.25), 2.0, 1.0, 1.0).Shape()

    from brepkernel.pipeline import boolean_brep
    best = None
    for k in range(3):
        t0 = time.perf_counter()
        _shape, report = boolean_brep(drilled, slot, "difference")
        dt = time.perf_counter() - t0
        rs = report["stages"]["assembly"].get("region_stats")
        vol = report["stages"]["assembly"].get("volume")
        print(f"run {k}: {dt:.2f}s volume={vol:.4f} "
              f"region_stats={rs}", flush=True)
        best = dt if best is None else min(best, dt)
    print(f"BEST: {best:.2f}s", flush=True)


if __name__ == "__main__":
    main()
