#!/usr/bin/env python3
"""S8B: assembly sub-stage profile for the two large benchmarks.

Extends the S8 decision rule one level down: splits the assembly stage
into classification vs topology-commit, and classification further into
witness generation / OCCT classifier / ray classifier / boundary
distance / region build / propagation / coincidence, so the S8C
parallelism question can be decided on measured shares.

Usage: profile_s8b.py [--many-n 8]
"""
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from bench_plate256 import build_drilled, CACHE, count_faces  # noqa

CLASSIFY_KEYS = [
    "classify.region_build",
    "classify.witness_gen",
    "classify.occt_classifier",
    "classify.ray_classifier",
    "classify.boundary_distance",
    "classify.coincident",
    "classify.propagation",
]
COMMIT_KEYS = [
    "commit.face_select",
    "commit.sewing",
    "commit.shell_extract",
    "commit.shell_records",
    "commit.solid_build",
    "commit.topo_checks",
    "commit.volume",
    "commit.lineage",
]


def show(title, report, wall_s):
    tm = report.get("timings_ms", {})
    asm_ms = float(tm.get("assembly", 0.0))
    sub = (report.get("stages", {}).get("assembly", {})
           .get("substage_ms", {}))
    print(f"=== {title} ===")
    print(f"wall {wall_s:.2f}s  assembly {asm_ms/1000:.3f}s "
          f"({100*asm_ms/1000/wall_s:.1f}% of wall)")
    if not sub:
        print("  (no substage_ms recorded)")
        return
    classify_ms = sub.get("classify", 0.0)
    commit_ms = sum(sub.get(k, 0.0) for k in COMMIT_KEYS)
    leaves_ms = sum(sub.get(k, 0.0) for k in CLASSIFY_KEYS)
    print(f"  CLASSIFY {classify_ms/1000:.3f}s "
          f"({100*classify_ms/1000/wall_s:.1f}% of wall, "
          f"{100*classify_ms/max(asm_ms,1e-9):.1f}% of assembly; "
          f"leaves sum {leaves_ms/1000:.3f}s)")
    for k in CLASSIFY_KEYS:
        v = sub.get(k, 0.0)
        if v > 0:
            print(f"    {k:32s} {v/1000:8.3f}s "
                  f"({100*v/1000/wall_s:5.1f}% wall)")
    print(f"  COMMIT   {commit_ms/1000:.3f}s "
          f"({100*commit_ms/1000/wall_s:.1f}% of wall, "
          f"{100*commit_ms/max(asm_ms,1e-9):.1f}% of assembly)")
    for k in COMMIT_KEYS:
        v = sub.get(k, 0.0)
        if v > 0:
            print(f"    {k:32s} {v/1000:8.3f}s "
                  f"({100*v/1000/wall_s:5.1f}% wall)")
    accounted = classify_ms + commit_ms
    print(f"  accounted {accounted/1000:.3f}s vs assembly "
          f"{asm_ms/1000:.3f}s "
          f"({100*accounted/max(asm_ms,1e-9):.1f}%)")
    # Decision-relevant rollups
    indep = (sub.get("classify.occt_classifier", 0.0)
             + sub.get("classify.ray_classifier", 0.0)
             + sub.get("classify.boundary_distance", 0.0))
    print(f"  point-classification core (occt+ray+bdistance): "
          f"{indep/1000:.3f}s ({100*indep/1000/wall_s:.1f}% of wall)")
    perf = report.get("performance", {})
    counts = perf.get("counters", {}) if isinstance(perf, dict) else {}
    interesting = ("ray_intersector_perform", "ray_aabb_prune",
                   "solid_classifier_eval", "occt_box_prune",
                   "ray_pair_cast")
    shown = {k: counts[k] for k in interesting if k in counts}
    if shown:
        print(f"  counters: " + ", ".join(
            f"{k}={v}" for k, v in shown.items()))


def plate():
    from OCP.TopoDS import TopoDS_Shape
    from OCP.BRepTools import BRepTools
    from OCP.BRep import BRep_Builder
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt
    from brepkernel.pipeline import boolean_brep

    if os.path.exists(CACHE):
        shp = TopoDS_Shape()
        BRepTools.Read_s(shp, CACHE, BRep_Builder())
        drilled = shp
    else:
        drilled = build_drilled()
    slot = BRepPrimAPI_MakeBox(gp_Pnt(3.0, 3.5, -0.25), 2.0, 1.0, 1.0).Shape()
    t0 = time.perf_counter()
    _shape, report = boolean_brep(drilled, slot, "difference")
    wall = time.perf_counter() - t0
    assert report["accepted"]
    show(f"262-face plate minus slot ({count_faces(drilled)} faces)",
         report, wall)


def many(n):
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt
    from brepkernel.many import boolean_brep_many

    base = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 8.0, 8.0, 0.5).Shape()
    step = 8.0 / n
    tools = [BRepPrimAPI_MakeBox(
        gp_Pnt(step/2 + i*step - 0.15, step/2 + j*step - 0.15, -0.25),
        0.3, 0.3, 1.0).Shape()
        for i in range(n) for j in range(n)]
    t0 = time.perf_counter()
    _out, summary = boolean_brep_many(base, tools, "difference",
                                      optimized=True)
    wall = time.perf_counter() - t0
    assert summary["optimized"] is True
    main_rep = summary["steps"][-1]
    assert main_rep["accepted"]
    show(f"{n*n}-hole optimized many", main_rep, wall)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    plate()
    print(flush=True)
    many(n)


if __name__ == "__main__":
    main()
