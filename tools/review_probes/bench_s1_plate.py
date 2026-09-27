#!/usr/bin/env python3
"""S1 benchmark: raw vs prepared on the 262-face plate minus slot.

Measures:
- prepare_brep time for plate and slot;
- raw boolean_brep (internal prepare) best-of-3 with counters;
- prepared boolean_brep (prepare once outside) best-of-3 with counters.

Pass criterion: prepared runs show zero face-box / edge-index rebuilds
during the Boolean (the one build happens once in prepare_brep).
"""
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                           "..", "..", "src"))

CACHE = os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "goals",
    "robust-brep-booleans-prototype", "hidden_files", "plate256.brep"))


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


def load_plate():
    from OCP.TopoDS import TopoDS_Shape
    from OCP.BRepTools import BRepTools
    from OCP.BRep import BRep_Builder
    if os.path.exists(CACHE):
        shp = TopoDS_Shape()
        BRepTools.Read_s(shp, CACHE, BRep_Builder())
        # Guard against a stale/empty cache (seen once: overlayfs lost
        # the file mid-session and Read_s silently yields an empty shape).
        from OCP.TopAbs import TopAbs_FACE
        from OCP.TopExp import TopExp_Explorer
        ex = TopExp_Explorer(shp, TopAbs_FACE)
        if ex.More():
            print(f"loaded cached plate from {CACHE}", flush=True)
            return shp
        print("cached plate unreadable/empty; re-drilling", flush=True)
    else:
        print("no cached plate; drilling with OCCT ...", flush=True)
    t0 = time.perf_counter()
    drilled = build_drilled()
    print(f"drilled in {time.perf_counter() - t0:.1f}s", flush=True)
    try:
        os.makedirs(os.path.dirname(CACHE), exist_ok=True)
        BRepTools.Write_s(drilled, CACHE)
        print(f"cached to {CACHE}", flush=True)
    except OSError as e:
        print(f"cache write failed ({e}); continuing", flush=True)
    return drilled


def main():
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt
    from brepkernel.pipeline import boolean_brep
    from brepkernel.prepared import prepare_brep

    drilled = load_plate()
    slot = BRepPrimAPI_MakeBox(gp_Pnt(3.0, 3.5, -0.25), 2.0, 1.0, 1.0).Shape()

    t0 = time.perf_counter()
    p_plate = prepare_brep(drilled)
    t_plate = time.perf_counter() - t0
    t0 = time.perf_counter()
    p_slot = prepare_brep(slot)
    t_slot = time.perf_counter() - t0
    print(f"prepare_brep: plate={t_plate:.2f}s slot={t_slot:.2f}s "
          f"(faces={p_plate.n_faces} edges={p_plate.n_edges})", flush=True)

    def counters_of(rep):
        return rep["performance"]["counters"]

    # Raw: internal prepare each call.
    best_raw, rep_raw = None, None
    for k in range(3):
        t0 = time.perf_counter()
        _s, rep = boolean_brep(drilled, slot, "difference",
                               collect_perf=True)
        dt = time.perf_counter() - t0
        print(f"raw run {k}: {dt:.2f}s", flush=True)
        if best_raw is None or dt < best_raw:
            best_raw, rep_raw = dt, rep
    c = counters_of(rep_raw)
    print(f"RAW best={best_raw:.2f}s face_box_build={c['face_box_build']} "
          f"edge_box_build={c['edge_box_build']} "
          f"prepared_face_box_hit={c['prepared_face_box_hit']} "
          f"prepared_edge_index_hit={c['prepared_edge_index_hit']}",
          flush=True)

    # Prepared: prepare once outside.
    best_pre, rep_pre = None, None
    for k in range(3):
        t0 = time.perf_counter()
        _s, rep = boolean_brep(p_plate, p_slot, "difference",
                               collect_perf=True)
        dt = time.perf_counter() - t0
        print(f"prepared run {k}: {dt:.2f}s", flush=True)
        if best_pre is None or dt < best_pre:
            best_pre, rep_pre = dt, rep
    c = counters_of(rep_pre)
    vol = rep_pre["stages"]["assembly"].get("volume")
    print(f"PREPARED best={best_pre:.2f}s volume={vol:.4f} "
          f"face_box_build={c['face_box_build']} "
          f"edge_box_build={c['edge_box_build']} "
          f"prepared_face_box_hit={c['prepared_face_box_hit']} "
          f"prepared_edge_index_hit={c['prepared_edge_index_hit']}",
          flush=True)
    print(f"verdict identical: "
          f"{rep_raw['accepted'] == rep_pre['accepted']}", flush=True)


if __name__ == "__main__":
    main()
