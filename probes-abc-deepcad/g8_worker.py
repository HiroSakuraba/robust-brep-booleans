#!/usr/bin/env python3
"""G8 single-case worker: one (modelA, modelB, op) through boolean_brep().

Reads a JSON case file:
  {"case_id": str, "corpus": str, "model_a": str, "model_b": str,
   "path_a": str, "path_b": str, "op": "union|difference|intersection",
   "out": str}
Runs it in this process (the driver isolates each case in its own process
with a 120 s join so a hang or crash cannot wedge the sweep) and writes a
JSON result to "out".

B is translated by a rigid motion so its bbox center lands near A's bbox
center: guaranteed overlap, no tangency-by-construction, same geometry.
"""

import json
import os
import sys
import time
import traceback

sys.path.insert(0, os.environ["G8_SRC"])


def load_shape(path):
    from OCP.STEPControl import STEPControl_Reader
    from OCP.IFSelect import IFSelect_RetDone
    r = STEPControl_Reader()
    if r.ReadFile(path) != IFSelect_RetDone:
        raise RuntimeError(f"STEP read failed: {path}")
    if r.TransferRoots() == 0:
        raise RuntimeError(f"STEP has no transferable roots: {path}")
    return r.OneShape()


def bbox_center_diag(shape):
    import numpy as np
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib
    box = Bnd_Box()
    BRepBndLib.Add_s(shape, box)
    xmin, ymin, zmin = box.GetXMin(), box.GetYMin(), box.GetZMin()
    xmax, ymax, zmax = box.GetXMax(), box.GetYMax(), box.GetZMax()
    lo = np.array([xmin, ymin, zmin], dtype=float)
    hi = np.array([xmax, ymax, zmax], dtype=float)
    return (lo + hi) / 2.0, float(np.linalg.norm(hi - lo))


def translate(shape, vec):
    from OCP.gp import gp_Trsf, gp_Vec
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
    t = gp_Trsf()
    t.SetTranslation(gp_Vec(float(vec[0]), float(vec[1]), float(vec[2])))
    return BRepBuilderAPI_Transform(shape, t, True).Shape()


def occt_volume(shape):
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    from OCP.TopAbs import TopAbs_SOLID
    from OCP.TopExp import TopExp_Explorer
    if shape.IsNull():
        return 0.0
    if not TopExp_Explorer(shape, TopAbs_SOLID).More():
        return 0.0
    p = GProp_GProps()
    err = BRepGProp.VolumePropertiesGK_s(
        shape, p, 1e-10, True, True, False, False, False)
    if float(err) < 0.0:
        raise RuntimeError("BRepGProp reported a negative error")
    return float(p.Mass())


def solid_count(shape):
    from OCP.TopAbs import TopAbs_SOLID
    from OCP.TopExp import TopExp_Explorer
    n = 0
    ex = TopExp_Explorer(shape, TopAbs_SOLID)
    while ex.More():
        n += 1
        ex.Next()
    return n


def prune(obj, depth=0, cap=200):
    """Make a report JSON-serializable and small."""
    if depth > 4:
        return "..."
    if isinstance(obj, dict):
        out = {}
        for k, v in obj.items():
            if isinstance(v, (list, tuple)) and len(v) > cap:
                out[k] = f"[{len(v)} items]"
            else:
                out[k] = prune(v, depth + 1, cap)
        return out
    if isinstance(obj, (list, tuple)):
        if len(obj) > cap:
            return f"[{len(obj)} items]"
        return [prune(v, depth + 1, cap) for v in obj]
    if isinstance(obj, float):
        if obj != obj or obj in (float("inf"), float("-inf")):
            return str(obj)
        return obj
    if isinstance(obj, (str, int, bool)) or obj is None:
        return obj
    return str(obj)


def main():
    case = json.load(open(sys.argv[1]))
    out = {"case_id": case["case_id"], "corpus": case["corpus"],
           "pair": case["pair"], "model_a": case["model_a"],
           "model_b": case["model_b"],
           "op": case["op"], "outcome": None}
    t0 = time.perf_counter()
    try:
        shape_a = load_shape(case["path_a"])
        shape_b = load_shape(case["path_b"])
        out["load_s"] = time.perf_counter() - t0
        ca, da = bbox_center_diag(shape_a)
        cb, _db = bbox_center_diag(shape_b)
        shift = (ca - cb) + 0.25 * da * __import__("numpy").array(
            [1.0, 0.0, 0.0])
        shape_b = translate(shape_b, shift)
        out["shift_applied"] = [float(x) for x in shift]
        out["diag_a"] = da

        from brepkernel import boolean_brep, BRepAmbiguousResult
        t1 = time.perf_counter()
        try:
            result, report = boolean_brep(shape_a, shape_b, case["op"])
        finally:
            out["boolean_s"] = time.perf_counter() - t1
        out["outcome"] = "accepted"
        out["result_solids"] = solid_count(result)
        out["result_volume"] = occt_volume(result)
        out["volume_a"] = occt_volume(shape_a)
        out["volume_b"] = occt_volume(shape_b)
        stages = report.get("stages", {})
        out["stage_keys"] = sorted(stages.keys())
        out["broadphase_max_tolerance"] = report.get(
            "broadphase_max_tolerance")
        out["timings_ms"] = report.get("timings_ms")
        ix = stages.get("intersection", {})
        out["intersection_summary"] = prune({
            k: ix.get(k) for k in (
                "n_sections", "section_tol", "max_section_tol_used",
                "accepted_section_tol", "section_tolerance",
                "n_contact_faces") if k in ix})
        asm = stages.get("assembly", {})
        out["assembly_summary"] = prune({
            k: asm.get(k) for k in (
                "n_result_faces", "n_result_solids", "volume",
                "volume_tolerance", "keep_provenance") if k in asm})
    except Exception as e:  # noqa: BLE001 - record everything
        out["outcome"] = "crashed" if not hasattr(e, "report") else None
        tb = traceback.format_exc(limit=30)
        try:
            from brepkernel import BRepAmbiguousResult
            if isinstance(e, BRepAmbiguousResult):
                ref = (e.report or {}).get("refusal", {})
                out["outcome"] = "refused"
                out["refusal_stage"] = ref.get("stage")
                out["refusal_kind"] = ref.get("kind")
                out["refusal_type"] = ref.get("type")
                out["refusal_message"] = str(ref.get("message"))[:500]
                out["refusal_evidence"] = prune(
                    ref.get("evidence")) if ref.get(
                    "evidence") is not None else None
                out["broadphase_max_tolerance"] = (e.report or {}).get(
                    "broadphase_max_tolerance")
                out["stage_keys"] = sorted(
                    (e.report or {}).get("stages", {}).keys())
        except Exception:  # noqa: BLE001
            pass
        if out["outcome"] is None:
            out["outcome"] = "exception"
        out["exception_type"] = type(e).__name__
        out["exception_message"] = str(e)[:500]
        out["traceback"] = tb
    out["wall_s"] = time.perf_counter() - t0
    with open(case["out"], "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
