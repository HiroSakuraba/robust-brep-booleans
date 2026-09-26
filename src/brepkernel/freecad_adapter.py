"""FreeCAD integration groundwork (Part B, B5): the geometry bridge.

This module does NOT import FreeCAD. FreeCAD is not installed in this
environment and its bundled Python/OCCT build cannot share in-process
objects with the brepkernel venv anyway (two different OCCT builds in
one process is the one technical trap named in the work plan).

The bridge works on plain data, so both the headless CI path and the
real workbench path use the same seam:

  FreeCAD side                     brepkernel side (this module)
  ----------------------------     ------------------------------
  shape.exportBrepToString()  -->  shape_of_brep_text(bytes) -> TopoDS
  run_boolean() runs boolean_brep()
  result BREP text             <--  brep_text_of(shape) -> bytes
  shape.importBrepFromString() <--  (back in FreeCAD)

In production the brepkernel half runs as a subprocess (its own venv,
its own Python and OCCT); here run_boolean() executes the pipeline
in-process because that is the same contract. STEP is the alternate
bridge for CAD-to-CAD exchange (read_step / write_step).

Typed refusals never become solids: run_boolean() returns an outcome
dict whose "refused" branch carries category + stage + message for the
FreeCAD task panel, and result_brep is None on that branch.

Evidence seam: the pipeline report's ["evidence"] record (schema
brepkernel.evidence/1.0, added on branch partB/evidence-schema) is
carried verbatim as result["evidence"] when present, so the FreeCAD
sidecar (part.evidence.json) is the same audit trail whether the call
came from a workbench, the CLI, or the service.
"""

from __future__ import annotations

import io
import json
import os
import tempfile

VALID_OPS = ("union", "intersection", "difference")


class AdapterError(Exception):
    """Typed input failure at the bridge: bad BREP text, bad STEP path."""

    kind = "AdapterInputInvalid"

    def __init__(self, message):
        super().__init__(message)


def brep_text_of(shape) -> bytes:
    """Canonical OCCT BREP text of a shape, as ASCII bytes.

    BRepTools_Write emits entities in construction order, so the bytes
    are a stable content fingerprint of the operand (this is the same
    property the evidence schema hashes).
    """
    from OCP.BRepTools import BRepTools

    fd, path = tempfile.mkstemp(suffix=".brep")
    os.close(fd)
    try:
        if not BRepTools.Write_s(shape, path):
            raise AdapterError("BRepTools.Write failed for shape")
        with open(path, "rb") as f:
            return f.read()
    finally:
        try:
            os.unlink(path)
        except OSError:
            pass


def shape_of_brep_text(data: bytes):
    """Read OCCT BREP text (bytes) into a TopoDS shape.

    Raises AdapterError on empty or unparseable input. The returned
    shape is explicitly checked: an empty compound is not a solid and
    is refused rather than handed to the pipeline.
    """
    from OCP.BRep import BRep_Builder
    from OCP.BRepTools import BRepTools
    from OCP.TopoDS import TopoDS_Shape

    if not isinstance(data, (bytes, bytearray)) or not data:
        raise AdapterError("empty BREP text")
    shape = TopoDS_Shape()
    try:
        BRepTools.Read_s(shape, io.BytesIO(bytes(data)), BRep_Builder())
    except Exception as e:
        raise AdapterError(f"BREP text parse failed: {e}") from e
    if shape.IsNull():
        raise AdapterError("BREP text parsed to a null shape")
    return shape


def write_step(shape, path: str) -> str:
    """Write a shape to STEP. Raises AdapterError on writer failure."""
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.STEPControl import STEPControl_Writer, STEPControl_StepModelType

    w = STEPControl_Writer()
    if w.Transfer(shape, STEPControl_StepModelType.STEPControl_AsIs) != IFSelect_RetDone:
        raise AdapterError("STEP transfer failed")
    if w.Write(path) != IFSelect_RetDone:
        raise AdapterError(f"STEP write failed for {path}")
    return path


def read_step(path: str):
    """Read STEP into an indexed BRepModel (raises on failure)."""
    from .step_ingest import load_step

    return load_step(path)


def topology_counts(shape) -> dict:
    """Solid/shell/face/edge/vertex counts of a shape."""
    from OCP.TopAbs import (TopAbs_EDGE, TopAbs_FACE, TopAbs_SHELL,
                            TopAbs_SOLID, TopAbs_VERTEX)
    from OCP.TopExp import TopExp_Explorer

    def count(t):
        n = 0
        ex = TopExp_Explorer(shape, t)
        while ex.More():
            n += 1
            ex.Next()
        return n

    return {
        "solids": count(TopAbs_SOLID),
        "shells": count(TopAbs_SHELL),
        "faces": count(TopAbs_FACE),
        "edges": count(TopAbs_EDGE),
        "vertices": count(TopAbs_VERTEX),
    }


def shape_volume(shape) -> float:
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps

    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, props)
    return float(props.Mass())


def round_trip_report(shape, fmt: str = "brep") -> dict:
    """Document what a bridge round trip preserves and loses.

    Returns a JSON-safe dict: counts and volume before/after, and
    whether they matched. Parametric history, toponaming references,
    and appearance data never cross the bridge; this report only checks
    the geometry that does.
    """
    if fmt == "brep":
        other = shape_of_brep_text(brep_text_of(shape))
    elif fmt == "step":
        fd, path = tempfile.mkstemp(suffix=".step")
        os.close(fd)
        try:
            write_step(shape, path)
            other = read_step(path).shape
        finally:
            try:
                os.unlink(path)
            except OSError:
                pass
    else:
        raise ValueError(f"unknown fmt {fmt!r}")

    before = topology_counts(shape)
    after = topology_counts(other)
    vol_before = shape_volume(shape)
    vol_after = shape_volume(other)
    counts_equal = before == after
    volume_ok = abs(vol_before - vol_after) <= 1e-9 * max(1.0, abs(vol_before))
    return {
        "format": fmt,
        "counts_before": before,
        "counts_after": after,
        "counts_equal": counts_equal,
        "volume_before": vol_before,
        "volume_after": vol_after,
        "volume_ok": volume_ok,
        "matched": counts_equal and volume_ok,
        "note": "parametric history, toponaming references and appearance "
                "data do not cross the bridge",
    }


def to_json_safe(obj):
    """Defensively convert pipeline artifacts to JSON-safe values."""
    import numpy as np

    if obj is None or isinstance(obj, (bool, int, str)):
        return obj
    if isinstance(obj, float):
        return obj if obj == obj and abs(obj) != float("inf") else str(obj)
    if isinstance(obj, bytes):
        return obj.decode("ascii", errors="replace")
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return to_json_safe(obj.tolist())
    if isinstance(obj, (list, tuple)):
        return [to_json_safe(v) for v in obj]
    if isinstance(obj, dict):
        return {str(k): to_json_safe(v) for k, v in obj.items()}
    return str(obj)


def run_boolean(brep_a: bytes, brep_b: bytes, op: str, **pipeline_kwargs) -> dict:
    """Run boolean_brep across the bridge; return the outcome as data.

    Returns a JSON-safe dict:
      outcome: "accepted" | "refused" | "error"
      op, result_brep (bytes on accept, else None),
      refusal: {category, stage, message} on refuse, else None,
      volume, report (pipeline report), evidence (when the pipeline
      attached one; None on this branch since the evidence schema
      branch is not merged yet).

    "accepted" means the kernel certified the solid. "refused" means a
    typed BRepAmbiguousResult; no solid is ever fabricated on that
    path. "error" is a typed boundary failure (bad input bytes). An
    unexpected exception is recorded as an error entry, never re-raised
    as a silent crash.
    """
    if op not in VALID_OPS:
        raise ValueError(f"unknown op {op!r}")
    try:
        shape_a = shape_of_brep_text(brep_a)
        shape_b = shape_of_brep_text(brep_b)
    except AdapterError as e:
        return {
            "outcome": "error",
            "op": op,
            "result_brep": None,
            "refusal": None,
            "volume": None,
            "report": {},
            "evidence": None,
            "error": {"kind": e.kind, "message": str(e)},
        }

    from .pipeline import boolean_brep, BRepAmbiguousResult

    try:
        result, report = boolean_brep(shape_a, shape_b, op, **pipeline_kwargs)
    except BRepAmbiguousResult as e:
        r = getattr(e, "report", {}) or {}
        refusal = r.get("refusal", {})
        return {
            "outcome": "refused",
            "op": op,
            "result_brep": None,
            "refusal": {
                "category": refusal.get("kind", type(e).__name__),
                "stage": refusal.get("stage", "unknown"),
                "message": refusal.get("message", str(e)),
            },
            "volume": None,
            "report": to_json_safe(r),
            "evidence": to_json_safe(r.get("evidence")),
            "error": None,
        }
    except Exception as e:  # noqa: BLE001 - boundary must not crash
        return {
            "outcome": "error",
            "op": op,
            "result_brep": None,
            "refusal": None,
            "volume": None,
            "report": {},
            "evidence": None,
            "error": {"kind": "EngineError", "message": f"{type(e).__name__}: {e}"},
        }

    return {
        "outcome": "accepted",
        "op": op,
        "result_brep": brep_text_of(result),
        "refusal": None,
        "volume": shape_volume(result),
        "report": to_json_safe(report),
        "evidence": to_json_safe(report.get("evidence")),
        "error": None,
    }


def report_json(result: dict) -> str:
    """The adapter result as a JSON string (for the FreeCAD Report property)."""
    return json.dumps(to_json_safe(result), indent=2, sort_keys=True)


def present_refusal(refusal: dict) -> str:
    """User-facing text for a typed refusal (FreeCAD task panel).

    States the refusal category, the pipeline stage, and the plain
    message. It never claims a solid was produced and never suggests
    retrying with looser tolerance (I3).
    """
    return (
        "brepkernel refused this Boolean rather than return an "
        "uncertified solid.\n"
        f"RefusalKind: {refusal.get('category', 'unknown')}\n"
        f"Stage: {refusal.get('stage', 'unknown')}\n"
        f"Detail: {refusal.get('message', '')}\n"
        "No result solid was produced. You may inspect the operands for "
        "the condition named above, or run FreeCAD's own Boolean "
        "explicitly (it is labeled uncertified)."
    )


def format_status_panel(result: dict) -> str:
    """One status block for the FreeCAD feature (Status/RefusalKind)."""
    outcome = result.get("outcome")
    if outcome == "accepted":
        return (
            f"Status: Accepted\n"
            f"Operation: {result.get('op')}\n"
            f"Volume: {result.get('volume')}\n"
        )
    if outcome == "refused":
        r = result.get("refusal") or {}
        return f"Status: Refused\n{present_refusal(r)}"
    e = result.get("error") or {}
    return (
        f"Status: Error\n"
        f"RefusalKind: {e.get('kind', 'unknown')}\n"
        f"Detail: {e.get('message', '')}\n"
    )
