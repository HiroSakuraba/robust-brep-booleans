"""Placed parts -> one named STEP file, checked before it replaces anything (``brepkernel export``).

The parts list is JSON::

    {"parts": [{"path": "part.brep", "name": "Bracket",
                "transform": [m00, m01, m02, m03, m10, ..., m23]}, ...]}

Each ``transform`` is a 3x4 row-major affine map from the part's own frame (a project ``.brep``:
meters, +Y up) into the scene (meters, +Y up). Moves, turns, mirrors and uniform scales are
applied exactly (the surfaces keep their type). Any other map (a stretch) is applied as a general
affine map: the part is still a valid solid of the right size, but a stretched cylinder is no
longer a cylinder, so its surfaces become B-splines and the part is reported ``exact: false``.

Checks, all before the output is replaced:

* every placed part is a valid solid (OCCT's checker) whose volume is ``|det|`` times the
  source's volume (relative 1e-6);
* the file is written next to the output under a temporary name, read back, and must hold the
  same number of solids with the same total volume (relative 1e-6);
* only then is it moved over the output. A failed check removes the temporary file and leaves
  any existing output untouched.

Units and axes: ``units`` "mm" or "m" (the STEP length unit), ``up`` "z" (most CAD programs) or
"y" (as in the editor). Every solid is named after its part (the STEP product name).
"""
from __future__ import annotations

import json
import math
import os

from . import _occt_compat  # noqa: F401

VOLUME_RTOL = 1e-6


class StepExportError(Exception):
    pass


def _read_parts(manifest_path: str) -> list[dict]:
    try:
        with open(manifest_path) as fh:
            data = json.load(fh)
    except (OSError, ValueError) as exc:
        raise StepExportError(f"cannot read the parts list {manifest_path}: {exc}") from exc
    parts = data.get("parts") if isinstance(data, dict) else None
    if not isinstance(parts, list) or not parts:
        raise StepExportError("the parts list has no parts")
    out = []
    for k, p in enumerate(parts):
        if not isinstance(p, dict) or not isinstance(p.get("path"), str):
            raise StepExportError(f"part {k}: needs a path")
        name = str(p.get("name") or os.path.splitext(os.path.basename(p["path"]))[0])
        t = p.get("transform", [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0])
        if (not isinstance(t, list) or len(t) != 12
                or not all(isinstance(v, (int, float)) and math.isfinite(v) for v in t)):
            raise StepExportError(f"{name}: transform needs 12 finite numbers (3x4 row-major)")
        out.append({"path": p["path"], "name": name, "transform": [float(v) for v in t]})
    return out


def _check_solid(shape, name: str) -> None:
    from OCP.BRepCheck import BRepCheck_Analyzer
    from .cli import _solid_count
    if shape is None or shape.IsNull() or _solid_count(shape) == 0:
        raise StepExportError(f"{name}: not a solid")
    if not BRepCheck_Analyzer(shape).IsValid():
        raise StepExportError(f"{name}: the placed part does not pass OCCT's validity check")


def _to_export_frame(shape, up: str):
    """Scene frame (+Y up) to the file's axes: for Z up, +Y becomes +Z (exact rotation)."""
    if up == "y":
        return shape
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
    from OCP.gp import gp_Ax1, gp_Dir, gp_Pnt, gp_Trsf
    t = gp_Trsf()
    t.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)), math.pi / 2)
    return BRepBuilderAPI_Transform(shape, t, True).Shape()


def _write_named_step(named_shapes, path: str, units: str) -> None:
    """Write [(name, shape in meters)] as one STEP file in `units`, one named product each."""
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.STEPCAFControl import STEPCAFControl_Writer
    from OCP.STEPControl import STEPControl_AsIs
    from OCP.TCollection import TCollection_ExtendedString
    from OCP.TDataStd import TDataStd_Name
    from OCP.TDocStd import TDocStd_Document
    from OCP.XCAFDoc import XCAFDoc_DocumentTool
    from ._occt_compat import occt_major_minor
    from .cli import frame_shape
    # Shapes are scaled to the file's unit here (exact), and the document unit says so, so OCCT
    # has nothing to convert. OCP 7.8's XCAF writer converts to millimetres whatever the STEP
    # unit and only relabels: for meters it gets a document declared in millimetres, which it
    # then copies unchanged. The read-back check (volume and declared unit) guards both.
    factor = 1000.0 if units == "mm" else 1.0
    doc_unit = 0.001 if (units == "m" and occt_major_minor() < (8, 0)) else 1.0 / factor
    doc = TDocStd_Document(TCollection_ExtendedString("XmlXCAF"))
    XCAFDoc_DocumentTool.SetLengthUnit_s(doc, doc_unit)          # meters per document unit
    tool = XCAFDoc_DocumentTool.ShapeTool_s(doc.Main())
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Copy
    for name, shape in named_shapes:
        # a copy per part: XCAF would merge two placements that share one shape into one product
        part = BRepBuilderAPI_Copy(shape, True, False).Shape()
        label = tool.AddShape(frame_shape(part, factor, False), False)
        TDataStd_Name.Set_s(label, TCollection_ExtendedString(name, True))
    w = STEPCAFControl_Writer()      # first: it registers the STEP settings (OCP 7.8 ignores them before)
    with _settings({"write.step.unit": "MM" if units == "mm" else "M", "write.step.schema": "AP214IS"}):
        w.SetNameMode(True)
        if not w.Transfer(doc, STEPControl_AsIs):
            raise StepExportError("STEP transfer failed")
        if w.Write(path) != IFSelect_RetDone:
            raise StepExportError(f"cannot write {path}")
    _ascii_escape_file(path)


class _settings:
    """Set OCCT's process-wide data-exchange settings for one write or read, then restore them,
    so a later STEP write in the same process (brepkernel or its caller) is not affected."""

    def __init__(self, values: dict):
        self.values, self.saved = values, {}

    def __enter__(self):
        from OCP.Interface import Interface_Static
        for key, value in self.values.items():
            self.saved[key] = Interface_Static.CVal_s(key) if Interface_Static.IsPresent_s(key) else None
            if not Interface_Static.SetCVal_s(key, value):
                raise StepExportError(f"cannot set {key}")
        return self

    def __exit__(self, *exc):
        from OCP.Interface import Interface_Static
        for key, value in self.saved.items():
            if value is not None:
                Interface_Static.SetCVal_s(key, value)
        return False


def _ascii_escape(text: str) -> str:
    """ISO 10303-21 control directives for every non-ASCII run: \\X2\\hhhh..\\X0\\ (UTF-16),
    or \\X4\\hhhhhhhh..\\X0\\ when the run holds characters beyond U+FFFF."""
    out, i = [], 0
    while i < len(text):
        if ord(text[i]) < 128:
            out.append(text[i])
            i += 1
            continue
        j = i
        while j < len(text) and ord(text[j]) >= 128:
            j += 1
        run = text[i:j]
        if all(ord(c) <= 0xFFFF for c in run):
            out.append("\\X2\\" + "".join(f"{ord(c):04X}" for c in run) + "\\X0\\")
        else:
            out.append("\\X4\\" + "".join(f"{ord(c):08X}" for c in run) + "\\X0\\")
        i = j
    return "".join(out)


def _ascii_escape_file(path: str) -> None:
    """OCCT writes non-ASCII names as raw UTF-8; rewrite them in the portable escaped form."""
    with open(path, "rb") as fh:
        data = fh.read()
    try:
        data.decode("ascii")
        return
    except UnicodeDecodeError:
        pass
    text = _ascii_escape(data.decode("utf-8"))
    with open(path, "wb") as fh:
        fh.write(text.encode("ascii"))


_UNIT_DECLARATION = {"mm": "SI_UNIT(.MILLI.,.METRE.)", "m": "SI_UNIT($,.METRE.)"}


def _read_back(path: str):
    """(solids, total volume in m^3, part names, declared length unit) of a written STEP file.

    Geometry is read with the plain reader set to millimetres; names with the XCAF reader (as a
    CAD program would see them). Both settings are made after the reader exists, because OCCT
    registers them only then."""
    import re
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.STEPCAFControl import STEPCAFControl_Reader
    from OCP.STEPControl import STEPControl_Reader
    from OCP.TCollection import TCollection_ExtendedString
    from OCP.TDataStd import TDataStd_Name
    from OCP.TDF import TDF_ChildIterator
    from OCP.TDocStd import TDocStd_Document
    from OCP.XCAFDoc import XCAFDoc_DocumentTool
    from .cli import _solid_count, _volume
    r = STEPControl_Reader()
    with _settings({"xstep.cascade.unit": "MM"}):
        if r.ReadFile(path) != IFSelect_RetDone:
            raise StepExportError("read-back check failed: the written file cannot be read")
        r.TransferRoots()
        shape = r.OneShape()
    solids, volume = _solid_count(shape), _volume(shape) * 1e-9

    x = STEPCAFControl_Reader()
    x.SetNameMode(True)
    names = []
    doc = TDocStd_Document(TCollection_ExtendedString("XmlXCAF"))
    if x.ReadFile(path) == IFSelect_RetDone and x.Transfer(doc):
        tool = XCAFDoc_DocumentTool.ShapeTool_s(doc.Main())
        it = TDF_ChildIterator(tool.Label(), False)
        while it.More():
            label = it.Value()
            it.Next()
            attr = TDataStd_Name()
            if label.FindAttribute(TDataStd_Name.GetID_s(), attr):
                names.append(attr.Get().ToExtString())

    with open(path, "rb") as fh:
        text = fh.read().decode("ascii", "replace")
    units = sorted(k for k, decl in _UNIT_DECLARATION.items()
                   if re.search(r"LENGTH_UNIT\(\)\s*NAMED_UNIT\(\*\)\s*" + re.escape(decl), text))
    return solids, volume, names, units


def export_parts(manifest_path: str, output: str, *, units: str = "mm", up: str = "z") -> dict:
    """Place, check, write and read back. Returns the report (see the module docstring)."""
    import numpy as np
    from .cli import CliError, _volume, frame_shape, place_shape, read_shape
    if units not in ("mm", "m") or up not in ("z", "y"):
        raise StepExportError("units must be mm or m, up must be z or y")
    if not output.lower().endswith((".step", ".stp")):
        raise StepExportError("the output must end in .step or .stp")
    parts = _read_parts(manifest_path)
    placed, records = [], []
    for p in parts:
        try:
            shape, fmt = read_shape(p["path"])
        except CliError as exc:
            raise StepExportError(f"{p['name']}: {exc}") from exc
        if fmt == "step":                                  # like `import`: mm, Z up -> m, Y up
            shape = frame_shape(shape, 0.001, True)
        _check_solid(shape, p["name"])
        source_volume = _volume(shape)
        try:
            shape, kind = place_shape(shape, p["transform"])
        except CliError as exc:
            raise StepExportError(f"{p['name']}: {str(exc).replace('--b-transform', 'transform')}") from exc
        _check_solid(shape, p["name"])
        det = abs(float(np.linalg.det(np.array(p["transform"]).reshape(3, 4)[:, :3])))
        volume = _volume(shape)
        expected = det * source_volume
        if not abs(volume - expected) <= VOLUME_RTOL * max(expected, 1e-300):
            raise StepExportError(f"{p['name']}: placed volume {volume:.9g} m^3, expected {expected:.9g} m^3")
        placed.append((p["name"], _to_export_frame(shape, up)))
        records.append({"name": p["name"], "path": os.path.abspath(p["path"]), "placement": kind,
                        "exact": kind != "general_affine", "volume_m3": volume,
                        "volume_mm3": volume * 1e9})
    total = sum(r["volume_m3"] for r in records)

    out_dir = os.path.dirname(os.path.abspath(output))
    os.makedirs(out_dir, exist_ok=True)
    stem, ext = os.path.splitext(os.path.abspath(output))
    tmp = f"{stem}.brepkernel-{os.getpid()}{ext}"
    try:
        _write_named_step(placed, tmp, units)
        solids, back, names, declared = _read_back(tmp)
        if declared != [units]:
            raise StepExportError(f"read-back check failed: the file declares {declared or 'no'} length unit, "
                                  f"not {units}")
        rel = abs(back - total) / max(total, 1e-300)
        if solids != len(placed):
            raise StepExportError(f"read-back check failed: {solids} solids in the file, {len(placed)} parts")
        if not rel <= VOLUME_RTOL:
            raise StepExportError(f"read-back check failed: volume differs by {rel:.3g} (relative)")
        missing = [r["name"] for r in records if r["name"] not in names]
        if missing:
            raise StepExportError(f"read-back check failed: names missing from the file: {', '.join(missing)}")
        os.replace(tmp, output)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    return {"parts": records, "units": units, "up": up, "solids": len(placed),
            "volume_mm3": total * 1e9,
            "readback": {"solids": solids, "volume_mm3": back * 1e9, "relative_error": rel}}

