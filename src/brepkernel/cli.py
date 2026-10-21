"""brepkernel command line: exact CAD import, certified Booleans, glTF meshes.

Built for game-engine editors (DVE) that keep an exact B-rep next to each
polygon object and call this tool as an external process. Every command
writes a JSON report (--report PATH, else stdout) and uses stable exit
codes:

    0  done (Boolean accepted)
    2  Boolean refused: no output file is written; the report carries the
       refusal and, with --advisory, the uncertified second opinions
    1  usage, input or output error (report has "error")

Commands

  import IN -o OUT.brep [--mesh OUT.glb] [--scale S] [--up z|y]
      Read STEP (.step/.stp), OCCT (.brep) or OpenUSD B-rep (.usda) and
      write a project-owned .brep in the engine frame (meters, +Y up). STEP
      arrives from OCCT in millimetres and is assumed +Z up, so its defaults
      are --scale 0.001 --up z; .brep and .usda default to --scale 1 --up y.

  boolean A B --op union|intersection|difference -o OUT.brep
          [--b-transform m00,m01,m02,m03,m10,...,m23] [--reverse] [--mesh OUT.glb]
          [--advisory] [--independent-engine smlib]
      Certified Boolean. --b-transform places B in A's frame (3x4 row-major
      affine); a rigid or uniformly scaled map stays exact, anything else is
      applied as a general affine map (surfaces become B-splines; reported).
      --reverse (difference only) computes B - A, still in A's frame.

  mesh IN -o OUT.glb [--deflection D] [--angular A] [--scale S] [--up z|y]
      Watertight glTF 2.0 binary with hard B-rep edges and exact normals.

Mesh options for import/boolean/mesh: --deflection (meters in the output
frame; default 1/500 of the part's diagonal), --angular (radians, 0.35).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

EXIT_OK, EXIT_ERROR, EXIT_REFUSED = 0, 1, 2
REPORT_SCHEMA = "brepkernel.cli/1"


class CliError(Exception):
    pass


# ---------------------------------------------------------------------------
# shapes in and out
# ---------------------------------------------------------------------------

def read_shape(path: str):
    """(shape, format) for .brep, .step/.stp or .usda."""
    from OCP.BRep import BRep_Builder
    from OCP.BRepTools import BRepTools
    from OCP.TopoDS import TopoDS_Shape
    if not os.path.isfile(path):
        raise CliError(f"input not found: {path}")
    low = path.lower()
    if low.endswith(".brep"):
        s = TopoDS_Shape()
        if not BRepTools.Read_s(s, path, BRep_Builder()) or s.IsNull():
            raise CliError(f"cannot read OCCT BREP {path}")
        return s, "brep"
    if low.endswith((".step", ".stp")):
        from .step_ingest import load_step
        try:
            return load_step(path, build_freeform=False).shape, "step"
        except Exception as exc:
            raise CliError(f"cannot read STEP {path}: {exc}") from exc
    if low.endswith(".usda"):
        from .usd import read_omnisolid_usda
        try:
            return read_omnisolid_usda(path), "usda"
        except Exception as exc:
            raise CliError(f"cannot read OpenUSD B-rep {path}: {exc}") from exc
    raise CliError(f"unsupported input format (use .brep, .step, .stp, .usda): {path}")


def write_shape(shape, path: str) -> None:
    from OCP.BRepTools import BRepTools
    low = path.lower()
    tmp = path + ".part"
    if low.endswith(".brep"):
        if not BRepTools.Write_s(shape, tmp):
            raise CliError(f"cannot write {path}")
    elif low.endswith(".usda"):
        from .usd import write_omnisolid_usda
        write_omnisolid_usda(shape, tmp)
    elif low.endswith((".step", ".stp")):
        from OCP.IFSelect import IFSelect_RetDone
        from OCP.STEPControl import STEPControl_AsIs, STEPControl_Writer
        w = STEPControl_Writer()
        w.Transfer(shape, STEPControl_AsIs)
        if w.Write(tmp) != IFSelect_RetDone:
            raise CliError(f"cannot write {path}")
    else:
        raise CliError(f"unsupported output format (use .brep, .usda, .step): {path}")
    os.replace(tmp, path)              # never leave a half-written result


def _solid_count(shape) -> int:
    from OCP.TopAbs import TopAbs_SOLID
    from OCP.TopExp import TopExp_Explorer
    n, ex = 0, TopExp_Explorer(shape, TopAbs_SOLID)
    while ex.More():
        n += 1
        ex.Next()
    return n


def _volume(shape) -> float:
    if shape is None or shape.IsNull() or _solid_count(shape) == 0:
        return 0.0
    from .assembly import _all_faces_analytic, _shape_volume
    return abs(float(_shape_volume(shape, 1e-12 if _all_faces_analytic(shape) else 1e-10)))


def _sha256(path: str) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def frame_shape(shape, scale: float, from_z_up: bool):
    """Shape converted to the engine frame (exact: similarity transform)."""
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
    from OCP.gp import gp_Ax1, gp_Dir, gp_Pnt, gp_Trsf
    out = shape
    if from_z_up:
        t = gp_Trsf()
        t.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)), -1.5707963267948966)
        out = BRepBuilderAPI_Transform(out, t, True).Shape()
    if scale != 1.0:
        t = gp_Trsf()
        t.SetScale(gp_Pnt(0, 0, 0), float(scale))
        out = BRepBuilderAPI_Transform(out, t, True).Shape()
    return out


def place_shape(shape, m12):
    """Apply a 3x4 row-major affine map. Returns (shape, kind)."""
    import numpy as np
    from OCP.BRepBuilderAPI import BRepBuilderAPI_GTransform, BRepBuilderAPI_Transform
    from OCP.gp import gp_GTrsf, gp_Mat, gp_Trsf, gp_XYZ
    if len(m12) != 12 or not all(np.isfinite(m12)):
        raise CliError("--b-transform needs 12 finite numbers (3x4 row-major)")
    M = np.array(m12, dtype=float).reshape(3, 4)
    L, t = M[:, :3], M[:, 3]
    det = float(np.linalg.det(L))
    if abs(det) < 1e-300:
        raise CliError("--b-transform is singular")
    s = abs(det) ** (1.0 / 3.0)
    R = L / s
    ortho_err = float(np.abs(R @ R.T - np.eye(3)).max())
    if np.allclose(L, np.eye(3), atol=0.0) and not t.any():
        return shape, "identity"
    if ortho_err <= 1e-9:
        # rigid, uniformly scaled or mirrored: gp_Trsf keeps it exact (a
        # mirror becomes a negative scale factor; checked in the tests)
        tr = gp_Trsf()
        tr.SetValues(*(float(x) for x in M.reshape(-1)))
        return BRepBuilderAPI_Transform(shape, tr, True).Shape(), (
            "similarity" if abs(s - 1.0) > 1e-12 else "rigid") + ("_mirror" if det < 0 else "")
    g = gp_GTrsf()
    g.SetVectorialPart(gp_Mat(*(float(x) for x in L.reshape(-1))))
    g.SetTranslationPart(gp_XYZ(*(float(x) for x in t)))
    op = BRepBuilderAPI_GTransform(shape, g, True)
    if not op.IsDone():
        raise CliError("general affine placement of B failed")
    return op.Shape(), "general_affine"


def _mesh(shape, path, args, report, *, scale=1.0, from_z_up=False):
    from .mesh_export import MeshExportError, write_glb
    try:
        m = write_glb(shape, path, linear_deflection=args.deflection,
                      angular_deflection=args.angular, scale=scale, from_z_up=from_z_up,
                      exact_volume=_volume(shape),
                      name=os.path.splitext(os.path.basename(path))[0])
    except MeshExportError as exc:
        raise CliError(f"mesh export failed: {exc}") from exc
    report["mesh"] = {"path": os.path.abspath(path), "sha256": _sha256(path), **m.summary()}


# ---------------------------------------------------------------------------
# commands
# ---------------------------------------------------------------------------

def _frame_defaults(fmt, args):
    scale = args.scale if args.scale is not None else (0.001 if fmt == "step" else 1.0)
    up = args.up if args.up is not None else ("z" if fmt == "step" else "y")
    return float(scale), up


def cmd_import(args, report) -> int:
    shape, fmt = read_shape(args.input)
    scale, up = _frame_defaults(fmt, args)
    out = frame_shape(shape, scale, up == "z")
    n = _solid_count(out)
    if n == 0:
        raise CliError("input contains no solid")
    write_shape(out, args.output)
    report.update({
        "input": {"path": os.path.abspath(args.input), "format": fmt,
                  "sha256": _sha256(args.input)},
        "frame": {"scale": scale, "up": up, "target": "meters, +Y up"},
        "output": {"path": os.path.abspath(args.output), "sha256": _sha256(args.output),
                   "solids": n, "volume": _volume(out)},
    })
    if args.mesh:
        _mesh(out, args.mesh, args, report)
    return EXIT_OK


def cmd_mesh(args, report) -> int:
    shape, fmt = read_shape(args.input)
    scale = args.scale if args.scale is not None else 1.0
    up = args.up or "y"
    _mesh(shape, args.output, args, report, scale=float(scale), from_z_up=(up == "z"))
    report["input"] = {"path": os.path.abspath(args.input), "format": fmt,
                       "sha256": _sha256(args.input)}
    return EXIT_OK


def cmd_boolean(args, report) -> int:
    from . import BRepAmbiguousResult, boolean_brep
    a, _fa = read_shape(args.a)
    b, _fb = read_shape(args.b)
    placement = "identity"
    if args.b_transform:
        try:
            m12 = [float(x) for x in args.b_transform.split(",")]
        except ValueError as exc:
            raise CliError("--b-transform must be 12 comma-separated numbers") from exc
        b, placement = place_shape(b, m12)
    report["inputs"] = {
        "a": {"path": os.path.abspath(args.a), "sha256": _sha256(args.a), "volume": _volume(a)},
        "b": {"path": os.path.abspath(args.b), "sha256": _sha256(args.b), "volume": _volume(b),
              "placement": placement},
    }
    report["op"] = args.op
    kw = {}
    if args.independent_engine:
        kw["independent_engine"] = args.independent_engine
    if args.advisory:
        kw["advisory_on_refusal"] = True
    if args.reverse and args.op != "difference":
        raise CliError("--reverse is only meaningful with --op difference")
    report["reverse"] = bool(args.reverse)
    try:
        # --reverse computes B - A (B placed by --b-transform); the result is
        # still expressed in A's frame
        out, rep = (boolean_brep(b, a, "difference", **kw) if args.reverse
                    else boolean_brep(a, b, args.op, **kw))
    except BRepAmbiguousResult as exc:
        r = exc.report
        report["status"] = "refused"
        report["refusal"] = r.get("refusal")
        if "advisory" in r:
            report["advisory"] = r["advisory"]
        if "evidence" in r:
            report["evidence_name"] = r["evidence"].get("name")
        return EXIT_REFUSED
    if _solid_count(out) == 0:
        report["status"] = "accepted_empty"
        report["empty"] = True
        report["result"] = {"volume": 0.0, "solids": 0}
        return EXIT_OK
    write_shape(out, args.output)
    report["status"] = "accepted"
    report["result"] = {"path": os.path.abspath(args.output), "sha256": _sha256(args.output),
                        "solids": _solid_count(out), "volume": _volume(out)}
    ver = rep.get("stages", {}).get("verification", {})
    report["certification"] = {
        "evidence_name": rep.get("evidence", {}).get("name"),
        "tolerance_ratio": ver.get("tolerance_ratio"),
        "independent_engine": rep.get("stages", {}).get("independent_engine", {}).get("status"),
    }
    if args.mesh:
        _mesh(out, args.mesh, args, report)
    return EXIT_OK


def cmd_primitive(args, report) -> int:
    """Analytic test solids in meters, +Y up (for trying the pipeline and for tests)."""
    from OCP.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder,
                                 BRepPrimAPI_MakeSphere)
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt
    at = [float(v) for v in args.at]
    if args.kind == "box":
        size = [float(v) for v in args.size]
        if min(size) <= 0:
            raise CliError("--size must be positive")
        lo = gp_Pnt(*(c - 0.5 * d for c, d in zip(at, size)))     # centred on --at
        hi = gp_Pnt(*(c + 0.5 * d for c, d in zip(at, size)))
        shape = BRepPrimAPI_MakeBox(lo, hi).Shape()
    elif args.kind == "cylinder":
        if args.radius <= 0 or args.height <= 0:
            raise CliError("--radius and --height must be positive")
        axis = {"x": (1, 0, 0), "y": (0, 1, 0), "z": (0, 0, 1)}[args.axis]
        base = gp_Pnt(*(c - 0.5 * args.height * a for c, a in zip(at, axis)))   # centred
        shape = BRepPrimAPI_MakeCylinder(gp_Ax2(base, gp_Dir(*axis)), args.radius,
                                         args.height).Shape()
    else:
        if args.radius <= 0:
            raise CliError("--radius must be positive")
        shape = BRepPrimAPI_MakeSphere(gp_Pnt(*at), args.radius).Shape()
    write_shape(shape, args.output)
    report["output"] = {"path": os.path.abspath(args.output), "sha256": _sha256(args.output),
                        "solids": _solid_count(shape), "volume": _volume(shape)}
    if args.mesh:
        _mesh(shape, args.mesh, args, report)
    return EXIT_OK


def _mjfit_result(xml, rep, args, report) -> int:
    report["certificate"] = rep
    report["status"] = rep["status"]
    if rep["status"] != "certified" or xml is None:
        return EXIT_REFUSED
    tmp = args.output + ".part"
    with open(tmp, "w") as fh:
        fh.write(xml)
    os.replace(tmp, args.output)
    report["output"] = {"path": os.path.abspath(args.output), "sha256": _sha256(args.output)}
    return EXIT_OK


def _peg_options(args) -> dict:
    kw = {}
    for name in ("peg_length", "handle_height", "spawn_height"):
        v = getattr(args, name, None)
        if v is not None:
            kw[name] = v
    return kw


def cmd_convex_pieces(args, report) -> int:
    """Certified convex collision pieces for an extruded part (JSON)."""
    from .convex_pieces import convex_pieces
    shape, fmt = read_shape(args.input)
    scale, up = _frame_defaults(fmt, args)
    shape = frame_shape(shape, scale, up == "z")
    rep = convex_pieces(shape, max_error=args.max_error, sample_spacing=args.sample_spacing)
    report["input"] = {"path": os.path.abspath(args.input), "format": fmt, "sha256": _sha256(args.input)}
    report["frame"] = {"scale": scale, "up": up, "target": "meters, +Y up"}
    report["status"] = rep["status"]
    pieces = rep.pop("pieces", None)
    report["certificate"] = rep
    if rep["status"] != "certified":
        return EXIT_REFUSED
    tmp = args.output + ".part"
    with open(tmp, "w") as fh:
        json.dump({"schema": rep["schema"], "units": "meters", "pieces": pieces}, fh)
    os.replace(tmp, args.output)
    report["output"] = {"path": os.path.abspath(args.output), "sha256": _sha256(args.output),
                        "pieces": len(pieces)}
    return EXIT_OK


def cmd_mjfit_certify(args, report) -> int:
    """Mode A: an mjfit peg config, checked against its true outlines."""
    import yaml
    from .mujoco_fit import certify_mjfit_config
    overrides = {}
    for item in args.set or []:
        key, _, value = item.partition("=")
        if not key or not _:
            raise CliError(f"--set needs key=value, got {item!r}")
        try:
            overrides[key] = float(value)            # YAML would read 8e-5 as text
        except ValueError:
            overrides[key] = yaml.safe_load(value)
    try:
        xml, rep = certify_mjfit_config(args.config, sample_spacing=args.sample_spacing, **overrides)
    except (FileNotFoundError, ValueError) as exc:
        raise CliError(str(exc)) from exc
    report["config"] = args.config
    return _mjfit_result(xml, rep, args, report)


def cmd_mjfit_from_cad(args, report) -> int:
    """Mode C: an mjfit peg-in-hole model for a CAD plate, checked against it."""
    from .mujoco_fit import mjfit_peg_from_cad
    xml, rep = mjfit_peg_from_cad(args.input, args.clearance, max_error=args.max_error, scale=args.scale,
                                  sample_spacing=args.sample_spacing, **_peg_options(args))
    return _mjfit_result(xml, rep, args, report)


# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="brepkernel", description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p, mesh_out=True):
        p.add_argument("--report", default=None, help="JSON report path (default: stdout)")
        if mesh_out:
            p.add_argument("--mesh", default=None, help="also write a .glb mesh")
        p.add_argument("--deflection", type=float, default=None,
                       help="mesh chord deflection in meters (output frame)")
        p.add_argument("--angular", type=float, default=0.35,
                       help="mesh angular deflection in radians")

    p = sub.add_parser("import", help="CAD file -> project .brep (meters, +Y up)")
    p.add_argument("input")
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--scale", type=float, default=None)
    p.add_argument("--up", choices=["y", "z"], default=None)
    common(p)

    p = sub.add_parser("boolean", help="certified Boolean of two solids")
    p.add_argument("a")
    p.add_argument("b")
    p.add_argument("--op", required=True, choices=["union", "intersection", "difference"])
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--b-transform", default=None)
    p.add_argument("--reverse", action="store_true", help="with difference: B - A")
    p.add_argument("--advisory", action="store_true")
    p.add_argument("--independent-engine", choices=["smlib"], default=None)
    common(p)

    p = sub.add_parser("mesh", help="B-rep -> watertight .glb")
    p.add_argument("input")
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--scale", type=float, default=None)
    p.add_argument("--up", choices=["y", "z"], default=None)
    common(p, mesh_out=False)

    p = sub.add_parser("primitive", help="analytic test solid (box, cylinder, sphere), centred on --at")
    p.add_argument("kind", choices=["box", "cylinder", "sphere"])
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--at", nargs=3, type=float, default=[0.0, 0.0, 0.0], metavar=("X", "Y", "Z"))
    p.add_argument("--size", nargs=3, type=float, default=[1.0, 1.0, 1.0], metavar=("DX", "DY", "DZ"))
    p.add_argument("--radius", type=float, default=0.5)
    p.add_argument("--height", type=float, default=1.0)
    p.add_argument("--axis", choices=["x", "y", "z"], default="y")
    common(p)

    p = sub.add_parser("convex-pieces", help="extruded part -> certified convex collision pieces (JSON)")
    p.add_argument("input")
    p.add_argument("-o", "--output", required=True, help="pieces JSON (written only when certified)")
    p.add_argument("--max-error", type=float, default=1e-4, help="polygon error for arcs, meters")
    p.add_argument("--scale", type=float, default=None)
    p.add_argument("--up", choices=["y", "z"], default=None)
    p.add_argument("--sample-spacing", type=float, default=None)
    p.add_argument("--report", default=None)

    p = sub.add_parser("mjfit-certify", help="mjfit peg config -> MJCF, checked against the true outlines")
    p.add_argument("config", help="bundled mjfit config name, YAML path")
    p.add_argument("-o", "--output", required=True, help="MJCF file (written only when certified)")
    p.add_argument("--set", action="append", metavar="KEY=VALUE", help="override a config value (YAML)")
    p.add_argument("--sample-spacing", type=float, default=None)
    p.add_argument("--report", default=None)

    p = sub.add_parser("mjfit-from-cad", help="CAD plate with a hole -> mjfit MJCF peg-in-hole, checked")
    p.add_argument("input", help=".step/.stp (millimeters), .brep or .usda (meters)")
    p.add_argument("-o", "--output", required=True, help="MJCF file (written only when certified)")
    p.add_argument("--clearance", type=float, required=True, help="peg to hole gap, meters")
    p.add_argument("--max-error", type=float, default=1e-4, help="polygon error for arcs, meters")
    p.add_argument("--scale", type=float, default=None, help="override the unit scale to meters")
    p.add_argument("--peg-length", type=float, default=None)
    p.add_argument("--handle-height", type=float, default=None)
    p.add_argument("--spawn-height", type=float, default=None)
    p.add_argument("--sample-spacing", type=float, default=None)
    p.add_argument("--report", default=None)
    return ap


def main(argv=None) -> int:
    ap = build_parser()
    try:
        args = ap.parse_args(argv)
    except SystemExit as exc:            # argparse usage error
        return EXIT_ERROR if exc.code else EXIT_OK
    report = {"schema": REPORT_SCHEMA, "command": args.cmd}
    t0 = time.perf_counter()
    try:
        from importlib.metadata import version
        report["brepkernel_version"] = version("brepkernel")
    except Exception:
        report["brepkernel_version"] = None
    # OCCT and usd-brep print to the C-level stdout/stderr (STEP transfer
    # statistics, pxr import warnings). Capture both while the command runs
    # so stdout carries only the JSON report; keep the tail on errors.
    import tempfile
    sys.stdout.flush()
    sys.stderr.flush()
    saved = os.dup(1), os.dup(2)
    log = tempfile.TemporaryFile()
    os.dup2(log.fileno(), 1)
    os.dup2(log.fileno(), 2)
    try:
        code = {"import": cmd_import, "boolean": cmd_boolean, "mesh": cmd_mesh,
                "primitive": cmd_primitive, "mjfit-certify": cmd_mjfit_certify,
                "convex-pieces": cmd_convex_pieces,
                "mjfit-from-cad": cmd_mjfit_from_cad}[args.cmd](args, report)
    except CliError as exc:
        report["status"] = "error"
        report["error"] = str(exc)
        code = EXIT_ERROR
    except Exception as exc:              # unexpected: still a structured report
        report["status"] = "error"
        report["error"] = f"{type(exc).__name__}: {exc}"
        code = EXIT_ERROR
    finally:
        sys.stdout.flush()
        sys.stderr.flush()
        os.dup2(saved[0], 1)
        os.dup2(saved[1], 2)
        os.close(saved[0])
        os.close(saved[1])
    if code == EXIT_ERROR:
        log.seek(0)
        tail = log.read()[-4000:].decode("utf-8", "replace")
        lines = [ln for ln in tail.splitlines()
                 if ln.strip() and "pxr" not in ln and "_LoadLibModules" not in ln]
        if lines:
            report["log_tail"] = lines[-20:]
    log.close()
    report.setdefault("status", "ok" if code == EXIT_OK else report.get("status"))
    report["exit_code"] = code
    report["seconds"] = round(time.perf_counter() - t0, 3)
    text = json.dumps(report, indent=1, default=str)
    if args.report:
        tmp = args.report + ".part"
        with open(tmp, "w") as fh:
            fh.write(text)
        os.replace(tmp, args.report)
    else:
        print(text)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
