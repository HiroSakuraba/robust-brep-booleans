"""Certified MuJoCo collision geometry with mjfit (optional: pip install mjfit).

`mjfit <https://github.com/younghyopark/mjfit>`_ (MIT) builds MuJoCo
peg-in-hole models from convex prisms whose union is the part, with the peg
exactly one clearance from the hole. This module checks that claim against
the true parts and feeds mjfit from CAD.

* :func:`certify_mjfit_peg` (mode A): run ``mjfit.peg`` and check the model
  against the true outlines rebuilt from the same shape specs (real circles,
  stadiums and Bezier curves).
* :func:`mjfit_peg_from_cad` (mode C): read a CAD plate (STEP, BREP, USD),
  confirm it is an extrusion with one hole, turn its outlines into polygons
  that never make the hole smaller, run ``mjfit.peg`` and check the result
  against the CAD part.

Both return ``(mjcf_xml, report)``; ``report["status"]`` is ``certified`` or
``refused`` (with ``report["refusal"]``), and ``xml`` is ``None`` when mode C
refuses before mjfit runs. See ``docs/mjfit_bridge.md`` for what each number
means.
"""
from __future__ import annotations

import math

import numpy as np

from .certify import SCHEMA, certify_model, default_resolution
from .geom2d import DifferenceRegion, LoopRegion, OffsetRegion, deviation, polygon_loop
from .reference import (Refusal, cap_levels, classify_loops, extrusion_axis, mjfit_centroid, orient,
                        polygonize, section_loops, spec_loop, to_axis_frame)

__all__ = ["certify_mjfit_peg", "certify_mjfit_config", "mjfit_peg_from_cad", "SCHEMA"]


def _refused(kind, message, **extra):
    return {"schema": SCHEMA, "status": "refused", "refusal": {"kind": kind, "message": message},
            "warnings": [], **extra}


def _require_mjfit():
    try:
        import mjfit
    except ImportError as exc:
        raise ImportError("this needs mjfit: pip install -e .[mjfit]") from exc
    return mjfit


# ---------------------------------------------------------------------------
# mode A
# ---------------------------------------------------------------------------


def _allowed_error(spec, scale):
    """How far mjfit's own polygon may sit from the true curve for this spec."""
    from mjfit import shapes
    if not isinstance(spec, dict):
        return 0.0
    kind = spec.get("type")
    if kind in ("circle", "stadium"):
        r = float(spec["radius"]) if kind == "circle" else min(float(spec["width"]), float(spec["height"])) / 2
        if spec.get("n"):
            n = int(spec["n"])
            n += n % 2 if kind == "stadium" else 0
            return scale * r * (1 / math.cos(math.pi / n) - 1)
        return scale * float(spec.get("max_error", shapes.MAX_ERROR))
    if kind == "bezier":
        return scale * float(spec.get("max_error", shapes.MAX_ERROR))
    return 0.0


def certify_mjfit_peg(sample_spacing: float | None = None, **peg_kwargs):
    """(xml, report) for ``mjfit.peg(**peg_kwargs)`` checked against true outlines."""
    mjfit = _require_mjfit()
    from mjfit import shapes

    xml = mjfit.peg(**peg_kwargs)
    from pathlib import Path
    base = Path(peg_kwargs["mesh_dir"]).parent if peg_kwargs.get("mesh_dir") else None
    scale = float(peg_kwargs.get("scale", 1.0))
    c = float(peg_kwargs.get("clearance", 1e-4)) * scale
    hole, peg = peg_kwargs.get("hole"), peg_kwargs.get("peg")
    source = shapes.from_spec(hole if hole is not None else peg) * scale
    shift = mjfit_centroid(source) if peg_kwargs.get("center", True) else np.zeros(2)
    size = float(np.abs(source - shift).max()) * 2
    specs = [s for s in [hole, peg, *peg_kwargs.get("islands", ())] if s is not None]
    max_err = max([_allowed_error(s, scale) for s in specs] + [0.0])
    h = sample_spacing or default_resolution(c, max_err or c, size)

    def loop(spec):
        lp, kind = spec_loop(spec, scale, shift)
        allowed = _allowed_error(spec, scale)
        own = LoopRegion([polygon_loop(shapes.from_spec(spec) * scale - shift)], h)
        dv = deviation(own, LoopRegion([lp], h))
        worst = max(dv.outside, dv.inside)
        if worst > allowed + 2 * h or dv.grid_mismatch_far:
            raise Refusal("ReferenceMismatch", f"could not rebuild the true {kind} outline: mjfit's "
                          f"polygon is {worst:.3g} m from it (allowed {allowed:.3g} m)")
        return lp, kind

    try:
        kinds = {}
        if hole is not None:
            hl, kinds["hole"] = loop(hole)
            isl = []
            for k, s in enumerate(peg_kwargs.get("islands", ())):
                lp, kinds[f"island_{k}"] = loop(s)
                isl.append(lp)
            free = LoopRegion([hl, *isl], h)
            hole_only = LoopRegion([hl], h)
            peg_ref = OffsetRegion(free, c, "erode")
            back_ref = OffsetRegion(hole_only, c, "erode") if isl else None
        else:
            pl, kinds["peg"] = loop(peg)
            peg_ref = LoopRegion([pl], h)
            free = OffsetRegion(peg_ref, c, "dilate")
            back_ref = None
        if peg_kwargs.get("plate") is not None:
            outer_lp, kinds["plate"] = loop(peg_kwargs["plate"])
        else:
            outer_lp = None
    except Refusal as exc:
        return xml, _refused(exc.kind, str(exc), sample_spacing=h)

    if outer_lp is None:
        # mjfit's square plate: the hull of the plate pieces is its exact outline
        from .geom2d import exact_hull
        from .mjcf import read_pieces
        pts = [(v[0], v[1]) for p in read_pieces(xml, base)["plate"] for v in p.verts]
        outer_lp = polygon_loop(np.array([[float(x), float(y)] for x, y in exact_hull(pts)]))
    outer = LoopRegion([outer_lp], h)
    refs = {"plate_hole": DifferenceRegion(outer, free), "plate_bottom": outer,
            "peg_insert": peg_ref, "peg_back": back_ref, "handle": None}
    n_isl = len(peg_kwargs.get("islands", ()) or ())
    loops = {"plate_hole": (1 + n_isl, 1), "plate_bottom": (1, 0), "peg_insert": (1, n_isl),
             "peg_back": (1, 0), "handle": (1, 0)}
    # with a hole spec, the plate is the outlines' polygons and must match them within
    # their polygon error; with a peg spec the real hole is a design choice (mjfit's
    # sharp-cornered offset or a round one), so its difference is reported, not refused
    budget = None if hole is None else max(
        [_allowed_error(s, scale) for s in [hole, peg_kwargs.get("plate"), *peg_kwargs.get("islands", ())]
         if s is not None] + [0.0])
    rep = certify_model(xml, refs, c, h, base_dir=base, loops=loops, plate_budget=budget)
    rep["mode"] = "mjfit_spec"
    rep["reference"] = {"outlines": kinds, "note": "true curves rebuilt from the mjfit shape specs; "
                        "the ideal peg is the true hole shrunk by the clearance"}
    return xml, rep


def certify_mjfit_config(config, sample_spacing: float | None = None, **overrides):
    """Mode A for an mjfit config (bundled name, YAML path or dict) of kind peg."""
    mjfit = _require_mjfit()
    cfg = {**mjfit.load(config), **overrides}
    kind = cfg.pop("kind", None)
    if kind != "peg":
        return None, _refused("UnsupportedGenerator", f"only peg configs are checked (this is {kind!r}); "
                              "threads are helical and gears need an involute reference")
    return certify_mjfit_peg(sample_spacing=sample_spacing, **cfg)


# ---------------------------------------------------------------------------
# mode C
# ---------------------------------------------------------------------------


def mjfit_peg_from_cad(source, clearance: float, *, max_error: float = 1e-4, scale: float | None = None,
                       sample_spacing: float | None = None, **peg_kwargs):
    """(xml or None, report): an mjfit peg-in-hole model for a CAD plate.

    ``source`` is a path (.step/.stp read as millimeters, .brep/.usda as
    meters, unless ``scale`` is given) or a TopoDS_Shape in meters. Extra
    keyword arguments go to ``mjfit.peg`` (peg_length, handle_height, ...).
    """
    mjfit = _require_mjfit()
    from ..cli import CliError, frame_shape, read_shape

    info: dict = {"mode": "cad"}
    try:
        if isinstance(source, str):
            try:
                shape, fmt = read_shape(source)
            except CliError as exc:
                return None, _refused("ReadFailed", str(exc))
            s = scale if scale is not None else (0.001 if fmt == "step" else 1.0)
            info["input"] = {"path": source, "format": fmt, "scale": s}
            shape = frame_shape(shape, s, False)
        else:
            shape = source
        axis = extrusion_axis(shape)
        from OCP.Bnd import Bnd_Box
        from OCP.BRepBndLib import BRepBndLib
        box = Bnd_Box()
        BRepBndLib.Add_s(shape, box)
        size = box.CornerMin().Distance(box.CornerMax())
        tol = 1e-9 * size
        hc = size / 4000

        def slabs_of(flip):
            sz, m = to_axis_frame(shape, axis, flip)
            lv = cap_levels(sz, tol)
            out = []
            for z0, z1 in zip(lv, lv[1:]):
                loops = section_loops(sz, (z0 + z1) / 2, 1e-7 * size)
                out.append((z0, z1, loops, classify_loops(loops, hc)))
            return sz, m, out

        sz, m, slabs = slabs_of(False)
        if len(slabs) == 2 and not slabs[0][3][1] and slabs[1][3][1]:
            pass
        elif len(slabs) == 2 and slabs[0][3][1] and not slabs[1][3][1]:
            sz, m, slabs = slabs_of(True)          # blind hole opening downward: turn it over
        elif len(slabs) != 1:
            raise Refusal("UnsupportedPlate", f"the plate has {len(slabs)} different cross-sections; "
                          "this version handles a through hole or one blind hole with a flat bottom")
        z0, z1, loops, (io, holes, isl) = slabs[-1]
        if len(holes) != 1:
            raise Refusal("UnsupportedPlate", f"the plate has {len(holes)} holes; mjfit models one hole "
                          "(with optional islands inside it)")
        hole_lp = orient(loops[holes[0]], True, hc)
        outer_lp = orient(loops[io], False, hc)
        isl_lp = [orient(loops[i], False, hc) for i in isl]
        hole_poly = polygonize(hole_lp, max_error, hc)
        outer_poly = polygonize(outer_lp, max_error, hc)[::-1]
        isl_poly = [polygonize(lp, max_error, hc) for lp in isl_lp]
        thick = z1 - z0
        bottom = slabs[0][1] - slabs[0][0] if len(slabs) == 2 else 0.0
        kw = dict(peg_kwargs)
        kw.setdefault("peg_length", max(0.03, 2 * thick))
        info.update({"axis": list(map(float, axis)), "frame": m.tolist(),
                     "plate_thickness": thick, "bottom_thickness": bottom,
                     "polygon_vertices": {"hole": len(hole_poly), "plate": len(outer_poly),
                                          "islands": [len(p) for p in isl_poly]},
                     "max_error": max_error})
        try:
            xml = mjfit.peg(hole=hole_poly, islands=isl_poly, plate=outer_poly, plate_thickness=thick,
                            bottom_thickness=bottom, clearance=clearance, center=False, **kw)
        except ValueError as exc:
            raise Refusal("MjfitRejected", f"mjfit rejected the outline: {exc}") from exc
    except Refusal as exc:
        return None, _refused(exc.kind, str(exc), **info)

    h = sample_spacing or default_resolution(clearance, max_error, size)
    plate_ref = LoopRegion([loops[io], loops[holes[0]], *[loops[i] for i in isl]], h)
    free = LoopRegion([loops[holes[0]], *[loops[i] for i in isl]], h)
    hole_only = LoopRegion([loops[holes[0]]], h)
    bottom_ref = LoopRegion(slabs[0][2], h) if len(slabs) == 2 else None
    try:
        refs = {"plate_hole": plate_ref, "plate_bottom": bottom_ref,
                "peg_insert": OffsetRegion(free, clearance, "erode"),
                "peg_back": OffsetRegion(hole_only, clearance, "erode") if isl else None, "handle": None}
    except ValueError as exc:
        return xml, _refused("ClearanceTooLarge", f"the hole shrunk by the clearance is empty: {exc}", **info)
    loops = {"plate_hole": (1 + len(isl), 1), "plate_bottom": (1, 0), "peg_insert": (1, len(isl)),
             "peg_back": (1, 0), "handle": (1, 0)}
    rep = certify_model(xml, refs, clearance, h, loops=loops, plate_budget=max_error)
    rep.update(info)
    return xml, rep
