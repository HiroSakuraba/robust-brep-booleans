"""NVIDIA SMLib (pip package ``usd-brep``) as an optional second engine.

brepkernel certifies its own results, but every check it runs is built on
OCCT. This module adds an engine that shares no code with OCCT: SMLib, the
NURBS kernel inside NVIDIA's usd-brep (Apache-2.0). Shapes cross between
the two through omniSolid BrepArray .usda files written and read by
brepkernel itself (omnisolid_write / omnisolid_read), so nothing here
needs a C++ build.

Install with ``pip install usd-brep`` (CPython 3.12 wheels) or
``pip install brepkernel[usd]``. Without it, ``available()`` is False and
every entry point raises SmlibUnavailable; the rest of brepkernel,
including USD export and import, works unchanged.

Public entry points:
    available()                          -> bool
    to_smlib(shape)                      -> usd_brep.Brep
    from_smlib(brep)                     -> TopoDS_Shape
    smlib_boolean(shapeA, shapeB, op)    -> TopoDS_Shape
    cross_check(shapeA, shapeB, op, result_shape, ...) -> dict

How the cross-check decides (see cross_check): SMLib computes the same
Boolean, its result is read back into OCCT, and the two results are
compared in one place with one integrator: solid count, volume, and a
deterministic set of point-membership probes. SMLib's own volume
integration is not used for the verdict because it loses accuracy on
B-spline faces (7e-4 relative on a lofted solid whose geometry is
correct; see docs/usd_brep_backend.md).
"""
from __future__ import annotations

import os
import tempfile

from .. import _occt_compat  # noqa: F401
from .omnisolid_read import OmniSolidReadError, read_omnisolid_usda
from .omnisolid_write import OmniSolidWriteError, solids_of, write_omnisolid_usda

ENGINE = "smlib"
_OPS = ("union", "intersection", "difference")


class SmlibUnavailable(RuntimeError):
    """usd-brep is not installed in this interpreter."""


class SmlibError(RuntimeError):
    """SMLib, or the conversion to or from it, did not produce a usable solid."""

    def __init__(self, message: str, kind: str):
        super().__init__(message)
        self.kind = kind


_U = None


def _load():
    global _U
    if _U is None:
        try:
            # OCP first: both packages ship their own shared libraries and
            # this order is the one tested.
            import OCP.TopoDS  # noqa: F401
            import usd_brep
        except ImportError as exc:
            raise SmlibUnavailable(
                "usd-brep is not installed; pip install usd-brep "
                "(or brepkernel[usd]) on CPython 3.12") from exc
        _U = usd_brep
    return _U


def available() -> bool:
    try:
        _load()
        return True
    except SmlibUnavailable:
        return False


def smlib_version() -> str | None:
    try:
        from importlib.metadata import version
        return version("usd-brep")
    except Exception:
        return None


# ---------------------------------------------------------------------------
# conversion
# ---------------------------------------------------------------------------

def to_smlib(shape, *, workdir: str | None = None):
    """OCCT solid(s) -> one SMLib Brep, through an omniSolid .usda file."""
    U = _load()
    with tempfile.TemporaryDirectory(dir=workdir) as tmp:
        path = os.path.join(tmp, "in.usda")
        try:
            prim = write_omnisolid_usda(shape, path)
        except OmniSolidWriteError as exc:
            raise SmlibError(f"cannot write omniSolid: {exc}", "ConversionFailed") from exc
        try:
            brep = U.usd.import_brep(path, prim, heal=False)
        except RuntimeError as exc:
            raise SmlibError(f"SMLib rejected the converted solid: {exc}",
                             "ConversionFailed") from exc
    if not brep.is_manifold_solid():
        raise SmlibError("converted solid is not a closed manifold solid in SMLib",
                         "ConversionFailed")
    census = U.material_census(brep)
    if census.get("solid_count") != len(solids_of(shape)):
        raise SmlibError(f"solid count changed in conversion: OCCT {len(solids_of(shape))}, "
                         f"SMLib {census}", "ConversionFailed")
    return brep


def from_smlib(brep, *, workdir: str | None = None):
    """SMLib Brep -> OCCT shape, through an omniSolid .usda file."""
    U = _load()
    with tempfile.TemporaryDirectory(dir=workdir) as tmp:
        path = os.path.join(tmp, "out.usda")
        try:
            U.usd.export_brep(brep, path)
        except RuntimeError as exc:
            raise SmlibError(f"SMLib could not export its result: {exc}",
                             "ResultUnreadable") from exc
        try:
            shape = read_omnisolid_usda(path)
        except OmniSolidReadError as exc:
            raise SmlibError(f"SMLib result could not be rebuilt in OCCT: {exc}",
                             "ResultUnreadable") from exc
    # Safety net against a converter error passing as a disagreement: the
    # rebuilt solid must match SMLib's own view of its result. The bound is
    # loose on purpose (SMLib's integration was seen 7e-4 off on B-spline
    # faces); it catches a wrong face region, not small deviations.
    census = U.material_census(brep)
    if census.get("solid_count") != len(solids_of(shape)):
        raise SmlibError(f"rebuilt result has {len(solids_of(shape))} solids, SMLib "
                         f"reports {census}", "ResultUnreadable")
    try:
        v_sm = float(U.compute_volume(brep, relative_accuracy=1e-6))
    except RuntimeError:
        v_sm = None
    if v_sm is not None:
        v_oc = _volume(shape)
        if abs(v_oc - v_sm) > 5e-3 * max(abs(v_sm), abs(v_oc), 1e-300):
            raise SmlibError(f"rebuilt result volume {v_oc:.9g} differs from SMLib's "
                             f"own {v_sm:.9g}", "ResultUnreadable")
    return shape


def smlib_boolean(shapeA, shapeB, op: str, *, workdir: str | None = None):
    """The same Boolean computed by SMLib, returned as an OCCT shape.

    No certification is attached: this is the raw second opinion. An
    empty result comes back as an empty compound.
    """
    U = _load()
    if op not in _OPS:
        raise ValueError(f"op must be one of {_OPS}")
    a = to_smlib(shapeA, workdir=workdir)
    b = to_smlib(shapeB, workdir=workdir)
    code = {"union": U.BooleanOp.UNION, "intersection": U.BooleanOp.INTERSECTION,
            "difference": U.BooleanOp.DIFFERENCE}[op]
    try:
        r = U.boolean(a, b, operation=code)      # consumes a and b
    except RuntimeError as exc:
        raise SmlibError(f"SMLib boolean failed: {exc}", "EngineRefused") from exc
    if r is None or r.face_count() == 0:
        return _empty()
    if not r.is_manifold_solid():
        # SMLib sometimes returns a broken body instead of raising (seen on
        # a sphere-cylinder intersection: 6 faces, one void region, no
        # solid). That is no answer, not a different answer.
        raise SmlibError(f"SMLib returned a body that is not a closed solid "
                         f"({r.face_count()} faces, {U.material_census(r)})",
                         "EngineBadResult")
    return from_smlib(r, workdir=workdir)


def _empty():
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Compound
    c = TopoDS_Compound()
    BRep_Builder().MakeCompound(c)
    return c


# ---------------------------------------------------------------------------
# comparison
# ---------------------------------------------------------------------------

def _volume(shape) -> float:
    """brepkernel's own volume routing (assembly._shape_volume), tightened
    to 1e-12 on all-analytic shapes: shapes rebuilt from USD carry fitted
    UV curves, and OCCT's default integration at 1e-9 lost 4e-7 relative
    on such a cylinder (1e-12: 6e-11)."""
    from ..assembly import _all_faces_analytic, _shape_volume
    if shape is None or shape.IsNull() or not solids_of(shape):
        return 0.0
    tol = 1e-12 if _all_faces_analytic(shape) else 1e-10
    return abs(float(_shape_volume(shape, tol)))


def _area(shape) -> float:
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    if shape is None or shape.IsNull() or not solids_of(shape):
        return 0.0
    g = GProp_GProps()
    BRepGProp.SurfaceProperties_s(shape, g, 1e-9, True)
    return abs(float(g.Mass()))


def _bbox(shapes):
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib
    box = Bnd_Box()
    for s in shapes:
        if s is not None and not s.IsNull():
            BRepBndLib.Add_s(s, box, False)
    if box.IsVoid():
        return None
    lo, hi = box.CornerMin(), box.CornerMax()
    return (lo.X(), lo.Y(), lo.Z()), (hi.X(), hi.Y(), hi.Z())


def _classifier(shape, band):
    """point -> 'in' | 'out' | 'on' (within band of the boundary)."""
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.TopAbs import TopAbs_IN, TopAbs_ON, TopAbs_OUT
    from OCP.gp import gp_Pnt
    solids = solids_of(shape) if shape is not None and not shape.IsNull() else []
    clfs = [BRepClass3d_SolidClassifier(s) for s in solids]

    def classify(p):
        on = False
        for c in clfs:
            c.Perform(gp_Pnt(*p), band)
            st = c.State()
            if st == TopAbs_IN:
                return "in"
            if st == TopAbs_ON:
                on = True
            elif st != TopAbs_OUT:
                on = True                       # UNKNOWN: never counted as evidence
        return "on" if on else "out"
    return classify


def _halton(i: int, base: int) -> float:
    f, r = 1.0, 0.0
    while i > 0:
        f /= base
        r += f * (i % base)
        i //= base
    return r


def compare_results(result_shape, other_shape, *, volume_tol: float | None = None,
                    n_points: int = 512, band: float | None = None) -> dict:
    """Compare two OCCT results of the same Boolean. Pure measurement.

    Points are a deterministic Halton set over the union bounding box; a
    point within `band` of either boundary is skipped, so only clear
    in/out contradictions count.
    """
    import numpy as np
    vr, vo = _volume(result_shape), _volume(other_shape)
    sr, so = len(solids_of(result_shape)), len(solids_of(other_shape))
    box = _bbox([result_shape, other_shape])
    if box is None:
        return {"volume": {"brepkernel": vr, "other": vo, "difference": 0.0,
                           "tolerance": 0.0, "ok": True},
                "solids": {"brepkernel": sr, "other": so, "ok": sr == so},
                "points": {"tested": 0, "skipped_near_boundary": 0, "disagreements": 0,
                           "ok": True}}
    lo, hi = np.array(box[0]), np.array(box[1])
    diag = float(np.linalg.norm(hi - lo))
    if band is None:
        band = max(1e-6, 1e-5 * diag)
    if volume_tol is None:
        # a boundary offset of `band` over the larger surface, floored at
        # 1e-7 of the volume scale
        area = max(_area(result_shape), _area(other_shape))
        volume_tol = max(band * area, 1e-7 * max(vr, vo, 1e-30))
    dv = abs(vr - vo)
    cr, co = _classifier(result_shape, band), _classifier(other_shape, band)
    tested = skipped = 0
    bad = []
    pad = 1e-3 * (hi - lo) + 1e-9
    lo2, hi2 = lo - pad, hi + pad
    for i in range(1, n_points + 1):
        p = lo2 + (hi2 - lo2) * np.array([_halton(i, 2), _halton(i, 3), _halton(i, 5)])
        a, b = cr(p), co(p)
        if a == "on" or b == "on":
            skipped += 1
            continue
        tested += 1
        if a != b:
            bad.append({"point": [float(x) for x in p], "brepkernel": a, "other": b})
    return {
        "volume": {"brepkernel": vr, "other": vo, "difference": dv,
                   "tolerance": volume_tol, "ok": dv <= volume_tol},
        "solids": {"brepkernel": sr, "other": so, "ok": sr == so},
        "points": {"tested": tested, "skipped_near_boundary": skipped,
                   "band": band, "disagreements": len(bad), "examples": bad[:5],
                   "ok": not bad},
    }


def cross_check(shapeA, shapeB, op: str, result_shape, *, volume_tol: float | None = None,
                n_points: int = 512, workdir: str | None = None) -> dict:
    """Recompute `op` in SMLib and compare it with brepkernel's result.

    Returns a JSON-ready dict whose "status" is one of
        "agree"        SMLib's result matches on solids, volume and points
        "disagree"     a clear contradiction (see "comparison")
        "unavailable"  usd-brep is not installed
        "inconclusive" SMLib could not give an answer: conversion failed,
                       SMLib refused the Boolean, or its result could not
                       be rebuilt ("kind" says which)
    "inconclusive" is not evidence either way; the caller decides whether
    that is acceptable (pipeline: independent_engine_strict).
    """
    out = {"engine": ENGINE, "engine_version": smlib_version(), "op": op}
    if not available():
        out.update(status="unavailable",
                   message="usd-brep is not installed (pip install usd-brep)")
        return out
    try:
        other = smlib_boolean(shapeA, shapeB, op, workdir=workdir)
    except SmlibError as exc:
        out.update(status="inconclusive", kind=exc.kind, message=str(exc))
        return out
    cmp = compare_results(result_shape, other, volume_tol=volume_tol, n_points=n_points)
    agree = cmp["volume"]["ok"] and cmp["solids"]["ok"] and cmp["points"]["ok"]
    out.update(status="agree" if agree else "disagree", comparison=cmp)
    return out
