"""DEFERRED: curved analytic and NURBS coincidence machinery.

Status: set aside by scope decision, 2026-09-25 (G2 narrowed to planar
coincident faces; see docs/REVIEW_LEDGER.md, G2 closeout). This module
is NOT imported by the active pipeline. It preserves the pre-narrowing
G2.3 coincidence machinery verbatim so a follow-up gate can re-enable
it deliberately instead of re-deriving it.

What was moved here:
  - exact Fraction parameter extraction for cylinder/sphere/cone/torus
    supports (_cylinder_params, _sphere_params, _cone_params,
    _torus_params);
  - the curved branches of the exact coincidence test
    (_exact_coincident for GeomAbs_Cylinder/Sphere/Cone/Torus);
  - the dense two-direction sampling coincidence *acceptance* path for
    NURBS/mixed pairs (the active pipeline now only uses sampling to
    prove "distinct", never to accept).

Gate: ENABLE_DEFERRED_COINCIDENCE defaults to False. The active
pipeline (brepkernel.coincidence.classify_support_pair) never consults
this module; any curved pair that could be coincident is case (c) and
the caller refuses with kind CoincidenceUndecidable.

What a follow-up gate must do before re-enabling:
  1. Recover canonical analytic parameters through
     ShapeAnalysis_CanonicalRecognition (bindings exist on
     cadquery-ocp 8.0.1.0.0) instead of trusting float adaptor
     parameters; hand-fitting from samples is not acceptable.
  2. Apply the same three-case rule as the planar path
     (brepkernel.coincidence): EXACT on canonical parameters,
     TOLERANCE-CERTIFIED within the entities' actual tolerances,
     UNDECIDABLE otherwise, with deviation + tolerances on the record.
  3. Add I4-first tests: coincident cylinder/sphere/cone/torus pairs
     that the old rule accepted, near-coincident curved negatives,
     and the NURBS sampling-acceptance cases with per-entity
     tolerance certification.
  4. Re-run the full suite and the review probes; curved acceptance
     must not regress the planar guarantees.
"""
from __future__ import annotations

from fractions import Fraction

# Gate: default off. Nothing in the active pipeline reads this flag;
# it exists so a follow-up gate has an explicit switch to flip after
# completing the work listed above.
ENABLE_DEFERRED_COINCIDENCE = False


def _fr(v: float) -> Fraction:
    """Exact Fraction of a float (the exact binary value)."""
    return Fraction(v)


def _as_fracs(xs) -> tuple:
    return tuple(_fr(float(x)) for x in xs)


def _adaptor(face):
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    return BRepAdaptor_Surface(face)


def _cylinder_params(face):
    s = _adaptor(face)
    cyl = s.Cylinder()
    ax = cyl.Axis()
    loc = ax.Location()
    d = ax.Direction()
    return (_as_fracs((loc.X(), loc.Y(), loc.Z())),
            _as_fracs((d.X(), d.Y(), d.Z())),
            _fr(float(cyl.Radius())))


def _sphere_params(face):
    s = _adaptor(face)
    sph = s.Sphere()
    loc = sph.Location()
    return (_as_fracs((loc.X(), loc.Y(), loc.Z())),
            _fr(float(sph.Radius())))


def _cone_params(face):
    s = _adaptor(face)
    cone = s.Cone()
    ax = cone.Axis()
    loc = ax.Location()
    d = ax.Direction()
    apex = cone.Apex()
    return (_as_fracs((loc.X(), loc.Y(), loc.Z())),
            _as_fracs((d.X(), d.Y(), d.Z())),
            _as_fracs((apex.X(), apex.Y(), apex.Z())),
            abs(_fr(float(cone.SemiAngle()))))


def _torus_params(face):
    s = _adaptor(face)
    tor = s.Torus()
    ax = tor.Axis()
    loc = tor.Location()
    d = ax.Direction()
    return (_as_fracs((loc.X(), loc.Y(), loc.Z())),
            _as_fracs((d.X(), d.Y(), d.Z())),
            _fr(float(tor.MajorRadius())),
            _fr(float(tor.MinorRadius())))


def _neg(v):
    return tuple(-x for x in v)


def _sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _is_zero(v):
    return all(x == 0 for x in v)


def _parallel_mod_sign(d1, d2):
    """Exact check that two direction tuples are parallel (modulo sign)."""
    return d1 == d2 or d1 == _neg(d2)


def _same_line(l1, d1, l2, d2):
    """Exact check that (l1,d1) and (l2,d2) define the same 3D line."""
    return _parallel_mod_sign(d1, d2) and _is_zero(_cross(_sub(l2, l1), d1))


def _exact_coincident_curved(kind: str, pa, pb) -> bool:
    """Exact Fraction coincidence test for curved analytic parameters.

    Moved verbatim from the pre-narrowing G2.3 implementation. Only
    usable after a follow-up gate completes the parameter-recovery work
    described in the module docstring.
    """
    if kind == "GeomAbs_Cylinder":
        (la, da, ra), (lb, db, rb) = pa, pb
        return ra == rb and _same_line(la, da, lb, db)
    if kind == "GeomAbs_Sphere":
        (ca, ra), (cb, rb) = pa, pb
        return ra == rb and ca == cb
    if kind == "GeomAbs_Cone":
        (la, da, apa, ha), (lb, db, apb, hb) = pa, pb
        return ha == hb and apa == apb and _same_line(la, da, lb, db)
    if kind == "GeomAbs_Torus":
        (la, da, ra, na), (lb, db, rb, nb) = pa, pb
        return ra == rb and na == nb and la == lb \
            and _parallel_mod_sign(da, db)
    return False


_CURVED_PARAM_FNS = {
    "GeomAbs_Cylinder": _cylinder_params,
    "GeomAbs_Sphere": _sphere_params,
    "GeomAbs_Cone": _cone_params,
    "GeomAbs_Torus": _torus_params,
}

# The pre-narrowing analytic set (plane included), for fidelity of the
# reference implementation below.
_OLD_ANALYTIC_KINDS = frozenset(
    {"GeomAbs_Plane"} | set(_CURVED_PARAM_FNS))


def deferred_classify_pair(fa, fb, tol: float):
    """Reference re-implementation of the pre-narrowing curved/NURBS path.

    NOT called by the active pipeline (see module docstring). Provided
    so a follow-up gate can diff its certified re-implementation against
    the old behavior. Returns "distinct", ("coincident", "same"|"opp"),
    or "undecidable", exactly as G2.3 did before the planar narrowing.
    """
    from .coincidence import (coincidence_band, _sample_points,
                             _max_distance_to_support, _recognized_plane,
                             _plane_params_from_pln, _planes_exactly_equal,
                             _sense)
    band = coincidence_band(tol)
    ka, kb = fa.surface_type, fb.surface_type

    pa_fn = _CURVED_PARAM_FNS.get(ka)
    pb_fn = _CURVED_PARAM_FNS.get(kb)
    if pa_fn is not None and pb_fn is not None and ka == kb:
        try:
            pa = pa_fn(fa.face)
            pb = pb_fn(fb.face)
        except Exception:
            return "undecidable"
        if _exact_coincident_curved(ka, pa, pb):
            sense = _sense(fa, fb, band)
            if sense == "undecidable":
                return "undecidable"
            return ("coincident", sense)
        pts = _sample_points(fa.face, fa.uv_bounds, 3)
        pts += _sample_points(fb.face, fb.uv_bounds, 3)
        if not pts:
            return "undecidable"
        d1, ok1 = _max_distance_to_support(pts, fb.face, band)
        d2, ok2 = _max_distance_to_support(pts, fa.face, band)
        if (ok1 and d1 > band) or (ok2 and d2 > band):
            return "distinct"
        return "undecidable"

    pln_a = _recognized_plane(fa.face) \
        if ka not in _OLD_ANALYTIC_KINDS else None
    pln_b = _recognized_plane(fb.face) \
        if kb not in _OLD_ANALYTIC_KINDS else None
    if pln_a is not None and pln_b is not None:
        pa = _plane_params_from_pln(pln_a)
        pb = _plane_params_from_pln(pln_b)
        if _planes_exactly_equal(pa, pb):
            sense = _sense(fa, fb, band)
            if sense == "undecidable":
                return "undecidable"
            return ("coincident", sense)
        return "undecidable"

    pts_a = _sample_points(fa.face, fa.uv_bounds, 7)
    pts_b = _sample_points(fb.face, fb.uv_bounds, 7)
    if not pts_a or not pts_b:
        return "undecidable"
    d_ab, ok_ab = _max_distance_to_support(pts_a, fb.face, band)
    d_ba, ok_ba = _max_distance_to_support(pts_b, fa.face, band)
    if not ok_ab or not ok_ba:
        return "undecidable"
    if d_ab > band or d_ba > band:
        return "distinct"
    sense = _sense(fa, fb, band)
    if sense == "undecidable":
        return "undecidable"
    return ("coincident", sense)
