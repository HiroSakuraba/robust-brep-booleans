"""G2-planar: three-case coincidence classification of face-support pairs.

Scope: PLANAR supports only. For a candidate face pair, decide whether
the two faces lie on the same geometric support surface. Returns one of:
  "distinct"                       the supports are provably different;
  ("coincident", sense, cert)       same support; sense is "same"/"opp";
                                   cert is "exact" or "tolerance_certified";
  CoincidenceUndecidable            cannot tell apart safely; the caller
                                   refuses with a typed kind.

The three cases:
  (a) EXACT: canonical analytic parameters agree bit-identically
      (direction modulo sign, offset exactly zero, decided with
      fractions.Fraction on the float inputs). For a recognized
      (non-analytic) face this means the parameters recovered through
      the same canonical recognizer agree with the other face's
      canonical parameters.
  (b) TOLERANCE-CERTIFIED: the parameters are not bit-identical, but at
      least one face's support parameters carry recovery uncertainty
      (planar recognition of a non-analytic face), and the measured
      support deviation is within the entities' ACTUAL tolerances: each
      face's own OCCT face/edge/vertex tolerances, and the pipeline
      contact_tol. The test uses the tolerance each entity reports for
      itself, never a global constant.
  (c) UNDECIDABLE: anything else. The record carries the measured
      deviation and every tolerance consulted; the caller refuses with
      NearCoincidentFaces (planar, in band but not exact) or
      CoincidenceUndecidable (deferred curved/NURBS, failed sense).

Analytic/analytic planar pairs never take case (b): their parameters
are exact, so a nonzero deviation is a genuine geometric difference,
not measurement noise. Per-entity tolerances make such a pair
indistinguishable, not identical; fusing on indistinguishability would
be a guess, so it is case (c). Near-coincidence inside the band but not
exact is case (c), never silently accepted.

Curved analytic (cylinder/sphere/cone/torus) and NURBS coincidence
machinery is deferred to src/brepkernel/coincidence_deferred.py. Here,
a curved pair is only ever called "distinct" when sampling proves the
supports differ beyond every consulted tolerance; any curved pair that
could be coincident is case (c).

Every tolerance used is recorded ToleranceLedger-style on the returned
record.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from typing import Optional, Union

import math as _math
import numpy as np

# Coincidence band = 4x base tolerance (matches the contact tolerance).
COINCIDENCE_BAND_RATIO = 4.0

# Normal-dot threshold for the same/opp sense decision. For truly
# coincident supports the outward normals at a common point are exactly
# +1 or -1; anything else means the geometry is not clean and the pair
# is undecidable.
_SENSE_DOT_MIN = 0.9

CoincidenceResult = Union[str, tuple, "CoincidenceUndecidable"]


def coincidence_band(base_tol: float) -> float:
    """Return the coincidence band for a base tolerance."""
    return COINCIDENCE_BAND_RATIO * float(base_tol)


@dataclass
class CoincidentPairRecord:
    """Provenance for one coincident face pair (G2.3/2.4)."""
    face_a: int
    face_b: int
    sense: str  # "same" or "opp"
    # ToleranceLedger-style entries: {"name","value","role"}.
    tolerances: list = field(default_factory=list)
    # How coincidence was established: "exact" or "tolerance_certified".
    certification: str = "exact"
    # Filled by the overlap split (G2.4).
    common_area: float = 0.0
    verified: bool = False


@dataclass
class CoincidenceUndecidable:
    """Case (c): the pair cannot be classified safely.

    Carries the measured support deviation and every tolerance that was
    consulted, so the refusal names its evidence instead of guessing.
    """
    face_a: int
    face_b: int
    # "near_coincident_planar" | "recognized_beyond_tolerance" |
    # "curved_deferred" | "nurbs_deferred" | "sense_failed" |
    # "plane_params_failed" | "support_sampling_failed"
    reason: str
    deviation: float
    deviation_kind: str
    # ToleranceLedger-style entries: {"name","value","role"}.
    tolerances: list = field(default_factory=list)


def _fr(v: float) -> Fraction:
    """Exact Fraction of a float (the exact binary value)."""
    return Fraction(v)


def _as_fracs(xs) -> tuple:
    return tuple(_fr(float(x)) for x in xs)


# ---------------------------------------------------------------------------
# Surface parameter extraction (planar only; curved is deferred).
# ---------------------------------------------------------------------------

def _adaptor(face):
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    return BRepAdaptor_Surface(face)


def _plane_params(face):
    """Return (loc, dir) as Fraction tuples for an analytic planar face."""
    s = _adaptor(face)
    pln = s.Plane()
    loc = pln.Location()
    d = pln.Axis().Direction()
    return (_as_fracs((loc.X(), loc.Y(), loc.Z())),
            _as_fracs((d.X(), d.Y(), d.Z())))



def _neg(v):
    return tuple(-x for x in v)


def _sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))



# Representation-rounding bound for plane coincidence. A rigid motion
# (rotation + translation) applied in float64 re-rounds each plane
# location component; the rounding it can introduce is a few ULPs of the
# component magnitude. 64 ULPs of the largest location component is a
# generous bound on what the motion itself can introduce. It is a
# property of the float64 representation, not of the geometry: no
# geometric epsilon is tuned here, and the bound scales with the model
# (it is ~1e-14 at unit scale, ~1e-11 at 1e3 offset).
_REPRESENTATION_ULPS = 64


def _representation_bound(la, lb) -> float:
    """Largest-tolerable plane offset attributable to float64 rounding."""
    m = 0.0
    for v in la + lb:
        c = abs(float(v))
        if c > m:
            m = c
    return _REPRESENTATION_ULPS * _math.ulp(m)


def _planes_equal_up_to_rounding(pa, pb) -> bool:
    """Case (a): canonical plane coincidence up to float64 rounding.

    Directions must agree modulo sign (bit-identical: the same rigid
    motion applied to identical direction doubles yields identical
    results) and the offset of one location from the other's plane must
    be within what float64 rounding of a rigid motion can introduce
    (see _REPRESENTATION_ULPS). Anything larger is genuine geometric
    separation, handled by the caller's band logic.
    """
    (la, da), (lb, db) = pa, pb
    if not (da == db or da == _neg(db)):
        return False
    off = _dot(_sub(lb, la), da)
    if off == 0:
        return True
    return abs(float(off)) <= _representation_bound(la, lb)



# ---------------------------------------------------------------------------
# Per-entity tolerances: what each entity reports for itself.
# ---------------------------------------------------------------------------

def _entity_tolerances(face) -> dict:
    """Face/edge/vertex tolerances reported by the entities themselves."""
    from OCP.BRep import BRep_Tool
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopAbs import TopAbs_EDGE, TopAbs_VERTEX
    from OCP.TopoDS import TopoDS
    f = TopoDS.Face(face)
    ft = float(BRep_Tool.Tolerance_s(f))
    et = ft
    ex = TopExp_Explorer(f, TopAbs_EDGE)
    while ex.More():
        et = max(et, float(BRep_Tool.Tolerance_s(
            TopoDS.Edge(ex.Current()))))
        ex.Next()
    vt = ft
    ex = TopExp_Explorer(f, TopAbs_VERTEX)
    while ex.More():
        vt = max(vt, float(BRep_Tool.Tolerance_s(
            TopoDS.Vertex(ex.Current()))))
        ex.Next()
    return {"face": ft, "edge_max": et, "vertex_max": vt}


def _tolerance_entries(tag: str, tols: dict, role: str) -> list:
    return [
        {"name": f"{k}_tolerance_{tag}", "value": v,
         "role": f"{role}: {k} tolerance reported by the entity itself"}
        for k, v in tols.items()
    ]


# ---------------------------------------------------------------------------
# Outward normals and the same/opp sense.
# ---------------------------------------------------------------------------

def _outward_normal(face, orientation: str, uv) -> Optional[np.ndarray]:
    """Outward unit normal at a surface (u, v), or None if unavailable."""
    try:
        from OCP.gp import gp_Pnt, gp_Vec
        s = _adaptor(face)
        p = gp_Pnt()
        du = gp_Vec()
        dv = gp_Vec()
        s.D1(float(uv[0]), float(uv[1]), p, du, dv)
        n = du.Crossed(dv)
        m = n.Magnitude()
        if not np.isfinite(m) or m <= 0.0:
            return None
        nv = np.array([n.X() / m, n.Y() / m, n.Z() / m])
        if orientation == "TopAbs_REVERSED":
            nv = -nv
        elif orientation != "TopAbs_FORWARD":
            return None
        return nv
    except Exception:
        return None


def _interior_uv(face, uv_bounds) -> Optional[tuple]:
    """A (u, v) strictly inside the face trim, or None."""
    try:
        from OCP.BRepClass3d import BRepClass3d_SolidExplorer
        from OCP.gp import gp_Pnt
        from OCP.TopoDS import TopoDS
        f = TopoDS.Face(face)
        explorer = BRepClass3d_SolidExplorer()
        p = gp_Pnt()
        if explorer.FindAPointInTheFace_s(f, p):
            # Project back to (u, v) on this face's support.
            from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf
            from OCP.BRepAdaptor import BRepAdaptor_Surface
            proj = GeomAPI_ProjectPointOnSurf(
                p, BRepAdaptor_Surface(f).Surface().Surface())
            if proj.NbPoints() > 0:
                up, vp = proj.LowerDistanceParameters()
                return (float(up), float(vp))
    except Exception:
        pass
    # Fallback: center of the UV bounds.
    try:
        u0, u1, v0, v1 = (float(x) for x in uv_bounds)
        if np.isfinite(u0) and np.isfinite(u1) \
                and np.isfinite(v0) and np.isfinite(v1):
            return ((u0 + u1) / 2.0, (v0 + v1) / 2.0)
    except Exception:
        pass
    return None


def _sense(fa, fb, band: float) -> str:
    """Return "same", "opp", or "undecidable" for coincident supports."""
    from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.TopoDS import TopoDS
    try:
        uv_a = _interior_uv(TopoDS.Face(fa.face), fa.uv_bounds)
        if uv_a is None:
            return "undecidable"
        na = _outward_normal(TopoDS.Face(fa.face), fa.orientation, uv_a)
        if na is None:
            return "undecidable"
        # 3D point on fa's support.
        s_a = _adaptor(TopoDS.Face(fa.face))
        pa = s_a.Value(uv_a[0], uv_a[1])
        # Corresponding point on fb's support (distance ~0 if coincident).
        s_b = _adaptor(TopoDS.Face(fb.face))
        proj = GeomAPI_ProjectPointOnSurf(pa, s_b.Surface().Surface())
        if proj.NbPoints() == 0 or proj.LowerDistance() > band:
            return "undecidable"
        ub, vb = proj.LowerDistanceParameters()
        nb = _outward_normal(TopoDS.Face(fb.face), fb.orientation,
                             (float(ub), float(vb)))
        if nb is None:
            return "undecidable"
        dot = float(np.dot(na, nb))
        if dot > _SENSE_DOT_MIN:
            return "same"
        if dot < -_SENSE_DOT_MIN:
            return "opp"
        return "undecidable"
    except Exception:
        return "undecidable"


# ---------------------------------------------------------------------------
# Sampling support for the conservative curved distinctness check.
# ---------------------------------------------------------------------------

def _sample_points(face, uv_bounds, n: int) -> list:
    """Deterministic n x n grid of 3D points on the face support."""
    from OCP.TopoDS import TopoDS
    pts = []
    try:
        s = _adaptor(TopoDS.Face(face))
        u0, u1, v0, v1 = (float(x) for x in uv_bounds)
        if not all(np.isfinite(x) for x in (u0, u1, v0, v1)):
            return pts
        for i in range(n):
            u = u0 + (u1 - u0) * (i + 0.5) / n
            for j in range(n):
                v = v0 + (v1 - v0) * (j + 0.5) / n
                try:
                    p = s.Value(u, v)
                    pts.append((p.X(), p.Y(), p.Z(), u, v))
                except Exception:
                    continue
    except Exception:
        pass
    return pts


def _max_distance_to_support(pts_xyz, face, band: float) -> tuple:
    """Max over pts of distance to the face's support surface.

    Returns (max_dist, ok) where ok=False means a projection failed
    (inconclusive). Short-circuits above the band.
    """
    from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf
    from OCP.gp import gp_Pnt
    from OCP.TopoDS import TopoDS
    try:
        surf = _adaptor(TopoDS.Face(face)).Surface().Surface()
    except Exception:
        return (float("inf"), False)
    m = 0.0
    for (x, y, z, _u, _v) in pts_xyz:
        try:
            proj = GeomAPI_ProjectPointOnSurf(gp_Pnt(x, y, z), surf)
            if proj.NbPoints() == 0:
                return (m, False)
            d = float(proj.LowerDistance())
        except Exception:
            return (m, False)
        if d > m:
            m = d
            if m > band:
                return (m, True)
    return (m, True)


def _recognized_plane(face):
    """Return the gp_Pln if a face's support is recognized as planar.

    Uses GeomLib_IsPlanarSurface and returns the fitted plane itself;
    callers must NOT use BRepAdaptor_Surface(face).Plane() here, which
    raises on non-analytic surfaces even when the support is planar.
    Returns None when the support is not planar or recognition fails.
    """
    try:
        from OCP.GeomLib import GeomLib_IsPlanarSurface
        from OCP.BRepAdaptor import BRepAdaptor_Surface
        from OCP.TopoDS import TopoDS
        surf = BRepAdaptor_Surface(TopoDS.Face(face)).Surface().Surface()
        checker = GeomLib_IsPlanarSurface(surf)
        if checker.IsPlanar():
            return checker.Plan()
    except Exception:
        pass
    return None


def _plane_params_from_pln(pln) -> tuple:
    """Canonical exact (loc, dir) Fraction params from a gp_Pln."""
    loc = pln.Location()
    d = pln.Axis().Direction()
    return (_as_fracs((loc.X(), loc.Y(), loc.Z())),
            _as_fracs((d.X(), d.Y(), d.Z())))



# ---------------------------------------------------------------------------
# Main entry point.
# ---------------------------------------------------------------------------

# Curved analytic kinds: the coincidence *machinery* for these is deferred
# (coincidence_deferred.py). They are named here only so the "different
# kinds cannot share a support" short-circuit stays provably correct.
_CURVED_ANALYTIC_KINDS = frozenset({
    "GeomAbs_Cylinder", "GeomAbs_Sphere", "GeomAbs_Cone", "GeomAbs_Torus",
})


def _base_entries(tol: float, band: float, contact_tol: float,
                  extra: list) -> list:
    """Shared ToleranceLedger entries for a pair decision."""
    entries = [
        {"name": "base_tol", "value": float(tol),
         "role": "input base tolerance"},
        {"name": "coincidence_band", "value": band,
         "role": "coincident support gate (4x base_tol)"},
        {"name": "contact_tol", "value": float(contact_tol),
         "role": "pipeline contact threshold consulted for certification"},
    ]
    entries.extend(extra)
    return entries


def classify_support_pair(fa, fb, tol: float,
                          contact_tol: Optional[float] = None
                          ) -> CoincidenceResult:
    """Classify the support-surface pair of two FaceRecords.

    Planar scope only. Returns "distinct", ("coincident", sense,
    certification), or a CoincidenceUndecidable record. `tol` is the
    base tolerance; the coincidence band is 4x tol. `contact_tol`
    defaults to 4x tol (the pipeline default).
    """
    band = coincidence_band(tol)
    ctol = float(contact_tol) if contact_tol is not None else band
    ka, kb = fa.surface_type, fb.surface_type

    # Different analytic kinds cannot share a support (a plane is never
    # a cylinder, etc.). Provably correct, no sampling needed.
    analytic_a = ka == "GeomAbs_Plane" or ka in _CURVED_ANALYTIC_KINDS
    analytic_b = kb == "GeomAbs_Plane" or kb in _CURVED_ANALYTIC_KINDS
    if analytic_a and analytic_b and ka != kb:
        return "distinct"

    planar_a = (ka == "GeomAbs_Plane")
    planar_b = (kb == "GeomAbs_Plane")

    # ---- Planar / planar --------------------------------------------
    # Recognized-planar faces go through the same canonical recognizer;
    # analytic planes use their exact parameters.
    pln_a = _recognized_plane(fa.face) if not planar_a else None
    pln_b = _recognized_plane(fb.face) if not planar_b else None
    if (planar_a or pln_a is not None) and (planar_b or pln_b is not None):
        # Canonical exact params (case a): analytic params for analytic
        # faces, recognizer-recovered params otherwise.
        try:
            pa = _plane_params(fa.face) if planar_a \
                else _plane_params_from_pln(pln_a)
            pb = _plane_params(fb.face) if planar_b \
                else _plane_params_from_pln(pln_b)
        except Exception:
            tols_a = _entity_tolerances(fa.face)
            tols_b = _entity_tolerances(fb.face)
            return CoincidenceUndecidable(
                fa.face_id, fb.face_id, "plane_params_failed",
                float("nan"), "plane_offset", _base_entries(
                    tol, band, ctol,
                    _tolerance_entries("a", tols_a, "face A")
                    + _tolerance_entries("b", tols_b, "face B")))

        tols_a = _entity_tolerances(fa.face)
        tols_b = _entity_tolerances(fb.face)
        t_entity = max(max(tols_a.values()), max(tols_b.values()))
        entries = _base_entries(
            tol, band, ctol,
            _tolerance_entries("a", tols_a, "face A")
            + _tolerance_entries("b", tols_b, "face B"))
        recognized = (pln_a is not None) or (pln_b is not None)
        outer = max(band, ctol, t_entity)

        if _planes_equal_up_to_rounding(pa, pb):
            sense = _sense(fa, fb, band)
            if sense == "undecidable":
                return CoincidenceUndecidable(
                    fa.face_id, fb.face_id, "sense_failed", 0.0,
                    "plane_offset", entries)
            # Case (a): canonical params agree up to float64 rounding.
            # This includes recognizer-recovered params, which come through
            # the same canonical recognizer on both sides.
            return ("coincident", sense, "exact")

        if not recognized:
            # Both faces analytic: exact geometry decides, no
            # sampling. Non-parallel planes are provably different
            # supports; parallel planes sit apart by the exact offset.
            (la, da), (lb, db) = pa, pb
            if not (da == db or da == _neg(db)):
                return "distinct"
            dev = float(abs(_dot(_sub(lb, la), da)))
            if dev > outer:
                return "distinct"
            # Case (c): parallel, in the band, but not exact. Never
            # silently accepted: analytic params are exact, so this is
            # genuine near-coincidence, not measurement noise.
            return CoincidenceUndecidable(
                fa.face_id, fb.face_id, "near_coincident_planar",
                dev, "plane_offset", entries)

        # A recognizer was involved, so the recovered params are
        # approximate: measure the support deviation by sampling.
        pts_a = _sample_points(fa.face, fa.uv_bounds, 7)
        pts_b = _sample_points(fb.face, fb.uv_bounds, 7)
        d_ab, ok_ab = _max_distance_to_support(
            pts_a, fb.face, outer) if pts_a else (0.0, False)
        d_ba, ok_ba = _max_distance_to_support(
            pts_b, fa.face, outer) if pts_b else (0.0, False)
        if d_ab > outer or d_ba > outer:
            return "distinct"
        if not (ok_ab and ok_ba):
            return CoincidenceUndecidable(
                fa.face_id, fb.face_id, "support_sampling_failed",
                float("nan"), "sampled_support_distance", entries)
        dev = max(d_ab, d_ba)

        # Case (b), tolerance-certified: the measured deviation is
        # within what the entities themselves claim AND the pipeline
        # contact threshold. Never a global constant.
        if dev <= t_entity and dev <= ctol:
            sense = _sense(fa, fb, band)
            if sense == "undecidable":
                return CoincidenceUndecidable(
                    fa.face_id, fb.face_id, "sense_failed", dev,
                    "sampled_support_distance", entries)
            return ("coincident", sense, "tolerance_certified")

        # Case (c): in the band but not certifiable.
        return CoincidenceUndecidable(
            fa.face_id, fb.face_id, "recognized_beyond_tolerance",
            dev, "sampled_support_distance", entries)

    # ---- Deferred: curved analytic or unrecognized NURBS ------------
    # Coincidence acceptance for these is deferred to a follow-up gate.
    # We only ever claim "distinct", and only when sampling proves the
    # supports differ beyond the outer band; anything that could be
    # coincident is case (c).
    nurbs_involved = ("BSpline" in ka) or ("BSpline" in kb) \
        or ("Bezier" in ka) or ("Bezier" in kb)
    n = 7 if nurbs_involved else 3
    pts_a = _sample_points(fa.face, fa.uv_bounds, n)
    pts_b = _sample_points(fb.face, fb.uv_bounds, n)
    tols_a = _entity_tolerances(fa.face)
    tols_b = _entity_tolerances(fb.face)
    entries = _base_entries(
        tol, band, ctol,
        _tolerance_entries("a", tols_a, "face A")
        + _tolerance_entries("b", tols_b, "face B"))
    dev = float("inf")
    if pts_a and pts_b:
        d_ab, ok_ab = _max_distance_to_support(pts_a, fb.face, band)
        d_ba, ok_ba = _max_distance_to_support(pts_b, fa.face, band)
        dev = max(d_ab, d_ba)
        if ok_ab and ok_ba and d_ab > band and d_ba > band:
            return "distinct"
    return CoincidenceUndecidable(
        fa.face_id, fb.face_id,
        "nurbs_deferred" if nurbs_involved else "curved_deferred",
        dev, "sampled_support_distance", entries)
