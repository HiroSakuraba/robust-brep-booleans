"""Trim-aware surface/surface intersection workset for Tier B/C.

This layer consumes the conservative candidate pairs from step_ingest.py and
runs expensive OCCT section work only where faces can actually interact.

For every section edge it requires:
- a 3D edge curve;
- p-curves on both originating faces;
- SameParameter correspondence (repaired on an edge copy if needed);
- adaptive samples whose 3D point agrees with both surface evaluations;
- all sampled UVs classified IN/ON the actual trimmed faces.

A no-edge result is never automatically treated as disjoint: point contacts
and near-zero face/face distance are reported as contact/ambiguous-contact.
That prevents tangent or almost-tangent freeform geometry from silently
falling through a broad-phase optimization.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Optional

import numpy as np

from .freeform import FreeformError
from .step_ingest import BRepModel, FaceRecord, candidate_face_pairs


class IntersectionError(FreeformError):
    pass


def _p3(p) -> np.ndarray:
    return np.array([p.X(), p.Y(), p.Z()], dtype=np.float64)


def _p2(p) -> np.ndarray:
    return np.array([p.X(), p.Y()], dtype=np.float64)


def _point_segment_distance(p: np.ndarray, a: np.ndarray,
                            b: np.ndarray) -> float:
    ab = b - a
    d = float(ab @ ab)
    if d <= 1e-300:
        return float(np.linalg.norm(p - a))
    t = float(np.clip(((p - a) @ ab) / d, 0.0, 1.0))
    return float(np.linalg.norm(p - (a + t * ab)))


def _adaptive_edge_samples(edge, chord_tol: float, *,
                           max_depth: int = 12
                           ) -> tuple[np.ndarray, np.ndarray]:
    """Hybrid deflection sampler with independent local verification.

    OCCT first proposes a non-uniform parameter distribution using its
    QuasiUniformDeflection algorithm. That avoids the dyadic over-refinement
    pattern seen on smoothly parameterized NURBS sections.

    We do *not* trust the proposal blindly. Every proposed interval is checked
    again at quarter/mid/three-quarter parameters with the same chord-deviation
    criterion used by the previous recursive sampler; failing intervals are
    recursively bisected. Global quartile parameters are always included, so
    straight/simple curves still receive at least five verification samples.

    Like the previous implementation, this controls verification density and
    is not claimed as a certified Hausdorff bound.
    """
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    from OCP.GCPnts import GCPnts_QuasiUniformDeflection

    if not chord_tol > 0:
        raise ValueError("chord_tol must be positive")
    c = BRepAdaptor_Curve(edge)
    t0 = float(c.FirstParameter())
    t1 = float(c.LastParameter())
    if not (np.isfinite(t0) and np.isfinite(t1) and t1 > t0):
        raise IntersectionError(
            f"invalid section-edge parameter range [{t0}, {t1}]",
            kind="InvalidSectionEdge")

    cache: dict[float, np.ndarray] = {}

    def point(t: float) -> np.ndarray:
        t = float(t)
        if t not in cache:
            cache[t] = _p3(c.Value(t))
        return cache[t]

    # Mandatory global coverage protects against a sampler returning only
    # endpoints on a geometrically simple curve and makes p-curve/surface
    # agreement checks less dependent on the deflection sampler's placement.
    # 33 points is intentionally modest: curvature-driven OCCT proposals and
    # recursive chord checks may add more where geometry requires it.
    seeds = [t0 + (t1 - t0) * (i / 32.0) for i in range(33)]

    try:
        proposed = GCPnts_QuasiUniformDeflection(
            c, float(chord_tol), t0, t1)
        if proposed.IsDone():
            n = int(proposed.NbPoints())
            if 2 <= n <= 100000:
                seeds.extend(float(proposed.Parameter(i))
                             for i in range(1, n + 1))
    except Exception:
        # Optimization only. The independently checked quartile seed set
        # remains sufficient to fall back to recursive refinement.
        pass

    # Collapse parameter duplicates introduced by combining mandatory and OCCT
    # seed points. Use a scale-aware parameter tolerance rather than exact
    # equality because periodic/trimmed curves can return numerically adjacent
    # endpoints.
    pscale = max(abs(t0), abs(t1), abs(t1 - t0), 1.0)
    ptol = 32.0 * np.finfo(np.float64).eps * pscale
    raw = sorted(t for t in seeds
                 if np.isfinite(t) and t0 - ptol <= t <= t1 + ptol)
    seed_params = []
    for t in raw:
        t = min(max(float(t), t0), t1)
        if not seed_params or abs(t - seed_params[-1]) > ptol:
            seed_params.append(t)
    if not seed_params or abs(seed_params[0] - t0) > ptol:
        seed_params.insert(0, t0)
    else:
        seed_params[0] = t0
    if abs(seed_params[-1] - t1) > ptol:
        seed_params.append(t1)
    else:
        seed_params[-1] = t1

    leaves: list[tuple[float, float]] = []

    def refine(a: float, b: float, depth: int):
        if b - a <= ptol:
            leaves.append((a, b))
            return
        pa, pb = point(a), point(b)
        probes = [a + (b - a) * q for q in (0.25, 0.5, 0.75)]
        dev = max(_point_segment_distance(point(t), pa, pb)
                  for t in probes)
        if dev > chord_tol and depth < max_depth:
            m = 0.5 * (a + b)
            refine(a, m, depth + 1)
            refine(m, b, depth + 1)
        else:
            leaves.append((a, b))

    for a, b in zip(seed_params[:-1], seed_params[1:]):
        refine(float(a), float(b), 0)

    ts = sorted(set([leaves[0][0]] + [b for _, b in leaves]))
    T = np.asarray(ts, dtype=np.float64)
    P = np.vstack([point(float(t)) for t in T])
    return T, P


def _surface_d1(adaptor, u: float, v: float
                ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    from OCP.gp import gp_Pnt, gp_Vec

    p = gp_Pnt()
    du = gp_Vec()
    dv = gp_Vec()
    adaptor.D1(float(u), float(v), p, du, dv)
    return (_p3(p),
            np.array([du.X(), du.Y(), du.Z()], dtype=np.float64),
            np.array([dv.X(), dv.Y(), dv.Z()], dtype=np.float64))


@dataclass
class SectionEdgeRecord:
    edge_index: int
    face_a: int
    face_b: int
    edge: object
    first: float
    last: float
    edge_tolerance: float
    verify_tolerance: float
    parameters: np.ndarray
    xyz: np.ndarray
    uv_a: np.ndarray
    uv_b: np.ndarray
    max_surface_error_a: float
    max_surface_error_b: float
    max_cross_surface_error: float
    trim_ok: bool
    min_transversality: float
    max_transversality: float
    risk_flags: tuple[str, ...] = ()
    repaired_same_parameter: bool = False
    exact_curve_on_surface_checked: bool = False
    exact_surface_error_a: Optional[float] = None
    exact_surface_error_b: Optional[float] = None
    shadow_crosschecked: bool = False
    shadow_max_distance: Optional[float] = None
    shadow_length_rel_error: Optional[float] = None
    # G2.5: True when every sample of this edge classifies TopAbs_ON on
    # the trim boundary of BOTH faces. Such an edge is a boundary
    # contact: it creates no split.
    is_boundary_contact: bool = False


@dataclass
class FaceIntersectionResult:
    face_a: int
    face_b: int
    status: str
    edges: list[SectionEdgeRecord] = field(default_factory=list)
    point_contacts: list[np.ndarray] = field(default_factory=list)
    min_distance: Optional[float] = None
    section_done: bool = True
    raw_curve_count: int = 0
    raw_trimmed_components: int = 0
    raw_unmatched_components: int = 0
    completeness_max_distance: float = 0.0
    completeness_components: list = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


@dataclass
class ModelIntersectionResult:
    pairs: list[FaceIntersectionResult]
    candidate_pairs: int
    section_calls: int
    verified_edges: int
    point_contacts: int
    ambiguous_contacts: int
    skipped_by_broadphase: int
    # G2.3: coincident face pairs (skipped the section; handled by the
    # overlap split). Each entry is a CoincidentPairRecord.
    coincident_pairs: list = field(default_factory=list)
    # G2.5: boundary-contact section edges (ON both trim boundaries).
    # They create no split; they resolve (or fail to resolve) contacts.
    boundary_edges: list = field(default_factory=list)
    shadow_section_calls: int = 0
    shadow_verified_edges: int = 0
    max_shadow_distance: float = 0.0
    completeness_probes: int = 0
    raw_curve_count: int = 0
    raw_trimmed_components: int = 0
    raw_unmatched_components: int = 0
    completeness_max_distance: float = 0.0
    completeness_components: list = field(default_factory=list)

    @property
    def has_ambiguous_contact(self) -> bool:
        return self.ambiguous_contacts > 0


def _repair_same_parameter(edge, tol: float):
    """Return (edge, repaired). Does not mutate the section result."""
    from OCP.BRep import BRep_Tool
    from OCP.BRepLib import BRepLib

    if BRep_Tool.SameParameter_s(edge):
        return edge, False
    new_edge = BRepLib.SameParameter_s(
        edge, float(tol), float(tol), False)
    if not BRep_Tool.SameParameter_s(new_edge):
        raise IntersectionError(
            "section edge could not be made SameParameter",
            kind="SectionNotSameParameter")
    return new_edge, True


def _exact_curve_on_surface_distance(edge, face) -> float:
    """Use OCCT's exact SameParameter curve-on-surface distance validator."""
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    from OCP.BRepLib import BRepLib_ValidateEdge

    try:
        c3 = BRepAdaptor_Curve(edge)
        cf = BRepAdaptor_Curve(edge, face)
        cos = cf.CurveOnSurface()
        v = BRepLib_ValidateEdge(c3, cos, True)
        v.SetExactMethod(True)
        v.SetParallel(False)
        v.Process()
        if not v.IsDone():
            raise IntersectionError(
                "OCCT exact curve-on-surface validation did not complete",
                kind="ExactCurveOnSurfaceValidationFailed")
        d = float(v.GetMaxDistance())
        if not np.isfinite(d) or d < 0.0:
            raise IntersectionError(
                f"invalid exact curve-on-surface distance {d}",
                kind="ExactCurveOnSurfaceValidationFailed")
        return d
    except IntersectionError:
        raise
    except Exception as exc:
        raise IntersectionError(
            f"OCCT exact curve-on-surface validation failed: "
            f"{type(exc).__name__}: {exc}",
            kind="ExactCurveOnSurfaceValidationFailed") from exc


def _section_tolerance_ceiling(fa: FaceRecord, fb: FaceRecord, *,
                               base_tol: float,
                               max_section_tol: Optional[float]) -> float:
    """Maximum acceptable tolerance for one section curve of this pair.

    Shared by _verify_section_edge (SectionToleranceTooLoose) and the
    completeness probe (RawIntersectionToleranceTooLoose) so the two
    refusals cannot drift apart (G3 item 1). The default ceiling is
    max(128 * base_tol, 1e-8 * model scale): 128x the base tolerance is
    the same headroom _verify_section_edge historically allowed, and the
    1e-8 * scale term keeps the ceiling meaningful on large models where
    128 * base_tol alone would be tighter than OCCT's own confusion.
    An explicit max_section_tol overrides the default; it must be
    positive. All inputs are OCCT-reported or existing pipeline
    quantities (I3).
    """
    scale = max(
        float(np.linalg.norm(fa.bbox_hi - fa.bbox_lo)),
        float(np.linalg.norm(fb.bbox_hi - fb.bbox_lo)),
        1.0)
    tol_limit = (max(128.0 * float(base_tol), 1e-8 * scale)
                 if max_section_tol is None
                 else float(max_section_tol))
    if not tol_limit > 0.0:
        raise ValueError("max_section_tol must be positive")
    return tol_limit


def _raw_curve_tolerance(ic) -> float:
    """OCCT-reported tolerance of one IntTools raw intersection curve.

    Separate helper (rather than an inline ic.Tolerance() call) so tests
    can inflate it artificially and check the RawIntersectionToleranceTooLoose
    refusal (G3 item 1, I4).
    """
    return float(ic.Tolerance())


def _verify_section_edge(edge, fa: FaceRecord, fb: FaceRecord,
                         edge_index: int, *,
                         base_tol: float,
                         chord_tol: Optional[float],
                         tangent_sin_tol: float,
                         max_section_tol: Optional[float]
                         ) -> SectionEdgeRecord:
    from OCP.BRep import BRep_Tool
    from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
    from OCP.BRepClass import BRepClass_FaceClassifier
    from OCP.TopAbs import TopAbs_IN, TopAbs_ON
    from OCP.gp import gp_Pnt2d

    et = float(BRep_Tool.Tolerance_s(edge))
    ft = max(float(BRep_Tool.Tolerance_s(fa.face)),
             float(BRep_Tool.Tolerance_s(fb.face)))

    # BRepAlgoAPI_Section may approximate the intersection curve.  Its edge
    # tolerance is useful evidence, but the verifier must not simply widen
    # itself to whatever tolerance the approximation produced.  Cap the
    # acceptable B-rep tolerance relative to the requested kernel tolerance
    # and local face scale; anything looser is a typed refusal.
    tol_limit = _section_tolerance_ceiling(
        fa, fb, base_tol=base_tol, max_section_tol=max_section_tol)
    if et > tol_limit or ft > tol_limit:
        raise IntersectionError(
            f"section edge {edge_index}: OCCT edge/face tolerance "
            f"{max(et, ft):.6g} exceeds acceptance ceiling "
            f"{tol_limit:.6g}",
            kind="SectionToleranceTooLoose")

    verify_tol = max(float(base_tol), 2.0 * et, 2.0 * ft)
    edge, repaired = _repair_same_parameter(edge, verify_tol)

    # SameParameter repair may legitimately update the edge tolerance.  Apply
    # the same ceiling after repair so the repair cannot silently loosen the
    # acceptance contract.
    et = float(BRep_Tool.Tolerance_s(edge))
    if et > tol_limit:
        raise IntersectionError(
            f"section edge {edge_index}: SameParameter repair tolerance "
            f"{et:.6g} exceeds acceptance ceiling {tol_limit:.6g}",
            kind="SectionToleranceTooLoose")
    verify_tol = max(float(base_tol), 2.0 * et, 2.0 * ft)
    pc_a = BRep_Tool.CurveOnSurface_s(edge, fa.face, 0.0, 0.0)
    pc_b = BRep_Tool.CurveOnSurface_s(edge, fb.face, 0.0, 0.0)
    if pc_a is None or pc_b is None:
        raise IntersectionError(
            f"section edge {edge_index} lacks a p-curve on an input face",
            kind="MissingPCurve")

    # This is stronger than our adaptive sample check: OCCT's exact
    # ValidateEdge mode computes the maximum distance between the 3D edge
    # curve and its SameParameter curve-on-surface representation.  Require it
    # independently on both authoritative input faces.
    exact_a = _exact_curve_on_surface_distance(edge, fa.face)
    exact_b = _exact_curve_on_surface_distance(edge, fb.face)
    exact_worst = max(exact_a, exact_b)
    if exact_worst > verify_tol:
        raise IntersectionError(
            f"section edge {edge_index}: exact curve-on-surface distance "
            f"{exact_worst:.6g} exceeds verification tolerance "
            f"{verify_tol:.6g}",
            kind="ExactCurveOnSurfaceMismatch")

    c3 = BRepAdaptor_Curve(edge)
    first = float(c3.FirstParameter())
    last = float(c3.LastParameter())
    # For the default path, verification sampling should be commensurate
    # with the geometric tolerance we are actually willing to accept.  Using
    # base_tol here used to force thousands of samples even when OCCT had
    # already bounded the section at ~1e-5.  A 2x verify-tolerance chord target
    # remains stricter than the acceptance ceiling while avoiding that
    # accidental over-sampling. Explicit caller chord_tol is never relaxed.
    local_chord = (max(2.0 * verify_tol, 1e-9)
                   if chord_tol is None else float(chord_tol))
    ts, xyz = _adaptive_edge_samples(edge, local_chord)

    sa = BRepAdaptor_Surface(fa.face)
    sb = BRepAdaptor_Surface(fb.face)
    uva: list[np.ndarray] = []
    uvb: list[np.ndarray] = []
    err_a: list[float] = []
    err_b: list[float] = []
    err_cross: list[float] = []
    trans: list[float] = []
    trim_ok = True

    for t, p in zip(ts, xyz):
        qa = pc_a.Value(float(t))
        qb = pc_b.Value(float(t))
        ua, ub = _p2(qa), _p2(qb)
        pa, sua, sva = _surface_d1(sa, *ua)
        pb, sub, svb = _surface_d1(sb, *ub)

        err_a.append(float(np.linalg.norm(p - pa)))
        err_b.append(float(np.linalg.norm(p - pb)))
        err_cross.append(float(np.linalg.norm(pa - pb)))

        na = np.cross(sua, sva)
        nb = np.cross(sub, svb)
        den = float(np.linalg.norm(na) * np.linalg.norm(nb))
        sinang = (0.0 if den <= 1e-300 else
                  float(np.linalg.norm(np.cross(na, nb)) / den))
        trans.append(sinang)

        ca = BRepClass_FaceClassifier(
            fa.face, gp_Pnt2d(float(ua[0]), float(ua[1])),
            verify_tol, True)
        cb = BRepClass_FaceClassifier(
            fb.face, gp_Pnt2d(float(ub[0]), float(ub[1])),
            verify_tol, True)
        trim_ok = trim_ok and (
            ca.State() in (TopAbs_IN, TopAbs_ON)
            and cb.State() in (TopAbs_IN, TopAbs_ON))
        uva.append(ua)
        uvb.append(ub)

    max_a = max(err_a, default=0.0)
    max_b = max(err_b, default=0.0)
    max_cross = max(err_cross, default=0.0)
    worst = max(max_a, max_b, max_cross)
    if worst > verify_tol:
        raise IntersectionError(
            f"section edge {edge_index}: 3D/p-curve surface mismatch "
            f"{worst:.6g} exceeds verification tolerance "
            f"{verify_tol:.6g}",
            kind="SectionGeometryMismatch")
    if not trim_ok:
        raise IntersectionError(
            f"section edge {edge_index}: p-curve left a trimmed input face",
            kind="SectionOutsideTrim")

    risk: list[str] = []
    min_tr = min(trans, default=1.0)
    max_tr = max(trans, default=1.0)
    if min_tr < tangent_sin_tol:
        risk.append("near_tangent")
    if BRep_Tool.IsClosed_s(edge, fa.face):
        risk.append("seam_on_a")
    if BRep_Tool.IsClosed_s(edge, fb.face):
        risk.append("seam_on_b")

    return SectionEdgeRecord(
        edge_index=edge_index,
        face_a=fa.face_id,
        face_b=fb.face_id,
        edge=edge,
        first=first,
        last=last,
        edge_tolerance=et,
        verify_tolerance=verify_tol,
        parameters=ts,
        xyz=xyz,
        uv_a=np.vstack(uva) if uva else np.zeros((0, 2)),
        uv_b=np.vstack(uvb) if uvb else np.zeros((0, 2)),
        max_surface_error_a=max_a,
        max_surface_error_b=max_b,
        max_cross_surface_error=max_cross,
        trim_ok=trim_ok,
        min_transversality=min_tr,
        max_transversality=max_tr,
        risk_flags=tuple(risk),
        repaired_same_parameter=repaired,
        exact_curve_on_surface_checked=True,
        exact_surface_error_a=exact_a,
        exact_surface_error_b=exact_b,
    )


def _shape_distance(a, b) -> Optional[float]:
    from OCP.BRepExtrema import BRepExtrema_DistShapeShape

    d = BRepExtrema_DistShapeShape(a, b)
    if not d.IsDone():
        d.Perform()
    if not d.IsDone():
        return None
    return float(d.Value())



def _run_section_engine(fa: FaceRecord, fb: FaceRecord, *,
                        approximation: bool,
                        fuzzy: float,
                        parallel: bool,
                        use_obb: bool):
    """Run one OCCT Section construction mode and return edges/vertices."""
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
    from OCP.BRep import BRep_Tool
    from OCP.TopAbs import TopAbs_EDGE, TopAbs_VERTEX
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS

    sec = BRepAlgoAPI_Section(fa.face, fb.face, False)
    sec.SetNonDestructive(True)
    sec.SetRunParallel(bool(parallel))
    sec.SetUseOBB(bool(use_obb))
    if fuzzy > 0.0:
        sec.SetFuzzyValue(float(fuzzy))
    sec.Approximation(bool(approximation))
    sec.ComputePCurveOn1(True)
    sec.ComputePCurveOn2(True)
    sec.Build()
    if not sec.IsDone():
        mode = "approx" if approximation else "nonapprox-shadow"
        raise IntersectionError(
            f"OCCT {mode} section failed for faces "
            f"{fa.face_id}/{fb.face_id}",
            kind="SectionEngineFailure")

    edges = []
    ex = TopExp_Explorer(sec.Shape(), TopAbs_EDGE)
    while ex.More():
        edges.append(TopoDS.Edge(ex.Current()))
        ex.Next()

    vertices = []
    ex = TopExp_Explorer(sec.Shape(), TopAbs_VERTEX)
    while ex.More():
        v = TopoDS.Vertex(ex.Current())
        vertices.append(_p3(BRep_Tool.Pnt_s(v)))
        ex.Next()
    return edges, vertices


def _edge_length(edge) -> float:
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps

    g = GProp_GProps()
    BRepGProp.LinearProperties_s(edge, g, False, False)
    return float(g.Mass())


def _edge_compound(records: list[SectionEdgeRecord]):
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Compound

    b = BRep_Builder()
    out = TopoDS_Compound()
    b.MakeCompound(out)
    for r in records:
        b.Add(out, r.edge)
    return out


def _point_shape_distance(point: np.ndarray, shape) -> float:
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
    from OCP.BRepExtrema import BRepExtrema_DistShapeShape
    from OCP.gp import gp_Pnt

    v = BRepBuilderAPI_MakeVertex(
        gp_Pnt(float(point[0]), float(point[1]), float(point[2]))).Vertex()
    d = BRepExtrema_DistShapeShape(v, shape)
    if not d.IsDone():
        d.Perform()
    if not d.IsDone():
        return float("inf")
    return float(d.Value())


def _subsample_xyz(records: list[SectionEdgeRecord],
                   max_per_edge: int = 33) -> list[np.ndarray]:
    pts = []
    for r in records:
        n = len(r.xyz)
        if not n:
            continue
        if n <= max_per_edge:
            idx = range(n)
        else:
            idx = np.unique(np.linspace(
                0, n - 1, max_per_edge, dtype=np.int64))
        pts.extend(r.xyz[int(i)] for i in idx)
    return pts


def _crosscheck_section_modes(primary: list[SectionEdgeRecord],
                              shadow: list[SectionEdgeRecord],
                              *,
                              base_tol: float
                              ) -> tuple[float, float]:
    """Require two OCCT section construction modes to describe one curve set.

    This is deliberately called a cross-check, not an exact Hausdorff proof.
    Both edge sets are independently verified against both trimmed input faces
    before this function is called. We then compare them bidirectionally in 3D
    and compare total curve length so a missing branch/segment cannot pass just
    because the surviving branch is close.
    """
    if not primary or not shadow:
        raise IntersectionError(
            "primary/shadow section edge set cardinality disagrees",
            kind="SectionConstructionDisagreement")

    pshape = _edge_compound(primary)
    sshape = _edge_compound(shadow)
    ppts = _subsample_xyz(primary)
    spts = _subsample_xyz(shadow)
    if not ppts or not spts:
        raise IntersectionError(
            "primary/shadow section has no comparable samples",
            kind="SectionConstructionDisagreement")

    d_ps = max((_point_shape_distance(p, sshape) for p in ppts),
               default=float("inf"))
    d_sp = max((_point_shape_distance(p, pshape) for p in spts),
               default=float("inf"))
    max_dist = max(d_ps, d_sp)

    tol = max(
        8.0 * float(base_tol),
        2.0 * max([r.verify_tolerance for r in primary + shadow],
                  default=float(base_tol)))
    if not np.isfinite(max_dist) or max_dist > tol:
        raise IntersectionError(
            f"approx/nonapprox section disagreement {max_dist:.6g} "
            f"exceeds cross-check tolerance {tol:.6g}",
            kind="SectionConstructionDisagreement")

    lp = sum(_edge_length(r.edge) for r in primary)
    ls = sum(_edge_length(r.edge) for r in shadow)
    rel = abs(lp - ls) / max(abs(lp), abs(ls), tol, 1e-300)
    len_tol = max(32.0 * tol, 2e-5 * max(lp, ls, 1.0))
    if abs(lp - ls) > len_tol:
        raise IntersectionError(
            f"approx/nonapprox section total length differs by "
            f"{abs(lp-ls):.6g} (rel {rel:.6g})",
            kind="SectionConstructionDisagreement")

    return float(max_dist), float(rel)



@dataclass
class CompletenessProbeReport:
    """Outcome of the raw-intersector completeness probe for one face pair.

    components holds one dict per material raw component with its
    per-interval matching numbers (see _match_raw_component); a component
    whose intervals cannot all be accounted for raises
    SectionCompletenessMismatch instead of appearing here.
    """
    raw_curve_count: int = 0
    trimmed_components: int = 0
    unmatched_components: int = 0
    components: list = field(default_factory=list)
    max_distance: float = 0.0


# ---------------------------------------------------------------------------
# G4 rework: completeness-probe scope. The probe runs for EVERY candidate
# face pair that produces section curves, not just freeform pairs. The
# earlier G4 classification (IntAna closed-form pairs skip the probe)
# rested on the assumption that OCCT's closed-form path cannot drop a
# branch its own lower-level intersector saw. That assumption is exactly
# what the probe is meant to check, so the classification was circular:
# it used a claim about OCCT's exactness to skip verifying OCCT. The probe
# now runs unconditionally; the per-interval matching rule (below) keeps
# analytic pairs cheap because their raw curves match at depth 0.
#
# The probe can still be disabled process-wide with the environment
# variable BREPKERNEL_COMPLETENESS_PROBE=0 (also accepts off/false/no);
# the default is on. This is the escape hatch used to measure the
# probe's cost (G4 item 2) and to recover if a future OCCT version ever
# makes the probe itself unreliable. Both section_face_pair and
# intersect_models honor it so direct and boolean_brep paths agree.
# ---------------------------------------------------------------------------

def _completeness_probe_enabled() -> bool:
    """Process-wide kill switch for the completeness probe (G4 item 2)."""
    import os
    return os.environ.get(
        "BREPKERNEL_COMPLETENESS_PROBE", "1").strip().lower() not in (
            "0", "off", "false", "no", "n")


# ---------------------------------------------------------------------------
# G3 rework: per-interval completeness matching.
#
# The probe checks each bounded raw IntTools_FaceFace curve component
# against the verified Section edges with an adaptive per-interval rule,
# replacing the old global coverage/approximation-gap rule. Rationale:
# a dropped MIDDLE of a branch is a completeness failure even when the
# kept ends match (the old rule's approximation_gap blanket-accepted
# any component with one near sample), while a raw tail that overshoots
# a trim boundary is trim noise, not a dropped branch (the old 95%
# coverage rule refused real cases like the F2 cone/sphere pair, whose
# raw curve overshoots the trim by 9.4e-5 past the Section edge end).
# ---------------------------------------------------------------------------

#: Initial uniform partition count for per-interval matching. With the
#: 3-point stencil (endpoints + midpoint, deduped) this gives 41 distinct
#: samples at depth 0.
_PROBE_INTERVALS0 = 20
#: Subdivision depth cap: a violation must persist down to ~1/640 of the
#: raw range before the probe refuses.
_PROBE_MAX_DEPTH = 5


def _face_wires_compound(fa_face, fb_face):
    """Compound of all trim wires of two faces (boundary proximity checks)."""
    from OCP.BRep import BRep_Builder
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopAbs import TopAbs_WIRE
    from OCP.TopoDS import TopoDS_Compound, TopoDS

    comp = TopoDS_Compound()
    BRep_Builder().MakeCompound(comp)
    for f in (fa_face, fb_face):
        ex = TopExp_Explorer(f, TopAbs_WIRE)
        while ex.More():
            BRep_Builder().Add(comp, TopoDS.Wire(ex.Current()))
            ex.Next()
    return comp


def _match_raw_component(ic_index, c3, c2a, c2b, fa, fb, *,
                         target, edge_ends, wires_thunk,
                         t0, t1, trim_tol, tight_tol, tol_i, scale):
    """Match one raw curve component against the verified edge set.

    Adaptive per-interval rule (G3 rework): the raw range [t0, t1] is
    split into _PROBE_INTERVALS0 uniform intervals; each interval is
    sampled at its endpoints and midpoint (in-trim classification at
    trim_tol). Every in-trim sample must lie within tol_i of the verified
    edge set, subject to the boundary-tail provision below. An interval
    with a violating sample subdivides (binary) until it matches or
    _PROBE_MAX_DEPTH is reached; every leaf interval must be accounted
    for (matched, boundary tail, or out of trim) or the component is
    refused.

    Boundary-tail provision: in an end interval (touching t0 or t1), a
    sample beyond tol_i is still accounted for when (a) its nearest point
    on the edge set is an edge END (the sample sits past the edge's end,
    not beside the edge), (b) the sample lies in the trim-boundary band
    (OUT at tight_tol but in/on at trim_tol), and (c) that edge end sits
    within trim_tol of the trim wires. This is the F2 situation: the raw
    IntTools curve overshoots the trim slightly past where the Section
    edge ends; the overshoot is trim-boundary noise, not a dropped
    branch. All three quantities are OCCT-reported or existing pipeline
    quantities (I3); no tolerance is widened.

    Returns (verdict, leaf_intervals, stats). verdict is "skipped"
    (fewer than 3 in-trim samples: point-contact evidence, same as the
    pre-G3 rule), "matched", or "refused".
    """
    from OCP.BRepClass import BRepClass_FaceClassifier
    from OCP.TopAbs import TopAbs_IN, TopAbs_ON
    from OCP.gp import gp_Pnt2d

    def classify(t, tol):
        ua = c2a.Value(float(t))
        ub = c2b.Value(float(t))
        ca = BRepClass_FaceClassifier(
            fa.face, gp_Pnt2d(float(ua.X()), float(ua.Y())), tol, True)
        cb = BRepClass_FaceClassifier(
            fb.face, gp_Pnt2d(float(ub.X()), float(ub.Y())), tol, True)
        return (ca.State() in (TopAbs_IN, TopAbs_ON)
                and cb.State() in (TopAbs_IN, TopAbs_ON))

    # 21-point materiality pre-scan, same as the pre-G3 rule: isolated
    # endpoint hits are point-contact evidence, not a trimmed component.
    # These are exactly the depth-0 interval endpoints, so the recursion
    # below is guaranteed to revisit every in-trim pre-scan point.
    n_probe = 21
    in_trim_pre = sum(
        1 for j in range(n_probe)
        if classify(t0 + (t1 - t0) * (j / (n_probe - 1.0)), trim_tol))
    stats = {"in_trim_samples": in_trim_pre, "nearest": 0.0,
             "max_sample_distance": 0.0}
    if in_trim_pre < 3:
        return "skipped", [], stats
    if target is None:
        # Material raw component but no verified edges at all: nothing
        # can be "the same branch".
        return "refused", [], stats

    endpoint_eps = 64.0 * 2.220446049250313e-16 * max(1.0, float(scale))
    point_cache: dict = {}
    trim_cache: dict = {}
    dist_cache: dict = {}
    band_cache: dict = {}
    wires_box: list = []

    def point_at(t):
        p = point_cache.get(t)
        if p is None:
            p = _p3(c3.Value(float(t)))
            point_cache[t] = p
        return p

    def in_trim(t):
        v = trim_cache.get(t)
        if v is None:
            v = classify(t, trim_tol)
            trim_cache[t] = v
        return v

    def nearest_info(t):
        """(distance to edge set, edge-end index or None).

        The end index is set when the nearest edge end is essentially as
        close as the nearest edge point, i.e. the sample sits past the
        edge's end rather than beside the edge.
        """
        v = dist_cache.get(t)
        if v is None:
            p = point_at(t)
            d = _point_shape_distance(p, target)
            best_k = None
            best_dv = float("inf")
            for k, (ep, _ei) in enumerate(edge_ends):
                dv = float(np.linalg.norm(p - ep))
                if dv < best_dv:
                    best_dv = dv
                    best_k = k
            v = (d, best_k if best_dv <= d + endpoint_eps else None)
            dist_cache[t] = v
        return v

    def in_boundary_band(t):
        v = band_cache.get(t)
        if v is None:
            # Samples reaching here are in-trim at trim_tol; the band is
            # OUT at the tight tolerance.
            v = not classify(t, tight_tol)
            band_cache[t] = v
        return v

    def edge_end_at_boundary(k):
        if not wires_box:
            wires_box.append(wires_thunk())
        return _point_shape_distance(edge_ends[k][0], wires_box[0]) <= trim_tol

    leaf_intervals: list = []

    def match_interval(a, b, depth, is_end_interval):
        stencil = [a, 0.5 * (a + b), b]
        seen = set()
        samples = []
        for t in stencil:
            if t in seen:
                continue
            seen.add(t)
            if in_trim(t):
                samples.append(t)
        rec = {"t_start": float(a), "t_end": float(b), "depth": depth,
               "n_in_trim": len(samples)}
        if not samples:
            rec["verdict"] = "out_of_trim"
            rec["max_distance"] = 0.0
            rec["boundary_tail_samples"] = 0
            leaf_intervals.append(rec)
            return True
        violating = []
        maxd = 0.0
        n_tail = 0
        for t in samples:
            d, end_k = nearest_info(t)
            if not np.isfinite(d):
                violating.append((t, d, "non-finite distance"))
                continue
            stats["max_sample_distance"] = max(
                stats["max_sample_distance"], d)
            maxd = max(maxd, d)
            if d <= tol_i:
                continue
            if (is_end_interval and end_k is not None
                    and in_boundary_band(t)
                    and edge_end_at_boundary(end_k)):
                n_tail += 1
                continue
            violating.append((t, d, "exceeds match tolerance"))
        # nearest is recomputed from the distance cache at the end.
        rec["max_distance"] = float(maxd)
        rec["boundary_tail_samples"] = n_tail
        if not violating:
            rec["verdict"] = ("matched_with_boundary_tails"
                              if n_tail else "matched")
            leaf_intervals.append(rec)
            return True
        if depth >= _PROBE_MAX_DEPTH:
            rec["verdict"] = "violated"
            rec["violations"] = [
                {"t": float(t), "distance": float(d), "reason": r}
                for (t, d, r) in violating]
            leaf_intervals.append(rec)
            return False
        mid = 0.5 * (a + b)
        ok_left = match_interval(a, mid, depth + 1,
                                 is_end_interval and a == t0)
        ok_right = match_interval(mid, b, depth + 1,
                                  is_end_interval and b == t1)
        return ok_left and ok_right

    h = (t1 - t0) / _PROBE_INTERVALS0
    ok = True
    for k in range(_PROBE_INTERVALS0):
        a = t0 + k * h
        b = t0 + (k + 1) * h if k + 1 < _PROBE_INTERVALS0 else t1
        if not match_interval(a, b, 0,
                              k == 0 or k == _PROBE_INTERVALS0 - 1):
            ok = False
    finite_ds = [v[0] for v in dist_cache.values()
                 if np.isfinite(v[0])]
    stats["nearest"] = float(min(finite_ds)) if finite_ds else 0.0
    stats["in_trim_samples"] = sum(1 for v in trim_cache.values() if v)
    return ("matched" if ok else "refused"), leaf_intervals, stats


def _raw_intersector_completeness_probe(
        fa: FaceRecord, fb: FaceRecord,
        verified: list[SectionEdgeRecord], *,
        base_tol: float, fuzzy: float, parallel: bool,
        max_section_tol: Optional[float] = None
        ) -> CompletenessProbeReport:
    """Probe Section post-processing against lower-level face/face curves.

    OCCT's BRepAlgoAPI_Section ultimately consumes IntTools_FaceFace
    curves. This probe reruns that lower-level intersector with the same
    approximation tolerance used by BOPAlgo_PaveFiller (1e-7) and requires
    every bounded raw curve component that actually lies in/on both
    trimmed faces to be represented by the final verified Section edges.

    Matching is per raw curve with the adaptive per-interval rule in
    _match_raw_component: the raw range is partitioned, every in-trim
    sample must lie within tol_i of the verified edge set (tol_i =
    max(16*base_tol, 4*max_verify_tol, 2*ic.Tolerance()), the pre-G3
    formula, unchanged), violating intervals subdivide to a fixed depth,
    and every leaf must be accounted for or the probe refuses with
    SectionCompletenessMismatch. End intervals additionally admit the
    boundary-tail provision (raw overshoot past a trim boundary where
    the Section edge ends). Per-interval numbers are recorded in each
    component record.

    A raw curve whose OCCT-reported tolerance exceeds the section
    tolerance ceiling refuses outright with RawIntersectionToleranceTooLoose
    (G3 item 1): the old code folded that tolerance into the match
    window, silently widening it.

    This is a post-processing completeness check, not a mathematical proof
    that IntTools_FaceFace itself discovered every true intersection
    branch.
    """
    from OCP.BRep import BRep_Tool
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    from OCP.IntTools import IntTools_FaceFace

    raw = IntTools_FaceFace()
    raw.SetParameters(True, True, True, 1e-7)
    if fuzzy > 0.0:
        raw.SetFuzzyValue(float(fuzzy))
    raw.Perform(fa.face, fb.face, bool(parallel))
    if not raw.IsDone():
        raise IntersectionError(
            f"lower-level face/face intersector failed for "
            f"{fa.face_id}/{fb.face_id}",
            kind="IntersectionCompletenessProbeFailed")
    # Mirror BOPAlgo_PaveFiller: it prepares the IntTools curves with
    # bSplitCurve=False before storing them in the face/face interference DS.
    raw.PrepareLines3D(False)

    lines = raw.Lines()
    raw_count = int(lines.Length())
    report = CompletenessProbeReport(raw_curve_count=raw_count)
    if raw_count == 0:
        return report

    target = _edge_compound(verified) if verified else None
    # Positional accuracy the Section verifier itself certifies for the
    # verified edges of this pair. A pipeline quantity, not a tuned number.
    max_verify_tol = max(
        (r.verify_tolerance for r in verified), default=float(base_tol))
    # Trim-classifier tolerance: unchanged pre-G3 floor.
    trim_tol = max(16.0 * float(base_tol), 4.0 * max_verify_tol)
    # Section tolerance ceiling, shared with _verify_section_edge (G3.1).
    tol_ceiling = _section_tolerance_ceiling(
        fa, fb, base_tol=base_tol, max_section_tol=max_section_tol)
    # Tight classifier tolerance for the trim-boundary band: the faces'
    # own OCCT tolerances (or base_tol), not the inflated trim_tol.
    tight_tol = max(float(base_tol),
                    float(BRep_Tool.Tolerance_s(fa.face)),
                    float(BRep_Tool.Tolerance_s(fb.face)))
    scale = max(
        float(np.linalg.norm(fa.bbox_hi - fa.bbox_lo)),
        float(np.linalg.norm(fb.bbox_hi - fb.bbox_lo)),
        1.0)

    # Verified edge end points (curve parameter ends: the ground truth
    # for "the sample sits past the edge's end").
    edge_ends: list = []
    for ei, rec in enumerate(verified):
        c = BRepAdaptor_Curve(rec.edge)
        edge_ends.append((_p3(c.Value(c.FirstParameter())), ei))
        edge_ends.append((_p3(c.Value(c.LastParameter())), ei))
    wires_box: list = []

    def wires_thunk():
        if not wires_box:
            wires_box.append(_face_wires_compound(fa.face, fb.face))
        return wires_box[0]

    trimmed_components = 0
    refused: list = []

    for i in range(1, raw_count + 1):
        ic = lines.Value(i)
        c3 = ic.Curve()
        if c3 is None:
            continue
        t0 = float(c3.FirstParameter())
        t1 = float(c3.LastParameter())
        if not (np.isfinite(t0) and np.isfinite(t1) and t1 > t0):
            # A face-trimmed FF component should normally be bounded. An
            # unbounded raw line cannot safely certify completeness.
            raise IntersectionError(
                f"raw intersection curve {i} is unbounded",
                kind="IntersectionCompletenessProbeFailed")

        c2first = ic.FirstCurve2d()
        c2second = ic.SecondCurve2d()
        if raw.Face1().IsSame(fa.face):
            c2a, c2b = c2first, c2second
        else:
            c2a, c2b = c2second, c2first
        if c2a is None or c2b is None:
            raise IntersectionError(
                f"raw intersection curve {i} lacks bilateral p-curves",
                kind="IntersectionCompletenessProbeFailed")

        # G3 item 1: refuse when OCCT's reported raw tolerance exceeds the
        # section tolerance ceiling instead of widening the window.
        # ic.TangentialTolerance() is recorded below but not folded in: it
        # governs OCCT's tangential-contact classification, not 3D
        # positional deviation.
        curve_tol = _raw_curve_tolerance(ic)
        if curve_tol > tol_ceiling:
            raise IntersectionError(
                f"raw intersection curve {i}: OCCT-reported tolerance "
                f"{curve_tol:.6g} exceeds the maximum acceptable section "
                f"tolerance {tol_ceiling:.6g}; refusing rather than "
                f"widening the completeness window",
                kind="RawIntersectionToleranceTooLoose")
        tol_i = max(16.0 * float(base_tol),
                    4.0 * max_verify_tol,
                    2.0 * curve_tol)

        verdict, leaf_intervals, stats = _match_raw_component(
            i, c3, c2a, c2b, fa, fb,
            target=target, edge_ends=edge_ends, wires_thunk=wires_thunk,
            t0=t0, t1=t1, trim_tol=trim_tol, tight_tol=tight_tol,
            tol_i=tol_i, scale=scale)
        if verdict == "skipped":
            continue
        trimmed_components += 1
        comp_record = {
            "component_index": i,
            "face_a": fa.face_id,
            "face_b": fb.face_id,
            "curve_tolerance": curve_tol,
            "tangential_tolerance": float(ic.TangentialTolerance()),
            "match_tolerance": float(tol_i),
            "tolerance_ceiling": float(tol_ceiling),
            "in_trim_samples": stats["in_trim_samples"],
            "nearest_distance": float(stats["nearest"]),
            "max_sample_distance": float(stats["max_sample_distance"]),
            "initial_intervals": _PROBE_INTERVALS0,
            "max_subdivision_depth": _PROBE_MAX_DEPTH,
            "leaf_interval_count": len(leaf_intervals),
            "max_depth_reached": max(
                (r["depth"] for r in leaf_intervals), default=0),
            "boundary_tail_intervals": sum(
                1 for r in leaf_intervals
                if r["verdict"] == "matched_with_boundary_tails"),
            "leaf_intervals": leaf_intervals,
            "verdict": verdict,
        }
        report.components.append(comp_record)
        report.max_distance = max(report.max_distance,
                                  float(stats["max_sample_distance"]))
        if verdict == "refused":
            bad = [r for r in leaf_intervals
                   if r["verdict"] == "violated"]
            if bad:
                worst = max(bad, key=lambda r: r["max_distance"])
                leaf_detail = (
                    f"{len(bad)} interval(s) violate the per-interval "
                    f"match at max depth {_PROBE_MAX_DEPTH}; worst leaf "
                    f"[{worst['t_start']:.6g}, {worst['t_end']:.6g}] "
                    f"max_distance={worst['max_distance']:.6g}; ")
            else:
                # No verified edges at all (target is None).
                leaf_detail = "no verified Section edges to match; "
            refused.append(
                f"component {i}: {leaf_detail}"
                f"(match_tolerance={tol_i:.6g}, "
                f"curve_tolerance={curve_tol:.6g}, "
                f"tolerance_ceiling={tol_ceiling:.6g}); "
                f"in_trim_samples={stats['in_trim_samples']}, "
                f"nearest={stats['nearest']:.6g}")

    report.trimmed_components = trimmed_components
    report.unmatched_components = len(refused)
    if refused:
        raise IntersectionError(
            f"lower-level intersector exposes {len(refused)} trimmed curve "
            f"component(s) not represented by verified Section edges "
            f"(raw={raw_count}, trimmed={trimmed_components}, "
            f"max_distance={report.max_distance:.6g}); "
            + "; ".join(refused),
            kind="SectionCompletenessMismatch")
    return report


def _edge_is_boundary_contact(edge, fa, fb, tol: float) -> bool:
    """True when the edge lies ON the trim boundary of both faces (G2.5).

    Every (u, v) sample of the edge must classify TopAbs_ON against both
    face trims. A boundary contact creates no split; it is resolved (or
    left unresolved) by _resolve_contacts in split.py.
    """
    from OCP.BRepClass import BRepClass_FaceClassifier
    from OCP.gp import gp_Pnt2d
    from OCP.TopAbs import TopAbs_ON
    from OCP.TopoDS import TopoDS
    try:
        uvs_a = np.asarray(edge.uv_a, dtype=float).reshape(-1, 2)
        uvs_b = np.asarray(edge.uv_b, dtype=float).reshape(-1, 2)
    except Exception:
        return False
    if len(uvs_a) == 0 or len(uvs_b) != len(uvs_a):
        return False
    try:
        ffa = TopoDS.Face(fa.face)
        ffb = TopoDS.Face(fb.face)
    except Exception:
        return False
    for (ua, va), (ub, vb) in zip(uvs_a, uvs_b):
        try:
            ca = BRepClass_FaceClassifier(
                ffa, gp_Pnt2d(float(ua), float(va)), float(tol))
            cb = BRepClass_FaceClassifier(
                ffb, gp_Pnt2d(float(ub), float(vb)), float(tol))
        except Exception:
            return False
        if ca.State() != TopAbs_ON or cb.State() != TopAbs_ON:
            return False
    return True


def _mark_boundary_contacts(result, fa, fb, tol: float) -> None:
    """Flag boundary-contact section edges on a pair result (G2.5)."""
    if result.status not in ("curve", "curve_near_tangent"):
        return
    for e in result.edges:
        edge_tol = max(float(tol), float(e.verify_tolerance),
                       float(e.edge_tolerance))
        e.is_boundary_contact = _edge_is_boundary_contact(e, fa, fb, edge_tol)
    if result.edges and all(e.is_boundary_contact for e in result.edges):
        result.status = "boundary_contact"
        result.notes.append(
            "all section edges lie on both trim boundaries; "
            "no split created")


def section_face_pair(fa: FaceRecord, fb: FaceRecord, *,
                      base_tol: float = 1e-7,
                      chord_tol: Optional[float] = None,
                      contact_tol: Optional[float] = None,
                      fuzzy: float = 0.0,
                      parallel: bool = True,
                      use_obb: bool = True,
                      tangent_sin_tol: float = 1e-4,
                      max_section_tol: Optional[float] = None,
                      crosscheck_nonapprox: bool = False,
                      completeness_probe: bool = True
                      ) -> FaceIntersectionResult:
    """Intersect one pair of *trimmed* faces and verify all section curves."""
    from OCP.BRep import BRep_Tool

    if not base_tol > 0:
        raise ValueError("base_tol must be positive")
    if fuzzy < 0:
        raise ValueError("fuzzy must be >= 0")
    # G6: the near-contact band is tolerance-aware in both directions.
    # A pair of faces whose OCCT tolerances sum to 2t is ambiguous at
    # distances up to about 2t, so the pipeline's explicit contact_tol
    # must never suppress that band: take the max. This rule was already
    # the default when contact_tol was left unset; the change only stops
    # an explicit small contact_tol from hiding a high-tolerance face.
    # Widening the band can only add typed refusals, never acceptances.
    contact_tol = max(
        0.0 if contact_tol is None else float(contact_tol),
        4.0 * base_tol,
        2.0 * float(BRep_Tool.Tolerance_s(fa.face)),
        2.0 * float(BRep_Tool.Tolerance_s(fb.face)))

    edges, vertices = _run_section_engine(
        fa, fb, approximation=True, fuzzy=float(fuzzy),
        parallel=bool(parallel), use_obb=bool(use_obb))

    verified = [
        _verify_section_edge(
            edge, fa, fb, i, base_tol=base_tol,
            chord_tol=chord_tol, tangent_sin_tol=tangent_sin_tol,
            max_section_tol=max_section_tol)
        for i, edge in enumerate(edges)
    ]

    if verified and crosscheck_nonapprox:
        shadow_edges, _ = _run_section_engine(
            fa, fb, approximation=False, fuzzy=float(fuzzy),
            parallel=bool(parallel), use_obb=bool(use_obb))
        shadow_verified = [
            _verify_section_edge(
                edge, fa, fb, i, base_tol=base_tol,
                chord_tol=chord_tol, tangent_sin_tol=tangent_sin_tol,
                max_section_tol=max_section_tol)
            for i, edge in enumerate(shadow_edges)
        ]
        max_shadow_distance, shadow_length_rel = _crosscheck_section_modes(
            verified, shadow_verified, base_tol=base_tol)
        for e in verified:
            e.shadow_crosschecked = True
            e.shadow_max_distance = max_shadow_distance
            e.shadow_length_rel_error = shadow_length_rel

    raw_count = raw_trimmed = raw_unmatched = 0
    completeness_distance = 0.0
    completeness_components: list = []
    # G4 rework: the probe runs for every candidate face pair that
    # produces section curves (no closed-form exemption); it can be
    # disabled process-wide with BREPKERNEL_COMPLETENESS_PROBE=0.
    if completeness_probe and _completeness_probe_enabled():
        probe = _raw_intersector_completeness_probe(
            fa, fb, verified, base_tol=base_tol,
            fuzzy=float(fuzzy), parallel=bool(parallel),
            max_section_tol=max_section_tol)
        raw_count = probe.raw_curve_count
        raw_trimmed = probe.trimmed_components
        raw_unmatched = probe.unmatched_components
        completeness_distance = probe.max_distance
        completeness_components = probe.components

    # Vertices include edge endpoints. They are useful diagnostic data, but
    # only a no-edge vertex set is a pure point-contact result.
    if verified:
        status = ("curve_near_tangent"
                  if any("near_tangent" in e.risk_flags for e in verified)
                  else "curve")
        return FaceIntersectionResult(
            fa.face_id, fb.face_id, status, verified, vertices,
            min_distance=0.0,
            raw_curve_count=raw_count,
            raw_trimmed_components=raw_trimmed,
            raw_unmatched_components=raw_unmatched,
            completeness_max_distance=completeness_distance,
            completeness_components=completeness_components)

    if vertices:
        return FaceIntersectionResult(
            fa.face_id, fb.face_id, "point_contact", [], vertices,
            min_distance=0.0,
            raw_curve_count=raw_count,
            raw_trimmed_components=raw_trimmed,
            raw_unmatched_components=raw_unmatched,
            completeness_max_distance=completeness_distance,
            completeness_components=completeness_components,
            notes=["zero-dimensional contact: topology-sensitive"])

    d = _shape_distance(fa.face, fb.face)
    if d is None:
        return FaceIntersectionResult(
            fa.face_id, fb.face_id, "distance_unknown", [],
            min_distance=None, section_done=True,
            notes=["section returned no curve and exact distance failed"])
    if d <= float(contact_tol):
        return FaceIntersectionResult(
            fa.face_id, fb.face_id, "ambiguous_contact", [],
            min_distance=d, section_done=True,
            notes=[f"no section curve but face distance {d:.6g} "
                   f"<= contact tolerance {contact_tol:.6g}"])
    return FaceIntersectionResult(
        fa.face_id, fb.face_id, "disjoint", [],
        min_distance=d)


def intersect_models(a: BRepModel, b: BRepModel, *,
                     broadphase_pad: float = 0.0,
                     base_tol: float = 1e-7,
                     chord_tol: Optional[float] = None,
                     contact_tol: Optional[float] = None,
                     fuzzy: float = 0.0,
                     parallel: bool = True,
                     use_obb: bool = True,
                     tangent_sin_tol: float = 1e-4,
                     max_section_tol: Optional[float] = None,
                     crosscheck_nonapprox: bool = False,
                     completeness_probe: bool = True
                     ) -> ModelIntersectionResult:
    """Run verified section work only for conservative candidate face pairs."""
    candidates = candidate_face_pairs(a, b, pad=float(broadphase_pad))
    results: list[FaceIntersectionResult] = []
    section_calls = 0
    verified_edges = 0
    point_contacts = 0
    ambiguous = 0
    shadow_calls = 0
    shadow_edges = 0
    max_shadow_distance = 0.0
    completeness_probes = 0
    raw_curve_count = 0
    raw_trimmed_components = 0
    raw_unmatched_components = 0
    completeness_max_distance = 0.0
    completeness_components: list = []

    by_a = {f.face_id: f for f in a.faces}
    by_b = {f.face_id: f for f in b.faces}
    # G2.3: coincidence stage runs between candidate selection and the
    # section. Coincident pairs skip the section (it would return
    # degenerate geometry); undecidable pairs refuse loudly.
    # G2-planar: three-case rule; the undecidable record carries the
    # measured deviation and every tolerance consulted.
    from .coincidence import (classify_support_pair, CoincidentPairRecord,
                              CoincidenceUndecidable, coincidence_band)
    coincident_pairs: list = []
    band = coincidence_band(base_tol)
    for pair in candidates:
        fa = by_a[pair["face_a"]]
        fb = by_b[pair["face_b"]]
        coinc = classify_support_pair(fa, fb, base_tol,
                                      contact_tol=contact_tol)
        if isinstance(coinc, CoincidenceUndecidable):
            if coinc.reason in ("near_coincident_planar",
                                "recognized_beyond_tolerance"):
                kind = "NearCoincidentFaces"
            else:
                kind = "CoincidenceUndecidable"
            tol_txt = "; ".join(
                f"{t['name']}={t['value']:.3g}" for t in coinc.tolerances)
            exc = IntersectionError(
                f"face pair A{fa.face_id}/B{fb.face_id}: support "
                f"coincidence undecidable ({coinc.reason}); measured "
                f"deviation {coinc.deviation:.3g} ({coinc.deviation_kind}); "
                f"tolerances consulted: {tol_txt}; refusing rather than "
                f"guessing",
                kind=kind)
            exc.evidence = {
                "reason": coinc.reason,
                "deviation": coinc.deviation,
                "deviation_kind": coinc.deviation_kind,
                "tolerances": coinc.tolerances,
            }
            raise exc
        if isinstance(coinc, tuple) and coinc[0] == "coincident":
            rec = CoincidentPairRecord(
                face_a=fa.face_id, face_b=fb.face_id, sense=coinc[1],
                certification=coinc[2] if len(coinc) > 2 else "exact")
            rec.tolerances.append({
                "name": "coincidence_band", "value": band,
                "role": "coincident support gate (4x base_tol)"})
            rec.tolerances.append({
                "name": "base_tol", "value": float(base_tol),
                "role": "input base tolerance"})
            rec.tolerances.append({
                "name": "certification", "value": rec.certification,
                "role": "how coincidence was established: "
                        "exact or tolerance_certified"})
            coincident_pairs.append(rec)
            results.append(FaceIntersectionResult(
                fa.face_id, fb.face_id, "coincident", [],
                min_distance=0.0, section_done=False,
                notes=[f"coincident supports, sense={coinc[1]}; "
                       f"handled by overlap split"]))
            continue
        # Distinct supports: run the section as before.
        section_calls += 1
        r = section_face_pair(
            fa, fb, base_tol=base_tol, chord_tol=chord_tol,
            contact_tol=contact_tol, fuzzy=fuzzy, parallel=parallel,
            use_obb=use_obb, tangent_sin_tol=tangent_sin_tol,
            max_section_tol=max_section_tol,
            crosscheck_nonapprox=bool(crosscheck_nonapprox),
            completeness_probe=bool(completeness_probe))
        # G2.5: flag boundary-contact edges (ON both trim boundaries);
        # they create no split.
        _mark_boundary_contacts(r, fa, fb, base_tol)
        results.append(r)
        verified_edges += len(r.edges)
        if any(e.shadow_crosschecked for e in r.edges):
            shadow_calls += 1
            shadow_edges += len(r.edges)
            max_shadow_distance = max(
                max_shadow_distance,
                max((float(e.shadow_max_distance or 0.0)
                     for e in r.edges), default=0.0))
        if completeness_probe and _completeness_probe_enabled():
            completeness_probes += 1
            raw_curve_count += int(r.raw_curve_count)
            raw_trimmed_components += int(r.raw_trimmed_components)
            raw_unmatched_components += int(r.raw_unmatched_components)
            completeness_max_distance = max(
                completeness_max_distance,
                float(r.completeness_max_distance))
            completeness_components.extend(r.completeness_components)
        point_contacts += int(r.status == "point_contact")
        ambiguous += int(r.status in ("ambiguous_contact",
                                      "distance_unknown",
                                      "curve_near_tangent"))

    return ModelIntersectionResult(
        pairs=results,
        candidate_pairs=len(candidates),
        section_calls=section_calls,
        verified_edges=verified_edges,
        point_contacts=point_contacts,
        ambiguous_contacts=ambiguous,
        skipped_by_broadphase=max(0, len(a.faces) * len(b.faces)
                                  - len(candidates)),
        coincident_pairs=coincident_pairs,
        boundary_edges=[
            e for r in results for e in r.edges
            if getattr(e, "is_boundary_contact", False)],
        shadow_section_calls=shadow_calls,
        shadow_verified_edges=shadow_edges,
        max_shadow_distance=max_shadow_distance,
        completeness_probes=completeness_probes,
        raw_curve_count=raw_curve_count,
        raw_trimmed_components=raw_trimmed_components,
        raw_unmatched_components=raw_unmatched_components,
        completeness_max_distance=completeness_max_distance,
        completeness_components=completeness_components,
    )
