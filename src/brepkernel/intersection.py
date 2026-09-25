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
    shadow_section_calls: int = 0
    shadow_verified_edges: int = 0
    max_shadow_distance: float = 0.0
    completeness_probes: int = 0
    raw_curve_count: int = 0
    raw_trimmed_components: int = 0
    raw_unmatched_components: int = 0
    completeness_max_distance: float = 0.0

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
    scale = max(
        float(np.linalg.norm(fa.bbox_hi - fa.bbox_lo)),
        float(np.linalg.norm(fb.bbox_hi - fb.bbox_lo)),
        1.0)
    tol_limit = (max(128.0 * float(base_tol), 1e-8 * scale)
                 if max_section_tol is None
                 else float(max_section_tol))
    if not tol_limit > 0.0:
        raise ValueError("max_section_tol must be positive")
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



def _raw_intersector_completeness_probe(
        fa: FaceRecord, fb: FaceRecord,
        verified: list[SectionEdgeRecord], *,
        base_tol: float, fuzzy: float, parallel: bool
        ) -> tuple[int, int, int, float]:
    """Probe Section post-processing against lower-level face/face curves.

    OCCT's BRepAlgoAPI_Section ultimately consumes IntTools_FaceFace curves.
    This probe reruns that lower-level intersector with the same approximation
    tolerance used by BOPAlgo_PaveFiller (1e-7), samples every bounded raw
    curve component that actually lies in/on both trimmed faces, and requires
    those components to be represented by the final verified Section edges.

    This is a post-processing completeness check, not a mathematical proof
    that IntTools_FaceFace itself discovered every true intersection branch.
    """
    from OCP.BRepClass import BRepClass_FaceClassifier
    from OCP.IntTools import IntTools_FaceFace
    from OCP.TopAbs import TopAbs_IN, TopAbs_ON
    from OCP.gp import gp_Pnt2d

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
    if raw_count == 0:
        return 0, 0, 0, 0.0

    target = _edge_compound(verified) if verified else None
    tol = max(
        16.0 * float(base_tol),
        4.0 * max((r.verify_tolerance for r in verified),
                  default=float(base_tol)))
    trimmed_components = 0
    unmatched = 0
    max_dist = 0.0

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

        in_points = []
        # Enough global coverage to recognize separate trimmed components
        # without turning this diagnostic into another high-density sampler.
        for j in range(17):
            t = t0 + (t1 - t0) * (j / 16.0)
            ua = c2a.Value(float(t))
            ub = c2b.Value(float(t))
            ca = BRepClass_FaceClassifier(
                fa.face, gp_Pnt2d(float(ua.X()), float(ua.Y())),
                tol, True)
            cb = BRepClass_FaceClassifier(
                fb.face, gp_Pnt2d(float(ub.X()), float(ub.Y())),
                tol, True)
            if (ca.State() in (TopAbs_IN, TopAbs_ON)
                    and cb.State() in (TopAbs_IN, TopAbs_ON)):
                in_points.append(_p3(c3.Value(float(t))))

        # Require several points before treating a raw line as a material
        # trimmed component; isolated endpoint hits are point-contact evidence.
        if len(in_points) < 3:
            continue
        trimmed_components += 1

        if target is None:
            unmatched += 1
            continue

        d = max((_point_shape_distance(p, target) for p in in_points),
                default=float("inf"))
        max_dist = max(max_dist, float(d))
        if not np.isfinite(d) or d > tol:
            unmatched += 1

    if unmatched:
        raise IntersectionError(
            f"lower-level intersector exposes {unmatched} trimmed curve "
            f"component(s) not represented by verified Section edges "
            f"(raw={raw_count}, trimmed={trimmed_components}, "
            f"max_distance={max_dist:.6g}, tol={tol:.6g})",
            kind="SectionCompletenessMismatch")
    return raw_count, trimmed_components, unmatched, float(max_dist)


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
    if contact_tol is None:
        contact_tol = max(4.0 * base_tol,
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
    # Restrict the extra intersector run to pairs involving freeform faces.
    # Analytic-only pairs already have exact closed-form/OCCT paths and do not
    # justify doubling FF work in the Tier B/C freeform route.
    if completeness_probe and (fa.freeform is not None or fb.freeform is not None):
        (raw_count, raw_trimmed, raw_unmatched,
         completeness_distance) = _raw_intersector_completeness_probe(
            fa, fb, verified, base_tol=base_tol,
            fuzzy=float(fuzzy), parallel=bool(parallel))

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
            completeness_max_distance=completeness_distance)

    if vertices:
        return FaceIntersectionResult(
            fa.face_id, fb.face_id, "point_contact", [], vertices,
            min_distance=0.0,
            raw_curve_count=raw_count,
            raw_trimmed_components=raw_trimmed,
            raw_unmatched_components=raw_unmatched,
            completeness_max_distance=completeness_distance,
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

    by_a = {f.face_id: f for f in a.faces}
    by_b = {f.face_id: f for f in b.faces}
    for pair in candidates:
        fa = by_a[pair["face_a"]]
        fb = by_b[pair["face_b"]]
        section_calls += 1
        r = section_face_pair(
            fa, fb, base_tol=base_tol, chord_tol=chord_tol,
            contact_tol=contact_tol, fuzzy=fuzzy, parallel=parallel,
            use_obb=use_obb, tangent_sin_tol=tangent_sin_tol,
            max_section_tol=max_section_tol,
            crosscheck_nonapprox=bool(crosscheck_nonapprox),
            completeness_probe=bool(completeness_probe))
        results.append(r)
        verified_edges += len(r.edges)
        if any(e.shadow_crosschecked for e in r.edges):
            shadow_calls += 1
            shadow_edges += len(r.edges)
            max_shadow_distance = max(
                max_shadow_distance,
                max((float(e.shadow_max_distance or 0.0)
                     for e in r.edges), default=0.0))
        if (fa.freeform is not None or fb.freeform is not None):
            completeness_probes += int(bool(completeness_probe))
            raw_curve_count += int(r.raw_curve_count)
            raw_trimmed_components += int(r.raw_trimmed_components)
            raw_unmatched_components += int(r.raw_unmatched_components)
            completeness_max_distance = max(
                completeness_max_distance,
                float(r.completeness_max_distance))
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
        shadow_section_calls=shadow_calls,
        shadow_verified_edges=shadow_edges,
        max_shadow_distance=max_shadow_distance,
        completeness_probes=completeness_probes,
        raw_curve_count=raw_curve_count,
        raw_trimmed_components=raw_trimmed_components,
        raw_unmatched_components=raw_unmatched_components,
        completeness_max_distance=completeness_max_distance,
    )
