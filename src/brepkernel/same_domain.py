"""Strict same-domain equivalence checks for Tier B/C fast paths.

OCCT exposes BOPTools_AlgoTools.AreFacesSameDomain(), but its implementation
uses an interior point from one face and checks that point on the other face.
That is useful evidence, not by itself a proof that two arbitrary trimmed
faces are identical.

This module deliberately layers several independent checks and prefers false
negatives over false positives.  It is used only to skip a Boolean when two
closed B-reps are independently constructed representations of the same
material boundary.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .freeform import FreeformError
from .step_ingest import BRepModel


class SameDomainError(FreeformError):
    pass


@dataclass
class FaceSameDomainEvidence:
    face_a: int
    face_b: int
    bbox_error: float
    area_a: float
    area_b: float
    area_rel_error: float
    perimeter_a: float
    perimeter_b: float
    perimeter_rel_error: float
    edge_count: int
    wire_count: int
    reverse_normals: Optional[bool]


@dataclass
class SameDomainResult:
    equivalent: bool
    reason: str
    matches: list[FaceSameDomainEvidence] = field(default_factory=list)
    candidate_counts: list[int] = field(default_factory=list)
    signed_volume_a: Optional[float] = None
    signed_volume_b: Optional[float] = None
    bbox_error: Optional[float] = None


def _bbox(shape) -> tuple[np.ndarray, np.ndarray]:
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib

    b = Bnd_Box()
    BRepBndLib.AddOptimal_s(shape, b, False, False)
    if b.IsVoid():
        z = np.zeros(3, dtype=np.float64)
        return z, z
    p0, p1 = b.CornerMin(), b.CornerMax()
    return (
        np.array([p0.X(), p0.Y(), p0.Z()], dtype=np.float64),
        np.array([p1.X(), p1.Y(), p1.Z()], dtype=np.float64),
    )


def _shape_scale(shape) -> float:
    lo, hi = _bbox(shape)
    return max(float(np.linalg.norm(hi - lo)),
               float(np.max(np.abs(lo))),
               float(np.max(np.abs(hi))), 1.0)


def _face_area(face) -> float:
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps

    g = GProp_GProps()
    err = BRepGProp.SurfaceProperties_s(face, g, 1e-10, False)
    if float(err) < 0.0:
        raise SameDomainError("adaptive face-area integration failed",
                              kind="SameDomainAreaFailed")
    return float(g.Mass())


def _perimeter(face) -> float:
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps

    g = GProp_GProps()
    BRepGProp.LinearProperties_s(face, g, False, False)
    return float(g.Mass())


def _signed_volume(shape) -> float:
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps

    g = GProp_GProps()
    err = BRepGProp.VolumePropertiesGK_s(
        shape, g, 1e-10, True, True, False, False, False)
    if float(err) < 0.0:
        raise SameDomainError("adaptive volume integration failed",
                              kind="SameDomainVolumeFailed")
    return float(g.Mass())


def _topology_signature(face) -> tuple[int, int]:
    from OCP.TopAbs import TopAbs_EDGE, TopAbs_WIRE
    from OCP.TopExp import TopExp_Explorer

    ne = 0
    ex = TopExp_Explorer(face, TopAbs_EDGE)
    while ex.More():
        ne += 1
        ex.Next()
    nw = 0
    ex = TopExp_Explorer(face, TopAbs_WIRE)
    while ex.More():
        nw += 1
        ex.Next()
    return ne, nw


def _rel_error(a: float, b: float, floor: float = 1e-300) -> float:
    return abs(a - b) / max(abs(a), abs(b), floor)


def _face_candidate(fa, fb, ida: int, idb: int, context, *,
                    bbox_tol: float, area_rel_tol: float,
                    perimeter_rel_tol: float, fuzz: float
                    ) -> Optional[FaceSameDomainEvidence]:
    from OCP.BOPTools import BOPTools_AlgoTools

    loa, hia = _bbox(fa)
    lob, hib = _bbox(fb)
    bbox_error = max(
        float(np.max(np.abs(loa - lob))),
        float(np.max(np.abs(hia - hib))))
    if bbox_error > bbox_tol:
        return None

    siga = _topology_signature(fa)
    sigb = _topology_signature(fb)
    if siga != sigb:
        return None

    aa, ab = _face_area(fa), _face_area(fb)
    are = _rel_error(aa, ab)
    if are > area_rel_tol:
        return None

    pa, pb = _perimeter(fa), _perimeter(fb)
    pre = _rel_error(pa, pb)
    if pre > perimeter_rel_tol:
        return None

    # OCCT's predicate is intentionally checked in both directions.  A one-way
    # interior-point test can accept a strict subset of another trimmed face.
    if not BOPTools_AlgoTools.AreFacesSameDomain_s(
            fa, fb, context, float(fuzz)):
        return None
    if not BOPTools_AlgoTools.AreFacesSameDomain_s(
            fb, fa, context, float(fuzz)):
        return None

    # Keep local normal-sense information for diagnostics, but do not use it
    # as a cross-model equivalence gate. Independently constructed coincident
    # faces can use different support-surface parameter frames even when their
    # enclosing closed solids represent the same material. Global material
    # orientation is checked from the signed closed-solid volume instead.
    try:
        # OCP 8 currently requires the nominally optional error argument.
        reverse = bool(BOPTools_AlgoTools.IsSplitToReverse_s(
            fa, fb, context, 0))
    except (TypeError, RuntimeError):
        # Diagnostic only. Closed-solid signed-volume orientation below is the
        # actual material-side invariant.
        reverse = None

    return FaceSameDomainEvidence(
        face_a=ida,
        face_b=idb,
        bbox_error=bbox_error,
        area_a=aa,
        area_b=ab,
        area_rel_error=are,
        perimeter_a=pa,
        perimeter_b=pb,
        perimeter_rel_error=pre,
        edge_count=siga[0],
        wire_count=siga[1],
        reverse_normals=reverse,
    )


def _perfect_matching(candidates: list[list[FaceSameDomainEvidence]],
                      n_b: int
                      ) -> Optional[list[FaceSameDomainEvidence]]:
    """Deterministic bipartite perfect matching, smallest candidate sets first."""
    order = sorted(range(len(candidates)),
                   key=lambda i: (len(candidates[i]), i))
    match_b: dict[int, FaceSameDomainEvidence] = {}

    def visit(ia: int, seen: set[int]) -> bool:
        rows = sorted(candidates[ia],
                      key=lambda e: (e.bbox_error,
                                     e.area_rel_error,
                                     e.perimeter_rel_error,
                                     e.face_b))
        for ev in rows:
            jb = ev.face_b
            if jb in seen:
                continue
            seen.add(jb)
            old = match_b.get(jb)
            if old is None:
                match_b[jb] = ev
                return True
            if visit(old.face_a, seen):
                match_b[jb] = ev
                return True
        return False

    for ia in order:
        if not candidates[ia] or not visit(ia, set()):
            return None
    if len(match_b) != n_b:
        return None

    by_a = {ev.face_a: ev for ev in match_b.values()}
    if len(by_a) != len(candidates):
        return None
    return [by_a[i] for i in range(len(candidates))]


def same_domain_models(a: BRepModel, b: BRepModel, *,
                       base_tol: float = 1e-7,
                       fuzz: Optional[float] = None,
                       area_rel_tol: float = 5e-8,
                       perimeter_rel_tol: float = 5e-8,
                       volume_rel_tol: float = 5e-8
                       ) -> SameDomainResult:
    """Strictly recognize independently built but equivalent closed B-reps.

    The recognizer intentionally handles only matching decompositions:
    same counts of solids/shells/faces and a one-to-one same-domain face map.
    Different valid decompositions of the same material are left to the normal
    Boolean/refusal path rather than guessed equivalent.
    """
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.IntTools import IntTools_Context

    if not base_tol > 0:
        raise ValueError("base_tol must be positive")
    if fuzz is None:
        fuzz = float(base_tol)
    if fuzz < 0:
        raise ValueError("fuzz must be >= 0")

    if not a.solids or not b.solids:
        return SameDomainResult(False, "closed solids required")
    if (len(a.solids) != len(b.solids)
            or len(a.shells) != len(b.shells)
            or len(a.faces) != len(b.faces)):
        return SameDomainResult(False, "topology counts differ")
    if not BRepCheck_Analyzer(a.shape, True).IsValid():
        return SameDomainResult(False, "operand A is B-rep invalid")
    if not BRepCheck_Analyzer(b.shape, True).IsValid():
        return SameDomainResult(False, "operand B is B-rep invalid")

    scale = max(_shape_scale(a.shape), _shape_scale(b.shape), 1.0)
    bbox_tol = max(8.0 * float(base_tol), 2e-10 * scale)

    loa, hia = _bbox(a.shape)
    lob, hib = _bbox(b.shape)
    model_bbox_error = max(
        float(np.max(np.abs(loa - lob))),
        float(np.max(np.abs(hia - hib))))
    if model_bbox_error > bbox_tol:
        return SameDomainResult(
            False, "model bounding boxes differ",
            bbox_error=model_bbox_error)

    try:
        va, vb = _signed_volume(a.shape), _signed_volume(b.shape)
    except SameDomainError as exc:
        return SameDomainResult(False, str(exc),
                                bbox_error=model_bbox_error)

    if _rel_error(abs(va), abs(vb)) > volume_rel_tol:
        return SameDomainResult(
            False, "model volumes differ",
            signed_volume_a=va, signed_volume_b=vb,
            bbox_error=model_bbox_error)
    # A negative/positive sign disagreement means the same support boundary is
    # being interpreted with opposite global orientation.
    if (va < 0.0) != (vb < 0.0):
        return SameDomainResult(
            False, "global material orientation differs",
            signed_volume_a=va, signed_volume_b=vb,
            bbox_error=model_bbox_error)

    context = IntTools_Context()
    candidates: list[list[FaceSameDomainEvidence]] = []
    try:
        for fa in a.faces:
            row = []
            for fb in b.faces:
                ev = _face_candidate(
                    fa.face, fb.face, fa.face_id, fb.face_id, context,
                    bbox_tol=bbox_tol,
                    area_rel_tol=float(area_rel_tol),
                    perimeter_rel_tol=float(perimeter_rel_tol),
                    fuzz=float(fuzz))
                if ev is not None:
                    row.append(ev)
            candidates.append(row)
    except Exception as exc:
        # A same-domain fast path is optional. Binding/projector failure is a
        # false negative, not permission to guess.
        return SameDomainResult(
            False, f"same-domain predicate failed: {type(exc).__name__}: {exc}",
            signed_volume_a=va, signed_volume_b=vb,
            bbox_error=model_bbox_error)

    matches = _perfect_matching(candidates, len(b.faces))
    if matches is None:
        return SameDomainResult(
            False, "no one-to-one bidirectional same-domain face matching",
            candidate_counts=[len(row) for row in candidates],
            signed_volume_a=va, signed_volume_b=vb,
            bbox_error=model_bbox_error)

    return SameDomainResult(
        True, "strict same-domain boundary match",
        matches=matches,
        candidate_counts=[len(row) for row in candidates],
        signed_volume_a=va,
        signed_volume_b=vb,
        bbox_error=model_bbox_error)


def same_domain_shapes(shape_a, shape_b, **kwargs) -> SameDomainResult:
    """Convenience wrapper for raw OCCT shapes."""
    from .step_ingest import index_shape
    return same_domain_models(index_shape(shape_a), index_shape(shape_b),
                              **kwargs)
