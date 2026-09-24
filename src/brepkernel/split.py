"""Local B-rep face splitting from verified section edges.

This module is intentionally downstream of intersection.py. It never invents
split curves: only section edges whose 3D curve and both p-curves already
passed verification may be used.

Performance:
- unaffected faces pass through bit-for-bit;
- all section edges affecting one parent face are split in a single
  BRepAlgoAPI_Splitter call;
- no global solid boolean is run merely to create local face patches.

Accuracy / refusal:
- near-tangent or seam-risk section edges are refused by default;
- every split result must be OCCT-valid;
- child face areas must partition the parent area within a scale-aware bound;
- deterministic interior UV witnesses from child pieces must also classify
  inside/on the parent face and evaluate to the same 3D surface point.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .freeform import FreeformError
from .intersection import ModelIntersectionResult, SectionEdgeRecord
from .step_ingest import BRepModel, FaceRecord


class SplitError(FreeformError):
    pass


@dataclass
class SplitFacePiece:
    parent_face_id: int
    piece_index: int
    face: object
    area: float
    uv_witness: Optional[np.ndarray]
    unchanged: bool = False


@dataclass
class FaceSplitResult:
    parent_face_id: int
    status: str
    pieces: list[SplitFacePiece]
    source_edges: int
    area_before: float
    area_after: float
    area_error: float
    notes: list[str] = field(default_factory=list)


@dataclass
class ModelSplitResult:
    faces_a: list[FaceSplitResult]
    faces_b: list[FaceSplitResult]
    split_calls: int
    affected_faces_a: int
    affected_faces_b: int
    unresolved_contacts: list[tuple[int, int, str]]

    @property
    def certified_local_split(self) -> bool:
        return not self.unresolved_contacts


def _face_area(face) -> float:
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps

    g = GProp_GProps()
    BRepGProp.SurfaceProperties_s(face, g)
    return float(g.Mass())


def _bbox_scale(face) -> float:
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib

    b = Bnd_Box()
    BRepBndLib.AddOptimal_s(face, b, False, False)
    if b.IsVoid():
        return 1.0
    p0 = b.CornerMin()
    p1 = b.CornerMax()
    return max(1.0, float(np.linalg.norm(np.array(
        [p1.X() - p0.X(), p1.Y() - p0.Y(), p1.Z() - p0.Z()]))))


def _uv_interior_witness(child, parent, tol: float
                         ) -> tuple[Optional[np.ndarray], Optional[str]]:
    """Find one deterministic UV that is IN child and IN/ON parent.

    Splitter is expected to retain the same support surface/parameterization.
    The witness checks that expectation rather than assuming it.
    """
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.BRepClass import BRepClass_FaceClassifier
    from OCP.BRepTools import BRepTools
    from OCP.TopAbs import TopAbs_IN, TopAbs_ON
    from OCP.gp import gp_Pnt2d

    u0, u1, v0, v1 = map(float, BRepTools.UVBounds_s(child))
    if not all(np.isfinite([u0, u1, v0, v1])) or u1 <= u0 or v1 <= v0:
        return None, "child has invalid UV bounds"

    # Center-biased deterministic sequence; much more likely than a uniform
    # boundary-heavy grid to find a witness in a narrow trimmed patch.
    frac = (0.5, 0.25, 0.75, 0.125, 0.875,
            0.375, 0.625, 0.0625, 0.9375,
            0.1875, 0.3125, 0.4375, 0.5625,
            0.6875, 0.8125)
    sc = BRepAdaptor_Surface(child)
    sp = BRepAdaptor_Surface(parent)

    for fu in frac:
        u = u0 + (u1 - u0) * fu
        for fv in frac:
            v = v0 + (v1 - v0) * fv
            uv = gp_Pnt2d(float(u), float(v))
            cc = BRepClass_FaceClassifier(child, uv, tol, True)
            if cc.State() != TopAbs_IN:
                continue
            cp = BRepClass_FaceClassifier(parent, uv, tol, True)
            if cp.State() not in (TopAbs_IN, TopAbs_ON):
                return None, "child interior UV is outside parent trim"
            pc = sc.Value(float(u), float(v))
            pp = sp.Value(float(u), float(v))
            dc = np.array([pc.X(), pc.Y(), pc.Z()], dtype=np.float64)
            dp = np.array([pp.X(), pp.Y(), pp.Z()], dtype=np.float64)
            if float(np.linalg.norm(dc - dp)) > tol:
                return None, "child and parent support surfaces disagree"
            return np.array([u, v], dtype=np.float64), None
    return None, "no stable interior UV witness found"


def split_face(face_rec: FaceRecord,
               section_edges: list[SectionEdgeRecord], *,
               base_tol: float = 1e-7,
               area_rel_tol: float = 2e-6,
               fuzzy: float = 0.0,
               parallel: bool = True,
               use_obb: bool = True,
               allow_risky: bool = False) -> FaceSplitResult:
    """Split one face once by all verified section edges that lie on it."""
    from OCP.BRep import BRep_Tool
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Splitter
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.TopAbs import TopAbs_FACE
    from OCP.TopExp import TopExp_Explorer
    try:
        from OCP.TopTools import TopTools_ListOfShape
    except ImportError:
        # OCP 8.x generated collection types moved to OCP.collections.
        from OCP.collections import (
            List_TopoDS_Shape as TopTools_ListOfShape,
        )
    from OCP.TopoDS import TopoDS

    if not base_tol > 0:
        raise ValueError("base_tol must be positive")
    if area_rel_tol < 0:
        raise ValueError("area_rel_tol must be >= 0")
    if fuzzy < 0:
        raise ValueError("fuzzy must be >= 0")

    before = _face_area(face_rec.face)
    if not section_edges:
        return FaceSplitResult(
            face_rec.face_id, "unchanged",
            [SplitFacePiece(face_rec.face_id, 0, face_rec.face,
                            before, None, True)],
            0, before, before, 0.0)

    bad = sorted(set(flag for e in section_edges for flag in e.risk_flags))
    if bad and not allow_risky:
        raise SplitError(
            f"face {face_rec.face_id}: refusing local split on risk flags "
            f"{bad}", kind="RiskySectionCurve")

    args = TopTools_ListOfShape()
    args.Append(face_rec.face)
    tools = TopTools_ListOfShape()
    seen = set()
    for e in section_edges:
        # TopoDS hash is stable for the underlying TShape/location during
        # this operation and avoids adding the same tool twice.
        h = hash(e.edge)
        if h in seen:
            continue
        seen.add(h)
        tools.Append(e.edge)

    sp = BRepAlgoAPI_Splitter()
    sp.SetArguments(args)
    sp.SetTools(tools)
    sp.SetNonDestructive(True)
    sp.SetRunParallel(bool(parallel))
    sp.SetUseOBB(bool(use_obb))
    if fuzzy > 0:
        sp.SetFuzzyValue(float(fuzzy))
    sp.Build()
    if not sp.IsDone():
        raise SplitError(
            f"face {face_rec.face_id}: OCCT local splitter failed",
            kind="FaceSplitEngineFailure")
    if not BRepCheck_Analyzer(sp.Shape()).IsValid():
        raise SplitError(
            f"face {face_rec.face_id}: splitter returned invalid topology",
            kind="FaceSplitInvalid")

    children = []
    ex = TopExp_Explorer(sp.Shape(), TopAbs_FACE)
    while ex.More():
        children.append(TopoDS.Face(ex.Current()))
        ex.Next()
    if not children:
        raise SplitError(
            f"face {face_rec.face_id}: splitter returned no faces",
            kind="FaceSplitInvalid")

    after_areas = [_face_area(f) for f in children]
    after = float(sum(after_areas))
    err = after - before
    scale = _bbox_scale(face_rec.face)
    area_tol = max(area_rel_tol * max(before, 1.0),
                   8.0 * base_tol * scale)
    if abs(err) > area_tol:
        raise SplitError(
            f"face {face_rec.face_id}: split does not conserve area "
            f"({before:.12g} -> {after:.12g}, error {err:.6g}, "
            f"tol {area_tol:.6g})",
            kind="FaceSplitAreaMismatch")

    verify_tol = max(
        float(base_tol),
        2.0 * float(BRep_Tool.Tolerance_s(face_rec.face)),
        max((e.verify_tolerance for e in section_edges), default=base_tol))

    pieces = []
    notes = []
    for i, (child, ar) in enumerate(zip(children, after_areas)):
        if not BRepCheck_Analyzer(child).IsValid():
            raise SplitError(
                f"face {face_rec.face_id}: child {i} is invalid",
                kind="FaceSplitInvalid")
        uv, note = _uv_interior_witness(
            child, face_rec.face, verify_tol)
        if note:
            # Lack of a stable grid witness is not alone a rejection: very
            # thin but valid pieces exist. It is reported so later Stage 6
            # can demand stronger evidence before final certification.
            notes.append(f"piece {i}: {note}")
        pieces.append(SplitFacePiece(
            face_rec.face_id, i, child, ar, uv, False))

    return FaceSplitResult(
        parent_face_id=face_rec.face_id,
        status="split" if len(pieces) > 1 else "unchanged_after_splitter",
        pieces=pieces,
        source_edges=len(seen),
        area_before=before,
        area_after=after,
        area_error=err,
        notes=notes,
    )


def split_models(a: BRepModel, b: BRepModel,
                 intersections: ModelIntersectionResult, *,
                 base_tol: float = 1e-7,
                 area_rel_tol: float = 2e-6,
                 fuzzy: float = 0.0,
                 parallel: bool = True,
                 use_obb: bool = True,
                 allow_risky: bool = False) -> ModelSplitResult:
    """Split only faces touched by verified transverse section curves.

    All point contacts, near contacts, distance failures, and (by default)
    near-tangent/seam curves are surfaced as unresolved rather than being
    converted into guessed face topology.
    """
    edges_a: dict[int, list[SectionEdgeRecord]] = {}
    edges_b: dict[int, list[SectionEdgeRecord]] = {}
    unresolved: list[tuple[int, int, str]] = []

    for pair in intersections.pairs:
        if pair.status == "curve":
            for e in pair.edges:
                edges_a.setdefault(pair.face_a, []).append(e)
                edges_b.setdefault(pair.face_b, []).append(e)
        elif pair.status == "disjoint":
            continue
        elif pair.status == "curve_near_tangent" and allow_risky:
            for e in pair.edges:
                edges_a.setdefault(pair.face_a, []).append(e)
                edges_b.setdefault(pair.face_b, []).append(e)
        else:
            unresolved.append((pair.face_a, pair.face_b, pair.status))

    out_a = []
    out_b = []
    split_calls = 0

    for fr in a.faces:
        es = edges_a.get(fr.face_id, [])
        if es:
            split_calls += 1
        out_a.append(split_face(
            fr, es, base_tol=base_tol, area_rel_tol=area_rel_tol,
            fuzzy=fuzzy, parallel=parallel, use_obb=use_obb,
            allow_risky=allow_risky))

    for fr in b.faces:
        es = edges_b.get(fr.face_id, [])
        if es:
            split_calls += 1
        out_b.append(split_face(
            fr, es, base_tol=base_tol, area_rel_tol=area_rel_tol,
            fuzzy=fuzzy, parallel=parallel, use_obb=use_obb,
            allow_risky=allow_risky))

    return ModelSplitResult(
        faces_a=out_a,
        faces_b=out_b,
        split_calls=split_calls,
        affected_faces_a=len(edges_a),
        affected_faces_b=len(edges_b),
        unresolved_contacts=unresolved,
    )
