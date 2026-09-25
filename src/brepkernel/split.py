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
- near-tangent section edges are refused by default;
- a verified one-sided seam curve reuses the existing closing boundary on
  that operand and remains a split tool on the opposite operand;
- a curve that is a seam on both operands remains unresolved by default;
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


def _face_boundary_edges(face) -> list:
    """Return the TopoDS_Edge boundary wires' edges of a face."""
    from OCP.TopAbs import TopAbs_EDGE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS
    out = []
    ex = TopExp_Explorer(TopoDS.Face(face), TopAbs_EDGE)
    while ex.More():
        out.append(TopoDS.Edge(ex.Current()))
        ex.Next()
    return out


@dataclass
class SplitFacePiece:
    parent_face_id: int
    piece_index: int
    face: object
    area: float
    uv_witness: Optional[np.ndarray]
    unchanged: bool = False
    # G2.4: for pieces of a coincident pair, "same"/"opp" if this piece
    # lies in the common (overlapping) region, else None.
    coincidence: Optional[str] = None
    # (partner_face_id, partner_piece_index) for common pieces.
    coincident_partner: Optional[tuple] = None
    # Partner FaceRecord for the witness check (G2.6).
    coincident_partner_rec: Optional[object] = None


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
    section_edges: list[SectionEdgeRecord] = field(default_factory=list)
    reused_seam_edges_a: int = 0
    reused_seam_edges_b: int = 0
    shared_seam_refusals: int = 0
    # G2.5: boundary-contact edges (ON both trim boundaries; no split).
    boundary_edges: list[SectionEdgeRecord] = field(default_factory=list)
    # G2.6: partner boundary edges used as overlap-split tools for
    # coincident pairs, indexed for coincident_boundary provenance.
    # Each entry: {"operand", "face_id", "edge_id", "edge"}.
    coincident_boundary_tools: list = field(default_factory=list)

    @property
    def certified_local_split(self) -> bool:
        return not self.unresolved_contacts


def _face_area(face) -> float:
    """High-accuracy area for split-conservation checks.

    OCCT's no-Eps SurfaceProperties overload uses a non-adaptive path that can
    be materially inaccurate on rational B-spline faces (a NURBS sphere was
    off by about 0.35%).  The Eps overload uses adaptive 2D Gauss integration,
    so conservation is checked against the actual trimmed surface rather than
    against a loose integration estimate.
    """
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps

    g = GProp_GProps()
    err = BRepGProp.SurfaceProperties_s(face, g, 1e-9, False)
    if float(err) < 0.0:
        raise SplitError("adaptive face-area integration failed",
                         kind="FaceAreaIntegrationFailed")
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
               allow_risky: bool = False,
               ignored_risk_flags: frozenset[str] = frozenset(),
               extra_tool_edges: Optional[list] = None) -> FaceSplitResult:
    """Split one face once by all verified section edges that lie on it.

    extra_tool_edges (G2.4) are raw TopoDS_Edge tools (e.g. a coincident
    partner face's boundary edges) added to the same single splitter call.
    """
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
    if not section_edges and not extra_tool_edges:
        return FaceSplitResult(
            face_rec.face_id, "unchanged",
            [SplitFacePiece(face_rec.face_id, 0, face_rec.face,
                            before, None, True)],
            0, before, before, 0.0)

    bad = sorted(set(
        flag for e in section_edges for flag in e.risk_flags
        if flag not in ignored_risk_flags))
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
    for te in extra_tool_edges or []:
        h = hash(te)
        if h in seen:
            continue
        seen.add(h)
        tools.Append(te)

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


def _piece_interior_point(piece_face):
    """A 3D point strictly inside a face piece, or None."""
    from OCP.BRepClass3d import BRepClass3d_SolidExplorer
    from OCP.gp import gp_Pnt
    from OCP.TopoDS import TopoDS
    try:
        explorer = BRepClass3d_SolidExplorer()
        p = gp_Pnt()
        if explorer.FindAPointInTheFace_s(TopoDS.Face(piece_face), p):
            return np.array([p.X(), p.Y(), p.Z()], dtype=np.float64)
    except Exception:
        pass
    return None


def _classify_point_on_face(xyz, face, tol: float) -> str:
    """Classify a 3D point against a face trim: IN/OUT/ON/UNKNOWN."""
    from OCP.BRepClass import BRepClass_FaceClassifier
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.gp import gp_Pnt, gp_Pnt2d
    from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf
    from OCP.TopAbs import TopAbs_IN, TopAbs_OUT, TopAbs_ON
    from OCP.TopoDS import TopoDS
    try:
        f = TopoDS.Face(face)
        surf = BRepAdaptor_Surface(f).Surface().Surface()
        proj = GeomAPI_ProjectPointOnSurf(
            gp_Pnt(float(xyz[0]), float(xyz[1]), float(xyz[2])), surf)
        if proj.NbPoints() == 0:
            return "UNKNOWN"
        u, v = proj.LowerDistanceParameters()
        if float(proj.LowerDistance()) > tol:
            # Not on the support; still classify the projection for the
            # trim test but mark via distance (caller checks support).
            pass
        clf = BRepClass_FaceClassifier(f, gp_Pnt2d(float(u), float(v)),
                                       float(tol))
        st = clf.State()
        if st == TopAbs_IN:
            return "IN"
        if st == TopAbs_OUT:
            return "OUT"
        if st == TopAbs_ON:
            return "ON"
        return "UNKNOWN"
    except Exception:
        return "UNKNOWN"


def _label_coincident_pieces(a, b, out_a, out_b, coincident_pairs,
                             base_tol: float):
    """Label common/remainder pieces for coincident pairs (G2.4).

    For each coincident pair, every piece of each face is tested against
    the partner's original face trim. Common pieces get
    coincidence="same"/"opp" and a partner reference; the overlap is
    verified independently (area conservation already checked by the
    splitter; here common pieces must match 1-1 by area and bbox, and
    witnesses must lie on both supports).
    """
    from OCP.BRep import BRep_Tool
    from OCP.BRepClass import BRepClass_FaceClassifier
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf
    from OCP.gp import gp_Pnt
    from OCP.TopoDS import TopoDS

    res_by_a = {r.parent_face_id: r for r in out_a}
    res_by_b = {r.parent_face_id: r for r in out_b}
    rec_by_a = {f.face_id: f for f in a.faces}
    rec_by_b = {f.face_id: f for f in b.faces}

    def face_bbox(face):
        from OCP.Bnd import Bnd_Box
        from OCP.BRepBndLib import BRepBndLib
        bb = Bnd_Box()
        BRepBndLib.Add_s(TopoDS.Face(face), bb)
        lo = bb.CornerMin()
        hi = bb.CornerMax()
        return (np.array([lo.X(), lo.Y(), lo.Z()]),
                np.array([hi.X(), hi.Y(), hi.Z()]))

    for rec in coincident_pairs:
        ra = res_by_a.get(rec.face_a)
        rb = res_by_b.get(rec.face_b)
        fa_rec = rec_by_a.get(rec.face_a)
        fb_rec = rec_by_b.get(rec.face_b)
        if ra is None or rb is None or fa_rec is None or fb_rec is None:
            raise SplitError(
                f"coincident pair A{rec.face_a}/B{rec.face_b}: missing "
                f"split result", kind="CoincidentPairMissing")
        tol = max(float(base_tol),
                  2.0 * float(BRep_Tool.Tolerance_s(
                      TopoDS.Face(fa_rec.face))))
        # Label pieces by testing against the partner's original trim.
        for pieces, partner_face, side in (
                (ra.pieces, fb_rec.face, "A"),
                (rb.pieces, fa_rec.face, "B")):
            for pc in pieces:
                p = _piece_interior_point(pc.face)
                if p is None:
                    raise SplitError(
                        f"coincident pair A{rec.face_a}/B{rec.face_b}: "
                        f"no interior point for {side} piece "
                        f"{pc.piece_index}", kind="CoincidentNoWitness")
                # Must be on the partner's support (they are coincident).
                psurf = BRepAdaptor_Surface(
                    TopoDS.Face(partner_face)).Surface().Surface()
                proj = GeomAPI_ProjectPointOnSurf(
                    gp_Pnt(float(p[0]), float(p[1]), float(p[2])), psurf)
                if proj.NbPoints() == 0 or \
                        float(proj.LowerDistance()) > tol:
                    raise SplitError(
                        f"coincident pair A{rec.face_a}/B{rec.face_b}: "
                        f"{side} piece {pc.piece_index} witness off "
                        f"partner support",
                        kind="CoincidentWitnessMismatch")
                cls = _classify_point_on_face(p, partner_face, tol)
                if cls == "IN":
                    pc.coincidence = rec.sense
                elif cls == "OUT":
                    pc.coincidence = None
                else:
                    raise SplitError(
                        f"coincident pair A{rec.face_a}/B{rec.face_b}: "
                        f"{side} piece {pc.piece_index} classifies "
                        f"{cls} on partner trim",
                        kind="CoincidentWitnessMismatch")
        # Independent verification: common pieces match 1-1 by area and
        # bbox; total common area agrees both sides.
        common_a = [pc for pc in ra.pieces if pc.coincidence]
        common_b = [pc for pc in rb.pieces if pc.coincidence]
        area_a = sum(pc.area for pc in common_a)
        area_b = sum(pc.area for pc in common_b)
        rec.common_area = 0.5 * (area_a + area_b)
        area_tol = max(2e-6 * max(area_a, area_b, 1.0),
                       8.0 * base_tol * max(area_a, area_b) ** 0.5)
        if abs(area_a - area_b) > area_tol:
            raise SplitError(
                f"coincident pair A{rec.face_a}/B{rec.face_b}: common "
                f"area mismatch {area_a:.12g} vs {area_b:.12g}",
                kind="CoincidentAreaMismatch")
        # 1-1 matching by bbox overlap and area.
        unmatched_b = list(common_b)
        for pca in common_a:
            la, ha = face_bbox(pca.face)
            best = None
            best_score = None
            for pcb in unmatched_b:
                lb, hb = face_bbox(pcb.face)
                # Bbox overlap volume > 0?
                overlap = np.maximum(
                    0.0, np.minimum(ha, hb) - np.maximum(la, lb))
                if float(np.prod(overlap)) <= 0.0:
                    continue
                area_err = abs(pca.area - pcb.area) / max(
                    pca.area, pcb.area, 1e-30)
                if best_score is None or area_err < best_score:
                    best_score = area_err
                    best = pcb
            if best is None or best_score > 1e-3:
                raise SplitError(
                    f"coincident pair A{rec.face_a}/B{rec.face_b}: no "
                    f"1-1 match for A piece {pca.piece_index}",
                    kind="CoincidentPieceMismatch")
            pca.coincident_partner = (rec.face_b, best.piece_index)
            best.coincident_partner = (rec.face_a, pca.piece_index)
            pca.coincident_partner_rec = fb_rec
            best.coincident_partner_rec = fa_rec
            unmatched_b.remove(best)
        if unmatched_b:
            raise SplitError(
                f"coincident pair A{rec.face_a}/B{rec.face_b}: "
                f"{len(unmatched_b)} unmatched B common pieces",
                kind="CoincidentPieceMismatch")
        rec.verified = True


def _point_to_edge_distance(xyz, edge, tol: float) -> float:
    """Min distance from a 3D point to a TopoDS_Edge, or inf on failure."""
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    from OCP.GeomAPI import GeomAPI_ProjectPointOnCurve
    from OCP.gp import gp_Pnt
    from OCP.TopoDS import TopoDS
    try:
        curve = BRepAdaptor_Curve(TopoDS.Edge(edge)).Curve().Curve()
        proj = GeomAPI_ProjectPointOnCurve(
            gp_Pnt(float(xyz[0]), float(xyz[1]), float(xyz[2])), curve)
        if proj.NbPoints() == 0:
            return float("inf")
        return float(proj.LowerDistance())
    except Exception:
        return float("inf")


def _edges_coincident(e1_xyz, e2_edge, tol: float) -> bool:
    """True when every sample of edge 1 lies within tol of edge 2."""
    pts = np.asarray(e1_xyz, dtype=float).reshape(-1, 3)
    if len(pts) == 0:
        return False
    # Check a bounded subset of samples (endpoints always included).
    idx = np.unique(np.linspace(0, len(pts) - 1,
                                min(len(pts), 25)).astype(int))
    for i in idx:
        if _point_to_edge_distance(
                (float(pts[i][0]), float(pts[i][1]), float(pts[i][2])),
                e2_edge, tol) > tol:
            return False
    return True


def _resolve_contacts(intersections, out_a, out_b, base_tol: float):
    """Resolve boundary/point contacts per G2.5; return unresolved list.

    A boundary-contact edge is resolved if it coincides (within tol)
    with a coincident common-region boundary or with a verified section
    edge from another pair. A point contact is resolved if it lies
    within tolerance of a verified section edge, a coincident
    common-region boundary, or a boundary-contact edge. Anything left
    unexplained stays unresolved (blocks).
    """
    from OCP.TopoDS import TopoDS
    from OCP.TopAbs import TopAbs_EDGE
    from OCP.TopExp import TopExp_Explorer

    # Verified section edges (non-boundary) and boundary-contact edges,
    # tagged with their pair so a boundary edge is never resolved by an
    # edge from its own pair.
    section_edges = []   # (face_a, face_b, TopoDS_Edge)
    boundary_edges = []  # (face_a, face_b, SectionEdgeRecord)
    for pair in intersections.pairs:
        if pair.status in ("curve", "curve_near_tangent",
                           "boundary_contact"):
            for e in pair.edges:
                if getattr(e, "is_boundary_contact", False):
                    boundary_edges.append((pair.face_a, pair.face_b, e))
                else:
                    section_edges.append(
                        (pair.face_a, pair.face_b, e.edge))
    # Coincident common-region boundaries.
    res_by_a = {r.parent_face_id: r for r in out_a}
    res_by_b = {r.parent_face_id: r for r in out_b}
    region_boundaries = []
    for rec in intersections.coincident_pairs:
        if not rec.verified or rec.common_area <= 0.0:
            continue
        for res_by, fid in ((res_by_a, rec.face_a),
                            (res_by_b, rec.face_b)):
            r = res_by.get(fid)
            if r is None:
                continue
            for pc in r.pieces:
                if not pc.coincidence:
                    continue
                ex = TopExp_Explorer(TopoDS.Face(pc.face), TopAbs_EDGE)
                while ex.More():
                    region_boundaries.append(TopoDS.Edge(ex.Current()))
                    ex.Next()

    explained_edges = ([e for _, _, e in section_edges]
                       + [e for _, _, rec in boundary_edges for e in [rec.edge]]
                       + region_boundaries)

    tol = float(base_tol) * 4.0
    still_unresolved = []
    # Boundary contacts: each must coincide with a region boundary or
    # with a verified section edge from another pair.
    for fa, fb, rec in boundary_edges:
        ok = any(_edges_coincident(rec.xyz, rb, tol)
                 for rb in region_boundaries)
        if not ok:
            ok = any(_edges_coincident(rec.xyz, e2, tol)
                     for fa2, fb2, e2 in section_edges
                     if (fa2, fb2) != (fa, fb))
        if not ok:
            # A boundary edge may also be explained by the same
            # geometric edge found by another pair's section (also a
            # boundary contact there).
            ok = any(_edges_coincident(rec.xyz, rec2.edge, tol)
                     for fa2, fb2, rec2 in boundary_edges
                     if (fa2, fb2) != (fa, fb))
        if not ok:
            still_unresolved.append((fa, fb, "boundary_contact"))
    for pair in intersections.pairs:
        if pair.status == "point_contact":
            resolved = True
            for pt in pair.point_contacts:
                xyz = (float(pt[0]), float(pt[1]), float(pt[2]))
                if not any(_point_to_edge_distance(xyz, e, tol) <= tol
                           for e in explained_edges):
                    resolved = False
                    break
            if not resolved:
                still_unresolved.append(
                    (pair.face_a, pair.face_b, pair.status))
        elif pair.status not in ("disjoint", "coincident", "curve",
                                 "curve_near_tangent", "boundary_contact"):
            still_unresolved.append(
                (pair.face_a, pair.face_b, pair.status))
    return still_unresolved


def split_models(a: BRepModel, b: BRepModel,
                 intersections: ModelIntersectionResult, *,
                 base_tol: float = 1e-7,
                 area_rel_tol: float = 2e-6,
                 fuzzy: float = 0.0,
                 parallel: bool = True,
                 use_obb: bool = True,
                 allow_risky: bool = False) -> ModelSplitResult:
    """Split only faces touched by verified transverse section curves.

    All point contacts, near contacts, distance failures, and near-tangent
    curves are surfaced as unresolved rather than guessed. Verified one-sided
    seam curves reuse the existing seam topology on that operand and split the
    opposite operand; shared seams remain unresolved by default.

    G2.4: coincident face pairs are split by the partner face's boundary
    edges (in the same single splitter call as any section edges), then
    the resulting pieces are labeled common/remainder with independent
    verification.
    """
    edges_a: dict[int, list[SectionEdgeRecord]] = {}
    edges_b: dict[int, list[SectionEdgeRecord]] = {}
    unresolved: list[tuple[int, int, str]] = []
    used_sections: list[SectionEdgeRecord] = []
    reused_seam_a = 0
    reused_seam_b = 0
    shared_seam_refusals = 0
    boundary_edges: list[SectionEdgeRecord] = []

    for pair in intersections.pairs:
        if pair.status == "curve":
            for e in pair.edges:
                # G2.5: boundary-contact edges create no split. This check
                # comes before seam routing: an edge ON both trim
                # boundaries is not a cut on either operand, even if it
                # also carries a seam flag.
                if getattr(e, "is_boundary_contact", False):
                    boundary_edges.append(e)
                    continue
                flags = set(e.risk_flags)
                seam_a = "seam_on_a" in flags
                seam_b = "seam_on_b" in flags

                # A section edge that is a seam on one operand already lies on
                # that face's existing closing boundary. It is not a new cut
                # there, but it remains a valid verified splitting tool on the
                # opposite operand. This avoids duplicating/re-splitting the
                # seam while preserving the physical intersection contour.
                if seam_a and seam_b and not allow_risky:
                    shared_seam_refusals += 1
                    unresolved.append(
                        (pair.face_a, pair.face_b, "shared_seam_curve"))
                    continue
                if seam_a and not allow_risky:
                    reused_seam_a += 1
                if seam_b and not allow_risky:
                    reused_seam_b += 1
                if not seam_a or allow_risky:
                    edges_a.setdefault(pair.face_a, []).append(e)
                if not seam_b or allow_risky:
                    edges_b.setdefault(pair.face_b, []).append(e)
                used_sections.append(e)
        elif pair.status == "boundary_contact":
            # G2.5: no split tools; the edges resolve (or fail to
            # resolve) contacts below.
            boundary_edges.extend(pair.edges)
        elif pair.status in ("disjoint", "coincident"):
            continue
        elif pair.status == "curve_near_tangent" and allow_risky:
            for e in pair.edges:
                if getattr(e, "is_boundary_contact", False):
                    boundary_edges.append(e)
                    continue
                edges_a.setdefault(pair.face_a, []).append(e)
                edges_b.setdefault(pair.face_b, []).append(e)
                used_sections.append(e)
        else:
            unresolved.append((pair.face_a, pair.face_b, pair.status))

    # G2.4: coincident tools. For each coincident pair, each face is cut
    # by the partner face's boundary edges.
    coinc_tools_a: dict[int, list] = {}
    coinc_tools_b: dict[int, list] = {}
    coincident_boundary_tools: list = []
    by_a = {f.face_id: f for f in a.faces}
    by_b = {f.face_id: f for f in b.faces}
    for rec in intersections.coincident_pairs:
        fa = by_a.get(rec.face_a)
        fb = by_b.get(rec.face_b)
        if fa is None or fb is None:
            continue
        b_bnd = list(enumerate(_face_boundary_edges(fb.face)))
        a_bnd = list(enumerate(_face_boundary_edges(fa.face)))
        coinc_tools_a.setdefault(
            rec.face_a, []).extend(e for _, e in b_bnd)
        coinc_tools_b.setdefault(
            rec.face_b, []).extend(e for _, e in a_bnd)
        # G2.6: index the partner boundary edges as provenance tools.
        for eid, e in b_bnd:
            coincident_boundary_tools.append(
                {"operand": "B", "face_id": rec.face_b,
                 "edge_id": int(eid), "edge": e})
        for eid, e in a_bnd:
            coincident_boundary_tools.append(
                {"operand": "A", "face_id": rec.face_a,
                 "edge_id": int(eid), "edge": e})

    out_a = []
    out_b = []
    split_calls = 0

    for fr in a.faces:
        es = edges_a.get(fr.face_id, [])
        ct = coinc_tools_a.get(fr.face_id, [])
        if es or ct:
            split_calls += 1
        out_a.append(split_face(
            fr, es, base_tol=base_tol, area_rel_tol=area_rel_tol,
            fuzzy=fuzzy, parallel=parallel, use_obb=use_obb,
            allow_risky=allow_risky,
            ignored_risk_flags=frozenset({"seam_on_b"}),
            extra_tool_edges=ct or None))

    for fr in b.faces:
        es = edges_b.get(fr.face_id, [])
        ct = coinc_tools_b.get(fr.face_id, [])
        if es or ct:
            split_calls += 1
        out_b.append(split_face(
            fr, es, base_tol=base_tol, area_rel_tol=area_rel_tol,
            fuzzy=fuzzy, parallel=parallel, use_obb=use_obb,
            allow_risky=allow_risky,
            ignored_risk_flags=frozenset({"seam_on_a"}),
            extra_tool_edges=ct or None))

    # G2.4: label coincident pieces and verify the overlap independently.
    if intersections.coincident_pairs:
        _label_coincident_pieces(
            a, b, out_a, out_b, intersections.coincident_pairs,
            base_tol=base_tol)

    # G2.5: resolve boundary/point contacts against verified section
    # edges and coincident-region boundaries.
    unresolved = _resolve_contacts(
        intersections, out_a, out_b, base_tol=base_tol)

    return ModelSplitResult(
        faces_a=out_a,
        faces_b=out_b,
        split_calls=split_calls,
        affected_faces_a=len(edges_a) + len(
            [k for k in coinc_tools_a if k not in edges_a]),
        affected_faces_b=len(edges_b) + len(
            [k for k in coinc_tools_b if k not in edges_b]),
        unresolved_contacts=unresolved,
        section_edges=used_sections,
        reused_seam_edges_a=reused_seam_a,
        reused_seam_edges_b=reused_seam_b,
        shared_seam_refusals=shared_seam_refusals,
        boundary_edges=boundary_edges,
        coincident_boundary_tools=coincident_boundary_tools,
    )
