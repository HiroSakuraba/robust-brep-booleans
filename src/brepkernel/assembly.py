"""Global B-rep patch classification and shell/solid assembly.

This stage turns locally split trimmed faces back into a material result.

Rules:
- classify each patch from an interior point on the exact B-rep face;
- classify against the other operand's original OCCT solids;
- refuse boundary/unknown witnesses rather than majority-voting;
- unresolved tangencies / near contacts block assembly;
- sew selected faces with OCCT and require closed manifold shells;
- normalize and geometrically nest shells so cavities and disconnected
  components are represented explicitly;
- require OCCT-valid, positively oriented result solids.

The mesh path may later accelerate classification, but it is not the source of
truth here.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .freeform import FreeformError
from .split import ModelSplitResult
from .step_ingest import BRepModel


class AssemblyError(FreeformError):
    pass


@dataclass
class PatchDecision:
    operand: str
    parent_face_id: int
    piece_index: int
    classification: str
    keep: bool
    reverse_for_difference: bool
    witness_xyz: np.ndarray
    witness_xyz_all: np.ndarray
    witness_classifications: tuple[str, ...]
    source_face: object
    selected_face: Optional[object] = None
    sewed_face: Optional[object] = None
    # G12b: if set, this decision was propagated from another face's
    # classification (the representative of an untouched region).
    # The value is the parent_face_id of the representative.
    propagated_from: Optional[int] = None


@dataclass
class ShellAssemblyRecord:
    shell_index: int
    shell: object
    normalized_shell: object
    temp_solid: object
    volume: float
    interior_point: np.ndarray
    parent_shell: Optional[int]
    depth: int


@dataclass
class SolidAssemblyRecord:
    solid_index: int
    solid: object
    outer_shell: int
    cavity_shells: tuple[int, ...]
    volume: float


@dataclass
class EdgeLineageRecord:
    result_edge_index: int
    piece_refs: tuple[tuple[str, int, int], ...]
    parent_faces: tuple[tuple[str, int], ...]
    operands: tuple[str, ...]
    intersection_refs: tuple[tuple[int, int, int], ...]
    source_boundary_refs: tuple[tuple[str, int], ...]
    provenance_kind: str
    verified_pcurves: bool


@dataclass
class SectionPayloadRecord:
    face_a: int
    face_b: int
    section_edge_index: int
    parameters: np.ndarray
    xyz: np.ndarray
    uv_a: np.ndarray
    uv_b: np.ndarray
    edge_tolerance: float
    verify_tolerance: float
    max_surface_error_a: float
    max_surface_error_b: float
    max_cross_surface_error: float
    min_transversality: float
    max_transversality: float
    risk_flags: tuple[str, ...]
    repaired_same_parameter: bool
    exact_curve_on_surface_checked: bool
    exact_surface_error_a: Optional[float]
    exact_surface_error_b: Optional[float]
    shadow_crosschecked: bool
    shadow_max_distance: Optional[float]
    shadow_length_rel_error: Optional[float]
    result_edge_indices: tuple[int, ...]


@dataclass
class BooleanAssemblyResult:
    operation: str
    decisions: list[PatchDecision]
    selected_faces: int
    sewed_shape: object
    shells: list[ShellAssemblyRecord]
    solids: list[SolidAssemblyRecord]
    shape: object
    volume: float
    free_edges: int
    multiple_edges: int
    edge_lineage: list[EdgeLineageRecord] = field(default_factory=list)
    section_payloads: list[SectionPayloadRecord] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    # G12b: untouched-region classification stats
    region_stats: dict = field(default_factory=dict)

    @property
    def is_empty(self) -> bool:
        return self.selected_faces == 0


def _p3(p) -> np.ndarray:
    return np.array([p.X(), p.Y(), p.Z()], dtype=np.float64)


def _bbox(shape) -> tuple[np.ndarray, np.ndarray]:
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib

    b = Bnd_Box()
    BRepBndLib.AddOptimal_s(shape, b, False, False)
    if b.IsVoid():
        z = np.zeros(3, dtype=np.float64)
        return z, z
    p0, p1 = b.CornerMin(), b.CornerMax()
    return _p3(p0), _p3(p1)


def _face_points(face, tol: float, *, max_points: int = 7,
                 min_points: int = 3) -> np.ndarray:
    """Find several deterministic points strictly inside a trimmed face.

    One centroid-like witness can silently misclassify a patch that still
    straddles the other solid because of a missed or degenerate split.  The
    Tier B/C path therefore samples several well-separated *face-interior*
    points and requires them all to have the same material classification.

    Failure to find enough stable interior witnesses is a refusal, not
    permission to fall back to one-point classification.
    """
    from OCP.BRepClass3d import BRepClass3d_SolidExplorer
    from OCP.gp import gp_Pnt

    if max_points < min_points or min_points < 1:
        raise ValueError("invalid face witness count")

    candidates: list[np.ndarray] = []
    p = gp_Pnt()
    try:
        if BRepClass3d_SolidExplorer.FindAPointInTheFace_s(face, p):
            candidates.append(_p3(p))
    except Exception:
        pass

    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.BRepClass import BRepClass_FaceClassifier
    from OCP.BRepTools import BRepTools
    from OCP.TopAbs import TopAbs_IN
    from OCP.gp import gp_Pnt2d

    u0, u1, v0, v1 = map(float, BRepTools.UVBounds_s(face))
    if not all(np.isfinite([u0, u1, v0, v1])) or u1 <= u0 or v1 <= v0:
        raise AssemblyError("cannot obtain finite face UV bounds",
                            kind="NoInteriorFaceWitness")
    surf = BRepAdaptor_Surface(face)

    # Center-biased low-discrepancy-ish sequence.  Pairing the U/V fractions
    # through two permutations avoids filling only one diagonal of UV space.
    fu = (0.5, 0.25, 0.75, 0.125, 0.875,
          0.375, 0.625, 0.0625, 0.9375,
          0.1875, 0.6875, 0.4375, 0.8125,
          0.3125, 0.5625)
    fv = (0.5, 0.75, 0.25, 0.375, 0.625,
          0.875, 0.125, 0.6875, 0.3125,
          0.9375, 0.4375, 0.1875, 0.5625,
          0.8125, 0.0625)

    lo, hi = _bbox(face)
    scale = max(float(np.linalg.norm(hi - lo)), 1.0)
    sep = max(8.0 * float(tol), 1e-10 * scale)

    def add_point(x: np.ndarray):
        if all(float(np.linalg.norm(x - q)) > sep for q in candidates):
            candidates.append(x)

    # First use the paired sequence, then a small Cartesian fallback for thin
    # or oddly trimmed regions.
    for a, b in zip(fu, fv):
        u = u0 + (u1 - u0) * a
        v = v0 + (v1 - v0) * b
        cl = BRepClass_FaceClassifier(
            face, gp_Pnt2d(float(u), float(v)), float(tol), True)
        if cl.State() == TopAbs_IN:
            add_point(_p3(surf.Value(float(u), float(v))))
            if len(candidates) >= max_points:
                break

    if len(candidates) < min_points:
        grid = (0.2, 0.4, 0.6, 0.8)
        for a in grid:
            if len(candidates) >= max_points:
                break
            for b in grid:
                u = u0 + (u1 - u0) * a
                v = v0 + (v1 - v0) * b
                cl = BRepClass_FaceClassifier(
                    face, gp_Pnt2d(float(u), float(v)), float(tol), True)
                if cl.State() == TopAbs_IN:
                    add_point(_p3(surf.Value(float(u), float(v))))
                    if len(candidates) >= max_points:
                        break

    if len(candidates) < min_points:
        raise AssemblyError(
            f"only {len(candidates)} stable interior witness(es) found; "
            f"{min_points} required",
            kind="InsufficientPatchWitnesses")
    return np.vstack(candidates[:max_points])

def _occt_point_verdict(point: np.ndarray, solids: list,
                        tol: float) -> str:
    """Raw BRepClass3d_SolidClassifier verdict against a list of solids."""
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.TopAbs import TopAbs_IN, TopAbs_ON, TopAbs_OUT
    from OCP.gp import gp_Pnt

    p = gp_Pnt(float(point[0]), float(point[1]), float(point[2]))
    saw_on = False
    for solid in solids:
        c = BRepClass3d_SolidClassifier(solid)
        c.Perform(p, float(tol))
        st = c.State()
        if st == TopAbs_IN:
            return "inside"
        if st == TopAbs_ON:
            saw_on = True
        elif st != TopAbs_OUT:
            return "unknown"
    return "boundary" if saw_on else "outside"


def _classify_point_in_model(point: np.ndarray, model: BRepModel,
                             tol: float) -> str:
    """Classify a point against the union of the model's OCCT solids."""
    if not model.solids:
        raise AssemblyError(
            "material patch classification requires closed OCCT solids; "
            "shell-only input is unsupported at this stage",
            kind="ShellOnlyClassificationUnsupported")
    return _occt_point_verdict(
        point, [sr.solid for sr in model.solids], tol)


# ---------------------------------------------------------------------------
# G5: second independent point classifier (multi-ray parity).
#
# Review finding F4: BRepClass3d_SolidClassifier returned a false IN for a
# point more than 1.0 from both input surfaces, and the kernel's patch and
# shell classification depended on that same classifier.  The second
# classifier below never calls the OCCT solid classifier; it casts rays with
# IntCurvesFace_ShapeIntersector (the face/ray intersector, a different OCCT
# subsystem with a different failure mode) and counts transverse surface
# crossings per solid: odd = inside, even = outside.  Per-solid parity
# (rather than one count over the union) keeps cavity shells correct: a
# point in a cavity sees an even crossing count and classifies outside.
# ---------------------------------------------------------------------------

# Fixed, deterministic, non-axis-aligned, well-spread ray directions.
# Direction choice, documented:
# - deterministic (no RNG): identical inputs give identical verdicts, and
#   fault-injection tests are reproducible;
# - no zero components: a ray can never run parallel to a coordinate plane,
#   the common degenerate case for axis-aligned CAD geometry;
# - components drawn from {1, 2, 3} with varied signs: keeps every direction
#   away from face diagonals of axis-aligned boxes while spreading the set
#   across all octants;
# - 12 candidates: up to 7 degenerate rays (edge/vertex hits, grazing hits,
#   near-origin hits) can be discarded while still seating the 5 valid rays
#   the verdict requires (minimum 3 to decide at all).
_RAY_DIRECTIONS = tuple(
    _d / np.linalg.norm(_d)
    for _d in (
        (1, 2, 3), (-2, 1, 3), (3, -1, 2), (1, -3, -2),
        (-3, 2, -1), (2, 3, -1), (-1, -2, 3), (3, 1, -2),
        (-1, 3, -2), (2, -1, -3), (-3, -2, 1), (1, -2, -3),
    )
)

_N_RAYS_TARGET = 5   # valid rays seated before deciding
_N_RAYS_MIN = 3      # minimum valid agreeing rays for a verdict
# |n.d| below this: grazing hit, ray discarded.  A near-tangent pass is
# the one configuration where the intersector demonstrably drops a
# crossing: observed on a NURBS sphere, a ray with a single reported hit
# at |n.d| = 0.023 whose paired near-tangent hit (0.045 further along the
# ray) was missed, flipping parity from even to odd.  The found hit of a
# near-tangent pass is itself near-tangent, so discarding the ray on the
# found hit's |n.d| removes the whole failure mode; 0.05 (~3 degrees)
# keeps about 2x margin over the observed case.  Discarding is always the
# safe direction: 12 fixed directions feed a 3-ray minimum.
_TANGENT_COS = 0.05
_RAY_PMAX = 1e100


class _MultiRayClassifier:
    """Independent point-in-solid classifier via ray parity.

    Crossing-count rule: for one ray, each solid contributes its transverse
    surface-crossing count; the point is inside a solid iff that count is
    odd, and inside the model iff inside at least one solid.  A ray is
    discarded (never counted) when any hit is degenerate:
    - the hit parameter is within tol of the ray origin (too close to call);
    - the hit is within tol of an edge or vertex of any solid;
    - the ray grazes the face (|oriented-normal . direction| < 0.05).
    A discarded ray is replaced with the next fixed direction.  Fewer than
    3 valid rays, or disagreement among the seated rays, yields "unknown"
    rather than a guessed verdict.

    Rays are cast bidirectionally: for each direction d, both +d and -d
    are cast, and the pair is seated only if every solid's (+d count +
    -d count) is even.  A line meets a closed solid in an even number of
    transverse crossings, so an odd sum proves the intersector missed or
    added a crossing on that line (observed failure: a near-tangent ray
    whose missed partner crossing flips the parity).  The pair's verdict
    comes from the +d counts.  Either side degenerate, or any solid with
    an odd bidirectional sum, discards the pair.
    """

    def __init__(self, solids: list, tol: float):
        from OCP.IntCurvesFace import IntCurvesFace_ShapeIntersector

        self._tol = float(tol)
        self._solids = list(solids)
        self._intersectors = []
        for solid in self._solids:
            inter = IntCurvesFace_ShapeIntersector()
            inter.Load(solid, self._tol)
            self._intersectors.append(inter)
        self._edges = self._edge_compound(self._solids)

    @staticmethod
    def _edge_compound(solids):
        """One compound holding every edge, for edge-proximity checks."""
        from OCP.BRep import BRep_Builder
        from OCP.TopAbs import TopAbs_EDGE
        from OCP.TopExp import TopExp_Explorer
        from OCP.TopoDS import TopoDS, TopoDS_Compound

        comp = TopoDS_Compound()
        builder = BRep_Builder()
        builder.MakeCompound(comp)
        for solid in solids:
            ex = TopExp_Explorer(solid, TopAbs_EDGE)
            while ex.More():
                builder.Add(comp, TopoDS.Edge(ex.Current()))
                ex.Next()
        return comp

    @staticmethod
    def _face_normal(face, u: float, v: float):
        from OCP.BRepAdaptor import BRepAdaptor_Surface
        from OCP.gp import gp_Pnt, gp_Vec

        try:
            surf = BRepAdaptor_Surface(face)
            p = gp_Pnt()
            du = gp_Vec()
            dv = gp_Vec()
            surf.D1(float(u), float(v), p, du, dv)
        except Exception:
            return None
        n = np.cross(
            np.array([du.X(), du.Y(), du.Z()]),
            np.array([dv.X(), dv.Y(), dv.Z()]))
        nn = float(np.linalg.norm(n))
        if nn <= 1e-300:
            return None
        return n / nn

    def _dist_to_edges(self, point: np.ndarray) -> float:
        from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
        from OCP.BRepExtrema import BRepExtrema_DistShapeShape
        from OCP.gp import gp_Pnt

        v = BRepBuilderAPI_MakeVertex(
            gp_Pnt(float(point[0]), float(point[1]), float(point[2]))).Vertex()
        d = BRepExtrema_DistShapeShape(v, self._edges)
        if not d.IsDone():
            d.Perform()
        if not d.IsDone():
            return 0.0  # fail safe: treat as degenerate
        return float(d.Value())

    def _cast_ray(self, point: np.ndarray,
                  direction: np.ndarray) -> Optional[list[int]]:
        """Crossing counts per solid, or None if the ray is degenerate."""
        from OCP.TopoDS import TopoDS
        from OCP.gp import gp_Ax1, gp_Dir, gp_Lin, gp_Pnt

        lin = gp_Lin(gp_Ax1(
            gp_Pnt(float(point[0]), float(point[1]), float(point[2])),
            gp_Dir(float(direction[0]), float(direction[1]),
                   float(direction[2]))))
        counts = []
        for inter in self._intersectors:
            inter.Perform(lin, 0.0, _RAY_PMAX)
            if not inter.IsDone():
                return None
            n = 0
            for i in range(1, inter.NbPnt() + 1):
                if float(inter.WParameter(i)) < self._tol:
                    return None  # hit at/behind the origin
                face = TopoDS.Face(inter.Face(i))
                nrm = self._face_normal(
                    face, float(inter.UParameter(i)),
                    float(inter.VParameter(i)))
                if nrm is None:
                    return None
                if abs(float(np.dot(nrm, direction))) < _TANGENT_COS:
                    return None  # grazing hit
                hp = inter.Pnt(i)
                if self._dist_to_edges(
                        np.array([hp.X(), hp.Y(), hp.Z()])) < self._tol:
                    return None  # within tol of an edge or vertex
                n += 1
            counts.append(n)
        return counts

    def _cast_bidirectional(self, point: np.ndarray,
                            direction: np.ndarray) -> Optional[list[int]]:
        """+d crossing counts per solid, or None if the pair is unusable.

        Casts both +d and -d.  Returns None when either side is degenerate
        (see _cast_ray) or when any solid's bidirectional crossing sum is
        odd, which proves the intersector missed or added a crossing on
        that line.
        """
        plus = self._cast_ray(point, direction)
        if plus is None:
            return None
        minus = self._cast_ray(point, -direction)
        if minus is None:
            return None
        for cp, cm in zip(plus, minus):
            if (cp + cm) % 2 == 1:
                return None
        return plus

    def classify(self, point: np.ndarray) -> str:
        """Return 'inside', 'outside', or 'unknown'.

        Seats up to 5 valid bidirectional ray pairs from the fixed
        direction list; all seated pairs must agree.  "unknown" is
        returned (never a guess) when fewer than 3 pairs validate or the
        pairs disagree with each other.
        """
        p = np.asarray(point, dtype=np.float64).reshape(3)
        verdicts = []
        for direction in _RAY_DIRECTIONS:
            counts = self._cast_bidirectional(p, direction)
            if counts is None:
                continue
            verdicts.append(
                "inside" if any(c % 2 == 1 for c in counts) else "outside")
            if len(verdicts) >= _N_RAYS_TARGET:
                break
        if len(verdicts) < _N_RAYS_MIN:
            return "unknown"
        if all(v == verdicts[0] for v in verdicts):
            return verdicts[0]
        return "unknown"


def _raise_classifier_disagreement(point: np.ndarray, occt_verdict: str,
                                   independent_verdict: str) -> None:
    """Refuse a decision the two classifiers do not agree on."""
    pt = (float(point[0]), float(point[1]), float(point[2]))
    exc = AssemblyError(
        f"point classifiers disagree at {pt}: "
        f"BRepClass3d_SolidClassifier={occt_verdict}, "
        f"multi-ray parity={independent_verdict}; refusing rather than "
        f"deciding on a single classifier",
        kind="ClassifierDisagreement")
    exc.occt_verdict = occt_verdict
    exc.independent_verdict = independent_verdict
    exc.point = pt
    raise exc


def _agreed_point_verdict(point: np.ndarray, model: BRepModel, tol: float,
                          ray: _MultiRayClassifier) -> str:
    """Classify a decision witness with both classifiers.

    Returns the agreed 'inside'/'outside' verdict, or the OCCT
    'boundary'/'unknown' state for the caller's existing refusal path.
    Any OCCT inside/outside verdict the multi-ray classifier does not
    confirm raises ClassifierDisagreement carrying both verdicts and the
    point.
    """
    occt = _classify_point_in_model(point, model, tol)
    if occt not in ("inside", "outside"):
        return occt
    independent = ray.classify(point)
    if independent == occt:
        return occt
    _raise_classifier_disagreement(point, occt, independent)


def classify_point_two_classifier(point, model, tol: float) -> dict:
    """Two-classifier point verdict (G5).

    `model` is a BRepModel or a sequence of TopoDS solids.  Returns a dict
    with keys: point, occt (BRepClass3d_SolidClassifier verdict),
    independent (multi-ray parity verdict), agreed (bool), decision (the
    agreed inside/outside verdict, or None when the classifiers do not
    agree).  The winding/tessellation arbiter in tools/review_probes stays
    the TESTING arbiter only; it is not used as a production classifier.
    """
    p = np.asarray(point, dtype=np.float64).reshape(-1)
    if p.shape != (3,):
        raise ValueError("point must have exactly 3 coordinates")
    if isinstance(model, BRepModel):
        if not model.solids:
            raise AssemblyError(
                "material patch classification requires closed OCCT "
                "solids; shell-only input is unsupported at this stage",
                kind="ShellOnlyClassificationUnsupported")
        solids = [sr.solid for sr in model.solids]
    else:
        solids = list(model)
    occt = _occt_point_verdict(p, solids, tol)
    independent = _MultiRayClassifier(solids, tol).classify(p)
    agreed = occt in ("inside", "outside") and independent == occt
    return {
        "point": (float(p[0]), float(p[1]), float(p[2])),
        "occt": occt,
        "independent": independent,
        "agreed": agreed,
        "decision": occt if agreed else None,
    }


def _point_boundary_distances(points: np.ndarray,
                              model: BRepModel) -> list[float]:
    """Distance from each point to the boundary of the model's solids.

    Measured against the model's FACES, not the solids: OCCT reports
    distance 0 for a vertex strictly inside a solid (containment), while
    the face distance is the true distance to the boundary.
    """
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
    from OCP.BRepExtrema import BRepExtrema_DistShapeShape
    from OCP.gp import gp_Pnt

    dists = []
    for p in points:
        v = BRepBuilderAPI_MakeVertex(
            gp_Pnt(float(p[0]), float(p[1]), float(p[2]))).Vertex()
        best = float("inf")
        for fr in model.faces:
            d = BRepExtrema_DistShapeShape(v, fr.face)
            if not d.IsDone():
                d.Perform()
            if d.IsDone():
                best = min(best, float(d.Value()))
        dists.append(best)
    return dists


def _witness_material_verdict(points: np.ndarray, classes: tuple[str, ...],
                              model: BRepModel, tol: float, *,
                              operand: str, parent_face_id: int,
                              piece_index: int,
                              min_points: int = 3) -> str:
    """Decide one patch's material state from dual-classified witnesses.

    G5 witness preference: witnesses at distance >= 10x tol from the other
    operand's boundary (BRepExtrema_DistShapeShape against the other model's
    faces) drive the decision, because the classifier confusion zone hugs
    the boundary.  Where at least `min_points` witnesses clear the band:
    - the far witnesses must be unanimous (else PatchClassificationInconsistent);
    - a far boundary/unknown verdict refuses (BoundaryOrUnknownPatch);
    - near-boundary boundary/unknown verdicts are confusion-zone noise and
      are ignored;
    - but a near witness that materially contradicts the far verdict (with
      both classifiers agreeing) still blocks with
      PatchClassificationInconsistent: that is the straddling-patch guard,
      and dropping it would let a thin sliver hide inside the far set.
    Where fewer than `min_points` witnesses clear the band, the preference
    is not feasible and the pre-G5 full-set unanimity rule applies
    unchanged.

    NOT applied in _solid_interior_points: those witnesses are deliberately
    near-boundary inward offsets, and thin shells would lose every witness
    under a 10x tol band.  NOT applied in _shell_records nesting either,
    which reuses those same near-boundary witnesses.
    """
    dists = _point_boundary_distances(points, model)
    band = 10.0 * float(tol)
    far = [c for c, d in zip(classes, dists) if d >= band]
    near = [c for c, d in zip(classes, dists) if d < band]

    def refuse(kind: str, detail: str) -> None:
        raise AssemblyError(
            f"{operand} face {parent_face_id} piece {piece_index}: "
            f"{detail}",
            kind=kind)

    def unanimous_among(which: list[str], what: str) -> str:
        invalid = sorted(set(
            x for x in which if x not in ("inside", "outside")))
        if invalid:
            refuse("BoundaryOrUnknownPatch",
                   f"{what} witness classification contains {invalid}")
        unique = set(which)
        if len(unique) != 1:
            refuse("PatchClassificationInconsistent",
                   f"supposedly split patch straddles material states "
                   f"{sorted(unique)} ({what} witnesses)")
        return which[0]

    if len(far) >= min_points:
        verdict = unanimous_among(far, "far-from-boundary")
        contra = sorted(set(
            x for x in near
            if x in ("inside", "outside") and x != verdict))
        if contra:
            refuse("PatchClassificationInconsistent",
                   f"near-boundary witness(es) contradict the "
                   f"far-witness verdict {verdict}: {contra}")
        return verdict
    return unanimous_among(list(classes), "witness")


def _decision_rule(operation: str, operand: str,
                   classification: str) -> tuple[bool, bool]:
    """Return (keep, reverse) for a patch in one of the four G2 states.

    The keep decision comes from the keep table (section 2.2) via
    keep_patch(); the table is the only place these rules live. The
    reverse flag is an orientation concern, not a keep rule: for
    difference, kept patches of B must be flipped to face outward from
    the A-B material.
    """
    keep = keep_patch(operation, operand, classification)
    reverse = (operation == "difference" and operand == "B"
               and classification == "IN" and keep)
    return keep, reverse


# ---------------------------------------------------------------------------
# G2.2 keep table: the core keep rules for regularized Booleans.
#
# A boundary patch of one operand is classified into one of four states
# relative to the other solid (G2.1):
#   IN      patch interior lies strictly inside the other solid;
#   OUT     patch interior lies strictly outside the other solid;
#   ON_SAME patch lies on the other solid's boundary and both outward
#           normals point the same way (materials on the same side);
#   ON_OPP  patch lies on the other solid's boundary and the outward
#           normals are opposite (materials touch from opposite sides).
# The ON states are decided from face-pair geometry, never from the 3D
# point classifier (which would just say ON).
#
# Why: for ON_SAME exactly one copy of the shared boundary survives in
# union and intersection (A's copy, for deterministic provenance). For
# ON_OPP the two materials meet from opposite sides, so the shared wall
# is interior to the union, zero-thickness for the intersection, and
# remains A's outer wall after subtracting B.
#
# This table is the ONLY place these rules live (pinned by
# tests/test_g2_keep_table.py, one test per cell).
# ---------------------------------------------------------------------------
_KEEP_TABLE = {
    "union": {
        "A": {"IN": False, "OUT": True, "ON_SAME": True, "ON_OPP": False},
        "B": {"IN": False, "OUT": True, "ON_SAME": False, "ON_OPP": False},
    },
    "intersection": {
        "A": {"IN": True, "OUT": False, "ON_SAME": True, "ON_OPP": False},
        "B": {"IN": True, "OUT": False, "ON_SAME": False, "ON_OPP": False},
    },
    "difference": {
        "A": {"IN": False, "OUT": True, "ON_SAME": False, "ON_OPP": True},
        "B": {"IN": True, "OUT": False, "ON_SAME": False, "ON_OPP": False},
    },
}


def keep_patch(operation: str, operand: str, state: str) -> bool:
    """Look up whether a patch in `state` is kept, from the keep table.

    Raises ValueError on any unknown operation, operand, or state rather
    than guessing a keep decision.
    """
    try:
        return _KEEP_TABLE[operation][operand][state]
    except KeyError as exc:
        raise ValueError(
            f"keep_patch: unknown operation/operand/state "
            f"({operation!r}, {operand!r}, {state!r})") from exc


def _reverse_face(face):
    from OCP.TopoDS import TopoDS
    return TopoDS.Face(face.Reversed())


def _classify_coincident_piece(piece, operand: str,
                               base_tol: float) -> str:
    """Return ON_SAME/ON_OPP for a coincident piece (G2.6).

    The state comes from the pair's sense relation. The witness must be
    ON the partner's support and the outward-normal dot sign must match
    the sense, else CoincidenceWitnessMismatch.
    """
    from OCP.BRep import BRep_Tool
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf
    from OCP.gp import gp_Pnt
    from OCP.TopoDS import TopoDS

    sense = piece.coincidence
    if sense not in ("same", "opp"):
        raise AssemblyError(
            f"{operand} face {piece.parent_face_id} piece "
            f"{piece.piece_index}: bad coincidence sense {sense!r}",
            kind="CoincidenceWitnessMismatch")
    partner_rec = piece.coincident_partner_rec
    if partner_rec is None:
        raise AssemblyError(
            f"{operand} face {piece.parent_face_id} piece "
            f"{piece.piece_index}: missing coincident partner",
            kind="CoincidenceWitnessMismatch")
    tol = max(float(base_tol),
              2.0 * float(BRep_Tool.Tolerance_s(piece.face)))
    points = _face_points(piece.face, tol)
    if len(points) == 0:
        raise AssemblyError(
            f"{operand} face {piece.parent_face_id} piece "
            f"{piece.piece_index}: no witnesses",
            kind="CoincidenceWitnessMismatch")
    partner_face = TopoDS.Face(partner_rec.face)
    psurf = BRepAdaptor_Surface(partner_face).Surface().Surface()
    for p in points:
        proj = GeomAPI_ProjectPointOnSurf(
            gp_Pnt(float(p[0]), float(p[1]), float(p[2])), psurf)
        if proj.NbPoints() == 0 or \
                float(proj.LowerDistance()) > tol:
            raise AssemblyError(
                f"{operand} face {piece.parent_face_id} piece "
                f"{piece.piece_index}: witness off partner support",
                kind="CoincidenceWitnessMismatch")
    # Normal-dot sign check at the first witness.
    from .coincidence import _outward_normal
    from OCP.BRepTools import BRepTools
    import numpy as np
    piece_face = TopoDS.Face(piece.face)
    try:
        u0, u1, v0, v1 = BRepTools.UVBounds_s(piece_face)
        uv = (float((u0 + u1) / 2.0), float((v0 + v1) / 2.0))
    except Exception:
        uv = None
    na = _outward_normal(
        piece_face,
        getattr(piece_face.Orientation(), "name", ""), uv) if uv else None
    if na is None:
        raise AssemblyError(
            f"{operand} face {piece.parent_face_id} piece "
            f"{piece.piece_index}: cannot compute outward normal",
            kind="CoincidenceWitnessMismatch")
    p0 = points[0]
    proj = GeomAPI_ProjectPointOnSurf(
        gp_Pnt(float(p0[0]), float(p0[1]), float(p0[2])), psurf)
    ub, vb = proj.LowerDistanceParameters()
    nb = _outward_normal(partner_face, partner_rec.orientation,
                         (float(ub), float(vb)))
    if nb is None:
        raise AssemblyError(
            f"{operand} face {piece.parent_face_id} piece "
            f"{piece.piece_index}: cannot compute partner normal",
            kind="CoincidenceWitnessMismatch")
    dot = float(np.dot(na, nb))
    if sense == "same" and dot < 0.9:
        raise AssemblyError(
            f"{operand} face {piece.parent_face_id} piece "
            f"{piece.piece_index}: normal-dot {dot:.3f} disagrees with "
            f"sense 'same'", kind="CoincidenceWitnessMismatch")
    if sense == "opp" and dot > -0.9:
        raise AssemblyError(
            f"{operand} face {piece.parent_face_id} piece "
            f"{piece.piece_index}: normal-dot {dot:.3f} disagrees with "
            f"sense 'opp'", kind="CoincidenceWitnessMismatch")
    return "ON_SAME" if sense == "same" else "ON_OPP"


# G12b: untouched-region classification.
#
# If a face has no section edges (status "unchanged"), it did not
# interact with the other model. Adjacent untouched faces connected
# through "clean" edges (edges not near any section vertex) share the
# same IN/OUT classification, so we classify one representative per
# region and propagate.


def _get_section_vertices(section_edges) -> list:
    """Extract vertex points from section edges."""
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopAbs import TopAbs_VERTEX
    from OCP.TopoDS import TopoDS
    from OCP.BRep import BRep_Tool

    verts = []
    for sec in section_edges:
        ex = TopExp_Explorer(sec.edge, TopAbs_VERTEX)
        while ex.More():
            v = TopoDS.Vertex(ex.Current())
            p = BRep_Tool.Pnt_s(v)
            verts.append((p.X(), p.Y(), p.Z()))
            ex.Next()
    return verts


def _edge_bbox_6tuple(edge):
    """Conservative bounding box of an edge as (xmin, ymin, zmin, xmax, ymax, zmax)."""
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib

    box = Bnd_Box()
    BRepBndLib.AddOptimal_s(edge, box, False, False)
    if box.IsVoid():
        return None
    lo, hi = box.CornerMin(), box.CornerMax()   # OCP 7.8 and 8.x
    return (lo.X(), lo.Y(), lo.Z(), hi.X(), hi.Y(), hi.Z())


def _bbox_near_point(bbox, pt, tol) -> bool:
    """Check if bbox is within tol of point."""
    xmin, ymin, zmin, xmax, ymax, zmax = bbox
    px, py, pz = pt
    # Clamp point to bbox, compute distance
    cx = max(xmin, min(px, xmax))
    cy = max(ymin, min(py, ymax))
    cz = max(zmin, min(pz, zmax))
    dx, dy, dz = px - cx, py - cy, pz - cz
    return (dx*dx + dy*dy + dz*dz) <= tol*tol


def _edge_is_clean(edge, section_vertices, tol) -> bool:
    """An edge is clean if its bbox is not within tol of any section vertex."""
    bbox = _edge_bbox_6tuple(edge)
    if bbox is None:
        return False  # Cannot determine, treat as not clean
    for v in section_vertices:
        if _bbox_near_point(bbox, v, tol):
            return False
    return True


def _build_untouched_regions(groups, section_edges, base_tol):
    """Build regions of untouched faces connected via clean edges.

    Returns (region_map, regions):
    - region_map: dict face_id -> region_id
    - regions: dict region_id -> list of face_ids
    """
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopAbs import TopAbs_EDGE
    from OCP.TopoDS import TopoDS

    # Untouched faces: status "unchanged" (no section edges)
    untouched = [fr for fr in groups if fr.status == "unchanged"]
    if not untouched:
        return {}, {}

    # Get section vertices
    section_vertices = _get_section_vertices(section_edges)

    # Get edges for each untouched face
    face_edges = {}
    for fr in untouched:
        # FaceSplitResult has pieces; for "unchanged", there's 1 piece
        # with the original face. We need the face object.
        # The piece has .face attribute.
        if not fr.pieces:
            continue
        face = fr.pieces[0].face
        edges = []
        ex = TopExp_Explorer(face, TopAbs_EDGE)
        while ex.More():
            edges.append(TopoDS.Edge(ex.Current()))
            ex.Next()
        face_edges[fr.parent_face_id] = edges

    # Union-find
    parent = {fr.parent_face_id: fr.parent_face_id for fr in untouched}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            # Deterministic: smaller face_id becomes parent
            if ra < rb:
                parent[rb] = ra
            else:
                parent[ra] = rb

    # Check all pairs of untouched faces for shared clean edges
    face_ids = list(face_edges.keys())
    for i in range(len(face_ids)):
        for j in range(i + 1, len(face_ids)):
            fid_a, fid_b = face_ids[i], face_ids[j]
            # Check if they share an edge
            shared_clean = False
            for ea in face_edges[fid_a]:
                for eb in face_edges[fid_b]:
                    if ea.IsSame(eb):
                        # Shared edge found, check if clean
                        if _edge_is_clean(ea, section_vertices, base_tol):
                            shared_clean = True
                            break
                if shared_clean:
                    break
            if shared_clean:
                union(fid_a, fid_b)

    # Build regions
    regions = {}
    region_map = {}
    for fid in face_ids:
        r = find(fid)
        region_map[fid] = r
        if r not in regions:
            regions[r] = []
        regions[r].append(fid)

    return region_map, regions


def _choose_representative(face_ids, groups):
    """Choose representative face for a region: largest area, tie-break by face_id.

    Deterministic: uses face_id for tie-breaking, not process-random ordering.
    """
    # Build map from face_id to FaceSplitResult
    fr_map = {fr.parent_face_id: fr for fr in groups}
    
    best = None
    best_area = -1.0
    for fid in sorted(face_ids):  # Sorted for determinism
        fr = fr_map.get(fid)
        if fr is None:
            continue
        # Use area_before as the area metric
        area = float(fr.area_before)
        if area > best_area or (area == best_area and (best is None or fid < best)):
            best = fid
            best_area = area
    return best


def _classify_pieces(model_a: BRepModel, model_b: BRepModel,
                     split: ModelSplitResult, operation: str,
                     base_tol: float) -> list[PatchDecision]:
    from OCP.BRep import BRep_Tool

    if split.unresolved_contacts:
        desc = ", ".join(
            f"A{a}/B{b}:{status}" for a, b, status
            in split.unresolved_contacts[:8])
        raise AssemblyError(
            "global assembly blocked by unresolved contact(s): " + desc,
            kind="UnresolvedContact")

    out: list[PatchDecision] = []

    def one_side(operand: str, groups, other: BRepModel):
        # One multi-ray classifier per side, built at the loosest piece
        # tolerance: a wider edge/near-origin discard band is the
        # conservative choice, and the per-piece OCCT tolerance still
        # governs the primary classifier.
        jobs = []
        for fr in groups:
            for piece in fr.pieces:
                tol = max(
                    float(base_tol),
                    2.0 * float(BRep_Tool.Tolerance_s(piece.face)))
    def one_side(operand: str, groups, other: BRepModel):
        # One multi-ray classifier per side, built at the loosest piece
        # tolerance: a wider edge/near-origin discard band is the
        # conservative choice, and the per-piece OCCT tolerance still
        # governs the primary classifier.
        #
        # G12b: Build untouched-face regions. Only region representatives
        # (and touched faces) are classified; results propagate to the
        # rest of the region.
        region_map, regions = _build_untouched_regions(
            groups, split.section_edges, base_tol)
        # Map face_id -> representative face_id (for non-representatives)
        # Representative maps to itself.
        rep_for = {}
        for region_id, face_ids in regions.items():
            rep = _choose_representative(face_ids, groups)
            for fid in face_ids:
                rep_for[fid] = rep
        # Faces to actually classify: touched faces + representatives
        # A face is classified if it's not in a region, or if it's the rep.
        def should_classify(face_id):
            if face_id not in rep_for:
                return True  # Touched face, or untouched but not in region
            return rep_for[face_id] == face_id

        jobs = []
        skipped = []  # (face_id, piece) for propagation
        for fr in groups:
            for piece in fr.pieces:
                tol = max(
                    float(base_tol),
                    2.0 * float(BRep_Tool.Tolerance_s(piece.face)))
                if should_classify(piece.parent_face_id):
                    jobs.append((piece, tol))
                else:
                    skipped.append((piece.parent_face_id, piece, tol))
        ray_tol = max([t for _, t in jobs], default=float(base_tol))
        ray = _MultiRayClassifier(
            [sr.solid for sr in other.solids], ray_tol)
        # Store decisions by (face_id, piece_index) for propagation
        decisions_by_key = {}
        for piece, tol in jobs:
            # G2.6: coincident pieces get their ON state from the pair
            # relation, not from 3D witnesses. The witness still has
            # to be ON the partner support with a matching normal-dot
            # sign, else CoincidenceWitnessMismatch.
            if piece.coincidence is not None:
                cls = _classify_coincident_piece(
                    piece, operand, base_tol)
                keep, rev = _decision_rule(operation, operand, cls)
                source = piece.face
                selected = _reverse_face(source) if keep and rev else (
                    source if keep else None)
                points = _face_points(source, tol)
                dec = PatchDecision(
                    operand=operand,
                    parent_face_id=piece.parent_face_id,
                    piece_index=piece.piece_index,
                    classification=cls,
                    keep=keep,
                    reverse_for_difference=rev,
                    witness_xyz=points[0],
                    witness_xyz_all=points,
                    witness_classifications=tuple([cls] * len(points)),
                    source_face=source,
                    selected_face=selected,
                )
                out.append(dec)
                decisions_by_key[(piece.parent_face_id, piece.piece_index)] = dec
                continue
            points = _face_points(piece.face, tol)
            classes = tuple(
                _agreed_point_verdict(p, other, tol, ray)
                for p in points)
            cls = _witness_material_verdict(
                points, classes, other, tol,
                operand=operand,
                parent_face_id=piece.parent_face_id,
                piece_index=piece.piece_index)
            # G2.1 canonical states: map the dual-classified
            # inside/outside verdict onto the four-state model before the
            # keep table.
            cls = {"inside": "IN", "outside": "OUT"}[cls]
            keep, rev = _decision_rule(operation, operand, cls)
            source = piece.face
            selected = _reverse_face(source) if keep and rev else (
                source if keep else None)
            dec = PatchDecision(
                operand=operand,
                parent_face_id=piece.parent_face_id,
                piece_index=piece.piece_index,
                classification=cls,
                keep=keep,
                reverse_for_difference=rev,
                witness_xyz=points[0],
                witness_xyz_all=points,
                witness_classifications=classes,
                source_face=source,
                selected_face=selected,
            )
            out.append(dec)
            decisions_by_key[(piece.parent_face_id, piece.piece_index)] = dec

        # G12b: Propagate representative decisions to skipped faces.
        for face_id, piece, tol in skipped:
            rep_id = rep_for[face_id]
            # Find the representative's decision (piece_index 0 for untouched)
            rep_key = (rep_id, 0)
            if rep_key not in decisions_by_key:
                # Fallback: should not happen for untouched faces
                continue
            rep_dec = decisions_by_key[rep_key]
            # Create propagated decision with this face's geometry
            # but the representative's classification.
            points = _face_points(piece.face, tol)
            keep, rev = _decision_rule(operation, operand, rep_dec.classification)
            source = piece.face
            selected = _reverse_face(source) if keep and rev else (
                source if keep else None)
            out.append(PatchDecision(
                operand=operand,
                parent_face_id=face_id,
                piece_index=piece.piece_index,
                classification=rep_dec.classification,
                keep=keep,
                reverse_for_difference=rev,
                witness_xyz=points[0],
                witness_xyz_all=points,
                witness_classifications=tuple(
                    [rep_dec.classification] * len(points)),
                source_face=source,
                selected_face=selected,
                propagated_from=rep_id,
            ))

        # Return region stats for the report
        return {
            "n_regions": len(regions),
            "n_propagated": len(skipped),
            "n_classified": len(jobs),
        }


    stats_a = one_side("A", split.faces_a, model_b)
    stats_b = one_side("B", split.faces_b, model_a)
    # Attach region stats to the output for the report
    # (stored on the function for access by caller)
    _classify_pieces.region_stats = {
        "A": stats_a,
        "B": stats_b,
    }
    return out


def _empty_compound():
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Compound
    c = TopoDS_Compound()
    b = BRep_Builder()
    b.MakeCompound(c)
    return c


# G13: per-Boolean-call volume cache.  Input volumes are measured several
# times during one Boolean call (operation invariants, evidence, shell
# checks); recomputing the adaptive 1e-10 integration each time is the
# dominant cost on rotated/titled geometry.  The cache is keyed by
# (id(shape), tol) and holds a strong reference to the shape so id() reuse
# after GC cannot alias a dead entry.  It is cleared at the start of every
# top-level Boolean call (see clear_volume_cache callers); entries never
# outlive one call, so stale results from mutated shapes are impossible.
_VOLUME_CACHE: dict = {}


def clear_volume_cache() -> None:
    """Drop all cached volume measurements.  Called once per Boolean call."""
    _VOLUME_CACHE.clear()


# Error budget for lower-precision internal volume checks (G13).  The
# coarse 1e-4 tolerance is used only for sign checks (vol > 0) on result
# shells and for operation-level volume invariants, where a relative
# error of 1e-4 cannot flip the verdict: volumes entering these checks
# are either exactly zero (degenerate, refused elsewhere) or bounded
# away from zero by construction tolerances >= 1e-7.  Any volume that is
# reported as evidence keeps the full 1e-10 adaptive integration.
_COARSE_VOLUME_TOL = 1e-4


def _shape_volume(shape, tol: float = 1e-10) -> float:
    """Adaptive volume measurement, including B-spline span integration.

    Results are cached for the duration of one Boolean call.  The default
    tol=1e-10 is the high-accuracy path used wherever a volume is reported
    as evidence.  Internal sign/coarse-bound checks may pass
    tol=_COARSE_VOLUME_TOL (1e-4); the error budget is documented above.
    """
    key = (id(shape), float(tol))
    hit = _VOLUME_CACHE.get(key)
    if hit is not None and hit[0] is shape:
        return hit[1]
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    g = GProp_GProps()
    err = BRepGProp.VolumePropertiesGK_s(
        shape, g, tol, True, True, False, False, False)
    if float(err) < 0.0:
        raise AssemblyError("adaptive volume integration failed",
                            kind="VolumeIntegrationFailed")
    vol = float(g.Mass())
    _VOLUME_CACHE[key] = (shape, vol)
    return vol


def _extract_shells(shape) -> list[object]:
    """Extract sewn shells and promote standalone closed periodic faces.

    OCCT sewing may legitimately return a FACE (or a compound of FACE
    objects) when a single periodic face already forms a closed skin, e.g. a
    full sphere.  That representation is geometrically closed but contains no
    explicit TopoDS_Shell.  We wrap only faces that are not already owned by a
    returned shell, and only when OCCT itself reports the one-face shell
    closed.
    """
    from OCP.BRep import BRep_Builder, BRep_Tool
    from OCP.TopAbs import TopAbs_FACE, TopAbs_SHELL
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS, TopoDS_Shell

    out = []
    claimed_faces = []

    ex = TopExp_Explorer(shape, TopAbs_SHELL)
    while ex.More():
        sh = TopoDS.Shell(ex.Current())
        if not any(sh.IsSame(x) for x in out):
            out.append(sh)
            ef = TopExp_Explorer(sh, TopAbs_FACE)
            while ef.More():
                f = TopoDS.Face(ef.Current())
                if not any(f.IsSame(x) for x in claimed_faces):
                    claimed_faces.append(f)
                ef.Next()
        ex.Next()

    builder = BRep_Builder()
    ef = TopExp_Explorer(shape, TopAbs_FACE)
    while ef.More():
        face = TopoDS.Face(ef.Current())
        if any(face.IsSame(x) for x in claimed_faces):
            ef.Next()
            continue
        sh = TopoDS_Shell()
        builder.MakeShell(sh)
        builder.Add(sh, face)
        if BRep_Tool.IsClosed_s(sh):
            sh.Closed(True)
            out.append(sh)
            claimed_faces.append(face)
        ef.Next()
    return out


def _make_outward_solid(shell):
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeSolid
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.BRepLib import BRepLib
    from OCP.TopAbs import TopAbs_SHELL
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS

    mk = BRepBuilderAPI_MakeSolid(shell)
    if not mk.IsDone():
        raise AssemblyError("could not convert closed shell to solid",
                            kind="SolidBuildFailed")
    solid = mk.Solid()
    if solid.IsNull() or not BRepLib.OrientClosedSolid_s(solid):
        raise AssemblyError("OCCT could not orient assembled closed shell",
                            kind="SolidOrientationFailed")
    if not BRepCheck_Analyzer(solid, True).IsValid():
        raise AssemblyError("assembled solid is not B-rep valid",
                            kind="SolidInvalid")
    ex = TopExp_Explorer(solid, TopAbs_SHELL)
    if not ex.More():
        raise AssemblyError("oriented solid lost its shell",
                            kind="SolidBuildFailed")
    outward_shell = TopoDS.Shell(ex.Current())
    return solid, outward_shell


def _solid_interior_points(solid, tol: float, *,
                           max_points: int = 7,
                           min_points: int = 3) -> np.ndarray:
    """Find several points just inside one closed shell.

    Shell nesting must not depend on a single center/grid witness: a point in
    the center of a large outer shell may also lie inside a nested cavity.  We
    therefore prefer several near-boundary inward offsets distributed across
    the shell and require containment decisions to agree for all of them.
    """
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    from OCP.TopAbs import TopAbs_FACE, TopAbs_IN
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.BRepClass import BRepClass_FaceClassifier
    from OCP.BRepTools import BRepTools
    from OCP.gp import gp_Pnt, gp_Pnt2d, gp_Vec

    if max_points < min_points or min_points < 1:
        raise ValueError("invalid solid witness count")

    clf = BRepClass3d_SolidClassifier(solid)
    # G5: witness search must not trust the OCCT classifier alone either;
    # a false IN here (finding F4) would plant a corrupt nesting witness.
    ray = _MultiRayClassifier([solid], tol)

    def is_in(x: np.ndarray) -> bool:
        clf.Perform(gp_Pnt(float(x[0]), float(x[1]), float(x[2])),
                    float(tol))
        if clf.State() != TopAbs_IN:
            return False
        return ray.classify(x) == "inside"

    lo, hi = _bbox(solid)
    scale = max(float(np.linalg.norm(hi - lo)), 1.0)
    sep = max(16.0 * float(tol), 1e-9 * scale)
    points: list[np.ndarray] = []

    def add(q: np.ndarray):
        if all(float(np.linalg.norm(q - p)) > sep for p in points):
            points.append(q)

    frac = (0.5, 0.25, 0.75, 0.125, 0.875,
            0.375, 0.625, 0.0625, 0.9375)

    # Prefer points just inside different boundary locations.  Test both
    # normal directions rather than trusting imported face orientation.
    ex = TopExp_Explorer(solid, TopAbs_FACE)
    while ex.More() and len(points) < max_points:
        face = TopoDS.Face(ex.Current())
        try:
            u0, u1, v0, v1 = map(float, BRepTools.UVBounds_s(face))
            if not all(np.isfinite([u0, u1, v0, v1])):
                ex.Next()
                continue
            surf = BRepAdaptor_Surface(face)
            for fu in frac:
                if len(points) >= max_points:
                    break
                u = u0 + (u1 - u0) * fu
                for fv in frac:
                    v = v0 + (v1 - v0) * fv
                    fc = BRepClass_FaceClassifier(
                        face, gp_Pnt2d(float(u), float(v)),
                        float(tol), True)
                    if fc.State() != TopAbs_IN:
                        continue
                    p = gp_Pnt()
                    du, dv = gp_Vec(), gp_Vec()
                    surf.D1(float(u), float(v), p, du, dv)
                    n = np.cross(
                        np.array([du.X(), du.Y(), du.Z()]),
                        np.array([dv.X(), dv.Y(), dv.Z()]))
                    nn = float(np.linalg.norm(n))
                    if nn <= 1e-300:
                        continue
                    n /= nn
                    x = _p3(p)
                    accepted = False
                    for eps in (1e-7, 1e-6, 1e-5, 1e-4):
                        d = max(16.0 * tol, eps * scale)
                        for sign in (-1.0, 1.0):
                            q = x + sign * d * n
                            if is_in(q):
                                add(q)
                                accepted = True
                                break
                        if accepted:
                            break
                    if len(points) >= max_points:
                        break
        except Exception:
            pass
        ex.Next()

    # Fallback points are allowed to fill out the witness set, but never
    # replace the requirement for multiple consistent witnesses.
    #
    # G10 Problem B: the center-of-mass witness is optional. The adaptive
    # volume integration behind it is expensive (it dominates NURBS
    # witness time), so it runs only when the ordinary face-based
    # strategy failed to produce enough points, at reduced accuracy, and
    # a failed integration only skips this candidate: it must never turn
    # an otherwise certifiable operation into a refusal.
    if len(points) < min_points:
        props = GProp_GProps()
        err = BRepGProp.VolumePropertiesGK_s(
            solid, props, 1e-4, True, True, True, False, False)
        if float(err) >= 0.0:
            cm = _p3(props.CentreOfMass())
            if is_in(cm):
                add(cm)

    if len(points) < min_points:
        for fx in frac:
            x = lo[0] + (hi[0] - lo[0]) * fx
            for fy in frac:
                y = lo[1] + (hi[1] - lo[1]) * fy
                for fz in frac:
                    z = lo[2] + (hi[2] - lo[2]) * fz
                    q = np.array([x, y, z], dtype=np.float64)
                    if is_in(q):
                        add(q)
                    if len(points) >= max_points:
                        break
                if len(points) >= max_points:
                    break
            if len(points) >= max_points:
                break

    if len(points) < min_points:
        raise AssemblyError(
            f"only {len(points)} stable shell interior witness(es) found; "
            f"{min_points} required",
            kind="InsufficientShellWitnesses")
    return np.vstack(points[:max_points])

def _shell_records(shells: list[object], tol: float
                   ) -> list[ShellAssemblyRecord]:
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.TopAbs import TopAbs_IN
    from OCP.gp import gp_Pnt

    tmp = []
    for i, sh in enumerate(shells):
        solid, outward = _make_outward_solid(sh)
        # G13: sign check only; coarse precision is safe per the documented
        # error budget above.
        vol = abs(_shape_volume(solid, tol=_COARSE_VOLUME_TOL))
        if not vol > 0:
            raise AssemblyError("assembled shell has non-positive volume",
                                kind="ZeroVolumeShell")
        points = _solid_interior_points(solid, tol)
        tmp.append({
            "index": i, "shell": sh, "outward": outward,
            "solid": solid, "volume": vol, "points": points,
            "point": points[0], "containers": []
        })

    from OCP.TopAbs import TopAbs_IN, TopAbs_ON, TopAbs_OUT
    # G5: one independent multi-ray classifier per candidate parent solid;
    # containment of a child shell in a parent is a material decision, so
    # the OCCT verdict on every witness must be confirmed by ray parity.
    ray_by_parent = {p["index"]: _MultiRayClassifier([p["solid"]], tol)
                     for p in tmp}
    for child in tmp:
        for parent in tmp:
            if parent is child:
                continue
            classifier = BRepClass3d_SolidClassifier(parent["solid"])
            ray = ray_by_parent[parent["index"]]
            states = []
            for q in child["points"]:
                classifier.Perform(
                    gp_Pnt(float(q[0]), float(q[1]), float(q[2])),
                    float(tol))
                st = classifier.State()
                if st == TopAbs_IN:
                    occt_name = "inside"
                elif st == TopAbs_OUT:
                    occt_name = "outside"
                else:
                    occt_name = None
                if occt_name is not None:
                    independent = ray.classify(q)
                    if independent != occt_name:
                        _raise_classifier_disagreement(
                            q, occt_name, independent)
                states.append(st)

            if any(st == TopAbs_ON for st in states):
                raise AssemblyError(
                    f"shell {child['index']} has a nesting witness on "
                    f"candidate parent {parent['index']} boundary",
                    kind="ShellContainmentAmbiguous")
            if all(st == TopAbs_IN for st in states):
                child["containers"].append(parent["index"])
            elif all(st == TopAbs_OUT for st in states):
                pass
            else:
                raise AssemblyError(
                    f"shell {child['index']} has inconsistent containment "
                    f"against shell {parent['index']}",
                    kind="ShellContainmentAmbiguous")

    parents: dict[int, Optional[int]] = {}
    for child in tmp:
        cs = child["containers"]
        if not cs:
            parents[child["index"]] = None
        else:
            parents[child["index"]] = min(
                cs, key=lambda j: tmp[j]["volume"])

    def depth(i: int) -> int:
        seen = set()
        d = 0
        p = parents[i]
        while p is not None:
            if p in seen:
                raise AssemblyError("cycle in shell containment hierarchy",
                                    kind="ShellNestingCycle")
            seen.add(p)
            d += 1
            p = parents[p]
        return d

    return [
        ShellAssemblyRecord(
            shell_index=x["index"],
            shell=x["shell"],
            normalized_shell=x["outward"],
            temp_solid=x["solid"],
            volume=x["volume"],
            interior_point=x["point"],
            parent_shell=parents[x["index"]],
            depth=depth(x["index"]),
        )
        for x in tmp
    ]


def _build_nested_solids(records: list[ShellAssemblyRecord]
                         ) -> list[SolidAssemblyRecord]:
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeSolid
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.BRepLib import BRepLib
    from OCP.TopoDS import TopoDS

    by_id = {r.shell_index: r for r in records}
    out = []
    for r in records:
        if r.depth % 2:
            continue
        mk = BRepBuilderAPI_MakeSolid()
        mk.Add(r.normalized_shell)
        cavities = sorted(
            x.shell_index for x in records
            if x.parent_shell == r.shell_index and x.depth == r.depth + 1)
        for cid in cavities:
            cavity = by_id[cid].normalized_shell
            mk.Add(TopoDS.Shell(cavity.Reversed()))
        if not mk.IsDone():
            raise AssemblyError("nested shell set could not build a solid",
                                kind="SolidBuildFailed")
        solid = mk.Solid()
        if solid.IsNull() or not BRepLib.OrientClosedSolid_s(solid):
            raise AssemblyError("nested result solid is not orientable",
                                kind="SolidOrientationFailed")
        if not BRepCheck_Analyzer(solid, True).IsValid():
            raise AssemblyError("nested result solid is B-rep invalid",
                                kind="SolidInvalid")
        # G13: sign check only; coarse precision is safe per the documented
        # error budget above.
        vol = _shape_volume(solid, tol=_COARSE_VOLUME_TOL)
        if not vol > 0:
            raise AssemblyError("result solid has non-positive volume",
                                kind="SolidInvalid")
        out.append(SolidAssemblyRecord(
            solid_index=len(out), solid=solid,
            outer_shell=r.shell_index,
            cavity_shells=tuple(cavities), volume=vol))
    return out


def _unique_edges(shape) -> list[object]:
    from OCP.TopAbs import TopAbs_EDGE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS

    out = []
    ex = TopExp_Explorer(shape, TopAbs_EDGE)
    while ex.More():
        e = TopoDS.Edge(ex.Current())
        if not any(e.IsSame(x) for x in out):
            out.append(e)
        ex.Next()
    return out


def _shape_has_edge(shape, edge) -> bool:
    from OCP.TopAbs import TopAbs_EDGE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS

    ex = TopExp_Explorer(shape, TopAbs_EDGE)
    while ex.More():
        if TopoDS.Edge(ex.Current()).IsSame(edge):
            return True
        ex.Next()
    return False


def _edge_length(edge) -> float:
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps

    g = GProp_GProps()
    BRepGProp.LinearProperties_s(edge, g, False, False)
    return float(g.Mass())


def _point_edge_distance(point: np.ndarray, edge) -> float:
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
    from OCP.BRepExtrema import BRepExtrema_DistShapeShape
    from OCP.gp import gp_Pnt

    v = BRepBuilderAPI_MakeVertex(
        gp_Pnt(float(point[0]), float(point[1]), float(point[2]))).Vertex()
    d = BRepExtrema_DistShapeShape(v, edge)
    if not d.IsDone():
        d.Perform()
    if not d.IsDone():
        return float("inf")
    return float(d.Value())


def _edge_on_face_boundary(edge, face, base_tol: float) -> bool:
    """Geometric check: does the edge lie on the face's boundary wire.

    Used as a fallback when the splitter rebuilt edges (TShape identity
    lost). A result edge that is a sub-edge of an original boundary
    edge is a source_boundary, not a coincident_boundary.
    """
    from OCP.TopAbs import TopAbs_EDGE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS
    try:
        ex = TopExp_Explorer(TopoDS.Face(face), TopAbs_EDGE)
        while ex.More():
            bnd = TopoDS.Edge(ex.Current())
            if _edge_matches_ref_edge(edge, bnd, float(base_tol),
                                      float(base_tol), float(base_tol)):
                return True
            ex.Next()
    except Exception:
        return False
    return False


def _edge_matches_ref_edge(edge, ref_edge, verify_tolerance: float,
                           edge_tolerance: float,
                           base_tol: float) -> bool:
    """Recognize a result edge as lying on a reference edge.

    Splitter/sewing may preserve the edge, copy it, or keep only a
    sub-edge. Therefore TShape identity is tried first, then a
    conservative geometric sub-edge check: the result edge may not
    exceed the reference length and several points along it must lie on
    the reference curve.
    """
    if edge.IsSame(ref_edge):
        return True

    from OCP.BRepAdaptor import BRepAdaptor_Curve

    le = _edge_length(edge)
    ls = _edge_length(ref_edge)
    tol = max(float(base_tol), float(verify_tolerance),
              float(edge_tolerance))
    len_tol = max(16.0 * tol, 1e-8 * max(le, ls, 1.0))
    if le > ls + len_tol:
        return False

    ce = BRepAdaptor_Curve(edge)
    t0 = float(ce.FirstParameter())
    t1 = float(ce.LastParameter())
    if not (np.isfinite(t0) and np.isfinite(t1) and t1 >= t0):
        return False
    match_tol = max(8.0 * tol, 1e-9 * max(le, ls, 1.0))
    for q in (0.0, 0.125, 0.25, 0.5, 0.75, 0.875, 1.0):
        t = t0 + (t1 - t0) * q
        p = ce.Value(float(t))
        x = np.array([p.X(), p.Y(), p.Z()], dtype=np.float64)
        if _point_edge_distance(x, ref_edge) > match_tol:
            return False
    return True


def _edge_matches_section(edge, sec, base_tol: float) -> bool:
    """Recognize a result edge as lying on a verified section edge."""
    return _edge_matches_ref_edge(
        edge, sec.edge, float(sec.verify_tolerance),
        float(sec.edge_tolerance), float(base_tol))


# ---------------------------------------------------------------------------
# G12a: call-scoped assembly indexes.
#
# _build_edge_lineage used to re-explore every face's edges and recompute
# every reference edge's length for each (result edge, candidate) pair.
# The classes below hoist that repeated work into per-call indexes built
# once at the top of _build_edge_lineage. Every predicate below is
# pairwise-equal to the linear scan it replaces:
#   - _EdgeSet.__contains__ uses the same TShape-identity (IsSame)
#     predicate as _shape_has_edge over the same edge multiset.
#   - _AssemblyIndexes._match_ref_edge runs the original
#     _edge_matches_ref_edge for every pair the conservative box
#     pre-filter cannot rule out. The pre-filter only skips pairs whose
#     bounding boxes are farther apart than the match tolerance, in
#     which case no sample point of the edge could be within the
#     tolerance of the reference edge and the original would return
#     False. Skips are therefore provably verdict-neutral.
# All state is created fresh per call. There is no module-level cache,
# no mutable default argument, and nothing is shared between calls.
# ---------------------------------------------------------------------------


def _edge_bbox(edge) -> tuple[np.ndarray, np.ndarray]:
    """Axis-aligned bounding box of an edge as (lo, hi)."""
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib

    box = Bnd_Box()
    BRepBndLib.Add_s(edge, box)
    if box.IsVoid():
        inf = float("inf")
        return (np.full(3, -inf), np.full(3, inf))
    _lo, _hi = box.CornerMin(), box.CornerMax()   # OCP 7.8 and 8.x
    lo = np.array([_lo.X(), _lo.Y(), _lo.Z()], dtype=np.float64)
    hi = np.array([_hi.X(), _hi.Y(), _hi.Z()], dtype=np.float64)
    return (lo, hi)


def _bbox_separated(bb1, bb2, margin: float) -> bool:
    """True only if every point of box 1 is farther than margin from
    every point of box 2.

    If the boxes are separated by more than margin on any axis, the
    Euclidean distance between any two points (one per box) exceeds
    margin on that axis alone. An edge is contained in its box, so a
    True result proves the geometric edge matcher (which needs sample
    points within margin of the reference edge) would return False.
    """
    (lo1, hi1), (lo2, hi2) = bb1, bb2
    m = float(margin)
    return bool(np.any(hi1 + m < lo2) or np.any(hi2 + m < lo1))


class _EdgeSet:
    """Call-scoped unique-edge set for one shape.

    Membership uses TShape identity (IsSame), exactly like
    _shape_has_edge; the explorer scan happens once here instead of
    once per query.
    """

    __slots__ = ("_edges",)

    def __init__(self, shape):
        self._edges = _unique_edges(shape)

    def __contains__(self, edge) -> bool:
        return any(e.IsSame(edge) for e in self._edges)

    def __len__(self) -> int:
        return len(self._edges)


class _RefEdgeData:
    """Call-scoped cached data for one reference edge (section edge,
    coincident-boundary tool edge, or parent-face boundary edge)."""

    __slots__ = ("edge", "bbox", "length", "verify_tolerance",
                 "edge_tolerance")

    def __init__(self, edge, verify_tolerance: float,
                 edge_tolerance: float):
        self.edge = edge
        self.bbox = _edge_bbox(edge)
        self.length = _edge_length(edge)
        self.verify_tolerance = float(verify_tolerance)
        self.edge_tolerance = float(edge_tolerance)


class _AssemblyIndexes:
    """Call-scoped indexes for one _build_edge_lineage call.

    Built once per call from (result_shape, selected, sections, tools,
    models). Sections/tools are duck-typed: section-likes expose
    .edge/.verify_tolerance/.edge_tolerance/.face_a/.face_b/.edge_index,
    tool-likes are dicts with "edge"/"operand"/"face_id"/"edge_id",
    models maps operand -> BRepModel.
    """

    def __init__(self, result_shape, selected, sections, tools,
                 models, base_tol: float):
        self.base_tol = float(base_tol)
        self.result_edges = _unique_edges(result_shape)
        self.result_bboxes = [_edge_bbox(e) for e in self.result_edges]
        self.result_lengths = [_edge_length(e) for e in self.result_edges]
        self.n_result_edges = len(self.result_edges)
        self.selected = list(selected)
        self.sewed_sets = [self._sewed_edge_set(d) for d in self.selected]
        self.sections = [_RefEdgeData(s.edge, s.verify_tolerance,
                                      s.edge_tolerance) for s in sections]
        self.section_keys = [(s.face_a, s.face_b, s.edge_index)
                             for s in sections]
        self.n_sections = len(self.sections)
        self.tools = [_RefEdgeData(t["edge"], self.base_tol,
                                   self.base_tol) for t in tools]
        self.tool_keys = [(t["operand"], t["face_id"], t["edge_id"])
                          for t in tools]
        self.n_tools = len(self.tools)
        self.face_map = {op: {f.face_id: f.face for f in m.faces}
                         for op, m in models.items()}
        self._parent_sets: dict = {}
        self._parent_boundary: dict = {}

    @staticmethod
    def _sewed_edge_set(d):
        sf = d.sewed_face if d.sewed_face is not None else d.selected_face
        return _EdgeSet(sf) if sf is not None else None

    def register_parent_face(self, key, face) -> None:
        """Register a parent face under an arbitrary key (used by the
        lazy parent lookup and directly by tests)."""
        self._parent_sets[key] = _EdgeSet(face)
        self._parent_boundary[key] = [_RefEdgeData(b, self.base_tol,
                                                   self.base_tol)
                                      for b in _unique_edges(face)]

    def _parent_key(self, operand, face_id):
        return (operand, face_id)

    def parent_set(self, operand, face_id):
        key = self._parent_key(operand, face_id)
        if key not in self._parent_sets:
            face = self.face_map[operand].get(face_id)
            if face is None:
                return None
            self.register_parent_face(key, face)
        return self._parent_sets[key]

    def sewed_contains(self, di: int, edge) -> bool:
        s = self.sewed_sets[di]
        return s is not None and edge in s

    def parent_contains(self, operand, face_id, edge) -> bool:
        s = self.parent_set(operand, face_id)
        return s is not None and edge in s

    def _match_ref_edge(self, edge_idx: int, ref: _RefEdgeData,
                        base_tol: float) -> bool:
        """Indexed form of _edge_matches_ref_edge.

        Pairwise-equal to the original: the IsSame fast path and the
        length gate run first on cached values; the box pre-filter only
        skips pairs the original provably answers False for; otherwise
        the original function decides.
        """
        edge = self.result_edges[edge_idx]
        if edge.IsSame(ref.edge):
            return True
        le = self.result_lengths[edge_idx]
        ls = ref.length
        tol = max(float(base_tol), ref.verify_tolerance,
                  ref.edge_tolerance)
        len_tol = max(16.0 * tol, 1e-8 * max(le, ls, 1.0))
        if le > ls + len_tol:
            return False
        match_tol = max(8.0 * tol, 1e-9 * max(le, ls, 1.0))
        if _bbox_separated(self.result_bboxes[edge_idx], ref.bbox,
                           match_tol):
            return False
        return _edge_matches_ref_edge(edge, ref.edge,
                                      ref.verify_tolerance,
                                      ref.edge_tolerance,
                                      float(base_tol))

    def edge_matches_section(self, edge_idx: int, sec_idx: int,
                             base_tol: float) -> bool:
        """Indexed form of _edge_matches_section (pairwise-equal)."""
        return self._match_ref_edge(edge_idx, self.sections[sec_idx],
                                    base_tol)

    def edge_matches_tool(self, edge_idx: int, tool_idx: int,
                          base_tol: float) -> bool:
        """Indexed form of the coincident-tool _edge_matches_ref_edge
        call (pairwise-equal)."""
        return self._match_ref_edge(edge_idx, self.tools[tool_idx],
                                    base_tol)

    def edge_on_face_boundary(self, edge_idx: int, key,
                              base_tol: float) -> bool:
        """Indexed form of _edge_on_face_boundary (pairwise-equal)."""
        if key not in self._parent_boundary:
            return False
        for ref in self._parent_boundary[key]:
            if self._match_ref_edge(edge_idx, ref, base_tol):
                return True
        return False


def _build_edge_lineage(result_shape, selected: list[PatchDecision],
                        split: ModelSplitResult,
                        model_a: BRepModel, model_b: BRepModel,
                        base_tol: float) -> list[EdgeLineageRecord]:
    """Attach final result edges to selected patches and section evidence.

    G12a: all repeated scans go through the call-scoped _AssemblyIndexes
    built below (one per call). The lineage records produced are
    identical to the old linear scans by construction: every predicate
    is pairwise-equal to the scan it replaces.
    """
    btol = float(base_tol)
    idx = _AssemblyIndexes(
        result_shape, selected,
        list(split.section_edges), list(split.coincident_boundary_tools),
        {"A": model_a, "B": model_b}, btol)
    out = []
    for i in range(idx.n_result_edges):
        edge = idx.result_edges[i]
        refs = []
        parents = []
        operands = []
        source_boundary = []
        for di, d in enumerate(selected):
            if idx.sewed_contains(di, edge):
                ref = (d.operand, d.parent_face_id, d.piece_index)
                if ref not in refs:
                    refs.append(ref)
                pf = (d.operand, d.parent_face_id)
                if pf not in parents:
                    parents.append(pf)
                if d.operand not in operands:
                    operands.append(d.operand)
                if (idx.parent_contains(d.operand, d.parent_face_id, edge)
                        or idx.edge_on_face_boundary(
                            i, (d.operand, d.parent_face_id), btol)):
                    if pf not in source_boundary:
                        source_boundary.append(pf)

        intersections = []
        for sj in range(idx.n_sections):
            if idx.edge_matches_section(i, sj, btol):
                key = idx.section_keys[sj]
                if key not in intersections:
                    intersections.append(key)

        # G2.6: coincident_boundary. A result edge that matches a
        # partner-face boundary edge used as an overlap-split tool, and
        # that lies in the interior of at least one parent face (i.e.
        # it was created by the overlap cut, not merely coincident with
        # an original boundary), came from the other operand's boundary
        # via the overlap split.
        coincident_sources = []
        interior_to_a_parent = False
        for tj in range(idx.n_tools):
            if idx.edge_matches_tool(i, tj, btol):
                key = idx.tool_keys[tj]
                if key not in coincident_sources:
                    coincident_sources.append(key)
        if coincident_sources:
            for operand, fid in parents:
                if idx.parent_set(operand, fid) is not None and not (
                        idx.edge_on_face_boundary(i, (operand, fid), btol)):
                    interior_to_a_parent = True
                    break
            source_boundary.extend(
                s for s in coincident_sources if s not in source_boundary)

        if intersections:
            kind = "boolean_section"
        elif coincident_sources and interior_to_a_parent:
            kind = "coincident_boundary"
        elif source_boundary:
            kind = "source_boundary"
        elif refs:
            kind = "split_or_sewn_boundary"
        else:
            kind = "unattributed"

        out.append(EdgeLineageRecord(
            result_edge_index=i,
            piece_refs=tuple(sorted(refs)),
            parent_faces=tuple(sorted(parents)),
            operands=tuple(sorted(operands)),
            intersection_refs=tuple(sorted(intersections)),
            source_boundary_refs=tuple(sorted(source_boundary)),
            provenance_kind=kind,
            verified_pcurves=bool(intersections),
        ))
    return out


def _build_section_payloads(split: ModelSplitResult,
                            edge_lineage: list[EdgeLineageRecord]
                            ) -> list[SectionPayloadRecord]:
    """Persist verified section evidence and link it to final result edges."""
    result_edges: dict[tuple[int, int, int], list[int]] = {}
    for lin in edge_lineage:
        for ref in lin.intersection_refs:
            result_edges.setdefault(ref, []).append(lin.result_edge_index)

    out = []
    seen = set()
    for sec in split.section_edges:
        key = (sec.face_a, sec.face_b, sec.edge_index)
        if key in seen:
            continue
        seen.add(key)
        out.append(SectionPayloadRecord(
            face_a=sec.face_a,
            face_b=sec.face_b,
            section_edge_index=sec.edge_index,
            parameters=np.array(sec.parameters, dtype=np.float64, copy=True),
            xyz=np.array(sec.xyz, dtype=np.float64, copy=True),
            uv_a=np.array(sec.uv_a, dtype=np.float64, copy=True),
            uv_b=np.array(sec.uv_b, dtype=np.float64, copy=True),
            edge_tolerance=float(sec.edge_tolerance),
            verify_tolerance=float(sec.verify_tolerance),
            max_surface_error_a=float(sec.max_surface_error_a),
            max_surface_error_b=float(sec.max_surface_error_b),
            max_cross_surface_error=float(sec.max_cross_surface_error),
            min_transversality=float(sec.min_transversality),
            max_transversality=float(sec.max_transversality),
            risk_flags=tuple(sec.risk_flags),
            repaired_same_parameter=bool(sec.repaired_same_parameter),
            exact_curve_on_surface_checked=bool(
                sec.exact_curve_on_surface_checked),
            exact_surface_error_a=(
                None if sec.exact_surface_error_a is None
                else float(sec.exact_surface_error_a)),
            exact_surface_error_b=(
                None if sec.exact_surface_error_b is None
                else float(sec.exact_surface_error_b)),
            shadow_crosschecked=bool(sec.shadow_crosschecked),
            shadow_max_distance=(None if sec.shadow_max_distance is None
                                 else float(sec.shadow_max_distance)),
            shadow_length_rel_error=(
                None if sec.shadow_length_rel_error is None
                else float(sec.shadow_length_rel_error)),
            result_edge_indices=tuple(sorted(set(result_edges.get(key, [])))),
        ))
    return out


def _compound_solids(solids: list[SolidAssemblyRecord]):
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Compound
    c = TopoDS_Compound()
    b = BRep_Builder()
    b.MakeCompound(c)
    for s in solids:
        b.Add(c, s.solid)
    return c


def _solids_touch(solids, tol: float) -> bool:
    """True when any two result solids touch geometrically (G2.5).

    Touching (distance ~0) outer solids make a union non-manifold;
    disjoint solids are a legitimate multi-solid result.
    """
    from OCP.BRepExtrema import BRepExtrema_DistShapeShape
    shapes = [s.solid for s in solids]
    for i in range(len(shapes)):
        for j in range(i + 1, len(shapes)):
            try:
                d = BRepExtrema_DistShapeShape(shapes[i], shapes[j])
                if d.Value() <= tol:
                    return True
            except Exception:
                # An unreadable distance is not evidence of touching.
                continue
    return False


def assemble_boolean(model_a: BRepModel, model_b: BRepModel,
                     split: ModelSplitResult, operation: str, *,
                     base_tol: float = 1e-7,
                     sew_tol: Optional[float] = None,
                     allow_nonmanifold: bool = False
                     ) -> BooleanAssemblyResult:
    """Classify exact B-rep patches and assemble union/intersection/A-B.

    allow_nonmanifold=True returns a compound of touching solids for a
    touching-only union instead of refusing NonManifoldResult.
    """
    # G13: one Boolean call, one volume-cache lifetime.
    clear_volume_cache()
    from OCP.BRep import BRep_Tool
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Sewing
    from OCP.BRepCheck import BRepCheck_Analyzer

    operation = str(operation).lower()
    if operation not in ("union", "intersection", "difference"):
        raise ValueError("operation must be union/intersection/difference")
    if not base_tol > 0:
        raise ValueError("base_tol must be positive")

    decisions = _classify_pieces(
        model_a, model_b, split, operation, float(base_tol))
    selected = [d for d in decisions if d.keep]
    # G12b: capture region stats from _classify_pieces
    region_stats = getattr(_classify_pieces, "region_stats", {})

    if not selected:
        empty = _empty_compound()
        return BooleanAssemblyResult(
            operation=operation, decisions=decisions, selected_faces=0,
            sewed_shape=empty, shells=[], solids=[], shape=empty,
            volume=0.0, free_edges=0, multiple_edges=0,
            edge_lineage=[],
            section_payloads=[],
            notes=["empty material result"],
            region_stats=region_stats)

    if sew_tol is None:
        face_tols = [
            float(BRep_Tool.Tolerance_s(d.selected_face))
            for d in selected if d.selected_face is not None]
        sew_tol = max(float(base_tol),
                      2.0 * max(face_tols, default=base_tol))
    if not sew_tol > 0:
        raise ValueError("sew_tol must be positive")

    sew = BRepBuilderAPI_Sewing(
        float(sew_tol), True, True, True, False)
    sew.SetSameParameterMode(True)
    sew.SetLocalTolerancesMode(True)
    for d in selected:
        sew.Add(d.selected_face)
    sew.Perform()

    free = int(sew.NbFreeEdges())
    multi = int(sew.NbMultipleEdges())
    if multi:
        raise AssemblyError(
            f"sewing produced {multi} non-manifold multiple edge(s)",
            kind="NonManifoldAssembly")
    if free:
        raise AssemblyError(
            f"sewing left {free} free boundary edge(s)",
            kind="OpenAssembly")

    sewed = sew.SewedShape()
    if sewed.IsNull():
        raise AssemblyError("OCCT sewing returned a null result",
                            kind="SewingFailed")
    if not BRepCheck_Analyzer(sewed, True).IsValid():
        raise AssemblyError("sewed patch complex is B-rep invalid",
                            kind="SewingInvalid")

    for d in selected:
        try:
            d.sewed_face = (sew.Modified(d.selected_face)
                            if sew.IsModified(d.selected_face)
                            else d.selected_face)
        except Exception:
            d.sewed_face = d.selected_face

    raw_shells = _extract_shells(sewed)
    if not raw_shells:
        raise AssemblyError("closed sewing result contains no shells",
                            kind="SewingFailed")

    shell_records = _shell_records(raw_shells, float(sew_tol))
    solids = _build_nested_solids(shell_records)
    if not solids:
        raise AssemblyError("closed shells produced no material solids",
                            kind="SolidBuildFailed")

    result_shape = (solids[0].solid if len(solids) == 1
                    else _compound_solids(solids))
    # G2.5: a union of solids that only touch (along an edge or at a
    # point) is non-manifold. Refuse typed by default; the optional
    # allow_nonmanifold flag returns the compound with the contact
    # recorded in the report notes.
    touching_solids = (operation == "union" and len(solids) > 1
                       and _solids_touch(solids, max(float(base_tol),
                                                    float(sew_tol)) * 4.0))
    if touching_solids and not allow_nonmanifold:
        raise AssemblyError(
            f"union of {len(solids)} solids that touch geometrically: "
            f"non-manifold result refused",
            kind="NonManifoldResult")
    notes = []
    if touching_solids:
        notes.append(
            f"non-manifold contact accepted: {len(solids)} touching "
            f"solids returned as a compound (allow_nonmanifold=True)")
    if not BRepCheck_Analyzer(result_shape, True).IsValid():
        raise AssemblyError("final assembled result is B-rep invalid",
                            kind="SolidInvalid")
    volume = float(sum(s.volume for s in solids))
    edge_lineage = _build_edge_lineage(
        result_shape, selected, split, model_a, model_b, float(base_tol))
    section_payloads = _build_section_payloads(split, edge_lineage)

    return BooleanAssemblyResult(
        operation=operation,
        decisions=decisions,
        selected_faces=len(selected),
        sewed_shape=sewed,
        shells=shell_records,
        solids=solids,
        shape=result_shape,
        volume=volume,
        free_edges=free,
        multiple_edges=multi,
        edge_lineage=edge_lineage,
        section_payloads=section_payloads,
        notes=notes,
        region_stats=region_stats,
    )
