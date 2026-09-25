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

def _classify_point_in_model(point: np.ndarray, model: BRepModel,
                             tol: float) -> str:
    """Classify a point against the union of the model's OCCT solids."""
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.TopAbs import TopAbs_IN, TopAbs_ON, TopAbs_OUT
    from OCP.gp import gp_Pnt

    if not model.solids:
        raise AssemblyError(
            "material patch classification requires closed OCCT solids; "
            "shell-only input is unsupported at this stage",
            kind="ShellOnlyClassificationUnsupported")

    p = gp_Pnt(float(point[0]), float(point[1]), float(point[2]))
    saw_on = False
    for sr in model.solids:
        c = BRepClass3d_SolidClassifier(sr.solid)
        c.Perform(p, float(tol))
        st = c.State()
        if st == TopAbs_IN:
            return "inside"
        if st == TopAbs_ON:
            saw_on = True
        elif st != TopAbs_OUT:
            return "unknown"
    return "boundary" if saw_on else "outside"


def _decision_rule(operation: str, operand: str,
                   classification: str) -> tuple[bool, bool]:
    """Return (keep, reverse) for a non-boundary patch."""
    if classification not in ("inside", "outside"):
        raise AssemblyError(
            f"cannot decide {operation} patch on {classification} witness",
            kind="BoundaryOrUnknownPatch")
    if operation == "union":
        return classification == "outside", False
    if operation == "intersection":
        return classification == "inside", False
    if operation == "difference":
        if operand == "A":
            return classification == "outside", False
        return classification == "inside", True
    raise ValueError(f"unsupported boolean operation {operation!r}")


def _reverse_face(face):
    from OCP.TopoDS import TopoDS
    return TopoDS.Face(face.Reversed())


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
        for fr in groups:
            for piece in fr.pieces:
                tol = max(
                    float(base_tol),
                    2.0 * float(BRep_Tool.Tolerance_s(piece.face)))
                points = _face_points(piece.face, tol)
                classes = tuple(
                    _classify_point_in_model(p, other, tol)
                    for p in points)
                invalid = sorted(set(
                    x for x in classes
                    if x not in ("inside", "outside")))
                if invalid:
                    raise AssemblyError(
                        f"{operand} face {piece.parent_face_id} piece "
                        f"{piece.piece_index}: witness classification "
                        f"contains {invalid}",
                        kind="BoundaryOrUnknownPatch")
                unique = set(classes)
                if len(unique) != 1:
                    raise AssemblyError(
                        f"{operand} face {piece.parent_face_id} piece "
                        f"{piece.piece_index}: supposedly split patch "
                        f"straddles material states {sorted(unique)}",
                        kind="PatchClassificationInconsistent")
                cls = classes[0]
                keep, rev = _decision_rule(operation, operand, cls)
                source = piece.face
                selected = _reverse_face(source) if keep and rev else (
                    source if keep else None)
                out.append(PatchDecision(
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
                ))

    one_side("A", split.faces_a, model_b)
    one_side("B", split.faces_b, model_a)
    return out


def _empty_compound():
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Compound
    c = TopoDS_Compound()
    b = BRep_Builder()
    b.MakeCompound(c)
    return c


def _shape_volume(shape) -> float:
    """Adaptive volume measurement, including B-spline span integration."""
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    g = GProp_GProps()
    err = BRepGProp.VolumePropertiesGK_s(
        shape, g, 1e-10, True, True, False, False, False)
    if float(err) < 0.0:
        raise AssemblyError("adaptive volume integration failed",
                            kind="VolumeIntegrationFailed")
    return float(g.Mass())


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

    def is_in(x: np.ndarray) -> bool:
        clf.Perform(gp_Pnt(float(x[0]), float(x[1]), float(x[2])),
                    float(tol))
        return clf.State() == TopAbs_IN

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
    props = GProp_GProps()
    err = BRepGProp.VolumePropertiesGK_s(
        solid, props, 1e-9, True, True, True, False, False)
    if float(err) < 0.0:
        raise AssemblyError("adaptive center-of-mass integration failed",
                            kind="VolumeIntegrationFailed")
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
        vol = abs(_shape_volume(solid))
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
    for child in tmp:
        for parent in tmp:
            if parent is child:
                continue
            classifier = BRepClass3d_SolidClassifier(parent["solid"])
            states = []
            for q in child["points"]:
                classifier.Perform(
                    gp_Pnt(float(q[0]), float(q[1]), float(q[2])),
                    float(tol))
                states.append(classifier.State())

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
        vol = _shape_volume(solid)
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


def _edge_matches_section(edge, sec, base_tol: float) -> bool:
    """Recognize a result edge as lying on a verified section edge.

    Splitter/sewing may preserve the edge, copy it, or keep only a sub-edge.
    Therefore TShape identity is tried first, then a conservative geometric
    sub-edge check: the result edge may not exceed the section length and
    several points along it must lie on the verified section curve.
    """
    if edge.IsSame(sec.edge):
        return True

    from OCP.BRepAdaptor import BRepAdaptor_Curve

    le = _edge_length(edge)
    ls = _edge_length(sec.edge)
    tol = max(float(base_tol), float(sec.verify_tolerance),
              float(sec.edge_tolerance))
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
        if _point_edge_distance(x, sec.edge) > match_tol:
            return False
    return True


def _build_edge_lineage(result_shape, selected: list[PatchDecision],
                        split: ModelSplitResult,
                        model_a: BRepModel, model_b: BRepModel,
                        base_tol: float) -> list[EdgeLineageRecord]:
    """Attach final result edges to selected patches and section evidence."""
    original = {
        "A": {f.face_id: f.face for f in model_a.faces},
        "B": {f.face_id: f.face for f in model_b.faces},
    }
    out = []
    for i, edge in enumerate(_unique_edges(result_shape)):
        refs = []
        parents = []
        operands = []
        source_boundary = []
        for d in selected:
            sf = d.sewed_face if d.sewed_face is not None else d.selected_face
            if sf is not None and _shape_has_edge(sf, edge):
                ref = (d.operand, d.parent_face_id, d.piece_index)
                if ref not in refs:
                    refs.append(ref)
                pf = (d.operand, d.parent_face_id)
                if pf not in parents:
                    parents.append(pf)
                if d.operand not in operands:
                    operands.append(d.operand)
                parent_face = original[d.operand].get(d.parent_face_id)
                if parent_face is not None and _shape_has_edge(parent_face, edge):
                    if pf not in source_boundary:
                        source_boundary.append(pf)

        intersections = []
        for sec in split.section_edges:
            if _edge_matches_section(edge, sec, float(base_tol)):
                key = (sec.face_a, sec.face_b, sec.edge_index)
                if key not in intersections:
                    intersections.append(key)

        if intersections:
            kind = "boolean_section"
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


def assemble_boolean(model_a: BRepModel, model_b: BRepModel,
                     split: ModelSplitResult, operation: str, *,
                     base_tol: float = 1e-7,
                     sew_tol: Optional[float] = None
                     ) -> BooleanAssemblyResult:
    """Classify exact B-rep patches and assemble union/intersection/A-B."""
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

    if not selected:
        empty = _empty_compound()
        return BooleanAssemblyResult(
            operation=operation, decisions=decisions, selected_faces=0,
            sewed_shape=empty, shells=[], solids=[], shape=empty,
            volume=0.0, free_edges=0, multiple_edges=0,
            edge_lineage=[],
            section_payloads=[],
            notes=["empty material result"])

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
    )
