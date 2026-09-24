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


def _face_point(face, tol: float) -> np.ndarray:
    """Find a deterministic point strictly in a trimmed face."""
    from OCP.BRepClass3d import BRepClass3d_SolidExplorer
    from OCP.gp import gp_Pnt

    p = gp_Pnt()
    try:
        if BRepClass3d_SolidExplorer.FindAPointInTheFace_s(face, p):
            return _p3(p)
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
    frac = (0.5, 0.25, 0.75, 0.125, 0.875,
            0.375, 0.625, 0.0625, 0.9375,
            0.1875, 0.3125, 0.4375, 0.5625,
            0.6875, 0.8125)
    for fu in frac:
        u = u0 + (u1 - u0) * fu
        for fv in frac:
            v = v0 + (v1 - v0) * fv
            c = BRepClass_FaceClassifier(
                face, gp_Pnt2d(float(u), float(v)), float(tol), True)
            if c.State() == TopAbs_IN:
                return _p3(surf.Value(float(u), float(v)))
    raise AssemblyError("could not find a stable point inside trimmed face",
                        kind="NoInteriorFaceWitness")


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
                p = _face_point(piece.face, tol)
                cls = _classify_point_in_model(p, other, tol)
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
                    witness_xyz=p,
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
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    g = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, g)
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


def _solid_interior_point(solid, tol: float) -> np.ndarray:
    """Find a deterministic point classified IN a valid closed solid."""
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    from OCP.TopAbs import TopAbs_FACE, TopAbs_IN
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.BRepTools import BRepTools
    from OCP.gp import gp_Pnt, gp_Vec

    clf = BRepClass3d_SolidClassifier(solid)

    def is_in(x: np.ndarray) -> bool:
        clf.Perform(gp_Pnt(float(x[0]), float(x[1]), float(x[2])),
                    float(tol))
        return clf.State() == TopAbs_IN

    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(solid, props)
    cm = _p3(props.CentreOfMass())
    if is_in(cm):
        return cm

    lo, hi = _bbox(solid)
    scale = max(float(np.linalg.norm(hi - lo)), 1.0)

    ex = TopExp_Explorer(solid, TopAbs_FACE)
    while ex.More():
        face = TopoDS.Face(ex.Current())
        try:
            from OCP.BRepClass import BRepClass_FaceClassifier
            from OCP.gp import gp_Pnt2d
            u0, u1, v0, v1 = map(float, BRepTools.UVBounds_s(face))
            u, v = 0.5 * (u0 + u1), 0.5 * (v0 + v1)
            fc = BRepClass_FaceClassifier(
                face, gp_Pnt2d(u, v), float(tol), True)
            if fc.State() == TopAbs_IN:
                s = BRepAdaptor_Surface(face)
                p = gp_Pnt()
                du, dv = gp_Vec(), gp_Vec()
                s.D1(u, v, p, du, dv)
                n = np.cross(
                    np.array([du.X(), du.Y(), du.Z()]),
                    np.array([dv.X(), dv.Y(), dv.Z()]))
                nn = float(np.linalg.norm(n))
                if nn > 1e-300:
                    n /= nn
                    x = _p3(p)
                    for eps in (1e-4, 1e-5, 1e-6, 1e-7):
                        d = max(8.0 * tol, eps * scale)
                        for sign in (-1.0, 1.0):
                            q = x + sign * d * n
                            if is_in(q):
                                return q
        except Exception:
            pass
        ex.Next()

    frac = (0.5, 0.25, 0.75, 0.125, 0.875,
            0.375, 0.625, 0.0625, 0.9375)
    for fx in frac:
        x = lo[0] + (hi[0] - lo[0]) * fx
        for fy in frac:
            y = lo[1] + (hi[1] - lo[1]) * fy
            for fz in frac:
                z = lo[2] + (hi[2] - lo[2]) * fz
                q = np.array([x, y, z], dtype=np.float64)
                if is_in(q):
                    return q
    raise AssemblyError("could not find interior point of assembled shell",
                        kind="NoSolidInteriorWitness")


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
        p = _solid_interior_point(solid, tol)
        tmp.append({
            "index": i, "shell": sh, "outward": outward,
            "solid": solid, "volume": vol, "point": p,
            "containers": []
        })

    for child in tmp:
        p = gp_Pnt(*map(float, child["point"]))
        for parent in tmp:
            if parent is child:
                continue
            c = BRepClass3d_SolidClassifier(parent["solid"])
            c.Perform(p, float(tol))
            if c.State() == TopAbs_IN:
                child["containers"].append(parent["index"])

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
    )
