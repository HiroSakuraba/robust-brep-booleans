"""STEP/B-rep ingest that preserves Solid -> Shell -> Face provenance.

This is the bridge from the current anonymous-mesh prototype to Tier B/C.
It deliberately keeps OCCT topology instead of flattening STEP into one
triangle soup and trying to rediscover cavities/provenance later.

The face broad phase is two-tier:
  1. precise OCCT face AABBs remove obviously disjoint face pairs;
  2. BSpline/NURBS pairs are refined by conservative knot-span control-hull
     AABBs from freeform.NurbsPatchIndex.
Only surviving pairs need expensive surface/surface intersection or fine
meshing.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .freeform import FreeformFaceAccel, FreeformError, candidate_patch_pairs


@dataclass
class FaceRecord:
    face_id: int
    solid_id: int
    shell_id: int
    face: object
    orientation: str
    surface_type: str
    uv_bounds: tuple[float, float, float, float]
    bbox_lo: np.ndarray
    bbox_hi: np.ndarray
    freeform: Optional[FreeformFaceAccel] = None


@dataclass
class ShellRecord:
    shell_id: int
    solid_id: int
    shell: object
    orientation: str
    face_ids: list[int] = field(default_factory=list)


@dataclass
class SolidRecord:
    solid_id: int
    solid: object
    shell_ids: list[int] = field(default_factory=list)


@dataclass
class BRepModel:
    shape: object
    solids: list[SolidRecord]
    shells: list[ShellRecord]
    faces: list[FaceRecord]

    @property
    def nurbs_faces(self) -> list[FaceRecord]:
        return [f for f in self.faces if f.freeform is not None]


def _shape_bbox(shape) -> tuple[np.ndarray, np.ndarray]:
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib
    b = Bnd_Box()
    # Geometry-aware bounding; does not require us to create a fine
    # triangulation merely to decide that two faces are far apart.
    BRepBndLib.AddOptimal_s(shape, b, False, False)
    if b.IsVoid():
        z = np.zeros(3, dtype=np.float64)
        return z, z
    # OCCT 8 changed Bnd_Box.Get() to return Bnd_Box::Limits, which the
    # Python binding does not currently convert. CornerMin/CornerMax are
    # stable gp_Pnt accessors in both OCP 7.x and 8.x.
    p0 = b.CornerMin()
    p1 = b.CornerMax()
    return (np.array([p0.X(), p0.Y(), p0.Z()], dtype=np.float64),
            np.array([p1.X(), p1.Y(), p1.Z()], dtype=np.float64))


def _surface_type_name(face) -> str:
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    t = BRepAdaptor_Surface(face).GetType()
    return getattr(t, "name", str(t))


def index_shape(shape, *, build_freeform: bool = True,
                trim_tol: float = 1e-8) -> BRepModel:
    """Preserve OCCT Solid->Shell->Face hierarchy and per-face provenance."""
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopAbs import TopAbs_SOLID, TopAbs_SHELL, TopAbs_FACE
    from OCP.TopoDS import TopoDS
    from OCP.BRepTools import BRepTools

    solids: list[SolidRecord] = []
    shells: list[ShellRecord] = []
    faces: list[FaceRecord] = []

    es = TopExp_Explorer(shape, TopAbs_SOLID)
    solid_id = shell_id = face_id = 0
    while es.More():
        solid = TopoDS.Solid(es.Current())
        sr = SolidRecord(solid_id, solid)
        esh = TopExp_Explorer(solid, TopAbs_SHELL)
        while esh.More():
            shell = TopoDS.Shell(esh.Current())
            shr = ShellRecord(shell_id, solid_id, shell,
                              getattr(shell.Orientation(), "name",
                                      str(shell.Orientation())))
            ef = TopExp_Explorer(shell, TopAbs_FACE)
            while ef.More():
                face = TopoDS.Face(ef.Current())
                lo, hi = _shape_bbox(face)
                st = _surface_type_name(face)
                ff = None
                if build_freeform and "BSpline" in st:
                    try:
                        ff = FreeformFaceAccel.from_occt_face(
                            face, trim_tol=trim_tol)
                    except FreeformError:
                        # Provenance remains even if acceleration cannot be
                        # built; callers can fall back to OCCT directly.
                        pass
                faces.append(FaceRecord(
                    face_id, solid_id, shell_id, face,
                    getattr(face.Orientation(), "name", str(face.Orientation())),
                    st, tuple(float(x) for x in BRepTools.UVBounds_s(face)),
                    lo, hi, ff))
                shr.face_ids.append(face_id)
                face_id += 1
                ef.Next()
            shells.append(shr)
            sr.shell_ids.append(shell_id)
            shell_id += 1
            esh.Next()
        solids.append(sr)
        solid_id += 1
        es.Next()

    # Do not silently discard shell-only STEP content.
    if not solids:
        esh = TopExp_Explorer(shape, TopAbs_SHELL)
        while esh.More():
            shell = TopoDS.Shell(esh.Current())
            shr = ShellRecord(shell_id, -1, shell,
                              getattr(shell.Orientation(), "name",
                                      str(shell.Orientation())))
            ef = TopExp_Explorer(shell, TopAbs_FACE)
            while ef.More():
                face = TopoDS.Face(ef.Current())
                lo, hi = _shape_bbox(face)
                st = _surface_type_name(face)
                ff = None
                if build_freeform and "BSpline" in st:
                    try:
                        ff = FreeformFaceAccel.from_occt_face(
                            face, trim_tol=trim_tol)
                    except FreeformError:
                        pass
                faces.append(FaceRecord(
                    face_id, -1, shell_id, face,
                    getattr(face.Orientation(), "name", str(face.Orientation())),
                    st, tuple(float(x) for x in BRepTools.UVBounds_s(face)),
                    lo, hi, ff))
                shr.face_ids.append(face_id)
                face_id += 1
                ef.Next()
            shells.append(shr)
            shell_id += 1
            esh.Next()

    return BRepModel(shape, solids, shells, faces)


def load_step(path: str, *, build_freeform: bool = True,
              trim_tol: float = 1e-8) -> BRepModel:
    """Read STEP with OCCT and preserve its B-rep hierarchy."""
    from OCP.STEPControl import STEPControl_Reader
    from OCP.IFSelect import IFSelect_RetDone

    r = STEPControl_Reader()
    stat = r.ReadFile(path)
    if stat != IFSelect_RetDone:
        raise FreeformError(f"STEP read failed with status {stat}",
                            "StepReadFailed")
    if r.TransferRoots() == 0:
        raise FreeformError("STEP contains no transferable roots",
                            "StepReadFailed")
    return index_shape(r.OneShape(), build_freeform=build_freeform,
                       trim_tol=trim_tol)


def _face_pairs_aabb(a: BRepModel, b: BRepModel,
                     pad: float) -> list[tuple[int, int]]:
    """Sweep-and-prune on precise face AABBs."""
    if not a.faces or not b.faces:
        return []
    alo = np.vstack([f.bbox_lo for f in a.faces]) - pad
    ahi = np.vstack([f.bbox_hi for f in a.faces]) + pad
    blo = np.vstack([f.bbox_lo for f in b.faces]) - pad
    bhi = np.vstack([f.bbox_hi for f in b.faces]) + pad
    oa = np.argsort(alo[:, 0], kind="mergesort")
    ob = np.argsort(blo[:, 0], kind="mergesort")
    active: list[int] = []
    jb = 0
    pairs = []
    for ia in oa:
        xmin, xmax = alo[ia, 0], ahi[ia, 0]
        while jb < len(ob) and blo[ob[jb], 0] <= xmax:
            active.append(int(ob[jb]))
            jb += 1
        # A minima are monotone, but A maxima are not.  Therefore
        # only expire B intervals that end before the current xmin.  A B
        # interval whose start is beyond this *particular* xmax may still
        # overlap a later, wider A interval and must remain active.
        active = [j for j in active if bhi[j, 0] >= xmin]
        for j in active:
            if (blo[j, 0] <= xmax and bhi[j, 0] >= xmin
                    and alo[ia, 1] <= bhi[j, 1] and ahi[ia, 1] >= blo[j, 1]
                    and alo[ia, 2] <= bhi[j, 2]
                    and ahi[ia, 2] >= blo[j, 2]):
                pairs.append((int(ia), int(j)))
    return pairs


def candidate_face_pairs(a: BRepModel, b: BRepModel,
                         pad: float = 0.0) -> list[dict]:
    """Conservative face/patch interaction candidates.

    NURBS pairs get a second conservative control-hull filter. Analytic or
    unsupported freeform pairs remain face-level candidates, so accuracy is
    never traded away for speed.
    """
    out = []
    for ia, ib in _face_pairs_aabb(a, b, float(pad)):
        fa, fb = a.faces[ia], b.faces[ib]
        patch_pairs = None
        if fa.freeform is not None and fb.freeform is not None:
            pp = candidate_patch_pairs(
                fa.freeform.index, fb.freeform.index, pad=float(pad))
            if not pp:
                continue
            patch_pairs = pp
        out.append({"face_a": fa.face_id, "face_b": fb.face_id,
                    "surface_a": fa.surface_type,
                    "surface_b": fb.surface_type,
                    "patch_pairs": patch_pairs})
    return out
