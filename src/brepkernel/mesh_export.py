"""B-rep -> triangle mesh -> glTF 2.0 binary (.glb), for game engines.

The mesh keeps what a B-rep knows and a triangle soup does not:

* hard edges where B-rep faces meet: every face gets its own vertices, so
  shading is smooth across a curved face and sharp across a face boundary;
* normals from the exact surfaces (not averaged triangle normals), flipped
  for reversed faces so they point out of the material;
* watertightness: OCCT's mesher discretizes each edge once and both faces
  use the same points, so after welding positions every edge is shared by
  exactly two triangles. tessellate() checks this and reports it;
* face identity: triangles are grouped by B-rep face, and the .glb records
  each face's index range (mesh.extras.brepkernel.faces) so an editor
  selection can be mapped back to the exact face.

Output follows the glTF specification: meters, +Y up, counter-clockwise
front faces. Shapes are expected to be in that frame already; pass
scale / from_z_up to convert CAD data (millimetres, +Z up) on the way.
"""
from __future__ import annotations

import json
import math
import struct
from dataclasses import dataclass, field

import numpy as np

from . import _occt_compat  # noqa: F401


class MeshExportError(ValueError):
    """The shape could not be meshed faithfully."""


@dataclass
class MeshData:
    positions: np.ndarray            # (n, 3) float64, output frame
    normals: np.ndarray              # (n, 3) float64, unit
    indices: np.ndarray              # (m, 3) uint32, counter-clockwise outward
    face_ranges: list = field(default_factory=list)   # per B-rep face: (first_tri, tri_count)
    linear_deflection: float = 0.0
    angular_deflection: float = 0.0
    watertight: bool = False
    boundary_edges: int = 0
    nonmanifold_edges: int = 0
    inconsistent_edges: int = 0
    mesh_volume: float = 0.0
    exact_volume: float | None = None

    @property
    def triangle_count(self) -> int:
        return int(self.indices.shape[0])

    @property
    def vertex_count(self) -> int:
        return int(self.positions.shape[0])

    def summary(self) -> dict:
        out = {
            "vertices": self.vertex_count,
            "triangles": self.triangle_count,
            "brep_faces": len(self.face_ranges),
            "linear_deflection": self.linear_deflection,
            "angular_deflection": self.angular_deflection,
            "watertight": self.watertight,
            "boundary_edges": self.boundary_edges,
            "nonmanifold_edges": self.nonmanifold_edges,
            "inconsistent_edges": self.inconsistent_edges,
            "mesh_volume": self.mesh_volume,
        }
        if self.exact_volume is not None:
            out["exact_volume"] = self.exact_volume
            out["volume_relative_error"] = (
                abs(self.mesh_volume - self.exact_volume) / max(abs(self.exact_volume), 1e-300))
        return out


def frame_transform(scale: float = 1.0, from_z_up: bool = False) -> np.ndarray:
    """3x3 linear map from shape coordinates to glTF (+Y up) coordinates."""
    m = np.eye(3) * float(scale)
    if from_z_up:
        # +Z up -> +Y up, right-handed: (x, y, z) -> (x, z, -y)
        m = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]]) @ m
    return m


def default_deflection(shape) -> float:
    """1/500 of the bounding-box diagonal, floored at 1e-6."""
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib
    box = Bnd_Box()
    BRepBndLib.Add_s(shape, box, False)
    if box.IsVoid():
        raise MeshExportError("shape is empty")
    lo, hi = box.CornerMin(), box.CornerMax()
    diag = math.dist((lo.X(), lo.Y(), lo.Z()), (hi.X(), hi.Y(), hi.Z()))
    return max(1e-6, diag / 500.0)


def _faces(shape):
    from OCP.TopAbs import TopAbs_FACE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS
    out, seen = [], {}
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        f = TopoDS.Face(ex.Current())
        ex.Next()
        bucket = seen.setdefault(hash(f), [])     # OCP hashes ignore orientation
        if any(f.IsSame(g) and f.Orientation() == g.Orientation() for g in bucket):
            continue
        bucket.append(f)
        out.append(f)
    return out


def tessellate(shape, *, linear_deflection: float | None = None,
               angular_deflection: float = 0.35, scale: float = 1.0,
               from_z_up: bool = False, exact_volume: float | None = None) -> MeshData:
    """Triangulate every face of `shape` (meshed on a copy; the input is untouched)."""
    from OCP.BRep import BRep_Tool
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Copy
    from OCP.BRepLProp import BRepLProp_SLProps
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.BRepMesh import BRepMesh_IncrementalMesh
    from OCP.TopAbs import TopAbs_REVERSED
    from OCP.TopLoc import TopLoc_Location

    if shape is None or shape.IsNull():
        raise MeshExportError("shape is null")
    lin = default_deflection(shape) if linear_deflection is None else float(linear_deflection)
    if not lin > 0 or not angular_deflection > 0:
        raise MeshExportError("deflections must be positive")
    work = BRepBuilderAPI_Copy(shape, True, False).Shape()
    mesher = BRepMesh_IncrementalMesh(work, lin, False, float(angular_deflection), True)
    if not mesher.IsDone():
        raise MeshExportError("OCCT mesher failed")
    M = frame_transform(scale, from_z_up)
    normal_M = np.linalg.inv(M).T

    positions, normals, tris, ranges = [], [], [], []
    for face in _faces(work):
        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation_s(face, loc)
        if tri is None:
            raise MeshExportError("a face has no triangulation")
        trsf = loc.Transformation()
        reversed_ = face.Orientation() == TopAbs_REVERSED
        n = tri.NbNodes()
        base = len(positions)
        surf = BRepAdaptor_Surface(face, True)
        props = BRepLProp_SLProps(surf, 1, 1e-9)
        has_uv = tri.HasUVNodes()
        pts = []
        for i in range(1, n + 1):
            p = tri.Node(i).Transformed(trsf)
            pts.append((p.X(), p.Y(), p.Z()))
        face_normals = [None] * n
        if has_uv:
            for i in range(1, n + 1):
                uv = tri.UVNode(i)
                props.SetParameters(uv.X(), uv.Y())
                if props.IsNormalDefined():
                    d = props.Normal()
                    v = np.array([d.X(), d.Y(), d.Z()])
                    face_normals[i - 1] = -v if reversed_ else v
        first_tri = len(tris)
        for k in range(1, tri.NbTriangles() + 1):
            a, b, c = tri.Triangle(k).Get()
            if reversed_:
                b, c = c, b
            pa, pb, pc = (np.array(pts[j - 1]) for j in (a, b, c))
            cr = np.cross(pb - pa, pc - pa)
            if not np.linalg.norm(cr) > 0:
                continue                      # degenerate (pole) triangle
            tris.append((base + a - 1, base + b - 1, base + c - 1))
            for j in (a, b, c):                 # fill singular-point normals
                if face_normals[j - 1] is None:
                    face_normals[j - 1] = cr / np.linalg.norm(cr)
        for i in range(n):
            positions.append(pts[i])
            normals.append(face_normals[i] if face_normals[i] is not None else np.zeros(3))
        ranges.append((first_tri, len(tris) - first_tri))

    if not tris:
        raise MeshExportError("no triangles produced")
    P = np.asarray(positions, dtype=np.float64) @ M.T
    N = np.asarray(normals, dtype=np.float64) @ normal_M.T
    lens = np.linalg.norm(N, axis=1)
    bad = lens < 1e-12
    if bad.any():
        raise MeshExportError(f"{int(bad.sum())} vertices without a defined normal")
    N = N / lens[:, None]
    T = np.asarray(tris, dtype=np.uint32)
    if np.linalg.det(M) < 0:                       # mirrored frame keeps CCW outward
        T = T[:, [0, 2, 1]]

    data = MeshData(positions=P, normals=N, indices=T, face_ranges=ranges,
                    linear_deflection=lin, angular_deflection=float(angular_deflection))
    _topology_check(data)
    data.mesh_volume = _signed_volume(P, T)
    if exact_volume is not None:
        data.exact_volume = float(exact_volume) * abs(np.linalg.det(M))
    return data


def _weld_keys(P: np.ndarray) -> np.ndarray:
    """Vertex -> welded id. OCCT gives both faces of an edge the same 3D
    points; the quantum only absorbs the last-bit noise of the transform."""
    span = float(np.ptp(P, axis=0).max()) if len(P) else 1.0
    q = max(span, 1.0) * 1e-9
    keys = np.round(P / q).astype(np.int64)
    _, inv = np.unique(keys, axis=0, return_inverse=True)
    return inv.reshape(-1)


def _topology_check(d: MeshData) -> None:
    w = _weld_keys(d.positions)
    T = w[d.indices.astype(np.int64)]
    counts: dict = {}
    for a, b, c in T:
        for u, v in ((a, b), (b, c), (c, a)):
            key = (u, v) if u < v else (v, u)
            sgn = 1 if u < v else -1
            rec = counts.setdefault(key, [0, 0])
            rec[0] += 1
            rec[1] += sgn
    boundary = sum(1 for n, _s in counts.values() if n == 1)
    nonman = sum(1 for n, _s in counts.values() if n > 2)
    inconsistent = sum(1 for n, s in counts.values() if n == 2 and s != 0)
    d.boundary_edges, d.nonmanifold_edges, d.inconsistent_edges = boundary, nonman, inconsistent
    d.watertight = boundary == 0 and nonman == 0 and inconsistent == 0


def _signed_volume(P: np.ndarray, T: np.ndarray) -> float:
    a, b, c = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


# ---------------------------------------------------------------------------
# glTF 2.0 binary
# ---------------------------------------------------------------------------

def _pad(b: bytes, fill: bytes) -> bytes:
    return b + fill * ((4 - len(b) % 4) % 4)


def to_glb(mesh: MeshData, *, name: str = "brepkernel", base_color=(0.72, 0.72, 0.74, 1.0),
           metallic: float = 0.0, roughness: float = 0.5, extras: dict | None = None) -> bytes:
    """One node, one mesh, one primitive (triangles grouped by B-rep face)."""
    pos = np.ascontiguousarray(mesh.positions, dtype="<f4")
    nor = np.ascontiguousarray(mesh.normals, dtype="<f4")
    idx = np.ascontiguousarray(mesh.indices.reshape(-1), dtype="<u4")
    blobs = [pos.tobytes(), nor.tobytes(), idx.tobytes()]
    views, offset = [], 0
    for blob, target in zip(blobs, (34962, 34962, 34963)):
        views.append({"buffer": 0, "byteOffset": offset, "byteLength": len(blob),
                      "target": target})
        offset += len(blob)
        offset += (4 - offset % 4) % 4
    binary = b"".join(_pad(b, b"\x00") for b in blobs)
    faces = [[int(first * 3), int(count * 3)] for first, count in mesh.face_ranges]
    meta = {"version": 1, "faces": faces, **mesh.summary(), **(extras or {})}
    doc = {
        "asset": {"version": "2.0", "generator": "brepkernel mesh_export"},
        "scene": 0,
        "scenes": [{"nodes": [0]}],
        "nodes": [{"mesh": 0, "name": name}],
        "meshes": [{"name": name, "primitives": [{
            "attributes": {"POSITION": 0, "NORMAL": 1},
            "indices": 2, "material": 0, "mode": 4}],
            "extras": {"brepkernel": meta}}],
        "materials": [{"name": f"{name}_material", "pbrMetallicRoughness": {
            "baseColorFactor": [float(x) for x in base_color],
            "metallicFactor": float(metallic), "roughnessFactor": float(roughness)}}],
        "buffers": [{"byteLength": len(binary)}],
        "bufferViews": views,
        "accessors": [
            {"bufferView": 0, "componentType": 5126, "count": int(pos.shape[0]), "type": "VEC3",
             "min": [float(x) for x in pos.min(axis=0)], "max": [float(x) for x in pos.max(axis=0)]},
            {"bufferView": 1, "componentType": 5126, "count": int(nor.shape[0]), "type": "VEC3"},
            {"bufferView": 2, "componentType": 5125, "count": int(idx.shape[0]), "type": "SCALAR"},
        ],
    }
    js = _pad(json.dumps(doc, separators=(",", ":")).encode("utf-8"), b" ")
    total = 12 + 8 + len(js) + 8 + len(binary)
    return (struct.pack("<III", 0x46546C67, 2, total)
            + struct.pack("<II", len(js), 0x4E4F534A) + js
            + struct.pack("<II", len(binary), 0x004E4942) + binary)


def write_glb(shape, path: str, *, require_watertight: bool = True, **kw) -> MeshData:
    """Tessellate `shape` and write `path` (.glb). Returns the MeshData.

    With require_watertight (default) a mesh with open, non-manifold or
    inconsistently wound edges raises MeshExportError and nothing is
    written: game engines voxelize and collide against these meshes.
    """
    glb_kw = {k: kw.pop(k) for k in ("name", "base_color", "metallic", "roughness", "extras")
              if k in kw}
    mesh = tessellate(shape, **kw)
    if require_watertight and not mesh.watertight:
        raise MeshExportError(
            f"mesh is not watertight (boundary={mesh.boundary_edges}, "
            f"nonmanifold={mesh.nonmanifold_edges}, inconsistent={mesh.inconsistent_edges})")
    with open(path, "wb") as fh:
        fh.write(to_glb(mesh, **glb_kw))
    return mesh


def read_glb(path: str):
    """(json document, positions, normals, indices) of a .glb this module wrote."""
    raw = open(path, "rb").read()
    magic, version, total = struct.unpack_from("<III", raw, 0)
    if magic != 0x46546C67 or version != 2 or total != len(raw):
        raise MeshExportError("not a glTF 2.0 binary")
    jl, jt = struct.unpack_from("<II", raw, 12)
    doc = json.loads(raw[20:20 + jl])
    bl, bt = struct.unpack_from("<II", raw, 20 + jl)
    binary = raw[28 + jl:28 + jl + bl]

    def acc(i, dtype, width):
        a = doc["accessors"][i]
        v = doc["bufferViews"][a["bufferView"]]
        arr = np.frombuffer(binary, dtype=dtype, count=a["count"] * width,
                            offset=v["byteOffset"])
        return arr.reshape(a["count"], width) if width > 1 else arr
    prim = doc["meshes"][0]["primitives"][0]
    return (doc, acc(prim["attributes"]["POSITION"], "<f4", 3),
            acc(prim["attributes"]["NORMAL"], "<f4", 3),
            acc(prim["indices"], "<u4", 1).reshape(-1, 3))
