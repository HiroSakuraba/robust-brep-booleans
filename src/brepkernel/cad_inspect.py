"""Exact answers about one B-rep for editors: what it is made of, and where a point lands on it.

Two read-only queries, used by DVE for Measure and exact snapping (``brepkernel describe`` and
``brepkernel query``). Everything is in the shape's own frame and units (for a project ``.brep``
written by ``brepkernel import``: meters, +Y up). Nothing is approximated by triangles.

describe_shape(shape)
    Faces in mesh order (the order :mod:`brepkernel.mesh_export` writes them, so a face index
    means the same face in both), each with its surface type, area and the defining numbers of
    that surface (plane origin and outward normal; axis, radius and whether a cylinder or cone
    is a hole; sphere and torus centres and radii; B-spline degrees). Plus the part's true
    corners, its edges (type, ends, arc-length midpoint, length) and its circles (centre, axis,
    radius, and whether the arcs on that circle close it), with the total volume and area.

    Seam edges (where a closed surface such as a cylinder meets itself) and degenerate edges
    (a sphere pole) are not features, so they are left out of ``edges``; a vertex is a corner
    only if an open, non-seam edge ends there (a cylinder has no corners).

query_points(shape, points)
    For each point: the closest point on the boundary, the distance, the face it lies on (mesh
    order) and the outward unit normal there. On an edge or a corner the face is the nearest of
    the faces that meet there, ties broken by mesh order.
"""
from __future__ import annotations

import math

from . import _occt_compat  # noqa: F401  (TopoDS cast spellings on OCP 7.8 and 8.0)

_SURFACE_TYPES = {
    "GeomAbs_Plane": "plane", "GeomAbs_Cylinder": "cylinder", "GeomAbs_Cone": "cone",
    "GeomAbs_Sphere": "sphere", "GeomAbs_Torus": "torus", "GeomAbs_BezierSurface": "bezier",
    "GeomAbs_BSplineSurface": "bspline", "GeomAbs_SurfaceOfRevolution": "revolution",
    "GeomAbs_SurfaceOfExtrusion": "extrusion", "GeomAbs_OffsetSurface": "offset",
}
_CURVE_TYPES = {"GeomAbs_Line": "line", "GeomAbs_Circle": "circle"}


def _xyz(p):
    return [float(p.X()), float(p.Y()), float(p.Z())]


def _type_name(enum_value, table):
    return table.get(getattr(enum_value, "name", str(enum_value).split(".")[-1]), None)


def _unit(v):
    n = math.sqrt(sum(c * c for c in v))
    return [c / n for c in v] if n > 0 else [0.0, 0.0, 0.0]


class _ShapeIndex:
    """Unique sub-shapes (orientation ignored) with the faces, in mesh order, that contain them.

    Built with explorers and IsSame so it works on every supported OCP (8.0 no longer exports
    the TopTools indexed-map typedefs)."""

    def __init__(self, faces, kind):
        from OCP.TopExp import TopExp_Explorer
        self.items, self.owners, self._buckets = [], [], {}
        for i, f in enumerate(faces):
            ex = TopExp_Explorer(f, kind)
            while ex.More():
                k = self._add(ex.Current())
                if i not in self.owners[k]:
                    self.owners[k].append(i)
                ex.Next()

    def _add(self, s):
        bucket = self._buckets.setdefault(hash(s), [])
        for k in bucket:
            if self.items[k].IsSame(s):
                return k
        self.items.append(s)
        self.owners.append([])
        bucket.append(len(self.items) - 1)
        return len(self.items) - 1

    def find(self, s):
        for k in self._buckets.get(hash(s), []):
            if self.items[k].IsSame(s):
                return k
        return -1


def _faces(shape):
    from .mesh_export import _faces as mesh_faces
    return mesh_faces(shape)


def _outward_normal(face, u, v):
    """Unit normal of the face at (u, v), pointing out of the material (face orientation applied)."""
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.gp import gp_Pnt, gp_Vec
    from OCP.TopAbs import TopAbs_REVERSED
    s = BRepAdaptor_Surface(face, False)
    p, du, dv = gp_Pnt(), gp_Vec(), gp_Vec()
    s.D1(u, v, p, du, dv)
    n = du.Crossed(dv)
    if n.Magnitude() <= 1e-300:
        # singular point (a cone apex, a sphere pole): step a little into the face
        u0, u1, v0, v1 = s.FirstUParameter(), s.LastUParameter(), s.FirstVParameter(), s.LastVParameter()
        um = u + 1e-6 * ((u0 + u1) / 2 - u) if math.isfinite(u0 + u1) else u
        vm = v + 1e-6 * ((v0 + v1) / 2 - v) if math.isfinite(v0 + v1) else v
        s.D1(um, vm, p, du, dv)
        n = du.Crossed(dv)
    out = _unit([n.X(), n.Y(), n.Z()])
    if face.Orientation() == TopAbs_REVERSED:
        out = [-c for c in out]
    return out


def _face_sample(face):
    """A (u, v) inside the face: the parameters of its centroid's projection."""
    from OCP.BRep import BRep_Tool
    from OCP.BRepGProp import BRepGProp
    from OCP.BRepTools import BRepTools
    from OCP.GProp import GProp_GProps
    from OCP.ShapeAnalysis import ShapeAnalysis_Surface
    props = GProp_GProps()
    BRepGProp.SurfaceProperties_s(face, props)
    c = props.CentreOfMass()
    sas = ShapeAnalysis_Surface(BRep_Tool.Surface_s(face))
    uv = sas.ValueOfUV(c, 1e-9)
    u, v = uv.X(), uv.Y()
    umin, umax, vmin, vmax = BRepTools.UVBounds_s(face)
    if not (umin <= u <= umax and vmin <= v <= vmax):        # centroid off the face (an annulus)
        u, v = 0.5 * (umin + umax), 0.5 * (vmin + vmax)
    return u, v, props


def _radial_hole(face, axis_point, axis_dir):
    """True when the material is on the axis side: the outward normal points toward the axis."""
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    u, v, _ = _face_sample(face)
    p = BRepAdaptor_Surface(face, False).Value(u, v)
    q = [p.X() - axis_point[0], p.Y() - axis_point[1], p.Z() - axis_point[2]]
    along = sum(a * b for a, b in zip(q, axis_dir))
    radial = [q[i] - along * axis_dir[i] for i in range(3)]
    n = _outward_normal(face, u, v)
    return sum(a * b for a, b in zip(n, radial)) < 0.0


def describe_face(face, index: int) -> dict:
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.TopAbs import TopAbs_REVERSED
    s = BRepAdaptor_Surface(face, True)
    kind = _type_name(s.GetType(), _SURFACE_TYPES) or "other"
    u, v, props = _face_sample(face)
    d = {"index": index, "type": kind, "area": float(props.Mass()),
         "reversed": face.Orientation() == TopAbs_REVERSED,
         "centroid": _xyz(props.CentreOfMass())}
    if kind == "plane":
        d["origin"] = d["centroid"]                        # on the plane, inside the face's span
        d["normal"] = _outward_normal(face, u, v)
    elif kind in ("cylinder", "cone"):
        g = s.Cylinder() if kind == "cylinder" else s.Cone()
        ax = g.Axis()
        d["axis_origin"] = _xyz(ax.Location())
        d["axis"] = _xyz(ax.Direction())
        d["radius"] = float(g.Radius() if kind == "cylinder" else g.RefRadius())
        if kind == "cone":
            d["half_angle"] = abs(float(g.SemiAngle()))
            d["apex"] = _xyz(g.Apex())
        d["hole"] = _radial_hole(face, d["axis_origin"], d["axis"])
    elif kind == "sphere":
        g = s.Sphere()
        d["centre"] = _xyz(g.Location())
        d["radius"] = float(g.Radius())
        d["hole"] = _radial_hole(face, d["centre"], [0.0, 0.0, 0.0])    # a spherical cavity
    elif kind == "torus":
        g = s.Torus()
        d["centre"] = _xyz(g.Location())
        d["axis"] = _xyz(g.Axis().Direction())
        d["major_radius"] = float(g.MajorRadius())
        d["minor_radius"] = float(g.MinorRadius())
    elif kind in ("bspline", "bezier"):
        d["u_degree"] = int(s.UDegree())
        d["v_degree"] = int(s.VDegree())
        d["rational"] = bool(s.IsURational() or s.IsVRational())
        d["poles"] = [int(s.NbUPoles()), int(s.NbVPoles())]
    return d


def _edge_record(edge):
    """(record, circle) for a feature edge; circle is (centre, axis, radius, length) or None."""
    from OCP.BRep import BRep_Tool
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    from OCP.GCPnts import GCPnts_AbscissaPoint
    from OCP.TopExp import TopExp
    c = BRepAdaptor_Curve(edge)
    kind = _type_name(c.GetType(), _CURVE_TYPES) or "curve"
    t0, t1 = c.FirstParameter(), c.LastParameter()
    length = float(GCPnts_AbscissaPoint.Length_s(c))
    if kind in ("line", "circle"):
        tm = 0.5 * (t0 + t1)
    else:
        tm = GCPnts_AbscissaPoint(c, 0.5 * length, t0).Parameter()
    v0, v1 = TopExp.FirstVertex_s(edge, True), TopExp.LastVertex_s(edge, True)
    rec = {"type": kind, "start": _xyz(BRep_Tool.Pnt_s(v0)), "end": _xyz(BRep_Tool.Pnt_s(v1)),
           "mid": _xyz(c.Value(tm)), "length": length}
    circle = None
    if kind == "circle":
        g = c.Circle()
        circle = (_xyz(g.Location()), _xyz(g.Axis().Direction()), float(g.Radius()), length)
    return rec, circle


def _merge_circles(found, size):
    """One entry per geometric circle; ``full`` when its arcs add up to the whole circle."""
    out = []
    tol = 1e-9 * max(size, 1e-300)
    for centre, axis, radius, length in found:
        for c in out:
            if (math.dist(c["centre"], centre) <= tol and abs(c["radius"] - radius) <= tol
                    and abs(abs(sum(a * b for a, b in zip(c["axis"], axis))) - 1.0) <= 1e-9):
                c["_length"] += length
                break
        else:
            out.append({"centre": centre, "axis": axis, "radius": radius, "_length": length})
    for c in out:
        c["full"] = c.pop("_length") >= 2.0 * math.pi * c["radius"] * (1.0 - 1e-9)
    return out


def describe_shape(shape) -> dict:
    from OCP.BRep import BRep_Tool
    from OCP.BRepBndLib import BRepBndLib
    from OCP.Bnd import Bnd_Box
    from OCP.TopAbs import TopAbs_EDGE, TopAbs_SOLID, TopAbs_VERTEX
    from OCP.TopExp import TopExp, TopExp_Explorer
    from OCP.TopoDS import TopoDS
    from .cli import _volume
    faces = _faces(shape)
    box = Bnd_Box()
    BRepBndLib.Add_s(shape, box)
    size = 1.0
    if not box.IsVoid():
        lo, hi = box.CornerMin(), box.CornerMax()
        size = math.dist(_xyz(lo), _xyz(hi)) or 1.0

    edge_index = _ShapeIndex(faces, TopAbs_EDGE)
    edges, circles, corners = [], [], _ShapeIndex([], TopAbs_VERTEX)
    seams = degenerate = 0
    for e, owners in zip(edge_index.items, edge_index.owners):
        e = TopoDS.Edge(e)
        if BRep_Tool.Degenerated_s(e):
            degenerate += 1
            continue
        if any(BRep_Tool.IsClosed_s(e, faces[i]) for i in owners):
            seams += 1
            continue
        rec, circle = _edge_record(e)
        edges.append(rec)
        if circle:
            circles.append(circle)
        v0, v1 = TopExp.FirstVertex_s(e), TopExp.LastVertex_s(e)
        if not v0.IsSame(v1):                              # an open edge: its ends are corners
            corners._add(v0)
            corners._add(v1)

    n_solids, ex = 0, TopExp_Explorer(shape, TopAbs_SOLID)
    while ex.More():
        n_solids += 1
        ex.Next()
    face_records = [describe_face(f, i) for i, f in enumerate(faces)]
    all_vertices = _ShapeIndex(faces, TopAbs_VERTEX)
    return {
        "solids": n_solids,
        "volume": _volume(shape),
        "area": float(sum(f["area"] for f in face_records)),
        "faces": face_records,
        "vertices": [_xyz(BRep_Tool.Pnt_s(TopoDS.Vertex(v))) for v in corners.items],
        "edges": edges,
        "circles": _merge_circles(circles, size),
        "skipped": {"seam_edges": seams, "degenerate_edges": degenerate,
                    "vertices_not_corners": len(all_vertices.items) - len(corners.items)},
        "frame": "the shape's own (a project .brep: meters, +Y up)",
    }


def query_points(shape, points) -> list[dict]:
    from OCP.BRep import BRep_Tool
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
    from OCP.BRepExtrema import BRepExtrema_DistShapeShape, BRepExtrema_SupportType
    from OCP.ShapeAnalysis import ShapeAnalysis_Surface
    from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_VERTEX
    from OCP.gp import gp_Pnt

    faces = _faces(shape)
    face_index = _ShapeIndex([], TopAbs_FACE)
    for f in faces:
        face_index._add(f)
    ancestors = {k: _ShapeIndex(faces, k) for k in (TopAbs_EDGE, TopAbs_VERTEX)}

    def index_of(face):
        return face_index.find(face)

    def nearest_face(vertex, p, candidates):
        """(index, distance, closest, uv, normal) of the nearest candidate face. Where several
        are equally near (an edge or a corner), the one whose outward normal points most
        toward p wins, then mesh order."""
        found = []
        for i in sorted(set(candidates)):
            d = BRepExtrema_DistShapeShape(vertex, faces[i])
            if not (d.IsDone() and d.NbSolution() > 0):
                continue
            q = d.PointOnShape2(1)
            uv = ShapeAnalysis_Surface(BRep_Tool.Surface_s(faces[i])).ValueOfUV(q, 1e-12)
            found.append((i, float(d.Value()), q, uv, _outward_normal(faces[i], uv.X(), uv.Y())))
        if not found:
            return None
        dmin = min(f[1] for f in found)
        near = [f for f in found if f[1] <= dmin + 1e-12 * max(1.0, dmin)]
        if len(near) == 1 or dmin <= 0.0:
            return near[0]
        def facing(f):
            q = f[2]
            r = _unit([p[0] - q.X(), p[1] - q.Y(), p[2] - q.Z()])
            return sum(a * b for a, b in zip(r, f[4]))
        return max(near, key=lambda f: (facing(f), -f[0]))

    out = []
    for p in points:
        rec = {"point": [float(c) for c in p], "found": False}
        out.append(rec)
        if not all(math.isfinite(c) for c in p) or not faces:
            continue
        vertex = BRepBuilderAPI_MakeVertex(gp_Pnt(*(float(c) for c in p))).Vertex()
        dist = BRepExtrema_DistShapeShape(vertex, shape)
        if not dist.IsDone() or dist.NbSolution() == 0:
            continue
        # the solution closest to the point, and among equals the first face in mesh order
        best = min(range(1, dist.NbSolution() + 1),
                   key=lambda k: dist.PointOnShape1(k).Distance(dist.PointOnShape2(k)))
        support = dist.SupportOnShape2(best)
        stype = dist.SupportTypeShape2(best)
        if stype == BRepExtrema_SupportType.BRepExtrema_IsInFace:
            cands = [index_of(support)]
        else:
            m = ancestors[TopAbs_EDGE if stype == BRepExtrema_SupportType.BRepExtrema_IsOnEdge else TopAbs_VERTEX]
            k = m.find(support)
            cands = list(m.owners[k]) if k >= 0 else []
        cands = [i for i in cands if i >= 0]
        chosen = nearest_face(vertex, [float(c) for c in p], cands) if cands else None
        if chosen is None:
            continue
        i, d, closest, uv, normal = chosen
        rec.update({"found": True, "closest": _xyz(closest), "distance": d, "face": i,
                    "normal": normal, "uv": [uv.X(), uv.Y()]})
    return out
