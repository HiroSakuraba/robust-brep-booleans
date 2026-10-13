"""OCCT TopoDS solids -> omniSolid BrepArray (.usda text).

Pure OCP + text output: no OpenUSD or usd-brep import is needed to write
the file. The layout follows the omniSolid schema shipped with NVIDIA's
usd-brep (source/schema/omniSolid/resources/schema.usda) and the files
usd-brep itself exports:

* one BrepArray prim per call; by default one Brep holding every solid,
  laid out as SMLib lays out separate lumps: the infinite void region
  (bounded by every outer shell), one solid region per solid, then one
  void region per cavity. Every face has two faceuses, one per side;
* faceuse "same" is the side the face's surface normal points to;
* loops are packed per face (outer loop first) and run counter-clockwise
  about the surface normal; edgeuses are packed per loop and "same" means
  the loop runs along the edge's curve direction. OCCT degenerate (pole)
  edges are dropped, as in usd-brep's own sphere export;
* planes, cylinders and cones (as BrepSurfaceConeAPI), spheres and tori
  are written analytically in the schema's parameterization, checked
  numerically against OCCT before writing; a left-handed OCCT frame, whose
  normal points inward, is written with sides swapped and loops reversed;
  every other surface becomes an exact rational B-spline;
* every edge is a B-spline (BrepCurve3dNurbAPI), circles and ellipses with
  at most a quarter turn per rational span, as usd-brep writes them;
* no UV trim curves are written; usd-brep rebuilds trims from the 3D edges.

Non-manifold input (an edge with other than two face uses), sheets and
wire edges are refused with OmniSolidWriteError rather than guessed.
See docs/usd_brep_backend.md for why each convention is what it is.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .. import _occt_compat  # noqa: F401  (TopoDS cast spellings on OCP 7.8)


class OmniSolidWriteError(ValueError):
    """The shape cannot be represented faithfully; nothing was written."""


# ---------------------------------------------------------------------------
# small identity helpers (OCP hashes ignore orientation, like IsSame)
# ---------------------------------------------------------------------------

class _ShapeIndex:
    def __init__(self):
        self._buckets: dict = {}
        self.items: list = []

    def index(self, shape, add: bool = True) -> int:
        bucket = self._buckets.setdefault(hash(shape), [])
        for s, i in bucket:
            if s.IsSame(shape):
                return i
        if not add:
            return -1
        i = len(self.items)
        bucket.append((shape, i))
        self.items.append(shape)
        return i


def _xyz(p) -> tuple:
    return (float(p.X()), float(p.Y()), float(p.Z()))


# ---------------------------------------------------------------------------
# geometry extraction
# ---------------------------------------------------------------------------

@dataclass
class _Nurb:
    cvs: list            # list of (x, y, z)
    weights: list
    knots: list          # full multiplicity-expanded knot vector
    order: int


@dataclass
class _NurbSurf:
    cvs: list            # row-major, U rows, V columns
    weights: list
    u_knots: list
    v_knots: list
    u_order: int
    v_order: int
    nu: int
    nv: int


def _bspline_curve_data(bs) -> _Nurb:
    if bs.IsPeriodic():
        bs = bs.Copy()
        bs.SetNotPeriodic()
    n = bs.NbPoles()
    cvs = [_xyz(bs.Pole(i)) for i in range(1, n + 1)]
    w = [float(bs.Weight(i)) for i in range(1, n + 1)]
    knots = []
    for i in range(1, bs.NbKnots() + 1):
        knots.extend([float(bs.Knot(i))] * int(bs.Multiplicity(i)))
    return _Nurb(cvs, w, knots, int(bs.Degree()) + 1)


def _edge_bspline(edge):
    """Exact B-spline of an edge's 3D curve over the edge's own range."""
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    from OCP.Geom import Geom_TrimmedCurve, Geom_BSplineCurve
    from OCP.GeomConvert import GeomConvert
    ad = BRepAdaptor_Curve(edge)
    first, last = float(ad.FirstParameter()), float(ad.LastParameter())
    geom = ad.Curve().Curve()
    geom = geom.Transformed(ad.Trsf())
    if isinstance(geom, Geom_BSplineCurve) and not geom.IsPeriodic():
        bs = geom.Copy()
        if abs(bs.FirstParameter() - first) > 0 or abs(bs.LastParameter() - last) > 0:
            bs.Segment(first, last)
    else:
        from OCP.Geom import Geom_Circle, Geom_Ellipse
        trimmed = Geom_TrimmedCurve(geom, first, last)
        if isinstance(geom, (Geom_Circle, Geom_Ellipse)):
            # One rational span per quarter turn at most, as usd-brep writes
            # circles itself. OCCT's default packs a full circle into three
            # 120-degree spans (weight 0.5), and SMLib's area and volume
            # integration loses about 2e-6 relative accuracy on those.
            from OCP.Convert import Convert_ParameterisationType as P
            spans = max(1, min(4, math.ceil((last - first) / (0.5 * math.pi) - 1e-9)))
            kind = {1: P.Convert_TgtThetaOver2_1, 2: P.Convert_TgtThetaOver2_2,
                    3: P.Convert_TgtThetaOver2_3, 4: P.Convert_TgtThetaOver2_4}[spans]
            bs = GeomConvert.CurveToBSplineCurve_s(trimmed, kind)
        else:
            bs = GeomConvert.CurveToBSplineCurve_s(trimmed)
    return bs


def _plane_frame(surface):
    """(origin, axis, refDirection, direct) of a Geom_Plane, as plain tuples.

    The schema plane normal is its axis and its v direction is axis x ref.
    OCCT's normal is X x Y, which is -Z for a left-handed frame, so the
    axis written is the true normal; then axis x ref equals OCCT's Y and
    the surface keeps OCCT's orientation and parameterization.
    """
    pln = surface.Pln()
    ax = pln.Position()
    z = _xyz(ax.Direction())
    if not ax.Direct():
        z = tuple(-c for c in z)
    return _xyz(ax.Location()), z, _xyz(ax.XDirection()), ax.Direct()


def _schema_point(kind, frame, u, v):
    """Evaluate an analytic omniSolid surface S(u, v) (schema.usda formulas)."""
    o, z, x, p = frame["origin"], frame["axis"], frame["ref"], frame
    y = np.cross(z, x)
    if kind == "cone":
        r = p["radius"] + v * math.tan(p["semiAngle"])
        return o + r * math.cos(u) * x + r * math.sin(u) * y + v * z
    if kind == "sphere":
        r = p["radius"]
        return o + r * (math.cos(v) * math.cos(u) * x + math.cos(v) * math.sin(u) * y
                        + math.sin(v) * z)
    if kind == "torus":
        a = p["majorRadius"] + p["minorRadius"] * math.cos(v)
        return o + a * math.cos(u) * x + a * math.sin(u) * y + p["minorRadius"] * math.sin(v) * z
    raise ValueError(kind)


def _analytic_surface(base, face, u0, u1, v0, v1):
    """Map an OCCT cylinder/cone/sphere/torus onto the schema's analytic API.

    Returns ('cone'|'sphere'|'torus', frame dict, (su0, sv0, su1, sv1)) or None
    for other surface types. OCCT and the schema share the angular u. They
    differ in two places, both handled here and then checked numerically:
      * a left-handed OCCT frame (Y = -(Z x X)): flip the axis, which makes
        axis x ref equal OCCT's Y and negates the v parameter;
      * OCCT cone v runs along the slant line; the schema's v is axial, so
        v_schema = v_occt * cos(semiAngle).
    A cylinder is written as a cone with semiAngle 0, as usd-brep's own
    exporter does.
    """
    from OCP.Geom import Geom_CylindricalSurface, Geom_ConicalSurface, \
        Geom_SphericalSurface, Geom_ToroidalSurface
    if isinstance(base, Geom_CylindricalSurface):
        kind, extra, vscale = "cone", {"radius": base.Radius(), "semiAngle": 0.0}, 1.0
    elif isinstance(base, Geom_ConicalSurface):
        a = float(base.SemiAngle())
        kind, extra, vscale = "cone", {"radius": base.RefRadius(), "semiAngle": a}, math.cos(a)
    elif isinstance(base, Geom_SphericalSurface):
        kind, extra, vscale = "sphere", {"radius": base.Radius()}, 1.0
    elif isinstance(base, Geom_ToroidalSurface):
        kind, extra, vscale = "torus", {"majorRadius": base.MajorRadius(),
                                        "minorRadius": base.MinorRadius()}, 1.0
    else:
        return None
    ax3 = base.Position()
    o = np.array(_xyz(ax3.Location()))
    z = np.array(_xyz(ax3.Direction()))
    x = np.array(_xyz(ax3.XDirection()))
    sign = 1.0
    if not ax3.Direct():
        z, sign = -z, -1.0
        if kind == "cone":
            extra["semiAngle"] = -extra["semiAngle"]
    frame = {"origin": o, "axis": z, "ref": x, **{k: float(v) for k, v in extra.items()}}
    # The schema normal of these surfaces always points away from the axis
    # or centre. OCCT's normal (dU x dV) points inward on a left-handed
    # frame, so such faces are written with their sides swapped and their
    # loops reversed (see _solid_tables).
    flipped = sign < 0
    sv = sorted((sign * vscale * v0, sign * vscale * v1))

    def primary(lo, hi):
        # shift an angular interval by whole turns so it starts in [0, 2pi);
        # S(u + 2k pi) == S(u), and SMLib rebuilds generating circles from
        # these intervals (a torus v of [-2pi, 0] fails AdjustSTEPUVDomain)
        k = math.floor(lo / (2 * math.pi) + 1e-12)
        return lo - 2 * math.pi * k, hi - 2 * math.pi * k

    su0, su1 = primary(float(u0), float(u1))
    if kind == "torus":
        sv = list(primary(sv[0], sv[1]))
    rng = (su0, sv[0], su1, sv[1])
    # numeric guard: the mapped schema surface must reproduce OCCT's points
    scale = max(1.0, float(np.abs(o).max()), abs(v0), abs(v1),
                max(abs(float(v)) for v in extra.values()))  # extra holds radii/angle only
    for fu in (0.0, 0.37, 1.0):
        for fv in (0.0, 0.61, 1.0):
            uo, vo = u0 + fu * (u1 - u0), v0 + fv * (v1 - v0)
            want = np.array(_xyz(base.Value(uo, vo)))
            got = _schema_point(kind, frame, uo, sign * vscale * vo)
            if np.linalg.norm(got - want) > 1e-9 * scale:
                raise OmniSolidWriteError(
                    f"analytic {kind} mapping check failed by {np.linalg.norm(got - want):.3g}")
    frame["flipped"] = flipped
    return kind, frame, rng


def _face_surface(face):
    """('plane', origin, axis, ref, direct), ('cone'|'sphere'|'torus', frame, range)
    or ('nurb', _NurbSurf, (u0,u1,v0,v1))."""
    from OCP.BRep import BRep_Tool
    from OCP.BRepTools import BRepTools
    from OCP.Geom import Geom_Plane, Geom_RectangularTrimmedSurface, Geom_BSplineSurface
    from OCP.GeomConvert import GeomConvert
    surf = BRep_Tool.Surface_s(face)
    base = surf
    while isinstance(base, Geom_RectangularTrimmedSurface):
        base = base.BasisSurface()
    if isinstance(base, Geom_Plane):
        return ("plane",) + _plane_frame(base)
    u0, u1, v0, v1 = (float(x) for x in BRepTools.UVBounds_s(face))
    analytic = _analytic_surface(base, face, u0, u1, v0, v1)
    if analytic is not None:
        return analytic
    if isinstance(base, Geom_BSplineSurface):
        bs = base.Copy()
        if bs.IsUPeriodic():
            bs.SetUNotPeriodic()
        if bs.IsVPeriodic():
            bs.SetVNotPeriodic()
    else:
        bs = GeomConvert.SurfaceToBSplineSurface_s(
            Geom_RectangularTrimmedSurface(base, u0, u1, v0, v1))
        if bs.IsUPeriodic():
            bs.SetUNotPeriodic()
        if bs.IsVPeriodic():
            bs.SetVNotPeriodic()
    nu, nv = bs.NbUPoles(), bs.NbVPoles()
    cvs, w = [], []
    for i in range(1, nu + 1):          # U rows
        for j in range(1, nv + 1):      # V columns
            cvs.append(_xyz(bs.Pole(i, j)))
            w.append(float(bs.Weight(i, j)))
    uk, vk = [], []
    for i in range(1, bs.NbUKnots() + 1):
        uk.extend([float(bs.UKnot(i))] * int(bs.UMultiplicity(i)))
    for i in range(1, bs.NbVKnots() + 1):
        vk.extend([float(bs.VKnot(i))] * int(bs.VMultiplicity(i)))
    data = _NurbSurf(cvs, w, uk, vk, int(bs.UDegree()) + 1, int(bs.VDegree()) + 1, nu, nv)
    b = bs.Bounds()
    return ("nurb", data, tuple(float(x) for x in b))


# ---------------------------------------------------------------------------
# topology walk
# ---------------------------------------------------------------------------

@dataclass
class _BrepTables:
    """Flattened arrays for one Brep, before packing into the BrepArray."""
    vertices: list = field(default_factory=list)          # positions
    edges: list = field(default_factory=list)             # (_Nurb, (t0, t1), (vs, ve))
    faces: list = field(default_factory=list)             # dict per face
    regions: list = field(default_factory=list)           # (type, [shell faceuse lists])
    extent: tuple = ((0, 0, 0), (0, 0, 0))
    tol: float = 2e-5


def _ordered_loop(wire, face):
    """[(edge, forward_bool)] in traversal order; degenerate edges dropped."""
    from OCP.BRep import BRep_Tool
    from OCP.BRepTools import BRepTools_WireExplorer
    from OCP.TopAbs import TopAbs_FORWARD
    out = []
    ex = BRepTools_WireExplorer(wire, face)
    while ex.More():
        e = ex.Current()
        if not BRep_Tool.Degenerated_s(e):
            out.append((e, e.Orientation() == TopAbs_FORWARD))
        ex.Next()
    return out


def _brep_tables(solids: list, tol: float) -> _BrepTables:
    """One Brep holding every solid in `solids` (SMLib's lump layout)."""
    from OCP.BRep import BRep_Tool
    from OCP.BRepTools import BRepTools
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib
    from OCP.TopAbs import (TopAbs_FACE, TopAbs_FORWARD, TopAbs_REVERSED,
                            TopAbs_SHELL, TopAbs_WIRE)
    from OCP.TopExp import TopExp, TopExp_Explorer
    from OCP.TopoDS import TopoDS

    T = _BrepTables(tol=tol)
    vidx, eidx, fidx = _ShapeIndex(), _ShapeIndex(), _ShapeIndex()
    per_solid = []      # per solid: per shell: [(face index, material-side token)]
    for solid in solids:
        shells = []
        ex = TopExp_Explorer(solid, TopAbs_SHELL)
        while ex.More():
            shells.append(TopoDS.Shell(ex.Current()))
            ex.Next()
        if not shells:
            raise OmniSolidWriteError("solid has no shell")
        if len(shells) > 1:
            from OCP.BRepClass3d import BRepClass3d
            outer = BRepClass3d.OuterShell_s(TopoDS.Solid(solid))
            if outer.IsNull():
                raise OmniSolidWriteError("cannot identify the outer shell")
            shells.sort(key=lambda s: 0 if s.IsSame(outer) else 1)

        # faces, loops, edgeuses
        shell_faces = []    # per shell: [(face index, material-side token)]
        for sh in shells:
            entries = []
            fx = TopExp_Explorer(sh, TopAbs_FACE)
            while fx.More():
                f_in_solid = TopoDS.Face(fx.Current())
                fi = fidx.index(f_in_solid)
                # material side of this face: OCCT oriented normal points away from material
                if f_in_solid.Orientation() == TopAbs_FORWARD:
                    material = "opposite"
                elif f_in_solid.Orientation() == TopAbs_REVERSED:
                    material = "same"
                else:
                    raise OmniSolidWriteError("face with INTERNAL/EXTERNAL orientation")
                if fi == len(T.faces):
                    fwd = TopoDS.Face(f_in_solid.Oriented(TopAbs_FORWARD))
                    surf = _face_surface(fwd)
                    flipped = bool(surf[0] in _ANALYTIC_API and surf[1].get("flipped"))
                    loops = []
                    outer_w = BRepTools.OuterWire_s(fwd)
                    wires = [outer_w]
                    wx = TopExp_Explorer(fwd, TopAbs_WIRE)
                    while wx.More():
                        w = TopoDS.Wire(wx.Current())
                        if not w.IsSame(outer_w):
                            wires.append(w)
                        wx.Next()
                    for w in wires:
                        uses = []
                        for e, forward in _ordered_loop(w, fwd):
                            ei = eidx.index(e)
                            if ei == len(T.edges):
                                ef = TopoDS.Edge(e.Oriented(TopAbs_FORWARD))
                                bs = _edge_bspline(ef)
                                v0 = TopExp.FirstVertex_s(ef)
                                v1 = TopExp.LastVertex_s(ef)
                                iv = []
                                for v in (v0, v1):
                                    k = vidx.index(v)
                                    if k == len(T.vertices):
                                        T.vertices.append(_xyz(BRep_Tool.Pnt_s(v)))
                                    iv.append(k)
                                T.edges.append((_bspline_curve_data(bs),
                                                (float(bs.FirstParameter()), float(bs.LastParameter())),
                                                tuple(iv), bs))
                            uses.append((ei, "same" if forward else "opposite"))
                        if uses:
                            loops.append(uses)
                    if not loops:
                        raise OmniSolidWriteError("face without a usable loop")
                    if flipped:
                        # CCW about the schema normal = reversed OCCT traversal
                        loops = [[(ei, "opposite" if o == "same" else "same")
                                  for ei, o in reversed(lp)] for lp in loops]
                    T.faces.append({"surface": surf, "loops": loops, "face": fwd,
                                    "flipped": flipped})
                if T.faces[fi]["flipped"]:
                    material = "same" if material == "opposite" else "opposite"
                entries.append((fi, material))
                fx.Next()
            shell_faces.append(entries)
        per_solid.append(shell_faces)

    # Regions, as SMLib lays out separate lumps: the infinite void bounded
    # by every outer shell, one solidRegion per solid (outer shell then its
    # cavities), then one voidRegion per cavity.
    flip = {"same": "opposite", "opposite": "same"}
    regions = [("voidRegion", [[(fi, flip[m]) for fi, m in sf[0]] for sf in per_solid])]
    for sf in per_solid:
        regions.append(("solidRegion", [list(sh) for sh in sf]))
    for sf in per_solid:
        for cav in sf[1:]:
            regions.append(("voidRegion", [[(fi, flip[m]) for fi, m in cav]]))
    T.regions = regions

    box = Bnd_Box()
    for solid in solids:
        BRepBndLib.Add_s(solid, box, False)
    lo, hi = box.CornerMin(), box.CornerMax()
    T.extent = (_xyz(lo), _xyz(hi))
    return T


# ---------------------------------------------------------------------------
# radial edge data
# ---------------------------------------------------------------------------

def _radial(T: _BrepTables):
    """(edgeuse list, nextRadial list, entry list) for one Brep.

    Edgeuse order = faces in order, loops in order, uses in order.
    Entry type: right-hand-rule traversal about the edge tangent t enters
    the face from its top (surface-normal) side iff n . (t x d) < 0, where d
    points from the edge into the face. With loops counter-clockwise about
    n, d = n x tau (tau = traversal direction), which reduces to
    "topEntry iff the use is opposite". Two-use edges are linked to each
    other; anything else is refused.
    """
    uses = []
    for fi, f in enumerate(T.faces):
        for loop in f["loops"]:
            for ei, ori in loop:
                uses.append((ei, ori, fi))
    by_edge: dict = {}
    for k, (ei, _o, _f) in enumerate(uses):
        by_edge.setdefault(ei, []).append(k)
    nxt = [0] * len(uses)
    for ei, ks in by_edge.items():
        if len(ks) != 2:
            raise OmniSolidWriteError(
                f"edge {ei} has {len(ks)} face uses; only manifold edges are supported")
        a, b = ks
        nxt[a], nxt[b] = b, a
    entry = ["topEntry" if o == "opposite" else "bottomEntry" for _e, o, _f in uses]
    return uses, nxt, entry


# ---------------------------------------------------------------------------
# .usda text
# ---------------------------------------------------------------------------

def _f(x: float) -> str:
    r = repr(float(x))
    return r if ("e" in r or "." in r or "n" in r) else r + ".0"


def _arr(values, kind: str) -> str:
    if kind == "token":
        return "[" + ", ".join(f'"{v}"' for v in values) + "]"
    if kind in ("double", "float"):
        return "[" + ", ".join(_f(v) for v in values) + "]"
    if kind == "uint":
        return "[" + ", ".join(str(int(v)) for v in values) + "]"
    if kind == "tuple":
        return "[" + ", ".join("(" + ", ".join(_f(c) for c in v) + ")" for v in values) + "]"
    if kind == "int2":
        return "[" + ", ".join(f"({int(a)}, {int(b)})" for a, b in values) + "]"
    raise ValueError(kind)


def _pack(tables: list) -> dict:
    """Concatenate per-Brep tables into BrepArray attribute arrays."""
    A = {k: [] for k in (
        "brep:intersectTol3d", "brep:extent", "brep:regionCount",
        "region:shellCount", "region:type", "shell:faceuseCount",
        "shell:wireEdgeCount", "shell:pointType", "faceuse:faceIndex",
        "faceuse:orientationType", "face:loopCount", "face:surfaceType",
        "face:trimType", "face:range", "loop:edgeuseCount", "loop:vertexIndex",
        "edgeuse:edgeIndex", "edgeuse:orientationType", "edgeuse:nextRadialEUIndex",
        "edgeuse:thisRadialEntryType", "edge:curveType", "edge:range",
        "edge:vertexIndices", "vertex:pointType",
        "brep:vertexPoint:point:position",
        "brep:edge3dNurb:curve3d:nurb:controlVertices",
        "brep:edge3dNurb:curve3d:nurb:knots", "brep:edge3dNurb:curve3d:nurb:order",
        "brep:edge3dNurb:curve3d:nurb:vertexCount", "brep:edge3dNurb:curve3d:nurb:weights",
        "brep:surface:plane:origin", "brep:surface:plane:axis",
        "brep:surface:plane:refDirection",
        "brep:surface:nurb:controlVertices", "brep:surface:nurb:uVertexCount",
        "brep:surface:nurb:vVertexCount", "brep:surface:nurb:uOrder",
        "brep:surface:nurb:vOrder", "brep:surface:nurb:uKnots",
        "brep:surface:nurb:vKnots", "brep:surface:nurb:weights") + _ANALYTIC_ATTRS}
    v_off = e_off = f_off = eu_off = 0
    for T in tables:
        uses, nxt, entry = _radial(T)
        A["brep:intersectTol3d"].append(T.tol)
        pad = 0.5 * T.tol
        lo, hi = T.extent
        A["brep:extent"] += [tuple(c - pad for c in lo), tuple(c + pad for c in hi)]
        A["brep:regionCount"].append(len(T.regions))
        for rtype, shells in T.regions:
            A["region:shellCount"].append(len(shells))
            A["region:type"].append(rtype)
            for sh in shells:
                A["shell:faceuseCount"].append(len(sh))
                A["shell:wireEdgeCount"].append(0)
                A["shell:pointType"].append("none")
                for fi, side in sh:
                    A["faceuse:faceIndex"].append(fi + f_off)
                    A["faceuse:orientationType"].append(side)
        for f in T.faces:
            A["face:loopCount"].append(len(f["loops"]))
            A["face:trimType"].append("general")
            s = f["surface"]
            if s[0] == "plane":
                _kind, origin, axis, ref, direct = s
                A["face:surfaceType"].append("BrepSurfacePlaneAPI")
                A["brep:surface:plane:origin"].append(origin)
                A["brep:surface:plane:axis"].append(axis)
                A["brep:surface:plane:refDirection"].append(ref)
                A["face:range"] += list(_plane_range(f["face"], origin, axis, ref))
            elif s[0] in _ANALYTIC_API:
                kind, frame, (su0, sv0, su1, sv1) = s
                A["face:surfaceType"].append(_ANALYTIC_API[kind])
                for key in _ANALYTIC_FIELDS[kind]:
                    val = frame[key]
                    if key in ("origin", "axis", "ref"):
                        val = tuple(float(c) for c in val)
                    A[_analytic_attr(kind, key)].append(val)
                A["face:range"] += [(su0, sv0), (su1, sv1)]
            else:
                _kind, d, bounds = s
                A["face:surfaceType"].append("BrepSurfaceNurbAPI")
                A["brep:surface:nurb:controlVertices"] += d.cvs
                A["brep:surface:nurb:weights"] += d.weights
                A["brep:surface:nurb:uVertexCount"].append(d.nu)
                A["brep:surface:nurb:vVertexCount"].append(d.nv)
                A["brep:surface:nurb:uOrder"].append(d.u_order)
                A["brep:surface:nurb:vOrder"].append(d.v_order)
                A["brep:surface:nurb:uKnots"] += d.u_knots
                A["brep:surface:nurb:vKnots"] += d.v_knots
                u0, u1, v0, v1 = bounds
                A["face:range"] += [(u0, v0), (u1, v1)]
            for loop in f["loops"]:
                A["loop:edgeuseCount"].append(len(loop))
                A["loop:vertexIndex"].append(9999999)
        for k, (ei, ori, _fi) in enumerate(uses):
            A["edgeuse:edgeIndex"].append(ei + e_off)
            A["edgeuse:orientationType"].append(ori)
            A["edgeuse:nextRadialEUIndex"].append(nxt[k] + eu_off)
            A["edgeuse:thisRadialEntryType"].append(entry[k])
        for nurb, (t0, t1), (vs, ve), _bs in T.edges:
            A["edge:curveType"].append("BrepCurve3dNurbAPI")
            A["edge:range"] += [t0, t1]
            A["edge:vertexIndices"].append((vs + v_off, ve + v_off))
            A["brep:edge3dNurb:curve3d:nurb:controlVertices"] += nurb.cvs
            A["brep:edge3dNurb:curve3d:nurb:weights"] += nurb.weights
            A["brep:edge3dNurb:curve3d:nurb:knots"] += nurb.knots
            A["brep:edge3dNurb:curve3d:nurb:order"].append(nurb.order)
            A["brep:edge3dNurb:curve3d:nurb:vertexCount"].append(len(nurb.cvs))
        for p in T.vertices:
            A["vertex:pointType"].append("BrepPointAPI")
            A["brep:vertexPoint:point:position"].append(p)
        v_off += len(T.vertices)
        e_off += len(T.edges)
        f_off += len(T.faces)
        eu_off += len(uses)
    return A


def _plane_range(face, origin, axis, ref):
    """UV box of the face in the (ref, axis x ref) frame of the plane."""
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    from OCP.TopAbs import TopAbs_EDGE
    o, n, x = np.array(origin), np.array(axis), np.array(ref)
    y = np.cross(n, x)
    pts = []
    ex = TopExp_Explorer(face, TopAbs_EDGE)
    while ex.More():
        c = BRepAdaptor_Curve(TopoDS.Edge(ex.Current()))
        a, b = c.FirstParameter(), c.LastParameter()
        for k in range(17):
            pts.append(np.array(_xyz(c.Value(a + (b - a) * k / 16.0))))
        ex.Next()
    if not pts:
        raise OmniSolidWriteError("planar face without edges")
    P = np.array(pts) - o
    u, v = P @ x, P @ y
    return (float(u.min()), float(v.min())), (float(u.max()), float(v.max()))


_ANALYTIC_API = {"cone": "BrepSurfaceConeAPI", "sphere": "BrepSurfaceSphereAPI",
                 "torus": "BrepSurfaceTorusAPI"}
_ANALYTIC_FIELDS = {
    "cone": ("origin", "axis", "ref", "radius", "semiAngle"),
    "sphere": ("origin", "axis", "ref", "radius"),
    "torus": ("origin", "axis", "ref", "majorRadius", "minorRadius"),
}


def _analytic_attr(kind: str, key: str) -> str:
    if key == "origin" and kind == "sphere":
        key = "center"
    if key == "ref":
        key = "refDirection"
    return f"brep:surface:{kind}:{key}"


_ANALYTIC_ATTRS = tuple(_analytic_attr(k, f) for k, fs in _ANALYTIC_FIELDS.items() for f in fs)


def _analytic_type(name: str):
    key = name.rsplit(":", 1)[1]
    if key in ("origin", "center"):
        return ("point3d[]", "tuple")
    if key in ("axis", "refDirection"):
        return ("vector3d[]", "tuple")
    return ("double[]", "double")


_TYPES = {
    "brep:intersectTol3d": ("double[]", "double"),
    "brep:extent": ("double3[]", "tuple"),
    "brep:regionCount": ("uint[]", "uint"),
    "region:shellCount": ("uint[]", "uint"),
    "region:type": ("token[]", "token"),
    "shell:faceuseCount": ("uint[]", "uint"),
    "shell:wireEdgeCount": ("uint[]", "uint"),
    "shell:pointType": ("token[]", "token"),
    "faceuse:faceIndex": ("uint[]", "uint"),
    "faceuse:orientationType": ("token[]", "token"),
    "face:loopCount": ("uint[]", "uint"),
    "face:surfaceType": ("token[]", "token"),
    "face:trimType": ("token[]", "token"),
    "face:range": ("double2[]", "tuple"),
    "loop:edgeuseCount": ("uint[]", "uint"),
    "loop:vertexIndex": ("uint[]", "uint"),
    "edgeuse:edgeIndex": ("uint[]", "uint"),
    "edgeuse:orientationType": ("token[]", "token"),
    "edgeuse:nextRadialEUIndex": ("uint[]", "uint"),
    "edgeuse:thisRadialEntryType": ("token[]", "token"),
    "edge:curveType": ("token[]", "token"),
    "edge:range": ("double[]", "double"),
    "edge:vertexIndices": ("int2[]", "int2"),
    "vertex:pointType": ("token[]", "token"),
    "brep:vertexPoint:point:position": ("point3d[]", "tuple"),
    "brep:edge3dNurb:curve3d:nurb:controlVertices": ("point3d[]", "tuple"),
    "brep:edge3dNurb:curve3d:nurb:knots": ("double[]", "double"),
    "brep:edge3dNurb:curve3d:nurb:order": ("uint[]", "uint"),
    "brep:edge3dNurb:curve3d:nurb:vertexCount": ("uint[]", "uint"),
    "brep:edge3dNurb:curve3d:nurb:weights": ("double[]", "double"),
    "brep:surface:plane:origin": ("point3d[]", "tuple"),
    "brep:surface:plane:axis": ("vector3d[]", "tuple"),
    "brep:surface:plane:refDirection": ("vector3d[]", "tuple"),
    "brep:surface:nurb:controlVertices": ("point3d[]", "tuple"),
    "brep:surface:nurb:uVertexCount": ("uint[]", "uint"),
    "brep:surface:nurb:vVertexCount": ("uint[]", "uint"),
    "brep:surface:nurb:uOrder": ("uint[]", "uint"),
    "brep:surface:nurb:vOrder": ("uint[]", "uint"),
    "brep:surface:nurb:uKnots": ("double[]", "double"),
    "brep:surface:nurb:vKnots": ("double[]", "double"),
    "brep:surface:nurb:weights": ("double[]", "double"),
    **{n: _analytic_type(n) for n in _ANALYTIC_ATTRS},
}


def solids_of(shape) -> list:
    from OCP.TopAbs import TopAbs_SOLID
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS
    out = []
    ex = TopExp_Explorer(shape, TopAbs_SOLID)
    while ex.More():
        out.append(TopoDS.Solid(ex.Current()))
        ex.Next()
    return out


def write_omnisolid_usda(shape, path: str, *, prim_name: str = "Brep0",
                         intersect_tol: float | None = None,
                         one_brep_per_solid: bool = False) -> str:
    """Write every solid of `shape` into one BrepArray prim. Returns the prim path.

    By default all solids go into a single Brep (one solidRegion per solid,
    the layout SMLib itself uses for separate lumps), so usd-brep's
    import_brep returns them as one body. one_brep_per_solid=True packs one
    Brep per solid instead; usd-brep then needs import_breps to see them
    all. brepkernel's reader accepts both.
    """
    from OCP.BRep import BRep_Tool
    from OCP.TopAbs import TopAbs_EDGE, TopAbs_VERTEX
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS
    solids = solids_of(shape)
    if not solids:
        raise OmniSolidWriteError("shape contains no solid")
    def tol_of(group):
        if intersect_tol is not None:
            return float(intersect_tol)
        mx = 0.0
        for s in group:
            for kind, cast in ((TopAbs_VERTEX, TopoDS.Vertex), (TopAbs_EDGE, TopoDS.Edge)):
                ex = TopExp_Explorer(s, kind)
                while ex.More():
                    mx = max(mx, float(BRep_Tool.Tolerance_s(cast(ex.Current()))))
                    ex.Next()
        return max(2e-5, 2.0 * mx)      # usd-brep's own exports use 2e-5
    groups = [[s] for s in solids] if one_brep_per_solid else [solids]
    tables = [_brep_tables(g, tol_of(g)) for g in groups]
    A = _pack(tables)
    schemas = ["BrepPointAPI:vertexPoint", "BrepCurve3dNurbAPI:edge3dNurb"]
    if A["brep:surface:plane:origin"]:
        schemas.append("BrepSurfacePlaneAPI")
    if A["brep:surface:nurb:uOrder"]:
        schemas.append("BrepSurfaceNurbAPI")
    for kind, api in _ANALYTIC_API.items():
        if A[_analytic_attr(kind, "axis")]:
            schemas.append(api)
    lines = ["#usda 1.0", "", 'def Xform "World"', "{",
             f'    def BrepArray "{prim_name}" (',
             "        prepend apiSchemas = [" + ", ".join(f'"{s}"' for s in schemas) + "]",
             '        customData = {', '            string source = "brepkernel"', "        }",
             "    )", "    {"]
    for name in sorted(A):
        values = A[name]
        if not values and (name.startswith("brep:surface") or name.startswith("brep:edge3dNurb")):
            continue
        usd_t, kind = _TYPES[name]
        lines.append(f"        uniform {usd_t} {name} = {_arr(values, kind)}")
    ext = A["brep:extent"]
    lo = tuple(min(p[i] for p in ext[0::2]) for i in range(3))
    hi = tuple(max(p[i] for p in ext[1::2]) for i in range(3))
    lines.append(f"        uniform float3[] extent = {_arr([lo, hi], 'tuple')}")
    lines.append("        uniform token[] wireEdge:curveType = []")
    lines.append("        uniform double[] wireEdge:range = []")
    lines.append("        uniform int2[] wireEdge:vertexIndices = []")
    lines += ["    }", "}", ""]
    with open(path, "w") as fh:
        fh.write("\n".join(lines))
    return f"/World/{prim_name}"
