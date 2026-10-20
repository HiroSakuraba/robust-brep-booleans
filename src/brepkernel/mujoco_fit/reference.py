"""True outlines: rebuilt from mjfit shape specs, or read from a CAD plate.

Mode A (mjfit specs): circles, stadiums and Bezier-rounded polygons are
rebuilt as true curves with the same placement mjfit applies (rotate, then
center, then scale, then the centring shift). Polygon shapes (rect,
trapezoid, ngon, star, polygon, csv, vertex lists) are their own truth.

Mode C (CAD): a plate is accepted when it is an extrusion (every face is a
cap perpendicular to one axis or a wall parallel to it), with one through
hole or one blind hole, optional islands inside the hole, and a convex
outer outline. Its cross-sections give the true loops, and
:func:`polygonize` turns lines and circle arcs into the polygons mjfit
needs, always keeping the polygon inside the material: the hole and the
free space around islands are never made smaller.
"""
from __future__ import annotations

import math

import numpy as np

from .geom2d import Arc, Curve, QuadBezier, SampledCurve, Segment, polygon_loop, sample_loop, signed_area


class Refusal(Exception):
    def __init__(self, kind: str, message: str):
        super().__init__(message)
        self.kind = kind


# ---------------------------------------------------------------------------
# mode A: mjfit shape specs
# ---------------------------------------------------------------------------


def mjfit_centroid(poly: np.ndarray) -> np.ndarray:
    """Same formula as mjfit's centring (area centroid of the polygon)."""
    following = np.roll(poly, -1, axis=0)
    cross = poly[:, 0] * following[:, 1] - following[:, 0] * poly[:, 1]
    return ((poly + following) * cross[:, None]).sum(axis=0) / (3.0 * cross.sum())


def spec_loop(spec, scale: float, shift: np.ndarray) -> tuple[list[Curve], str]:
    """(true loop, kind) for one mjfit outline spec, after mjfit's placement."""
    from mjfit import shapes

    poly = shapes.from_spec(spec) * scale - shift
    if not isinstance(spec, dict):
        return polygon_loop(poly), "polygon"
    kind = spec.get("type")
    center = np.asarray(spec.get("center", (0.0, 0.0)), dtype=float)
    rot = math.radians(float(spec.get("rotate", 0.0)))
    place = scale * center - shift
    if kind == "circle":
        r = float(spec["radius"])
        return [Arc(np.zeros(2), r, 0.0, 2 * math.pi).mapped(rot, scale, place)], "circle"
    if kind == "stadium":
        w, hgt = float(spec["width"]), float(spec["height"])
        long, short = max(w, hgt), min(w, hgt)
        r, l2 = short / 2, (long - short) / 2
        local: list[Curve] = [Arc((l2, 0.0), r, -math.pi / 2, math.pi)]
        if l2 > 0:
            local.append(Segment((l2, r), (-l2, r)))
        local.append(Arc((-l2, 0.0), r, math.pi / 2, math.pi))
        if l2 > 0:
            local.append(Segment((-l2, -r), (l2, -r)))
        extra = math.pi / 2 if hgt > w else 0.0          # mjfit turns tall stadiums by +90 deg
        return [c.mapped(rot + extra, scale, place) for c in local], "stadium"
    if kind == "bezier":
        return _bezier_loop(spec, scale, place, rot), "bezier"
    return polygon_loop(poly), "polygon"


def _bezier_loop(spec, scale, place, rot) -> list[Curve]:
    original = np.asarray(spec["vertices"], dtype=float)
    a = np.asarray(spec.get("rounding", 0.5), dtype=float)
    control = original if signed_area(original) > 0 else original[::-1]
    alphas = np.full(len(control), float(a)) if a.ndim == 0 else np.asarray(a, dtype=float)
    if signed_area(original) < 0 and a.ndim:
        alphas = alphas[::-1]
    n = len(control)
    ends = []
    for i in range(n):
        v, al = control[i], alphas[i]
        if al == 0:
            ends.append((v, v, None))
        else:
            entry = v + al * (control[i - 1] - v)
            exit_ = v + al * (control[(i + 1) % n] - v)
            ends.append((entry, exit_, QuadBezier(entry, v, exit_)))
    loop: list[Curve] = []
    for i in range(n):
        entry, exit_, curve = ends[i]
        if curve is not None:
            loop.append(curve)
        nxt = ends[(i + 1) % n][0]
        if np.linalg.norm(nxt - exit_) > 1e-15:
            loop.append(Segment(exit_, nxt))
    return [c.mapped(rot, scale, place) for c in loop]


# ---------------------------------------------------------------------------
# mode C: CAD plates
# ---------------------------------------------------------------------------


def _edge_curve(edge) -> Curve:
    """An OCCT edge in a z = const plane as a traversal-ordered curve."""
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    from OCP.GeomAbs import GeomAbs_Circle, GeomAbs_Line
    from OCP.gp import gp_Pnt, gp_Vec
    from OCP.TopAbs import TopAbs_REVERSED

    ad = BRepAdaptor_Curve(edge)
    u0, u1 = ad.FirstParameter(), ad.LastParameter()
    rev = edge.Orientation() == TopAbs_REVERSED

    def xy(u):
        p = ad.Value(u)
        return (p.X(), p.Y())

    a_pt, b_pt = (np.array(xy(u1)), np.array(xy(u0))) if rev else (np.array(xy(u0)), np.array(xy(u1)))
    t = ad.GetType()
    if t == GeomAbs_Line:
        return Segment(a_pt, b_pt)
    if t == GeomAbs_Circle:
        circ = ad.Circle()
        c = circ.Location()
        center = np.array([c.X(), c.Y()])
        sgn = 1.0 if circ.Axis().Direction().Z() > 0 else -1.0
        sweep = (u1 - u0) * sgn * (-1.0 if rev else 1.0)
        a0 = math.atan2(a_pt[1] - center[1], a_pt[0] - center[0])
        return Arc(center, circ.Radius(), a0, sweep)

    def d1(u):
        p, v = gp_Pnt(), gp_Vec()
        ad.D1(u, p, v)
        return (v.X(), v.Y())

    from OCP.GCPnts import GCPnts_AbscissaPoint
    length = GCPnts_AbscissaPoint.Length_s(ad)
    vec = dvec = None
    from OCP.GeomAbs import GeomAbs_BSplineCurve, GeomAbs_Ellipse
    if t == GeomAbs_BSplineCurve:
        vec, dvec = _bspline_eval(ad.BSpline())
    elif t == GeomAbs_Ellipse:
        vec, dvec = _ellipse_eval(ad.Ellipse())
    if vec is not None:
        if rev:
            return SampledCurve(lambda u: tuple(vec(np.array([u0 + u1 - u]))[0]),
                                lambda u: tuple(-dvec(np.array([u0 + u1 - u]))[0]), u0, u1, length,
                                kind=str(t).split(".")[-1], vec=lambda u: vec(u0 + u1 - u),
                                dvec=lambda u: -dvec(u0 + u1 - u))
        return SampledCurve(lambda u: tuple(vec(np.array([u]))[0]), lambda u: tuple(dvec(np.array([u]))[0]),
                            u0, u1, length, kind=str(t).split(".")[-1], vec=vec, dvec=dvec)
    if rev:
        def back(u):
            x, y = d1(u0 + u1 - u)
            return (-x, -y)
        return SampledCurve(lambda u: xy(u0 + u1 - u), back, u0, u1, length, kind=str(t).split(".")[-1])
    return SampledCurve(xy, d1, u0, u1, length, kind=str(t).split(".")[-1])


def _bspline_eval(curve):
    """Vectorized (rational) B-spline in xy from OCCT poles, knots and weights."""
    from scipy.interpolate import BSpline
    deg = curve.Degree()
    poles = np.array([[curve.Pole(i).X(), curve.Pole(i).Y()] for i in range(1, curve.NbPoles() + 1)])
    w = np.array([curve.Weight(i) for i in range(1, curve.NbPoles() + 1)])
    knots = []
    for i in range(1, curve.NbKnots() + 1):
        knots += [curve.Knot(i)] * curve.Multiplicity(i)
    knots = np.array(knots)
    if curve.IsPeriodic():
        raise Refusal("UnsupportedOutlineCurve", "a periodic B-spline outline edge")
    num = BSpline(knots, poles * w[:, None], deg, extrapolate=False)
    den = BSpline(knots, w, deg, extrapolate=False)
    dnum, dden = num.derivative(), den.derivative()
    lo, hi = knots[deg], knots[-deg - 1]

    def vec(u):
        u = np.clip(np.asarray(u, dtype=float), lo, hi)
        return num(u) / den(u)[:, None]

    def dvec(u):
        u = np.clip(np.asarray(u, dtype=float), lo, hi)
        d = den(u)[:, None]
        return (dnum(u) * d - num(u) * dden(u)[:, None]) / (d * d)
    return vec, dvec


def _ellipse_eval(e):
    pos = e.Position()
    c = np.array([pos.Location().X(), pos.Location().Y()])
    xd = np.array([pos.XDirection().X(), pos.XDirection().Y()])
    yd = np.array([pos.YDirection().X(), pos.YDirection().Y()])
    a, b = e.MajorRadius(), e.MinorRadius()

    def vec(u):
        u = np.asarray(u, dtype=float)[:, None]
        return c + a * np.cos(u) * xd + b * np.sin(u) * yd

    def dvec(u):
        u = np.asarray(u, dtype=float)[:, None]
        return -a * np.sin(u) * xd + b * np.cos(u) * yd
    return vec, dvec


def section_loops(shape, z: float, tol: float) -> list[list[Curve]]:
    """Closed loops of the plane section z = const, chained by endpoints."""
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
    from OCP.gp import gp_Dir, gp_Pln, gp_Pnt
    from OCP.TopAbs import TopAbs_EDGE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS

    sec = BRepAlgoAPI_Section(shape, gp_Pln(gp_Pnt(0, 0, z), gp_Dir(0, 0, 1)), False)
    sec.ComputePCurveOn1(False)
    sec.Approximation(False)
    sec.Build()
    if not sec.IsDone():
        raise Refusal("SectionFailed", f"cannot section the part at z = {z:g}")
    curves = []
    ex = TopExp_Explorer(sec.Shape(), TopAbs_EDGE)
    while ex.More():
        curves.append(_edge_curve(TopoDS.Edge(ex.Current())))
        ex.Next()
    loops, used = [], [False] * len(curves)
    for i in range(len(curves)):
        if used[i]:
            continue
        used[i] = True
        loop = [curves[i]]
        while np.linalg.norm(loop[-1].end() - loop[0].start()) > tol:
            tail = loop[-1].end()
            for j in range(len(curves)):
                if used[j]:
                    continue
                c = curves[j]
                if np.linalg.norm(c.start() - tail) <= tol:
                    loop.append(c)
                    used[j] = True
                    break
                if np.linalg.norm(c.end() - tail) <= tol:
                    loop.append(_reversed(c))
                    used[j] = True
                    break
            else:
                raise Refusal("OpenSection", f"the section at z = {z:g} does not close into loops")
        loops.append(loop)
    return loops


def _reversed(c: Curve) -> Curve:
    if isinstance(c, Segment):
        return Segment(c.p1, c.p0)
    if isinstance(c, Arc):
        return Arc(c.center, c.radius, c.a0 + c.sweep, -c.sweep)
    if isinstance(c, QuadBezier):
        return QuadBezier(c.p2, c.p1, c.p0)
    f, d, v, dv = c.fn, c.dfn, c.vec, c.dvec
    a, b = c.u0, c.u1

    def back(u):
        x, y = d(a + b - u)
        return (-x, -y)
    return SampledCurve(lambda u: f(a + b - u), back, a, b, c.len, c.kind,
                        vec=lambda u: v(a + b - np.asarray(u, dtype=float)),
                        dvec=lambda u: -dv(a + b - np.asarray(u, dtype=float)))


def reverse_loop(loop: list[Curve]) -> list[Curve]:
    return [_reversed(c) for c in reversed(loop)]


def _faces(shape):
    from OCP.TopAbs import TopAbs_FACE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        yield TopoDS.Face(ex.Current())
        ex.Next()


def _face_area(face) -> float:
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    g = GProp_GProps()
    BRepGProp.SurfaceProperties_s(face, g)
    return g.Mass()


def _sampled_normals(face):
    """(normals, points) on a 5 x 5 grid inside the face's parameter box."""
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.BRepLProp import BRepLProp_SLProps
    s = BRepAdaptor_Surface(face)
    normals, points = [], []
    props = BRepLProp_SLProps(s, 1, 1e-9)
    for a in np.linspace(0.1, 0.9, 5):
        for b in np.linspace(0.1, 0.9, 5):
            u = s.FirstUParameter() + a * (s.LastUParameter() - s.FirstUParameter())
            v = s.FirstVParameter() + b * (s.LastVParameter() - s.FirstVParameter())
            props.SetParameters(u, v)
            if props.IsNormalDefined():
                n, q = props.Normal(), props.Value()
                normals.append([n.X(), n.Y(), n.Z()])
                points.append([q.X(), q.Y(), q.Z()])
    return np.array(normals), np.array(points)


def _flat(normals, points, tol=1e-9):
    """Normal of a face whose samples are all coplanar (a B-spline plane), else None."""
    if len(normals) < 3:
        return None
    n0 = normals[0]
    if not np.all(np.abs(normals @ n0) > 1 - tol):
        return None
    span = float(np.ptp(points, axis=0).max()) or 1.0
    if float(np.abs((points - points[0]) @ n0).max()) > tol * span:
        return None
    return n0


def extrusion_axis(shape) -> np.ndarray:
    """The axis the part is extruded along, or Refusal NotAnExtrusion."""
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.GeomAbs import GeomAbs_Cylinder, GeomAbs_Plane

    info = []
    for f in _faces(shape):
        s = BRepAdaptor_Surface(f)
        t = s.GetType()
        if t == GeomAbs_Plane:
            d = s.Plane().Axis().Direction()
            info.append(("plane", np.array([d.X(), d.Y(), d.Z()]), _face_area(f)))
        elif t == GeomAbs_Cylinder:
            d = s.Cylinder().Axis().Direction()
            info.append(("cyl", np.array([d.X(), d.Y(), d.Z()]), 0.0))
        else:
            normals, points = _sampled_normals(f)
            flat = _flat(normals, points)
            if flat is not None:                     # a plane stored as a B-spline
                info.append(("plane", flat, _face_area(f)))
            else:
                info.append(("other", normals, 0.0))
    candidates = []
    for kind, d, area in info:
        if kind == "plane" and not any(abs(float(np.dot(d, c))) > 1 - 1e-12 for c in candidates):
            candidates.append(d)
    best, best_area = None, -1.0
    for c in candidates:
        cap_area, ok = 0.0, True
        for kind, d, area in info:
            if kind == "plane":
                dot = abs(float(np.dot(d, c)))
                if dot > 1 - 1e-9:
                    cap_area += area
                elif dot > 1e-9:
                    ok = False
            elif kind == "cyl":
                ok = ok and abs(float(np.dot(d, c))) > 1 - 1e-9
            else:
                ok = ok and len(d) > 0 and bool(np.all(np.abs(d @ c) < 1e-7))
            if not ok:
                break
        if ok and cap_area > best_area:
            best, best_area = c, cap_area
    if best is not None:
        k = int(np.argmax(np.abs(best)))
        if abs(best[k]) > 1 - 1e-12:              # rounding from an earlier frame turn: snap
            best = np.zeros(3)
            best[k] = 1.0
    if best is None:
        raise Refusal("NotAnExtrusion",
                      "the part is not a straight extrusion along any axis (every face must be "
                      "a flat cap across the axis or a wall parallel to it)")
    return best


def to_axis_frame(shape, axis: np.ndarray, flip: bool = False, shift_to_zero: bool = True):
    """Shape turned so `axis` is +z (or -z with flip), lowest point moved to z = 0.

    Returns (shape, m) with m the 3x4 map from the input frame. Without the
    shift, an axis-aligned part maps by a signed permutation, which is exact.
    """
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
    from OCP.gp import gp_Ax3, gp_Dir, gp_Pnt, gp_Trsf, gp_Vec

    z = axis * (-1.0 if flip else 1.0)
    ref = np.array([1.0, 0.0, 0.0]) if abs(z[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    x = ref - np.dot(ref, z) * z
    x /= np.linalg.norm(x)
    t = gp_Trsf()
    t.SetTransformation(gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(*z), gp_Dir(*x)))
    out = BRepBuilderAPI_Transform(shape, t, True).Shape()
    m = np.array([[t.Value(i, j) for j in range(1, 5)] for i in range(1, 4)])
    if not shift_to_zero:
        return out, m
    box = Bnd_Box()
    BRepBndLib.Add_s(out, box)
    zmin = box.CornerMin().Z()
    t2 = gp_Trsf()
    t2.SetTranslation(gp_Vec(0, 0, -zmin))
    out = BRepBuilderAPI_Transform(out, t2, True).Shape()
    m = np.array([[t.Value(i, j) for j in range(1, 5)] for i in range(1, 4)])
    m[2, 3] -= zmin
    return out, m


def cap_levels(shape, tol: float) -> list[float]:
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.GeomAbs import GeomAbs_Plane
    zs = []
    for f in _faces(shape):
        s = BRepAdaptor_Surface(f)
        if s.GetType() == GeomAbs_Plane:
            d = s.Plane().Axis().Direction()
            if abs(d.Z()) > 1 - 1e-9:
                zs.append(s.Plane().Location().Z())
        else:
            normals, points = _sampled_normals(f)
            flat = _flat(normals, points)
            if flat is not None and abs(flat[2]) > 1 - 1e-9:
                zs.append(float(points[:, 2].mean()))
    zs.sort()
    out = []
    for z in zs:
        if not out or z - out[-1] > tol:
            out.append(z)
    return out


def classify_loops(loops: list[list[Curve]], h: float):
    """Split loops into (outer, hole, islands) by nesting depth."""
    import shapely
    polys = [sample_loop(lp, h).points for lp in loops]
    geoms = [shapely.Polygon(p) for p in polys]
    depth = []
    for i, p in enumerate(polys):
        probe = p[0]
        depth.append(sum(1 for j, g in enumerate(geoms) if j != i and g.contains(shapely.Point(*probe))))
    by = {}
    for i, d in enumerate(depth):
        by.setdefault(d, []).append(i)
    if len(by.get(0, [])) != 1:
        raise Refusal("UnsupportedPlate", f"the plate has {len(by.get(0, []))} separate outer outlines; "
                      "mjfit needs one plate")
    if max(by) > 2:
        raise Refusal("UnsupportedPlate", "outlines are nested more than three deep (outer, hole, island)")
    holes = by.get(1, [])
    return by[0][0], holes, by.get(2, [])


def orient(loop: list[Curve], free_on_left_ccw: bool, h: float) -> list[Curve]:
    """Loop traversed CCW (free_on_left_ccw) or CW."""
    area = signed_area(sample_loop(loop, h).points)
    if (area > 0) != free_on_left_ccw:
        return reverse_loop(loop)
    return loop


def _polygonize_curve(c: Curve, max_error: float, h: float) -> list:
    """Material-side polygon for any smooth curve (B-splines, ellipses, Bezier).

    Free space is on the left. Where the curve turns right (toward the
    material) its chords lie in the material; where it turns left, the
    corners where neighboring tangents meet do. The curve is split where it
    changes direction of turn, and each run is covered greedily with the
    longest steps that stay within max_error of the dense samples. The
    certificate then measures the result independently.
    """
    step = min(max_error / 4, h) if h > 0 else max_error / 4
    p, t = c.sample(step)
    p = np.vstack([p, c.end()[None]])
    t = np.vstack([t, c.tangent_at_end()[None]])
    n = len(p)
    if n < 3:
        return [p[0]]
    turn = t[:-1, 0] * t[1:, 1] - t[:-1, 1] * t[1:, 0]      # > 0: turning left (toward free space)
    sign = np.sign(np.where(np.abs(turn) < 1e-14, 0.0, turn))
    # runs of constant turning sign (straight stretches join the run they touch)
    runs, start = [], 0
    cur = 0.0
    for i, sg in enumerate(sign):
        if sg == 0:
            continue
        if cur == 0:
            cur = sg
        elif sg != cur:
            runs.append((start, i, cur))
            start, cur = i, sg
    runs.append((start, n - 1, cur))
    out = []

    def chord_ok(i, j):
        a, b = p[i], p[j]
        d = b - a
        L = np.linalg.norm(d)
        if L == 0:
            return True
        q = p[i:j + 1] - a
        dist = np.abs(q[:, 0] * d[1] - q[:, 1] * d[0]) / L
        return float(dist.max()) <= max_error

    def corner(i, j):
        a, b, ta, tb = p[i], p[j], t[i], t[j]
        den = ta[0] * tb[1] - ta[1] * tb[0]
        if abs(den) < 1e-15:
            return None
        s_ = ((b - a)[0] * tb[1] - (b - a)[1] * tb[0]) / den
        return a + s_ * ta

    def corner_ok(i, j):
        x = corner(i, j)
        if x is None:
            return None
        return x if float(np.linalg.norm(p[i:j + 1] - x, axis=1).min()) <= max_error else None

    for a_i, b_i, sg in runs:
        out.append(p[a_i])
        i = a_i
        while i < b_i:
            # longest step from i that stays within max_error (doubling, then bisection)
            ok = (lambda j: chord_ok(i, j)) if sg <= 0 else (lambda j: corner_ok(i, j) is not None)
            hi = i + 1
            while hi < b_i and ok(min(b_i, i + 2 * (hi - i))):
                hi = min(b_i, i + 2 * (hi - i))
            top = min(b_i, i + 2 * (hi - i))
            while top - hi > 1:
                mid = (hi + top) // 2
                if ok(mid):
                    hi = mid
                else:
                    top = mid
            j = hi
            if sg > 0:
                x = corner_ok(i, j)
                if x is not None:
                    out.append(x)
            if j < b_i:
                out.append(p[j])
            i = j
    return out


def polygonize(loop: list[Curve], max_error: float, h: float) -> np.ndarray:
    """Polygon for a loop oriented with free space on its left.

    Arcs that curve around the free space (center on the left) get tangent
    polygons, so their corners sit in the material; arcs that curve around
    the material get chords, which also lie in the material. Either way the
    polygon never takes free space away, and it deviates by at most
    max_error.
    """
    pts: list[np.ndarray] = []
    for c in loop:
        if isinstance(c, Segment):
            pts.append(c.p0)
        elif isinstance(c, Arc):
            r, sw = c.radius, c.sweep
            if sw > 0:            # CCW about its center: center on the left (free side)
                half = math.acos(r / (r + max_error))
                k = max(math.ceil(abs(sw) / (2 * half)), 6 if abs(sw) > 2 * math.pi - 1e-9 else 1)
                th = sw / (2 * k)
                pts.append(c.start())
                for j in range(k):
                    a = c.a0 + (2 * j + 1) * th
                    pts.append(c.center + (r / math.cos(th)) * np.array([math.cos(a), math.sin(a)]))
            else:                 # chords
                half = math.acos(max(-1.0, 1 - max_error / r))
                k = max(math.ceil(abs(sw) / (2 * half)), 6 if abs(sw) > 2 * math.pi - 1e-9 else 1)
                for j in range(k):
                    a = c.a0 + sw * j / k
                    pts.append(c.center + r * np.array([math.cos(a), math.sin(a)]))
        else:
            pts.extend(_polygonize_curve(c, max_error, h))
    poly = np.array(pts)
    # drop repeated and collinear vertices (a tangent polygon starts on a flat)
    keep = []
    n = len(poly)
    scale = float(np.abs(poly).max()) or 1.0
    for i in range(n):
        a, b, cc = poly[i - 1], poly[i], poly[(i + 1) % n]
        if np.linalg.norm(b - a) <= 1e-15 * scale:
            continue
        cr = (b[0] - a[0]) * (cc[1] - b[1]) - (b[1] - a[1]) * (cc[0] - b[0])
        if abs(cr) <= 1e-13 * scale * scale and np.dot(b - a, cc - b) > 0:
            continue
        keep.append(b)
    return np.array(keep)
