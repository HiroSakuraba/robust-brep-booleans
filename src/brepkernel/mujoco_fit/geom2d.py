"""Planar geometry for checking convex MuJoCo pieces against true outlines.

Two kinds of check live here.

Exact checks run on the piece coordinates as read from the MJCF file
(``repr`` floats, so every value is an exact binary rational) using
:class:`fractions.Fraction`: convex hulls, whether two convex polygons
overlap, and areas. Their answers have no rounding error.

Measured checks compare a union of pieces with a true region bounded by
lines, circle arcs, Bezier curves or any OCCT edge. Both boundaries are
sampled with arc-length spacing at most ``h``; distances are taken to the
samples (k-d tree). Distance to a set is 1-Lipschitz, so a maximum taken
over samples is within half a spacing of the true one, and a distance to
samples is within half a spacing of the distance to the curve. Offset
boundaries keep candidates up to 1.5 spacings past a corner. Every
reported number is therefore within ``2 h`` of the true value; that bound
travels with it as ``resolution``.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from fractions import Fraction

import numpy as np

# ---------------------------------------------------------------------------
# curves and loops
# ---------------------------------------------------------------------------


class Curve:
    """One boundary piece, traversed from start to end."""

    def sample(self, h: float) -> tuple[np.ndarray, np.ndarray]:
        """Points (start included, end excluded) and unit tangents, spacing <= h."""
        raise NotImplementedError

    def start(self) -> np.ndarray:
        raise NotImplementedError

    def end(self) -> np.ndarray:
        raise NotImplementedError

    def tangent_at_start(self) -> np.ndarray:
        raise NotImplementedError

    def tangent_at_end(self) -> np.ndarray:
        raise NotImplementedError

    def mapped(self, rot: float, scale: float, shift) -> "Curve":
        """p -> scale * R(rot) p + shift."""
        raise NotImplementedError


def _rot(angle: float) -> np.ndarray:
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, -s], [s, c]])


def _unit(v: np.ndarray) -> np.ndarray:
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


@dataclass
class Segment(Curve):
    p0: np.ndarray
    p1: np.ndarray

    def __post_init__(self):
        self.p0 = np.asarray(self.p0, dtype=float)
        self.p1 = np.asarray(self.p1, dtype=float)

    def length(self) -> float:
        return float(np.linalg.norm(self.p1 - self.p0))

    def sample(self, h):
        n = max(1, math.ceil(self.length() / h))
        t = np.arange(n)[:, None] / n
        pts = self.p0 + t * (self.p1 - self.p0)
        tan = np.repeat(_unit(self.p1 - self.p0)[None], n, axis=0)
        return pts, tan

    def start(self):
        return self.p0

    def end(self):
        return self.p1

    def tangent_at_start(self):
        return _unit(self.p1 - self.p0)

    tangent_at_end = tangent_at_start

    def mapped(self, rot, scale, shift):
        r = _rot(rot)
        return Segment(scale * (r @ self.p0) + shift, scale * (r @ self.p1) + shift)


@dataclass
class Arc(Curve):
    """Circle arc: center, radius, start angle a0, signed sweep (CCW > 0)."""
    center: np.ndarray
    radius: float
    a0: float
    sweep: float

    def __post_init__(self):
        self.center = np.asarray(self.center, dtype=float)

    def point(self, a):
        return self.center + self.radius * np.stack([np.cos(a), np.sin(a)], axis=-1)

    def sample(self, h):
        n = max(1, math.ceil(abs(self.sweep) * self.radius / h))
        a = self.a0 + self.sweep * np.arange(n) / n
        pts = self.point(a)
        sgn = 1.0 if self.sweep > 0 else -1.0
        tan = sgn * np.stack([-np.sin(a), np.cos(a)], axis=-1)
        return pts, tan

    def start(self):
        return self.point(self.a0)

    def end(self):
        return self.point(self.a0 + self.sweep)

    def _tan(self, a):
        sgn = 1.0 if self.sweep > 0 else -1.0
        return sgn * np.array([-math.sin(a), math.cos(a)])

    def tangent_at_start(self):
        return self._tan(self.a0)

    def tangent_at_end(self):
        return self._tan(self.a0 + self.sweep)

    def mapped(self, rot, scale, shift):
        return Arc(scale * (_rot(rot) @ self.center) + shift, scale * self.radius,
                   self.a0 + rot, self.sweep)


@dataclass
class QuadBezier(Curve):
    p0: np.ndarray
    p1: np.ndarray
    p2: np.ndarray

    def __post_init__(self):
        self.p0, self.p1, self.p2 = (np.asarray(p, dtype=float) for p in (self.p0, self.p1, self.p2))

    def _at(self, t):
        t = np.asarray(t)[:, None]
        return (1 - t) ** 2 * self.p0 + 2 * t * (1 - t) * self.p1 + t ** 2 * self.p2

    def _d(self, t):
        t = np.asarray(t)[:, None]
        return 2 * (1 - t) * (self.p1 - self.p0) + 2 * t * (self.p2 - self.p1)

    def sample(self, h):
        # speed <= 2 * max leg, so this parameter step moves at most h
        legs = max(np.linalg.norm(self.p1 - self.p0), np.linalg.norm(self.p2 - self.p1))
        n = max(1, math.ceil(2 * legs / h))
        t = np.arange(n) / n
        return self._at(t), _unit(self._d(t))

    def start(self):
        return self.p0

    def end(self):
        return self.p2

    def tangent_at_start(self):
        return _unit(self.p1 - self.p0)

    def tangent_at_end(self):
        return _unit(self.p2 - self.p1)

    def mapped(self, rot, scale, shift):
        r = _rot(rot)
        return QuadBezier(*(scale * (r @ p) + shift for p in (self.p0, self.p1, self.p2)))


class SampledCurve(Curve):
    """Any curve given as a function of a parameter (used for OCCT edges)."""

    def __init__(self, fn, dfn, u0, u1, length, kind="curve"):
        self.fn, self.dfn, self.u0, self.u1, self.len, self.kind = fn, dfn, u0, u1, length, kind

    def sample(self, h):
        # uniform in arc length via a dense parameter table
        m = max(64, math.ceil(4 * self.len / h))
        u = np.linspace(self.u0, self.u1, m + 1)
        pts = np.array([self.fn(x) for x in u])
        s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])
        n = max(1, math.ceil(s[-1] / h))
        target = s[-1] * np.arange(n) / n
        uu = np.interp(target, s, u)
        p = np.array([self.fn(x) for x in uu])
        t = _unit(np.array([self.dfn(x) for x in uu]))
        return p, t

    def start(self):
        return np.asarray(self.fn(self.u0))

    def end(self):
        return np.asarray(self.fn(self.u1))

    def tangent_at_start(self):
        return _unit(np.asarray(self.dfn(self.u0)))

    def tangent_at_end(self):
        return _unit(np.asarray(self.dfn(self.u1)))


def polygon_loop(poly) -> list[Curve]:
    poly = np.asarray(poly, dtype=float)
    return [Segment(poly[i], poly[(i + 1) % len(poly)]) for i in range(len(poly))]


@dataclass
class LoopSamples:
    points: np.ndarray            # (n, 2), closed loop order
    tangents: np.ndarray          # (n, 2)
    corners: list                 # (vertex, tangent_in, tangent_out)
    probe: int = 0                # a sample in the middle of the longest curve


def sample_loop(loop: list[Curve], h: float) -> LoopSamples:
    pts, tans, corners = [], [], []
    best, probe, start = -1, 0, 0
    for i, c in enumerate(loop):
        p, t = c.sample(h)
        if len(p) > best:
            best, probe = len(p), start + len(p) // 2
        start += len(p)
        pts.append(p)
        tans.append(t)
        nxt = loop[(i + 1) % len(loop)]
        t_in, t_out = c.tangent_at_end(), nxt.tangent_at_start()
        if float(np.dot(t_in, t_out)) < 1 - 1e-12:
            corners.append((np.asarray(nxt.start()), t_in, t_out))
    return LoopSamples(np.vstack(pts), np.vstack(tans), corners, probe)


def signed_area(poly: np.ndarray) -> float:
    x, y = poly[:, 0], poly[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))


# ---------------------------------------------------------------------------
# regions
# ---------------------------------------------------------------------------


class Region:
    """A closed planar region known through membership and boundary samples."""
    h: float

    def contains(self, pts: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    def boundary(self) -> np.ndarray:
        raise NotImplementedError

    def tree(self):
        if getattr(self, "_tree", None) is None:
            from scipy.spatial import cKDTree
            self._tree = cKDTree(self.boundary())
        return self._tree

    def dist_to_boundary(self, pts: np.ndarray, upper: float = np.inf) -> np.ndarray:
        """Distance to the boundary samples; inf where it exceeds `upper`.

        Far queries against a dense curve are slow (the search ball meets many
        leaves), so callers that only need "within d or not" pass `upper`.
        """
        if len(pts) == 0:
            return np.zeros(0)
        d, _ = self.tree().query(pts, distance_upper_bound=upper)
        return d

    def bbox(self):
        b = self.boundary()
        return b.min(axis=0), b.max(axis=0)


class LoopRegion(Region):
    """Even-odd region of closed curve loops (outer, holes, islands...)."""

    def __init__(self, loops: list[list[Curve]], h: float):
        import shapely
        self.loops, self.h = loops, h
        self.samples = [sample_loop(lp, h) for lp in loops]
        geom = None
        for s in self.samples:
            # membership only needs the outline within h/10, not every sample
            ring = shapely.simplify(shapely.LinearRing(s.points), h / 10, preserve_topology=True)
            poly = shapely.Polygon(ring)
            if not poly.is_valid:
                poly = shapely.make_valid(poly)
            geom = poly if geom is None else shapely.symmetric_difference(geom, poly)
        self.geom = geom
        shapely.prepare(self.geom)
        self._b = np.vstack([s.points for s in self.samples])

    def contains(self, pts):
        import shapely
        return shapely.contains_xy(self.geom, pts[:, 0], pts[:, 1])

    def boundary(self):
        return self._b

    def sides(self):
        """Per loop: +1 if the region is on the left of the traversal, else -1."""
        out = []
        for smp in self.samples:
            k = smp.probe
            left = np.array([-smp.tangents[k, 1], smp.tangents[k, 0]])
            pair = np.array([smp.points[k] + 4 * self.h * left, smp.points[k] - 4 * self.h * left])
            a, b = self.contains(pair)
            if a == b:
                raise ValueError("cannot tell which side of an outline is inside")
            out.append(1.0 if a else -1.0)
        return out


class OffsetRegion(Region):
    """Erosion (shrink) or dilation (grow) of a LoopRegion by distance d.

    The boundary is built from the base boundary moved d along its normals
    plus circle arcs of radius d at the corners that need them; candidates
    closer than d - h to the base boundary are dropped (they belong to a
    different part of the offset).
    """

    def __init__(self, base: LoopRegion, d: float, mode: str):
        assert mode in ("erode", "dilate")
        self.base, self.d, self.mode, self.h = base, d, mode, base.h
        h = self.h
        sign = 1.0 if mode == "erode" else -1.0      # move into or out of the base
        cands = []
        for smp, side in zip(base.samples, base.sides()):
            left = np.column_stack([-smp.tangents[:, 1], smp.tangents[:, 0]])
            cands.append(smp.points + (sign * side * d) * left)
            # circle arcs of radius d around corners, the short way between the two
            # moved normals; at corners where the offset edges cross instead, the arc
            # points are closer than d to an edge and the filter below drops them
            for vtx, t_in, t_out in smp.corners:
                m_in = sign * side * np.array([-t_in[1], t_in[0]])
                m_out = sign * side * np.array([-t_out[1], t_out[0]])
                a0 = math.atan2(m_in[1], m_in[0])
                da = (math.atan2(m_out[1], m_out[0]) - a0 + math.pi) % (2 * math.pi) - math.pi
                m = max(2, math.ceil(abs(da) * d / h) + 1)
                a = a0 + da * np.linspace(0, 1, m)
                cands.append(vtx + d * np.column_stack([np.cos(a), np.sin(a)]))
        cand = np.vstack(cands)
        dist = base.dist_to_boundary(cand, upper=2 * d + 4 * h)
        keep = dist >= d - 1.5 * h
        if mode == "erode":
            keep &= base.contains(cand)
        else:
            keep &= ~base.contains(cand)
        self._b = cand[keep]
        if len(self._b) == 0:
            raise ValueError(f"offset by {d:g} leaves nothing")

    def contains(self, pts):
        inside = self.base.contains(pts)
        dist = self.base.dist_to_boundary(pts, upper=self.d * (1 + 1e-9) + self.h)
        if self.mode == "erode":
            return inside & (dist >= self.d)
        return inside | (dist <= self.d)

    def boundary(self):
        return self._b


class DifferenceRegion(Region):
    """A minus B (B's boundary inside A becomes boundary of the difference)."""

    def __init__(self, a: Region, b: Region):
        self.a, self.b, self.h = a, b, max(a.h, b.h)
        ba, bb = a.boundary(), b.boundary()
        keep_a = ~b.contains(ba) | (b.dist_to_boundary(ba, upper=2 * self.h) <= self.h)
        keep_b = a.contains(bb) | (a.dist_to_boundary(bb, upper=2 * self.h) <= self.h)
        self._b = np.vstack([ba[keep_a], bb[keep_b]])

    def contains(self, pts):
        return self.a.contains(pts) & ~self.b.contains(pts)

    def boundary(self):
        return self._b


class PieceRegion(Region):
    """Union of convex polygons (the MuJoCo pieces of one slab)."""

    def __init__(self, polys: list[np.ndarray], h: float, eps: float):
        self.polys = [np.asarray(p, dtype=float) for p in polys]
        self.h, self.eps = h, eps
        pts = []
        for i, p in enumerate(self.polys):
            for k in range(len(p)):
                seg = Segment(p[k], p[(k + 1) % len(p)])
                q, _ = seg.sample(h)
                pts.append(np.vstack([q, seg.p1[None]]))
        pts = np.vstack(pts)
        # an edge sample covered by another piece (closed, within eps) is interior
        owner_count = self._count(pts, closed=True)
        self._b = pts[owner_count <= 1]

    def _count(self, pts, closed):
        count = np.zeros(len(pts), dtype=int)
        for p in self.polys:
            e = np.roll(p, -1, axis=0) - p
            ok = np.ones(len(pts), dtype=bool)
            for k in range(len(p)):
                cr = e[k, 0] * (pts[:, 1] - p[k, 1]) - e[k, 1] * (pts[:, 0] - p[k, 0])
                lim = -self.eps * np.linalg.norm(e[k]) if closed else self.eps * np.linalg.norm(e[k])
                ok &= cr >= lim
            count += ok
        return count

    def contains(self, pts):
        return self._count(pts, closed=True) > 0

    def boundary(self):
        return self._b


# ---------------------------------------------------------------------------
# comparisons
# ---------------------------------------------------------------------------


@dataclass
class Deviation:
    outside: float                    # sim sticks out of the reference by this much
    outside_at: list | None
    inside: float                     # reference sticks out of sim by this much
    inside_at: list | None
    resolution: float                 # each number is within this of the true value
    grid_mismatch_far: int            # grid points in one region only, far from both boundaries
    grid_spacing: float
    extra: dict = field(default_factory=dict)

    def as_dict(self):
        return {"sim_outside_reference": self.outside, "sim_outside_reference_at": self.outside_at,
                "reference_outside_sim": self.inside, "reference_outside_sim_at": self.inside_at,
                "resolution": self.resolution, "grid_mismatch_far": self.grid_mismatch_far,
                "grid_spacing": self.grid_spacing, **self.extra}


def deviation(sim: Region, ref: Region, grid: int = 300) -> Deviation:
    h = max(sim.h, ref.h)
    bs, br = sim.boundary(), ref.boundary()
    out_mask = ~ref.contains(bs)
    d_out = np.zeros(len(bs))
    d_out[out_mask] = ref.dist_to_boundary(bs[out_mask])
    in_mask = ~sim.contains(br)
    d_in = np.zeros(len(br))
    d_in[in_mask] = sim.dist_to_boundary(br[in_mask])
    io, ii = int(np.argmax(d_out)), int(np.argmax(d_in))
    outside, inside = float(d_out[io]), float(d_in[ii])
    # every simulated edge should run close to a real one; an edge far from all of
    # them (a notch or slit cut into the material) shows up here and nowhere else
    budget = 2 * max(outside, inside) + 4 * h
    far_edge = ref.dist_to_boundary(bs, upper=4 * budget)
    ie = int(np.argmax(far_edge))
    edge_far = float(far_edge[ie])
    # region-level check: a mismatch far from both boundaries means the regions
    # differ by a whole component (e.g. a covered hole), which boundary samples miss
    lo = np.minimum(bs.min(axis=0), br.min(axis=0))
    hi = np.maximum(bs.max(axis=0), br.max(axis=0))
    g = float((hi - lo).max()) / grid
    xs = np.arange(lo[0] + g / 2, hi[0], g)
    ys = np.arange(lo[1] + g / 2, hi[1], g)
    gx, gy = np.meshgrid(xs, ys)
    pts = np.column_stack([gx.ravel(), gy.ravel()])
    mism = sim.contains(pts) != ref.contains(pts)
    far = 0
    if mism.any():
        q = pts[mism]
        lim = max(outside, inside) + 2 * h + g
        near = np.minimum(sim.dist_to_boundary(q, upper=2 * lim), ref.dist_to_boundary(q, upper=2 * lim))
        far = int((near > lim).sum())
    return Deviation(outside, bs[io].tolist() if outside > 0 else None, inside,
                     br[ii].tolist() if inside > 0 else None, 2 * h, far, g,
                     {"sim_edge_to_reference_max": edge_far, "sim_edge_budget": budget,
                      "sim_edge_far_at": bs[ie].tolist() if edge_far > budget else None})


def min_distance_between(polys_a: list[np.ndarray], polys_b: list[np.ndarray]) -> tuple[float, list]:
    """Smallest distance between two sets of polygons' edges (float, ~1e-18 error)."""
    def edges(polys):
        a, b = [], []
        for p in polys:
            a.append(p)
            b.append(np.roll(p, -1, axis=0))
        return np.vstack(a), np.vstack(b)

    a0, a1 = edges(polys_a)
    b0, b1 = edges(polys_b)

    def pt_seg(p, s0, s1):
        d = s1 - s0
        t = np.clip(np.einsum("...j,...j", p - s0, d) / np.einsum("...j,...j", d, d), 0, 1)
        return np.linalg.norm(p - (s0 + t[..., None] * d), axis=-1)

    best, where = math.inf, None
    for i in range(0, len(a0), 256):
        A0, A1 = a0[i:i + 256, None], a1[i:i + 256, None]
        B0, B1 = b0[None], b1[None]
        cand = np.minimum.reduce([pt_seg(A0, B0, B1), pt_seg(A1, B0, B1),
                                  pt_seg(B0, A0, A1), pt_seg(B1, A0, A1)])
        k = np.unravel_index(int(np.argmin(cand)), cand.shape)
        if cand[k] < best:
            best = float(cand[k])
            where = ((a0[i + k[0]] + a1[i + k[0]]) / 2).tolist()
    return best, where


# ---------------------------------------------------------------------------
# exact checks
# ---------------------------------------------------------------------------


def frac_points(pts) -> list[tuple[Fraction, Fraction]]:
    return [(Fraction(float(x)), Fraction(float(y))) for x, y in pts]


def exact_hull(pts: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    """Strict convex hull, CCW (collinear points dropped), exact."""
    p = sorted(set(pts))
    if len(p) < 3:
        return p

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower, upper = [], []
    for q in p:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], q) <= 0:
            lower.pop()
        lower.append(q)
    for q in reversed(p):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], q) <= 0:
            upper.pop()
        upper.append(q)
    return lower[:-1] + upper[:-1]


def exact_area(hull) -> Fraction:
    s = Fraction(0)
    for i in range(len(hull)):
        x0, y0 = hull[i]
        x1, y1 = hull[(i + 1) % len(hull)]
        s += x0 * y1 - x1 * y0
    return s / 2


def exact_overlap(a, b) -> bool:
    """Whether two CCW convex polygons share interior points (exact SAT).

    They do not when some edge line of either has every vertex of the other
    on or outside it (touching along an edge or at a point is not overlap).
    """
    for poly, other in ((a, b), (b, a)):
        n = len(poly)
        for i in range(n):
            x0, y0 = poly[i]
            x1, y1 = poly[(i + 1) % n]
            ex, ey = x1 - x0, y1 - y0
            if all(ex * (y - y0) - ey * (x - x0) <= 0 for x, y in other):
                return False
    return True


def loop_signature(loops) -> tuple[int, int]:
    """(outer boundaries, holes) = (CCW loops, CW loops) of exact union loops."""
    ccw = sum(1 for lp in loops if exact_area(lp) > 0)
    return ccw, len(loops) - ccw


def exact_union_loops(hulls: list[list[tuple[Fraction, Fraction]]]) -> list[list[tuple[Fraction, Fraction]]]:
    """Boundary loops of a union of interior-disjoint CCW convex polygons, exactly.

    Every edge is split at the endpoints of the collinear edges that overlap
    it. A sub-edge covered once from each side lies between two pieces and
    cancels; the rest is the union's boundary, chained into loops (CCW for
    outer boundaries, CW for holes). A gap between pieces, however thin,
    shows up as an extra loop.
    """
    def line_key(p, q):
        a, b = q[1] - p[1], p[0] - q[0]           # a x + b y = c
        c = a * p[0] + b * p[1]
        g = a if a != 0 else b                    # signed: both directions share a key
        return (a / g, b / g, c / g)

    edges = []
    for hl in hulls:
        n = len(hl)
        for i in range(n):
            edges.append((hl[i], hl[(i + 1) % n]))
    groups: dict = {}
    for p, q in edges:
        groups.setdefault(line_key(p, q), []).append((p, q))

    boundary = []
    for key, es in groups.items():
        a, b, _ = key
        d = (-b, a)                                  # direction along the line

        def t(pt):
            return pt[0] * d[0] + pt[1] * d[1]

        cuts = sorted({t(x) for e in es for x in e})
        pts_at = {}
        for e in es:
            for x in e:
                pts_at[t(x)] = x
        cover = {}
        for p, q in es:
            tp, tq = t(p), t(q)
            sgn = 1 if tq > tp else -1
            lo, hi = min(tp, tq), max(tp, tq)
            for k in range(len(cuts) - 1):
                if cuts[k] >= lo and cuts[k + 1] <= hi:
                    cover[k] = cover.get(k, 0) + sgn
        for k, v in cover.items():
            if v == 0:
                continue
            u, w = pts_at[cuts[k]], pts_at[cuts[k + 1]]
            for _ in range(abs(v)):
                boundary.append((u, w) if v > 0 else (w, u))

    out_edges: dict = {}
    for u, w in boundary:
        out_edges.setdefault(u, []).append(w)
    loops = []
    while out_edges:
        start = next(iter(out_edges))
        loop, cur = [start], start
        while True:
            nxt = out_edges[cur].pop()
            if not out_edges[cur]:
                del out_edges[cur]
            if nxt == start:
                break
            loop.append(nxt)
            cur = nxt
            if cur not in out_edges:
                raise ValueError("union boundary does not close")
        loops.append(loop)
    return loops
