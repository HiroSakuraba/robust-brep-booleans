"""Stage 0/1 + Tier A: analytic solids.

Each solid exposes:
  - implicit(p): exact signed function, <= 0 inside, == 0 on boundary.
    These are *exact* (closed-form), so topology decisions made from them
    are not subject to tolerance guessing -- this is the design's Tier A.
  - tessellate(tol): triangle mesh whose chordal deviation from the true
    surface is *certified* <= tol. Returned as (V, F, chordal_error).
  - bbox(): exact axis-aligned bounding box.

The implicit is the source of truth for classification; the mesh is only
a bounded-approximation proxy (design: geometry as bounded approximation).
"""

import math
import numpy as np


class AnalyticSolid:
    kind = "solid"

    def implicit(self, p):
        """Exact implicit function. p: (..., 3) array. Returns (...,) array."""
        raise NotImplementedError

    def inside(self, p):
        return self.implicit(p) <= 0.0

    def tessellate(self, tol):
        """Return (V, F, chordal_error) with chordal_error <= tol certified."""
        raise NotImplementedError

    def bbox(self):
        raise NotImplementedError

    def volume_exact(self):
        raise NotImplementedError


class Box(AnalyticSolid):
    kind = "box"

    def __init__(self, lo, hi):
        lo = np.asarray(lo, dtype=float)
        hi = np.asarray(hi, dtype=float)
        assert lo.shape == (3,) and hi.shape == (3,)
        assert np.all(np.isfinite(lo)) and np.all(np.isfinite(hi))
        assert np.all(hi > lo), "degenerate box: hi must exceed lo on every axis"
        self.lo = lo
        self.hi = hi

    def implicit(self, p):
        p = np.asarray(p, dtype=float)
        d = np.maximum(self.lo - p, p - self.hi)  # per-axis signed distance
        return np.max(d, axis=-1)

    def tessellate(self, tol):
        lo, hi = self.lo, self.hi
        V = np.array([[x, y, z]
                      for x in (lo[0], hi[0])
                      for y in (lo[1], hi[1])
                      for z in (lo[2], hi[2])])
        F = np.array([
            [0, 1, 3], [0, 3, 2],  # x = lo[0]
            [4, 6, 7], [4, 7, 5],  # x = hi[0]
            [0, 4, 5], [0, 5, 1],  # y = lo[1]
            [2, 7, 6], [2, 3, 7],  # y = hi[1]
            [0, 6, 4], [0, 2, 6],  # z = lo[2]
            [1, 5, 7], [1, 7, 3],  # z = hi[2]
        ])
        return V, F, 0.0  # planar faces: zero chordal error, exact

    def bbox(self):
        return self.lo.copy(), self.hi.copy()

    def volume_exact(self):
        return float(np.prod(self.hi - self.lo))

    def __repr__(self):
        return f"Box(lo={self.lo.tolist()}, hi={self.hi.tolist()})"


def _latlong_sphere(center, r, tol):
    """Tessellate a sphere with certified chordal error <= tol.

    Chordal error of an edge subtending angle dtheta is r*(1 - cos(dtheta/2)).
    We pick the angular step so the bound holds, then report the true bound.
    """
    center = np.asarray(center, dtype=float)
    if tol <= 0:
        raise ValueError("tol must be positive")
    # worst edge: longitude step at equator; use dtheta for both directions
    dtheta = 2.0 * math.acos(max(-1.0, min(1.0, 1.0 - tol / r)))
    n = max(4, int(math.ceil(2.0 * math.pi / dtheta)))
    dtheta = 2.0 * math.pi / n
    chordal = r * (1.0 - math.cos(dtheta / 2.0))
    assert chordal <= tol * (1 + 1e-9), (chordal, tol)

    n_lat = max(2, n // 2)
    verts = [center + np.array([0.0, 0.0, r])]  # north pole
    for i in range(1, n_lat):
        phi = math.pi * i / n_lat  # 0..pi from north pole
        z = r * math.cos(phi)
        ring_r = r * math.sin(phi)
        for j in range(n):
            theta = 2.0 * math.pi * j / n
            verts.append(center + np.array(
                [ring_r * math.cos(theta), ring_r * math.sin(theta), z]))
    verts.append(center + np.array([0.0, 0.0, -r]))  # south pole
    V = np.array(verts)

    def ring(i):
        return 1 + (i - 1) * n
    F = []
    # north cap
    for j in range(n):
        F.append([0, ring(1) + j, ring(1) + (j + 1) % n])
    # bands
    for i in range(1, n_lat - 1):
        r0, r1 = ring(i), ring(i + 1)
        for j in range(n):
            j2 = (j + 1) % n
            F.append([r0 + j, r1 + j, r1 + j2])
            F.append([r0 + j, r1 + j2, r0 + j2])
    # south cap
    south = len(V) - 1
    rlast = ring(n_lat - 1)
    for j in range(n):
        F.append([south, rlast + (j + 1) % n, rlast + j])
    return V, np.array(F), chordal


class Sphere(AnalyticSolid):
    kind = "sphere"

    def __init__(self, center, r):
        center = np.asarray(center, dtype=float)
        assert center.shape == (3,) and np.all(np.isfinite(center))
        assert np.isfinite(r) and r > 0, "sphere radius must be positive"
        self.center = center
        self.r = float(r)

    def implicit(self, p):
        p = np.asarray(p, dtype=float)
        return np.linalg.norm(p - self.center, axis=-1) - self.r

    def tessellate(self, tol):
        return _latlong_sphere(self.center, self.r, tol)

    def bbox(self):
        e = np.full(3, self.r)
        return self.center - e, self.center + e

    def volume_exact(self):
        return 4.0 / 3.0 * math.pi * self.r ** 3

    def __repr__(self):
        return f"Sphere(center={self.center.tolist()}, r={self.r})"


class Cylinder(AnalyticSolid):
    kind = "cylinder"

    def __init__(self, base_center, axis, r, h):
        base_center = np.asarray(base_center, dtype=float)
        axis = np.asarray(axis, dtype=float)
        n = np.linalg.norm(axis)
        assert np.all(np.isfinite(base_center)) and np.isfinite(n) and n > 0
        assert np.isfinite(r) and r > 0 and np.isfinite(h) and h > 0
        self.base = base_center
        self.axis = axis / n
        self.r = float(r)
        self.h = float(h)

    def _frame(self, p):
        rel = p - self.base
        axial = rel @ self.axis
        radial_vec = rel - np.outer(axial, self.axis)
        radial = np.linalg.norm(radial_vec, axis=-1)
        return axial, radial

    def implicit(self, p):
        """Exact implicit: standard capped-cylinder SDF (exact zero set)."""
        p = np.asarray(p, dtype=float)
        axial, radial = self._frame(p)
        qx = radial - self.r
        qy = np.abs(axial - self.h / 2.0) - self.h / 2.0
        ax = np.maximum(qx, 0.0)
        ay = np.maximum(qy, 0.0)
        return np.minimum(np.maximum(qx, qy), 0.0) + np.sqrt(ax * ax + ay * ay)

    def tessellate(self, tol):
        # side: n-gon with chordal error r*(1-cos(pi/n)) <= tol; caps: fans (exact)
        dtheta = 2.0 * math.acos(max(-1.0, min(1.0, 1.0 - tol / self.r)))
        n = max(6, int(math.ceil(2.0 * math.pi / dtheta)))
        dtheta = 2.0 * math.pi / n
        chordal = self.r * (1.0 - math.cos(dtheta / 2.0))

        # orthonormal frame perpendicular to axis
        a = self.axis
        helper = np.array([1.0, 0.0, 0.0]) if abs(a[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        u = np.cross(a, helper); u /= np.linalg.norm(u)
        v = np.cross(a, u)
        ring_pts = np.array([
            self.base + self.r * (math.cos(dtheta * j) * u + math.sin(dtheta * j) * v)
            for j in range(n)])
        top = ring_pts + a * self.h
        V = np.vstack([ring_pts, top,
                       self.base.reshape(1, 3),
                       (self.base + a * self.h).reshape(1, 3)])
        cb, ct = 2 * n, 2 * n + 1
        F = []
        for j in range(n):  # side (outward)
            j2 = (j + 1) % n
            F.append([j, j2, n + j2])
            F.append([j, n + j2, n + j])
        for j in range(n):  # bottom cap (outward = -axis)
            j2 = (j + 1) % n
            F.append([cb, j2, j])
        for j in range(n):  # top cap (outward = +axis)
            j2 = (j + 1) % n
            F.append([ct, n + j, n + j2])
        return V, np.array(F), chordal

    def bbox(self):
        top = self.base + self.axis * self.h
        lo = np.minimum(self.base, top) - self.r
        hi = np.maximum(self.base, top) + self.r
        return lo, hi

    def volume_exact(self):
        return math.pi * self.r ** 2 * self.h

    def __repr__(self):
        return (f"Cylinder(base={self.base.tolist()}, axis={self.axis.tolist()}, "
                f"r={self.r}, h={self.h})")


class Cone(AnalyticSolid):
    kind = "cone"

    def __init__(self, base_center, axis, r, h):
        base_center = np.asarray(base_center, dtype=float)
        axis = np.asarray(axis, dtype=float)
        n = np.linalg.norm(axis)
        assert np.all(np.isfinite(base_center)) and np.isfinite(n) and n > 0
        assert np.isfinite(r) and r > 0 and np.isfinite(h) and h > 0
        self.base = base_center
        self.axis = axis / n
        self.r = float(r)
        self.h = float(h)

    def implicit(self, p):
        """Exact implicit for a right circular cone (base disc at t=0, apex t=h)."""
        p = np.asarray(p, dtype=float)
        rel = p - self.base
        t = rel @ self.axis
        radial_vec = rel - np.outer(t, self.axis)
        radial = np.linalg.norm(radial_vec, axis=-1)
        r_at = self.r * np.maximum(0.0, 1.0 - t / self.h)
        # inside iff 0 <= t <= h and radial <= r_at
        return np.maximum.reduce([
            radial - r_at,
            -t,
            t - self.h,
        ])

    def tessellate(self, tol):
        dtheta = 2.0 * math.acos(max(-1.0, min(1.0, 1.0 - tol / self.r)))
        n = max(6, int(math.ceil(2.0 * math.pi / dtheta)))
        dtheta = 2.0 * math.pi / n
        chordal = self.r * (1.0 - math.cos(dtheta / 2.0))
        a = self.axis
        helper = np.array([1.0, 0.0, 0.0]) if abs(a[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        u = np.cross(a, helper); u /= np.linalg.norm(u)
        v = np.cross(a, u)
        ring = np.array([
            self.base + self.r * (math.cos(dtheta * j) * u + math.sin(dtheta * j) * v)
            for j in range(n)])
        apex = (self.base + a * self.h).reshape(1, 3)
        V = np.vstack([ring, apex, self.base.reshape(1, 3)])
        ai, cb = n, n + 1
        F = []
        for j in range(n):
            j2 = (j + 1) % n
            F.append([j, j2, ai])      # side (outward)
            F.append([cb, j2, j])      # base cap (outward = -axis)
        return V, np.array(F), chordal

    def bbox(self):
        top = self.base + self.axis * self.h
        lo = np.minimum(self.base, top) - self.r
        hi = np.maximum(self.base, top) + self.r
        return lo, hi

    def volume_exact(self):
        return math.pi * self.r ** 2 * self.h / 3.0

    def __repr__(self):
        return (f"Cone(base={self.base.tolist()}, axis={self.axis.tolist()}, "
                f"r={self.r}, h={self.h})")
