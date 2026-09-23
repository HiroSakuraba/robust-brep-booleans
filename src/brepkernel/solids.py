"""Stage 0/1 + Tier A: analytic solids.

Each solid exposes:
  - implicit(p): exact signed function, <= 0 inside, == 0 on boundary.
    These are *exact* (closed-form), so topology decisions made from them
    are not subject to tolerance guessing -- this is the design's Tier A.
  - tessellate(tol): triangle mesh whose chordal deviation from the true
    surface is *certified* <= tol. The certificate is computed from the
    mesh actually produced (lattice sampling + Lipschitz bound), never
    from the formula that chose the mesh density.
  - bbox(): exact axis-aligned bounding box.

The implicit is the source of truth for classification; the mesh is only
a bounded-approximation proxy (design: geometry as bounded approximation).
"""

import math
import numpy as np


def _max_pair_angle(a, b, c):
    """Max pairwise angle between corresponding vectors (N,D); zero-safe."""
    def ang(x, y):
        nx = np.linalg.norm(x, axis=1)
        ny = np.linalg.norm(y, axis=1)
        denom = np.maximum(nx * ny, 1e-300)
        cosv = np.clip(np.einsum("ij,ij->i", x, y) / denom, -1.0, 1.0)
        # pairs involving a zero vector contribute 0
        return np.where((nx < 1e-300) | (ny < 1e-300), 0.0, np.arccos(cosv))
    return np.maximum.reduce([ang(a, b), ang(b, c), ang(a, c)])


def _perp_basis(axis):
    a = np.asarray(axis, dtype=float)
    helper = np.array([1.0, 0.0, 0.0]) if abs(a[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(a, helper); u /= np.linalg.norm(u)
    v = np.cross(a, u)
    return u, v


def analytic_certificate(V, F, solid):
    """Exact, rigorous chordal-error bound from the mesh actually produced.

    Wedge lemma: if a triangle's vertices lie on the surface within an
    angular spread dtheta of each other (as seen from the center/axis),
    every interior point p satisfies |p - c| >= R*cos(dtheta/2), so the
    chordal deviation is <= R*(1 - cos(dtheta/2)). Vertices are on the
    surface by construction; triangles are chords of convex solids, so
    interior points lie inside and the deviation is one-sided.
    Planar parts (box faces, caps): max |implicit| at vertices.
    Note |implicit| overestimates true distance for the cone (|grad|>1),
    which is the safe direction for an upper bound.
    """
    V = np.asarray(V, dtype=float)
    F = np.asarray(F)
    vert_dev = np.max(np.abs(solid.implicit(V)))
    n = len(F)
    tri_bound = np.zeros(n)

    kind = solid.kind
    if kind == "sphere":
        # Rigorous: |sum λw_i|^2 = r^2 Σλ_iλ_j cosθ_ij
        #   >= r^2 (1 - 2(1-cosΔθ) Σ_{i<j}λ_iλ_j) >= r^2 (1 - (2/3)(1-cosΔθ)).
        w0 = V[F[:, 0]] - solid.center
        w1 = V[F[:, 1]] - solid.center
        w2 = V[F[:, 2]] - solid.center
        dtheta = _max_pair_angle(w0, w1, w2)
        tri_bound = solid.r * (1.0 - np.sqrt(np.maximum(
            0.0, 1.0 - (2.0 / 3.0) * (1.0 - np.cos(dtheta)))))
    elif kind in ("cylinder", "cone"):
        # Rigorous: |Σλs_i|^2 = Σλ_iλ_j |s_i||s_j| cosθ_ij
        #   >= cosΔθ (Σλ|s_i|)^2, so the Jensen gap
        #   Σλ|s_i| - |Σλs_i| <= (1 - sqrt(cosΔθ)) R_max.
        u, v = _perp_basis(solid.axis)
        rel = V - solid.base
        axial = rel @ solid.axis
        # cap detection eps scales with coordinate magnitude (item 4):
        # axial = (V - base) @ axis rounds at ~eps64 * max|V|.
        coord_scale = max(1.0, float(np.max(np.abs(V)))) if n else 1.0
        cap_eps = 1e-12 * coord_scale
        q = np.stack([rel @ u, rel @ v], axis=1)
        q0, q1, q2 = q[F[:, 0]], q[F[:, 1]], q[F[:, 2]]
        dtheta = _max_pair_angle(q0, q1, q2)
        gap = 1.0 - np.sqrt(np.maximum(np.cos(dtheta), 0.0))
        a0, a1, a2 = axial[F[:, 0]], axial[F[:, 1]], axial[F[:, 2]]
        if kind == "cylinder":
            cap = ((np.abs(a0) < cap_eps) & (np.abs(a1) < cap_eps) & (np.abs(a2) < cap_eps)) | \
                  ((np.abs(a0 - solid.h) < cap_eps) & (np.abs(a1 - solid.h) < cap_eps) & (np.abs(a2 - solid.h) < cap_eps))
            tri_bound = np.where(~cap & (solid.h > 4 * cap_eps), solid.r * gap, 0.0)
        else:
            t_min = np.minimum.reduce([a0, a1, a2])
            cap = (np.abs(a0) < cap_eps) & (np.abs(a1) < cap_eps) & (np.abs(a2) < cap_eps)
            r_at = solid.r * (1.0 - t_min / solid.h)
            tri_bound = np.where(~cap & (solid.h > 4 * cap_eps), r_at * gap, 0.0)
    elif kind == "box":
        pass  # planar: vertex deviation only
    else:
        raise ValueError(f"no certificate rule for {kind}")
    return float(max(vert_dev, np.max(tri_bound) if n else 0.0))


def _certified_tessellate(solid, tol, build):
    """Build a mesh, certify it with the exact bound; refine until <= tol."""
    k = 0
    while True:
        V, F = build(tol, k)
        bound = analytic_certificate(np.asarray(V, float), np.asarray(F), solid)
        if bound <= tol:
            return np.asarray(V, dtype=float), np.asarray(F), float(bound)
        k += 1
        if k > 10:
            raise RuntimeError(
                f"could not certify tessellation to tol={tol} (bound={bound})")


def _angular_density(tol, r, safety=0.5):
    """Segment angle whose certified bound is about tol*safety.

    Certified bounds behave like r*dtheta^2/6 (sphere) and r*dtheta^2/4
    (cylinder/cone) for small angles; the refine loop guarantees the rest.
    """
    return math.sqrt(4.0 * tol * safety / r)


class AnalyticSolid:
    kind = "solid"
    lipschitz = 1.0  # Lipschitz constant of implicit(); used for certificates

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

    def params_equal(self, other):
        """Bitwise parameter equality (Tier A identical-input fast path)."""
        return False


class Box(AnalyticSolid):
    kind = "box"

    def __init__(self, lo, hi):
        lo = np.asarray(lo, dtype=float)
        hi = np.asarray(hi, dtype=float)
        if lo.shape != (3,) or hi.shape != (3,):
            raise ValueError("box corners must be 3-vectors")
        if not (np.all(np.isfinite(lo)) and np.all(np.isfinite(hi))):
            raise ValueError("box corners must be finite")
        if not np.all(hi > lo):
            raise ValueError("degenerate box: hi must exceed lo on every axis")
        self.lo = lo
        self.hi = hi

    def implicit(self, p):
        p = np.asarray(p, dtype=float)
        d = np.maximum(self.lo - p, p - self.hi)  # per-axis signed distance
        return np.max(d, axis=-1)

    def tessellate(self, tol):
        if not tol > 0:
            raise ValueError("tol must be positive")
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
        bound = analytic_certificate(V, F, self)
        if bound > tol:
            raise RuntimeError(
                f"box tessellation certificate {bound} exceeds tol {tol}")
        return V, F, bound

    def bbox(self):
        return self.lo.copy(), self.hi.copy()

    def volume_exact(self):
        return float(np.prod(self.hi - self.lo))

    def params_equal(self, other):
        return (isinstance(other, Box)
                and np.array_equal(self.lo, other.lo)
                and np.array_equal(self.hi, other.hi))

    def __repr__(self):
        return f"Box(lo={self.lo.tolist()}, hi={self.hi.tolist()})"


def _sphere_lattice_mesh(center, r, tol, k):
    """Lat-long sphere mesh; density from tol, doubles with k."""
    dtheta = _angular_density(tol, r)
    n = max(12, int(math.ceil(2.0 * math.pi / dtheta)) * (2 ** k))
    n_lat = max(6, n // 2)
    verts = [center + np.array([0.0, 0.0, r])]
    for i in range(1, n_lat):
        phi = math.pi * i / n_lat
        z = r * math.cos(phi)
        ring_r = r * math.sin(phi)
        th = np.arange(n) * (2.0 * math.pi / n)
        ring = np.stack([ring_r * np.cos(th), ring_r * np.sin(th),
                         np.full(n, z)], axis=1) + center
        verts.append(ring)
    verts.append(center + np.array([0.0, 0.0, -r]))
    V = np.vstack([v.reshape(-1, 3) if v.ndim > 1 else v.reshape(1, 3)
                   for v in verts])

    def ring(i):
        return 1 + (i - 1) * n
    F = []
    for j in range(n):
        F.append([0, ring(1) + j, ring(1) + (j + 1) % n])
    for i in range(1, n_lat - 1):
        r0, r1 = ring(i), ring(i + 1)
        for j in range(n):
            j2 = (j + 1) % n
            F.append([r0 + j, r1 + j, r1 + j2])
            F.append([r0 + j, r1 + j2, r0 + j2])
    south = len(V) - 1
    rlast = ring(n_lat - 1)
    for j in range(n):
        F.append([south, rlast + (j + 1) % n, rlast + j])
    return V, np.array(F)


class Sphere(AnalyticSolid):
    kind = "sphere"

    def __init__(self, center, r):
        center = np.asarray(center, dtype=float)
        if center.shape != (3,) or not np.all(np.isfinite(center)):
            raise ValueError("sphere center must be a finite 3-vector")
        if not (np.isfinite(r) and r > 0):
            raise ValueError("sphere radius must be positive")
        self.center = center
        self.r = float(r)

    def implicit(self, p):
        p = np.asarray(p, dtype=float)
        return np.linalg.norm(p - self.center, axis=-1) - self.r

    def tessellate(self, tol):
        if not tol > 0:
            raise ValueError("tol must be positive")
        V, F, bound = _certified_tessellate(
            self, tol,
            lambda t, k: _sphere_lattice_mesh(self.center, self.r, t, k))
        return V, F, bound

    def bbox(self):
        e = np.full(3, self.r)
        return self.center - e, self.center + e

    def volume_exact(self):
        return 4.0 / 3.0 * math.pi * self.r ** 3

    def params_equal(self, other):
        return (isinstance(other, Sphere)
                and np.array_equal(self.center, other.center)
                and self.r == other.r)

    def __repr__(self):
        return f"Sphere(center={self.center.tolist()}, r={self.r})"


def _cylinder_mesh(base, axis, r, h, tol, k):
    dtheta = _angular_density(tol, r)
    n = max(12, int(math.ceil(2.0 * math.pi / dtheta)) * (2 ** k))
    th = np.arange(n) * (2.0 * math.pi / n)
    a = axis
    helper = np.array([1.0, 0.0, 0.0]) if abs(a[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(a, helper); u /= np.linalg.norm(u)
    v = np.cross(a, u)
    ring_pts = (base + r * (np.cos(th)[:, None] * u + np.sin(th)[:, None] * v))
    top = ring_pts + a * h
    V = np.vstack([ring_pts, top,
                   base.reshape(1, 3), (base + a * h).reshape(1, 3)])
    cb, ct = 2 * n, 2 * n + 1
    F = []
    for j in range(n):
        j2 = (j + 1) % n
        F.append([j, j2, n + j2])      # side (outward)
        F.append([j, n + j2, n + j])
        F.append([cb, j2, j])          # bottom cap (outward = -axis)
        F.append([ct, n + j, n + j2])  # top cap (outward = +axis)
    return V, np.array(F)


class Cylinder(AnalyticSolid):
    kind = "cylinder"

    def __init__(self, base_center, axis, r, h):
        base_center = np.asarray(base_center, dtype=float)
        axis = np.asarray(axis, dtype=float)
        n = np.linalg.norm(axis)
        if not (np.all(np.isfinite(base_center)) and np.isfinite(n) and n > 0):
            raise ValueError("cylinder base/axis must be finite, axis nonzero")
        if not (np.isfinite(r) and r > 0 and np.isfinite(h) and h > 0):
            raise ValueError("cylinder radius and height must be positive")
        self.base = base_center
        self.raw_axis = axis  # as given, unnormalized: exact Tier A tests
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
        if not tol > 0:
            raise ValueError("tol must be positive")
        V, F, bound = _certified_tessellate(
            self, tol,
            lambda t, k: _cylinder_mesh(self.base, self.axis, self.r, self.h, t, k))
        return V, F, bound

    def bbox(self):
        top = self.base + self.axis * self.h
        lo = np.minimum(self.base, top) - self.r
        hi = np.maximum(self.base, top) + self.r
        return lo, hi

    def volume_exact(self):
        return math.pi * self.r ** 2 * self.h

    def params_equal(self, other):
        return (isinstance(other, Cylinder)
                and np.array_equal(self.base, other.base)
                and np.array_equal(self.axis, other.axis)
                and self.r == other.r and self.h == other.h)

    def __repr__(self):
        return (f"Cylinder(base={self.base.tolist()}, axis={self.axis.tolist()}, "
                f"r={self.r}, h={self.h})")


def _cone_mesh(base, axis, r, h, tol, k):
    dtheta = _angular_density(tol, r)
    n = max(12, int(math.ceil(2.0 * math.pi / dtheta)) * (2 ** k))
    th = np.arange(n) * (2.0 * math.pi / n)
    a = axis
    helper = np.array([1.0, 0.0, 0.0]) if abs(a[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(a, helper); u /= np.linalg.norm(u)
    v = np.cross(a, u)
    ring = (base + r * (np.cos(th)[:, None] * u + np.sin(th)[:, None] * v))
    apex = (base + a * h).reshape(1, 3)
    V = np.vstack([ring, apex, base.reshape(1, 3)])
    ai, cb = n, n + 1
    F = []
    for j in range(n):
        j2 = (j + 1) % n
        F.append([j, j2, ai])      # side (outward)
        F.append([cb, j2, j])      # base cap (outward = -axis)
    return V, np.array(F)


class Cone(AnalyticSolid):
    kind = "cone"
    # implicit = max(radial - r(1-t/h), -t, t-h). The first term has
    # gradient radial_unit + (r/h)*axis, of norm sqrt(1+(r/h)^2), which
    # dominates the other two (norm 1): that is the Lipschitz constant.
    # Set per-instance in __init__.

    def __init__(self, base_center, axis, r, h):
        base_center = np.asarray(base_center, dtype=float)
        axis = np.asarray(axis, dtype=float)
        n = np.linalg.norm(axis)
        if not (np.all(np.isfinite(base_center)) and np.isfinite(n) and n > 0):
            raise ValueError("cone base/axis must be finite, axis nonzero")
        if not (np.isfinite(r) and r > 0 and np.isfinite(h) and h > 0):
            raise ValueError("cone radius and height must be positive")
        self.base = base_center
        self.raw_axis = axis  # as given, unnormalized: exact Tier A tests
        self.axis = axis / n
        self.r = float(r)
        self.h = float(h)
        self.lipschitz = math.sqrt(1.0 + (self.r / self.h) ** 2)

    def implicit(self, p):
        """Exact implicit for a right circular cone (base disc at t=0, apex t=h)."""
        p = np.asarray(p, dtype=float)
        rel = p - self.base
        t = rel @ self.axis
        radial_vec = rel - np.outer(t, self.axis)
        radial = np.linalg.norm(radial_vec, axis=-1)
        r_at = self.r * np.maximum(0.0, 1.0 - t / self.h)
        return np.maximum.reduce([
            radial - r_at,
            -t,
            t - self.h,
        ])

    def tessellate(self, tol):
        if not tol > 0:
            raise ValueError("tol must be positive")
        V, F, bound = _certified_tessellate(
            self, tol,
            lambda t, k: _cone_mesh(self.base, self.axis, self.r, self.h, t, k))
        return V, F, bound

    def bbox(self):
        top = self.base + self.axis * self.h
        lo = np.minimum(self.base, top) - self.r
        hi = np.maximum(self.base, top) + self.r
        return lo, hi

    def volume_exact(self):
        return math.pi * self.r ** 2 * self.h / 3.0

    def params_equal(self, other):
        return (isinstance(other, Cone)
                and np.array_equal(self.base, other.base)
                and np.array_equal(self.axis, other.axis)
                and self.r == other.r and self.h == other.h)

    def __repr__(self):
        return (f"Cone(base={self.base.tolist()}, axis={self.axis.tolist()}, "
                f"r={self.r}, h={self.h})")
