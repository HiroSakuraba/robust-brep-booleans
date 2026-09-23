"""Stages 3-4: Tier A exact analysis + per-face audit of the engine.

Each result face lies on an input surface, so it is classified against the
*other* solid only (origin tracked per triangle from the engine), and the
engine's keep/discard decision is audited against the exact implicit:

  union:        kept A-face must satisfy B_implicit >= 0; kept B-face: A >= 0
  intersection: kept A-face must satisfy B_implicit <= 0; kept B-face: A <= 0
  difference:   kept A-face must satisfy B_implicit >= 0; kept B-face: A <= 0

Per-face verdict (s = required polarity, L = other solid's Lipschitz
constant, m = margin):
  * centroid strictly on the kept side (s*f > m) AND an adaptive
    subdivision finds no region of the face provably beyond the margin on
    the wrong side  -> VERIFIED
  * subdivision proves some sub-triangle lies wholly beyond the margin on
    the wrong side (s*f(c) + L*r < -m)  -> VIOLATION (engine is wrong;
    this catches "forgot to cut" cases like a bump poking through a face,
    where the centroid alone can sit outside the intruding solid)
  * otherwise -> AMBIGUOUS (inside the margin band; cannot verify)

Faces are then grouped into patches (edge-adjacent, same origin), split
wherever A-faces meet B-faces. A patch touching such a split carries a
transverse crossing, so its in-band slivers take the patch verdict: the
patch is verified if some face verifies decisively and none violate (the
residual error is geometric, bounded by the margin, and cannot change
topology). A patch with no crossing gets no such transfer: every face
must verify on its own, so tangencies (tangent cylinders/spheres, point
or edge box contacts) with undecidable faces still block.

Tier A (below) detects degeneracies from defining parameters with exact
rational predicates -- no tolerance voting.
"""

from fractions import Fraction

import numpy as np


def _Fr(x):
    """Exact rational for a float input."""
    return Fraction(x)


def _point_segment_distance(P, A, B):
    """Distance from P[i] to segment A[i]B[i]. P,A,B: (n,3)."""
    ab = B - A
    denom = np.maximum(np.einsum('ij,ij->i', ab, ab), 1e-300)
    t = np.clip(np.einsum('ij,ij->i', P - A, ab) / denom, 0.0, 1.0)
    Q = A + t[:, None] * ab
    return np.linalg.norm(P - Q, axis=1)


def _point_triangle_dist_lb(P, T):
    """Sound lower bound on distance from P[i] to triangle T[i].

    Uses the standard barycentric classification when the triangle is
    well-conditioned (exact up to a 1e-12 fp slop), and falls back to
    min(plane distance, edge distance) -- always sound, merely
    conservative -- for degenerate or needle triangles. Taking the
    minimum with the edge distance keeps the bound sound even if
    floating-point noise misclassifies a projection that falls almost
    exactly on the triangle boundary (the resulting overestimate is
    O(fp noise), covered by the slop).
    """
    A, B, C = T[:, 0], T[:, 1], T[:, 2]
    ab = B - A
    ac = C - A
    ap = P - A
    d00 = np.einsum('ij,ij->i', ab, ab)
    d01 = np.einsum('ij,ij->i', ab, ac)
    d11 = np.einsum('ij,ij->i', ac, ac)
    d20 = np.einsum('ij,ij->i', ap, ab)
    d21 = np.einsum('ij,ij->i', ap, ac)
    den = d00 * d11 - d01 * d01  # |ab x ac|^2 >= 0
    # well-conditioned iff the angle at A is not tiny
    well = den > 1e-8 * np.maximum(d00 * d11, 1e-300)
    den_safe = np.where(well, den, 1.0)
    u = (d11 * d20 - d01 * d21) / den_safe
    v = (d00 * d21 - d01 * d20) / den_safe
    # Scale-aware uncertainty band: fp error in (u,v) grows like
    # |ap|/|edge| / sin^2(theta); the band below covers it with wide
    # headroom, so a "clear" classification cannot be a fp artifact.
    band = 1e-7 * (1.0 + (np.abs(d20) + np.abs(d21))
                   / np.maximum(den, 1e-300))
    clearly_in = (well & (u > band) & (v > band)
                  & (u + v < 1.0 - band))
    clearly_out = (well & ((u < -band) | (v < -band)
                           | (u + v > 1.0 + band)))
    n = np.cross(ab, ac)
    n2 = np.einsum('ij,ij->i', n, n)
    t = np.einsum('ij,ij->i', ap, n) / np.maximum(n2, 1e-300)
    plane_dist = np.abs(t) * np.sqrt(n2)
    edge = _point_segment_distance(P, A, B)
    edge = np.minimum(edge, _point_segment_distance(P, B, C))
    edge = np.minimum(edge, _point_segment_distance(P, C, A))
    # inside -> plane distance is exact; outside -> edge distance is
    # exact; near the boundary or ill-conditioned -> min of the two,
    # which is always a sound lower bound
    d = np.where(clearly_in, plane_dist,
                 np.where(clearly_out, edge,
                          np.minimum(plane_dist, edge)))
    return d - 1e-12 * (1.0 + d)


def _face_lo(T, solid, s):
    """Rigorous lower bound on min(s * f) over each triangle.

    T: (n,3,3) triangles, f: solid's implicit, s: required polarity (+1/-1).
    A face with lo > margin is wholly on the kept side (VERIFIED); a face
    with lo >= -margin cannot contain a provable violation, so the costly
    subdivision only runs when lo < -margin. Bounds used:
      box:      f = max of linear pieces; min(max l_k) >= max_k min_vert l_k,
                max(max l_k) = max_k max_vert l_k (linear -> vertices exact)
      sphere:   min |p-c| = point-triangle distance (exact lb above);
                max |p-c| at a vertex (convex)
      cylinder: f = max(radial - r, -t, t - h); radial is convex, t linear
      cone:     f = max(|q| + (r/h) t - r, -t, t - h); the first piece is
                convex, so its max is at a vertex; its min >= dist_lb(|q|)
                + (r/h) min_vert(t) - r
    """
    from .solids import _perp_basis
    n = len(T)
    kind = solid.kind
    if kind == "box":
        p0 = solid.lo[None, None, :] - T
        p1 = T - solid.hi[None, None, :]
        pieces = np.concatenate([p0, p1], axis=2)  # (n,3,6)
        if s > 0:
            return np.max(np.min(pieces, axis=1), axis=1)
        return -np.max(np.max(pieces, axis=1), axis=1)
    if kind == "sphere":
        c = np.broadcast_to(solid.center, (n, 3))
        if s > 0:
            return _point_triangle_dist_lb(c, T) - solid.r
        maxd = np.max(np.linalg.norm(T - c[:, None, :], axis=2), axis=1)
        return solid.r - maxd
    if kind in ("cylinder", "cone"):
        u, v = _perp_basis(solid.axis)
        rel = T - solid.base[None, None, :]
        tax = rel @ solid.axis                        # (n,3) axial coords
        q = np.stack([rel @ u, rel @ v], axis=-1)     # (n,3,2)
        q3 = np.zeros((n, 3, 3))
        q3[..., :2] = q
        slope = solid.r / solid.h if kind == "cone" else 0.0
        if s > 0:
            d_lb = _point_triangle_dist_lb(np.zeros((n, 3)), q3)
            lo_a = d_lb - solid.r + slope * np.min(tax, axis=1)
            lo_b = np.min(-tax, axis=1)
            lo_c = np.min(tax - solid.h, axis=1)
            return np.maximum.reduce([lo_a, lo_b, lo_c])
        # s < 0: max_T f is exact at vertices (convex pieces)
        fv = solid.implicit(T.reshape(-1, 3)).reshape(n, 3)
        return -np.max(fv, axis=1)
    raise ValueError(f"no face-bound rule for {kind}")


def _face_has_violation(tri, implicit_other, s, L, margin, budget=20000):
    """Search the triangle for a region provably beyond the margin on the
    wrong side. Returns True only on proof (sound); False means no proof
    was found within the budget/margin-scale limits.

    Prunes a sub-triangle when s*f(c) - L*r >= -margin: no point of it can
    lie beyond the band on the wrong side. Proves violation when
    s*f(c) + L*r < -margin: every point of it does. Stops below the margin
    scale (r <= margin/L), where the audit cannot decide by construction.
    """
    r_min = margin / L
    stack = [(tri, 0)]
    nodes = 0
    while stack:
        t, depth = stack.pop()
        nodes += 1
        if nodes > budget:
            return False
        c = t.mean(axis=0)
        r = float(np.max(np.linalg.norm(t - c, axis=1)))
        g = s * float(implicit_other(c.reshape(1, 3))[0])
        if g + L * r < -margin:
            return True
        if g - L * r >= -margin:
            continue
        if r <= r_min or depth >= 14:
            continue
        m01 = (t[0] + t[1]) * 0.5
        m12 = (t[1] + t[2]) * 0.5
        m20 = (t[2] + t[0]) * 0.5
        stack.append((np.array([t[0], m01, m20]), depth + 1))
        stack.append((np.array([t[1], m12, m01]), depth + 1))
        stack.append((np.array([t[2], m20, m12]), depth + 1))
        stack.append((np.array([m01, m12, m20]), depth + 1))
    return False


def _patches(F, origins):
    """Group faces into patches: edge-adjacent and same origin.

    Returns (patch_id per face, touches_crossing per face): a patch
    touches a crossing if any of its faces shares an edge with a face of
    a different origin.
    """
    n = len(F)
    E = np.sort(np.concatenate(
        [F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]), axis=1)
    all_e = E
    all_f = np.tile(np.arange(n), 3)  # E stacks the 3 edge slots per face
    order = np.lexsort((all_e[:, 1], all_e[:, 0]))
    all_e = all_e[order]
    all_f = all_f[order]
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    touch = np.zeros(n, dtype=bool)
    m = len(all_e)
    i = 0
    while i < m:
        j = i + 1
        ei = all_e[i]
        while j < m and all_e[j, 0] == ei[0] and all_e[j, 1] == ei[1]:
            j += 1
        faces = all_f[i:j]
        if len(faces) > 1:
            for o in np.unique(origins[faces]):
                ff = faces[origins[faces] == o]
                for f in ff[1:]:
                    union(int(ff[0]), int(f))
            if len(np.unique(origins[faces])) > 1:
                touch[faces] = True
        i = j
    roots = np.array([find(f) for f in range(n)])
    _, pid = np.unique(roots, return_inverse=True)
    pt = np.zeros(int(pid.max()) + 1, dtype=bool)
    np.logical_or.at(pt, pid, touch)
    return pid, pt[pid]


def audit_faces(F, V, origins, solidA, solidB, margin, op):
    """Audit the engine's keep decision per face against exact implicits.

    Returns dict with boolean masks: verified, ambiguous, violation, plus
    the raw other-solid implicit values and the patch count.
    """
    F = np.asarray(F)
    n = len(F)
    z = np.zeros(n, dtype=bool)
    if n == 0:
        return {"verified": z, "ambiguous": z, "violation": z,
                "f_other": np.zeros(0), "origins": origins, "n_patches": 0}
    origins = np.asarray(origins)
    T = V[F]  # (n,3,3)
    fromA = origins == "A"
    L = np.where(fromA, solidB.lipschitz, solidA.lipschitz)
    # required polarity s: keep ⟹ s * f_other >= 0
    if op == "union":
        s = np.ones(n)
    elif op == "intersection":
        s = -np.ones(n)
    elif op == "difference":
        s = np.where(fromA, 1.0, -1.0)
    else:
        raise ValueError(f"unknown op {op!r}")

    # Rigorous lower bound on min(s*f_other) over each face. lo > margin:
    # the whole face is decisively on the kept side. lo >= -margin: no
    # sub-triangle can be provably wrong-side, so the subdivision is
    # skipped and the face is at worst ambiguous.
    lo = np.empty(n)
    if np.any(fromA):
        lo[fromA] = _face_lo(T[fromA], solidB, float(s[fromA][0]))
    if np.any(~fromA):
        lo[~fromA] = _face_lo(T[~fromA], solidA, float(s[~fromA][0]))
    f_other = np.empty(n)  # centroid values, for the report only
    if np.any(fromA):
        f_other[fromA] = solidB.implicit(T[fromA].mean(axis=1))
    if np.any(~fromA):
        f_other[~fromA] = solidA.implicit(T[~fromA].mean(axis=1))

    verified = lo > margin
    violated = np.zeros(n, dtype=bool)
    # Centroid itself wrong-side beyond the margin: that point of the face
    # is provably mis-kept.
    g = s * f_other
    violated[g < -margin] = True
    # Otherwise, a violation is provable only if lo < -margin; run the
    # subdivision solely in that case. If it cannot prove one (budget),
    # the face stays ambiguous -- never verified on an unproven dip.
    for i in np.nonzero((lo < -margin) & ~violated)[0]:
        i = int(i)
        imp = solidB.implicit if fromA[i] else solidA.implicit
        if _face_has_violation(T[i], imp, float(s[i]), float(L[i]), margin):
            violated[i] = True
    # No violation provable and centroid decisive: the engine's keep
    # decision is certified (this covers faces that touch the other
    # solid's surface at their boundary, which a whole-face bound
    # can never verify).
    verified[(lo >= -margin) & (g > margin)] = True
    ambiguous = ~(verified | violated)

    # Patch-level verdicts.
    pid, ptouch = _patches(F, origins)
    n_patches = int(pid.max()) + 1
    pv = np.zeros(n_patches, dtype=np.int8)  # 0 amb, 1 ok, 2 violated
    p_viol = np.zeros(n_patches, dtype=bool)
    p_has_ok = np.zeros(n_patches, dtype=bool)
    p_touch = np.zeros(n_patches, dtype=bool)
    np.logical_or.at(p_viol, pid, violated)
    np.logical_or.at(p_has_ok, pid, verified)
    np.logical_or.at(p_touch, pid, ptouch)
    n_per_patch = np.bincount(pid)
    n_ok_per_patch = np.bincount(pid, weights=verified.astype(int))
    p_all_ok = n_ok_per_patch == n_per_patch
    ok = ~p_viol
    pv[ok & p_touch & p_has_ok] = 1       # crossing patch: slivers rescued
    pv[ok & ~p_touch & p_all_ok] = 1      # no crossing: every face verified
    pv[p_viol] = 2
    verified = pv[pid] == 1
    violation = pv[pid] == 2
    ambiguous = ~(verified | violation)
    return {"verified": verified, "ambiguous": ambiguous,
            "violation": violation, "f_other": f_other, "origins": origins,
            "n_patches": n_patches}


def identical_inputs(solidA, solidB):
    """Exact fast path: bitwise-identical defining parameters."""
    return solidA.params_equal(solidB)


def _axis_relations(a_lo, a_hi, b_lo, b_hi):
    """Per-axis interval relation: 'disjoint', 'touch', or 'overlap'."""
    rel = []
    for ax in range(3):
        if a_hi[ax] < b_lo[ax] or b_hi[ax] < a_lo[ax]:
            rel.append("disjoint")
        elif a_hi[ax] == b_lo[ax] or b_hi[ax] == a_lo[ax]:
            rel.append("touch")
        else:
            rel.append("overlap")
    return rel


def tierA_degeneracy(solidA, solidB, margin):
    """Exact degeneracy detection from the *defining parameters*.

    Findings are informational: they explain why faces may be
    unverifiable, but they never block on their own. Blocking comes
    from the per-face audit and Stage 6. No tolerance voting: every
    predicate here is exact (rational arithmetic on the raw inputs;
    float == is already exact).
    """
    from .solids import Box, Sphere, Cylinder
    findings = []
    if isinstance(solidA, Box) and isinstance(solidB, Box):
        rel = _axis_relations(solidA.lo, solidA.hi, solidB.lo, solidB.hi)
        if "disjoint" in rel:
            findings.append("disjoint_boxes")
        else:
            touches = [ax for ax, r in enumerate(rel) if r == "touch"]
            if touches:
                kinds = {3: "point", 2: "edge", 1: "face"}[len(touches)]
                findings.append(
                    f"touching_boxes: {kinds} contact on axes {touches}")
            # coincident face planes, only meaningful when closures meet
            for ax in range(3):
                if solidA.lo[ax] == solidB.lo[ax]:
                    findings.append(f"coincident_planes: A.lo[{ax}]==B.lo[{ax}]")
                if solidA.hi[ax] == solidB.hi[ax]:
                    findings.append(f"coincident_planes: A.hi[{ax}]==B.hi[{ax}]")
    elif isinstance(solidA, Sphere) and isinstance(solidB, Sphere):
        d2 = sum((_Fr(a) - _Fr(b)) ** 2
                 for a, b in zip(solidA.center, solidB.center))
        rs = _Fr(solidA.r) + _Fr(solidB.r)
        rd = _Fr(solidA.r) - _Fr(solidB.r)
        if d2 == 0:
            findings.append("concentric_spheres")
        elif d2 == rs * rs:
            findings.append("tangent_spheres: externally tangent")
        elif d2 == rd * rd:
            findings.append("tangent_spheres: internally tangent")
    elif isinstance(solidA, Cylinder) and isinstance(solidB, Cylinder):
        a = [_Fr(x) for x in solidA.raw_axis]
        b = [_Fr(x) for x in solidB.raw_axis]
        cr = [a[1] * b[2] - a[2] * b[1],
              a[2] * b[0] - a[0] * b[2],
              a[0] * b[1] - a[1] * b[0]]
        if cr[0] == 0 and cr[1] == 0 and cr[2] == 0:  # exactly parallel
            w = [_Fr(q) - _Fr(p)
                 for p, q in zip(solidA.base, solidB.base)]
            aa = sum(x * x for x in a)
            wad = sum(x * y for x, y in zip(w, a))
            w2 = sum(x * x for x in w)
            # |w_perp|^2 = |w|^2 - (w.a)^2/(a.a), exact rational
            dist2 = w2 - wad * wad / aa
            rr = _Fr(solidA.r) + _Fr(solidB.r)
            if dist2 == rr * rr:
                findings.append("tangent_cylinders: axis distance == rA+rB")
            elif dist2 == 0 and solidA.r == solidB.r:
                findings.append("coincident_cylindrical_surfaces")
    return findings


def exact_separation(solidA, solidB):
    """Exactly decide disjointness where a closed-form predicate exists.

    Returns True (disjoint), False (not disjoint), or None (unknown).
    """
    from .solids import Box, Sphere
    if isinstance(solidA, Box) and isinstance(solidB, Box):
        rel = _axis_relations(solidA.lo, solidA.hi, solidB.lo, solidB.hi)
        if "disjoint" in rel:
            return True
        return False
    if isinstance(solidA, Sphere) and isinstance(solidB, Sphere):
        d2 = sum((_Fr(a) - _Fr(b)) ** 2
                 for a, b in zip(solidA.center, solidB.center))
        rs = _Fr(solidA.r) + _Fr(solidB.r)
        return d2 > rs * rs
    # bbox-level: disjoint bboxes are exactly disjoint
    alo, ahi = solidA.bbox()
    blo, bhi = solidB.bbox()
    if np.any(ahi < blo) or np.any(bhi < alo):
        return True
    return None


def _box_grid(solidA, solidB, op):
    """Exact (n_shells, chi) for box-box, from the cell grid.

    The 12 face planes cut space into a grid of cells; each cell is
    exactly inside/outside each closed box (no cell straddles a face).
    Keep cells per op, take boundary faces, count edge-connected
    boundary components (shells) and chi = V - E + F over the boundary.
    Returns (None, None) for non-box inputs.
    """
    from .solids import Box
    if not (isinstance(solidA, Box) and isinstance(solidB, Box)):
        return None, None
    alo, ahi = solidA.lo, solidA.hi
    blo, bhi = solidB.lo, solidB.hi
    px = sorted({alo[0], ahi[0], blo[0], bhi[0]})
    py = sorted({alo[1], ahi[1], blo[1], bhi[1]})
    pz = sorted({alo[2], ahi[2], blo[2], bhi[2]})
    xs = [(px[i], px[i + 1]) for i in range(len(px) - 1) if px[i] < px[i + 1]]
    ys = [(py[i], py[i + 1]) for i in range(len(py) - 1) if py[i] < py[i + 1]]
    zs = [(pz[i], pz[i + 1]) for i in range(len(pz) - 1) if pz[i] < pz[i + 1]]
    nx, ny, nz = len(xs), len(ys), len(zs)
    kept = np.zeros((nx, ny, nz), dtype=bool)
    for i, (x0, x1) in enumerate(xs):
        ain = alo[0] <= x0 and x1 <= ahi[0]
        bin = blo[0] <= x0 and x1 <= bhi[0]
        for j, (y0, y1) in enumerate(ys):
            ain_y = ain and alo[1] <= y0 and y1 <= ahi[1]
            bin_y = bin and blo[1] <= y0 and y1 <= bhi[1]
            for k, (z0, z1) in enumerate(zs):
                a = ain_y and alo[2] <= z0 and z1 <= ahi[2]
                b = bin_y and blo[2] <= z0 and z1 <= bhi[2]
                if op == "union":
                    kept[i, j, k] = a or b
                elif op == "intersection":
                    kept[i, j, k] = a and b
                elif op == "difference":
                    kept[i, j, k] = a and not b
                else:
                    raise ValueError(op)
    if not kept.any():
        return 0, 0
    # boundary faces: (axis, i, j, k); axis 0 -> x-plane px[i], etc.
    faces = []
    for i in range(nx + 1):
        for j in range(ny):
            for k in range(nz):
                left = kept[i - 1, j, k] if i > 0 else False
                right = kept[i, j, k] if i < nx else False
                if left != right:
                    faces.append((0, i, j, k))
    for i in range(nx):
        for j in range(ny + 1):
            for k in range(nz):
                dn = kept[i, j - 1, k] if j > 0 else False
                up = kept[i, j, k] if j < ny else False
                if dn != up:
                    faces.append((1, i, j, k))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz + 1):
                bk = kept[i, j, k - 1] if k > 0 else False
                fr = kept[i, j, k] if k < nz else False
                if bk != fr:
                    faces.append((2, i, j, k))

    def face_edges(ax, i, j, k):
        if ax == 0:
            return [("y", i, j, k), ("y", i, j, k + 1),
                    ("z", i, j, k), ("z", i, j + 1, k)]
        if ax == 1:
            return [("x", i, j, k), ("x", i, j, k + 1),
                    ("z", i, j, k), ("z", i + 1, j, k)]
        return [("x", i, j, k), ("x", i, j + 1, k),
                ("y", i, j, k), ("y", i + 1, j, k)]

    def face_verts(ax, i, j, k):
        if ax == 0:
            return [(i, j, k), (i, j + 1, k), (i, j, k + 1), (i, j + 1, k + 1)]
        if ax == 1:
            return [(i, j, k), (i + 1, j, k), (i, j, k + 1), (i + 1, j, k + 1)]
        return [(i, j, k), (i + 1, j, k), (i, j + 1, k), (i + 1, j + 1, k)]

    edge_to_faces = {}
    for fid, (ax, i, j, k) in enumerate(faces):
        for e in face_edges(ax, i, j, k):
            edge_to_faces.setdefault(e, []).append(fid)
    parent = list(range(len(faces)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for flist in edge_to_faces.values():
        for f in flist[1:]:
            ra, rb = find(flist[0]), find(f)
            if ra != rb:
                parent[ra] = rb
    n_shells = len({find(f) for f in range(len(faces))})
    verts, edges = set(), set()
    for (ax, i, j, k) in faces:
        verts.update(face_verts(ax, i, j, k))
        edges.update(face_edges(ax, i, j, k))
    chi = len(verts) - len(edges) + len(faces)
    return n_shells, chi


def expected_shells(solidA, solidB, op):
    """Expected number of connected triangle-shells, or None if unknown."""
    from .solids import Box
    if isinstance(solidA, Box) and isinstance(solidB, Box):
        n, _ = _box_grid(solidA, solidB, op)
        return n
    sep = exact_separation(solidA, solidB)
    if op == "union":
        if sep is True:
            return 2
        if sep is False:
            return 1
        return None
    if op == "intersection":
        return None  # decided by emptiness downstream
    if op == "difference":
        if sep is True:
            return 1  # A - B = A
        return None  # relation not exactly decidable: report only
    raise ValueError(op)


def expected_euler(solidA, solidB, op):
    """Predicted Euler characteristic, or None if not exactly decidable.

    Exact for box-box, derived from the cell grid (handles cavities,
    splits, and tunnels that interval relations get wrong).
    """
    from .solids import Box
    if isinstance(solidA, Box) and isinstance(solidB, Box):
        _, chi = _box_grid(solidA, solidB, op)
        return chi
    return None


def exact_op_volume(solidA, solidB, op):
    """Closed-form result volume, or None if not exactly decidable.

    Box-box: the intersection is a box (possibly empty/degenerate), so
    union / intersection / difference volumes are exact via inclusion.
    """
    from .solids import Box
    if not (isinstance(solidA, Box) and isinstance(solidB, Box)):
        return None
    lo = np.maximum(solidA.lo, solidB.lo)
    hi = np.minimum(solidA.hi, solidB.hi)
    v_inter = float(np.prod(np.maximum(hi - lo, 0.0)))
    vA, vB = solidA.volume_exact(), solidB.volume_exact()
    if op == "union":
        return vA + vB - v_inter
    if op == "intersection":
        return v_inter
    if op == "difference":
        return vA - v_inter
    raise ValueError(op)
