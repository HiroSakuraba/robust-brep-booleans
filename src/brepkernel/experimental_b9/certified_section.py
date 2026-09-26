"""Certified section-curve enclosure for analytic surface pairs (EXPERIMENTAL, B9).

Question under test: can we compute a section curve between two analytic
surfaces with a *certificate* instead of OCCT tolerances plus refusal
gates?

What this module does, for one pair of implicit surfaces (F, G) on a
bounded domain box D:

1. Exclusion subdivision (rigorous): recursively split D into boxes.
   A box is discarded only when interval evaluation proves 0 is not in
   F(box) or not in G(box). Every discarded box carries its exclusion
   proof. The union of the retained boxes is therefore a *rigorous
   superset* of the true intersection curve: no branch can be missed.
   This is the certified answer to the plan's "open completeness
   problem" for analytic pairs: a missed loop would have to hide in a
   discarded box, which the interval proof forbids.

2. Transversality flags: on a retained box, if the interval cross
   product grad F x grad G excludes the zero vector, the curve inside
   is regular (Jacobian rank 2 everywhere in the box). No singular
   point hides there.

3. Krawczyk existence/uniqueness per box: on a transverse box, one
   Krawczyk step on the square system (F, G, t.(x - m)) with t the
   approximate tangent at the box center proves a curve point in the
   slicing plane ("exists"), or, when the Krawczyk image lands
   strictly inside the box, exactly one ("unique"). This is the
   "interval Newton / Krawczyk steps that prove a unique curve piece
   in each box" step from the plan's B9 section. With the
   transversality flag, "unique" means a unique regular arc.

4. Refusal instead of guessing: boxes that reach the minimum size
   without being excludable or certifiable are reported as
   *unresolved* (singular or near-singular region). The enclosure never
   claims more than it proved.

What this module does NOT do: freeform NURBS (needs the plan's Bezier
decomposition first), coincident/overlapping surfaces (F and G
proportional: every box is retained, reported as degenerate input),
or trimming (the enclosure is in 3D; trimming to face domains is a
separate step).

Agreement protocol with OCCT (the plan's "run both intersectors"):
every 3D sample point of OCCT's section edges must lie inside the
retained enclosure; a sample outside is a soundness alarm (a bug in
this code or in OCCT, never silently ignored). Retained-box clusters
with no OCCT sample nearby are flagged as unexplained retained
regions: over-approximation the tracer may have missed, or a genuine
missed loop. Either way it is reported, not certified.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .intervals import Interval


# ---------------------------------------------------------------------------
# Box type

@dataclass
class Box3:
    lo: np.ndarray  # (3,) float
    hi: np.ndarray  # (3,) float

    def diag(self) -> float:
        return float(np.linalg.norm(self.hi - self.lo))

    def center(self) -> np.ndarray:
        return 0.5 * (self.lo + self.hi)

    def intervals(self) -> tuple[Interval, Interval, Interval]:
        return tuple(Interval(float(a), float(b))
                     for a, b in zip(self.lo, self.hi))

    def children(self) -> list["Box3"]:
        m = self.center()
        out = []
        for i in range(8):
            lo = self.lo.copy()
            hi = self.hi.copy()
            for ax in range(3):
                if i & (1 << ax):
                    lo[ax] = m[ax]
                else:
                    hi[ax] = m[ax]
            out.append(Box3(lo, hi))
        return out


# ---------------------------------------------------------------------------
# Interval linear algebra (3x3 / 3-vector, minimal)

def _imat_zeros():
    return [[Interval(0.0, 0.0) for _ in range(3)] for _ in range(3)]


def _imat_zeros2():
    return [[Interval(0.0, 0.0) for _ in range(2)] for _ in range(2)]


def _fmat_imat_mul(C: np.ndarray, J) -> list:
    """Float 3x3 times interval 3x3 -> interval 3x3."""
    R = _imat_zeros()
    for i in range(3):
        for j in range(3):
            s = Interval(0.0, 0.0)
            for k in range(3):
                s = s + J[k][j] * float(C[i, k])
            R[i][j] = s
    return R


def _fmat_imat_mul2(C: np.ndarray, J) -> list:
    """Float 2x2 times interval 2x2 -> interval 2x2."""
    R = _imat_zeros2()
    for i in range(2):
        for j in range(2):
            s = Interval(0.0, 0.0)
            for k in range(2):
                s = s + J[k][j] * float(C[i, k])
            R[i][j] = s
    return R


def _fmat_ivec_mul(C: np.ndarray, v) -> list:
    R = [Interval(0.0, 0.0) for _ in range(3)]
    for i in range(3):
        s = Interval(0.0, 0.0)
        for k in range(3):
            s = s + v[k] * float(C[i, k])
        R[i] = s
    return R


def _fmat_ivec_mul2(C: np.ndarray, v) -> list:
    R = [Interval(0.0, 0.0) for _ in range(2)]
    for i in range(2):
        s = Interval(0.0, 0.0)
        for k in range(2):
            s = s + v[k] * float(C[i, k])
        R[i] = s
    return R


def _imat_ivec_mul(M, v) -> list:
    R = [Interval(0.0, 0.0) for _ in range(3)]
    for i in range(3):
        s = Interval(0.0, 0.0)
        for k in range(3):
            s = s + M[i][k] * v[k]
        R[i] = s
    return R


def _imat_ivec_mul2(M, v) -> list:
    R = [Interval(0.0, 0.0) for _ in range(2)]
    for i in range(2):
        s = Interval(0.0, 0.0)
        for k in range(2):
            s = s + M[i][k] * v[k]
        R[i] = s
    return R


# ---------------------------------------------------------------------------
# Transversality and Krawczyk certificates

def transverse_certificate(f, g, box: Box3) -> tuple[bool, tuple]:
    """True when grad F x grad G provably avoids 0 on the whole box.

    Returns (ok, normal_intervals). A True result means the Jacobian of
    (F, G) has rank 2 at every common zero in the box: no singular
    point of the intersection curve lies inside.
    """
    X = box.intervals()
    gf = f.interval_grad(X)
    gg = g.interval_grad(X)
    n = (
        gf[1] * gg[2] - gf[2] * gg[1],
        gf[2] * gg[0] - gf[0] * gg[2],
        gf[0] * gg[1] - gf[1] * gg[0],
    )
    ok = any(not c.contains_zero() for c in n)
    return ok, n


def _newton_refine(f, g, m, t, steps=12):
    """Gauss-Newton closest point on {F = G = 0} to m.

    Solves min ||dx|| s.t. J dx = -(F, G), i.e.
    dx = -J^T (J J^T)^{-1} F. This converges to the NEAREST curve
    point, so the limit x1 lies in the box whenever the box
    genuinely contains curve near its center (the old slicing-plane
    Newton could converge to a far intersection outside the box).
    Returns (x1, t1) with the unit tangent at x1, or (None, None).
    """
    x = np.array(m, dtype=float)
    for _ in range(steps):
        Fv = np.array([f.float_value(x), g.float_value(x)])
        if float(np.linalg.norm(Fv)) < 1e-14:
            break
        J = np.stack([np.asarray(f.float_grad(x), dtype=float),
                      np.asarray(g.float_grad(x), dtype=float)])
        try:
            dx = -J.T @ np.linalg.solve(J @ J.T, Fv)
        except np.linalg.LinAlgError:
            return None, None
        if not np.all(np.isfinite(dx)):
            return None, None
        x = x + dx
        if float(np.linalg.norm(dx)) < 1e-14:
            break
    if float(np.linalg.norm([f.float_value(x),
                             g.float_value(x)])) > 1e-8:
        return None, None
    gf = np.asarray(f.float_grad(x), dtype=float)
    gg = np.asarray(g.float_grad(x), dtype=float)
    t1 = np.cross(gf, gg)
    n1 = float(np.linalg.norm(t1))
    if not np.isfinite(n1) or n1 < 1e-300:
        return None, None
    return x, t1 / n1


def _krawczyk_on_box(f, g, box: Box3, m, t) -> str:
    """Single Krawczyk step with given center m and tangent t.

    The step is centered at a Newton-refined point x1 (not the raw
    box center): K = x1 - C.H(x1) + (I - C.J(X))(X - x1). Centering
    at the approximate zero makes the enclosure dramatically
    tighter when the box center is off the curve.
    """
    x1, t1 = _newton_refine(f, g, m, t)
    if x1 is None:
        return "none"
    gf = np.asarray(f.float_grad(x1), dtype=float)
    gg = np.asarray(g.float_grad(x1), dtype=float)
    Jm = np.array([gf, gg, t1])
    try:
        C = np.linalg.inv(Jm)
    except np.linalg.LinAlgError:
        return "none"
    if not np.all(np.isfinite(C)):
        return "none"

    Hx = np.array([f.float_value(x1), g.float_value(x1), 0.0])
    xc = x1 - C @ Hx  # base point: x1 - C.H(x1)

    X = box.intervals()
    Jint = [list(f.interval_grad(X)), list(g.interval_grad(X)),
            [Interval.point(float(v)) for v in t1]]
    I_minus_CJ = _imat_zeros()
    CJ = _fmat_imat_mul(C, Jint)
    for i in range(3):
        for j in range(3):
            I_minus_CJ[i][j] = Interval.point(1.0 if i == j else 0.0) \
                - CJ[i][j]
    # Standard Krawczyk form: K = xc + (I - C.J(X)) (X - x1).
    d = [Interval(float(box.lo[k]) - float(x1[k]),
                  float(box.hi[k]) - float(x1[k])) for k in range(3)]
    K = [Interval.point(float(v)) for v in xc]
    corr = _imat_ivec_mul(I_minus_CJ, d)
    K = [K[k] + corr[k] for k in range(3)]
    inside = all(K[k].lo >= float(box.lo[k]) and
                 K[k].hi <= float(box.hi[k]) for k in range(3))
    if not inside:
        return "none"
    strict = all(K[k].lo > float(box.lo[k]) and
                 K[k].hi < float(box.hi[k]) for k in range(3))
    return "unique" if strict else "exists"


def krawczyk_prove(f, g, box: Box3) -> str:
    """Krawczyk certificate for the sliced system H(x) = 0 in the box.

    H(x) = (F(x), G(x), t.(x - m)) with m the box center and t the
    approximate curve tangent there. Returns:

    - "unique": K(box) lies strictly inside box. Existence AND
      uniqueness of the zero: exactly one curve point in the slicing
      plane, and with transversality a unique regular arc.
    - "exists": K(box) is contained in box but touches the boundary.
      At least one curve point in the slicing plane lies in the box
      (existence only).
    - "none": no certificate.

    Grid-alignment note: when the true curve runs exactly along a box
    face (common with axis-aligned test geometry), the sliced zero is
    on the boundary and strict containment is impossible. In that
    case the test is retried on a slightly expanded box (2% of the
    diagonal per side); a "unique" verdict then certifies a unique
    regular arc in the expanded box, which still covers this box.
    """
    m = box.center()
    gf = np.asarray(f.float_grad(m), dtype=float)
    gg = np.asarray(g.float_grad(m), dtype=float)
    t = np.cross(gf, gg)
    nt = float(np.linalg.norm(t))
    if not np.isfinite(nt) or nt < 1e-300:
        return "none"
    t = t / nt
    r = _krawczyk_on_box(f, g, box, m, t)
    if r != "none":
        return r
    # Retry on a slightly expanded box (see grid-alignment note).
    pad = 0.02 * box.diag()
    big = Box3(box.lo - pad, box.hi + pad)
    return _krawczyk_on_box(f, g, big, m, t)


def krawczyk_certificate(f, g, box: Box3) -> bool:
    """Boolean wrapper: True only for the uniqueness certificate."""
    return krawczyk_prove(f, g, box) == "unique"


def graph_verdict(f, g, box: Box3) -> str:
    """Interval implicit-function-theorem verdict: 'unique', 'excluded',
    or 'none'.

    On a transverse box, some tangent component t_k provably avoids 0,
    so the curve is a graph over coordinate k. Write u for the other
    two coordinates and run a *parametric* 2D Krawczyk step: with the
    interval residual over the full k-range folded in,
    K(U) = u - C.[]H(u, Kparam) + (I - C.[]J_u(U, Kparam)) (U - u).

    - K(U) strictly inside U: for EVERY parameter value in the box's
      k-range the 2x2 system has exactly one solution in U. The box
      holds a unique regular arc: existence and uniqueness, no
      slicing plane, independent of where the curve sits in the box.
    - K(U) disjoint from U: no parameter value admits a solution in
      U (every per-parameter zero must lie in K(U)). The box
      provably contains NO curve: a rigorous exclusion that kills
      the false positives of the plain F(B)/G(B) interval test.
    - otherwise 'none': subdivide.

    This is the per-box certificate a static (non-traced) subdivision
    can actually prove; the slicing-plane Krawczyk below is kept as a
    fallback for the boxes this misses.
    """
    X = box.intervals()
    gf = f.interval_grad(X)
    gg = g.interval_grad(X)
    n = (
        gf[1] * gg[2] - gf[2] * gg[1],
        gf[2] * gg[0] - gf[0] * gg[2],
        gf[0] * gg[1] - gf[1] * gg[0],
    )
    # Parameter coordinate: tangent component provably nonzero,
    # strongest (largest |mid|) wins.
    k = -1
    best = -1.0
    for i in range(3):
        if n[i].contains_zero():
            continue
        m = abs(n[i].mid())
        if m > best:
            best = m
            k = i
    if k < 0:
        return "none"
    u_idx = [i for i in range(3) if i != k]

    def H2(vals, kk):
        p = [0.0, 0.0, 0.0]
        p[u_idx[0]] = vals[0]
        p[u_idx[1]] = vals[1]
        p[k] = kk
        return np.array([f.float_value(p), g.float_value(p)])

    def J2(vals, kk):
        p = [0.0, 0.0, 0.0]
        p[u_idx[0]] = vals[0]
        p[u_idx[1]] = vals[1]
        p[k] = kk
        gu = f.float_grad(p)
        gv = g.float_grad(p)
        return np.array([[gu[u_idx[0]], gu[u_idx[1]]],
                         [gv[u_idx[0]], gv[u_idx[1]]]])

    k0 = float(box.lo[k] + box.hi[k]) / 2.0
    m = box.center()
    u = np.array([m[u_idx[0]], m[u_idx[1]]])
    for _ in range(5):
        J = J2(u, k0)
        try:
            du = np.linalg.solve(J, H2(u, k0))
        except np.linalg.LinAlgError:
            return "none"
        if not np.all(np.isfinite(du)):
            return "none"
        u = u - du
        if float(np.linalg.norm(du)) < 1e-13:
            break
    J = J2(u, k0)
    try:
        C = np.linalg.inv(J)
    except np.linalg.LinAlgError:
        return "none"
    if not np.all(np.isfinite(C)):
        return "none"

    # Interval residual over the full parameter range.
    Xu = [Interval.point(float(u[0])), Interval.point(float(u[1]))]
    Xfull = [None, None, None]
    Xfull[u_idx[0]] = Xu[0]
    Xfull[u_idx[1]] = Xu[1]
    Xfull[k] = X[k]
    R = [f.interval_eval(tuple(Xfull)), g.interval_eval(tuple(Xfull))]

    # 2x2 interval Jacobian w.r.t. u over the whole 3D box.
    gu = f.interval_grad(X)
    gv = g.interval_grad(X)
    Jint = [[gu[u_idx[0]], gu[u_idx[1]]],
            [gv[u_idx[0]], gv[u_idx[1]]]]

    I2 = _imat_zeros2()
    CJ = _fmat_imat_mul2(C, Jint)
    for i in range(2):
        for j in range(2):
            I2[i][j] = Interval.point(1.0 if i == j else 0.0) - CJ[i][j]
    U = [X[u_idx[0]], X[u_idx[1]]]
    d = [Interval(float(box.lo[u_idx[i]]) - float(u[i]),
                  float(box.hi[u_idx[i]]) - float(u[i]))
         for i in range(2)]
    # Parametric Krawczyk: K(U) = u - C.[]H(u, Kparam)
    #     + (I - C.[]J_u(U, Kparam)) (U - u).
    CR = _fmat_ivec_mul2(C, R)
    K = [Interval.point(float(u[i])) - CR[i] for i in range(2)]
    corr = _imat_ivec_mul2(I2, d)
    K = [K[i] + corr[i] for i in range(2)]
    if all(K[i].lo > float(box.lo[u_idx[i]]) and
           K[i].hi < float(box.hi[u_idx[i]]) for i in range(2)):
        return "unique"
    if any(K[i].hi < float(box.lo[u_idx[i]]) or
           K[i].lo > float(box.hi[u_idx[i]]) for i in range(2)):
        # Every per-parameter zero in U would have to lie in K(U);
        # K(U) misses U entirely, so the box holds no curve at all.
        return "excluded"
    return "none"


def graph_certificate(f, g, box: Box3) -> str:
    """Boolean-style wrapper: 'unique' or 'none' (exclusion folded in)."""
    v = graph_verdict(f, g, box)
    return "unique" if v == "unique" else "none"


# ---------------------------------------------------------------------------
# Exclusion subdivision

@dataclass
class RetainedBox:
    lo: tuple
    hi: tuple
    f_interval: tuple  # (lo, hi), both contain 0 by construction
    g_interval: tuple
    transverse: bool
    krawczyk: str  # "unique" | "exists" | "none"
    unresolved: bool  # kept without a full uniqueness certificate


@dataclass
class EnclosureResult:
    domain: Box3
    retained: list = field(default_factory=list)
    excluded_count: int = 0
    excluded_proofs: list = field(default_factory=list)
    unresolved_count: int = 0
    boxes_visited: int = 0
    min_size: float = 0.0

    def check_partition(self) -> tuple[bool, str]:
        """Certificate check: leaves must partition the domain volume."""
        vol = 0.0
        for r in self.retained:
            lo = np.array(r.lo)
            hi = np.array(r.hi)
            vol += float(np.prod(hi - lo))
        for p in self.excluded_proofs_all:
            lo = np.array(p["lo"])
            hi = np.array(p["hi"])
            vol += float(np.prod(hi - lo))
        d = self.domain
        want = float(np.prod(d.hi - d.lo))
        ok = abs(vol - want) <= 1e-9 * max(1.0, want)
        return ok, f"leaf_volume={vol:.6g} domain_volume={want:.6g}"

    # Filled by enclosure(): every excluded leaf (compact form).
    excluded_proofs_all: list = field(default_factory=list)


def enclosure(f, g, domain: Box3, *, min_size: float,
              max_boxes: int = 400000,
              prove_krawczyk: bool = True) -> EnclosureResult:
    """Rigorous superset enclosure of {F = 0} intersect {G = 0} in domain.

    A box is discarded only with an interval proof that F or G keeps a
    nonzero sign on it. Returns the enclosure plus per-box certificates.
    """
    res = EnclosureResult(domain=domain, min_size=min_size)
    stack = [domain]
    while stack:
        box = stack.pop()
        res.boxes_visited += 1
        if res.boxes_visited > max_boxes:
            # Honest stop: remaining work is unresolved, not certified.
            stack.append(box)
            for b in stack:
                res.retained.append(RetainedBox(
                    lo=tuple(b.lo), hi=tuple(b.hi),
                    f_interval=(float("nan"), float("nan")),
                    g_interval=(float("nan"), float("nan")),
                    transverse=False, krawczyk="none", unresolved=True))
                res.unresolved_count += 1
            break
        X = box.intervals()
        fi = f.interval_eval(X)
        gi = g.interval_eval(X)
        if not fi.contains_zero():
            res.excluded_count += 1
            proof = {"lo": tuple(box.lo), "hi": tuple(box.hi),
                     "surface": "F",
                     "interval": (fi.lo, fi.hi)}
            res.excluded_proofs_all.append(proof)
            if len(res.excluded_proofs) < 8:
                res.excluded_proofs.append(proof)
            continue
        if not gi.contains_zero():
            res.excluded_count += 1
            proof = {"lo": tuple(box.lo), "hi": tuple(box.hi),
                     "surface": "G",
                     "interval": (gi.lo, gi.hi)}
            res.excluded_proofs_all.append(proof)
            if len(res.excluded_proofs) < 8:
                res.excluded_proofs.append(proof)
            continue
        if box.diag() <= min_size:
            # Leaf: attempt the certificates; only boxes that resist
            # both a uniqueness proof and a rigorous exclusion are
            # honestly reported as unresolved.
            trans, _ = transverse_certificate(f, g, box)
            verdict = "none"
            if trans and prove_krawczyk:
                verdict = graph_verdict(f, g, box)
            if verdict == "excluded":
                res.excluded_count += 1
                proof = {"lo": tuple(box.lo), "hi": tuple(box.hi),
                         "surface": "graph-Krawczyk",
                         "interval": None}
                res.excluded_proofs_all.append(proof)
                if len(res.excluded_proofs) < 8:
                    res.excluded_proofs.append(proof)
                continue
            kra = "none"
            if verdict == "unique":
                kra = "unique"
            elif trans and prove_krawczyk:
                kra = krawczyk_prove(f, g, box)  # slicing fallback
            certified = trans and kra == "unique"
            if not certified:
                res.unresolved_count += 1
            res.retained.append(RetainedBox(
                lo=tuple(float(v) for v in box.lo),
                hi=tuple(float(v) for v in box.hi),
                f_interval=(fi.lo, fi.hi), g_interval=(gi.lo, gi.hi),
                transverse=trans, krawczyk=kra,
                unresolved=not certified))
            continue
        trans, _ = transverse_certificate(f, g, box)
        verdict = "none"
        if trans and prove_krawczyk:
            verdict = graph_verdict(f, g, box)
        if verdict == "excluded":
            # Rigorous exclusion: the parametric Krawczyk image misses
            # the box, so no curve point can lie inside it.
            res.excluded_count += 1
            proof = {"lo": tuple(box.lo), "hi": tuple(box.hi),
                     "surface": "graph-Krawczyk",
                     "interval": None}
            res.excluded_proofs_all.append(proof)
            if len(res.excluded_proofs) < 8:
                res.excluded_proofs.append(proof)
            continue
        kra = "none"
        if verdict == "unique":
            kra = "unique"
        elif trans and prove_krawczyk:
            kra = krawczyk_prove(f, g, box)  # slicing fallback
        if trans and kra == "unique":
            # Certified unique regular arc: keep the box as a certified
            # leaf, no need to subdivide further for the enclosure.
            res.retained.append(RetainedBox(
                lo=tuple(float(v) for v in box.lo),
                hi=tuple(float(v) for v in box.hi),
                f_interval=(fi.lo, fi.hi), g_interval=(gi.lo, gi.hi),
                transverse=True, krawczyk="unique", unresolved=False))
            continue
        stack.extend(box.children())
    return res


# ---------------------------------------------------------------------------
# Agreement with OCCT + retained-region clustering

def point_in_enclosure(pt: np.ndarray, res: EnclosureResult,
                       tol: float = 1e-9) -> bool:
    for r in res.retained:
        lo = np.array(r.lo)
        hi = np.array(r.hi)
        if np.all(pt >= lo - tol) and np.all(pt <= hi + tol):
            return True
    return False


def distance_to_enclosure(pt: np.ndarray, res: EnclosureResult) -> float:
    best = float("inf")
    for r in res.retained:
        lo = np.array(r.lo)
        hi = np.array(r.hi)
        d = np.maximum(np.maximum(lo - pt, pt - hi), 0.0)
        best = min(best, float(np.linalg.norm(d)))
    return best


def occt_agreement(res: EnclosureResult, points: np.ndarray,
                   tol: float = 1e-9) -> dict:
    """Every OCCT section sample must lie in the retained enclosure.

    Returns coverage stats. covered < 1.0 is a soundness alarm.
    """
    pts = np.asarray(points, dtype=float).reshape(-1, 3)
    if len(pts) == 0:
        return {"n_points": 0, "covered": 0, "coverage": 1.0,
                "max_gap": 0.0}
    covered = 0
    max_gap = 0.0
    for p in pts:
        if point_in_enclosure(p, res, tol):
            covered += 1
        else:
            max_gap = max(max_gap, distance_to_enclosure(p, res))
    return {"n_points": int(len(pts)), "covered": int(covered),
            "coverage": covered / len(pts), "max_gap": float(max_gap)}


def retained_clusters(res: EnclosureResult) -> list[list[int]]:
    """Connected components of retained boxes on a min_size voxel grid.

    26-connectivity: conservative, avoids splitting one curve into
    pieces over voxel corners.
    """
    s = res.min_size
    if s <= 0:
        return [[i] for i in range(len(res.retained))]
    d = res.domain
    origin = d.lo
    cells: dict[tuple[int, int, int], list[int]] = {}
    for i, r in enumerate(res.retained):
        lo = np.array(r.lo)
        hi = np.array(r.hi)
        c0 = np.floor((lo - origin) / s).astype(int)
        c1 = np.floor((hi - origin) / s).astype(int)
        for ix in range(c0[0], c1[0] + 1):
            for iy in range(c0[1], c1[1] + 1):
                for iz in range(c0[2], c1[2] + 1):
                    cells.setdefault((ix, iy, iz), []).append(i)
    parent = list(range(len(res.retained)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for (ix, iy, iz), members in cells.items():
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    nb = cells.get((ix + dx, iy + dy, iz + dz))
                    if nb:
                        for a in members:
                            for b in nb:
                                union(a, b)
    groups: dict[int, list[int]] = {}
    for i in range(len(res.retained)):
        groups.setdefault(find(i), []).append(i)
    return list(groups.values())


def unexplained_regions(res: EnclosureResult, points: np.ndarray,
                        tol: float = 1e-9) -> list[list[int]]:
    """Retained clusters containing no OCCT sample point.

    These are regions the certified enclosure could not rule out but
    the floating-point tracer did not find: either over-approximation
    or a genuinely missed loop. Reported, never certified either way.
    """
    pts = np.asarray(points, dtype=float).reshape(-1, 3)
    out = []
    for cl in retained_clusters(res):
        explained = False
        for p in pts:
            for i in cl:
                r = res.retained[i]
                lo = np.array(r.lo)
                hi = np.array(r.hi)
                if np.all(p >= lo - tol) and np.all(p <= hi + tol):
                    explained = True
                    break
            if explained:
                break
        if not explained:
            out.append(cl)
    return out
