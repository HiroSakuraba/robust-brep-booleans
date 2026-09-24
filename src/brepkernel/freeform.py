"""Freeform / NURBS acceleration helpers for the future Tier B/C path.

Design goals:
- preserve the current certify-or-refuse contract;
- make broad-phase work cheap by indexing NURBS knot spans with conservative
  control-hull AABBs;
- make point-to-surface verification fast with cached local span seeds and a
  Gauss-Newton projector;
- never silently treat an untrimmed-surface hit as a trimmed-face hit;
- fall back to OCCT's trimmed-face distance for hard trim/seam cases.

This module is optional: OCP/cadquery-ocp is imported lazily. The pure
NurbsSurface evaluator/index can also be used without OCCT.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional

import numpy as np


class FreeformError(Exception):
    """Typed freeform-path refusal / extraction error."""

    def __init__(self, message: str, kind: str = "FreeformError"):
        super().__init__(message)
        self.kind = kind


def _expanded_knots(unique: np.ndarray, mult: np.ndarray) -> np.ndarray:
    out = np.repeat(np.asarray(unique, dtype=np.float64),
                    np.asarray(mult, dtype=np.int64))
    if len(out) < 2 or np.any(np.diff(out) < 0):
        raise FreeformError("invalid knot sequence", "InvalidKnots")
    return out


def _find_span(n: int, p: int, u: float, U: np.ndarray) -> int:
    """The NURBS Book A2.1. n is highest control-point index."""
    if u >= U[n + 1]:
        return n
    if u <= U[p]:
        return p
    lo, hi = p, n + 1
    mid = (lo + hi) // 2
    while u < U[mid] or u >= U[mid + 1]:
        if u < U[mid]:
            hi = mid
        else:
            lo = mid
        mid = (lo + hi) // 2
    return mid


def _basis_funs(span: int, u: float, p: int, U: np.ndarray) -> np.ndarray:
    """Non-zero basis values N_{span-p..span,p}(u)."""
    N = np.zeros(p + 1, dtype=np.float64)
    left = np.zeros(p + 1, dtype=np.float64)
    right = np.zeros(p + 1, dtype=np.float64)
    N[0] = 1.0
    for j in range(1, p + 1):
        left[j] = u - U[span + 1 - j]
        right[j] = U[span + j] - u
        saved = 0.0
        for r in range(j):
            den = right[r + 1] + left[j - r]
            temp = 0.0 if abs(den) < 1e-300 else N[r] / den
            N[r] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        N[j] = saved
    return N


def _basis_and_derivative(span: int, u: float, p: int,
                          U: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Basis and first derivative for the p+1 active basis functions."""
    N = _basis_funs(span, u, p, U)
    dN = np.zeros(p + 1, dtype=np.float64)
    if p == 0:
        return N, dN
    Nm1 = _basis_funs(span, u, p - 1, U)
    first_p = span - p
    first_m = span - (p - 1)
    for local in range(p + 1):
        i = first_p + local
        a = 0.0
        b = 0.0
        j = i - first_m
        if 0 <= j < p:
            den = U[i + p] - U[i]
            if abs(den) > 1e-300:
                a = p * Nm1[j] / den
        j2 = (i + 1) - first_m
        if 0 <= j2 < p:
            den = U[i + p + 1] - U[i + 1]
            if abs(den) > 1e-300:
                b = p * Nm1[j2] / den
        dN[local] = a - b
    return N, dN


@dataclass(frozen=True)
class NurbsSurface:
    degree_u: int
    degree_v: int
    knots_u: np.ndarray
    knots_v: np.ndarray
    poles: np.ndarray
    weights: np.ndarray
    uv_bounds: tuple[float, float, float, float]
    periodic_u: bool = False
    periodic_v: bool = False

    def __post_init__(self):
        P = np.asarray(self.poles, dtype=np.float64)
        W = np.asarray(self.weights, dtype=np.float64)
        U = np.asarray(self.knots_u, dtype=np.float64)
        V = np.asarray(self.knots_v, dtype=np.float64)
        if P.ndim != 3 or P.shape[2] != 3:
            raise FreeformError("poles must have shape (nu,nv,3)",
                                "InvalidControlNet")
        if W.shape != P.shape[:2]:
            raise FreeformError("weight/control-net shape mismatch",
                                "InvalidControlNet")
        if not np.all(np.isfinite(P)) or not np.all(np.isfinite(W)):
            raise FreeformError("non-finite NURBS data", "InvalidControlNet")
        if len(U) != P.shape[0] + self.degree_u + 1:
            raise FreeformError("U knot count does not match degree/poles",
                                "InvalidKnots")
        if len(V) != P.shape[1] + self.degree_v + 1:
            raise FreeformError("V knot count does not match degree/poles",
                                "InvalidKnots")
        if np.any(np.diff(U) < 0) or np.any(np.diff(V) < 0):
            raise FreeformError("knots must be non-decreasing", "InvalidKnots")
        object.__setattr__(self, "poles", P)
        object.__setattr__(self, "weights", W)
        object.__setattr__(self, "knots_u", U)
        object.__setattr__(self, "knots_v", V)

    @property
    def has_positive_weights(self) -> bool:
        return bool(np.all(self.weights > 0.0))

    @property
    def domain(self) -> tuple[float, float, float, float]:
        nu, nv = self.poles.shape[:2]
        return (float(self.knots_u[self.degree_u]),
                float(self.knots_u[nu]),
                float(self.knots_v[self.degree_v]),
                float(self.knots_v[nv]))

    def clamp_uv(self, u: float, v: float) -> tuple[float, float]:
        u0, u1, v0, v1 = self.uv_bounds
        return (float(np.clip(u, u0, u1)), float(np.clip(v, v0, v1)))

    def eval_d1(self, u: float, v: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Evaluate S, dS/du, dS/dv using local-support rational basis."""
        u, v = self.clamp_uv(u, v)
        nu, nv = self.poles.shape[:2]
        su = _find_span(nu - 1, self.degree_u, u, self.knots_u)
        sv = _find_span(nv - 1, self.degree_v, v, self.knots_v)
        Nu, dNu = _basis_and_derivative(su, u, self.degree_u, self.knots_u)
        Nv, dNv = _basis_and_derivative(sv, v, self.degree_v, self.knots_v)
        iu = np.arange(su - self.degree_u, su + 1)
        iv = np.arange(sv - self.degree_v, sv + 1)
        P = self.poles[np.ix_(iu, iv)]
        W = self.weights[np.ix_(iu, iv)]

        B = Nu[:, None] * Nv[None, :] * W
        Bu = dNu[:, None] * Nv[None, :] * W
        Bv = Nu[:, None] * dNv[None, :] * W
        den = float(np.sum(B))
        if abs(den) < 1e-300:
            raise FreeformError("rational denominator vanished",
                                "DegenerateWeights")
        A = np.einsum("ij,ijk->k", B, P)
        Au = np.einsum("ij,ijk->k", Bu, P)
        Av = np.einsum("ij,ijk->k", Bv, P)
        Wu = float(np.sum(Bu))
        Wv = float(np.sum(Bv))
        S = A / den
        Su = (Au - S * Wu) / den
        Sv = (Av - S * Wv) / den
        return S, Su, Sv

    def value(self, u: float, v: float) -> np.ndarray:
        return self.eval_d1(u, v)[0]

    def nonzero_spans(self) -> tuple[list[tuple[float, float]],
                                     list[tuple[float, float]]]:
        """Unique non-zero knot spans clipped to the face UV rectangle."""
        u0, u1, v0, v1 = self.uv_bounds
        us = []
        for a, b in zip(self.knots_u[:-1], self.knots_u[1:]):
            a2, b2 = max(float(a), u0), min(float(b), u1)
            if b2 > a2 and (not us or us[-1] != (a2, b2)):
                us.append((a2, b2))
        vs = []
        for a, b in zip(self.knots_v[:-1], self.knots_v[1:]):
            a2, b2 = max(float(a), v0), min(float(b), v1)
            if b2 > a2 and (not vs or vs[-1] != (a2, b2)):
                vs.append((a2, b2))
        return us, vs

    def active_control_block(self, u: float, v: float) -> np.ndarray:
        nu, nv = self.poles.shape[:2]
        su = _find_span(nu - 1, self.degree_u, u, self.knots_u)
        sv = _find_span(nv - 1, self.degree_v, v, self.knots_v)
        iu = np.arange(su - self.degree_u, su + 1)
        iv = np.arange(sv - self.degree_v, sv + 1)
        return self.poles[np.ix_(iu, iv)]


@dataclass(frozen=True)
class PatchRecord:
    u0: float
    u1: float
    v0: float
    v1: float
    lo: np.ndarray
    hi: np.ndarray

    @property
    def uv_center(self) -> tuple[float, float]:
        return ((self.u0 + self.u1) * 0.5, (self.v0 + self.v1) * 0.5)


class NurbsPatchIndex:
    """Conservative knot-span broad phase for a positive-weight NURBS face."""

    def __init__(self, surface: NurbsSurface):
        self.surface = surface
        self.records: list[PatchRecord] = []
        safe_local = (surface.has_positive_weights and
                      not surface.periodic_u and not surface.periodic_v)
        if safe_local:
            us, vs = surface.nonzero_spans()
            for ua, ub in us:
                for va, vb in vs:
                    um = 0.5 * (ua + ub)
                    vm = 0.5 * (va + vb)
                    cp = surface.active_control_block(um, vm).reshape(-1, 3)
                    self.records.append(PatchRecord(
                        ua, ub, va, vb, cp.min(axis=0), cp.max(axis=0)))
        else:
            cp = surface.poles.reshape(-1, 3)
            u0, u1, v0, v1 = surface.uv_bounds
            self.records.append(PatchRecord(
                u0, u1, v0, v1, cp.min(axis=0), cp.max(axis=0)))
        self.lo = (np.vstack([r.lo for r in self.records])
                   if self.records else np.zeros((0, 3)))
        self.hi = (np.vstack([r.hi for r in self.records])
                   if self.records else np.zeros((0, 3)))

    def query_aabb(self, lo: Iterable[float], hi: Iterable[float],
                   pad: float = 0.0) -> np.ndarray:
        lo = np.asarray(lo, dtype=np.float64) - pad
        hi = np.asarray(hi, dtype=np.float64) + pad
        if not len(self.records):
            return np.zeros(0, dtype=np.int64)
        hit = (np.all(self.lo <= hi[None, :], axis=1)
               & np.all(self.hi >= lo[None, :], axis=1))
        return np.nonzero(hit)[0]

    def nearest_seed_indices(self, point: Iterable[float],
                             k: int = 6) -> np.ndarray:
        p = np.asarray(point, dtype=np.float64)
        if not len(self.records):
            return np.zeros(0, dtype=np.int64)
        q = np.maximum(np.maximum(self.lo - p, p - self.hi), 0.0)
        d2 = np.einsum("ij,ij->i", q, q)
        k = min(max(int(k), 1), len(d2))
        if k == len(d2):
            return np.argsort(d2)
        idx = np.argpartition(d2, k - 1)[:k]
        return idx[np.argsort(d2[idx])]

    def seed_uvs(self, point: Iterable[float],
                 k: int = 6) -> list[tuple[float, float]]:
        return [self.records[int(i)].uv_center
                for i in self.nearest_seed_indices(point, k)]


@dataclass
class ProjectionResult:
    distance: float
    point: np.ndarray
    uv: tuple[float, float]
    converged: bool
    trimmed_in: Optional[bool] = None
    used_exact_trim_fallback: bool = False


class FreeformFaceAccel:
    """Cached NURBS evaluator + patch broad phase for one OCCT face."""

    def __init__(self, face, nurbs: NurbsSurface, trim_tol: float = 1e-8):
        self.face = face
        self.nurbs = nurbs
        self.index = NurbsPatchIndex(nurbs)
        self.trim_tol = float(trim_tol)

    @classmethod
    def from_occt_face(cls, face, trim_tol: float = 1e-8):
        return cls(face, nurbs_from_occt_face(face), trim_tol=trim_tol)

    def _trim_contains(self, u: float, v: float) -> bool:
        from OCP.BRepClass import BRepClass_FaceClassifier
        from OCP.gp import gp_Pnt2d
        from OCP.TopAbs import TopAbs_IN, TopAbs_ON
        c = BRepClass_FaceClassifier(self.face,
                                     gp_Pnt2d(float(u), float(v)),
                                     self.trim_tol, True)
        return c.State() in (TopAbs_IN, TopAbs_ON)

    def _exact_trimmed_distance(self, point: np.ndarray) -> float:
        from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
        from OCP.BRepExtrema import BRepExtrema_DistShapeShape
        from OCP.gp import gp_Pnt
        v = BRepBuilderAPI_MakeVertex(
            gp_Pnt(*map(float, point))).Vertex()
        d = BRepExtrema_DistShapeShape(v, self.face)
        if not d.IsDone():
            d.Perform()
        if not d.IsDone():
            raise FreeformError("OCCT trimmed-face distance failed",
                                "ProjectionFailed")
        return float(d.Value())

    def project_point(self, point: Iterable[float], *, seeds: int = 6,
                      max_iter: int = 15, xyz_tol: float = 1e-10,
                      uv_tol: float = 1e-12,
                      exact_trim_fallback: bool = True) -> ProjectionResult:
        """Fast local projection, with exact trimmed-face fallback."""
        x = np.asarray(point, dtype=np.float64)
        best = None
        u0, u1, v0, v1 = self.nurbs.uv_bounds
        seed_uvs = self.index.seed_uvs(x, seeds)
        if not seed_uvs:
            seed_uvs = [((u0 + u1) * 0.5, (v0 + v1) * 0.5)]
        seed_uvs += [((u0 + u1) * 0.5, (v0 + v1) * 0.5),
                     (u0, v0), (u0, v1), (u1, v0), (u1, v1)]
        seen = set()
        for us, vs in seed_uvs:
            key = (round(float(us), 15), round(float(vs), 15))
            if key in seen:
                continue
            seen.add(key)
            u, v = float(us), float(vs)
            converged = False
            for _ in range(max_iter):
                S, Su, Sv = self.nurbs.eval_d1(u, v)
                r = S - x
                J = np.column_stack([Su, Sv])
                H = J.T @ J
                g = J.T @ r
                lam = max(1e-18, 1e-12 * float(np.trace(H) + 1.0))
                try:
                    step = np.linalg.solve(H + lam * np.eye(2), g)
                except np.linalg.LinAlgError:
                    break
                un = float(np.clip(u - step[0], u0, u1))
                vn = float(np.clip(v - step[1], v0, v1))
                if np.linalg.norm([un - u, vn - v]) <= uv_tol:
                    u, v = un, vn
                    converged = True
                    break
                u, v = un, vn
            S, _, _ = self.nurbs.eval_d1(u, v)
            dist = float(np.linalg.norm(S - x))
            if dist <= xyz_tol:
                converged = True
            cand = (dist, S, (u, v), converged)
            if best is None or dist < best[0]:
                best = cand

        assert best is not None
        dist, S, uv, converged = best
        trimmed_in = self._trim_contains(*uv)
        if trimmed_in:
            return ProjectionResult(dist, S, uv, converged, True, False)
        if not exact_trim_fallback:
            return ProjectionResult(dist, S, uv, converged, False, False)
        dtrim = self._exact_trimmed_distance(x)
        return ProjectionResult(dtrim, S, uv, True, False, True)

    def refinement_plan(self, target_sag: float, *,
                        max_edge: float | None = None) -> list[dict]:
        """Rank knot spans for adaptive meshing using sampled curvature.

        This is a scheduling hint, not a certificate. It may spend MORE
        triangles on difficult patches but never relaxes final verification.
        """
        if not target_sag > 0:
            raise ValueError("target_sag must be positive")
        from OCP.BRepAdaptor import BRepAdaptor_Surface
        from OCP.GeomLProp import GeomLProp_SLProps

        bs = BRepAdaptor_Surface(self.face).BSpline()
        rows = []
        for i, rec in enumerate(self.index.records):
            u, v = rec.uv_center
            props = GeomLProp_SLProps(bs, float(u), float(v), 2, 1e-12)
            if props.IsCurvatureDefined():
                kappa = max(abs(float(props.MaxCurvature())),
                            abs(float(props.MinCurvature())))
            else:
                kappa = float("inf")
            if not np.isfinite(kappa):
                edge = 0.0
            elif kappa <= 1e-15:
                edge = float("inf")
            else:
                edge = float(np.sqrt(8.0 * target_sag / kappa))
            if max_edge is not None:
                edge = min(edge, float(max_edge))
            diag = float(np.linalg.norm(rec.hi - rec.lo))
            priority = (diag / max(edge, 1e-15)
                        if np.isfinite(edge) else 0.0)
            rows.append({"patch": i,
                         "uv": (rec.u0, rec.u1, rec.v0, rec.v1),
                         "kappa_sample": kappa,
                         "recommended_edge": edge,
                         "bbox_diag": diag,
                         "priority": priority})
        rows.sort(key=lambda r: -r["priority"])
        return rows

    def verify_points(self, points: np.ndarray, tolerance: float,
                      *, exact_trim_fallback: bool = True) -> dict:
        """Verify result samples against the original trimmed face."""
        P = np.asarray(points, dtype=np.float64).reshape(-1, 3)
        rows = []
        worst = 0.0
        n_fallback = 0
        for i, p in enumerate(P):
            r = self.project_point(
                p, exact_trim_fallback=exact_trim_fallback)
            worst = max(worst, r.distance)
            n_fallback += int(r.used_exact_trim_fallback)
            rows.append({"index": i, "distance": r.distance,
                         "uv": r.uv, "trimmed_in": r.trimmed_in,
                         "fallback": r.used_exact_trim_fallback,
                         "pass": r.distance <= tolerance})
        return {"pass": all(r["pass"] for r in rows),
                "n": len(rows), "max_distance": worst,
                "n_exact_trim_fallback": n_fallback,
                "rows": rows}


def _shift_knots_to_face_window(knots: np.ndarray, degree: int,
                                n_poles: int, face0: float, face1: float,
                                period: float, axis: str) -> np.ndarray:
    """Translate a clamped periodic copy by a whole number of periods.

    The geometric surface is unchanged by a whole-period translation, while
    keeping the knot frame aligned with the original face is essential because
    its p-curves remain expressed in that original UV frame.
    """
    K = np.asarray(knots, dtype=np.float64).copy()
    d0 = float(K[degree])
    d1 = float(K[n_poles])
    scale = max(abs(face0), abs(face1), abs(d0), abs(d1), abs(period), 1.0)
    eps = 1e-10 * scale
    if not period > 0.0 or not np.isfinite(period):
        raise FreeformError(f"invalid {axis}-period {period}",
                            "InvalidPeriodicSurface")
    if face1 - face0 > period + eps:
        raise FreeformError(
            f"{axis}-face window spans more than one period; "
            "periodic acceleration refused",
            "PeriodicWindowTooWide")

    center_delta = 0.5 * ((face0 + face1) - (d0 + d1))
    k0 = int(round(center_delta / period))
    for k in (k0, k0 - 1, k0 + 1, k0 - 2, k0 + 2):
        s0, s1 = d0 + k * period, d1 + k * period
        if face0 >= s0 - eps and face1 <= s1 + eps:
            return K + k * period
    raise FreeformError(
        f"could not align clamped {axis}-periodic knot domain "
        f"[{d0}, {d1}] to face window [{face0}, {face1}]",
        "PeriodicWindowUncovered")


def _verify_nurbs_parameter_frame(face, nurbs: NurbsSurface,
                                  rel_tol: float = 5e-11) -> None:
    """Refuse an accelerator whose UV frame disagrees with OCCT.

    This is especially important after de-periodizing a B-spline: geometric
    equivalence alone is insufficient because trim p-curves use the original
    face parameters.
    """
    from OCP.BRepAdaptor import BRepAdaptor_Surface

    ad = BRepAdaptor_Surface(face)
    u0, u1, v0, v1 = nurbs.uv_bounds
    fractions = (0.07, 0.21, 0.43, 0.67, 0.89)
    worst = 0.0
    scale = 1.0
    for fu in fractions:
        u = u0 + (u1 - u0) * fu
        for fv in fractions:
            v = v0 + (v1 - v0) * fv
            po = ad.Value(float(u), float(v))
            a = np.array([po.X(), po.Y(), po.Z()], dtype=np.float64)
            b = nurbs.value(float(u), float(v))
            worst = max(worst, float(np.linalg.norm(a - b)))
            scale = max(scale, float(np.linalg.norm(a)), float(np.linalg.norm(b)))
    tol = rel_tol * scale
    if worst > tol:
        raise FreeformError(
            f"NURBS UV-frame verification failed: max error {worst:.6g} "
            f"> {tol:.6g}",
            "NurbsParameterFrameMismatch")


def nurbs_from_occt_face(face) -> NurbsSurface:
    """Extract transformed OCCT BSpline face into NumPy NURBS data."""
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.BRepTools import BRepTools
    from OCP.GeomAbs import GeomAbs_BSplineSurface

    ad = BRepAdaptor_Surface(face)
    if ad.GetType() != GeomAbs_BSplineSurface:
        raise FreeformError(
            f"face is {ad.GetType()}, not a BSpline surface",
            "UnsupportedSurfaceType")
    source = ad.BSpline()
    was_periodic_u = bool(source.IsUPeriodic())
    was_periodic_v = bool(source.IsVPeriodic())
    period_u = float(source.UPeriod()) if was_periodic_u else 0.0
    period_v = float(source.VPeriod()) if was_periodic_v else 0.0

    # Work on a copy: de-periodizing changes the pole/knot representation and
    # must never mutate the imported B-rep that owns the authoritative
    # p-curves.
    bs = source.Copy()
    if was_periodic_u:
        bs.SetUNotPeriodic()
    if was_periodic_v:
        bs.SetVNotPeriodic()

    pu, pv = int(bs.UDegree()), int(bs.VDegree())
    nu, nv = int(bs.NbUPoles()), int(bs.NbVPoles())
    poles = np.empty((nu, nv, 3), dtype=np.float64)
    weights = np.ones((nu, nv), dtype=np.float64)
    for i in range(1, nu + 1):
        for j in range(1, nv + 1):
            p = bs.Pole(i, j)
            poles[i - 1, j - 1] = (p.X(), p.Y(), p.Z())
            weights[i - 1, j - 1] = float(bs.Weight(i, j))
    ku = np.array([bs.UKnot(i)
                   for i in range(1, bs.NbUKnots() + 1)],
                  dtype=np.float64)
    kv = np.array([bs.VKnot(i)
                   for i in range(1, bs.NbVKnots() + 1)],
                  dtype=np.float64)
    mu = np.array([bs.UMultiplicity(i)
                   for i in range(1, bs.NbUKnots() + 1)],
                  dtype=np.int64)
    mv = np.array([bs.VMultiplicity(i)
                   for i in range(1, bs.NbVKnots() + 1)],
                  dtype=np.int64)
    U = _expanded_knots(ku, mu)
    V = _expanded_knots(kv, mv)
    bounds = tuple(float(x) for x in BRepTools.UVBounds_s(face))
    u0, u1, v0, v1 = bounds

    # Set*NotPeriodic preserves the surface but clamps its knot frame.  A face
    # produced by previous modeling operations may address the same periodic
    # surface one or more whole periods away, so translate the copied knot
    # vector back into the face's parameter frame before constructing the
    # NumPy evaluator.
    if was_periodic_u:
        U = _shift_knots_to_face_window(
            U, pu, nu, u0, u1, period_u, "U")
    if was_periodic_v:
        V = _shift_knots_to_face_window(
            V, pv, nv, v0, v1, period_v, "V")

    nurbs = NurbsSurface(
        pu, pv, U, V, poles, weights, bounds, False, False)
    _verify_nurbs_parameter_frame(face, nurbs)
    return nurbs


def candidate_patch_pairs(a: NurbsPatchIndex, b: NurbsPatchIndex,
                          pad: float = 0.0) -> list[tuple[int, int]]:
    """Sweep-and-prune conservative NURBS span-pair broad phase."""
    if not len(a.records) or not len(b.records):
        return []
    alo, ahi = a.lo - pad, a.hi + pad
    blo, bhi = b.lo - pad, b.hi + pad
    oa = np.argsort(alo[:, 0], kind="mergesort")
    ob = np.argsort(blo[:, 0], kind="mergesort")
    pairs = []
    active: list[int] = []
    jb = 0
    for ia in oa:
        xmin, xmax = alo[ia, 0], ahi[ia, 0]
        while jb < len(ob) and blo[ob[jb], 0] <= xmax:
            active.append(int(ob[jb]))
            jb += 1
        active = [j for j in active
                  if bhi[j, 0] >= xmin and blo[j, 0] <= xmax]
        for j in active:
            if (alo[ia, 1] <= bhi[j, 1]
                    and ahi[ia, 1] >= blo[j, 1]
                    and alo[ia, 2] <= bhi[j, 2]
                    and ahi[ia, 2] >= blo[j, 2]):
                pairs.append((int(ia), int(j)))
    return pairs
