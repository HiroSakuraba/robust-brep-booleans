"""Analytic implicit surfaces with interval and exact evaluation (EXPERIMENTAL, B9).

Each surface stores its parameters twice: as Fractions (for the exact
point classifier in exact_classify.py) and as floats (for interval
evaluation). The implicit function f is a polynomial in (x, y, z) with
f < 0 inside, f = 0 on, f > 0 outside (outward orientation where it
matters).

Surfaces covered: these are exactly the analytic surface kinds the
Tier B/C pipeline ingests from STEP (plane, cylinder, cone, sphere,
torus; see step_ingest.py). Freeform NURBS are out of scope for this
spike: the plan's Bezier-decomposition route would be the follow-up.

Interval evaluation is inclusion-monotone: for a box B,
f_interval(B) is guaranteed to contain { f(p) : p in B }. The gradient
is likewise evaluated over intervals, which the certified intersector
uses for the transversality test (grad F x grad G provably nonzero on
a box means the intersection inside that box is a regular curve with
no singular point).
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

import numpy as np

from .intervals import Interval


def _fr(x) -> Fraction:
    return x if isinstance(x, Fraction) else Fraction(x)


@dataclass(frozen=True)
class ImplicitSurface:
    """Polynomial implicit surface f(x, y, z) = 0.

    terms: tuple of (coeff: Fraction, (i, j, k)) meaning
        coeff * x**i * y**j * z**k.
    grad_terms: for each coordinate, the term list of the partial
        derivative (derived once, exactly, from terms).
    """

    name: str
    terms: tuple
    grad_terms: tuple

    def _eval_terms(self, terms, X, Y, Z):
        total = Interval.point(0.0)
        for coeff, (i, j, k) in terms:
            c = Interval.point(float(coeff))
            t = c * X.pow_int(i) * Y.pow_int(j) * Z.pow_int(k)
            total = total + t
        return total

    def interval_eval(self, box) -> Interval:
        """Rigorous enclosure of f over the box (Interval triple)."""
        X, Y, Z = box
        return self._eval_terms(self.terms, X, Y, Z)

    def interval_grad(self, box) -> tuple[Interval, Interval, Interval]:
        """Rigorous enclosure of grad f over the box."""
        X, Y, Z = box
        return tuple(
            self._eval_terms(g, X, Y, Z) for g in self.grad_terms
        )

    def exact_value(self, p) -> Fraction:
        """Exact f(p) for a rational point p = (Fraction, Fraction, Fraction)."""
        x, y, z = (_fr(v) for v in p)
        total = Fraction(0)
        for coeff, (i, j, k) in self.terms:
            total += coeff * x**i * y**j * z**k
        return total

    def exact_sign(self, p) -> int:
        v = self.exact_value(p)
        return (v > 0) - (v < 0)

    def float_value(self, p) -> float:
        x, y, z = (float(v) for v in p)
        total = 0.0
        for coeff, (i, j, k) in self.terms:
            total += float(coeff) * x**i * y**j * z**k
        return total

    def float_grad(self, p) -> tuple[float, float, float]:
        x, y, z = (float(v) for v in p)
        out = []
        for g in self.grad_terms:
            total = 0.0
            for coeff, (i, j, k) in g:
                total += float(coeff) * x**i * y**j * z**k
            out.append(total)
        return tuple(out)


def _diff_terms(terms, var: int):
    out = []
    for coeff, exps in terms:
        e = exps[var]
        if e:
            exps2 = list(exps)
            exps2[var] = e - 1
            out.append((coeff * e, tuple(exps2)))
    return tuple(out)


def _make(name, terms) -> ImplicitSurface:
    terms = tuple((_fr(c), tuple(e)) for c, e in terms)
    grad = tuple(_diff_terms(terms, v) for v in range(3))
    return ImplicitSurface(name=name, terms=terms, grad_terms=grad)


def _sub(p, c):
    return tuple(_fr(pi) - _fr(ci) for pi, ci in zip(p, c))


def _dot(a, b):
    return sum(_fr(x) * _fr(y) for x, y in zip(a, b))


def _norm2(a):
    return _dot(a, a)


def plane(point, normal) -> ImplicitSurface:
    """Plane through point with the given (outward) normal: n.(x - p) = 0."""
    n = tuple(_fr(v) for v in normal)
    p = tuple(_fr(v) for v in point)
    d = _dot(n, p)
    terms = [(n[0], (1, 0, 0)), (n[1], (0, 1, 0)), (n[2], (0, 0, 1)),
             (-d, (0, 0, 0))]
    return _make(f"plane(n={normal})", terms)


def sphere(center, radius) -> ImplicitSurface:
    """Sphere |x - c|^2 - r^2 = 0."""
    c = tuple(_fr(v) for v in center)
    r = _fr(radius)
    terms = []
    for i in range(3):
        e2 = [0, 0, 0]
        e2[i] = 2
        terms.append((Fraction(1), tuple(e2)))
        e1 = [0, 0, 0]
        e1[i] = 1
        terms.append((-2 * c[i], tuple(e1)))
    terms.append((c[0] * c[0] + c[1] * c[1] + c[2] * c[2] - r * r,
                  (0, 0, 0)))
    return _make(f"sphere(c={center}, r={radius})", terms)


def _axis_frame(axis):
    """Return the axis as a unit Fraction vector; requires a rational
    axis that is already unit length (the pipeline's analytic cases)."""
    a = tuple(_fr(v) for v in axis)
    n2 = _norm2(a)
    if n2 != 1:
        raise ValueError("axis must be a unit vector for exact mode")
    return a


def cylinder(point, axis, radius) -> ImplicitSurface:
    """Circular cylinder |x-p|^2 - ((x-p).a)^2 - r^2 = 0."""
    p = tuple(_fr(v) for v in point)
    a = _axis_frame(axis)
    r = _fr(radius)
    # |x|^2 - 2 p.x + |p|^2 - ((x.a)^2 - 2 (p.a)(x.a) + (p.a)^2) - r^2
    terms = []
    pa = _dot(p, a)
    for i in range(3):
        e2 = [0, 0, 0]
        e2[i] = 2
        terms.append((Fraction(1), tuple(e2)))
        e1 = [0, 0, 0]
        e1[i] = 1
        terms.append((-2 * p[i], tuple(e1)))
    # -(x.a)^2 = -sum_ij a_i a_j x_i x_j
    for i in range(3):
        for j in range(3):
            if a[i] * a[j] == 0:
                continue
            e = [0, 0, 0]
            e[i] += 1
            e[j] += 1
            terms.append((-a[i] * a[j], tuple(e)))
    # +2 (p.a)(x.a)
    for i in range(3):
        if 2 * pa * a[i] == 0:
            continue
        e = [0, 0, 0]
        e[i] = 1
        terms.append((2 * pa * a[i], tuple(e)))
    const = _norm2(p) - pa * pa - r * r
    terms.append((const, (0, 0, 0)))
    return _make(f"cylinder(p={point}, r={radius})", terms)


def cone(apex, axis, cos_half_angle) -> ImplicitSurface:
    """Double cone ((x-p).a)^2 - cos^2(h) |x-p|^2 = 0.

    The caller restricts the domain box to one nappe when needed."""
    p = tuple(_fr(v) for v in apex)
    a = _axis_frame(axis)
    ch = _fr(cos_half_angle)
    ch2 = ch * ch
    pa = _dot(p, a)
    terms = []
    # (x.a)^2 - 2(p.a)(x.a) + (p.a)^2 - ch^2(|x|^2 - 2 p.x + |p|^2)
    for i in range(3):
        for j in range(3):
            if a[i] * a[j] == 0:
                continue
            e = [0, 0, 0]
            e[i] += 1
            e[j] += 1
            terms.append((a[i] * a[j], tuple(e)))
    for i in range(3):
        if 2 * pa * a[i] == 0:
            continue
        e = [0, 0, 0]
        e[i] = 1
        terms.append((-2 * pa * a[i], tuple(e)))
    for i in range(3):
        e2 = [0, 0, 0]
        e2[i] = 2
        terms.append((-ch2, tuple(e2)))
        e1 = [0, 0, 0]
        e1[i] = 1
        terms.append((2 * ch2 * p[i], tuple(e1)))
    const = pa * pa - ch2 * _norm2(p)
    terms.append((const, (0, 0, 0)))
    return _make(f"cone(apex={apex})", terms)


def torus(center, axis, major_radius, minor_radius) -> ImplicitSurface:
    """Torus (|x-c|^2 - R^2 - r^2)^2 - 4 R^2 (|x-c|^2 - ((x-c).a)^2) = 0."""
    c = tuple(_fr(v) for v in center)
    a = _axis_frame(axis)
    R = _fr(major_radius)
    r = _fr(minor_radius)
    # Let u = x - c. s = |u|^2, t = u.a.
    # (sqrt(s - t^2) - R)^2 + t^2 = r^2  <=>  2 R sqrt(s - t^2) = s + R^2 - r^2
    # f = (s + R^2 - r^2)^2 - 4 R^2 (s - t^2)
    # Expand in x directly via helper: write polynomial in u, then shift.
    # Polynomial in u: coefficients keyed by exponent triple.
    from collections import defaultdict
    poly = defaultdict(Fraction)
    K = R * R - r * r
    # (s + K)^2 = s^2 + 2K s + K^2; then - 4R^2 (s - t^2)
    for i in range(3):
        e = [0, 0, 0]
        e[i] = 2
        poly[tuple(e)] += 2 * K  # +2K s term
        poly[tuple(e)] += -4 * R * R  # -4R^2 s term
    for i in range(3):
        for j in range(3):
            e = [0, 0, 0]
            e[i] += 2
            e[j] += 2
            poly[tuple(e)] += Fraction(1)  # s^2 term
    # +4 R^2 t^2 = 4R^2 sum_ij a_i a_j u_i u_j
    for i in range(3):
        for j in range(3):
            if a[i] * a[j] == 0:
                continue
            e = [0, 0, 0]
            e[i] += 1
            e[j] += 1
            poly[tuple(e)] += 4 * R * R * a[i] * a[j]
    poly[(0, 0, 0)] += K * K
    # Shift u = x - c: expand each monomial prod_i (x_i - c_i)^e_i.
    from math import comb
    terms = []
    for (i, j, k), coeff in poly.items():
        if coeff == 0:
            continue
        # expand (x-cx)^i (y-cy)^j (z-cz)^k
        for a1 in range(i + 1):
            for b1 in range(j + 1):
                for d1 in range(k + 1):
                    cc = (coeff * comb(i, a1) * comb(j, b1)
                          * comb(k, d1)
                          * ((-c[0]) ** (i - a1))
                          * ((-c[1]) ** (j - b1))
                          * ((-c[2]) ** (k - d1)))
                    if cc:
                        terms.append((cc, (a1, b1, d1)))
    return _make(f"torus(c={center}, R={major_radius}, r={minor_radius})",
                 terms)
