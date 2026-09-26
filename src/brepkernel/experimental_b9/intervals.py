"""Rigorous interval arithmetic over binary64 floats (EXPERIMENTAL, B9 research).

This module implements the interval primitive that the B9 certified
intersector spike builds on. Every operation widens the computed float
result outward by a fixed number of ulps, so the returned interval is
guaranteed to contain the true real result as long as:

- inputs are finite (no infinities or NaNs enter the computation),
- the correctly rounded IEEE-754 result is within 0.5 ulp of the true
  result (round-to-nearest, the default mode), and
- no overflow to infinity occurs (all our domains are bounded and the
  polynomials are low degree, so this is checked by assertion).

The widening margin is deliberately generous (4 ulps each side) rather
than the tightest possible, trading a little enclosure sharpness for a
clearer soundness argument. This is a research prototype, not a
production interval library; see the module README.
"""

from __future__ import annotations

import math

# Number of ulps of outward widening applied to every operation result.
_WIDEN_ULPS = 4


def _widen(lo: float, hi: float) -> tuple[float, float]:
    for _ in range(_WIDEN_ULPS):
        lo = math.nextafter(lo, -math.inf)
        hi = math.nextafter(hi, math.inf)
    return lo, hi


class Interval:
    """Closed interval [lo, hi] with lo <= hi, both finite floats."""

    __slots__ = ("lo", "hi")

    def __init__(self, lo: float, hi: float):
        lo = float(lo)
        hi = float(hi)
        if not (math.isfinite(lo) and math.isfinite(hi)):
            raise ValueError("interval endpoints must be finite")
        if lo > hi:
            raise ValueError("interval lower bound exceeds upper bound")
        self.lo = lo
        self.hi = hi

    @staticmethod
    def point(x: float) -> "Interval":
        """Degenerate interval holding exactly the float x, widened once.

        The widening covers the case where x itself is a rounded
        representation of a true rational value (e.g. float(Fraction)).
        """
        return Interval(*_widen(float(x), float(x)))

    def __repr__(self) -> str:
        return f"Interval({self.lo!r}, {self.hi!r})"

    def width(self) -> float:
        return self.hi - self.lo

    def mid(self) -> float:
        return 0.5 * (self.lo + self.hi)

    def contains_zero(self) -> bool:
        return self.lo <= 0.0 <= self.hi

    def __neg__(self) -> "Interval":
        return Interval(*_widen(-self.hi, -self.lo))

    def __add__(self, other: "Interval") -> "Interval":
        o = _as_interval(other)
        return Interval(*_widen(self.lo + o.lo, self.hi + o.hi))

    def __radd__(self, other: "Interval") -> "Interval":
        return self.__add__(other)

    def __sub__(self, other: "Interval") -> "Interval":
        o = _as_interval(other)
        return Interval(*_widen(self.lo - o.hi, self.hi - o.lo))

    def __rsub__(self, other: "Interval") -> "Interval":
        return _as_interval(other).__sub__(self)

    def __mul__(self, other: "Interval") -> "Interval":
        o = _as_interval(other)
        if self.lo == 0.0 and self.hi == 0.0:
            return Interval(*_widen(0.0, 0.0))
        if o.lo == 0.0 and o.hi == 0.0:
            return Interval(*_widen(0.0, 0.0))
        prods = (
            self.lo * o.lo,
            self.lo * o.hi,
            self.hi * o.lo,
            self.hi * o.hi,
        )
        for p in prods:
            if not math.isfinite(p):
                raise ArithmeticError(
                    "interval multiply overflowed; domain not bounded"
                )
        return Interval(*_widen(min(prods), max(prods)))

    def __rmul__(self, other: "Interval") -> "Interval":
        return self.__mul__(other)

    def __truediv__(self, other: "Interval") -> "Interval":
        o = _as_interval(other)
        if o.contains_zero():
            raise ZeroDivisionError(
                "interval division by an interval containing zero"
            )
        return self * Interval(*_widen(1.0 / o.hi, 1.0 / o.lo))

    def pow_int(self, n: int) -> "Interval":
        """Integer power via repeated squaring (n >= 0)."""
        if n < 0:
            raise ValueError("negative integer power not supported")
        if n == 0:
            return Interval(*_widen(1.0, 1.0))
        result = Interval(*_widen(1.0, 1.0))
        base = self
        while n:
            if n & 1:
                result = result * base
            base = base * base
            n >>= 1
        return result

    def __pow__(self, n: int) -> "Interval":
        return self.pow_int(int(n))


def _as_interval(x) -> Interval:
    if isinstance(x, Interval):
        return x
    return Interval.point(float(x))


def interval_pow_even_corrected(x: Interval, n: int) -> Interval:
    """x**n for even n, with the lower bound clamped to >= 0 when
    appropriate. Used by the implicit-surface evaluators."""
    r = x.pow_int(n)
    if n % 2 == 0 and x.contains_zero():
        return Interval(0.0 if r.lo < 0.0 else r.lo, r.hi)
    return r
