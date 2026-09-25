"""Thin test helper wrapping the G0 independent membership arbiter.

Tier B/C test scripts that ACCEPT a Boolean result call check_accepted()
with the two input shapes, the result shape, and the op. The arbiter judges
the result wrong at a point only when the generalized winding number of the
result disagrees with the set operation of the winding numbers of the
inputs, which is independent of OCCT's Boolean/intersector machinery.

Tests that expect a refusal must NOT call this helper.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_HERE)
for _p in (_REPO_ROOT, os.path.join(_REPO_ROOT, "tools", "review_probes")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np

# Note: `arbiter` (tools/review_probes/arbiter.py) is imported lazily inside
# raw_audit so this module stays importable even before sys.path is fixed.


# Pure memoization caches. Triangle soups and winding numbers are pure
# functions of (shape, deflection) and (soup, point); caching them only
# avoids recomputation when several audits in one test process share inputs
# (e.g. the same A/B pair audited for union, intersection and difference).
# Shapes/soups are held by strong references so id() keys stay valid.
_tri_cache: dict = {}
_wind_cache: dict = {}


def raw_audit(a, b, out, op, n=300, seed=20260925, lo=None, hi=None):
    """Run the membership audit with pure memoization; no assertions.

    Returns (result_dict, seconds). Triangle soups and winding numbers are
    cached per test process; the cached values are bit-identical to fresh
    computation.
    """
    from time import perf_counter
    import arbiter as arb

    real_triangles, real_winding = arb.triangles, arb.winding

    def cached_triangles(shape, deflection=2e-4):
        key = (id(shape), float(deflection))
        hit = _tri_cache.get(key)
        if hit is None or hit[0] is not shape:
            soup = real_triangles(shape, deflection)
            _tri_cache[key] = (shape, soup)
            return soup
        return hit[1]

    def cached_winding(tris, q):
        key = (id(tris),
               np.ascontiguousarray(q, dtype=np.float64).tobytes())
        hit = _wind_cache.get(key)
        if hit is None or hit[0] is not tris:
            w = real_winding(tris, q)
            _wind_cache[key] = (tris, w)
            return w
        return hit[1]

    orig_membership = arb.membership_audit

    def patched_membership(aa, bb, oo, opp, rng, n=300, lo=-2.2, hi=2.2,
                           band=2e-3, deflection=2e-4):
        arb.triangles, arb.winding = cached_triangles, cached_winding
        try:
            return orig_membership(aa, bb, oo, opp, rng, n=n, lo=lo,
                                   hi=hi, band=band, deflection=deflection)
        finally:
            arb.triangles, arb.winding = real_triangles, real_winding

    arb.membership_audit = patched_membership
    try:
        if lo is None or hi is None:
            lo, hi = audit_bounds(a, b)
        rng = np.random.default_rng(seed)
        t0 = perf_counter()
        res = arb.membership_audit(a, b, out, op, rng, n=n, lo=lo, hi=hi)
        dt = perf_counter() - t0
    finally:
        arb.membership_audit = orig_membership
    return res, dt


def check_accepted(name, check_fn, a, b, out, op, n=300, seed=20260925,
                   lo=None, hi=None):
    """Audit an accepted Boolean result; feed the verdict into check_fn.

    Asserts (through the test file's own check() reporter) that the arbiter
    found zero kernel errors over at least 200 sample points. Returns
    (ok, raw_audit_dict, seconds).
    """
    res, dt = raw_audit(a, b, out, op, n=n, seed=seed, lo=lo, hi=hi)
    detail = (f"checked={res['checked']} "
              f"kernel_errors={res['kernel_errors']} "
              f"audit_ms={dt * 1000.0:.0f}")
    ok = check_fn(
        f"{name} independent membership audit",
        res["kernel_errors"] == 0 and res["checked"] >= 200,
        detail)
    return bool(ok), res, dt


def audit_bounds(a, b, pad_frac=0.1):
    """Sampling box covering both input shapes, returned as (lo, hi)."""
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib
    box = Bnd_Box()
    BRepBndLib.Add_s(a, box)
    BRepBndLib.Add_s(b, box)
    lo = np.array(box.CornerMin().Coord(), dtype=float)
    hi = np.array(box.CornerMax().Coord(), dtype=float)
    pad = pad_frac * float(np.max(hi - lo))
    return lo - pad, hi + pad
