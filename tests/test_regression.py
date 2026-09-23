"""Regression tests for specific review findings.

Each test reproduces a probe that failed on the pre-fix code:
  1. 1e-9-apart boxes fused into one solid (silent topology change).
  2. Two-corner shape attack: volume preserved, shape corrupted.
  3. Non-round coordinates scored 0.0 engine/classifier agreement.
  4. Sphere certificate understated the true error ~2x.
  5. Flipped triangle winding passed the closure check.
"""

import sys
import numpy as np

sys.path.insert(0, "src")
from brepkernel import boolean, AmbiguousResult
from brepkernel.solids import Box, Sphere, Cylinder
from brepkernel.verify import signed_volume, verify, connected_shells


def vol(mesh):
    return abs(signed_volume(mesh["V"], mesh["F"]))


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return cond


def t1_gap_boxes():
    """Two boxes 1e-9 apart: union must be 2 shells, volume exactly 2."""
    A = Box([0, 0, 0], [1, 1, 1])
    B = Box([1 + 1e-9, 0, 0], [2 + 1e-9, 1, 1])
    try:
        m, r = boolean(A, B, "union")
    except AmbiguousResult as e:
        return check("1e-9 gap boxes", False,
                     f"rejected: {e.report.get('ambiguities')}")
    ok_v = abs(vol(m) - 2.0) < 1e-12
    ok_s = connected_shells(m["F"]) == 2
    return check("1e-9 gap boxes", r["accepted"] and ok_v and ok_s,
                 f"vol={vol(m)!r} shells={connected_shells(m['F'])}")


def _dVdx(V, F, i):
    """Analytic dV/dx_i. V is linear in each vertex x-coordinate."""
    tris = F[np.any(F == i, axis=1)]
    c = 0.0
    for a, b, cc in tris:
        x = V[[a, b, cc]]
        if i == a:
            c += x[1, 1] * x[2, 2] - x[2, 1] * x[1, 2]
        elif i == b:
            c += x[2, 1] * x[0, 2] - x[0, 1] * x[2, 2]
        else:
            c += x[0, 1] * x[1, 2] - x[1, 1] * x[0, 2]
    return c / 6.0


def t2_corner_push():
    """Push two corners +/-0.3 with volume preserved: verify must reject.

    The old checks (closure, even Euler, orientation, MC volume, same-
    engine x-check) all pass on this; only a shape-aware check catches it.
    """
    A = Box([0, 0, 0], [2, 2, 2])
    C = Box([1, 1, 1], [3, 3, 3])
    m, r = boolean(A, C, "union")
    assert r["accepted"]
    V = m["V"].copy()
    F = m["F"].copy()
    v0 = signed_volume(V, F)
    cands = np.where(V[:, 0] == V[:, 0].max())[0]
    pair = None
    for a in range(len(cands)):
        for b in range(a + 1, len(cands)):
            i, j = int(cands[a]), int(cands[b])
            Vp = V.copy()
            Vp[i, 0] += 0.3
            Vp[j, 0] -= 0.3
            if abs(signed_volume(Vp, F) - v0) < 1e-9:
                pair = (i, j, Vp)
                break
        if pair:
            break
    assert pair, "no volume-preserving corner pair found"
    i, j, Vp = pair
    from brepkernel.ingest import ToleranceLedger
    from brepkernel.proxy import certified_proxy
    ledger = ToleranceLedger(proxy_tol=1e-3)
    pA, pB = certified_proxy(A, 1e-3), certified_proxy(C, 1e-3)
    margin = ledger.degeneracy_margin(pA["chordal_error"], pB["chordal_error"])
    assembled = {"V": Vp, "F": F, "empty": False, "margin": margin,
                 "audit": None}
    accepted, vrep = verify(assembled, A, C, "union", pA, pB)
    vv = vrep["checks"]["verts_on_surface"]
    return check("corner-push attack rejected",
                 (not accepted) and (not vv["pass"]),
                 f"verts_on_surface pass={vv['pass']} dV={signed_volume(Vp,F)-v0:.1e}")


def t3_nonround_coords():
    """Fractional coordinates: every face must verify, none ambiguous."""
    A = Box([0.1, -0.2, 0.3], [1.7, 0.9, 2.1])
    B = Box([0.9, 0.4, -0.1], [2.3, 1.8, 1.2])
    try:
        m, r = boolean(A, B, "union")
    except AmbiguousResult as e:
        return check("non-round coords", False,
                     f"rejected: {e.report.get('ambiguities')}")
    asm = r["stages"]["assembly"]
    ok = (r["accepted"] and asm["n_ambiguous"] == 0
          and asm["n_violation"] == 0
          and asm["n_verified"] == asm["n_faces"] > 0)
    return check("non-round coords", ok, f"assembly={asm}")


def t4_sphere_certificate():
    """Measured max deviation must never exceed the claimed bound."""
    rng = np.random.default_rng(7)
    ok = True
    for tol in (1e-2, 1e-3):
        s = Sphere([0.1, -0.2, 0.3], 1.5)
        V, F, bound = s.tessellate(tol)
        N = 200000
        tr = rng.integers(0, len(F), N)
        r1 = np.sqrt(rng.random(N))
        r2 = rng.random(N)
        w = np.stack([1 - r1, r1 * (1 - r2), r1 * r2], 1)
        P = (w[:, 0:1] * V[F[tr, 0]] + w[:, 1:2] * V[F[tr, 1]]
             + w[:, 2:3] * V[F[tr, 2]])
        measured = float(np.max(np.abs(s.implicit(P))))
        good = measured <= bound
        ok &= good
        print(f"    tol={tol}: measured={measured:.3e} "
              f"bound={bound:.3e} ratio={bound/max(measured,1e-300):.2f}")
    return check("sphere certificate honest", ok)


def t5_flipped_winding():
    """One flipped triangle must fail the direction-aware closure check."""
    from brepkernel.verify import check_closure_directed
    A = Box([0, 0, 0], [2, 2, 2])
    C = Box([1, 1, 1], [3, 3, 3])
    m, r = boolean(A, C, "union")
    assert r["accepted"]
    V = m["V"].copy()
    F = m["F"].copy()
    F[0] = F[0][::-1]  # flip winding of one triangle
    ok, info = check_closure_directed(V, F)
    return check("flipped winding rejected", not ok, f"info={info}")


def main():
    ok = True
    ok &= t1_gap_boxes()
    ok &= t2_corner_push()
    ok &= t3_nonround_coords()
    ok &= t4_sphere_certificate()
    ok &= t5_flipped_winding()
    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
