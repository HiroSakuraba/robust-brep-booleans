"""Regression tests for specific review findings.

Each test reproduces a probe that failed on the pre-fix code:
  1. 1e-9-apart boxes fused into one solid (silent topology change).
  2. Two-corner shape attack: volume preserved, shape corrupted.
  3. Non-round coordinates scored 0.0 engine/classifier agreement.
  4. Sphere certificate understated the true error ~2x.
  5. Flipped triangle winding passed the closure check.
Second-review findings (23 Sept 2026):
  6. Inside-out shell passed all eight Stage 6 checks.
  7. Sphere bump through a face (r=0.1) verified 12/12 faces.
  8. Slab-split difference: predictor said 1 shell / chi=2 (answer: 2, 4).
  9. Tunnel difference: predictor said chi=2 (answer: chi=0).
  10. Boxes at 1e8 offset refused by the absolute 1e-9 margin.
Third-review probe (23 Sept 2026):
  11. Hairline bridge silently accepted as 1 shell at proxy_tol=3e-3.
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
    """Two boxes 1e-8 apart: union must be 2 shells, volume exactly 2.

    (The margin for boxes is 1e-9; the old 1e-9 gap sat exactly ON the
    margin and passed only by a float-rounding hair. 1e-8 is cleanly
    above it; the at/below-margin behavior is pinned by stress test A.)
    """
    A = Box([0, 0, 0], [1, 1, 1])
    B = Box([1 + 1e-8, 0, 0], [2 + 1e-8, 1, 1])
    try:
        m, r = boolean(A, B, "union")
    except AmbiguousResult as e:
        return check("1e-9 gap boxes", False,
                     f"rejected: {e.report.get('ambiguities')}")
    ok_v = abs(vol(m) - 2.0) < 1e-12
    ok_s = connected_shells(m["F"]) == 2
    return check("1e-9 gap boxes", r["accepted"] and ok_v and ok_s,
                 f"vol={vol(m)!r} shells={connected_shells(m['F'])}")


def t6_inside_out_shell():
    """A shell with flipped winding must fail winding AND per-shell orientation.

    Second-review attack: disjoint-sphere union with one shell reversed
    passed all eight old checks (directed closure is per-shell
    self-consistent, total-volume orientation is sign-blind, ray parity
    ignores face direction). The new V7 (generalized winding number) and
    V3 (per-shell signed volume vs the shell's role) must both fail.
    """
    from brepkernel.verify import verify
    from brepkernel.proxy import certified_proxy
    from brepkernel.ingest import ToleranceLedger
    A = Box([0, 0, 0], [1, 1, 1])
    B = Box([3, 3, 3], [4, 4, 4])
    m, r = boolean(A, B, "union")
    assert r["accepted"]
    V, F = m["V"], m["F"]
    # disjoint boxes: no cut faces, so centroid-in-bbox gives origins exactly
    c = V[F].mean(axis=1)
    inA = np.all((c >= 0.0) & (c <= 1.0), axis=1)
    origins = np.where(inA, "A", "B")
    assert np.any(~inA) and np.any(inA)
    F2 = F.copy()
    F2[~inA] = F2[~inA][:, ::-1]  # turn shell B inside-out
    ledger = ToleranceLedger(proxy_tol=1e-3)
    pA, pB = certified_proxy(A, 1e-3), certified_proxy(B, 1e-3)
    margin = ledger.degeneracy_margin(pA["chordal_error"], pB["chordal_error"])
    assembled = {"V": V, "F": F2, "origins": origins, "empty": False,
                 "margin": margin}
    accepted, vrep = verify(assembled, A, B, "union", pA, pB)
    w = vrep["checks"]["sample_membership"]["status"]
    o = vrep["checks"]["orientation"]["status"]
    return check("inside-out shell caught",
                 (not accepted) and w == "fail" and o == "fail",
                 f"accepted={accepted} winding={w} orientation={o}")


def _bump_verdict(R):
    """One big A-face, union, with a B-sphere poking through it away from
    the centroid: the centroid is decisively outside B, so only the
    subdivision can prove the engine forgot to cut. Returns the audit."""
    from brepkernel.classify import audit_faces
    V = np.array([[-1., -1., 0.], [1., -1., 0.], [0., 1., 0.]])
    F = np.array([[0, 1, 2]])
    A = Box([-2, -2, -2], [2, 2, 2])
    B = Sphere([0.3, 0.2, -0.02], R)
    return audit_faces(F, V, np.array(["A"]), A, B, 2e-9, "union")


def t7_bump_through_face():
    """Sphere bumps (r=0.1 and r=0.15) poking through a face away from its
    centroid must be audit VIOLATIONS, not verifications.

    Second-review attack: every face centroid and vertex lay outside the
    sphere, so the centroid audit verified 12/12 faces on the r=0.1 bump.
    """
    ok = True
    for R in (0.1, 0.15):
        a = _bump_verdict(R)
        good = bool(a["violation"][0]) and not bool(a["verified"][0])
        ok &= good
        print(f"    bump r={R}: violation={bool(a['violation'][0])} "
              f"verified={bool(a['verified'][0])}")
    return check("bump through face is a violation", ok)


def t8_slab_split_difference():
    """Slab cutting A in two: 2 shells, chi=4, vol exactly 7.2.

    The old interval-relations predictor said 1 shell / chi=2 here; the
    exact cell-grid predictor gets 2 shells and chi=4.
    """
    A = Box([0, 0, 0], [2, 2, 2])
    B = Box([-1, 0.9, -1], [3, 1.1, 3])
    try:
        m, r = boolean(A, B, "difference")
    except AmbiguousResult as e:
        return check("slab split", False, f"rejected: {e}")
    v = vol(m)
    nsh = connected_shells(m["F"])
    from brepkernel.verify import euler_chi
    chi = euler_chi(m["V"], m["F"])
    ok = (r["accepted"] and nsh == 2 and chi == 4
          and abs(v - 7.2) < 1e-9)
    return check("slab split", ok,
                 f"vol={v!r} shells={nsh} chi={chi}")


def t9_tunnel_difference():
    """Square tunnel through A: 1 shell, chi=0 (genus 1), vol exactly 60."""
    A = Box([0, 0, 0], [4, 4, 4])
    B = Box([1.5, 1.5, -1], [2.5, 2.5, 5])
    try:
        m, r = boolean(A, B, "difference")
    except AmbiguousResult as e:
        return check("tunnel", False, f"rejected: {e}")
    v = vol(m)
    nsh = connected_shells(m["F"])
    from brepkernel.verify import euler_chi
    chi = euler_chi(m["V"], m["F"])
    ok = (r["accepted"] and nsh == 1 and chi == 0
          and abs(v - 60.0) < 1e-9)
    return check("tunnel", ok, f"vol={v!r} shells={nsh} chi={chi}")


def t10_large_offset_union():
    """Boxes offset by 1e8: accepted, volume exactly 15.0.

    The old absolute 1e-9 margin was below float64 coordinate rounding at
    1e8 (~1.5e-8), so this was refused; the scale-aware margin
    (1e-9 * coord_scale) plus local-origin volume computation fix it.
    """
    E = 1e8
    A = Box([E, 0, 0], [E + 2, 2, 2])
    B = Box([E + 1, 1, 1], [E + 3, 3, 3])
    try:
        m, r = boolean(A, B, "union")
    except AmbiguousResult as e:
        return check("1e8 offset union", False, f"rejected: {e}")
    v = vol(m)
    ok = r["accepted"] and abs(v - 15.0) < 1e-6
    return check("1e8 offset union", ok, f"vol={v!r}")


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


def t11_thin_bridge():
    """Hairline bridge refused; puncture lips still accepted.

    Third-review probe: a bar (1.0 thick in x) minus a y-axis cylinder of
    radius 0.5+delta, delta below the proxy chordal error. The true surface
    pierces both x faces (2 pieces); at proxy_tol=3e-3 and delta=0.25*chordal
    the inscribed proxy falls short, leaving hairline bridges. The bridge
    config must raise AmbiguousResult (sub_margin_thin_feature) -- v0.3
    silently accepted 1 shell. A larger delta (0.9*chordal) has no bridge
    and must accept with 2 shells, chi=4.
    """
    from brepkernel.proxy import certified_proxy
    from brepkernel.verify import euler_chi
    A = Box([1.5, 0, -2], [2.5, 1, 3])
    tol = 3e-3
    chordal = certified_proxy(
        Cylinder([2, -1, 0.5], [0, 1, 0], 0.5, 3.0), tol)["chordal_error"]
    ok = True
    B = Cylinder([2, -1, 0.5], [0, 1, 0], 0.5 + 0.25 * chordal, 3.0)
    try:
        m, r = boolean(A, B, "difference", proxy_tol=tol)
        ok &= check("bridge refused", False,
                    f"silently accepted {connected_shells(m['F'])} shells")
    except AmbiguousResult as e:
        types = [a.get("check", a["type"]) for a in e.report["ambiguities"]]
        ok &= check("bridge refused", "sub_margin_thin_feature" in types,
                    f"ambiguities={types[:4]}")
    B2 = Cylinder([2, -1, 0.5], [0, 1, 0], 0.5 + 0.9 * chordal, 3.0)
    try:
        m, r = boolean(A, B2, "difference", proxy_tol=tol)
    except AmbiguousResult as e:
        return ok & check("no-bridge accepted", False, f"rejected: {e}")
    nsh = connected_shells(m["F"])
    chi = euler_chi(m["V"], m["F"])
    ok &= check("no-bridge accepted",
                r["accepted"] and nsh == 2 and chi == 4,
                f"shells={nsh} chi={chi}")
    return ok


def main():
    ok = True
    ok &= t1_gap_boxes()
    ok &= t2_corner_push()
    ok &= t3_nonround_coords()
    ok &= t4_sphere_certificate()
    ok &= t5_flipped_winding()
    ok &= t6_inside_out_shell()
    ok &= t7_bump_through_face()
    ok &= t8_slab_split_difference()
    ok &= t9_tunnel_difference()
    ok &= t10_large_offset_union()
    ok &= t11_thin_bridge()
    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
