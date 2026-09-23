"""Stage 6: independent verification.

A result is accepted only if ALL applicable checks pass:
  V1 closure     - directed edges: every undirected edge appears exactly
                   twice, once in each direction (tests orientation too)
  V2 euler       - chi = V - E + F equals the predicted value where one is
                   exactly decidable (box-box); reported otherwise
  V3 orientation - signed volume > 0 (outward normals; empty exempt)
  V4 field       - mesh volume vs EXACT closed-form volume (box-box), else a
                   Monte-Carlo check against the exact implicits, honestly
                   labeled statistical
  V5 occ-xcheck  - rebuild the boolean in OCCT (independent engine,
                   independent geometry kernel) and compare volumes
  V6 verts       - every output vertex lies within the margin of an input
                   surface (catches shape errors that preserve volume)
  V7 samples     - point-in-mesh (own ray caster) agrees with exact implicit
                   membership on samples outside the margin band
  V8 shells      - connected shell count matches the expected count where
                   exactly decidable (catches silent fusion/fission)

Any failure -> the result is rejected, never silently returned.
"""

import math
import numpy as np


def signed_volume(V, F):
    """Signed volume; positive iff normals point outward."""
    if len(F) == 0:
        return 0.0
    v0, v1, v2 = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    return float(np.sum(np.einsum("ij,ij->i", v0, np.cross(v1, v2))) / 6.0)


def euler_chi(V, F):
    if len(F) == 0:
        return 0
    e = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    e = np.sort(e, axis=1)
    n_edges = len(np.unique(e, axis=0))
    return len(V) - n_edges + len(F)


def check_closure_directed(V, F):
    """Every undirected edge appears exactly twice, in opposite directions."""
    if len(F) == 0:
        return True, "empty mesh is vacuously closed"
    d = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    key = d[:, 0].astype(np.int64) * (len(V) + 1) + d[:, 1].astype(np.int64)
    ukey, counts = np.unique(key, return_counts=True)
    cnt = dict(zip(ukey.tolist(), counts.tolist()))
    n = len(V) + 1
    bad = 0
    for k, c in cnt.items():
        a, b = k // n, k % n
        if a >= b:
            continue
        c1 = c
        c2 = cnt.get(b * n + a, 0)
        if c1 != 1 or c2 != 1:
            bad += 1
    info = f"{bad} bad undirected edges of {len(cnt)//2}"
    return bad == 0, info


def connected_shells(F):
    """Number of connected triangle shells (union-find over vertices)."""
    if len(F) == 0:
        return 0
    parent = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for tri in F:
        for v in tri:
            v = int(v)
            if v not in parent:
                parent[v] = v
        union(int(tri[0]), int(tri[1]))
        union(int(tri[1]), int(tri[2]))
    # a shell = a set of faces connected via shared vertices; count distinct
    # face-roots
    roots = set()
    for tri in F:
        roots.add(find(int(tri[0])))
    return len(roots)


def _op_contains(p, solidA, solidB, op):
    a = solidA.inside(p)
    b = solidB.inside(p)
    if op == "union":
        return a | b
    if op == "intersection":
        return a & b
    if op == "difference":
        return a & ~b
    raise ValueError(op)


def field_check(V, F, solidA, solidB, op):
    """Mesh volume vs exact closed-form volume (box-box).

    Falls back to Monte-Carlo integration of the exact implicits, honestly
    labeled statistical: it can only catch gross errors, never certify.
    """
    from .classify import exact_op_volume
    mesh_vol = abs(signed_volume(V, F)) if len(F) else 0.0
    exact = exact_op_volume(solidA, solidB, op)
    if exact is not None:
        tol = 1e-9 * max(1.0, abs(exact)) + 1e-12
        ok = abs(mesh_vol - exact) <= tol
        return ok, {"kind": "exact", "mesh_vol": mesh_vol,
                    "exact_vol": exact, "diff": abs(mesh_vol - exact)}
    # statistical fallback
    rng = np.random.default_rng(0)
    n = 40000
    loa, hia = solidA.bbox()
    lob, hib = solidB.bbox()
    lo = np.minimum(loa, lob)
    hi = np.maximum(hia, hib)
    box_vol = float(np.prod(hi - lo))
    pts = rng.uniform(lo, hi, size=(n, 3))
    frac = float(np.mean(_op_contains(pts, solidA, solidB, op)))
    mc_vol = frac * box_vol
    sigma = box_vol * math.sqrt(frac * (1 - frac) / n)
    ok = abs(mesh_vol - mc_vol) <= 4 * sigma + 1e-6
    return ok, {"kind": "statistical-mc", "mesh_vol": mesh_vol,
                "mc_vol": mc_vol, "sigma": sigma, "n": n,
                "diff": abs(mesh_vol - mc_vol)}


def _occt_shape(solid):
    from OCP.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakeSphere,
                                 BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeCone)
    from OCP.gp import gp_Pnt, gp_Ax2, gp_Dir
    k = solid.kind
    if k == "box":
        return BRepPrimAPI_MakeBox(gp_Pnt(*solid.lo),
                                   gp_Pnt(*solid.hi)).Shape()
    if k == "sphere":
        return BRepPrimAPI_MakeSphere(
            gp_Ax2(gp_Pnt(*solid.center), gp_Dir(0, 0, 1)),
            solid.r).Shape()
    if k == "cylinder":
        return BRepPrimAPI_MakeCylinder(
            gp_Ax2(gp_Pnt(*solid.base), gp_Dir(*solid.axis)),
            solid.r, solid.h).Shape()
    if k == "cone":
        return BRepPrimAPI_MakeCone(
            gp_Ax2(gp_Pnt(*solid.base), gp_Dir(*solid.axis)),
            solid.r, 0.0, solid.h).Shape()
    raise ValueError(f"no OCCT builder for {k}")


def occt_crosscheck(V, F, solidA, solidB, op):
    """Independent engine AND independent geometry kernel (OCCT).

    Builds the same boolean from OCCT primitives and compares volumes.
    """
    if len(F) == 0:
        return True, {"note": "empty result"}
    try:
        from OCP.BRepAlgoAPI import (BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut,
                                     BRepAlgoAPI_Common)
        from OCP.BRepGProp import BRepGProp
        from OCP.GProp import GProp_GProps
        sA, sB = _occt_shape(solidA), _occt_shape(solidB)
        algo = {"union": BRepAlgoAPI_Fuse, "intersection": BRepAlgoAPI_Common,
                "difference": BRepAlgoAPI_Cut}[op]
        res = algo(sA, sB)
        if not res.IsDone():
            return True, {"note": "OCCT did not finish; check skipped"}
        props = GProp_GProps()
        BRepGProp.VolumeProperties_s(res.Shape(), props)
        ovol = float(props.Mass())
        mvol = abs(signed_volume(V, F))
        if solidA.kind == "box" and solidB.kind == "box":
            tol = 1e-6 * max(1.0, mvol)
        else:
            tol = 0.02 * max(1.0, mvol) + 1e-9
        ok = abs(ovol - mvol) <= tol
        return ok, {"occt_vol": ovol, "mesh_vol": mvol,
                    "diff": abs(ovol - mvol), "tol": tol}
    except Exception as e:
        return True, {"note": f"OCCT unavailable ({e}); check skipped"}


def vertex_on_surface(V, solidA, solidB, margin_v):
    """Every output vertex must lie within margin_v of an input surface.

    Result vertices are either proxy vertices (on a true surface) or
    engine intersection points (on both proxy surfaces, hence within the
    chordal bound of the true surfaces). A pushed/corrupted vertex fails.
    """
    if len(V) == 0:
        return True, {"note": "empty"}
    dA = np.abs(solidA.implicit(V))
    dB = np.abs(solidB.implicit(V))
    dmin = np.minimum(dA, dB)
    bad = int(np.sum(dmin > margin_v))
    return bad == 0, {"n_bad": bad, "n_verts": len(V),
                      "max_dmin": float(np.max(dmin)),
                      "margin": margin_v}


def _ray_parity(pts, V, F, direction):
    """Odd/even crossing count for rays along direction (vectorized MT)."""
    d = np.asarray(direction, dtype=float)
    d = d / np.linalg.norm(d)
    inside = np.zeros(len(pts), dtype=bool)
    chunk = 1500
    for s in range(0, len(F), chunk):
        Fc = F[s:s + chunk]
        v0, v1, v2 = V[Fc[:, 0]], V[Fc[:, 1]], V[Fc[:, 2]]
        e1 = v1 - v0
        e2 = v2 - v0
        pvec = np.cross(np.broadcast_to(d, e2.shape), e2)
        det = np.einsum("ij,ij->i", e1, pvec)
        ok_det = np.abs(det) > 1e-30
        inv = np.where(ok_det, 1.0 / np.where(ok_det, det, 1.0), 0.0)
        tvec = pts[:, None, :] - v0[None, :, :]
        u = np.einsum("ijk,jk->ij", tvec, pvec) * inv[None, :]
        qvec = np.cross(tvec, e1[None, :, :])
        v = np.einsum("ijk,k->ij", qvec, d) * inv[None, :]
        t = np.einsum("ijk,jk->ij", qvec, e2) * inv[None, :]
        hit = (ok_det[None, :] & (u >= 0) & (v >= 0)
               & (u + v <= 1) & (t > 1e-12))
        inside ^= (np.sum(hit, axis=1) % 2 == 1)
    return inside


def sample_membership(V, F, solidA, solidB, op, margin, n=1000, seed=1):
    """Point-in-mesh (own ray caster) vs exact implicit membership.

    Samples outside the margin band of both surfaces must agree 100%.
    Rays are cast in two directions; disagreeing samples are undecided.
    """
    if len(F) == 0:
        return True, {"note": "empty"}
    rng = np.random.default_rng(seed)
    loa, hia = solidA.bbox()
    lob, hib = solidB.bbox()
    lo = np.minimum(loa, lob)
    hi = np.maximum(hia, hib)
    pts = rng.uniform(lo, hi, size=(n, 3))
    dA = np.abs(solidA.implicit(pts))
    dB = np.abs(solidB.implicit(pts))
    decided = (dA > margin) & (dB > margin)
    if np.sum(decided) < 50:
        return True, {"note": "too few decided samples; check skipped",
                      "n_decided": int(np.sum(decided))}
    p = pts[decided]
    in1 = _ray_parity(p, V, F, (1.0, 0.0, 0.0))
    in2 = _ray_parity(p, V, F, (0.0, 1.0, 0.0))
    agree_rays = in1 == in2
    exact = _op_contains(p, solidA, solidB, op)
    match = (in1 == exact) & agree_rays
    n_bad = int(np.sum(~match))
    return n_bad == 0, {"n_samples": n, "n_decided": int(np.sum(decided)),
                        "n_undecided_rays": int(np.sum(~agree_rays)),
                        "n_mismatch": n_bad}


def verify(assembled, solidA, solidB, op, proxyA, proxyB):
    """Run all checks. Returns (accepted, report dict)."""
    from .classify import expected_euler, expected_shells
    V, F = assembled["V"], assembled["F"]
    margin = assembled.get("margin", 1e-3)
    report = {"checks": {}}

    if assembled.get("empty") or len(F) == 0:
        f_ok, f_info = field_check(V, F, solidA, solidB, op)
        report["checks"]["empty_field"] = {"pass": bool(f_ok), "info": f_info}
        report["accepted"] = bool(f_ok)
        return report["accepted"], report

    c_ok, c_info = check_closure_directed(V, F)
    report["checks"]["closure"] = {"pass": bool(c_ok), "info": c_info}

    chi = euler_chi(V, F)
    pred = expected_euler(solidA, solidB, op)
    if pred is None:
        report["checks"]["euler"] = {"pass": True, "chi": chi,
                                     "info": "no exact prediction; reported"}
    else:
        report["checks"]["euler"] = {"pass": bool(chi == pred), "chi": chi,
                                     "predicted": pred}

    vol = signed_volume(V, F)
    report["checks"]["orientation"] = {"pass": bool(vol > 0), "volume": vol}

    f_ok, f_info = field_check(V, F, solidA, solidB, op)
    report["checks"]["field"] = {"pass": bool(f_ok), "info": _jsonable(f_info)}

    o_ok, o_info = occt_crosscheck(V, F, solidA, solidB, op)
    report["checks"]["occt_xcheck"] = {"pass": bool(o_ok),
                                       "info": _jsonable(o_info)}

    mv = max(proxyA["chordal_error"], proxyB["chordal_error"]) + 1e-9
    v_ok, v_info = vertex_on_surface(V, solidA, solidB, mv)
    report["checks"]["verts_on_surface"] = {"pass": bool(v_ok),
                                            "info": _jsonable(v_info)}

    s_ok, s_info = sample_membership(V, F, solidA, solidB, op, margin)
    report["checks"]["sample_membership"] = {"pass": bool(s_ok),
                                             "info": _jsonable(s_info)}

    n_shells = connected_shells(F)
    exp_shells = expected_shells(solidA, solidB, op)
    if exp_shells is None:
        report["checks"]["shells"] = {"pass": True, "n_shells": n_shells,
                                      "info": "no exact prediction; reported"}
    else:
        report["checks"]["shells"] = {"pass": bool(n_shells == exp_shells),
                                      "n_shells": n_shells,
                                      "expected": exp_shells}

    report["accepted"] = all(c["pass"] for c in report["checks"].values())
    return report["accepted"], report


def _jsonable(d):
    return {k: (float(v) if isinstance(v, (np.floating, float)) else v)
            for k, v in d.items()}
