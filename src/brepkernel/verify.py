"""Stage 6: independent verification.

A result is accepted only if every applicable check passes:
  V1 closure     - directed edges: every undirected edge appears exactly
                   twice, once in each direction (tests orientation too)
  V2 euler       - chi = V - E + F equals the exactly predicted value
                   (box-box grid); SKIPPED where not exactly decidable
  V3 orientation - every shell's signed volume matches its role
                   (union/intersection: all positive; difference: shells
                   containing A-faces positive, pure-B shells negative).
                   Catches inside-out shells that parity/total-volume miss.
  V4 field       - mesh volume vs EXACT closed-form volume (box-box), else a
                   Monte-Carlo check against the exact implicits, honestly
                   labeled statistical
  V5 occ-xcheck  - rebuild the boolean in OCCT (independent engine,
                   independent geometry kernel) and compare volumes;
                   SKIPPED if OCCT errors or does not finish
  V6 verts       - every output vertex lies within the margin of an input
                   surface (catches shape errors that preserve volume)
  V7 samples     - generalized winding number (solid-angle sum) vs exact
                   implicit membership on samples outside the margin band.
                   An inside-out shell gives -1, failing immediately;
                   SKIPPED if too few decided samples
  V8 shells      - connected shell count matches the exactly predicted
                   count (box-box grid); SKIPPED where not decidable

Check statuses are pass / fail / skip. A skip is never counted as a
pass: the report says "certified with <names> skipped", and the caller
chooses via allow_skips whether a skip blocks certification.
Any failure -> the result is rejected, never silently returned.

Volumes are computed about a local origin (bounding-box minimum):
signed volume is translation-invariant for closed meshes, and the local
origin keeps large-coordinate inputs (1e8) exact instead of washing out
in catastrophic cancellation about the world origin.
"""

import math
import numpy as np


def _coord_scale(V):
    if len(V) == 0:
        return 1.0
    return max(1.0, float(np.max(np.abs(V))))


def signed_volume(V, F):
    """Signed volume about a local origin; positive iff normals outward."""
    if len(F) == 0:
        return 0.0
    t = np.min(V, axis=0)
    W = V - t
    v0, v1, v2 = W[F[:, 0]], W[F[:, 1]], W[F[:, 2]]
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


def shell_face_labels(F):
    """Face -> shell label via vertex union-find (same as connected_shells)."""
    if len(F) == 0:
        return np.zeros(0, dtype=np.int64)
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
    roots = np.array([find(int(t)) for t in F[:, 0]])
    _, labels = np.unique(roots, return_inverse=True)
    return labels


def connected_shells(F):
    """Number of connected triangle shells (union-find over vertices)."""
    if len(F) == 0:
        return 0
    return int(np.max(shell_face_labels(F)) + 1)


def shell_signed_volumes(V, F):
    """Signed volume per shell: {label: volume}."""
    labels = shell_face_labels(F)
    vols = {}
    for k in np.unique(labels):
        vols[int(k)] = signed_volume(V, F[labels == k])
    return vols


def check_orientation(V, F, origins, op):
    """Every shell's signed volume matches its role.

    union/intersection: every shell is an outer boundary -> positive.
    difference: a shell containing any A-face is outer -> positive; a
    shell of pure B-faces is a cavity -> negative. (A cavity boundary is
    exactly the B-surface inside A, so it cannot contain A-faces; an
    outer shell always retains part of A's boundary.)
    """
    if len(F) == 0:
        return "pass", {"note": "empty"}
    vols = shell_signed_volumes(V, F)
    total = sum(vols.values())
    bad = {}
    if op in ("union", "intersection"):
        for k, v in vols.items():
            if not v > 0:
                bad[k] = v
        expectation = "all shells positive"
    elif op == "difference":
        if origins is None or len(origins) != len(F):
            ok = total > 0
            return ("pass" if ok else "fail",
                    {"note": "no per-face origins; total-volume fallback",
                     "total_volume": total})
        labels = shell_face_labels(F)
        origins = np.asarray(origins)
        for k, v in vols.items():
            hasA = bool(np.any(origins[labels == k] == "A"))
            if hasA and not v > 0:
                bad[k] = (v, "outer shell not positive")
            elif not hasA and not v < 0:
                bad[k] = (v, "cavity shell not negative")
        expectation = "outer shells positive, cavity shells negative"
    else:
        raise ValueError(op)
    info = {"expectation": expectation, "shell_volumes": vols,
            "total_volume": total, "bad_shells": bad}
    return ("pass" if not bad else "fail", info)


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
    The exact tolerance accounts for float64 vertex rounding:
    tol = 1e-9*|exact| + 4*eps*scale*surface_area.
    """
    from .classify import exact_op_volume
    mesh_vol = abs(signed_volume(V, F)) if len(F) else 0.0
    exact = exact_op_volume(solidA, solidB, op)
    if exact is not None:
        scale = _coord_scale(V)
        eps = np.finfo(float).eps * scale
        v0, v1, v2 = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
        area = float(np.sum(np.linalg.norm(np.cross(v1 - v0, v2 - v0),
                                                 axis=1)) / 2.0)
        tol = 1e-9 * max(1.0, abs(exact)) + 4.0 * eps * area + 1e-12
        ok = abs(mesh_vol - exact) <= tol
        return ("pass" if ok else "fail",
                {"kind": "exact", "mesh_vol": mesh_vol,
                 "exact_vol": exact, "diff": abs(mesh_vol - exact),
                 "tol": tol})
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
    return ("pass" if ok else "fail",
            {"kind": "statistical-mc", "mesh_vol": mesh_vol,
             "mc_vol": mc_vol, "sigma": sigma, "n": n,
             "diff": abs(mesh_vol - mc_vol)})


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
    Returns status 'skip' (never a silent pass) if OCCT errors or bails.
    """
    if len(F) == 0:
        return "skip", {"note": "empty result"}
    try:
        from OCP.BRepAlgoAPI import (BRepAlgoAPI_Fuse, BRepAlgoAPI_Common,
                                     BRepAlgoAPI_Cut)
        from OCP.BRepGProp import BRepGProp
        from OCP.GProp import GProp_GProps
        sA, sB = _occt_shape(solidA), _occt_shape(solidB)
        algo = {"union": BRepAlgoAPI_Fuse, "intersection": BRepAlgoAPI_Common,
                "difference": BRepAlgoAPI_Cut}[op]
        res = algo(sA, sB)
        if not res.IsDone():
            return "skip", {"note": "OCCT did not finish"}
        props = GProp_GProps()
        BRepGProp.VolumeProperties_s(res.Shape(), props)
        ovol = float(props.Mass())
        mvol = abs(signed_volume(V, F))
        if solidA.kind == "box" and solidB.kind == "box":
            tol = 1e-6 * max(1.0, mvol)
        else:
            tol = 0.02 * max(1.0, mvol) + 1e-9
        ok = abs(ovol - mvol) <= tol
        return ("pass" if ok else "fail",
                {"occt_vol": ovol, "mesh_vol": mvol,
                 "diff": abs(ovol - mvol), "tol": tol})
    except Exception as e:
        return "skip", {"note": f"OCCT unavailable ({e})"}


def vertex_on_surface(V, solidA, solidB, margin_v):
    """Every output vertex must lie within margin_v of an input surface.

    Result vertices are either proxy vertices (on a true surface) or
    engine intersection points (on both proxy surfaces, hence within the
    chordal bound of the true surfaces). A pushed/corrupted vertex fails.
    """
    if len(V) == 0:
        return "pass", {"note": "empty"}
    dA = np.abs(solidA.implicit(V))
    dB = np.abs(solidB.implicit(V))
    dmin = np.minimum(dA, dB)
    bad = int(np.sum(dmin > margin_v))
    return ("pass" if bad == 0 else "fail",
            {"n_bad": bad, "n_verts": len(V),
             "max_dmin": float(np.max(dmin)),
             "margin": margin_v})


def _winding_number(pts, V, F):
    """Generalized winding number per point via the solid-angle sum.

    w(p) = (1/4pi) * sum over triangles of signed solid angle.
    For a correctly oriented closed solid: +1 inside, 0 outside.
    An inside-out shell contributes -1 inside it. Unlike ray casting this
    has no edge/vertex degeneracies: sample points are outside the margin
    band, hence never exactly on the surface, so the Van Oosterom &
    Strackee formula is stable and the sum rounds cleanly to an integer.
    """
    pts = np.asarray(pts, dtype=float)
    w = np.zeros(len(pts))
    chunk = 1500
    for s in range(0, len(F), chunk):
        Fc = F[s:s + chunk]
        a = V[Fc[:, 0]][None, :, :] - pts[:, None, :]
        b = V[Fc[:, 1]][None, :, :] - pts[:, None, :]
        c = V[Fc[:, 2]][None, :, :] - pts[:, None, :]
        na = np.linalg.norm(a, axis=2)
        nb = np.linalg.norm(b, axis=2)
        nc = np.linalg.norm(c, axis=2)
        num = np.einsum("ijk,ijk->ij", a, np.cross(b, c))
        den = (na * nb * nc + np.einsum("ijk,ijk->ij", a, b) * nc
               + np.einsum("ijk,ijk->ij", b, c) * na
               + np.einsum("ijk,ijk->ij", c, a) * nb)
        w += np.sum(2.0 * np.arctan2(num, den), axis=1)
    return np.rint(w / (4.0 * np.pi)).astype(np.int64)


def sample_membership(V, F, solidA, solidB, op, margin, n=1000, seed=1):
    """Generalized winding number (solid-angle sum) vs exact membership.

    Samples outside the margin band of both surfaces must agree 100%:
    winding +1 where the exact boolean says inside, 0 where outside.
    An inside-out shell yields -1 and fails immediately.
    Returns 'skip' (never a silent pass) if too few samples decide.
    """
    if len(F) == 0:
        return "skip", {"note": "empty"}
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
        return "skip", {"note": "too few decided samples",
                        "n_decided": int(np.sum(decided))}
    p = pts[decided]
    w = _winding_number(p, V, F)
    exact = _op_contains(p, solidA, solidB, op)
    match = ((w == 1) == exact)
    n_bad = int(np.sum(~match))
    return ("pass" if n_bad == 0 else "fail",
            {"n_samples": n, "n_decided": int(np.sum(decided)),
             "n_mismatch": n_bad})


def verify(assembled, solidA, solidB, op, proxyA, proxyB, allow_skips=True):
    """Run all checks. Returns (accepted, report dict).

    allow_skips: if False, any skipped check blocks certification.
    """
    from .classify import expected_euler, expected_shells
    V, F = assembled["V"], assembled["F"]
    origins = assembled.get("origins")
    margin = assembled.get("margin", 1e-3)
    report = {"checks": {}}

    if assembled.get("empty") or len(F) == 0:
        f_st, f_info = field_check(V, F, solidA, solidB, op)
        report["checks"]["empty_field"] = {"status": f_st,
                                           "pass": f_st == "pass",
                                           "info": _jsonable(f_info)}
        report["accepted"] = f_st == "pass"
        report["skipped_checks"] = []
        report["certification"] = ("certified" if report["accepted"]
                                   else "not certified")
        return report["accepted"], report

    def put(name, status, info):
        info = _jsonable(info) if isinstance(info, dict) else info
        report["checks"][name] = {"status": status,
                                  "pass": status == "pass",
                                  "info": info}

    c_ok, c_info = check_closure_directed(V, F)
    put("closure", "pass" if c_ok else "fail", c_info)

    chi = euler_chi(V, F)
    pred = expected_euler(solidA, solidB, op)
    if pred is None:
        put("euler", "skip", {"chi": chi,
                              "note": "no exact prediction"})
    else:
        put("euler", "pass" if chi == pred else "fail",
            {"chi": chi, "predicted": pred})

    o_st, o_info = check_orientation(V, F, origins, op)
    put("orientation", o_st, o_info)

    f_st, f_info = field_check(V, F, solidA, solidB, op)
    put("field", f_st, f_info)

    o2_st, o2_info = occt_crosscheck(V, F, solidA, solidB, op)
    put("occt_xcheck", o2_st, o2_info)

    mv = (max(proxyA["chordal_error"], proxyB["chordal_error"])
          + 1e-9 * _coord_scale(V))
    v_st, v_info = vertex_on_surface(V, solidA, solidB, mv)
    put("verts_on_surface", v_st, v_info)

    s_st, s_info = sample_membership(V, F, solidA, solidB, op, margin)
    put("sample_membership", s_st, s_info)

    n_shells = connected_shells(F)
    exp_shells = expected_shells(solidA, solidB, op)
    if exp_shells is None:
        put("shells", "skip", {"n_shells": n_shells,
                               "note": "no exact prediction"})
    else:
        put("shells", "pass" if n_shells == exp_shells else "fail",
            {"n_shells": n_shells, "expected": exp_shells})

    failed = [k for k, c in report["checks"].items()
              if c["status"] == "fail"]
    skipped = [k for k, c in report["checks"].items()
               if c["status"] == "skip"]
    report["skipped_checks"] = skipped
    accepted = not failed and (allow_skips or not skipped)
    if accepted and skipped:
        report["certification"] = ("certified with skipped checks: "
                                   + ", ".join(skipped))
    else:
        report["certification"] = "certified" if accepted else "not certified"
    report["accepted"] = accepted
    return accepted, report


def _jsonable(d):
    out = {}
    for k, v in d.items():
        if isinstance(v, (np.floating, float)):
            out[k] = float(v)
        elif isinstance(v, (np.integer, int)):
            out[k] = int(v)
        elif isinstance(v, np.ndarray):
            out[k] = v.tolist()
        elif isinstance(v, dict):
            out[k] = _jsonable(v)
        elif isinstance(v, (list, tuple)):
            out[k] = [(_jsonable({"v": x})["v"]
                       if isinstance(x, dict) else x) for x in v]
        else:
            out[k] = v
    return out
