"""Stage 6: independent verification (design G6-G9).

A result is accepted only if ALL checks pass:
  V1 closure      -- every edge used exactly twice (no boundary, no fins)
  V2 euler        -- chi = V - E + F is an even integer (closed surface)
  V3 orientation  -- signed volume > 0 (outward normals; empty result exempt)
  V4 field        -- mesh volume agrees with Monte-Carlo integration of the
                    EXACT implicits (independent of the arrangement engine)
  V5 engine-diverse cross-check -- rebuild via manifold3d and compare volume

Any failure -> the result is rejected, never silently returned.
"""

import math
import numpy as np


def _edges(F):
    e = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    e = np.sort(e, axis=1)
    return e


def check_closure(V, F):
    """Every edge appears exactly twice."""
    if len(F) == 0:
        return True, "empty mesh is vacuously closed"
    e = _edges(F)
    _, counts = np.unique(e, axis=0, return_counts=True)
    bad = np.sum(counts != 2)
    return bad == 0, f"{bad} non-manifold edges of {len(e)//1}"


def euler_chi(V, F):
    """chi = V - E + F for a closed triangle mesh (E = 3F/2)."""
    if len(F) == 0:
        return 0
    return len(V) - len(F) // 2


def signed_volume(V, F):
    """Signed volume; positive iff normals point outward."""
    if len(F) == 0:
        return 0.0
    v0, v1, v2 = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    return float(np.sum(np.einsum("ij,ij->i", v0, np.cross(v1, v2))) / 6.0)


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


def field_consistency(V, F, solidA, solidB, op, n=160000, seed=0):
    """Monte-Carlo volume from EXACT implicits vs mesh volume.

    Independent of the arrangement engine: samples the true solids.
    Passes if |mesh_vol - mc_vol| <= 4*sigma + slack. The wide band and
    large n keep the flake rate negligible -- a verification gate must
    not roll dice on every run.
    """
    rng = np.random.default_rng(seed)
    if len(F) == 0:
        mesh_vol = 0.0
    else:
        mesh_vol = abs(signed_volume(V, F))
    loa, hia = solidA.bbox()
    lob, hib = solidB.bbox()
    lo = np.minimum(loa, lob)
    hi = np.maximum(hia, hib)
    box_vol = float(np.prod(hi - lo))
    if box_vol == 0:
        return True, {"mesh_vol": mesh_vol, "mc_vol": 0.0, "note": "zero box"}
    pts = rng.uniform(lo, hi, size=(n, 3))
    hits = _op_contains(pts, solidA, solidB, op)
    frac = float(np.mean(hits))
    mc_vol = frac * box_vol
    sigma = box_vol * math.sqrt(frac * (1 - frac) / n)
    ok = abs(mesh_vol - mc_vol) <= 4 * sigma + 1e-6
    return ok, {"mesh_vol": mesh_vol, "mc_vol": mc_vol, "sigma": sigma,
                "n": n, "diff": abs(mesh_vol - mc_vol)}


def manifold_crosscheck(V, F):
    """Rebuild through manifold3d and compare volume (diverse engine path)."""
    if len(F) == 0:
        return True, {"note": "empty result"}
    try:
        from manifold3d import Mesh, Manifold
        V = np.ascontiguousarray(V, dtype=np.float32)
        F = np.ascontiguousarray(F, dtype=np.uint32)
        m = Manifold(Mesh(vert_properties=V, tri_verts=F))
        ok = m.status().name == "NoError"
        vol = m.volume()
        return ok, {"volume": vol, "status_ok": ok}
    except Exception as e:
        return False, {"error": str(e)}


def verify(assembled, solidA, solidB, op):
    """Run all checks. Returns (accepted, report dict)."""
    V, F = assembled["V"], assembled["F"]
    report = {"checks": {}}
    if assembled.get("empty"):
        ok, info = field_consistency(V, F, solidA, solidB, op)
        report["checks"]["empty_field"] = {"pass": bool(ok), "info": info}
        report["accepted"] = bool(ok)
        return report["accepted"], report

    c_ok, c_info = check_closure(V, F)
    report["checks"]["closure"] = {"pass": bool(c_ok), "info": c_info}

    chi = euler_chi(V, F)
    chi_ok = (chi % 2 == 0)
    report["checks"]["euler"] = {"pass": bool(chi_ok), "chi": chi}

    vol = signed_volume(V, F)
    report["checks"]["orientation"] = {"pass": bool(vol > 0), "volume": vol}

    f_ok, f_info = field_consistency(V, F, solidA, solidB, op)
    f_info = {k: (float(v) if isinstance(v, (np.floating, float)) else v)
              for k, v in f_info.items()}
    report["checks"]["field"] = {"pass": bool(f_ok), "info": f_info}

    m_ok, m_info = manifold_crosscheck(V, F)
    m_info = {k: (float(v) if isinstance(v, (np.floating, float)) else v)
              for k, v in m_info.items()}
    if "volume" in m_info and vol:
        m_info["rel_volume_diff"] = abs(m_info["volume"] - abs(vol)) / abs(vol)
        m_ok = m_ok and m_info["rel_volume_diff"] < 1e-6
    report["checks"]["manifold_xcheck"] = {"pass": bool(m_ok), "info": m_info}

    report["accepted"] = all(c["pass"] for c in report["checks"].values())
    return report["accepted"], report
