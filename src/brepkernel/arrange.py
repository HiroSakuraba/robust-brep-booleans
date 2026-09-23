"""Stage 2: arrangement via an external boolean engine.

The engine computes the geometric arrangement of the two proxy meshes:
where the surfaces intersect and how the pieces connect. Topology verdicts
(keep/discard per face) are NOT trusted from the engine; they are audited
in classify.py against the exact implicits. The engine therefore affects
only intersection-curve geometry, never the final topology decision.

Currently manifold3d Mesh64: full float64, so the engine and the Tier A
float64 checks look at the same geometry. Per-triangle origin is tracked
through the boolean via face_id tags, so every result face is classified
against the *other* solid only.
"""

import numpy as np

ENGINE = "manifold3d-mesh64"


def _to_manifold(proxy, tag0):
    from manifold3d import Mesh64, Manifold
    V = np.ascontiguousarray(np.asarray(proxy["V"], dtype=np.float64))
    F = np.ascontiguousarray(np.asarray(proxy["F"], dtype=np.uint32))
    n = len(F)
    mesh = Mesh64(vert_properties=V, tri_verts=F,
                  face_id=np.arange(n, dtype=np.uint32) + tag0)
    return Manifold(mesh), n


def arrange(proxy_a, proxy_b, op):
    """Boolean the two proxy meshes.

    Returns {'V','F','origins','engine'} with V float64, origins in {'A','B'}
    per result face. Raises ArrangementError on engine failure.
    """
    if op not in ("union", "intersection", "difference"):
        raise ValueError(f"unknown op {op!r}")
    try:
        mA, nA = _to_manifold(proxy_a, 0)
        mB, nB = _to_manifold(proxy_b, nA)
        if op == "union":
            mR = mA + mB
        elif op == "intersection":
            mR = mA ^ mB
        else:
            mR = mA - mB
    except Exception as e:
        raise ArrangementError(f"engine failed on {op}: {e}")
    if mR.status().name != "NoError":
        raise ArrangementError(f"engine returned error status {mR.status()}")
    out = mR.to_mesh64()
    V = np.asarray(out.vert_properties, dtype=np.float64)
    F = np.asarray(out.tri_verts, dtype=np.int64).reshape(-1, 3)
    face_id = np.asarray(out.face_id, dtype=np.int64).reshape(-1)
    if len(face_id) != len(F):
        raise ArrangementError("engine did not return per-face origin ids")
    if len(F) and np.any(face_id >= nA + nB):
        raise ArrangementError("engine returned out-of-range face ids")
    origins = np.where(face_id < nA, "A", "B")
    return {"V": V, "F": F, "origins": origins, "engine": ENGINE,
            "n_faces_a": nA, "n_faces_b": nB,
            "empty": len(F) == 0}


class ArrangementError(Exception):
    pass
