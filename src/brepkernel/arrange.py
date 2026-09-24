"""Stage 2: arrangement via an external boolean engine (thin wrapper).

The engine computes the geometric arrangement of the two proxy meshes.
Topology verdicts (keep/discard per face) are NOT trusted from the engine;
they are audited in classify.py against the exact implicits. The engine
therefore affects only intersection-curve geometry, never the final
topology decision.

Solid-semantics normalization (validate -> nesting tree -> canonicalize)
lives in normalize.py; this module prepares both operands, runs the engine
boolean, and then runs the semantic verification leg (Boolean inequalities
+ adversarial membership witnesses) refuse-fast on the result.

Currently manifold3d Mesh64: full float64, so the engine and the Tier A
float64 checks look at the same geometry. Per-triangle origin is tracked
through the boolean via face_id tags, so every result face is classified
against the *other* solid only.
"""

import numpy as np

from .normalize import (prepare_operand, check_boolean_semantics,
                        ArrangementError)

ENGINE = "manifold3d-mesh64"


def arrange(proxy_a, proxy_b, op):
    """Boolean the two proxy meshes.

    Returns {'V','F','origins','engine'} with V float64, origins in {'A','B'}
    per result face. Raises ArrangementError on invalid input, engine
    failure, normalization failure, or semantic violation, or when the
    engine drops per-face origin ids.
    """
    if op not in ("union", "intersection", "difference"):
        raise ValueError(f"unknown op {op!r}")
    nA = len(np.asarray(proxy_a["F"]))
    nB = len(np.asarray(proxy_b["F"]))
    Va = np.ascontiguousarray(np.asarray(proxy_a["V"], dtype=np.float64))
    Vb = np.ascontiguousarray(np.asarray(proxy_b["V"], dtype=np.float64))
    scale = 1.0
    if len(Va):
        scale = max(scale, float(np.max(np.abs(Va))))
    if len(Vb):
        scale = max(scale, float(np.max(np.abs(Vb))))
    try:
        prepA = prepare_operand(proxy_a, 0, "input A")
        prepB = prepare_operand(proxy_b, nA, "input B")
        mA, mB = prepA["manifold"], prepB["manifold"]
        if op == "union":
            mR = mA + mB
        elif op == "intersection":
            mR = mA ^ mB
        else:
            mR = mA - mB
    except ArrangementError:
        raise
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
    check_boolean_semantics(V, F, prepA, prepB, op, scale)
    return {"V": V, "F": F, "origins": origins, "engine": ENGINE,
            "n_faces_a": nA, "n_faces_b": nB,
            "empty": len(F) == 0}
