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

Two hardening rules live here (both added 24 Sept 2026 after the ABC/DeepCAD
independent-oracle stress rounds found silent wrong-accepts):

1. Inputs must be closed 2-manifolds. manifold3d silently returns a wrong
   mesh on non-manifold input (observed: 2 non-manifold edges in, impossible
   volumes out -- negative difference, intersection bigger than the input --
   while every mesh verify function passed). Such inputs are refused with a
   typed ArrangementError instead.

2. Each input is first reduced to its set-union (its own shells folded with
   union). manifold3d classifies each input's faces only against the OTHER
   mesh, so shells of the SAME input that interpenetrate are kept verbatim
   and their overlap volume is double-counted in union/difference (observed:
   +18.4% union volume on a clean closed 3-shell ABC model whose own shells
   interpenetrate; all four mesh verify functions passed on the wrong mesh
   because it is a valid closed mesh -- just not the set-union). Folding
   each input's shells first gives the cross-mesh boolean proper inputs.
   Single-shell inputs skip the fold: bit-identical behavior.
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


def _check_input_manifold(V, F, name):
    """Refuse inputs that are not closed 2-manifolds (typed, loud).

    Directed closure requires every undirected edge to appear exactly twice,
    once in each direction: this rejects open meshes (boundary edges),
    non-manifold edges (shared by 3+ faces), and inconsistently oriented
    meshes alike. Such inputs have no defined solid semantics; the engine
    may silently return a geometrically impossible mesh.
    """
    if len(F) == 0:
        return
    from .verify import check_closure_directed
    ok, info = check_closure_directed(V, F)
    if not ok:
        raise ArrangementError(
            f"{name}: input is not a closed 2-manifold ({info})")


def _split_shells(V, F, tag0):
    """One Manifold per connected shell; faces keep global face_id tags."""
    from .verify import shell_face_labels
    labels = shell_face_labels(F)
    mans = []
    tag = tag0
    for lab in sorted(np.unique(labels).tolist()):
        m = labels == lab
        Fm = F[m]
        used = np.unique(Fm)
        remap = np.full(len(V), -1, dtype=np.int64)
        remap[used] = np.arange(len(used))
        man, _ = _to_manifold(
            {"V": np.ascontiguousarray(V[used]),
             "F": np.ascontiguousarray(remap[Fm])}, tag)
        mans.append(man)
        tag += int(np.sum(m))
    return mans


def _set_union(mans):
    """Fold shells with union: the input mesh interpreted as a geometric set.

    Every fold step is a boolean of proper closed meshes, so the
    accumulator never carries interpenetrating shells. Raises
    ArrangementError (typed refusal) if the engine balks.
    """
    acc = mans[0]
    for m in mans[1:]:
        acc = acc + m
        if acc.status().name != "NoError":
            raise ArrangementError(
                "engine failed while merging input shells")
    return acc


def _prepare(proxy, tag0, name):
    """Validate + reduce one input to its set-union Manifold.

    Single-shell (and empty) inputs take exactly the old construction path.
    """
    from .verify import connected_shells
    V = np.ascontiguousarray(np.asarray(proxy["V"], dtype=np.float64))
    F = np.ascontiguousarray(np.asarray(proxy["F"], dtype=np.int64))
    _check_input_manifold(V, F, name)
    if len(F) == 0 or connected_shells(F) <= 1:
        man, _ = _to_manifold(proxy, tag0)
        return man
    return _set_union(_split_shells(V, F, tag0))


def arrange(proxy_a, proxy_b, op):
    """Boolean the two proxy meshes.

    Returns {'V','F','origins','engine'} with V float64, origins in {'A','B'}
    per result face. Raises ArrangementError on engine failure, on
    non-manifold input, or when the engine drops per-face origin ids.
    """
    if op not in ("union", "intersection", "difference"):
        raise ValueError(f"unknown op {op!r}")
    nA = len(np.asarray(proxy_a["F"]))
    nB = len(np.asarray(proxy_b["F"]))
    try:
        mA = _prepare(proxy_a, 0, "input A")
        mB = _prepare(proxy_b, nA, "input B")
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
    return {"V": V, "F": F, "origins": origins, "engine": ENGINE,
            "n_faces_a": nA, "n_faces_b": nB,
            "empty": len(F) == 0}


class ArrangementError(Exception):
    pass
