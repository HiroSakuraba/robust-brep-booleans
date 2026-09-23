"""Stage 2: arrangement (proxy boolean).

The arrangement engine intersects the two proxy meshes. This prototype
uses manifold3d (robust-inexact: guaranteed manifold output, approximate
intersection curves). The exact core from the design (libigl/CGAL
MESH_BOOLEAN_TYPE_RESOLVE with the J birth-triangle map) slots in here
behind the same interface when the wheel is available.

Design rule honored: the arrangement is *geometry only*. Which faces
survive is decided in Stage 4 by exact classification, and the result
is never trusted without Stage 6 verification.
"""

import numpy as np

try:
    from manifold3d import Manifold, CrossSection  # noqa
    _HAS_MANIFOLD = True
except Exception:
    _HAS_MANIFOLD = False


class ArrangementError(Exception):
    pass


def _to_manifold(proxy):
    from manifold3d import Mesh, Manifold
    V = np.ascontiguousarray(proxy["V"], dtype=np.float32)
    F = np.ascontiguousarray(proxy["F"], dtype=np.uint32)
    mesh = Mesh(vert_properties=V, tri_verts=F)
    return Manifold(mesh)


def arrange(proxyA, proxyB, op):
    """Boolean the two proxy meshes. Returns {'V','F','engine'}.

    op in {'union','intersection','difference'} (difference = A - B).
    Raises ArrangementError if the engine reports failure.
    """
    if not _HAS_MANIFOLD:
        raise ArrangementError("manifold3d not installed")
    mA = _to_manifold(proxyA)
    mB = _to_manifold(proxyB)
    try:
        if op == "union":
            mR = mA + mB
        elif op == "intersection":
            mR = mA ^ mB
        elif op == "difference":
            mR = mA - mB
        else:
            raise ValueError(f"unknown op {op}")
    except Exception as e:
        raise ArrangementError(f"engine failed on {op}: {e}")
    if mR.status().name != "NoError":
        raise ArrangementError(f"engine returned error status {mR.status()}")
    mesh = mR.to_mesh()
    V = np.asarray(mesh.vert_properties[:, :3], dtype=float)
    F = np.asarray(mesh.tri_verts, dtype=np.int64)
    if len(F) == 0:
        # empty result is a legitimate outcome (e.g. A - A, disjoint intersect)
        return {"V": V.reshape(0, 3), "F": F.reshape(0, 3), "engine": "manifold3d",
                "empty": True}
    return {"V": V, "F": F, "engine": "manifold3d", "empty": False}
