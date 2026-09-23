"""Stage 1: certified proxy meshes.

A proxy is a triangle mesh with a *certificate*: a proven upper bound on
its chordal deviation from the true analytic surface. Downstream stages
may treat the mesh as geometry, but every topological decision must go
through the exact implicit (Tier A) or an explicit degeneracy analysis.
"""

import numpy as np
from . import solids


def certified_proxy(solid, tol):
    """Build a proxy mesh with chordal error certified <= tol.

    Returns dict with V (n,3), F (m,3), chordal_error, solid, tol.
    """
    V, F, chordal = solid.tessellate(tol)
    assert chordal <= tol * (1 + 1e-9), \
        f"certificate violated: {chordal} > {tol}"
    # sanity: all vertices lie on/near the true surface
    dev = np.abs(solid.implicit(V))
    return {
        "V": V, "F": F,
        "chordal_error": float(chordal),
        "max_vertex_deviation": float(np.max(dev)),
        "solid": solid,
        "tol": float(tol),
        "n_vertices": len(V), "n_faces": len(F),
    }
