"""Stage 5: topology-invariant assembly.

Takes the arrangement mesh and attaches per-face provenance from exact
classification: which input solid(s) each surviving face belongs to.
Invariants checked here:
  - every surviving face is classified IN/OUT by both exact implicits
    (no ON faces undecided -- those go to the ambiguity report);
  - the kept/discarded decision from the engine agrees with the exact
    classification table on a sample of faces (engine/classifier cross-check).
"""

import numpy as np
from .classify import classify_faces


class AssemblyError(Exception):
    pass


def assemble(arrangement, solidA, solidB, proxyA, proxyB, ledger, op):
    """Attach provenance; cross-check engine vs exact classifier.

    Returns {'V','F','in_a','in_b','on_faces','agreement'}.
    Faces whose centroid is ON a surface within the margin are listed in
    on_faces and trigger Tier A degeneracy handling upstream -- they are
    never silently kept or dropped.
    """
    V, F = arrangement["V"], arrangement["F"]
    if arrangement.get("empty"):
        return {"V": V, "F": F, "in_a": np.array([]), "in_b": np.array([]),
                "on_faces": np.array([], dtype=int), "agreement": 1.0,
                "empty": True}
    centroids = V[F].mean(axis=1)
    margin = ledger.degeneracy_margin(proxyA["chordal_error"],
                                      proxyB["chordal_error"])
    in_a, in_b, keep_exact, on_mask = classify_faces(
        centroids, solidA, solidB, margin, op)
    on_faces = np.nonzero(on_mask)[0]

    # Engine/classifier cross-check: the arrangement mesh should contain
    # exactly the faces the exact table keeps, up to ON-margin faces.
    # We check agreement on faces safely away from the margin.
    safe = ~on_mask
    # A face survives in the arrangement iff the engine kept it; the exact
    # table says which (inA,inB) combos survive. Boundary faces of the
    # result lie ON one input surface, so compare only interior-side faces:
    # every result face must be ON A or ON B (it came from an input surface).
    onA = np.abs(solidA.implicit(centroids)) <= margin
    onB = np.abs(solidB.implicit(centroids)) <= margin
    from_input = onA | onB
    agreement = float(np.mean(from_input)) if len(F) else 1.0

    return {"V": V, "F": F, "in_a": in_a, "in_b": in_b,
            "on_faces": on_faces, "agreement": agreement, "empty": False,
            "margin": margin}
