"""Stage 5: assembly with per-face audit.

Attaches the origin-aware audit from classify.py to the arrangement mesh.
Any ambiguous or violating face blocks the result upstream -- nothing is
silently kept or dropped.
"""

import numpy as np
from .classify import audit_faces


def assemble(arrangement, solidA, solidB, proxyA, proxyB, ledger, op,
             skip_audit=False):
    """Attach the per-face audit. Returns the assembled record.

    skip_audit: for the Tier A identical-input fast path, where there is
    no engine decision to audit.
    """
    V, F = arrangement["V"], arrangement["F"]
    margin = ledger.degeneracy_margin(proxyA["chordal_error"],
                                      proxyB["chordal_error"])
    if arrangement.get("empty") or len(F) == 0:
        return {"V": V, "F": F, "empty": True, "margin": margin,
                "audit": None}
    origins = arrangement.get("origins")
    if origins is None or len(origins) != len(F):
        raise RuntimeError("arrangement lacks per-face origins")
    audit = None if skip_audit else audit_faces(
        F, V, origins, solidA, solidB, margin, op)
    return {"V": V, "F": F, "empty": False, "margin": margin,
            "audit": audit, "origins": origins}
