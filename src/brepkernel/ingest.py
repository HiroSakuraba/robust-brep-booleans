"""Stage 0: ingest + audit. Margin conventions live here."""

import numpy as np


class IngestError(Exception):
    pass


def audit_solid(solid):
    """Validate a solid's defining parameters. Returns an audit record."""
    issues = []
    lo, hi = solid.bbox()
    if not (np.all(np.isfinite(lo)) and np.all(np.isfinite(hi))):
        issues.append("non-finite bounding box")
    if np.any(hi <= lo):
        issues.append("degenerate bounding box")
    diag = float(np.linalg.norm(hi - lo))
    if diag == 0.0:
        issues.append("zero-size solid")
    return {"solid": repr(solid), "bbox_diag": diag, "issues": issues,
            "ok": len(issues) == 0}


class ToleranceLedger:
    """Every epsilon used by the pipeline, in one place.

    margin: a face centroid (or vertex) whose |implicit| against the other
    solid is <= margin cannot be classified reliably, because the proxy
    surface can deviate from the true surface by up to the chordal bound.
    The chordal bounds are rigorous (certified in solids.py), so
    margin = chordal_A + chordal_B + numeric_eps is principled, not tuned.
    """

    def __init__(self, proxy_tol, numeric_eps=1e-9):
        if not proxy_tol > 0:
            raise IngestError("proxy_tol must be positive")
        self.entries = {
            "proxy_tol": float(proxy_tol),   # certified chordal bound
            "numeric_eps": float(numeric_eps),
        }

    def degeneracy_margin(self, chordal_a, chordal_b):
        m = chordal_a + chordal_b + self.entries["numeric_eps"]
        self.entries["degeneracy_margin"] = float(m)
        return m

    def record(self, key, value):
        self.entries[key] = value

    def as_dict(self):
        return dict(self.entries)
