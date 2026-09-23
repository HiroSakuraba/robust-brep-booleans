"""Stage 0: ingest audit + tolerance ledger.

Every input is audited *before* any geometry is built. The tolerance
ledger records every epsilon the pipeline will use, so a result can
always be traced back to the tolerances that produced it (design G1).
"""

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
    """Every epsilon used by the pipeline, in one place."""

    def __init__(self, proxy_tol, degeneracy_margin_scale=2.0):
        if not proxy_tol > 0:
            raise IngestError("proxy_tol must be positive")
        self.entries = {
            "proxy_tol": float(proxy_tol),          # certified chordal bound
            "degeneracy_margin_scale": float(degeneracy_margin_scale),
        }

    def degeneracy_margin(self, chordal_a, chordal_b):
        """Margin below which a query point counts as ON a surface.

        A proxy face can deviate from the true surface by up to its
        chordal error, so a point closer than chordal_a + chordal_b to
        both surfaces cannot be classified by the proxy alone -- it
        needs Tier A exact analysis (design Stage 3).
        """
        return (self.entries["degeneracy_margin_scale"]
                * (chordal_a + chordal_b) + 1e-12)

    def record(self, key, value):
        self.entries[key] = value

    def as_dict(self):
        return dict(self.entries)
