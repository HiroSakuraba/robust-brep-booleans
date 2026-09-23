"""The 7-stage pipeline (design Section 4), prototype slice.

  Stage 0  ingest audit + tolerance ledger
  Stage 1  certified proxies (chordal-error-bounded tessellation)
  Stage 2  arrangement (proxy boolean; exact CGAL core slots in here)
  Stage 3  degeneracy tiers: Tier A analytic detection (this slice);
           Tier B/C reserved for freeform surfaces
  Stage 4  winding-number classification with margins (exact implicits)
  Stage 5  topology-invariant assembly + engine/classifier cross-check
  Stage 6  independent verification (closure, Euler, orientation, field,
           manifold cross-check)

Contract: never silently return a broken solid. If verification fails,
or an undecided degeneracy remains, the result carries an explicit
ambiguity report (or raises) -- it is never presented as clean.
"""

import numpy as np
from .ingest import audit_solid, ToleranceLedger, IngestError
from .proxy import certified_proxy
from .arrange import arrange, ArrangementError
from .classify import tierA_degeneracy
from .assemble import assemble
from .verify import verify


class AmbiguousResult(Exception):
    """Raised when the pipeline cannot certify a result.

    Carries the partial mesh and the ambiguity report -- the caller
    decides, the pipeline never decides silently.
    """
    def __init__(self, message, mesh, report):
        super().__init__(message)
        self.mesh = mesh
        self.report = report


def boolean(solidA, solidB, op, proxy_tol=1e-3):
    """Run the full pipeline. Returns (result_mesh_dict, full_report).

    Raises IngestError / ArrangementError / AmbiguousResult on failure.
    """
    if op not in ("union", "intersection", "difference"):
        raise ValueError(f"unknown op {op!r}")
    report = {"op": op, "stages": {}}

    # Stage 0
    ledger = ToleranceLedger(proxy_tol=proxy_tol)
    audits = [audit_solid(solidA), audit_solid(solidB)]
    report["stages"]["ingest"] = {"audits": audits,
                                  "ledger": ledger.as_dict()}
    for a in audits:
        if not a["ok"]:
            raise IngestError(f"ingest failed: {a['issues']}")

    # Stage 1
    proxyA = certified_proxy(solidA, proxy_tol)
    proxyB = certified_proxy(solidB, proxy_tol)
    report["stages"]["proxy"] = {
        "A": {k: proxyA[k] for k in ("chordal_error", "n_vertices", "n_faces")},
        "B": {k: proxyB[k] for k in ("chordal_error", "n_vertices", "n_faces")},
    }

    # Stage 3 (Tier A) -- run BEFORE arrangement so degeneracies are known
    # before any geometric decisions are made.
    margin = ledger.degeneracy_margin(proxyA["chordal_error"],
                                      proxyB["chordal_error"])
    degeneracies = tierA_degeneracy(solidA, solidB, margin)
    report["stages"]["degeneracy"] = {"tier": "A (analytic)",
                                      "findings": degeneracies,
                                      "margin": margin}

    # Stage 2
    arrangement = arrange(proxyA, proxyB, op)
    report["stages"]["arrangement"] = {
        "engine": arrangement["engine"],
        "n_faces": len(arrangement["F"]),
        "empty": arrangement.get("empty", False),
    }

    # Stages 4+5
    assembled = assemble(arrangement, solidA, solidB, proxyA, proxyB,
                         ledger, op)
    report["stages"]["assembly"] = {
        "n_on_margin_faces": len(assembled["on_faces"]),
        "engine_classifier_agreement": assembled["agreement"],
        "empty": assembled.get("empty", False),
    }

    # Stage 6
    accepted, vreport = verify(assembled, solidA, solidB, op)
    report["stages"]["verification"] = vreport
    report["accepted"] = accepted

    mesh = {"V": assembled["V"], "F": assembled["F"],
            "empty": assembled.get("empty", False)}

    ambiguities = []
    if degeneracies:
        ambiguities.append({"type": "degeneracy",
                            "findings": degeneracies})
    if len(assembled["on_faces"]) > 0:
        ambiguities.append({"type": "on_margin_faces",
                            "count": len(assembled["on_faces"])})
    if assembled["agreement"] < 1.0 - 1e-9:
        ambiguities.append({"type": "engine_classifier_disagreement",
                            "agreement": assembled["agreement"]})
    report["ambiguities"] = ambiguities

    if not accepted:
        raise AmbiguousResult(
            "Stage 6 verification FAILED -- result rejected, not returned "
            "as a valid solid.", mesh, report)
    return mesh, report
