"""The 7-stage pipeline (design Section 4), prototype slice.

  Stage 0  ingest audit + tolerance ledger
  Stage 1  certified proxies (chordal-error-bounded tessellation)
  Stage 2  arrangement (proxy boolean in float64; exact CGAL core slots in)
  Stage 3  Tier A: exact degeneracy analysis of defining parameters
  Stage 4  per-face audit of the engine against the exact implicits
  Stage 5  assembly (provenance + audit attached)
  Stage 6  independent verification (8 checks)

Contract: never silently return a broken solid. Ambiguities BLOCK: any
unverifiable face, any engine-decision violation, or any failed check
raises AmbiguousResult carrying the partial mesh and the full report.
Tier A identical inputs (A op A) resolve exactly without the engine.
"""

import numpy as np
from .ingest import audit_solid, ToleranceLedger, IngestError
from .proxy import certified_proxy
from .arrange import arrange, ArrangementError
from .classify import tierA_degeneracy, identical_inputs
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
    margin = ledger.degeneracy_margin(proxyA["chordal_error"],
                                      proxyB["chordal_error"])

    # Stage 3 (Tier A), before any geometric decisions
    degeneracies = tierA_degeneracy(solidA, solidB, margin)
    report["stages"]["degeneracy"] = {"tier": "A (analytic)",
                                      "findings": degeneracies,
                                      "margin": margin}

    # Tier A exact resolution: identical inputs need no engine.
    fast_path = identical_inputs(solidA, solidB)
    if fast_path:
        report["stages"]["tierA_resolution"] = {
            "identical_inputs": True,
            "resolution": {"union": "A", "intersection": "A",
                           "difference": "empty"}[op],
        }

    # Stage 2
    if fast_path:
        if op == "difference":
            arrangement = {"V": np.zeros((0, 3)), "F": np.zeros((0, 3), int),
                           "engine": "tierA-exact", "empty": True}
        else:
            arrangement = {"V": proxyA["V"], "F": proxyA["F"],
                           "origins": np.array(["A"] * len(proxyA["F"])),
                           "engine": "tierA-exact", "empty": False}
    else:
        arrangement = arrange(proxyA, proxyB, op)
    report["stages"]["arrangement"] = {
        "engine": arrangement["engine"],
        "n_faces": len(arrangement["F"]),
        "empty": arrangement.get("empty", False),
    }

    # Stages 4+5
    assembled = assemble(arrangement, solidA, solidB, proxyA, proxyB,
                         ledger, op, skip_audit=fast_path)
    audit = assembled["audit"]
    if audit is not None:
        report["stages"]["assembly"] = {
            "n_faces": len(assembled["F"]),
            "n_verified": int(np.sum(audit["verified"])),
            "n_ambiguous": int(np.sum(audit["ambiguous"])),
            "n_violation": int(np.sum(audit["violation"])),
            "empty": False,
        }
    else:
        report["stages"]["assembly"] = {
            "audit": "skipped (Tier A exact resolution)",
            "empty": assembled.get("empty", False),
        }

    # Stage 6
    accepted, vreport = verify(assembled, solidA, solidB, op,
                               proxyA, proxyB)
    report["stages"]["verification"] = vreport

    mesh = {"V": assembled["V"], "F": assembled["F"],
            "empty": assembled.get("empty", False)}

    # Ambiguities BLOCK. Findings alone are notes; unverifiable faces,
    # engine violations, and failed checks raise.
    ambiguities = []
    if audit is not None:
        n_amb = int(np.sum(audit["ambiguous"]))
        n_vio = int(np.sum(audit["violation"]))
        if n_amb:
            ambiguities.append({"type": "unverifiable_faces",
                                "count": n_amb,
                                "faces": [int(i) for i in
                                          np.nonzero(audit["ambiguous"])[0][:25]]})
        if n_vio:
            ambiguities.append({"type": "engine_decision_violation",
                                "count": n_vio,
                                "faces": [int(i) for i in
                                          np.nonzero(audit["violation"])[0][:25]]})
    failed_checks = [k for k, c in vreport["checks"].items()
                     if not c["pass"]]
    for k in failed_checks:
        ambiguities.append({"type": "verification_failed", "check": k,
                            "info": vreport["checks"][k].get("info")})
    report["ambiguities"] = ambiguities
    report["degeneracy_findings"] = degeneracies  # informational
    report["accepted"] = accepted and not ambiguities

    if ambiguities or not accepted:
        raise AmbiguousResult(
            f"result not certified ({len(ambiguities)} blocking issue(s)); "
            "see report['ambiguities']", mesh, report)
    return mesh, report
