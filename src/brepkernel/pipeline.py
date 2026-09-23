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


def _boolean_once(solidA, solidB, op, proxy_tol=1e-3, allow_skips=True):
    """Single-tolerance core of the pipeline (stages 0-6).

    Raises IngestError / ArrangementError / AmbiguousResult on failure.
    allow_skips: when False, a verification check that could not run
    (status 'skip', e.g. no exact Euler prediction for non-box inputs)
    blocks certification instead of being reported as
    'certified with skipped checks: ...'.
    """
    if op not in ("union", "intersection", "difference"):
        raise ValueError(f"unknown op {op!r}")
    report = {"op": op, "stages": {}}

    # Stage 0
    audits = [audit_solid(solidA), audit_solid(solidB)]
    # coordinate scale for scale-aware tolerances (item 4): rounding of
    # vertex coordinates happens at ~eps64 * max|coord|
    coord_scale = 1.0
    for a_solid in (solidA, solidB):
        lo, hi = a_solid.bbox()
        coord_scale = max(coord_scale,
                          float(np.max(np.abs(lo))),
                          float(np.max(np.abs(hi))))
    ledger = ToleranceLedger(proxy_tol=proxy_tol, coord_scale=coord_scale)
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
            "n_rescued": int(audit.get("n_rescued", 0)),
            "n_thin_features": len(audit.get("thin_features", [])),
            "empty": False,
        }
    else:
        report["stages"]["assembly"] = {
            "audit": "skipped (Tier A exact resolution)",
            "empty": assembled.get("empty", False),
        }

    # Stage 6
    accepted, vreport = verify(assembled, solidA, solidB, op,
                               proxyA, proxyB, allow_skips=allow_skips)
    report["stages"]["verification"] = vreport

    mesh = {"V": assembled["V"], "F": assembled["F"],
            "empty": assembled.get("empty", False)}

    # Ambiguities BLOCK. Findings alone are notes; unverifiable faces,
    # engine violations, and failed checks raise.
    ambiguities = []
    if audit is not None:
        thin = [int(i) for i in audit.get("thin_features", [])]
        thin_set = set(thin)
        n_vio = int(np.sum(audit["violation"]))
        amb_idx = [int(i) for i in np.nonzero(audit["ambiguous"])[0]
                   if int(i) not in thin_set]
        if amb_idx:
            ambiguities.append({"type": "unverifiable_faces",
                                "count": len(amb_idx),
                                "faces": amb_idx[:25]})
        if thin:
            ambiguities.append({
                "type": "sub_margin_thin_feature",
                "count": len(thin),
                "faces": thin[:25],
                "note": ("rescued band faces within 2x margin of an "
                         "opposite-facing result surface; patch rescue "
                         "withheld (possible hairline wall/bridge)")})
        if n_vio:
            ambiguities.append({"type": "engine_decision_violation",
                                "count": n_vio,
                                "faces": [int(i) for i in
                                          np.nonzero(audit["violation"])[0][:25]]})
    failed_checks = [k for k, c in vreport["checks"].items()
                     if c["status"] == "fail"]
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


def boolean(solidA, solidB, op, proxy_tol=1e-3, allow_skips=True):
    """Run the full pipeline. Returns (result_mesh_dict, full_report).

    Raises IngestError / ArrangementError / AmbiguousResult on failure.
    allow_skips: when False, a verification check that could not run
    (status 'skip', e.g. no exact Euler prediction for non-box inputs)
    blocks certification instead of being reported as
    'certified with skipped checks: ...'.

    Topology stability: after the primary run at proxy_tol certifies, the
    boolean is re-run at 2*proxy_tol (a different tessellation phase for
    curved solids). Shell count and chi are topological invariants of the
    true result, so if both runs certify but disagree, the inputs are in
    a near-degenerate region at this resolution and the result is refused
    instead of silently returning one phase's answer. A check run that
    does not certify is recorded but does not block the primary result.
    """
    from .verify import connected_shells, euler_chi
    mesh, report = _boolean_once(solidA, solidB, op, proxy_tol,
                                 allow_skips=allow_skips)
    check_tol = 2.0 * proxy_tol
    stability = {"check_tol": check_tol, "ran": False}
    topo2 = None
    try:
        mesh2, _ = _boolean_once(solidA, solidB, op, check_tol,
                                 allow_skips=True)
        topo2 = (connected_shells(mesh2["F"]),
                 euler_chi(mesh2["V"], mesh2["F"]))
        stability["ran"] = True
    except (AmbiguousResult, IngestError, ArrangementError) as e:
        stability["note"] = (f"check run at tol={check_tol:g} did not "
                             f"certify ({type(e).__name__}); primary "
                             "result kept")
    if topo2 is not None:
        topo1 = (connected_shells(mesh["F"]),
                 euler_chi(mesh["V"], mesh["F"]))
        stability["topology_at_tol"] = list(topo1)
        stability["topology_at_check_tol"] = list(topo2)
        stability["agree"] = topo1 == topo2
    report["stages"]["topology_stability"] = stability
    if topo2 is not None and topo1 != topo2:
        raise AmbiguousResult(
            f"topology differs across tessellation phases "
            f"({tuple(topo1)} at tol={proxy_tol:g} vs {tuple(topo2)} at "
            f"tol={check_tol:g}); near-degenerate region, refusing",
            mesh, report)
    return mesh, report
