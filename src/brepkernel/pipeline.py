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


class BRepAmbiguousResult(Exception):
    """Tier B/C B-rep refusal with the completed stage report attached."""

    def __init__(self, message, report, cause=None):
        super().__init__(message)
        self.report = report
        self.cause = cause


def _empty_brep_compound():
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Compound
    c = TopoDS_Compound()
    b = BRep_Builder()
    b.MakeCompound(c)
    return c


def boolean_brep(shapeA, shapeB, op, *, base_tol=1e-7,
                 broadphase_pad=None, chord_tol=None, contact_tol=None,
                 fuzzy=0.0, parallel=True, use_obb=True,
                 tangent_sin_tol=1e-4, area_rel_tol=2e-6,
                 sew_tol=None):
    """Run the exact trimmed-B-rep Tier B/C pipeline.

    Returns (TopoDS_Shape, report) and leaves the existing mesh/proxy
    boolean() API unchanged.

    shapeA and shapeB may be OCCT TopoDS shapes or already-indexed BRepModel
    objects. Ambiguous contacts, failed p-curve verification, split
    inconsistencies, open/non-manifold sewing, and boundary-only patch
    classifications refuse through BRepAmbiguousResult rather than being
    converted into a guessed result.

    Difference means A - B.
    """
    from OCP.BRepCheck import BRepCheck_Analyzer

    from .freeform import FreeformError
    from .step_ingest import BRepModel, index_shape
    from .intersection import intersect_models
    from .split import split_models
    from .assembly import assemble_boolean
    from .same_domain import same_domain_models

    if op not in ("union", "intersection", "difference"):
        raise ValueError(f"unknown op {op!r}")
    if not base_tol > 0:
        raise ValueError("base_tol must be positive")
    if fuzzy < 0:
        raise ValueError("fuzzy must be >= 0")

    if contact_tol is None:
        contact_tol = 4.0 * float(base_tol)
    if broadphase_pad is None:
        # Preserve tolerance-near contacts for the exact contact classifier;
        # a fast AABB rejection must not make them disappear.
        broadphase_pad = max(float(contact_tol), 4.0 * float(base_tol))
    if chord_tol is None:
        chord_tol = max(4.0 * float(base_tol), 1e-9)

    report = {
        "op": op,
        "route": "Tier B/C exact trimmed B-rep",
        "stages": {},
        "accepted": False,
    }

    def refuse(stage, exc):
        report["refusal"] = {
            "stage": stage,
            "type": type(exc).__name__,
            "kind": getattr(exc, "kind", type(exc).__name__),
            "message": str(exc),
        }
        raise BRepAmbiguousResult(
            f"B-rep result refused in {stage}: {exc}", report, exc) from exc

    try:
        a = shapeA if isinstance(shapeA, BRepModel) else index_shape(shapeA)
        b = shapeB if isinstance(shapeB, BRepModel) else index_shape(shapeB)
    except FreeformError as exc:
        refuse("ingest", exc)

    report["stages"]["ingest"] = {
        "A": {
            "solids": len(a.solids), "shells": len(a.shells),
            "faces": len(a.faces), "freeform_accels": len(a.nurbs_faces),
            "local_patches": sum(
                len(fr.freeform.index.records)
                for fr in a.faces if fr.freeform is not None),
        },
        "B": {
            "solids": len(b.solids), "shells": len(b.shells),
            "faces": len(b.faces), "freeform_accels": len(b.nurbs_faces),
            "local_patches": sum(
                len(fr.freeform.index.records)
                for fr in b.faces if fr.freeform is not None),
        },
    }

    # Exact oriented topological equality is the cheapest identity path.
    # IsSame() is intentionally NOT used: OCCT documents that IsSame ignores
    # orientation, so a reversed view of the same TShape would otherwise be
    # accepted as identical material.
    if a.shape.IsEqual(b.shape):
        if op == "difference":
            out = _empty_brep_compound()
            resolution = "empty"
        else:
            out = a.shape
            resolution = "A"
        report["stages"]["identity"] = {
            "exact_topological_identity": True,
            "same_domain_equivalent": True,
            "resolution": resolution,
        }
        report["stages"]["verification"] = {
            "brep_valid": True if op == "difference"
            else bool(BRepCheck_Analyzer(out, True).IsValid()),
            "identity_exact": True,
        }
        if not report["stages"]["verification"]["brep_valid"]:
            exc = FreeformError("exact-identity result is not B-rep valid",
                                "IdentityResultInvalid")
            refuse("verification", exc)
        report["accepted"] = True
        return out, report

    # Independently constructed B-reps can represent the same material
    # boundary without sharing a TShape. Use the strict optional recognizer;
    # failure means "not proven equivalent", never "different".
    sd = same_domain_models(
        a, b, base_tol=float(base_tol),
        fuzz=max(float(base_tol), float(fuzzy)))
    def canonical_report(ev):
        if ev is None:
            return None
        return {
            "changed": ev.changed,
            "faces_before": ev.faces_before,
            "faces_after": ev.faces_after,
            "shells_before": ev.shells_before,
            "shells_after": ev.shells_after,
            "solids_before": ev.solids_before,
            "solids_after": ev.solids_after,
            "bbox_error": ev.bbox_error,
            "volume_before": ev.volume_before,
            "volume_after": ev.volume_after,
            "volume_rel_error": ev.volume_rel_error,
        }

    report["stages"]["same_domain"] = {
        "equivalent": sd.equivalent,
        "reason": sd.reason,
        "matched_faces": len(sd.matches),
        "candidate_counts": list(sd.candidate_counts),
        "signed_volume_A": sd.signed_volume_a,
        "signed_volume_B": sd.signed_volume_b,
        "bbox_error": sd.bbox_error,
        "canonicalized": sd.canonicalized,
        "canonical_A": canonical_report(sd.canonical_a),
        "canonical_B": canonical_report(sd.canonical_b),
    }
    if sd.equivalent:
        if op == "difference":
            out = _empty_brep_compound()
            resolution = "empty"
        else:
            out = a.shape
            resolution = "A"
        report["stages"]["same_domain"]["resolution"] = resolution
        report["stages"]["same_domain"]["face_matches"] = [
            {
                "A": ev.face_a,
                "B": ev.face_b,
                "bbox_error": ev.bbox_error,
                "area_rel_error": ev.area_rel_error,
                "perimeter_rel_error": ev.perimeter_rel_error,
                "edges": ev.edge_count,
                "wires": ev.wire_count,
            }
            for ev in sd.matches
        ]
        report["stages"]["verification"] = {
            "brep_valid": True if op == "difference"
            else bool(BRepCheck_Analyzer(out, True).IsValid()),
            "strict_same_domain": True,
        }
        if not report["stages"]["verification"]["brep_valid"]:
            exc = FreeformError("same-domain fast-path result is invalid",
                                "SameDomainResultInvalid")
            refuse("verification", exc)
        report["accepted"] = True
        return out, report

    try:
        ix = intersect_models(
            a, b, broadphase_pad=float(broadphase_pad),
            base_tol=float(base_tol), chord_tol=float(chord_tol),
            contact_tol=float(contact_tol), fuzzy=float(fuzzy),
            parallel=bool(parallel), use_obb=bool(use_obb),
            tangent_sin_tol=float(tangent_sin_tol))
    except FreeformError as exc:
        refuse("intersection", exc)

    report["stages"]["intersection"] = {
        "candidate_face_pairs": ix.candidate_pairs,
        "section_calls": ix.section_calls,
        "verified_edges": ix.verified_edges,
        "point_contacts": ix.point_contacts,
        "ambiguous_contacts": ix.ambiguous_contacts,
        "face_pairs_skipped": ix.skipped_by_broadphase,
    }

    try:
        sp = split_models(
            a, b, ix, base_tol=float(base_tol),
            area_rel_tol=float(area_rel_tol), fuzzy=float(fuzzy),
            parallel=bool(parallel), use_obb=bool(use_obb))
    except FreeformError as exc:
        refuse("split", exc)

    report["stages"]["split"] = {
        "split_calls": sp.split_calls,
        "affected_faces_A": sp.affected_faces_a,
        "affected_faces_B": sp.affected_faces_b,
        "unresolved_contacts": list(sp.unresolved_contacts),
    }

    try:
        assembled = assemble_boolean(
            a, b, sp, op, base_tol=float(base_tol), sew_tol=sew_tol)
    except FreeformError as exc:
        refuse("assembly", exc)

    report["stages"]["assembly"] = {
        "selected_faces": assembled.selected_faces,
        "shells": len(assembled.shells),
        "solids": len(assembled.solids),
        "free_edges": assembled.free_edges,
        "multiple_edges": assembled.multiple_edges,
        "volume": assembled.volume,
        "empty": assembled.is_empty,
        "section_payloads": [
            {
                "ref": [p.face_a, p.face_b, p.section_edge_index],
                "samples": int(len(p.parameters)),
                "result_edges": list(p.result_edge_indices),
                "edge_tolerance": p.edge_tolerance,
                "verify_tolerance": p.verify_tolerance,
                "max_surface_error_A": p.max_surface_error_a,
                "max_surface_error_B": p.max_surface_error_b,
                "max_cross_surface_error": p.max_cross_surface_error,
                "min_transversality": p.min_transversality,
                "max_transversality": p.max_transversality,
                "risk_flags": list(p.risk_flags),
                "repaired_same_parameter": p.repaired_same_parameter,
            }
            for p in assembled.section_payloads
        ],
        "edge_lineage": {
            "result_edges": len(assembled.edge_lineage),
            "boolean_section_edges": sum(
                1 for e in assembled.edge_lineage
                if e.provenance_kind == "boolean_section"),
            "source_boundary_edges": sum(
                1 for e in assembled.edge_lineage
                if e.provenance_kind == "source_boundary"),
            "unattributed_edges": sum(
                1 for e in assembled.edge_lineage
                if e.provenance_kind == "unattributed"),
            "records": [
                {
                    "edge": e.result_edge_index,
                    "kind": e.provenance_kind,
                    "operands": list(e.operands),
                    "parent_faces": [list(x) for x in e.parent_faces],
                    "piece_refs": [list(x) for x in e.piece_refs],
                    "intersection_refs": [list(x) for x in e.intersection_refs],
                    "verified_pcurves": e.verified_pcurves,
                }
                for e in assembled.edge_lineage
            ],
        },
        "decisions": [
            {
                "operand": d.operand,
                "parent_face_id": d.parent_face_id,
                "piece_index": d.piece_index,
                "classification": d.classification,
                "keep": d.keep,
                "reversed": d.reverse_for_difference,
            }
            for d in assembled.decisions
        ],
    }

    valid = (True if assembled.is_empty
             else bool(BRepCheck_Analyzer(assembled.shape, True).IsValid()))
    report["stages"]["verification"] = {
        "brep_valid": valid,
        "closed": assembled.free_edges == 0,
        "manifold_edges": assembled.multiple_edges == 0,
        "unresolved_contacts": len(sp.unresolved_contacts),
    }
    if not valid:
        exc = FreeformError("final B-rep validity check failed",
                            "FinalBRepInvalid")
        refuse("verification", exc)

    report["accepted"] = True
    return assembled.shape, report
