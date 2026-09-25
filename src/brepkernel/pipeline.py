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
                 tangent_sin_tol=1e-4, max_section_tol=None,
                 area_rel_tol=2e-6, sew_tol=None,
                 include_full_evidence=False,
                 shadow_section_crosscheck=False,
                 crosscheck_ops=False,
                 allow_nonmanifold=False):
    """Run the exact trimmed-B-rep Tier B/C pipeline.

    Returns (TopoDS_Shape, report) and leaves the existing mesh/proxy
    boolean() API unchanged.

    shapeA and shapeB may be OCCT TopoDS shapes or already-indexed BRepModel
    objects. Ambiguous contacts, failed p-curve verification, split
    inconsistencies, open/non-manifold sewing, and boundary-only patch
    classifications refuse through BRepAmbiguousResult rather than being
    converted into a guessed result.

    Difference means A - B.

    If include_full_evidence is True, the report also contains the complete
    verified section sample payloads (parameters, XYZ, and UV on both input
    faces) converted to ordinary Python lists so the report can be serialized
    directly to JSON. The default remains a compact summary.

    shadow_section_crosscheck=True additionally runs OCCT's alternative
    non-approximated section construction and requires geometric agreement.
    It is an optional stress mode, not part of the default certification path:
    on difficult freeform walking intersections that alternative construction
    can be less accurate than the primary representation.

    crosscheck_ops=True additionally runs the two other Boolean operations on
    the same inputs and requires the volume identities
    vol(AuB)+vol(AnB)=vol(A)+vol(B) and vol(A-B)=vol(A)-vol(AnB) within the
    scale-aware volume tolerance. A violation is a typed refusal with
    kind=OperationIdentityFailed. Off by default because it costs three
    pipeline runs; intended for CI-level tests, not every call.

    allow_nonmanifold=True accepts a touching-only union as a compound
    of the touching solids (the contact is recorded in the report)
    instead of refusing with kind=NonManifoldResult.
    """
    from time import perf_counter
    from OCP.BRepCheck import BRepCheck_Analyzer

    from .freeform import FreeformError
    from .step_ingest import BRepModel, index_shape, model_max_tolerance
    from .intersection import intersect_models
    from .split import split_models
    from .assembly import assemble_boolean, _shape_volume
    from .same_domain import same_domain_models

    if op not in ("union", "intersection", "difference"):
        raise ValueError(f"unknown op {op!r}")
    if not base_tol > 0:
        raise ValueError("base_tol must be positive")
    if fuzzy < 0:
        raise ValueError("fuzzy must be >= 0")

    if contact_tol is None:
        contact_tol = 4.0 * float(base_tol)
    _pad_explicit = broadphase_pad is not None
    if broadphase_pad is None:
        # Preserve tolerance-near contacts for the exact contact classifier;
        # a fast AABB rejection must not make them disappear.
        broadphase_pad = max(float(contact_tol), 4.0 * float(base_tol))
    # Keep an unspecified chord tolerance as None.  The intersection verifier
    # can then derive a local sampling tolerance from the *actual* accepted
    # OCCT edge/face tolerance.  Turning None into a global base_tol-derived
    # number here caused severe over-sampling on otherwise well-bounded NURBS
    # sections. Explicit caller values are still honored exactly.
    t_total = perf_counter()
    report = {
        "op": op,
        "route": "Tier B/C exact trimmed B-rep",
        "stages": {},
        "timings_ms": {},
        "accepted": False,
    }

    def refuse(stage, exc):
        report["refusal"] = {
            "stage": stage,
            "type": type(exc).__name__,
            "kind": getattr(exc, "kind", type(exc).__name__),
            "message": str(exc),
        }
        # Structured refusal evidence (e.g. the coincidence deviation
        # and the tolerances consulted) rides along when provided.
        evidence = getattr(exc, "evidence", None)
        if evidence is not None:
            report["refusal"]["evidence"] = evidence
        raise BRepAmbiguousResult(
            f"B-rep result refused in {stage}: {exc}", report, exc) from exc

    # Options forwarded verbatim to the two companion runs of the optional
    # cross-operation identity check. crosscheck_ops itself is forced off
    # there so the check never recurses.
    _pass_kwargs = dict(
        base_tol=base_tol, broadphase_pad=broadphase_pad,
        chord_tol=chord_tol, contact_tol=contact_tol, fuzzy=fuzzy,
        parallel=parallel, use_obb=use_obb,
        tangent_sin_tol=tangent_sin_tol, max_section_tol=max_section_tol,
        area_rel_tol=area_rel_tol, sew_tol=sew_tol,
        include_full_evidence=include_full_evidence,
        shadow_section_crosscheck=shadow_section_crosscheck)

    def _finalize(out, result_report):
        """Run the optional cross-operation identity check (G1).

        With crosscheck_ops=True the two other Boolean operations are
        computed on the same inputs and the volume identities
        vol(AuB)+vol(AnB)=vol(A)+vol(B) and vol(A-B)=vol(A)-vol(AnB) must
        hold within the scale-aware volume tolerance. A violation is a
        typed refusal with kind=OperationIdentityFailed. A companion run
        that itself refuses is recorded as skipped, never silently passed.
        """
        if not crosscheck_ops:
            return out, result_report
        others = [o for o in ("union", "intersection", "difference")
                  if o != op]
        vols = {}
        tols = []
        skipped = {}

        def _safe_volume(shape):
            # Empty results (no solids) have volume zero; the adaptive
            # integrator reports an error on them instead of 0.0.
            from OCP.TopAbs import TopAbs_SOLID
            from OCP.TopExp import TopExp_Explorer
            if shape.IsNull():
                return 0.0
            if not TopExp_Explorer(shape, TopAbs_SOLID).More():
                return 0.0
            return abs(float(_shape_volume(shape)))

        def _report_tol(rep):
            ver = rep.get("stages", {}).get("verification", {})
            t = ver.get("volume_tolerance")
            return float(t) if t is not None else None

        t_primary = _report_tol(result_report)
        if t_primary is not None:
            tols.append(t_primary)
        for other in others:
            try:
                other_out, other_report = boolean_brep(
                    shapeA, shapeB, other, crosscheck_ops=False,
                    **_pass_kwargs)
            except BRepAmbiguousResult as exc:
                skipped[other] = exc.report.get("refusal", {}).get(
                    "kind", type(exc).__name__)
                continue
            vols[other] = _safe_volume(other_out)
            t_other = _report_tol(other_report)
            if t_other is not None:
                tols.append(t_other)
        vols[op] = _safe_volume(out)
        va = _safe_volume(a.shape)
        vb = _safe_volume(b.shape)
        volume_tol = max(tols) if tols else 2e-8 * max(va, vb, 1.0)

        identities = {}
        failed = []
        if "union" in vols and "intersection" in vols:
            err = abs(vols["union"] + vols["intersection"] - va - vb)
            ok = err <= volume_tol
            identities["union_plus_intersection"] = {
                "error": err, "tolerance": volume_tol, "ok": ok}
            if not ok:
                failed.append(("union_plus_intersection", err))
        else:
            identities["union_plus_intersection"] = {
                "checked": False,
                "reason": "companion operation refused",
                "skipped": skipped,
            }
        if "difference" in vols and "intersection" in vols:
            err = abs(vols["difference"] - (va - vols["intersection"]))
            ok = err <= volume_tol
            identities["difference_from_intersection"] = {
                "error": err, "tolerance": volume_tol, "ok": ok}
            if not ok:
                failed.append(("difference_from_intersection", err))
        else:
            identities["difference_from_intersection"] = {
                "checked": False,
                "reason": "companion operation refused",
                "skipped": skipped,
            }
        result_report["stages"]["crosscheck"] = {
            "enabled": True,
            "volumes": vols,
            "input_volumes": {"A": va, "B": vb},
            "tolerance": volume_tol,
            "identities": identities,
            "skipped_operations": skipped,
        }
        if failed:
            names = ", ".join(f"{n} err={e:.6g}" for n, e in failed)
            exc = FreeformError(
                f"operation volume identities violated: {names} "
                f"(tolerance {volume_tol:.6g})",
                "OperationIdentityFailed")
            result_report["accepted"] = False
            refuse("crosscheck", exc)
        return out, result_report

    t_stage = perf_counter()
    try:
        a = shapeA if isinstance(shapeA, BRepModel) else index_shape(shapeA)
        b = shapeB if isinstance(shapeB, BRepModel) else index_shape(shapeB)
    except FreeformError as exc:
        report["timings_ms"]["ingest"] = (perf_counter() - t_stage) * 1000.0
        report["timings_ms"]["total"] = (perf_counter() - t_total) * 1000.0
        refuse("ingest", exc)
    report["timings_ms"]["ingest"] = (perf_counter() - t_stage) * 1000.0

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

    # G6: make the broad phase tolerance-aware. A face whose OCCT
    # tolerance is t is geometrically uncertain over a band of width
    # about 2t around it; a candidate test padded by less than that can
    # silently drop a genuinely ambiguous contact as "disjoint".
    # pad = max(contact_tol, 2 * max tolerance over all vertices, edges
    # and faces of both models). This is additive padding, never a
    # loosening of an acceptance check: it only widens the candidate
    # set, which can add typed refusals but never new acceptances.
    _tol_max = max(model_max_tolerance(a), model_max_tolerance(b))
    if not _pad_explicit:
        broadphase_pad = max(float(broadphase_pad), 2.0 * float(_tol_max))
    report["broadphase_pad"] = float(broadphase_pad)
    report["broadphase_max_tolerance"] = float(_tol_max)

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
        t_verify = perf_counter()
        report["stages"]["verification"] = {
            "brep_valid": True if op == "difference"
            else bool(BRepCheck_Analyzer(out, True).IsValid()),
            "identity_exact": True,
        }
        report["timings_ms"]["verification"] = (
            perf_counter() - t_verify) * 1000.0
        if not report["stages"]["verification"]["brep_valid"]:
            exc = FreeformError("exact-identity result is not B-rep valid",
                                "IdentityResultInvalid")
            refuse("verification", exc)
        report["accepted"] = True
        report["timings_ms"]["total"] = (perf_counter() - t_total) * 1000.0
        return _finalize(out, report)

    # Independently constructed B-reps can represent the same material
    # boundary without sharing a TShape. Use the strict optional recognizer;
    # failure means "not proven equivalent", never "different".
    t_stage = perf_counter()
    sd = same_domain_models(
        a, b, base_tol=float(base_tol),
        fuzz=max(float(base_tol), float(fuzzy)))
    report["timings_ms"]["same_domain"] = (
        perf_counter() - t_stage) * 1000.0
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
        t_verify = perf_counter()
        report["stages"]["verification"] = {
            "brep_valid": True if op == "difference"
            else bool(BRepCheck_Analyzer(out, True).IsValid()),
            "strict_same_domain": True,
        }
        report["timings_ms"]["verification"] = (
            perf_counter() - t_verify) * 1000.0
        if not report["stages"]["verification"]["brep_valid"]:
            exc = FreeformError("same-domain fast-path result is invalid",
                                "SameDomainResultInvalid")
            refuse("verification", exc)
        report["accepted"] = True
        report["timings_ms"]["total"] = (perf_counter() - t_total) * 1000.0
        return _finalize(out, report)

    t_stage = perf_counter()
    try:
        ix = intersect_models(
            a, b, broadphase_pad=float(broadphase_pad),
            base_tol=float(base_tol),
            chord_tol=(None if chord_tol is None else float(chord_tol)),
            contact_tol=float(contact_tol), fuzzy=float(fuzzy),
            parallel=bool(parallel), use_obb=bool(use_obb),
            tangent_sin_tol=float(tangent_sin_tol),
            max_section_tol=(None if max_section_tol is None
                             else float(max_section_tol)),
            crosscheck_nonapprox=bool(shadow_section_crosscheck))
    except FreeformError as exc:
        report["timings_ms"]["intersection"] = (
            perf_counter() - t_stage) * 1000.0
        report["timings_ms"]["total"] = (perf_counter() - t_total) * 1000.0
        refuse("intersection", exc)
    report["timings_ms"]["intersection"] = (
        perf_counter() - t_stage) * 1000.0

    report["stages"]["intersection"] = {
        "candidate_face_pairs": ix.candidate_pairs,
        "section_calls": ix.section_calls,
        "verified_edges": ix.verified_edges,
        "point_contacts": ix.point_contacts,
        "ambiguous_contacts": ix.ambiguous_contacts,
        "face_pairs_skipped": ix.skipped_by_broadphase,
        "shadow_section_calls": ix.shadow_section_calls,
        "shadow_verified_edges": ix.shadow_verified_edges,
        "max_shadow_distance": ix.max_shadow_distance,
        "completeness_probes": ix.completeness_probes,
        "raw_curve_count": ix.raw_curve_count,
        "raw_trimmed_components": ix.raw_trimmed_components,
        "raw_unmatched_components": ix.raw_unmatched_components,
        "completeness_max_distance": ix.completeness_max_distance,
        "approximation_gaps": ix.approximation_gaps,
        "approximation_gap_count": len(ix.approximation_gaps),
        "completeness_components": ix.completeness_components,
        "coincident_pairs": len(ix.coincident_pairs),
        "boundary_contact_edges": len(ix.boundary_edges),
    }

    t_stage = perf_counter()
    try:
        sp = split_models(
            a, b, ix, base_tol=float(base_tol),
            area_rel_tol=float(area_rel_tol), fuzzy=float(fuzzy),
            parallel=bool(parallel), use_obb=bool(use_obb))
    except FreeformError as exc:
        report["timings_ms"]["split"] = (
            perf_counter() - t_stage) * 1000.0
        report["timings_ms"]["total"] = (perf_counter() - t_total) * 1000.0
        refuse("split", exc)
    report["timings_ms"]["split"] = (
        perf_counter() - t_stage) * 1000.0

    report["stages"]["split"] = {
        "split_calls": sp.split_calls,
        "affected_faces_A": sp.affected_faces_a,
        "affected_faces_B": sp.affected_faces_b,
        "reused_seam_edges_A": sp.reused_seam_edges_a,
        "reused_seam_edges_B": sp.reused_seam_edges_b,
        "shared_seam_refusals": sp.shared_seam_refusals,
        "unresolved_contacts": list(sp.unresolved_contacts),
        "boundary_contact_edges": len(sp.boundary_edges),
    }

    t_stage = perf_counter()
    try:
        assembled = assemble_boolean(
            a, b, sp, op, base_tol=float(base_tol), sew_tol=sew_tol,
            allow_nonmanifold=bool(allow_nonmanifold))
    except FreeformError as exc:
        report["timings_ms"]["assembly"] = (
            perf_counter() - t_stage) * 1000.0
        report["timings_ms"]["total"] = (perf_counter() - t_total) * 1000.0
        refuse("assembly", exc)
    report["timings_ms"]["assembly"] = (
        perf_counter() - t_stage) * 1000.0

    section_sample_counts = [
        int(len(p.parameters)) for p in assembled.section_payloads]
    report["stages"]["assembly"] = {
        "selected_faces": assembled.selected_faces,
        "shells": len(assembled.shells),
        "solids": len(assembled.solids),
        "free_edges": assembled.free_edges,
        "multiple_edges": assembled.multiple_edges,
        "volume": assembled.volume,
        "empty": assembled.is_empty,
        "notes": list(assembled.notes),
        "section_sampling": {
            "sections": len(section_sample_counts),
            "total_samples": int(sum(section_sample_counts)),
            "max_samples": int(max(section_sample_counts, default=0)),
            "mean_samples": (float(sum(section_sample_counts))
                             / len(section_sample_counts)
                             if section_sample_counts else 0.0),
        },
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
                "exact_curve_on_surface_checked":
                    p.exact_curve_on_surface_checked,
                "exact_surface_error_A": p.exact_surface_error_a,
                "exact_surface_error_B": p.exact_surface_error_b,
                "shadow_crosschecked": p.shadow_crosschecked,
                "shadow_max_distance": p.shadow_max_distance,
                "shadow_length_rel_error": p.shadow_length_rel_error,
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
            "coincident_boundary_edges": sum(
                1 for e in assembled.edge_lineage
                if e.provenance_kind == "coincident_boundary"),
            "records": [
                {
                    "edge": e.result_edge_index,
                    "kind": e.provenance_kind,
                    "operands": list(e.operands),
                    "parent_faces": [list(x) for x in e.parent_faces],
                    "piece_refs": [list(x) for x in e.piece_refs],
                    "intersection_refs": [list(x) for x in e.intersection_refs],
                    "source_boundary_refs": [list(x) for x in
                                             e.source_boundary_refs],
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

    if include_full_evidence:
        report["stages"]["assembly"]["full_section_payloads"] = [
            {
                "ref": [p.face_a, p.face_b, p.section_edge_index],
                "parameters": p.parameters.tolist(),
                "xyz": p.xyz.tolist(),
                "uv_A": p.uv_a.tolist(),
                "uv_B": p.uv_b.tolist(),
                "edge_tolerance": p.edge_tolerance,
                "verify_tolerance": p.verify_tolerance,
                "max_surface_error_A": p.max_surface_error_a,
                "max_surface_error_B": p.max_surface_error_b,
                "max_cross_surface_error": p.max_cross_surface_error,
                "min_transversality": p.min_transversality,
                "max_transversality": p.max_transversality,
                "risk_flags": list(p.risk_flags),
                "repaired_same_parameter": p.repaired_same_parameter,
                "exact_curve_on_surface_checked":
                    p.exact_curve_on_surface_checked,
                "exact_surface_error_A": p.exact_surface_error_a,
                "exact_surface_error_B": p.exact_surface_error_b,
                "shadow_crosschecked": p.shadow_crosschecked,
                "shadow_max_distance": p.shadow_max_distance,
                "shadow_length_rel_error": p.shadow_length_rel_error,
                "result_edges": list(p.result_edge_indices),
            }
            for p in assembled.section_payloads
        ]

    t_stage = perf_counter()
    valid = (True if assembled.is_empty
             else bool(BRepCheck_Analyzer(assembled.shape, True).IsValid()))

    # Provenance is part of the Tier B/C acceptance contract, not optional
    # debugging metadata. A final edge that cannot be traced to a selected
    # patch/source boundary or a verified section curve means topology changed
    # somewhere in splitting/sewing without an auditable cause.
    unattributed = [
        e.result_edge_index for e in assembled.edge_lineage
        if e.provenance_kind == "unattributed"]
    bad_sections = [
        e.result_edge_index for e in assembled.edge_lineage
        if e.provenance_kind == "boolean_section"
        and (not e.verified_pcurves or not e.intersection_refs)]
    # G2.6: coincident_boundary edges must carry the other operand's
    # source edge IDs; a bare claim without refs is unaudited.
    bad_coincident = [
        e.result_edge_index for e in assembled.edge_lineage
        if e.provenance_kind == "coincident_boundary"
        and not e.source_boundary_refs]
    complete_lineage = (not unattributed and not bad_sections
                        and not bad_coincident)
    exact_curve_surface_complete = all(
        p.exact_curve_on_surface_checked
        and p.exact_surface_error_a is not None
        and p.exact_surface_error_b is not None
        and p.exact_surface_error_a <= p.verify_tolerance
        and p.exact_surface_error_b <= p.verify_tolerance
        for p in assembled.section_payloads)
    shadow_complete = (
        all(p.shadow_crosschecked for p in assembled.section_payloads)
        if shadow_section_crosscheck else None)

    # Cheap operation-level volume invariants catch catastrophic selection or
    # shell-orientation errors without using a second Boolean engine.
    va = abs(float(_shape_volume(a.shape)))
    vb = abs(float(_shape_volume(b.shape)))
    vr = abs(float(assembled.volume))
    all_faces = a.faces + b.faces
    if all_faces:
        lo = np.min(np.vstack([fr.bbox_lo for fr in all_faces]), axis=0)
        hi = np.max(np.vstack([fr.bbox_hi for fr in all_faces]), axis=0)
        volume_scale = max(float(np.linalg.norm(hi - lo)), 1.0)
    else:
        volume_scale = 1.0
    volume_tol = max(
        1e-10,
        256.0 * float(base_tol) * volume_scale * volume_scale,
        2e-8 * max(va, vb, vr, 1.0))
    if op == "union":
        volume_bounds_ok = (
            vr + volume_tol >= max(va, vb)
            and vr <= va + vb + volume_tol)
        volume_bounds = [max(va, vb), va + vb]
    elif op == "intersection":
        volume_bounds_ok = vr <= min(va, vb) + volume_tol
        volume_bounds = [0.0, min(va, vb)]
    elif op == "difference":
        volume_bounds_ok = (vr <= va + volume_tol) and (
            vr + volume_tol >= va - vb)
        volume_bounds = [max(0.0, va - vb), va]

    report["stages"]["verification"] = {
        "brep_valid": valid,
        "closed": assembled.free_edges == 0,
        "manifold_edges": assembled.multiple_edges == 0,
        "unresolved_contacts": len(sp.unresolved_contacts),
        "complete_edge_lineage": complete_lineage,
        "unattributed_edges": unattributed,
        "coincident_boundary_edges_missing_source": bad_coincident,
        "section_edges_missing_verified_pcurves": bad_sections,
        "exact_curve_on_surface_complete": exact_curve_surface_complete,
        "shadow_section_crosscheck_requested":
            bool(shadow_section_crosscheck),
        "shadow_section_crosscheck_complete": shadow_complete,
        "volume_bounds_ok": volume_bounds_ok,
        "volume_bounds": volume_bounds,
        "volume_tolerance": volume_tol,
        "input_volume_A": va,
        "input_volume_B": vb,
        "result_volume": vr,
    }
    report["timings_ms"]["verification"] = (
        perf_counter() - t_stage) * 1000.0
    if not valid:
        exc = FreeformError("final B-rep validity check failed",
                            "FinalBRepInvalid")
        refuse("verification", exc)
    if not complete_lineage:
        exc = FreeformError(
            f"final B-rep has unaudited edge lineage: "
            f"unattributed={unattributed}, bad_sections={bad_sections}",
            "IncompleteEdgeLineage")
        refuse("verification", exc)
    if not exact_curve_surface_complete:
        exc = FreeformError(
            "one or more accepted section edges lack exact bilateral "
            "curve-on-surface validation evidence",
            "ExactCurveOnSurfaceValidationMissing")
        refuse("verification", exc)
    if shadow_section_crosscheck and not shadow_complete:
        exc = FreeformError(
            "one or more accepted section edges lack requested "
            "approx/nonapprox construction cross-check evidence",
            "SectionConstructionCrosscheckMissing")
        refuse("verification", exc)
    if not volume_bounds_ok:
        exc = FreeformError(
            f"result volume {vr:.12g} violates {op} bounds "
            f"{volume_bounds} with tolerance {volume_tol:.6g}",
            "OperationVolumeInvariantFailed")
        refuse("verification", exc)

    report["accepted"] = True
    report["timings_ms"]["total"] = (perf_counter() - t_total) * 1000.0
    return _finalize(assembled.shape, report)
