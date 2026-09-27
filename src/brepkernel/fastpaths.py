"""S3: whole-operation fast paths.

Skip the section / split / assembly stages when conservative
model-level facts prove they cannot change the verdict.

Fast path A -- disjoint expanded component boxes
    Every solid of A is strictly separated from every solid of B after
    entity-derived padding.  Each solid's component box is the min/max
    over its OWN faces' prepared boxes (containment by construction, no
    reliance on OCCT box-nesting semantics) and the expansion is the max
    per-face broad-phase pad over the solid's faces -- the exact pads the
    sweep-and-prune broad phase uses.  A strictly separated component
    pair therefore implies every face pair across it is a broad-phase
    non-candidate, so no section curve, coincidence, or contact can
    exist:
      intersection -> canonical empty result
      difference   -> verified copy of A
      union        -> compound of the disjoint solids, then validity
                      and evidence checks

Fast path B -- zero broad-phase candidates, overlapping boxes
    When the face broad phase (the same candidate_face_pairs call
    intersect_models makes) finds zero candidates but component boxes
    overlap, the solids' boundaries cannot meet: each solid pair is
    disjoint or nested.  One clean interior witness per solid (>= 10x tol
    from its own boundary, confirmed inside by the G5 dual classifier)
    then decides the relation.  A UNIFORM relation resolves the operation
    without face splitting:
      A in B   : intersection -> A, union -> B, difference -> empty
      B in A   : intersection -> B, union -> A, difference -> FALLBACK
                 (a cavity needs real face splitting)
      disjoint : intersection -> empty, union -> compound, difference -> A
    Any unknown / near-boundary / disagreeing / non-uniform outcome falls
    back to the full pipeline.

Uncertainty never raises: every undecided case returns None ("no fast
path").  Unexpected exceptions are caught, counted (s3_internal_error),
and reported, so a bug in this speed code can only cost time, never
correctness (R3: every shortcut has a fallback).
"""

import numpy as np

from .perf import (
    count as _perf_count,
    S3_FAST_PATH_ATTEMPT,
    S3_FAST_PATH_HIT,
    S3_FAST_PATH_INTERNAL_ERROR,
)

# (u, v) fractions sampled per face when hunting a solid-interior witness.
_WITNESS_UV = ((0.5, 0.5), (0.25, 0.75), (0.75, 0.25))


def _solid_face_indices(model):
    """{solid_id: [face indices]} for a BRepModel."""
    groups = {}
    for i, fr in enumerate(model.faces):
        groups.setdefault(fr.solid_id, []).append(i)
    return groups


def _effective_pads(model, pads_array, scalar):
    """Per-face effective pad, mirroring step_ingest._effective_face_pads."""
    n = len(model.faces)
    if pads_array is None:
        return np.full(n, float(scalar), dtype=np.float64)
    return np.maximum(
        np.asarray(pads_array, dtype=np.float64).reshape(n), float(scalar))


def _component_boxes(prep, eff_pads):
    """Per-solid (lo, hi, max_pad) from the prepared face boxes.

    The component box is the min/max over the solid's own face boxes, so
    every expanded face box is contained in its solid's expanded
    component box by construction.
    """
    boxes = np.asarray(prep.face_boxes, dtype=np.float64)
    groups = _solid_face_indices(prep.model)
    out = []
    for sid in sorted(groups):
        idx = groups[sid]
        if not idx:
            continue
        sub = boxes[idx]
        lo = sub[:, 0:3].min(axis=0)
        hi = sub[:, 3:6].max(axis=0)
        out.append((sid, lo, hi, float(eff_pads[idx].max())))
    return out


def _all_components_separated(comp_a, comp_b):
    """True iff every (A solid, B solid) pair is strictly separated.

    Strict: expanded boxes that merely touch do NOT count as separated
    (the broad phase treats touching as a candidate via <= / >=).
    """
    if not comp_a or not comp_b:
        return False
    lo_a = np.array([c[1] for c in comp_a])
    hi_a = np.array([c[2] for c in comp_a])
    ma = np.array([c[3] for c in comp_a])
    lo_b = np.array([c[1] for c in comp_b])
    hi_b = np.array([c[2] for c in comp_b])
    mb = np.array([c[3] for c in comp_b])
    sep = np.zeros((len(comp_a), len(comp_b)), dtype=bool)
    for k in range(3):
        sep |= ((hi_a[:, k, None] + ma[:, None]
                 < lo_b[None, :, k] - mb[None, :])
                | (hi_b[None, :, k] + mb[None, :]
                   < lo_a[:, k, None] - ma[:, None]))
    return bool(sep.all())


def _dual_verdict(point, model, tol, ray):
    """G5 agreed verdict of a point vs a model, or None.

    None covers OCCT boundary/unknown states and classifier
    disagreement (which _agreed_point_verdict raises): both mean "not
    decisive", i.e. fall back to the full pipeline.
    """
    from .assembly import _agreed_point_verdict

    try:
        v = _agreed_point_verdict(np.asarray(point, dtype=np.float64),
                                  model, float(tol), ray)
    except Exception:
        return None
    return v if v in ("inside", "outside") else None


def _solid_interior_witness(prep, solid_id, *, tol, ray):
    """A clean interior witness point for one solid, or None.

    Returns (point, boundary_distance).  Samples face-interior UV
    points, offsets both ways along the surface normal by 11*tol, and
    keeps the offset farthest from the model's own boundary that (a)
    clears the 10*tol confusion band and (b) the dual classifier
    confirms as inside the model.  The dual check (not face orientation
    lore) decides which offset side is interior.
    """
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.BRepClass import BRepClass_FaceClassifier
    from OCP.gp import gp_Pnt2d
    from OCP.TopAbs import TopAbs_IN

    from .assembly import (
        _MultiRayClassifier,
        _point_boundary_distances,
    )

    model = prep.model
    k = 10.0 * float(tol)          # confusion-band requirement
    off = 11.0 * float(tol)        # witness offset (10% margin over k,
                                   # so the own face reads d = off > k
                                   # with no floating-point ambiguity)
    cap = 2.0 * off
    groups = _solid_face_indices(model)
    idx = groups.get(solid_id, [])
    best = None  # (dist, point)
    for i in idx:
        fr = model.faces[i]
        try:
            surf = BRepAdaptor_Surface(fr.face)
        except Exception:
            continue
        try:
            u0, u1, v0, v1 = (float(x) for x in fr.uv_bounds)
        except Exception:
            continue
        if not (np.isfinite([u0, u1, v0, v1]).all()
                and u1 > u0 and v1 > v0):
            continue
        for a, b in _WITNESS_UV:
            u = u0 + (u1 - u0) * a
            v = v0 + (v1 - v0) * b
            try:
                cl = BRepClass_FaceClassifier(
                    fr.face, gp_Pnt2d(float(u), float(v)),
                    float(tol), True)
                if cl.State() != TopAbs_IN:
                    continue
                n = _MultiRayClassifier._face_normal(fr.face, u, v)
            except Exception:
                continue
            if n is None:
                continue
            p = np.array([surf.Value(float(u), float(v)).X(),
                          surf.Value(float(u), float(v)).Y(),
                          surf.Value(float(u), float(v)).Z()],
                         dtype=np.float64)
            for sgn in (1.0, -1.0):
                w = p + sgn * n * off
                try:
                    d = _point_boundary_distances(
                        w.reshape(1, 3), model, cap=cap,
                        face_boxes=prep.face_boxes,
                        face_index=prep.face_index)[0]
                except Exception:
                    continue
                if not np.isfinite(d) or d < k:
                    continue
                if _dual_verdict(w, model, tol, ray) != "inside":
                    continue
                if best is None or d > best[0]:
                    best = (float(d), w)
    return None if best is None else best[::-1]  # (point, dist)


def _classifiers(prep, tol):
    """A _MultiRayClassifier over a prepared model's solids.

    Reuses the prepared deduplicated edges/boxes/index exactly as the
    assembly classifier does (S1/S2); any mismatch falls back inside
    the classifier itself.
    """
    from .assembly import _MultiRayClassifier

    return _MultiRayClassifier(
        [sr.solid for sr in prep.model.solids], float(tol),
        edges=prep.edges, edge_boxes=prep.edge_boxes,
        edge_index=prep.edge_index)


def _containment_relation(pa, pb, *, tol):
    """Classify the solid-level relation, or None when undecided.

    Returns (relation, detail) with relation in {"A_in_B", "B_in_A",
    "disjoint"}; None means fall back to the full pipeline.
    """
    sids_a = sorted(_solid_face_indices(pa.model))
    sids_b = sorted(_solid_face_indices(pb.model))
    if not sids_a or not sids_b:
        return None, {"reason": "no_solids"}
    try:
        ray_a = _classifiers(pa, tol)
        ray_b = _classifiers(pb, tol)
    except Exception:
        return None, {"reason": "classifier_build_failed"}

    wit_a, wit_b = {}, {}
    detail = {"witness_dist": {"A": {}, "B": {}}}
    for sid in sids_a:
        got = _solid_interior_witness(pa, sid, tol=tol, ray=ray_a)
        if got is None:
            return None, {"reason": f"no_clean_witness_A_solid_{sid}"}
        w, d = got
        wit_a[sid] = w
        detail["witness_dist"]["A"][str(sid)] = d
    for sid in sids_b:
        got = _solid_interior_witness(pb, sid, tol=tol, ray=ray_b)
        if got is None:
            return None, {"reason": f"no_clean_witness_B_solid_{sid}"}
        w, d = got
        wit_b[sid] = w
        detail["witness_dist"]["B"][str(sid)] = d

    v_a = {}
    for sid, w in wit_a.items():
        v = _dual_verdict(w, pb.model, tol, ray_b)
        if v is None:
            return None, {"reason": f"witness_A_solid_{sid}_undecided"}
        v_a[sid] = v
    v_b = {}
    for sid, w in wit_b.items():
        v = _dual_verdict(w, pa.model, tol, ray_a)
        if v is None:
            return None, {"reason": f"witness_B_solid_{sid}_undecided"}
        v_b[sid] = v

    a_in_b = (all(v == "inside" for v in v_a.values())
              and all(v == "outside" for v in v_b.values()))
    b_in_a = (all(v == "inside" for v in v_b.values())
              and all(v == "outside" for v in v_a.values()))
    disjoint = (all(v == "outside" for v in v_a.values())
                and all(v == "outside" for v in v_b.values()))
    if sum((a_in_b, b_in_a, disjoint)) != 1:
        return None, {"reason": "mixed_relation"}
    relation = "A_in_B" if a_in_b else ("B_in_A" if b_in_a else "disjoint")
    detail.update({
        "verdicts_A_vs_B": {str(s): v_a[s] for s in v_a},
        "verdicts_B_vs_A": {str(s): v_b[s] for s in v_b},
    })
    return relation, detail


def _copy_shape(shape):
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Copy

    return BRepBuilderAPI_Copy(shape).Shape()


def _empty_result():
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Compound

    c = TopoDS_Compound()
    BRep_Builder().MakeCompound(c)
    return c


def _disjoint_union_shape(pa, pb):
    """Compound of copies of every solid of A and B."""
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Compound

    builder = BRep_Builder()
    comp = TopoDS_Compound()
    builder.MakeCompound(comp)
    for prep in (pa, pb):
        for sr in prep.model.solids:
            builder.Add(comp, _copy_shape(sr.solid))
    return comp


def _resolve(pa, pb, op, path, relation):
    """Map (path, relation, op) to a result shape, or None to fall back."""
    if path == "A_disjoint":
        if op == "intersection":
            return _empty_result(), "empty"
        if op == "difference":
            return _copy_shape(pa.model.shape), "A"
        return _disjoint_union_shape(pa, pb), "union"
    # path == "B_containment"
    if relation == "A_in_B":
        if op == "intersection":
            return _copy_shape(pa.model.shape), "A"
        if op == "difference":
            return _empty_result(), "empty"
        return _copy_shape(pb.model.shape), "B"
    if relation == "B_in_A":
        if op == "intersection":
            return _copy_shape(pb.model.shape), "B"
        if op == "difference":
            return None, None  # cavity: needs real face splitting
        return _copy_shape(pa.model.shape), "A"
    # relation == "disjoint" with overlapping component boxes
    if op == "intersection":
        return _empty_result(), "empty"
    if op == "difference":
        return _copy_shape(pa.model.shape), "A"
    return _disjoint_union_shape(pa, pb), "union"


def try_fast_path(pa, pb, op, *, base_tol, contact_tol,
                  pad_scalar, pads):
    """Attempt an S3 whole-operation fast path.

    Returns (shape, info, candidates); shape None means "no fast path,
    run the full pipeline".  `candidates` is the broad-phase candidate
    list computed during the probe (None when the probe never ran the
    broad phase, e.g. path A fired): the caller should hand it to
    intersect_models as precomputed_candidates so the broad phase is
    not run twice.  `pad_scalar` / `pads` must be exactly what the
    caller will pass to intersect_models (broadphase_pad /
    broadphase_face_pads) so the candidate set agrees with the pipeline.
    `contact_tol` sets the witness confusion band (10x) and the dual
    classifier tolerance.
    """
    from .step_ingest import candidate_face_pairs

    _perf_count(S3_FAST_PATH_ATTEMPT)
    info = {"attempted": True, "path": None, "reason": None,
            "resolution": None}
    try:
        return _try_fast_path_inner(
            pa, pb, op, info, base_tol=base_tol,
            contact_tol=contact_tol, pad_scalar=pad_scalar, pads=pads,
            candidate_face_pairs=candidate_face_pairs)
    except Exception as exc:  # R3: a speed bug costs time, never correctness
        _perf_count(S3_FAST_PATH_INTERNAL_ERROR)
        info["reason"] = "s3_internal_error"
        info["s3_internal_error"] = f"{type(exc).__name__}: {exc}"
        return None, info, None


def _try_fast_path_inner(pa, pb, op, info, *, base_tol, contact_tol,
                         pad_scalar, pads, candidate_face_pairs):
    tol = float(contact_tol)
    pads_a, pads_b = pads if pads is not None else (None, None)
    eff_a = _effective_pads(pa.model, pads_a, pad_scalar)
    eff_b = _effective_pads(pb.model, pads_b, pad_scalar)
    comp_a = _component_boxes(pa, eff_a)
    comp_b = _component_boxes(pb, eff_b)
    info["component_pairs_checked"] = len(comp_a) * len(comp_b)

    # Fast path A: every component pair strictly separated.
    if _all_components_separated(comp_a, comp_b):
        shape, resolution = _resolve(pa, pb, op, "A_disjoint", None)
        _perf_count(S3_FAST_PATH_HIT)
        info.update({"path": "A_disjoint", "resolution": resolution})
        return shape, info, None

    # Fast path B: no broad-phase candidate at all, then containment.
    # The candidate list is returned so the caller can hand it to
    # intersect_models instead of running the broad phase twice.
    candidates = candidate_face_pairs(
        pa.model, pb.model, pad=float(pad_scalar),
        pads_a=pads_a, pads_b=pads_b)
    info["broadphase_candidates"] = len(candidates)
    if candidates:
        info["reason"] = "candidates_exist"
        return None, info, candidates
    relation, detail = _containment_relation(pa, pb, tol=tol)
    info.update(detail)
    if relation is None:
        info["reason"] = detail.get("reason", "containment_undecided")
        return None, info, candidates
    shape, resolution = _resolve(pa, pb, op, "B_containment", relation)
    if shape is None:
        info["reason"] = "unsupported_difference_cavity"
        return None, info, candidates
    _perf_count(S3_FAST_PATH_HIT)
    info.update({"path": "B_containment", "relation": relation,
                 "resolution": resolution})
    return shape, info, candidates
