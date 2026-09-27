"""S7: boolean_brep_many -- multi-tool Booleans.

Phase 1 (reference): plain sequential loop, each tool applied to the
previous result with ``boolean_brep``.  This is the semantics every
optimized path must reproduce exactly.

Phase 2 (optimized, ``optimized=True``): instead of N sequential
full-pipeline ops -- each re-splitting and re-assembling the growing
intermediate solid -- combine the tools once with the dual operation
and run a single pipeline op.  The identities are exact set theory,
not approximations::

    B - t1 - t2 - ... - tn  =  B - (t1 ∪ t2 ∪ ... ∪ tn)
    B ∪ t1 ∪ t2 ∪ ... ∪ tn  =  B ∪ (t1 ∪ t2 ∪ ... ∪ tn)

The per-tool section work is unchanged; what disappears is the
repeated split / classify / assembly / verification over the whole
solid on every step.  ``intersection`` is intentionally *not*
optimized: folding tools first would feed possibly-empty intermediate
shapes back into the pipeline, whose refusal behavior differs from the
sequential reference (probed 27 Sept 2026).  Intersection falls back to
the Phase 1 loop and reports ``optimized: False``.

Conservative rules (plan R1/R2/R3):
- A tool is skipped (difference only) only when its padded boxes prove
  it cannot touch any base face box.  The pad arithmetic mirrors the
  pipeline's own broad phase (``face_pads(contact_tol)``), so a skip
  fires exactly when the pipeline itself would find zero candidates.
- Tools are merged into a free compound only when their padded overall
  boxes are pairwise disjoint -- then the compound *is* the union.
- Tools whose padded boxes overlap are combined with real, certified
  ``boolean_brep`` union ops (small inputs, cheap).
- Every shortcut is observable in the returned summary (R4): skip
  counts, group counts, combine-op reports, and per-tool ancestry.

Provenance (pass criterion): the summary carries ``tool_ancestry`` --
one entry per input tool with its bbox, skip/group status -- plus
``face_attribution``, a per-result-face tool index (or -1 for
base-derived / ambiguous faces), so every result edge's ancestry is
identifiable by the faces that own it.
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# public API
# ---------------------------------------------------------------------------

def boolean_brep_many(base, tools, op="difference", *, optimized=False,
                      base_tol=1e-7, shadow_reference=False, **kw):
    """Apply a Boolean op between *base* and each of *tools* in order.

    Phase 1 (``optimized=False``, default): sequential reference --
    ``out = boolean_brep(out, tool, op)`` per tool.  Returns
    ``(final_shape, {"steps": [reports...], "optimized": False})``.

    Phase 2 (``optimized=True``): for ``difference``/``union`` the tools
    are combined once (free compound when pairwise disjoint, certified
    boolean unions otherwise) and a single pipeline op runs against the
    prepared base.  ``intersection`` falls back to the sequential loop.

    With ``shadow_reference=True`` (and ``optimized=True``) the
    sequential reference is also run; any mismatch in acceptance,
    volume, or face/solid counts raises ``AssertionError``, and the
    *reference* result is served (same pattern as the S3 shadow).
    """
    if op == "intersection" and optimized:
        # Folding tools first would feed possibly-empty intermediates
        # into the pipeline, whose refusal behavior differs from the
        # sequential reference.  Stay on the reference path (R1).
        optimized = False
    if not optimized:
        out, summary = _sequential_reference(base, tools, op,
                                             base_tol=base_tol, **kw)
        summary["shadow_reference"] = False
        return out, summary
    out, summary = _optimized_many(base, tools, op, base_tol=base_tol,
                                   shadow_reference=shadow_reference, **kw)
    if shadow_reference:
        ref_out, ref_summary = _sequential_reference(base, tools, op,
                                                     base_tol=base_tol, **kw)
        _assert_equivalent(summary, ref_summary, ref_out, out)
        summary["shadow_reference"] = True
        summary["shadow_matched"] = True
        return ref_out, summary
    summary["shadow_reference"] = False
    return out, summary


def _sequential_reference(base, tools, op, *, base_tol=1e-7, **kw):
    """S7 Phase 1: the exact sequential semantics (moved from session)."""
    from .pipeline import boolean_brep

    out = base
    reports = []
    tools = list(tools)
    for tool in tools:
        out, report = boolean_brep(out, tool, op, base_tol=base_tol, **kw)
        reports.append(report)
    return out, {"steps": reports, "optimized": False, "op": op,
                 "n_tools": len(tools)}


# ---------------------------------------------------------------------------
# Phase 2
# ---------------------------------------------------------------------------

def _optimized_many(base, tools, op, *, base_tol=1e-7,
                    shadow_reference=False, **kw):
    from .pipeline import boolean_brep
    from .prepared import ensure_prepared

    tools = list(tools)
    contact_tol = float(kw.get("contact_tol", 4.0 * float(base_tol)))
    prep_base = ensure_prepared(base, base_tol=base_tol)
    summary = {"optimized": True, "op": op, "n_tools": len(tools),
               "n_skipped": 0, "n_groups": 0, "n_combine_ops": 0,
               "combine_reports": [], "tool_ancestry": [],
               "face_attribution": []}

    if not tools:
        return prep_base.model.shape, {**summary, "steps": []}

    prep_tools = [ensure_prepared(t, base_tol=base_tol) for t in tools]

    # 1. Batched broad phase: which tools can touch the base at all?
    #    (difference only -- union needs every tool in the result.)
    contacting = _batched_contact(prep_base, prep_tools, contact_tol)
    if op == "difference":
        kept_idx = [i for i, c in enumerate(contacting) if c]
        summary["n_skipped"] = len(tools) - len(kept_idx)
    else:
        kept_idx = list(range(len(tools)))
    if not kept_idx:
        # Every tool provably misses the base: difference is a no-op.
        for i in range(len(tools)):
            summary["tool_ancestry"].append(_ancestry_entry(
                i, prep_tools[i], contact_tol, skipped=True, group_id=-1))
        return prep_base.model.shape, {**summary, "steps": []}
    kept = [prep_tools[i] for i in kept_idx]

    # 2. Partition kept tools into interaction clusters: connected
    #    components of the padded-box overlap graph.
    boxes = np.array([_padded_overall_box(p, contact_tol) for p in kept])
    groups = _overlap_groups(boxes)
    summary["n_groups"] = len(groups)

    # 3. Combine within each group with certified boolean unions
    #    (single-tool groups need no work).
    group_shapes = []
    for g in groups:
        if len(g) == 1:
            group_shapes.append(kept[g[0]].model.shape)
            continue
        acc = kept[g[0]]
        for j in g[1:]:
            acc_shape, rep = boolean_brep(acc, kept[j], "union",
                                          base_tol=base_tol, **kw)
            summary["combine_reports"].append(rep)
            summary["n_combine_ops"] += 1
            acc = ensure_prepared(acc_shape, base_tol=base_tol)
        group_shapes.append(acc.model.shape
                            if hasattr(acc, "model") else acc)

    # 4. Groups are pairwise disjoint, so compounding them is exactly
    #    their union -- one batched broad phase, one split, one assembly.
    combined = _make_compound(group_shapes) if len(group_shapes) > 1 \
        else group_shapes[0]

    # 5. The single main op, on the prepared base.
    out, main_rep = boolean_brep(prep_base, combined, op,
                                 base_tol=base_tol, **kw)

    for rank, i in enumerate(kept_idx):
        summary["tool_ancestry"].append(_ancestry_entry(
            i, prep_tools[i], contact_tol, skipped=False,
            group_id=_group_of(groups, rank)))
    for i in range(len(tools)):
        if i not in kept_idx:
            summary["tool_ancestry"].append(_ancestry_entry(
                i, prep_tools[i], contact_tol, skipped=True, group_id=-1))
    summary["tool_ancestry"].sort(key=lambda e: e["tool_index"])
    summary["face_attribution"] = _attribute_faces(
        out, [ _padded_overall_box(prep_tools[i], contact_tol)
               for i in kept_idx ])
    summary["steps"] = summary["combine_reports"] + [main_rep]
    return out, summary


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _overall_box(prep) -> np.ndarray:
    """Overall (6,) AABB of a prepared shape from its face boxes."""
    fb = np.asarray(prep.face_boxes, dtype=np.float64).reshape(-1, 6)
    if len(fb) == 0:
        z = np.zeros(6)
        return z
    lo = fb[:, :3].min(axis=0)
    hi = fb[:, 3:].max(axis=0)
    return np.concatenate([lo, hi])


def _padded_overall_box(prep, contact_tol: float) -> np.ndarray:
    """Overall box padded by the tool's own broad-phase pad (conservative:

    the pad covers every per-face pad, so overlap of padded overall boxes
    is a superset of any per-face-box overlap).
    """
    b = _overall_box(prep)
    ft = np.asarray(prep.face_tol, dtype=np.float64).reshape(-1)
    pad = float(contact_tol) + (float(ft.max()) if len(ft) else 0.0)
    return b + np.array([-pad, -pad, -pad, pad, pad, pad])


def _batched_contact(prep_base, prep_tools, contact_tol) -> np.ndarray:
    """One batched query: which tools can touch the base?

    Mirrors the pipeline's broad-phase pad arithmetic exactly
    (``face_pads(contact_tol)`` on both sides); a ``False`` here means
    the pipeline itself would find zero candidate pairs for that tool.
    """
    fb = np.asarray(prep_base.face_boxes, dtype=np.float64).reshape(-1, 6)
    if len(fb) == 0 or not prep_tools:
        return np.zeros(len(prep_tools), dtype=bool)
    pads = np.asarray(prep_base.face_pads(contact_tol),
                      dtype=np.float64).reshape(-1)
    fb_p = fb + np.stack([-pads, -pads, -pads, pads, pads, pads],
                         axis=1)
    tb = np.array([_padded_overall_box(p, contact_tol)
                   for p in prep_tools])  # (m, 6)
    over = ((tb[:, None, 0] <= fb_p[None, :, 3]) &
            (tb[:, None, 3] >= fb_p[None, :, 0]) &
            (tb[:, None, 1] <= fb_p[None, :, 4]) &
            (tb[:, None, 4] >= fb_p[None, :, 1]) &
            (tb[:, None, 2] <= fb_p[None, :, 5]) &
            (tb[:, None, 5] >= fb_p[None, :, 2]))
    return over.any(axis=1)


def _boxes_overlap(a: np.ndarray, b: np.ndarray) -> bool:
    return bool(a[0] <= b[3] and a[3] >= b[0] and
                a[1] <= b[4] and a[4] >= b[1] and
                a[2] <= b[5] and a[5] >= b[2])


def _overlap_groups(boxes: np.ndarray):
    """Connected components of the box-overlap graph (union-find)."""
    m = len(boxes)
    parent = list(range(m))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for i in range(m):
        for j in range(i + 1, m):
            if _boxes_overlap(boxes[i], boxes[j]):
                ri, rj = find(i), find(j)
                if ri != rj:
                    parent[ri] = rj
    comps = {}
    for i in range(m):
        comps.setdefault(find(i), []).append(i)
    return [sorted(v) for v in comps.values()]


def _group_of(groups, rank: int) -> int:
    for gi, g in enumerate(groups):
        if rank in g:
            return gi
    return -1


def _make_compound(shapes):
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Compound
    comp = TopoDS_Compound()
    bld = BRep_Builder()
    bld.MakeCompound(comp)
    for s in shapes:
        bld.Add(comp, s)
    return comp


def _ancestry_entry(i, prep, contact_tol, *, skipped, group_id):
    b = _overall_box(prep)
    return {"tool_index": i, "bbox": [float(v) for v in b],
            "skipped": bool(skipped), "group_id": int(group_id)}


def _face_box(face) -> np.ndarray:
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib
    b = Bnd_Box()
    BRepBndLib.AddOptimal_s(face, b, False, True)
    if b.IsVoid():
        z = np.zeros(6)
        return z
    lo = b.CornerMin()
    hi = b.CornerMax()
    return np.array([lo.X(), lo.Y(), lo.Z(), hi.X(), hi.Y(), hi.Z()])


def _attribute_faces(result_shape, tool_boxes):
    """Per-result-face tool attribution.

    Returns a list aligned with TopExp_Explorer face order: the index
    into the kept-tool list whose padded box uniquely contains the face
    box, or -1 for base-derived / ambiguous faces.  With pairwise
    disjoint tools the attribution is unambiguous: a tool-derived face
    must lie inside its own tool's box.
    """
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopAbs import TopAbs_FACE
    attr = []
    ex = TopExp_Explorer(result_shape, TopAbs_FACE)
    while ex.More():
        fb = _face_box(ex.Current())
        hits = [i for i, tb in enumerate(tool_boxes)
                if _boxes_overlap(tb, fb)
                and tb[0] <= fb[0] and tb[3] >= fb[3]
                and tb[1] <= fb[1] and tb[4] >= fb[4]
                and tb[2] <= fb[2] and tb[5] >= fb[5]]
        attr.append(hits[0] if len(hits) == 1 else -1)
        ex.Next()
    return attr


def _assert_equivalent(opt_summary, ref_summary, ref_out, opt_out):
    """R5 shadow check: optimized must match the reference exactly."""
    from .assembly import _shape_volume
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID

    def count(shape, what):
        ex = TopExp_Explorer(shape, what)
        n = 0
        while ex.More():
            n += 1
            ex.Next()
        return n

    ro = ref_summary["steps"][-1] if ref_summary["steps"] else None
    oo = opt_summary["steps"][-1] if opt_summary["steps"] else None
    ra = all(r["accepted"] for r in ref_summary["steps"])
    oa = all(r["accepted"] for r in opt_summary["steps"])
    assert ra == oa, f"acceptance mismatch: ref {ra} vs opt {oa}"
    if ra:
        rv = abs(float(_shape_volume(ref_out)))
        ov = abs(float(_shape_volume(opt_out)))
        # NOTE: exact bit equality cannot hold here.  The two paths run
        # different op sequences, so the volume integral sums the same
        # face contributions in a different order (1-ulp noise, ~1e-15
        # relative, observed 27 Sept 2026).  A 1e-9 relative band admits
        # only float noise, never a real geometric change (R1).
        assert abs(rv - ov) <= 1e-9 * max(1.0, abs(rv)), \
            f"volume mismatch: ref {rv} vs opt {ov}"
        assert count(ref_out, TopAbs_FACE) == count(opt_out, TopAbs_FACE), \
            "face count mismatch"
        assert count(ref_out, TopAbs_SOLID) == count(opt_out, TopAbs_SOLID), \
            "solid count mismatch"
