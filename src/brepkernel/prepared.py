"""S1: PreparedBRep - immutable prepared model owning all acceleration data.

A PreparedBRep wraps one BRepModel and precomputes, once, the data that
the pipeline stages otherwise rebuild on every Boolean call:

- face_boxes: (n, 6) conservative AABBs built by the exact same
  ``assembly._conservative_boxes`` the C6 distance path uses, so the
  boxes are bit-identical to a per-call rebuild;
- edges / edge_boxes: the deduplicated union of every solid's unique
  edges (IsSame identity) plus their conservative boxes, built by the
  same helpers ``_MultiRayClassifier`` uses;
- face_tol: per-face tol_face array; per-call broad-phase pads are
  ``contact_tol + face_tol`` (same arithmetic as
  ``step_ingest.face_broadphase_pads``);
- face_adjacency: index-based adjacency (faces sharing an IsSame edge),
  built with the C7 hash-bucket pattern;
- solid_boxes, analytic_mask / freeform_mask, all_faces_analytic,
  max_tolerance, topological counts.

Lifetime (plan rule R6): a PreparedBRep is immutable (frozen dataclass;
every numpy array is an owned copy marked read-only at construction)
and is valid for as long as the caller holds it. ``boolean_brep``
accepts a raw shape, a BRepModel, or a PreparedBRep; raw shapes and
plain BRepModels are prepared internally, so existing callers are
unaffected.

Explicit contract: the underlying TopoDS_Shape is shared by reference
and OCCT shapes are mutable in principle (a caller could move faces or
rewrite tolerances through handles).  A PreparedBRep assumes nobody
mutates the shape's topology, geometry, or tolerances after
preparation; violating that invalidates every cached box, index, and
mask.  Preparation itself never mutates the input shape.

Deliberately NOT stored here (per the plan): cached classification
verdicts, tolerance-dependent section results, and mutable OCCT
algorithm objects. Those belong to a per-call query context (S4).

Preparation performs no tessellation and does not mutate the input
shape: every builder below uses explorers, BRepBndLib box queries with
triangulation disabled, and BRep_Tool tolerance/surface-type reads.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Union

import numpy as np

from . import perf as _perf
from .spatial import HybridBoxIndex


@dataclass(frozen=True)
class PreparedBRep:
    """Immutable acceleration data for one indexed shape.

    All numpy arrays are owned, C-contiguous, read-only copies
    (write=False); the frozen dataclass additionally blocks attribute
    reassignment.  The wrapped ``model`` / TopoDS_Shape is shared by
    reference and must not be mutated after preparation (see the
    module docstring contract).
    """

    model: object                      # BRepModel from index_shape()
    face_boxes: np.ndarray              # (n, 6) conservative AABBs
    edges: tuple                        # unique edges (deduped by IsSame)
    edge_boxes: np.ndarray              # (m, 6) conservative AABBs
    face_tol: np.ndarray                # (n,) per-face tol_face
    face_adjacency: tuple               # tuple[tuple[int, ...]] by face index
    solid_boxes: np.ndarray             # (s, 6) conservative AABBs
    analytic_mask: np.ndarray           # (n,) bool, True = analytic surface
    freeform_mask: np.ndarray           # (n,) bool, True = freeform/NURBS
    all_faces_analytic: bool
    max_tolerance: float
    face_index: HybridBoxIndex             # S2: proximity index over face_boxes
    edge_index: HybridBoxIndex             # S2: proximity index over edge_boxes
    base_tol: float = 1e-7

    @property
    def n_faces(self) -> int:
        return len(self.model.faces)

    @property
    def n_edges(self) -> int:
        return len(self.edges)

    def face_pads(self, contact_tol: float) -> np.ndarray:
        """Per-face broad-phase pads: contact_tol + tol_face.

        Identical arithmetic to step_ingest.face_broadphase_pads; the
        contact_tol term is per-call, so it is added here rather than
        frozen at prepare time.
        """
        return np.asarray(self.face_tol, dtype=np.float64) + float(contact_tol)


def _freeze_array(a, dtype=None):
    """Owned, C-contiguous, read-only copy of a numeric array.

    frozen=True on the dataclass stops attribute reassignment but not
    in-place mutation (prepared.face_boxes[0, 0] = ...).  The prepared
    arrays back correctness arguments (and HybridBoxIndex builds BVH
    node bounds from the same boxes), so every array stored on the
    dataclass is frozen here at construction.
    """
    b = np.array(a, dtype=dtype, copy=True, order="C")
    b.setflags(write=False)
    return b


def _import_assembly():
    # Lazy: assembly.py must not import prepared.py at module level
    # (consumers there only duck-type on PreparedBRep), so the    # assembly -> prepared direction stays import-cycle free.
    from . import assembly as _a
    return _a


def _build_face_boxes(model) -> np.ndarray:
    asm = _import_assembly()
    faces = [fr.face for fr in model.faces]
    if not faces:
        return np.empty((0, 6), dtype=np.float64)
    _perf.count(_perf.FACE_BOX_BUILD)
    return asm._conservative_boxes(faces)


def _build_unique_edges(model) -> tuple:
    """Deduplicated union of every solid's unique edges.

    Built per solid (not from the whole shape) so the edge SET is
    exactly what _MultiRayClassifier sees when it explores the model's
    solids: substituting the prepared list cannot change its min-distance
    verdicts. Hash buckets keep the IsSame dedup linear; first-seen
    explorer order is preserved.
    """
    asm = _import_assembly()
    buckets: dict = {}
    out: list = []
    for sr in model.solids:
        for e in asm._unique_edges(sr.solid):
            bucket = buckets.setdefault(hash(e), [])
            if any(e.IsSame(x) for x in bucket):
                continue
            bucket.append(e)
            out.append(e)
    return tuple(out)


def _build_edge_boxes(edges: tuple) -> np.ndarray:
    asm = _import_assembly()
    if not edges:
        return np.empty((0, 6), dtype=np.float64)
    _perf.count(_perf.EDGE_BOX_BUILD)
    return asm._conservative_boxes(list(edges))


def _build_face_adjacency(model) -> tuple:
    """Index-based face adjacency via shared IsSame edges.

    adj[i] lists the indices j != i whose faces share at least one edge
    with face i. Symmetric by construction.
    """
    from OCP.TopAbs import TopAbs_EDGE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS

    owners: dict = {}
    for i, fr in enumerate(model.faces):
        ex = TopExp_Explorer(fr.face, TopAbs_EDGE)
        while ex.More():
            e = TopoDS.Edge(ex.Current())
            ex.Next()
            owners.setdefault(hash(e), []).append((e, i))
    adj = [set() for _ in model.faces]
    for bucket in owners.values():
        for k in range(len(bucket)):
            ea, ia = bucket[k]
            for m in range(k + 1, len(bucket)):
                eb, ib = bucket[m]
                if ia != ib and ea.IsSame(eb):
                    adj[ia].add(ib)
                    adj[ib].add(ia)
    return tuple(tuple(sorted(s)) for s in adj)


def _build_solid_boxes(model) -> np.ndarray:
    asm = _import_assembly()
    solids = [sr.solid for sr in model.solids]
    if not solids:
        return np.empty((0, 6), dtype=np.float64)
    return asm._conservative_boxes(solids)


def prepare_brep(shape, *, base_tol: float = 1e-7) -> PreparedBRep:
    """Index a raw TopoDS_Shape and precompute all acceleration data."""
    from .step_ingest import index_shape

    model = index_shape(shape)
    return _prepare_from_model(model, base_tol=float(base_tol))


def prepare_model(model, *, base_tol: float = 1e-7) -> PreparedBRep:
    """Precompute acceleration data for an existing BRepModel.

    The model itself is reused as-is (no re-indexing); the PreparedBRep
    shares it by reference.
    """
    return _prepare_from_model(model, base_tol=float(base_tol))


def _prepare_from_model(model, *, base_tol: float) -> PreparedBRep:
    from .step_ingest import model_max_tolerance

    asm = _import_assembly()
    edges = _build_unique_edges(model)
    face_boxes = _freeze_array(_build_face_boxes(model), dtype=np.float64)
    edge_boxes = _freeze_array(_build_edge_boxes(edges), dtype=np.float64)
    surface_kinds = [str(fr.surface_type) for fr in model.faces]
    # Canonical analytic predicate: the same _ANALYTIC_SURFACES set
    # assembly._all_faces_analytic uses (plane, cylinder, cone, sphere,
    # torus).  A "BSpline absent" test would wrongly mark Bezier,
    # offset, extrusion, and revolution surfaces as analytic.
    _analytic = asm._ANALYTIC_SURFACES
    analytic_mask = _freeze_array(
        [(k or "").split(".")[-1] in _analytic for k in surface_kinds],
        dtype=bool)
    freeform_mask = _freeze_array(
        [fr.freeform is not None for fr in model.faces], dtype=bool)
    return PreparedBRep(
        model=model,
        face_boxes=face_boxes,
        edges=edges,
        edge_boxes=edge_boxes,
        face_tol=_freeze_array(
            [float(fr.tol_face) for fr in model.faces], dtype=np.float64),
        face_adjacency=_build_face_adjacency(model),
        solid_boxes=_freeze_array(_build_solid_boxes(model),
                                  dtype=np.float64),
        analytic_mask=analytic_mask,
        freeform_mask=freeform_mask,
        all_faces_analytic=asm._all_faces_analytic(model.shape),
        max_tolerance=float(model_max_tolerance(model)),
        base_tol=float(base_tol),
        face_index=HybridBoxIndex(face_boxes),
        edge_index=HybridBoxIndex(edge_boxes),
    )


def ensure_prepared(shape_or_model, *, base_tol: float = 1e-7) -> PreparedBRep:
    """Return a PreparedBRep for a raw shape, BRepModel, or PreparedBRep.

    PreparedBRep passes through untouched; a BRepModel is wrapped without
    re-indexing; anything else is indexed via prepare_brep.
    """
    from .step_ingest import BRepModel

    if isinstance(shape_or_model, PreparedBRep):
        return shape_or_model
    if isinstance(shape_or_model, BRepModel):
        return prepare_model(shape_or_model, base_tol=base_tol)
    return prepare_brep(shape_or_model, base_tol=base_tol)
