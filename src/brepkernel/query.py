"""S4: per-Boolean-call shared query context.

This module consolidates every cache that is valid for exactly one
operand pair and one tolerance configuration.  It replaces two pieces
of hidden cross-call state:

* the module-global ``_VOLUME_CACHE`` dict in :mod:`brepkernel.assembly`
  (an R6 violation: two simultaneous ``boolean_brep`` calls could see
  each other's entries), and
* the ``_classify_pieces.region_stats`` function attribute,

and adds bounded per-call memos for repeated exact queries (point
verdicts, boundary distances) plus one shared multi-ray intersector per
operand/tolerance bucket.

Lifetime and threading
----------------------
A :class:`QueryContext` is created at the :func:`boolean_brep` public
boundary and threaded explicitly through ``_boolean_brep_impl``,
``assemble_boolean`` and ``_classify_pieces``.  The S3 shadow re-run
shares the outer call's context (same operands, same tolerances), so
its classifiers, volumes and point verdicts are not rebuilt.

Every memo is a pure function of its key within one call: the shapes
are never mutated (S1 immutability contract), so a cached answer is
always the answer a fresh computation would give.  A missing entry, a
key mismatch, or ``ctx=None`` falls back to the uncached computation;
nothing here can change a verdict, only skip recomputation.

Observability (R4): cache hits are counted with the S0 perf counters
(``ctx_*`` names) so the report shows what the context saved.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Hashable, Optional, Tuple


def _point_key(p) -> Tuple[float, float, float]:
    """Stable hashable key for a witness point (1 nm rounding)."""
    return (round(float(p[0]), 9), round(float(p[1]), 9), round(float(p[2]), 9))


def _conservative_model_bbox(model) -> Optional[Tuple[float, float, float,
                                                     float, float, float]]:
    """Conservative AABB of a BRepModel's solids, or None when it has none.

    Uses BRepBndLib.Add (not AddOptimal): the box always contains the
    shape, enlarged by entity tolerances.  A point strictly outside this
    box is strictly outside every solid of the model -- the geometric
    fact the SX point-verdict short-circuit rests on.
    """
    solids = getattr(model, "solids", None)
    if not solids:
        return None
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib
    b = Bnd_Box()
    for sr in solids:
        BRepBndLib.Add_s(sr.solid, b, False)
    if b.IsVoid():
        return None
    lo, hi = b.CornerMin(), b.CornerMax()
    return (float(lo.X()), float(lo.Y()), float(lo.Z()),
            float(hi.X()), float(hi.Y()), float(hi.Z()))


@dataclass
class QueryContext:
    """Shared per-Boolean-call query state (S4).

    All dicts below live exactly for one top-level ``boolean_brep``
    call.  Two simultaneous calls get separate contexts and cannot see
    each other's entries.
    """

    base_tol: float = 1e-7
    prepared_a: Any = None
    prepared_b: Any = None

    # (id(shape), tol) -> (shape, volume).  The strong shape reference
    # guards against id() reuse after GC: the dict keeps the original
    # wrapper alive, and volume_get verifies `stored is shape` before
    # trusting a hit.  (The old global cache stored only the int id.)
    #
    # CONTRACT: a cached volume is valid only while the shape is not
    # modified in place afterwards.  Callers must cache only shapes whose
    # geometry is finalized -- input solids (TopoDS non-mutation
    # contract) and freshly built result solids.  Never cache a shape
    # that will be mutated later in the call.
    volumes: Dict[tuple, tuple] = field(default_factory=dict)

    # (point_key, id(model), tol) -> agreed "inside"/"outside"/...
    # verdict from _agreed_point_verdict (dual classifiers).
    point_verdicts: Dict[tuple, str] = field(default_factory=dict)

    # (point_key, id(model), cap) -> boundary distance (capped).
    boundary_distances: Dict[tuple, float] = field(default_factory=dict)

    # id(model) -> conservative (xmin, ymin, zmin, xmax, ymax, zmax) of
    # the model's solids (SX).  The models are never mutated during a
    # call (S1 immutability contract), so the box is valid for the whole
    # call; a missing entry falls back to computing it.
    model_bboxes: Dict[int, tuple] = field(default_factory=dict)

    # (side, tol) -> _MultiRayClassifier.  One intersector set per
    # operand and tolerance bucket; built lazily via get_or_create.
    classifiers: Dict[tuple, Any] = field(default_factory=dict)

    # Replaces the _classify_pieces.region_stats function attribute.
    region_stats: Dict[str, Any] = field(default_factory=dict)

    # Candidate face ids from the broad phase, stored for
    # observability once known.
    candidate_face_ids_a: Optional[frozenset] = None
    candidate_face_ids_b: Optional[frozenset] = None

    def new_call(self) -> "QueryContext":
        """Return a fresh context for a new Boolean call.

        Carries over the tolerance policy and prepared references
        (they describe the inputs, not the query), but drops every
        memo, classifier, and region stat.  Used when one public call
        fans out into logically separate Boolean operations.
        """
        return QueryContext(
            base_tol=self.base_tol,
            prepared_a=self.prepared_a,
            prepared_b=self.prepared_b,
        )

    def get_or_create(self, key: Hashable, factory: Callable[[], Any]) -> Any:
        """Return the cached object for `key`, building it on first use."""
        try:
            return self.classifiers[key]
        except KeyError:
            pass
        obj = factory()
        self.classifiers[key] = obj
        return obj

    # -- volume cache (replaces assembly._VOLUME_CACHE) -----------------

    def volume_get(self, shape, tol: float) -> Optional[float]:
        """Cached volume, or None on miss.  Never raises.

        A hit requires the stored shape object to *be* the queried shape
        (identity, not just id equality), so a recycled id() can never
        serve a stale volume.
        """
        hit = self.volumes.get((id(shape), float(tol)))
        if hit is not None and hit[0] is shape:
            return hit[1]
        return None

    def volume_put(self, shape, tol: float, vol: float) -> None:
        """Cache a volume.  Only call for finalized shapes (see CONTRACT
        on `volumes` above)."""
        self.volumes[(id(shape), float(tol))] = (shape, float(vol))

    # -- point verdict memo -------------------------------------------

    def point_verdict_get(self, point, model, tol: float) -> Optional[str]:
        return self.point_verdicts.get(
            (_point_key(point), id(model), float(tol)))

    def point_verdict_put(self, point, model, tol: float, verdict: str) -> None:
        self.point_verdicts[(_point_key(point), id(model), float(tol))] = verdict

    # -- boundary distance memo ---------------------------------------

    def boundary_distance_get(self, point, model, cap: float) -> Optional[float]:
        return self.boundary_distances.get(
            (_point_key(point), id(model), float(cap)))

    def boundary_distance_put(self, point, model, cap: float, d: float) -> None:
        self.boundary_distances[(_point_key(point), id(model), float(cap))] = float(d)

    # -- model bbox (SX point-verdict short-circuit) --------------------

    def model_bbox(self, model) -> Optional[Tuple[float, float, float,
                                                 float, float, float]]:
        """Conservative bbox of the model's solids, cached for the call.

        Never raises: returns None when the model has no solids (the
        caller then falls back to the uncached full path, which raises
        the usual AssemblyError for shell-only input).
        """
        key = id(model)
        hit = self.model_bboxes.get(key)
        if hit is not None:
            return hit
        bb = _conservative_model_bbox(model)
        if bb is not None:
            self.model_bboxes[key] = bb
        return bb
