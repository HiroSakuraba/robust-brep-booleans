"""S0: stage and expensive-call profiler (zero cost when disabled).

This module owns the optional performance counters for the exact
trimmed-B-rep pipeline. The design follows the speed plan's rules:

* Verdict preservation: counting never changes control flow. Every
  instrumented site calls :func:`count`, which is a no-op unless a
  :class:`PerfCounters` is active for the current Boolean call.
* Zero cost when disabled: :func:`count` is one ``ContextVar.get()``
  plus a ``None`` check. No counters object is allocated and no dict
  is touched unless profiling was requested.
* Explicit bounded lifetime (R6): counters live exactly for one
  ``boolean_brep`` call. They are installed with :func:`scoped` at the
  public API boundary and reset on exit, so nested or concurrent calls
  cannot see each other's counts. No global cache keyed by object id.
* Observability (R4): when enabled, ``report["performance"]`` carries
  the counter snapshot. When disabled, the report has no
  ``"performance"`` key at all, so canonical evidence records are
  byte-identical between perf-on and perf-off runs apart from the
  intentionally added diagnostics section.

Enable per call with ``boolean_brep(..., collect_perf=True)`` or for
the whole process with the ``BREPKERNEL_PERF=1`` environment variable.
"""

from __future__ import annotations

import contextvars
import os
from typing import Dict, Optional

# Canonical counter names. Instrumentation sites must use these so that
# budgets and reports can rely on stable keys.
EXACT_FACE_DISTANCE = "exact_face_distance"      # BRepExtrema vertex-vs-face
EXACT_EDGE_DISTANCE = "exact_edge_distance"      # BRepExtrema vertex-vs-edge
EXACT_SOLID_DISTANCE = "exact_solid_distance"    # BRepExtrema solid-vs-solid
SOLID_CLASSIFIER_EVAL = "solid_classifier_eval"  # BRepClass3d_SolidClassifier
FACE_CLASSIFIER_EVAL = "face_classifier_eval"    # BRepClass_FaceClassifier
RAY_INTERSECTOR_PERFORM = "ray_intersector_perform"  # IntCurvesFace Perform
RAY_PAIR_CAST = "ray_pair_cast"                  # bidirectional ray pairs
SECTION_PAIR_ATTEMPT = "section_pair_attempt"    # per face-pair section
SECTION_ENGINE_CALL = "section_engine_call"      # intersect_models calls
RAW_INTERSECTOR_PERFORM = "raw_intersector_perform"  # IntTools_FaceFace
ADAPTIVE_EDGE_SAMPLES = "adaptive_edge_samples"  # verified section samples
COMPLETENESS_LEAF_INTERVALS = "completeness_leaf_intervals"
VOLUME_INTEGRATION = "volume_integration"        # _shape_volume calls
FACE_BOX_BUILD = "face_box_build"                # _conservative_boxes (faces)
EDGE_BOX_BUILD = "edge_box_build"                # _conservative_boxes (edges)
CURVE_ON_SURFACE_PROJECTION = "curve_on_surface_projection"
SINGLE_WITNESS_ATTEMPT = "single_witness_attempt"    # C9 shortcut tries
SINGLE_WITNESS_HIT = "single_witness_hit"            # shortcut applied
SINGLE_WITNESS_FALLBACK = "single_witness_fallback"  # full rule ran instead
REGION_PROPAGATED = "region_propagated"          # G12b propagated decisions
PREPARED_FACE_BOX_HIT = "prepared_face_box_hit"  # S1: reused prepared face boxes
PREPARED_EDGE_INDEX_HIT = "prepared_edge_index_hit"  # S1: reused prepared edges
BVH_QUERY = "bvh_query"                      # S2: HybridBoxIndex query via BVH
VECTOR_SCAN_QUERY = "vector_scan_query"      # S2: HybridBoxIndex query via flat scan

CANONICAL_COUNTERS = (
    EXACT_FACE_DISTANCE,
    EXACT_EDGE_DISTANCE,
    EXACT_SOLID_DISTANCE,
    SOLID_CLASSIFIER_EVAL,
    FACE_CLASSIFIER_EVAL,
    RAY_INTERSECTOR_PERFORM,
    RAY_PAIR_CAST,
    SECTION_PAIR_ATTEMPT,
    SECTION_ENGINE_CALL,
    RAW_INTERSECTOR_PERFORM,
    ADAPTIVE_EDGE_SAMPLES,
    COMPLETENESS_LEAF_INTERVALS,
    VOLUME_INTEGRATION,
    FACE_BOX_BUILD,
    EDGE_BOX_BUILD,
    CURVE_ON_SURFACE_PROJECTION,
    SINGLE_WITNESS_ATTEMPT,
    SINGLE_WITNESS_HIT,
    SINGLE_WITNESS_FALLBACK,
    REGION_PROPAGATED,
    PREPARED_FACE_BOX_HIT,
    PREPARED_EDGE_INDEX_HIT,
    BVH_QUERY,
    VECTOR_SCAN_QUERY,
)

ENV_VAR = "BREPKERNEL_PERF"


class PerfCounters:
    """Mutable per-call counters. Only touched while active."""

    __slots__ = ("counts",)

    def __init__(self) -> None:
        self.counts: Dict[str, int] = {}

    def bump(self, name: str, n: int = 1) -> None:
        d = self.counts
        d[name] = d.get(name, 0) + int(n)

    def snapshot(self) -> Dict[str, int]:
        # Canonical order, then any extra keys (forward compatibility).
        out = {k: self.counts.get(k, 0) for k in CANONICAL_COUNTERS}
        for k in sorted(self.counts):
            if k not in out:
                out[k] = self.counts[k]
        return out


# The active counters for this thread's current Boolean call, or None.
_active: contextvars.ContextVar[Optional[PerfCounters]] = (
    contextvars.ContextVar("brepkernel_perf", default=None)
)


def count(name: str, n: int = 1) -> None:
    """Bump a counter if profiling is active; otherwise do nothing."""
    p = _active.get()
    if p is not None:
        p.bump(name, n)


def current() -> Optional[PerfCounters]:
    """Return the active counters, or None when profiling is off."""
    return _active.get()


class scoped:
    """Install counters for one Boolean call (reentrant)."""

    def __init__(self, counters: PerfCounters) -> None:
        self._counters = counters
        self._token = None

    def __enter__(self) -> PerfCounters:
        self._token = _active.set(self._counters)
        return self._counters

    def __exit__(self, *exc) -> None:
        _active.reset(self._token)
        self._token = None


def wants_perf() -> bool:
    """True when the BREPKERNEL_PERF=1 environment override is set."""
    return os.environ.get(ENV_VAR, "") == "1"


def report_section(counters: PerfCounters) -> dict:
    """Build the report["performance"] payload for a finished call."""
    return {
        "enabled": True,
        "counters": counters.snapshot(),
    }
