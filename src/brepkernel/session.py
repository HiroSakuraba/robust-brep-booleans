"""S6: BooleanSession -- reuse one prepared base across repeated operations.

An interactive CAD user rarely computes one isolated Boolean and exits: a
base part is edited repeatedly while tools move, features are suppressed,
or dimensions change. Rebuilding the base's immutable indexes
(face/edge boxes, adjacency, analytic masks, spatial indexes) on every
call is pure setup overhead. A BooleanSession prepares the base once and
reuses it; every operation still gets its own QueryContext (S4), so no
per-operation mutable state is shared.

Usage::

    with BooleanSession(base_shape) as session:
        r1, rep1 = session.boolean(tool_a, "difference")
        r2, rep2 = session.boolean(tool_b, "difference")
        r3, rep3 = session.boolean(tool_c, "intersection")

The session owns the prepared *base*, not Boolean results. If the base
geometry changes, create a new session.

This module also hosts the S7 Phase 1 reference ``boolean_brep_many``:
plain sequential semantics (each tool applied to the previous result),
which later S7 phases may optimize. It is intentionally *not* session
based, because the base changes on every step.
"""

from __future__ import annotations


class BRepSessionClosed(Exception):
    """Raised when boolean() is called on a closed BooleanSession."""


class BooleanSession:
    """Owns one prepared base shape for repeated Boolean operations."""

    def __init__(self, base, *, base_tol: float = 1e-7, **pipeline_defaults):
        """Prepare *base* once; *pipeline_defaults* apply to every op.

        *base* may be a raw TopoDS_Shape, a BRepModel, or an already
        prepared PreparedBRep (all flow through ``ensure_prepared``).
        """
        from .prepared import ensure_prepared

        self._base_tol = float(base_tol)
        self._prepared_base = ensure_prepared(base, base_tol=self._base_tol)
        self._defaults = dict(pipeline_defaults)
        self._closed = False
        self._ops = 0

    @property
    def prepared_base(self):
        """The PreparedBRep shared by every operation in this session."""
        return self._prepared_base

    @property
    def base_tol(self) -> float:
        return self._base_tol

    @property
    def op_count(self) -> int:
        """Number of completed boolean() calls."""
        return self._ops

    @property
    def closed(self) -> bool:
        return self._closed

    def boolean(self, tool, op, **kw):
        """Run one Boolean of the session base against *tool*.

        *tool* may be raw or prepared. Every call gets a fresh
        QueryContext inside ``boolean_brep`` (S4); only the immutable
        prepared base is shared. Returns (shape, report), and raises the
        same typed refusals ``boolean_brep`` raises.
        """
        if self._closed:
            raise BRepSessionClosed(
                "BooleanSession is closed; create a new session for "
                "further operations.")
        if "base_tol" in kw and float(kw["base_tol"]) != self._base_tol:
            raise ValueError(
                f"session base was prepared with base_tol={self._base_tol}, "
                f"got base_tol={kw['base_tol']}; mixing tolerances across a "
                f"shared prepared base is not allowed.")
        from .pipeline import boolean_brep

        merged = dict(self._defaults)
        merged.update(kw)
        merged.setdefault("base_tol", self._base_tol)
        result = boolean_brep(self._prepared_base, tool, op, **merged)
        self._ops += 1
        return result

    def close(self):
        """Mark the session closed; the prepared base is simply dropped."""
        self._closed = True

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False


def boolean_brep_many(base, tools, op="difference", **kw):
    """S7 Phase 1 reference: sequential multi-tool Boolean.

    Applies each tool to the previous result with plain ``boolean_brep``
    calls (no cross-tool optimization). Returns
    ``(final_shape, {"steps": [reports...], "optimized": False})``.
    Later S7 phases may return ``optimized: True`` once equivalence with
    this reference is demonstrated.
    """
    from .pipeline import boolean_brep

    out = base
    reports = []
    for tool in tools:
        out, report = boolean_brep(out, tool, op, **kw)
        reports.append(report)
    return out, {"steps": reports, "optimized": False}
