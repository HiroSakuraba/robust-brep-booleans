# robust-brep-booleans

A prototype of **robust B-rep Boolean operations** for CAD: topology is decided
only by exact tools, geometry is treated as bounded approximation, and the
pipeline never silently returns a broken solid.

This implements the Tier A (analytic solids) slice of the design in
`docs/` — "B-reps that don't break: a design for robust CAD Booleans".

## The idea in one paragraph

Boolean operations on boundary representations fail in practice because
floating-point geometry is used to make *topological* decisions
(is this face inside or outside?). This prototype separates the two:
every inside/outside verdict comes from an **exact** implicit function,
while the triangle mesh is only a **certified proxy** — a tessellation with
a proven bound on its chordal deviation from the true surface. Anything the
exact tools can't decide within their margins is reported explicitly as an
ambiguity instead of being guessed.

## Pipeline (7 stages)

| Stage | Module | What it does |
|---|---|---|
| 0 | `ingest.py` | Audits input parameters; `ToleranceLedger` records every epsilon in one place |
| 1 | `proxy.py` | `certified_proxy()`: mesh with a *proven* chordal-error bound (asserted, not hoped) |
| 2 | `arrange.py` | Proxy boolean (manifold3d). The exact CGAL/libigl core plugs in behind this interface |
| 3 | `classify.py` | **Tier A**: exact degeneracy detection from defining parameters (coincident planes, tangent cylinders) — no tolerance voting |
| 4 | `classify.py` | Classification with margins from **exact** implicits; ON-margin faces are never decided by the proxy |
| 5 | `assemble.py` | Per-face provenance + engine/classifier cross-check |
| 6 | `verify.py` | Independent verification: closure, Euler characteristic, orientation, Monte-Carlo field check against the exact implicits, manifold3d cross-check |

## The contract

```python
from brepkernel import boolean, AmbiguousResult
from brepkernel.solids import Box, Sphere

try:
    mesh, report = boolean(Box([0,0,0],[2,2,2]), Sphere([1,0,0], 1.5), "union")
except AmbiguousResult as e:
    # e.mesh: partial result. e.report: full stage-by-stage ambiguity report.
    ...
```

`boolean()` either returns a mesh that passed **all** Stage 6 checks, or raises
`AmbiguousResult` carrying the partial mesh and the full report. It never
silently returns a broken solid.

## Quickstart

```bash
python3 -m venv .venv
.venv/bin/pip install manifold3d
.venv/bin/python tests/test_metamorphic.py
.venv/bin/python tests/test_degenerate.py
```

(`cadquery-ocp` and `rhino3dm` are needed only for the future STEP/freeform slice.)

## Test results

- `tests/test_metamorphic.py` — **13/13 pass**: A∪A=A, A−A=∅, disjoint
  union/intersection identities, commutativity of ∪/∩, translation invariance,
  inclusion–exclusion, cylinder/sphere/cone smoke tests.
- `tests/test_degenerate.py` — **8/8 pass, zero silent failures**: coincident
  faces, tangent cylinders, point-touching spheres, corner-touching boxes,
  nested boxes. Degeneracies are detected analytically (Tier A) and either
  resolved or explicitly reported.

## Honest limits

- **Arrangement engine**: manifold3d (robust-inexact) stands in for the exact
  libigl/CGAL core — no `pyigl` wheel exists on PyPI (conda-only). Topology
  verdicts don't depend on it (they use the exact implicits), so this affects
  only intersection-curve geometry.
- **Freeform (NURBS) faces** — Tier B/C — are not in this slice. That's the
  next build, via OCCT STEP ingest + openNURBS as an independent evaluator.

See `docs/PROTOTYPE.md` for the full write-up.

## License

MIT — see `LICENSE`.
