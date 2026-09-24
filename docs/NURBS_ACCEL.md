# v0.8 freeform / NURBS acceleration foundation

This branch adds the first Tier B/C geometry infrastructure without changing
the certified Tier A boolean path.

## Why this exists

The current prototype verifies analytic solids well, but imported STEP/NURBS
geometry creates two competing pressures:

1. **Accuracy:** trim loops, high curvature, seams, near tangencies, rational
   weights, and small local features make uniform coarse tessellation unsafe.
2. **Performance:** uniformly fine tessellation makes every downstream stage
   expensive even when 99% of two parts never come near one another.

The new path therefore separates **broad-phase rejection**, **local freeform
evaluation**, **trim validation**, and **final geometric verification**.

## New modules

### `src/brepkernel/freeform.py`

- NumPy rational B-spline evaluator with first derivatives.
- Expanded knot-vector extraction from OCCT `Geom_BSplineSurface`.
- Conservative knot-span AABBs from active control points for positive-weight,
  non-periodic NURBS.
- Sweep-and-prune knot-span pair culling.
- Fast local Gauss-Newton point projection using span-AABB seeds.
- Trim-aware acceptance. If the untrimmed closest point does not belong to the
  trimmed face, the fast path does **not** guess: it can fall back to OCCT
  `BRepExtrema_DistShapeShape`.
- Curvature-ranked refinement plan. This is explicitly a scheduling heuristic,
  not a correctness certificate: it can request more triangles in hard regions,
  never less verification.

### `src/brepkernel/step_ingest.py`

- Preserves STEP/OCCT `Solid -> Shell -> Face` hierarchy.
- Records per-face solid/shell IDs, orientation, exact OCCT UV bounds, surface
  type and precise geometry-aware AABB.
- Builds a `FreeformFaceAccel` for B-spline faces when supported.
- Two-tier candidate search:
  1. face-level OCCT AABB sweep-and-prune;
  2. NURBS span-level conservative control-hull culling.

Unsupported/analytic surface pairs remain face-level candidates. The optimizer
is therefore conservative: missing an acceleration path costs time, not
correctness.

## Local validation performed before push

The new regression test was run against the installed OCP build and passed:

- rational NURBS value agreement vs OCCT: ~1.9e-15 max error;
- first-derivative agreement vs OCCT: ~8.3e-15 max error;
- all sampled points stayed inside their control-hull span AABBs;
- known 0.003 normal offsets projected back to 0.003 distance;
- a separated NURBS pair produced zero candidate span pairs;
- a partially overlapping pair reduced 256 naive span comparisons to 160;
- OCCT box indexing preserved 1 solid / 1 shell / 6 faces;
- the face-level broad phase rejected far B-spline faces and retained the
  partially overlapping pair.

## What this does **not** claim yet

This is not a finished NURBS boolean kernel.

Still needed:

- exact/verified trimmed surface-surface intersection curves;
- p-curve/trim-loop provenance through boolean splitting;
- certified local tessellation error bounds for general NURBS rather than only
  a curvature scheduling heuristic;
- periodic-surface span subdivision with certified wrapped control hulls;
- integration of STEP face provenance into `classify.py` / Stage 6;
- a canonical `SolidComplex` rather than mesh-first result storage;
- freeform corpus regression against OCCT at several local feature scales.

The important boundary is maintained: optimization may reduce work, but
ambiguous geometry must still be refined, checked with OCCT, or refused.
