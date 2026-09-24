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
evaluation**, **trim validation**, **surface/surface intersection**, and
**final geometric verification**.

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

### `src/brepkernel/intersection.py`

The next stage consumes only the surviving face pairs.

- Runs `BRepAlgoAPI_Section` on **trimmed faces**, not infinite underlying
  surfaces.
- Requests 3D section curves and p-curves on both originating faces.
- Requires SameParameter correspondence. If OCCT returns an edge without it,
  the code repairs a copy and refuses if correspondence still cannot be
  established.
- Samples section curves adaptively according to 3D chord deviation.
- At every verification sample, independently evaluates:
  - the section edge's 3D curve;
  - face A at p-curve A's UV;
  - face B at p-curve B's UV.
  All three positions must agree within a tolerance derived from the input face
  and section-edge tolerances.
- Every sampled UV must classify IN/ON its actual trimmed face.
- Records surface transversality. Near-tangent intersections and periodic seam
  crossings are marked as risk conditions instead of being silently treated as
  ordinary transverse cuts.
- Distinguishes:
  - `curve`
  - `curve_near_tangent`
  - `point_contact`
  - `ambiguous_contact`
  - `disjoint`
  - `distance_unknown`
- A no-edge result is **not** automatically called disjoint. Exact tangencies
  are preserved as point contacts; a no-curve pair whose exact OCCT face
  distance lies inside the contact band becomes `ambiguous_contact`.
- Model-level work is driven by the conservative face/span broad phase, so
  distant faces cost zero section calls.

This moves the freeform contract from "triangle cuts look plausible" toward a
true B-rep statement:

```
section edge
    ├── 3D curve
    ├── p-curve on original face A
    └── p-curve on original face B

all three representations agree within tolerance
```


### `src/brepkernel/split.py`

The verified section curves are now usable as local B-rep split tools without
running a global solid Boolean merely to partition faces.

- Only faces touched by verified transverse section edges are sent to
  `BRepAlgoAPI_Splitter`; unaffected faces pass through unchanged.
- All section edges affecting one parent face are applied in one splitter call.
- Near-tangent, seam-risk, point-contact, and near-contact cases remain
  unresolved by default rather than being forced into guessed topology.
- Every split result is checked with OCCT's shape validator.
- The summed child-face area must reproduce the parent-face area within a
  scale-aware tolerance.
- Child pieces retain their parent face ID. A deterministic UV witness is
  sought on each piece and checked against the parent trimmed face and support
  surface; very thin pieces that do not yield a stable witness are flagged for
  stronger downstream verification rather than silently trusted.

This gives the freeform path a concrete local topology pipeline:

```
candidate faces
    ↓
verified 3D section + bilateral p-curves
    ↓
local face partition
    ↓
parent-face provenance retained
```

The next boundary is global assembly: joining split patches across adjacent
faces into coherent wires/shells and classifying which patches belong in the
requested Boolean result.

## Local validation performed before push

The NURBS evaluator/accelerator regressions were exercised against the installed
OCP build:

- rational NURBS value agreement vs OCCT: ~1.9e-15 max error;
- first-derivative agreement vs OCCT: ~8.3e-15 max error;
- all sampled points stayed inside their control-hull span AABBs;
- known 0.003 normal offsets projected back to 0.003 distance;
- a separated NURBS pair produced zero candidate span pairs;
- a partially overlapping pair reduced 256 naive span comparisons to 160;
- OCCT box indexing preserved 1 solid / 1 shell / 6 faces;
- the face-level broad phase rejected far B-spline faces and retained the
  partially overlapping pair.

The new intersection regression set pins:

- a transverse cubic B-spline / plane cut with verified p-curves;
- a UV-trimmed B-spline whose section is clipped to the actual face;
- a far face pair producing zero section calls;
- exact sphere/plane tangency surviving as a point contact;
- a near-tangent positive gap inside the contact band becoming
  `ambiguous_contact`, never "disjoint".

The local split regressions additionally pin:

- a verified NURBS section partitioning the affected trimmed face;
- child area conservation and parent-face provenance;
- incomplete open contours not being promoted into invented splits;
- far faces remaining bit-identical passthroughs;
- tangent and near-tangent contacts blocking speculative splitting.

## Why this should be faster

The intended work funnel is now:

```
all face pairs
   ↓  geometry-aware face AABBs
possible face pairs
   ↓  NURBS control-hull knot-span AABBs
possible local patch pairs
   ↓  OCCT trimmed-face section only here
verified section curves
```

Uniform fine tessellation is no longer the first response to difficult
freeform geometry.

## Why this should be more accurate

The original B-rep remains the source of truth. A tessellation can be used as
a computational proxy, but the emerging Tier B/C acceptance path can verify
geometry against:

- original trimmed faces;
- original NURBS surface evaluations;
- original p-curves;
- section-edge 3D geometry;
- exact OCCT trimmed-face distance in fallback cases.

In particular, being close to the *underlying untrimmed surface* is not enough:
a point or edge must also belong to the actual trimmed face.

## What this does **not** claim yet

This is not a finished NURBS boolean kernel.

Still needed:

- joining locally split patches across adjacent faces into robust global
  wires/shells;
- p-curve provenance through result assembly and classification;
- certified local tessellation error bounds for general NURBS rather than only
  a curvature scheduling heuristic;
- periodic-surface span subdivision with certified wrapped control hulls;
- same-domain/coincident-face handling beyond the current contact escalation;
- integration of STEP face provenance into `classify.py` / Stage 6;
- a canonical `SolidComplex` rather than mesh-first result storage;
- freeform corpus regression against OCCT at several local feature scales.

The important boundary is maintained: optimization may reduce work, but
ambiguous geometry must still be refined, checked with OCCT, or refused.
