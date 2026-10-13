# OpenUSD B-rep exchange and the SMLib second engine

This branch connects brepkernel to NVIDIA's
[usd-brep](https://github.com/NVIDIA-Omniverse/usd-brep) (Apache-2.0,
v1.0.0). It adds two things:

1. **OpenUSD B-rep files.** Solids can be written to and read from the
   omniSolid `BrepArray` schema in `.usda` files, the format Omniverse
   uses for exact B-reps. This is pure Python on OCCT and needs no
   extra package.
2. **A second, independent engine.** With `pip install usd-brep`, every
   accepted Boolean can be recomputed by SMLib, the NURBS kernel inside
   usd-brep, and compared with brepkernel's result. Until now every
   certification check in brepkernel ran on OCCT, so OCCT was checking
   OCCT. SMLib shares no code with OCCT.

Nothing is vendored. usd-brep is an optional dependency (`brepkernel[usd]`)
and no NVIDIA source is copied into this repository. Its wheels are built
for CPython 3.12.

## Quick use

```python
from brepkernel import boolean_brep
from brepkernel.usd import write_omnisolid_usda, read_omnisolid_usda

out, report = boolean_brep(a, b, "difference", independent_engine="smlib")
report["stages"]["independent_engine"]["status"]     # "agree"
write_omnisolid_usda(out, "part.usda")                # opens in Omniverse
back = read_omnisolid_usda("part.usda")               # back to OCCT
```

`independent_engine="smlib"` outcomes:

| status | meaning | effect |
|---|---|---|
| `agree` | same solid count, volume within tolerance, no point contradicts | result served |
| `disagree` | a clear contradiction | typed refusal, kind `IndependentEngineDisagreement` |
| `inconclusive` | SMLib gave no usable answer (conversion failed, SMLib raised, or SMLib returned a broken body) | result served; with `independent_engine_strict=True`, refusal kind `IndependentEngineInconclusive` |
| `unavailable` | usd-brep not installed | as inconclusive |

The outcome is stored at `report["stages"]["independent_engine"]` and
summarized in the evidence record under `certification.independent_engine`.
When the option is off, reports and evidence records are byte-identical to
a build without it.

Lower-level entry points live in `brepkernel.usd.smlib_backend`:
`to_smlib`, `from_smlib`, `smlib_boolean` (SMLib's own answer as an OCCT
shape, uncertified), `compare_results`, and `cross_check`.

## How the comparison decides

SMLib's result is exported, read back into OCCT with brepkernel's own
reader, and compared with brepkernel's result in one place:

* solid count;
* volume, both computed by OCCT's integrator. Tolerance: a boundary offset
  of `band` over the larger surface area, floored at 1e-7 of the volume,
  where `band = max(1e-6, 1e-5 x bounding-box diagonal)`;
* 512 deterministic Halton points over the joint bounding box, classified
  in or out of both solids. Points within `band` of either boundary are
  skipped, so only clear contradictions count.

SMLib's own `compute_volume` is deliberately not used for the verdict.
It loses accuracy on B-spline faces: on a lofted solid whose geometry
imports correctly (SMLib's own fine mesh of it gives 11.63366 against
OCCT's 11.63365) `compute_volume` reports 11.62550, 7e-4 low.

What the check catches: wrong patch selection, missing or extra lumps,
inverted cavities, anything that moves volume or flips clear interior
points. What it does not catch: deviations thinner than `band`.

## The converters

`brepkernel/usd/omnisolid_write.py` (OCCT to BrepArray) and
`omnisolid_read.py` (BrepArray to OCCT) are written from the schema in
`source/schema/omniSolid/resources/schema.usda` and checked against
SMLib's importer and exporter. Things that had to be right, each found by
a failing test:

* **Analytic surfaces stay analytic.** Planes, cylinders and cones
  (`BrepSurfaceConeAPI`, a cylinder is semi-angle 0, as SMLib writes
  them), spheres and tori are written in the schema's own
  parameterization. OCCT's cone `v` runs along the slant line and the
  schema's along the axis, so `v` is scaled by `cos(semiAngle)`. Each
  mapping is checked numerically against OCCT's surface at nine points
  before it is written.
* **Left-handed frames.** On a mirrored part, or OCCT fillet cylinders,
  the surface frame is left-handed and OCCT's normal points inward, while
  the schema's cylinder, cone, sphere and torus normals always point
  outward. Such faces are written with their sides swapped and their loops
  reversed. Planes take their axis from the true normal.
* **Circles as quarter turns.** Edges on circles and ellipses are written
  as rational B-splines with at most 90 degrees per span, as SMLib does.
  OCCT's default packs a full circle into three 120-degree spans, and
  SMLib's integration then lost 2e-6 relative volume on a plain cylinder.
* **Angular ranges in the primary period.** Angular face ranges are
  shifted by whole turns to start in [0, 2pi). SMLib rebuilds generating
  circles from them and rejected a mirrored torus whose `v` range came out
  as [-2pi, 0].
* **UV curves are built, not projected.** OCCT needs a UV curve for every
  edge on every face; the file has exact 3D edges. OCCT's projection
  (ShapeFix) treats each edge alone, so an edge starting on a periodic
  seam can land on the wrong copy of the period and cut back across the
  face. OCCT's same-parameter step then quietly widens tolerances to hide
  it (seen: 0.068 on a file whose tolerance is 2e-5), and the shape still
  passes `BRepCheck_Analyzer`. The reader therefore builds them itself:
  on planes the 3D B-spline's poles are mapped affinely (exact); on other
  surfaces the edge is sampled per knot span, inverted onto the surface,
  unwrapped, fitted per span (joined C0 at the edge's own knots, because
  a circle stored as rational quarter arcs is only C1 there) and checked
  in 3D, aiming for 1e-10 of the part size and never accepted above a
  quarter of the file tolerance.
* **Seams and parameter windows.** A seam (an edge used twice by one face)
  gets two curves one period apart; the face lies left of each loop in UV,
  so the use running +v along a u-seam is the high copy. The seam then fixes
  the face's window (exactly one period), every piece of a seam shares that
  window, and every other edge use is moved by whole periods into it. An
  edge lying exactly on the window edge (an edge along the seam line used
  once) takes the copy that connects to its neighbours. Where a loop passes
  through a pole, the next edge takes the copy nearest the previous one.
  Each of these rules came from a failing case: a 7/8 sphere whose
  boundary runs through the pole, full and cut tori with two seams split
  into pieces, a cylinder face with an edge on its seam line.
* **Separate lumps.** All solids of a shape go into one Brep with one
  `solidRegion` each, as SMLib lays out the lumps of a split result, so
  `usd_brep.usd.import_brep` returns them as one body.
  `one_brep_per_solid=True` packs one Brep per solid instead (SMLib then
  needs `import_breps`).
* Every rebuilt OCCT shape must pass `BRepCheck_Analyzer`, have a positive
  volume, and keep every edge and vertex tolerance within 10x the file
  tolerance; an edge ending further than that from its vertex is refused
  with the gap in the message. Otherwise the reader raises
  `OmniSolidReadError` (a rejected shape rides on `.shape` for diagnosis).
* `from_smlib` adds one more net: the rebuilt solid count must equal
  SMLib's and the rebuilt volume must be within 5e-3 of SMLib's own. That
  bound is loose on purpose; it stops a converter error from ever being
  reported as a disagreement (it becomes "inconclusive").
* Volumes of rebuilt shapes need a tight integrator setting: their UV
  curves are fitted B-splines, and OCCT's default adaptive integration at
  1e-9 was 4e-7 off on a plain cylinder (1e-12: 6e-11; Gauss-Kronrod:
  exact). The backend and tests use brepkernel's own `_shape_volume` at
  1e-12 for analytic shapes.

## Measured (2026-10-12, OCP 8.0.1, usd-brep 1.0.0)

Round trip OCCT to `.usda` to OCCT, relative volume change: below 3e-14
on box, hole, cylinder, sphere, torus, cone, filleted box (26 faces),
lofted B-spline solid, hollow box, far-off rotated part, mirrored fillet
and a three-solid compound; 8e-11 on a mirrored torus-and-cone union.

OCCT to SMLib (`import_brep`, no healing), SMLib volume at 1e-8 against
OCCT: 1e-10 to 1e-7 on all of the above except the loft (the
`compute_volume` issue above).

SMLib to OCCT: all SMLib exports tried rebuild as valid OCCT solids,
including Boolean results with spheres, tori and cylinders, a two-lump
split and a three-Brep file. Fixtures are in `tests/data/usd`.

Cross-check on brepkernel's accepted results:

* hand-picked everyday cases (holes, shared-face union, cylinder T,
  fillet minus cylinder, torus minus box, far rotated boxes, a split into
  two lumps, disjoint parts): all agree;
* `smlib_crosscheck_fuzz.py`, seeds 7 (30 trials), 20261012 (60) and
  42 (150), 240 trials in all: 227 accepted and agree, **0 disagree**,
  2 inconclusive, 11 refused by brepkernel. Both inconclusive cases are
  SMLib defects: a sphere-cylinder intersection returned as a zero-volume
  body with no solid region (no error raised), and an export whose edge
  ends 2.4e-3 from its vertex, 120 times the file's stated tolerance;
* a result deliberately swapped for the wrong operation: `disagree`,
  511 of 512 points contradict.

The fuzz earlier reported three disagreements; every one was a reader
bug (wrong-period UV curves hidden by tolerance growth, an edge along a
seam line), found because SMLib's own view of its result was right. They
are fixed and are now guarded against twice (tolerance-growth refusal in
the reader, volume and solid-count check against SMLib in from_smlib).

Cost: 0.01 s (disjoint) to about 4 s per check, typically under 1 s.

## Findings about brepkernel from the second engine

SMLib can answer cases brepkernel refuses. Its answers are never served as
certified results; they are leads. The first round of leads has been
worked through (docs/REVIEW_LEDGER.md, T1-T6):

* Box union a sphere centred on its corner refused `NonManifoldResult`.
  Cause: two distinct degenerated (pole) edges counted as a self-touching
  edge pair. Fixed (T4); accepted with the exact volume.
* `InsufficientPatchWitnesses` on thin slivers: the witness search missed
  pieces filling a few percent of their UV box. Fixed (T2).
* `SectionCompletenessMismatch`: the completeness probe paired raw
  p-curves with the wrong faces. Fixed (T6); a safety fix as well.
* Seed 42, 150 trials: refusal rate 4.0% -> 0.67%, suspected
  over-refusals 3.3% -> 0%. The one refusal left is justified: OCCT's own
  Boolean drops a torus/plane branch (0.4% volume error against SMLib)
  and brepkernel's probe sees it.
* Equal-radius crossing cylinders (Steinmetz solid) still refuse by
  design; OCCT and SMLib both answer 16/3 to 1e-9, so it stays in the
  queue for a tangency gate.

Run `smlib_crosscheck_fuzz.py --second-opinion --dump-dir DIR` to grow
the queue; `boolean_brep(advisory_on_refusal=True)` gives the same second
opinion for a single call.

## Limits

* `.usda` text only; binary `.usdc` is refused with a clear message.
  Export from usd-brep with a `.usda` file name.
* Wire edges, sheets, vertex-only loops and non-manifold edges are refused
  by both converters. Only closed solids cross.
* Edges are read from `BrepCurve3dNurbAPI` only (the only edge type SMLib
  writes). The line, circle and ellipse edge APIs are refused.
* Authored UV curves (`BrepCurveUvNurbAPI`) are not written and are
  ignored on read; UV curves are always recomputed from the exact 3D
  geometry. The schema marks them optional and SMLib rebuilds them.
* `thisRadialEntryType` follows a local rule (top entry iff the edgeuse is
  `opposite`). SMLib accepts it for two-face edges, which are the only
  edges written. Its own exports do not follow a purely local rule, so
  non-manifold edges would need the full radial ordering.
* SMLib occasionally returns a broken body instead of raising; the backend
  treats that as no answer. SMLib's `compute_volume` is not trusted for
  verdicts (see above).
* `boolean_brep_many` and `BooleanSession` do not take the option yet.

## Regenerating the fixtures

```python
import usd_brep as U
a = U.create_box((0, 0, 0), 2, 2, 2); b = U.create_sphere((2, 2, 2), 1.2)
U.usd.export_brep(U.boolean(a, b, operation=U.BooleanOp.DIFFERENCE), "sm_boxsph.usda")
a = U.create_box((0, 0, 0), 2, 2, 2); b = U.create_sphere((2, 2, 2), 1.2)
U.usd.export_brep(U.boolean(a, b, operation=U.BooleanOp.UNION), "sm_corner_union.usda")
c = U.create_cylinder((0, 0, 0), 1.0, 2.0); d = U.create_cylinder((0.5, 0, 1), 0.6, 2.0)
U.usd.export_brep(U.boolean(c, d, operation=U.BooleanOp.UNION), "sm_cylcyl.usda")
a = U.create_box((0, 0, 0), 4, 2, 2); b = U.create_box((1.5, -1, -1), 1, 4, 4)
U.usd.export_brep(U.boolean(a, b, operation=U.BooleanOp.DIFFERENCE), "sm_split.usda")
# multi_sm.usda: write the three-solid compound of tests/test_usd_omnisolid.py
# with one_brep_per_solid=True, load it with U.usd.import_breps, and write it
# back with U.usd.export_breps.
```
