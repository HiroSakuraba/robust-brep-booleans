# v0.8 freeform / NURBS Boolean path

This branch adds a conservative Tier B/C B-rep path without changing the
certified Tier A Boolean path.

## Design rule

The imported B-rep is the source of truth. Meshes and local NURBS acceleration
structures may reduce work, but they do not get to invent topology.

For difficult geometry the rule is:

> reject cheaply when we can prove two regions cannot interact; otherwise use
> the original trimmed surfaces, preserve provenance, verify the result, and
> refuse unresolved topology rather than guessing.

## Pipeline

```
STEP / OCCT B-rep
       |
       v
Solid -> Shell -> Face provenance
       |
       v
geometry-aware face AABB broad phase
       |
       v
NURBS knot-span control-hull broad phase
       |
       v
trimmed face / face section
       |
       v
verified 3D curve + p-curve on A + p-curve on B
       |
       v
local affected-face partition
       |
       v
exact patch classification against original solids
       |
       v
selected-face sewing
       |
       v
closed-shell nesting / cavity reconstruction
       |
       v
valid oriented OCCT solid(s)
```

## `src/brepkernel/freeform.py`

The freeform module provides a NumPy rational B-spline evaluator with first
derivatives and a conservative local acceleration index.

### Local knot-span acceleration

For positive-weight NURBS, each non-zero knot-span pair is associated with the
active control block. The convex-hull property gives a conservative local 3D
axis-aligned bounding box (AABB). Sweep-and-prune then removes local patch
pairs that cannot interact.

This is useful for:

- rejecting a face pair even when its coarse face AABBs overlap;
- seeding local point-to-surface projection;
- ranking local refinement work;
- avoiding uniform fine tessellation of an entire freeform face.

A partially overlapping synthetic regression reduces 256 naive span pairs to
160 conservative candidate pairs.

### Periodic NURBS

Periodic B-splines are no longer forced into a single whole-surface box.

The implementation:

1. copies the OCCT B-spline so the authoritative imported surface is never
   mutated;
2. records its exact U/V periods;
3. applies `SetUNotPeriodic()` / `SetVNotPeriodic()` to the copy;
4. translates the copied knot vector only by whole periods when necessary so
   its parameter domain covers the original face's UV window;
5. refuses acceleration if the face spans more than one representable period;
6. constructs the local knot-span control-hull index from the clamped copy;
7. independently samples the original OCCT face and the NumPy evaluator in
   the same UV frame and rejects the accelerator if they disagree.

This last check is important: geometric equivalence is not enough. Trim
p-curves are expressed in the original surface's parameter coordinates.

The current end-to-end periodic regression converts a sphere to a B-spline
before ingestion and produces **8 local patch records per spherical face**
instead of the previous single fallback box.

### Projection and refinement

- Local Gauss-Newton projection is seeded from nearby patch AABBs.
- A projected point is not accepted merely because it lies on the underlying
  surface; its UV must also lie in/on the actual trimmed face.
- Difficult trim cases can fall back to exact OCCT face distance.
- Curvature-ranked refinement is a scheduling heuristic only. It may request
  more work in difficult regions; it is not used as a correctness certificate.

**Current call-path boundary:** `project_point()`, `verify_points()`, and
`refinement_plan()` are tested freeform utilities, but they are not currently
called by `boolean_brep()` when deciding whether to accept a Boolean result.
The live acceptance path verifies section geometry against the original OCCT
trimmed faces and bilateral p-curves. The NumPy NURBS implementation currently
affects the live pipeline through conservative span AABBs / candidate culling
and accelerator metadata only. These utilities remain available for future
profiling-driven refinement work; their existence should not be read as part
of the present certification argument.

## `src/brepkernel/step_ingest.py`

STEP/OCCT ingestion preserves:

- solid IDs;
- shell IDs;
- face IDs;
- face orientation;
- exact OCCT UV bounds;
- surface type;
- geometry-aware face AABB;
- optional `FreeformFaceAccel`.

The candidate search is two-tier:

1. face-level geometry AABB sweep-and-prune;
2. NURBS local control-hull span rejection when both faces support it.

Unsupported accelerators cost performance, not correctness. The face remains a
candidate for the exact OCCT path.

## `src/brepkernel/intersection.py`

Only surviving face pairs receive expensive trimmed face/face section work.

For each section edge the kernel requires:

- a 3D curve;
- a p-curve on originating face A;
- a p-curve on originating face B;
- SameParameter correspondence;
- sampled agreement of all three geometric representations;
- UV membership in/on both actual trimmed faces.

The code also records transversality and seam risk.

Results distinguish:

- `curve`;
- `curve_near_tangent`;
- `point_contact`;
- `ambiguous_contact`;
- `disjoint`;
- `distance_unknown`.

A no-edge result is not automatically disjoint. Exact tangencies survive as
point contacts, and a no-curve pair inside the configured contact band becomes
ambiguous.

## `src/brepkernel/split.py`

Verified transverse section edges can partition only the faces they actually
touch.

- Unaffected faces pass through unchanged.
- All verified tools affecting one parent face are applied in one local
  `BRepAlgoAPI_Splitter` call.
- Near-tangent/seam-risk and unresolved contacts are not forced through by
  default.
- Every result is checked with OCCT validity.
- Child faces retain their parent-face ID.
- A child UV witness is checked against the original parent trim/support.

### High-accuracy area conservation

The split verifier uses the adaptive `BRepGProp.SurfaceProperties` overload,
not the loose default integration path. This mattered on rational B-spline
surfaces: the default property calculation was sufficiently inaccurate to
look like a split-area defect even when the split was correct.

The child areas must reproduce the parent area within the scale-aware split
tolerance.

## `src/brepkernel/assembly.py`

This branch now includes the first global B-rep material assembly stage.

### Exact patch classification

Each split patch gets a deterministic point in the actual trimmed face and is
classified against the other operand's original OCCT solids with
`BRepClass3d_SolidClassifier`.

Current material rules are:

| Operation | A patch | B patch |
| --- | --- | --- |
| Union | keep if outside B | keep if outside A |
| Intersection | keep if inside B | keep if inside A |
| A - B | keep if outside B | keep if inside A, reversed |

A witness classified ON/boundary or unknown is not majority-voted into a
result. Assembly refuses it.

Likewise, unresolved point/near-tangent contacts from the intersection stage
block assembly.

### Sewing and result construction

Selected exact B-rep faces are sewn with OCCT.

The result must have:

- zero free boundary edges;
- zero non-manifold multiple edges;
- a valid OCCT B-rep.

A full periodic face may itself be a closed skin (for example a one-face
sphere). If sewing returns such a standalone face rather than an explicit
shell, it is promoted to a one-face shell only after OCCT confirms that shell
is closed.

### Disconnected solids and cavities

Closed shells are normalized and geometrically nested.

A near-boundary interior witness is used for nesting rather than merely a
center-of-mass point. This avoids the concentric-shell error where the center
of an outer sphere is also inside its inner cavity shell.

Even-depth shells become material outers. Direct odd-depth children become
cavity shells. Disconnected outer shells become separate solids in a compound.

Final solids must be orientable, B-rep valid, and have positive volume.

### High-accuracy volume verification

B-spline result volume is measured with adaptive Gauss-Kronrod integration
(`VolumePropertiesGK`) with span-aware integration enabled. The default OCCT
volume-property path produced a visible integration error on the converted
NURBS sphere even though our assembled result and an independent OCCT Boolean
agreed exactly.

## Same-domain equivalence

`src/brepkernel/same_domain.py` adds a conservative fast path for
independently constructed closed B-reps that represent the same material
boundary.

OCCT's `AreFacesSameDomain` is useful evidence, but its implementation is
essentially an interior-point validity test and is not used alone. The kernel
requires:

- the same solid/shell/face counts;
- valid closed B-reps on both operands;
- model bounding boxes agreeing within a scale-aware tolerance;
- adaptive signed-volume agreement, including global material orientation;
- per-face bounding-box agreement;
- matching edge/wire topology signatures;
- adaptive face-area agreement;
- boundary-perimeter agreement;
- `AreFacesSameDomain(A,B)` **and** `AreFacesSameDomain(B,A)`;
- a deterministic one-to-one perfect matching of all faces.

This intentionally prefers false negatives. Equivalent solids with a different
face decomposition are not yet recognized by this fast path.

The exact identity shortcut was also tightened from `IsSame()` to
`IsEqual()`. OCCT documents that `IsSame()` ignores orientation, so a
reversed view of the same TShape must not be treated as identical material.

Current regressions pin:

- independently constructed identical boxes -> same-domain fast path;
- independently constructed identical NURBS spheres -> same-domain fast path;
- translated but otherwise identical boxes -> rejected;
- a reversed view of the same TShape -> rejected by global material
  orientation.

## Result edge / p-curve lineage

Verified section records now survive local splitting, and final assembled edges
carry lineage records.

For every unique result edge, the assembler records:

- adjacent selected patch references `(operand, parent_face, piece)`;
- parent input-face references;
- participating operands;
- whether the edge is still an original source boundary;
- matching verified Boolean-section references
  `(face_A, face_B, section_edge_index)`;
- whether bilateral verified p-curve evidence exists for that edge.

Splitter and sewing are allowed to copy or shorten section edges. Lineage first
uses topological identity, then a conservative geometric sub-edge check using
edge length and multiple point-to-section distance probes. This lets a final
edge remain traceable even if its TShape changed during local topology work.

The public `boolean_brep()` report summarizes result-edge lineage, including
Boolean-section edges, source-boundary edges and unattributed edges. Current
overlapping-sphere regressions require at least one final section edge with
bilateral operand provenance and verified p-curves; disjoint unions are checked
to ensure they do not invent section lineage.

## CI / regression status

The dedicated Freeform NURBS CI currently runs:

1. NURBS evaluator / broad-phase regressions;
2. trimmed surface-intersection regressions;
3. local face-split regressions;
4. global B-rep assembly regressions;
5. true NURBS end-to-end Boolean regression;
6. strict same-domain equivalence regressions;
7. public one-call B-rep pipeline regressions;
8. imported multi-face NURBS STEP regression;
9. adversarial periodic/sliver/oblique/imported NURBS corpus;
10. existing multi-shell semantic regressions.

The latest strict code-bearing run passes all ten groups with pinned
`cadquery-ocp 8.0.1.0.0`.

### Numerical checks

- NURBS value agreement vs OCCT: about `1.9e-15` max error.
- First derivatives vs OCCT: about `8.3e-15` max error.
- Known `0.003` normal offsets project back to `0.003`.
- Synthetic partial overlap: 160 conservative span pairs vs 256 naive.
- Far freeform faces can produce zero expensive section calls.
- Periodic converted sphere: 1/1 accelerators built on each operand,
  8 local patch records per face.
- Default converted-sphere section verification: 614 samples total across two
  section edges (max 377), down from 4,098 before local-tolerance sampling.
- Imported 26-face rounded STEP case: 392 verification samples total across
  eight section edges (33–65 each).

### End-to-end analytic B-rep assembly

For two radius-1 spheres with centers separated by 1, the assembled exact
B-rep path produces valid closed solids for union, intersection, and A-B and
agrees with independent OCCT Boolean results. Before the NURBS conversion
test was added, the analytic-surface regression errors relative to the
closed-form volumes were about `6e-9`.

The assembly suite also pins:

- disjoint union -> two material solids;
- disjoint intersection -> empty result;
- contained subtraction -> one material solid with a cavity shell;
- exact tangency -> refusal rather than speculative topology.

### True NURBS end-to-end regression

Both spheres are first converted to B-spline surfaces with
`BRepBuilderAPI_NurbsConvert` before this kernel sees them.

The strict run reports:

```
A surface: GeomAbs_BSplineSurface
B surface: GeomAbs_BSplineSurface
accelerators: 1/1 and 1/1
local patches: 8 and 8
verified section edges: 2
ambiguous contacts: 0
selected result faces: 2
result shells: 1
result solids: 1
```

The assembled union is B-rep valid, preserves selected-face provenance from
both operands, and matches the independent OCCT NURBS Boolean volume.


## Public Tier B/C API

The new path is exposed without changing the existing mesh/proxy
`boolean()` contract:

```python
from brepkernel import boolean_brep, BRepAmbiguousResult

try:
    result_shape, report = boolean_brep(shape_a, shape_b, "union")
except BRepAmbiguousResult as exc:
    report = exc.report
    # inspect report["refusal"] instead of accepting guessed topology
```

`boolean_brep()` accepts OCCT shapes or already indexed `BRepModel` objects
and runs ingestion, conservative candidate selection, verified intersection,
local splitting, exact patch classification, global assembly and final B-rep
validity checks.

Its report includes:

- solid/shell/face counts and accelerator/patch counts;
- candidate face pairs and expensive Section-call count;
- verified section-edge/contact counts;
- affected-face/split-call counts;
- unresolved contacts;
- per-patch material decisions and input provenance;
- selected faces, shells, solids, free/multiple edges and volume;
- complete result-edge lineage status;
- operation-level volume bounds and measured input/result volumes;
- final validity/closed/manifold status;
- per-stage wall-clock diagnostics.

Accepted ordinary results now require every final edge to have auditable
lineage. Boolean-section edges must retain verified bilateral p-curve
references. The final volume must also satisfy the basic mathematical bounds
implied by the operation (union, intersection, or A-B), using a
scale-aware tolerance. These checks are deliberately engine-independent and
can catch gross patch-selection or shell-orientation failures.

Passing `include_full_evidence=True` adds JSON-ready complete section
payloads to the report: curve parameters, XYZ samples, UV samples on both
input faces, tolerances, errors, transversality/risk data, SameParameter
repair state, and final result-edge links. The default report remains compact.

An exact oriented `TopoDS_Shape.IsEqual()` identity fast path handles
literal A op A without intersection. Independently constructed but equivalent
closed B-reps can also bypass intersection when the strict same-domain
recognizer proves a complete one-to-one boundary match. Union/intersection
return A and A-A returns an empty compound.

Refusals are surfaced as `BRepAmbiguousResult` with the partial stage report
and underlying typed cause attached.


## Different-decomposition canonicalization

Strict one-to-one face matching remains the preferred same-domain path, but
equivalent material can legitimately be represented with extra coplanar seams.

When strict matching fails, the optional second stage independently runs
`ShapeUpgrade_UnifySameDomain` on **deep geometry copies** of both operands.
The regression suite verifies that the original input keeps its pre-
canonicalization face decomposition. Before doing
that expensive normalization, invariant failures that canonicalization is
required to preserve now terminate the equivalence attempt immediately:
closed-solid availability, B-rep validity, model bounding box, adaptive volume,
and global material orientation. Topology-count or face-matching failures are
*not* terminal because canonicalization exists specifically to remove redundant
same-domain decomposition.

`ShapeUpgrade_UnifySameDomain` is then run only when it can plausibly help. A canonical
shape is accepted for comparison only if it preserves:

- B-rep validity;
- solid count;
- scale-aware bounding box;
- adaptive signed volume;
- global material orientation.

The canonicalized shapes then have to pass the **same strict face matcher**;
canonicalization is not itself an equivalence verdict.

The pinned regression compares:

- A: one ordinary `2 x 1 x 1` box;
- B: the same material produced by fusing two adjacent `1 x 1 x 1` boxes,
  retaining additional coplanar face decomposition before normalization.

The boundaries differ before canonicalization and resolve to the same six-face
boundary afterward. The public report records face/shell/solid counts before
and after normalization plus bbox and volume-preservation errors.

A separate adversarial regression uses two two-solid arrangements with the
same overall bounding box and the same total volume but material in different
corners. It remains non-equivalent, demonstrating that bbox/volume agreement
does not certify the fast path.

## Persistent verified section payloads

`BooleanAssemblyResult` now retains the complete verified section evidence,
not just section references.

For every verified Boolean section it copies and stores:

- curve parameters;
- sampled XYZ points;
- UV samples on operand A;
- UV samples on operand B;
- edge and verification tolerances;
- maximum surface-A, surface-B and cross-surface errors;
- minimum/maximum transversality;
- risk flags;
- whether SameParameter repair was needed;
- final result-edge indices derived from that section.

Final edge-lineage records reference these payloads by
`(face_A, face_B, section_edge_index)`. The ordinary public report stays
compact: it exposes sample counts, errors, tolerances and result-edge links
without dumping every XYZ/UV sample.

## Imported multi-face NURBS STEP corpus

The CI now contains a real STEP round-trip rather than only in-memory test
solids:

1. build a rectangular solid;
2. fillet all 12 edges;
3. convert the rounded B-rep to NURBS;
4. write it through `STEPControl_Writer`;
5. read it back through `STEPControl_Reader`;
6. run an A-B clip with a transverse box;
7. compare the assembled result against an independent OCCT cut oracle.

Latest measured run:

```
imported solids              1
imported faces              26
NURBS/freeform faces        26 / 26

naive face pairs           156
candidate face pairs         8
face pairs culled          148
OCCT Section calls           8
verified section edges       8
ambiguous contacts           0

local split calls             9
affected faces A              8
affected faces B              1
unresolved contacts           0

selected result faces        18
result shells                 1
result solids                 1
free edges                    0
multiple/non-manifold edges   0

final result edges           40
Boolean-section edges         8
source-boundary edges        24
unattributed edges            0
persisted section payloads    8

assembled volume      1.47334470035
OCCT oracle volume    1.47334470035
absolute error       4.44e-16
```

This case is useful because the imported shape is genuinely multi-face and
trimmed: all 26 imported faces are B-spline surfaces, eight separate source
faces intersect the cutter, and final edge lineage remains complete.

The default sampler has now been retuned so an unspecified public
`chord_tol` is derived from the section's accepted local OCCT tolerance
instead of the much smaller global base tolerance. Explicit caller
`chord_tol` values are still honored exactly, and every section retains 33
mandatory global verification parameters plus OCCT deflection proposals and
recursive local chord checks.

On the current imported rounded-box STEP regression the eight verified
sections use 33 or 65 samples each (392 total, maximum 65), replacing the older
513-sample curved-section pattern. On the converted-sphere union the two
sections dropped from 4,098 total samples to 614 total (maximum 377) while the
analytic/oracle, p-curve, lineage and validity checks remained green.

This is a verified reduction in work, not a claimed wall-clock speedup:
shared-runner timings remain noisy and OCCT Section / sewing still dominate
many runs.

## Exact bilateral curve-on-surface validation

The normal Tier B/C acceptance path now uses OCCT's
`BRepLib_ValidateEdge` in its exact mode after SameParameter repair.

For every accepted section edge, the kernel independently validates the 3D
edge curve against its p-curve lifted through **each** authoritative input
surface. OCP documents this mode as computing the maximum curve-on-surface
distance rather than checking only a finite control-point set.

Acceptance therefore requires, on both faces:

- a p-curve exists;
- SameParameter is true;
- exact curve-on-surface validation completes;
- the exact maximum distance is no larger than the section's verification
  tolerance;
- the underlying edge/face tolerance remains under the configured
  `max_section_tol` ceiling.

The previous adaptive XYZ/UV sampling remains in place for trim membership,
cross-surface comparison, transversality and risk classification. It is now a
complement to the exact curve-on-surface distance check rather than the only
3D/p-curve fidelity test.

On the current converted-sphere union, representative exact distances are
approximately:

```
section 0
  exact error on A  8.66003e-6
  exact error on B  8.66008e-6
  verification tol  1.73202e-5

section 1
  exact error on A  1.63152e-6
  exact error on B  1.63196e-6
  verification tol  1.0e-5
```

This materially strengthens the section-curve check, but it is not a proof
that OCCT discovered every connected component of the mathematical
surface/surface intersection. Completeness of the intersection set remains a
separate correctness problem.

### Alternative non-approximated shadow mode

`boolean_brep(..., shadow_section_crosscheck=True)` is an optional stress
mode. It builds the same face pair again with `Section.Approximation(False)`,
verifies that representation independently, and compares the two curve sets
bidirectionally in 3D plus total length.

This is **not** part of normal acceptance. OCCT's own source shows that
`Approximation(False)` can still turn a walking intersection into a B-spline;
it is not an exact analytic oracle. In the tiny NURBS-cap experiment that
alternative representation differed from the accurate primary section by
about `4.65e-5`, even though the primary result matched the independent OCCT
Boolean oracle. That experiment is why the shadow mode remains diagnostic
instead of vetoing ordinary accepted results.

## Intersection component completeness probe

Exact bilateral curve-on-surface validation certifies the geometric fidelity
of each **returned** section edge. A different failure mode is possible: a
lower-level surface/surface intersector can find several disconnected
components, while later Boolean post-processing accidentally drops one.

For every candidate pair involving at least one freeform face, the Tier B/C
path now reruns OCCT's lower-level `IntTools_FaceFace` intersector using the
same approximation tolerance used by `BOPAlgo_PaveFiller` (`1e-7`) and
calls `PrepareLines3D(False)`, matching the preparation step used by the
Boolean data structure.

Each bounded raw curve that has several samples lying IN/ON both trimmed input
faces is treated as a trimmed intersection component. Those raw components
must be geometrically represented by the verified final `Section` edge set.
If not, the operation refuses with `SectionCompletenessMismatch`.

This specifically protects against losing a curve branch between the
face/face intersector and final section-edge construction.

### Two-loop regression

A radius-`3` / tube-radius-`1` NURBS torus intersected by its equatorial
plane has two disconnected circular components from **one face pair**.

Current measured result:

```
candidate face pairs          1
final verified section edges  2
raw IntTools curves           2
raw trimmed components        2
raw unmatched components      0
max raw -> final distance      2.27e-15
```

Both final loops independently pass exact bilateral curve-on-surface
validation.

The torus loops touch the torus parameter seam. The intersection stage is
therefore accepted, but the full Boolean currently refuses later with
`RiskySectionCurve` / `seam_on_a` rather than guessing how to split the
periodic face. This is intentional: intersection completeness and safe seam
splitting are separate problems.

A negative regression deliberately gives the completeness checker only one of
the two verified loops. It refuses with:

```
SectionCompletenessMismatch
raw=2, trimmed=2, unmatched=1
```

The missing outer/inner loop is approximately distance `2` from the surviving
loop, so this is a real disconnected-component omission rather than a
parameterization difference.

### What this does and does not certify

OCCT's own `BOPAlgo_PaveFiller` contains special re-processing of
"potentially problematic Face-Face intersections" explicitly to avoid missing
section edges. The new probe adds another check around that post-processing.

It does **not** prove that `IntTools_FaceFace` itself found every connected
component of the mathematical intersection. The certification boundary is now:

1. conservative face-pair broad phase cannot drop a true candidate;
2. the lower-level OCCT face/face curve set is checked against final Section
   edges so post-processing cannot silently drop an observed branch;
3. each final edge gets exact bilateral curve-on-surface validation;
4. trims, transversality, lineage, splitting, sewing and final volume/topology
   are independently checked;
5. completeness of the *lower-level numerical surface intersector itself*
   remains an open certification problem.

## Adversarial NURBS corpus

CI now includes a separate corpus intended to exercise geometry classes that
were previously listed only as future work.

### Periodic NURBS cylinder

A cylinder is converted to NURBS and cut transversely by an analytic box.
The periodic B-spline support is preserved on the authoritative input face.

Measured run:

```
periodic input faces          1
NURBS accelerators            3 / 3
candidate face pairs          3
verified section edges        4
ambiguous contacts            0
free / multiple edges         0 / 0
unattributed result edges     0
assembled volume       4.13107608215
OCCT oracle volume     4.13107608215
absolute error         8.88e-16
```

### Sliver / tiny-feature sweep

A NURBS sphere is cut so only a spherical cap remains. Feature heights
`1e-3`, `1e-4`, and `1e-5` are exercised at `base_tol=1e-7`.
All three currently certify and agree with the independent OCCT oracle to the
precision printed by the tests. The smallest accepted cap has volume about
`3.14158e-10`.

The test contract still allows the smallest regime to become a typed refusal
if future tolerance changes make certification inappropriate; an inaccurate
accepted solid is never allowed.

### Oblique imported blend cut

The 26-face all-NURBS filleted STEP solid is cut with a box rotated 11 degrees
about the part center. The run keeps eight verified section edges, complete
edge provenance, one closed material solid, and oracle-identical volume to the
printed precision.

### Pre-trimmed periodic STEP round-trip

A cylinder is first notched by an OCCT Boolean, then the already-trimmed result
is converted to NURBS, written to STEP, read back, and cut again on the
opposite side. This exercises periodic parameter-frame recovery after previous
topology work and data exchange rather than only a pristine primitive.

Measured run:

```
imported NURBS faces          8
periodic imported faces       1
candidate face pairs          3
verified section edges        4
ambiguous contacts            0
Boolean-section result edges  4
free / multiple edges         0 / 0
unattributed result edges     0
assembled volume       4.20499408213
OCCT oracle volume     4.20499408213
absolute error         8.88e-16
```

## Runtime profiling

The public `boolean_brep()` report now records wall-clock timings (milliseconds)
for completed Tier B/C stages:

- ingest;
- same-domain/canonicalization check;
- intersection;
- local splitting;
- assembly;
- final B-rep verification;
- total call time.

Assembly reporting also summarizes section-verification sampling with section
count, total samples, maximum samples on one section, and mean samples. This is
intended to identify real production bottlenecks before introducing more
aggressive patch subdivision or approximation machinery.

Timings are diagnostics, not acceptance criteria, and can vary across machines.

The profiling also exposed avoidable same-domain work. Canonicalization is now
skipped when a preserved invariant already differs (for example model bounding
box, volume, or global material orientation). On one shared-runner
converted-sphere union this reduced the reported same-domain stage from roughly
374 ms before the short-circuit to roughly 30 ms afterward. That is a measured
example, not a stable benchmark.

## Performance boundary

The work funnel is now:

```
all face pairs
   ↓ face AABBs
candidate face pairs
   ↓ NURBS control-hull span rejection
candidate local interactions
   ↓ one OCCT trimmed-face Section call per surviving face pair
verified section geometry
   ↓ only affected faces split
classified exact patches
   ↓ sew only selected result patches
valid B-rep solid(s)
```

The span index currently rejects whole face interactions and accelerates local
projection/refinement. It does **not** yet subdivide a surviving face pair into
many small OCCT Section calls. That is intentional: doing so creates seam,
duplicate-curve, and curve-stitching problems. It should only be introduced if
profiling shows whole-face Section is the dominant cost on real STEP corpora.

## What remains

This is now materially beyond a mesh-first prototype, but it is not yet a
general certified NURBS Boolean kernel.

The main remaining work is:

- broader different-decomposition equivalence beyond cases reducible by
  same-domain face/edge unification;
- larger and genuinely industrial STEP/NURBS corpora: G2 blends, complex
  lofts/sweeps, pathological trim loops, imported tolerance damage, and
  assemblies with hundreds or thousands of interacting faces;
- more degenerate freeform contact classes: coincident seams, overlapping
  same-surface trims, multiple curves meeting at one vertex, and features at
  or below the configured tolerance scale;
- stronger guarantees on the **lower-level numerical intersector itself**:
  final Section post-processing is now checked against raw
  `IntTools_FaceFace` curve components and every returned edge gets exact
  bilateral curve-on-surface validation, but neither check proves that
  `IntTools_FaceFace` found every mathematical component;
- deciding whether/when the legacy mesh/proxy `boolean()` and the new
  `boolean_brep()` should share a common dispatch surface; they are currently
  separate public routes so the stable Tier A contract is not silently changed;
- a canonical `SolidComplex` result representation and durable on-disk
  evidence artifact format; full XYZ/UV evidence is now JSON-exportable from
  the public report but not yet versioned as a standalone interchange schema;
- certified local tessellation error bounds if tessellation is used for later
  acceleration or downstream consumers;
- profiling on production-scale STEP assemblies before adding more subdivision
  machinery.

The invariant remains unchanged: an optimization may remove work only when it
cannot remove a real geometric interaction. Ambiguous geometry is checked with
the exact B-rep, escalated, or refused.
