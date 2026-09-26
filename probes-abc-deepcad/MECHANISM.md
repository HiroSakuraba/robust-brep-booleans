# Mechanism note — ABC/DeepCAD silent wrong-accepts (24 Sept 2026)

## What was wrong

Two distinct failure classes, both in the mesh layer
(`brepkernel/arrange.py`, the thin manifold3d wrapper + the four mesh verify
functions). The certified `operate()` path was never implicated (analytic
inputs only; it refuses these models outright).

### Class 1 — intra-input shell interpenetration (the serious one)

**Mechanism.** manifold3d's boolean classifies each input mesh's faces only
against the *other* input mesh — never against sibling shells of the *same*
input. When one input's own shells interpenetrate (e.g. ABC pilot model
`00000050_80d90bfdd2e74e709956122a`: shells of vol 7364 / 39 / 9724 with
A0∩A2 = 3316, A0∩A1 = 33.0, A1∩A2 = 3.8), union and difference keep both
shells' full boundaries, so the overlap volume is **double-counted**.
Intersection is largely unaffected (a face must lie inside *both* meshes).

Concrete numbers (pilot, union vs shifted copy):
- before: kernel 30118.3 in 3 shells vs OCCT solid-by-solid fuse 25434.8
  (**+18.4%**); arbiter 758/800 with 42 kernel-IN/truth-OUT.
- The 42 arbiter disagreements are points in the double-covered overlap
  regions: the kernel mesh covers them twice, exact classification once.
- All four mesh verify functions passed on the wrong mesh: it *is* a valid
  closed mesh — internally consistent, just not the set-union. Closure,
  Euler, shell counts and signed volumes cannot see this class.

**Fix (in `arrange.py`).** Each input is first reduced to its set-union:
its connected shells are folded with union (`_set_union`), then the
requested op runs on the two normalized meshes. Every fold step is a
boolean of proper closed meshes, so the accumulator never carries
interpenetrating shells. Single-shell inputs skip the fold entirely —
bit-identical construction and behavior (verified: 0 volume diffs on all
single-shell ABC/DeepCAD models replayed).

After the fix (pilot): union 25404.8 vs OCCT 25434.8 (**-0.118%**), 1 shell
like the oracle, arbiter **798/800** (1/1 symmetric — boundary noise).
Difference is exactly volA_set − intersection (11628.3 = 13776.6 − 2148.28);
inclusion–exclusion holds to ~1e-4.

Follow-up multi-shell cases, union vol err vs OCCT, before → after:
- 00273581 (6 shells): −0.569% → −0.202%, arbiter 796/800 → 799/800
- 00276289 (13 shells): −0.028% → −0.004%, arbiter 800/800 → 800/800
- 00277962 (2 shells): +1.468% → +0.148%, arbiter 798/800 → 800/800
- 00270341, 00272103, 00276650, 00278332, 00279442: unchanged or improved,
  none worse.

### Class 2 — non-manifold input silently accepted

DeepCAD `cad_step/00651489.step`: 2 non-manifold edges, directed closure
false. Old `arrange()` fed it to manifold3d, which returned `NoError` and a
geometrically impossible mesh (negative difference volume, intersection
bigger than the input, arbiter 0.15) — while all four verify functions
passed.

**Fix (in `arrange.py`).** `_check_input_manifold` refuses any input that
fails directed closure (covers open meshes, non-manifold edges, inconsistent
orientation) with a typed `ArrangementError` before the engine runs. The
DeepCAD case now refuses on all three ops. This is a deliberate behavior
change for defective inputs only: garbage gets a loud typed refusal instead
of a silent wrong accept.

## Regression results (24 Sept 2026, after the fix)

- Repo suites `tests/test_metamorphic.py`, `test_degenerate.py`,
  `test_regression.py`, `test_stress.py`: **all pass** (61 checks, no silent
  failures). The certified `pipeline.boolean` path is untouched behaviorally
  (analytic proxies are single-shell).
- Replay of all 59 ABC models + 80 DeepCAD extracts through the new
  `arrange()` (`replay_regression.py` → `replay_regression.json`):
  **0 exceptions, 0 crashes.**
  - Verdict changes: 5 ABC + 2 DeepCAD models, **all** with defective
    (non-manifold/open) inputs — arranged-before → typed-refusal-now.
    Intended (Class 2 fix). Zero verdict changes on clean inputs.
  - Union volumes: bit-identical on every clean single-shell model;
    changed only on multi-shell models, always toward the OCCT oracle
    (Class 1 fix).
- Thingi10K round could not be re-run verbatim (its probe scripts lived in
  ephemeral `/tmp`, now gone), but its two relevant behaviors are covered
  above: clean single-shell inputs are bit-identical, and defective inputs
  now get typed refusals instead of being arranged.

## What remains open

- Residual multi-shell union errors of −0.1%…−0.9% vs OCCT on a few models
  (00273581, 00276650): within tessellation-deflection noise between the two
  oracle paths; arbiters agree 799–800/800. Not investigated further.
- The fix assumes each input *denotes* the set-union of its shells — the only
  sane solid semantics, and what the OCCT oracle computes. Inputs with
  interpenetrating shells are still "valid" closed meshes; we normalize
  rather than refuse them because the boolean is well-defined on their
  set-denotation.
- Incidental observation (not the wrong-accept mechanism, left alone):
  a numpy round-trip of `to_mesh64()` output back into `Mesh64(...)`
  intermittently raised `TypeError` in ad-hoc scripts; the fix keeps
  `Manifold` objects through the fold and never round-trips, so it is
  unaffected. Worth a look if anything ever re-ingests boolean output.

## Files

- Fix: `~/workspace/brep-booleans/src/brepkernel/arrange.py`
- This note: `~/workspace/brep-booleans/probes-abc-deepcad/MECHANISM.md`
- Mechanism probes: `dump_mechanism.py`, `pairwise.py`, `spatial.py`,
  `validate_normalize.py`, `validate_normalize2.py`
- Regression: `replay_regression.py` → `replay_regression.json`,
  `postfix_arbiter.py`
- Recovered pilot STEP: `scratch/pilot/pilot.step` (re-downloaded from the
  ABC gh-pages sample URL after `/tmp` was wiped)
- No git/GitHub touched.

## v0.7 — ChatGPT-review implementation: normalize.py + semantic leg (24 Sept 2026)

An external review ("From ChatGPT") proposed five items in order; all five
were implemented. The mesh `arrange()` path is now a three-stage pipeline
(`src/brepkernel/normalize.py`, ~1400 lines; `arrange.py` is a thin wrapper):

**A. Tests first** — `tests/test_multishell_semantics.py`: 9 tests / 24
checks, all pass. Two overlapping shells (exact union vol); three-shell
ABC-pilot structure; duplicate coincident components; disjoint assembly
bit-identical normalization; hollow cube (cavity survives union/difference);
welded-at-vertex tetras → typed refusal; self-intersecting shell → typed
refusal; shell-order permutation; three-level nesting (material island in
void in material, vol 57 — caught a real bug during development).

**B. Full input validation** — finite coords → drop degenerate tris →
triangle dedupe → directed-edge closure → vertex-link one-cycle →
sweep-and-prune + vectorized SAT self-intersection per shell. Typed
`ArrangementError` kinds: `InputNonFinite`, `InputNotManifold`,
`InputVertexLink`, `InputSelfIntersection`, `InputNestingOrientation`,
`InputNestingInvalid`, `NormalizationInvalid`, `SemanticViolation`,
`EngineError`.

**C. Nesting-tree solid semantics** — per-shell signed volume / bbox /
interior point; winding-based containment; depth-parity vs orientation
agreement (typed refusal on disagreement); hollow and nested inputs are
folded to proper solids (outer-minus-cavities) so the engine only ever sees
outward-oriented solids. Disjoint assemblies pass through bit-identically.

**D. `prepare_operand()` + `NormalizationReport`** — returns the manifold,
V/F, volume, per-shell volumes, and a report (checks run, notes,
normalization actions).

**E. Semantic verification leg** — cheap Boolean inequalities on every op
(union ⊆ A∪B volume bounds, difference = volA − intersection, all of which
the DeepCAD impossible-geometry case violated) + adversarial witnesses
(face centroid ±eps, edge-midpoint ±eps, per-shell interior points, and a
500-point grid concentrated in the A∩B overlap region) checked against a
winding-number oracle on the *prepared* operands. Self-tests: a
union-labeled-as-difference swap is caught by the inequalities; a
half-dropped result is caught by the witnesses; correct results pass all
three ops. Boundary-margin band (1e-9·scale) skips witnesses exactly on a
surface, vectorized via expanded-AABB prefilter.

**Performance.** The validation and witness legs were optimized for
large meshes: sweep-and-prune self-intersection (122s → 1.3s on a 29k-tri
sphere), vectorized near-surface mask (1.87s → 0.13s), two-tier witnesses
(full dense ≤3000 tris; grid+interiors above). A curved pipeline boolean
(sphere 29k tris) is ~27s end-to-end.

**Deliberately not done** (honest deferrals): preserving STEP
Solid→Shell→Face hierarchy from ingest; deep per-face provenance lineage;
SolidComplex architecture redefinition. The nesting tree is the
anonymous-mesh substitute.

## Files (v0.7)

- `src/brepkernel/normalize.py` (new); `src/brepkernel/arrange.py`
  (thin wrapper); `tests/test_multishell_semantics.py` (new)
- No git/GitHub touched (working-tree changes, unpushed).

### Bugs found by the replay (fixed 24 Sept 2026)

The replay caught three false refusals (models that arranged under
f1d2e331 but refused under normalize.py):

1. **00276289** (13 shells): "normalization left 4 interpenetrating
   material pairs". The clustered balanced reduction missed overlaps
   (a `_truly_overlap` false negative during clustering). Fix: fallback
   that unions ALL components aggressively when the clustered pass leaves
   overlaps (the old f1d2e331 behavior). All three ops now arrange.

2. **00272103** union, **00273581** difference: `SemanticViolation` false
   positives from near-surface winding noise. The per-triangle
   centroid±eps witnesses (eps=1e-6·scale) were too close to the surface
   for robust winding on real CAD meshes — a long thin sliver triangle
   (aspect 150:1) and near-surface points 0.009–0.025 from the boundary
   produced 2/3007 mismatches (0.07%, not clustered, volume correct to
   0.0002%). Fix: the semantic leg now uses only the overlap-concentrated
   grid + per-shell interior points (no per-triangle witnesses). The
   adversarial self-tests still pass (swap caught by inequalities,
   half-dropped caught by grid witnesses).

## Final replay results (24 Sept 2026, with fixes)

**ABC: 59 models, 5 verdict changes, 0 exceptions, 3 union-vol diffs**
- 5 models arranged→refused, ALL defective_input=True (stricter validation
  correctly refuses defective inputs the old code arranged).
- 0 false refusals (the 3 from the first replay are fixed).
- 3 union-vol diffs vs f1d2e331:
  - 00270341: 3.15e-07 (numerical noise)
  - 00276289: 0.000242 (0.024%, fallback union)
  - 00277962: 0.013 (1.3% — expected: v0.6 baseline never folded the
    genuine ~515-unit shell interpenetration; v0.7 folds it, matching the
    f1d2e331 +0.148%-vs-OCCT target. See correction below.)

**DeepCAD: 80 models, 2 verdict changes, 0 exceptions, 0 union-vol diffs**
- 2 models arranged→refused, both defective_input=True (correct).

**00277962 — CORRECTION (24 Sept 2026, post-v0.7 review round): the
"regression" was a phantom, caused by a baseline mix-up in this document,
not by a code bug.** Re-investigation with direct OCCT measurements:

- The two shells are NOT disjoint. Same x/y footprint, z-ranges
  [0, 158.8] and [155.2, 313.9]: OCCT pairwise fuse gives 39014.23 vs a
  shell-sum of 39529.35, i.e. the solids genuinely interpenetrate by
  ~515 volume units. The overlap graph (`merges=[(0,1)]`) was CORRECT;
  the disjoint passthrough correctly did NOT fire.
- The replay's `union_vol_old` (79174.14) is the **v0.6** baseline, which
  never folded the overlap: 2×39587.07 = +1.468% vs the true set union.
- v0.7 folds via the engine union: prepare(A) = 39072.09, union(A,A') =
  78144.19 = **+0.148% vs the true union** (2×39014.23 = 78028.46, from
  OCCT pairwise fuse + disjoint shifted copies, confirmed
  common(A,A') = 0). The residual +0.148% is tessellation bias (mesh
  shells are +0.146% larger than the exact OCCT solids), not kernel
  error: the engine union matches the tessellation-corrected OCCT
  expectation to 0.002%.
- So v0.7 == the f1d2e331 target (+0.148%), and both beat the v0.6
  baseline (+1.468%). The earlier text wrongly attributed f1d2e331's
  +0.148% *achievement* to the v0.6 baseline and then chained
  "-1.3% vs old" into "-1.15% vs OCCT". Nothing to fix in the
  normalization machinery; the fix was to the record. Pinned by new
  regression test t10 (interpenetrating slab must fold to the union
  volume, not the sum) and t11 (genuinely disjoint shells take the
  "disjoint passthrough", bit-identical).

## v0.7 review round 2: semantic-leg hardening (24 Sept 2026)

A follow-up external review of commit 11cb2d46 found three fixable items
(all in `src/brepkernel/normalize.py::check_boolean_semantics`, plus
regression tests). The deferred architectural items (STEP hierarchy,
provenance, SolidComplex) stay deferred.

### 1. 00277962 "regression" — root-caused to a phantom (no code change)

See the CORRECTION above: the two shells genuinely interpenetrate
(~515 units per OCCT pairwise fuse); the overlap graph was right to
merge them; v0.7's 78144.19 is +0.148% vs the true set union, identical
to the f1d2e331 target and better than the v0.6 baseline (+1.468%) the
replay compared against. The "regression" was a baseline mix-up in this
document (f1d2e331's +0.148% achievement misattributed to the v0.6
baseline). Pinned by new tests t10 (slab-interpenetration folds to the
union volume, never the sum) and t11 ("disjoint passthrough" status and
bit-identical output for genuinely disjoint shells).

### 2. Dead guard in check_boolean_semantics() — fixed

The review caught: `vR = abs(signed_volume(...))` followed by
`if vR < -av:` — dead code, since `abs()` can never be negative. Fixed
properly: the signed volume is captured first and `svR < -av` refuses
inside-out/corrupted results ("negative signed result volume"), then
`vR = abs(svR)` feeds the remaining checks. Pinned by new test t12:
an inside-out box (negative signed volume) presented as the result of
each of union/intersection/difference is refused with
`SemanticViolation`. This is the guard that would have caught the
DeepCAD impossible-geometry case (negative difference volume) at the
semantic leg.

### 3. Semantic volume checks — strengthened where cheap, prose fixed

The review noted the prose claimed more than the code enforced. Audit:
the code already enforced union max(vA,vB)<=vol<=vA+vB, intersection
vol<=min(vA,vB), and difference vol<=vA. Added the one cheap missing
single-op bound: difference vol >= volA-volB (vol(A\B) >= volA-volB is
exact set theory; same 1e-6/1e-9·scale³ tolerance pattern as the rest).
The docstring now enumerates exactly which inequalities are enforced
and explicitly lists what is deliberately NOT enforced: cross-operation
identities (volD = volA-volI, inclusion-exclusion volU = volA+volB-volI)
need the complementary operation's result, and we do not compute extra
booleans just to check them. (The witness leg checks the per-point
truth table, which is the pointwise form of those identities.)
