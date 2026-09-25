# Review Gate Ledger

One entry per gate. A gate is closed only when every criterion is PASS;
anything else stays open and is recorded honestly here (invariant I6).

---

## G0 - Tooling, baseline, CI, dependency pins

- Date: 2026-09-25
- Branch: gate/G0-tooling (local only; never pushed, never merged to main)
- Base: 5554713 "Merge v0.8 exact NURBS Boolean pipeline"
- Commits:
  - 6814862 "tools: add review probe scripts and review_20260925 brep fixtures"
  - 45fb234 "ci: add Tier A workflow and pin freeform dependencies"
- Environment: ~/workspace/brep-booleans/.venv (CPython 3.12.3;
  the work plan brief said 3.11 - deviation D1),
  numpy==2.5.3, manifold3d==3.5.3, cadquery-ocp==8.0.1.0.0 (OCCT 8.0.1).

### Deviations from the plan

- D1: venv Python is 3.12.3, not 3.11 as the brief stated. All runs used
  ~/workspace/brep-booleans/.venv/bin/python.
- D2: TASK 1 was not done by HTML extraction. Ben uploaded the appendix
  files as a zip; per parent instruction the files were copied verbatim
  from /tmp/review_probes_zip/ (repo-relative layout already correct).
  A diff against the plan-appendix HTML extraction showed the zip is
  byte-identical in content modulo trailing newlines; the zip contained
  zero U+2014 em dashes.
- D3: a VM/runtime restart mid-baseline killed the first fuzz run
  (probes.txt and repro.txt were already complete). The three fuzz runs
  were re-run from scratch; fuzz is seed-deterministic so results are
  identical to a fresh run.
- D4: OCP.Standard exposes no Standard_Version symbol in this build; the
  OCCT version (8.0.1) is inferred from cadquery-ocp==8.0.1.0.0 and recorded
  in docs/baseline_20260925/versions.txt.
- D5: arbiter.py produces no report output of its own, so the Task 4
  versions print is implemented as arbiter.version_banner() and is printed
  at the top of the report output of fuzz_brep.py, common_cad_probes.py,
  and repro_findings.py (and stored in the fuzz JSON "versions" field).

### Commands run

```
python tools/review_probes/common_cad_probes.py --baseline | tee docs/baseline_20260925/probes.txt
python tools/review_probes/repro_findings.py | tee docs/baseline_20260925/repro.txt
python tools/review_probes/fuzz_brep.py --trials 150 --seed 7 --out docs/baseline_20260925/fuzz_analytic.json
python tools/review_probes/fuzz_brep.py --trials 60 --seed 11 --nurbs --out docs/baseline_20260925/fuzz_nurbs.json
python tools/review_probes/fuzz_brep.py --trials 200 --seed 5 --snap 0.5 --kinds box,cyl --points 150 --out docs/baseline_20260925/fuzz_snapped.json
for t in tests/test_*.py; do python $t; done   # all 14 test files
python tools/review_probes/fuzz_brep.py --trials 30 --seed 7 --out /tmp/tiera_fuzz.json   # CI step validation
```

### Baseline results (trimmed)

Probes (17 cases, --baseline mode): wrong=0, unmet_accepts=11, exit 0.
4 accepted and correct (slot overshoot cut, through-hole cylinder,
crossing cylinders r1 != r2, sphere-from-corner cut); 10 refused as
expect=accept (all assembly/UnresolvedContact - finding F1, the G2 work);
the singular equal-radius crossing cylinders correctly refuse as
expect=refuse; touching-edge union refuses as expect=either.

repro_findings (exit 1, expected at baseline): F2 REFUSED with
SectionCompletenessMismatch (raw=3, trimmed=2, max_distance=9.41089e-05,
tol=4e-05), OCCT common volume 0.54081549348893. F4: kernel union correct
(winding +0.0000 at the probe point, volume == vA + vB - vCommon) while
BRepClass3d_SolidClassifier reports IN for a point >1.0 from both
surfaces - on the kernel union and on OCCT's own fuse alike (finding F4).

Fuzz tallies:
- analytic (150 trials, seed 7): accept 143, refuse:InsufficientPatchWitnesses 5,
  refuse:PatchClassificationInconsistent 1, refuse:SectionToleranceTooLoose 1,
  WRONG 0, CRASH 0
- nurbs (60 trials, seed 11): accept 59, refuse:InsufficientPatchWitnesses 1,
  WRONG 0, CRASH 0
- snapped (200 trials, seed 5, snap 0.5, kinds box,cyl): accept 93,
  refuse:UnresolvedContact 91, refuse:RiskySectionCurve 13,
  refuse:SectionOutsideTrim 2, refuse:FaceSplitAreaMismatch 1,
  WRONG 0, CRASH 0

### CI validation

- .github/workflows/tierA-ci.yml parses as valid YAML (name "Tier A CI",
  one job "tiera", 9 steps, push-to-main + pull_request path filters).
- Each step was run locally: pip install step maps to requirements-freeform.txt
  (venv already holds the pinned versions); the four Tier A test files all
  exit 0; common_cad_probes.py --baseline exits 0 (fails only on wrong
  answers); fuzz --trials 30 --seed 7 exits 0 with TALLY
  {accept: 29, refuse:InsufficientPatchWitnesses: 1} - no WRONG, no CRASH.
- freeform-ci.yml (pre-existing): all ten of its test steps were run as part
  of the full 14-test sweep below, all exit 0.
- Remote (GitHub) CI green is PENDING a push; local-only branch, so remote
  status is unverified by construction. Recorded honestly, not claimed.

### Full suite (TASK 5)

14/14 test files exit 0 (venv python):
test_brep_pipeline (11s), test_degenerate (47s), test_freeform_assembly (53s),
test_freeform_intersection (5s), test_freeform_nurbs (2s),
test_freeform_split (1s), test_metamorphic (60s), test_multishell_semantics (1s),
test_nurbs_adversarial_corpus (14s), test_nurbs_boolean_end_to_end (7s),
test_regression (26s), test_same_domain (1s), test_step_nurbs_multiface (3s),
test_stress (136s).

### Gate verdicts

| Criterion | Result |
|---|---|
| Both CI workflows defined and locally validated | PASS (remote green pending push, recorded as unverified) |
| Baseline files committed with tallies | PASS |
| Fuzz baseline 0 WRONG 0 CRASH | PASS (all three runs) |
| Full suite green (14/14) | PASS |
| I1-I9 invariants held | PASS (no crash, no tolerance loosened, no push/merge) |

G0 gate: CLOSED. Remote CI verification is the only deferred item and is
blocked on a push, which the plan forbids at this stage.

---

## G1 - Independent arbiter and volume identities

- Date: 2026-09-25
- Branch: gate/G1-arbiter (local only; never pushed, never merged to main)
- Base: 3c494ab "docs: add G0 gate entry to the review ledger"
- Commits:
  - b5c9fdd "tests: add G1 arbiter helper and negative corruption test"
  - 6dfeeb9 "pipeline: add difference volume lower bound and crosscheck_ops mode"
  - 0af4caa "tests: add membership-audit assertions to Tier B/C acceptance tests"
- Environment: ~/workspace/brep-booleans/.venv (CPython 3.12.3),
  numpy==2.5.3, manifold3d==3.5.3, cadquery-ocp==8.0.1.0.0 (OCCT 8.0.1).

### What was built

Task 1 (arbiter assertions on Tier B/C tests):
- New tests/_arbiter.py: thin wrapper over the G0 membership_audit
  (tools/review_probes/arbiter.py). check_accepted() asserts
  kernel_errors == 0 and checked >= 200 via the file's check();
  raw_audit() runs the audit with no assertions (for the negative test);
  audit_bounds() builds the sampling box from the input bboxes.
  Pure memoization of triangle soups and winding numbers per test process
  (bit-identical results) keeps repeat audits cheap.
- 23 audits added across 5 files, on every test that ACCEPTS a Boolean
  result: test_brep_pipeline (8: p1, p2 union/difference, p4
  union/difference, p5, p6 intersection/difference), test_freeform_assembly
  (6: t1 union/intersection/difference, t2, t3, t4),
  test_nurbs_boolean_end_to_end (2: n1, n2), test_step_nurbs_multiface (1),
  test_nurbs_adversarial_corpus (6: a1, a2 at h=1e-3/1e-4/1e-5, a3, a4).
  Refusal-expecting tests and internal-stage files
  (test_freeform_nurbs/intersection/split, test_same_domain,
  test_multishell_semantics) accept no Boolean result and get no assertions.

Task 2 (volume identities and cross-operation checks):
- difference lower bound in brepkernel/pipeline.py: volume_bounds =
  [max(0, va - vb), va] with the paired check
  (vr <= va + tol) and (vr + tol >= va - vb). Union upper and intersection
  lower/upper unchanged. No tolerance values changed (invariant I3).
- New optional boolean_brep(..., crosscheck_ops=False): after acceptance the
  two companion ops run on the same inputs (recursion forced off) and the
  identities vol(AuB) + vol(AnB) = vol(A) + vol(B) and
  vol(A-B) = vol(A) - vol(AnB) are checked within max volume_tol. Violation
  raises a typed refusal kind=OperationIdentityFailed (stage crosscheck,
  accepted=False, stage report attached). A companion op that itself refuses
  is recorded as skipped, never silently passed. _safe_volume maps empty
  (no-solid) results to 0.0. All three return sites route through _finalize.
- crosscheck_ops=True on correct results (overlapping/disjoint spheres x 3
  ops, identity fast path) produces no false refusals; identities ok=True,
  empty-intersection companion handled.

Task 3 (negative test):
- New tests/test_g1_negative_corruption.py, committed FIRST as failing per
  invariant I4. t1: disjoint-sphere union passes all pre-G1 checks silently
  (documents the blind spot). t2: monkeypatching _decision_rule to flip the
  first keep flag produces a corruption the pipeline accepts-but-wrong
  (result = sphere B alone, vol 4.188790 vs true 8.377580, closed, valid,
  volume_bounds_ok=True); the arbiter catches it (kernel_errors=9,
  checked=300, audit_ms=11865). t3: boolean_brep(crosscheck_ops=True) on the
  corrupted union raises typed OperationIdentityFailed
  (union_plus_intersection err=4.18879 vs tol=0.00112373). Note: the
  overlapping-sphere corruption variant was tried first and refused via
  OpenAssembly, so the disjoint variant is the one that is silent.
- Full file: ALL PASS (3/3) on unpatched code, exit 0, ~43s.

### Deviations from the plan

- D1 (G1): /tmp/plan/PartA.txt did not exist; the G1 tasks, four-step test
  protocol, and pass criteria were taken from the uploaded work-plan file
  workspace/user/files/brepkernel_muse_workplan_0_p3h8.html (G1 section),
  which lists the same items. Nothing in the task list was dropped.
- D2 (G1): the audit is slow on NURBS geometry (one tessellation at 2e-4
  deflection = 115k triangles; per-audit cost 4-42s, winding-number
  sampling dominates). Mitigated with bit-identical memoization inside
  tests/_arbiter.py; the five modified files take 8s-134s total and stay
  deterministic. Timing is recorded below, not hidden.
- D3 (G1): test_nurbs_adversarial_corpus t2's 1e-5 sliver was ACCEPTED (not
  refused) and audited cleanly, so that file carries 6 audits, not 5.

### Commands run

```
python tests/test_g1_negative_corruption.py          # negative test (committed failing first per I4)
python tests/test_brep_pipeline.py                    # 117s, 27 checks
python tests/test_freeform_assembly.py                # 134s, 39 checks
python tests/test_nurbs_boolean_end_to_end.py         # 64s, 10 checks
python tests/test_step_nurbs_multiface.py             # 8s, 8 checks
python tests/test_nurbs_adversarial_corpus.py         # 73s, 24 checks
for f in test_degenerate test_freeform_intersection test_freeform_nurbs \
    test_freeform_split test_metamorphic test_multishell_semantics \
    test_regression test_same_domain test_stress; do python tests/$f.py; done
git grep -P '\x{2014}' -- tests/_arbiter.py tests/test_g1_negative_corruption.py \
    src/brepkernel/pipeline.py tests/test_brep_pipeline.py \
    tests/test_freeform_assembly.py tests/test_nurbs_boolean_end_to_end.py \
    tests/test_step_nurbs_multiface.py tests/test_nurbs_adversarial_corpus.py \
    docs/REVIEW_LEDGER.md  # 0 matches: no em dashes
```

### Results (trimmed)

- test_brep_pipeline: ALL PASS, 27 checks, ~117s; 8 audits, all
  kernel_errors=0, checked>=291 (p1 41.9s, p2 4.7/3.0s, p4 143/54ms,
  p5 154ms, p6 30.1/26.5s).
- test_freeform_assembly: ALL PASS, 39 checks, ~134s; 6 audits, all
  kernel_errors=0, checked>=299.
- test_nurbs_boolean_end_to_end: ALL PASS, 10 checks, ~64s; 2 audits,
  kernel_errors=0, checked=300/299.
- test_step_nurbs_multiface: ALL PASS, 8 checks, ~8s; 1 audit,
  kernel_errors=0, checked=297.
- test_nurbs_adversarial_corpus: ALL PASS, 24 checks, ~73s; 6 audits, all
  kernel_errors=0, checked>=298.
- Negative test: t1 documents the silent blind spot; t2 arbiter catches the
  flipped decision rule (kernel_errors=9, checked=300); t3 crosscheck raises
  typed OperationIdentityFailed.
- Full sweep of the remaining 9 files (degenerate, freeform_intersection,
  freeform_nurbs, freeform_split, metamorphic, multishell_semantics,
  regression, same_domain, stress): all exit 0, 0 [FAIL] lines.
- Total: 15/15 test files exit 0 (14 existing + tests/test_g1_negative_corruption).
- No em dashes in any touched file, prose, or commit message (verified with
  git grep -P '\x{2014}': 0 matches).

### Gate verdicts

| Criterion | Result |
|---|---|
| Tier B/C tests pass with arbiter assertions added | PASS (23 audits, 0 kernel errors, all checked>=200) |
| Negative test catches the corruption (fails as intended pre-fix, catches post-fix) | PASS (t2 kernel_errors=9; t3 typed OperationIdentityFailed) |
| Full suite green, ledger updated | PASS (15/15 exit 0; this entry) |
| I1-I9 invariants held | PASS (no crash, no tolerance changed, no push/merge, local-only branch) |

G1 gate: CLOSED.

---

## G3 - Completeness probe tolerance

- Date: 2026-09-25
- Branch: gate/G3-probe-tol (local only; never pushed, never merged to main)
- Base: 6e35327 "docs: add G1 gate entry to the review ledger" (gate/G1-arbiter tip)
- Commits:
  - 0ec58ae "tests: add G3 F2 completeness-probe false-alarm regression (fails pre-fix)"
  - 5abd47f "intersection: per-curve tolerance and coverage matching in completeness probe"
- Environment: ~/workspace/brep-booleans/.venv (CPython 3.12.3),
  numpy==2.5.3, manifold3d==3.5.3, cadquery-ocp==8.0.1.0.0 (OCCT 8.0.1).

### I3 tolerance derivation (required by the plan; no constant tuned to the test)

Quantities from OCCT per entity: `ic.Tolerance()` on each IntTools_Curve
(the 3D deviation OCCT's intersector reports for that exact raw curve
approximation). `ic.TangentialTolerance()` is recorded per component but not
folded into the positional match tolerance: it governs OCCT's
tangential-contact classification, not 3D positional deviation.

Existing pipeline quantities: `base_tol` (pipeline base tolerance, 1e-7)
and `max_verify_tol` (max verify_tolerance over the verified Section edges
of the pair; the positional accuracy the Section verifier itself certifies).

Per-curve match tolerance: tol_i = max(16*base_tol, 4*max_verify_tol,
2*ic.Tolerance()).

- 16*base_tol: the pre-G3 probe floor (old code: max(16*base_tol,
  4*max_verify_tol)); kept unchanged, so no existing floor is loosened.
- 4*max_verify_tol: the raw curve and the verified edge are two independent
  approximations of the same true branch, each good to about
  max_verify_tol. The 4x covers the triangle inequality across the two (2x)
  plus slack for verify_tolerance bounding distance-to-surface (which
  amplifies into 3D position error at shallow intersection angles) and
  sampling slack. This term dominates in F2 (4*1e-5 = 4e-5) and is the same
  dominant term the old probe used, so tol_i's leading term is unchanged
  from pre-G3.
- 2*ic.Tolerance(): OCCT's per-curve positional deviation claim, doubled
  because the deviation is two-sided around the true curve. Minor term in
  practice (F2 gap component: 2e-7; t9 omitted loop: ~1e-14).

95% coverage with 21 samples: 20/21 = 95.24%, so exactly one borderline
sample (trim-classifier edge effect or a locally coarse verified-edge
polyline) does not condemn a component, while a systematic offset (2+ of 21
samples off) fails. The absolute cap (every sample within 8*tol_i) ensures a
matched component tracks the verified set along its whole length. The
sampler was bumped 17 -> 21 because with 17 samples 16/17 = 94.1% < 95%,
which would degenerate the coverage rule back into the worst-point rule the
task explicitly moves away from.

8x separation band: must sit strictly between the largest plausible
same-branch approximation offset and the smallest plausible distinct-branch
separation. Measured: same-branch worst offset 9.41089e-05 = 2.35*tol_i
(F2); distinct-branch separation 2.0 = 26,000*tol_i (t9 omitted torus loop).
The 8x band edge (3.2e-4 for F2, 6.15e-4 for t9) gives 3.4x headroom above
the largest observed same-branch deviation and sits 3+ orders of magnitude
below the genuine miss. Any band in (2.35x, 26000x) separates the evidence;
8x is a round conservative choice inside that wide interval, not edge-tuned
to make the test pass (the test passes for the whole interval).

Arithmetic:

- F2 (cone side face a2 vs sphere b0): tol_i = max(16*1e-7, 4*1e-5, 2*1e-7)
  = max(1.6e-6, 4e-5, 2e-7) = 4e-5; band = 3.2e-4. The component: 21 in-trim
  samples, coverage 20/21 = 0.9524 >= 0.95, max_sample_distance 9.41089e-05
  <= 3.2e-4, nearest 1.3e-15 -> matched. The old worst-point rule refused it
  (9.41e-5 > 4e-5).
- t9 (deliberately omitted torus loop, one of two edges passed):
  tol_i = max(1.6e-6, 4*1.92245e-5, 2*5.42e-15) = 7.68979e-5;
  band = 6.15183e-4. Omitted loop: coverage 0.0, nearest = max = 2.0, and
  2.0 > 6.15e-4 -> SectionCompletenessMismatch, refused.
- Case separation: 2.0 / 9.41089e-05 = ~21,000x (4.3 orders of magnitude,
  as the plan states).

Accepted residual (I1 analysis): the separation check uses nearest
distance, so a genuinely dropped branch passing within 8*tol_i of a
verified edge at one or more samples would be recorded as approximation_gap
and accepted. Both the Section edges and the raw curves derive from the
same IntTools curves, so a raw component inside the band of a verified edge
is geometrically the same branch through a coarser approximation, not a
distinct loop (a distinct loop, t9-style, sits orders of magnitude
farther). Near-coincident branches are the near-tangent regime, which the
Section verifier flags separately (near_tangent risk flags,
tangent_sin_tol). Every approximation_gap acceptance is logged with full
numbers (coverage, nearest, max, tol_i, band) in
stages.intersection.approximation_gaps for audit. The accepted geometry is
always the verified edge set, which already represents the branch.

### What was built

- `_raw_intersector_completeness_probe` (src/brepkernel/intersection.py)
  reworked; now returns CompletenessProbeReport with raw_curve_count,
  trimmed_components, unmatched_components (separation-failing only),
  approximation_gaps, per-component records (components), max_distance.
  Per-curve tol_i as derived above; match by coverage (>= 95% of in-trim
  samples within tol_i AND every sample within 8*tol_i); separation check
  (nearest < 8*tol_i -> approximation_gap, recorded and accepted; else
  SectionCompletenessMismatch whose message logs per-component
  nearest/max distances, tol_i, band, coverage, curve tolerance).
  Non-finite distance evaluations and the no-verified-edges case still
  refuse (conservative, as pre-G3).
- The F2 component takes the coverage-match path (verdict "matched"), not
  the approximation_gap path: 20/21 samples within tol_i is exactly the
  "match by coverage, not just worst point" behavior the task specifies.
  The gap path is the backstop for components failing coverage but inside
  the band.
- New report plumbing: FaceIntersectionResult / ModelIntersectionResult
  gain approximation_gaps and completeness_components; pipeline report
  stages.intersection gains approximation_gaps, approximation_gap_count,
  completeness_components. Refusals keep kind=SectionCompletenessMismatch
  with the stage report attached (I2).
- tests/test_g3_f2_regression.py: F2 repro committed FIRST as failing
  (I4), then updated to pin the post-fix contract (acceptance, OCCT volume
  agreement, arbiter, per-component coverage record).

### Deviations from the plan

- D1 (G3): sampler bumped 17 -> 21 in-trim probe samples (justified in the
  derivation above; cost is a few extra classifier/distance evaluations per
  freeform pair).
- D2 (G3): the per-component record also carries verdict, face ids, and
  the tangential tolerance, beyond the "both numbers" the task requires;
  strictly more evidence in the report, no behavior change.

### Commands run

```
git checkout -b gate/G3-probe-tol gate/G1-arbiter
python tests/test_g3_f2_regression.py                       # pre-fix: FAIL (refused SectionCompletenessMismatch), exit 1
python tools/review_probes/repro_findings.py                # pre-fix baseline: F2 REFUSED raw=3 trimmed=2 max_distance=9.41089e-05 tol=4e-05
python tests/test_g3_f2_regression.py                       # post-fix: ALL PASS, exit 0
python tests/test_freeform_intersection.py                  # t9 still refuses SectionCompletenessMismatch
python tools/review_probes/fuzz_brep.py --trials 60 --seed 11 --nurbs --out /tmp/g3_fuzz_nurbs.json
for t in tests/test_*.py; do python $t; done                # full suite
git grep -P '\x{2014}' -- src/brepkernel/intersection.py src/brepkernel/pipeline.py \
    tests/test_g3_f2_regression.py docs/REVIEW_LEDGER.md    # 0 matches: no em dashes
```

### Results (trimmed)

- tests/test_g3_f2_regression.py (new, committed failing first per I4):
  pre-fix [FAIL] f2 cone/sphere intersection accepted
  (refusal kind=SectionCompletenessMismatch, exit 1); post-fix ALL PASS:
  accepted; kernel vol 0.540815496 vs OCCT oracle 0.540815493,
  rel_err=4.74e-09 (<= 1e-6); arbiter checked=299 kernel_errors=0;
  recorded component: verdict=matched coverage=0.9524 nearest=1.30067e-15
  max=9.41089e-05 tol_i=4e-05 band=0.00032.
- tests/test_freeform_intersection.py: ALL PASS; t9 omitted torus loop
  still refuses kind=SectionCompletenessMismatch with the new message
  logging both distances:
  "component 1: nearest_distance=2, max_sample_distance=2,
  match_tolerance=7.68979e-05 (8x band 0.000615183), coverage=0.0000 over
  21 in-trim samples, curve_tolerance=5.42137e-15".
- NURBS fuzz (60 trials, seed 11): TALLY {accept: 59,
  refuse:InsufficientPatchWitnesses: 1}, WRONG 0, CRASH 0, exit 0.
  Zero SectionCompletenessMismatch refusals, so the
  "no refusal with max distance below 8*tol_i" criterion holds vacuously.
  Tally identical to the G0 baseline run (same seed).
- Full suite: 16/16 test files exit 0 (15 pre-existing + the new
  test_g3_f2_regression), 0 [FAIL] lines.
- No em dashes in any touched file, comment, docstring, commit message,
  or generated report field (git grep -P '\x{2014}': 0 matches).

### Gate verdicts

| Criterion | Result |
|---|---|
| F2 case accepted; OCCT volume agreement within 1e-6 relative; arbiter 0 errors | PASS (rel_err=4.74e-09; checked=299, kernel_errors=0) |
| Omitted-torus-loop negative test still refuses with SectionCompletenessMismatch | PASS (t9, both distances logged: nearest=max=2.0 vs band 6.15e-4) |
| NURBS fuzz (60 trials): 0 WRONG, no SectionCompletenessMismatch refusal with max distance below 8*tol_i | PASS (0 WRONG, 0 CRASH; zero such refusals) |
| Ledger documents the tolerance derivation (I3) | PASS (derivation above) |
| Full suite green (I5) | PASS (16/16 exit 0) |
| I1-I9 invariants held | PASS (no crash; refusal always typed; tolerance derived not tuned; test committed failing first; local branch only, no push/merge) |

G3 gate: CLOSED.

## G4 - Completeness probe coverage for analytic pairs

- Date: 2026-09-25
- Branch: gate/G4-probe-coverage (local only; never pushed, never merged to main)
- Base: ca95880 "docs: add G3 gate entry to the review ledger"
- Commits:
  - 372171e "intersection: run completeness probe for non-coaxial quadric and torus pairs" (implementation worker: `_pair_needs_completeness_probe` + coaxiality helpers, `section_face_pair` comment fix, tests/test_g4_probe_coverage.py, 37/37 checks)
  - 6d5ea06 "intersection: remove temporary G4 probe-scope measurement switch" (finishing: env switch removed after measurement, probe scope unconditional)
- Environment: ~/workspace/brep-booleans/.venv (CPython 3.12.3),
  numpy==2.5.3, manifold3d==3.5.3, cadquery-ocp==8.0.1.0.0 (OCCT 8.0.1).

### Pair-type classification (from the code comments / test CASES)

Probe ON wherever OCCT intersects the pair with the numerical walking
intersector (IntPatch); probe OFF where OCCT solves in closed form (IntAna).
In doubt the probe is ON (a false refusal is allowed, a missed dropped
branch is not, I1).

| Pair | Probe | Reason |
|---|---|---|
| plane/plane, plane/cylinder, plane/sphere, plane/cone | OFF | closed form via IntAna |
| cylinder/cylinder coaxial, cylinder/cone coaxial, cylinder/sphere coaxial | OFF | closed-form circle intersection |
| sphere/sphere (any offset) | OFF | trivially coaxial (radical plane always exists) |
| cylinder/cylinder offset or tilted | ON | IntPatch walking |
| cone/cone tilted, cone/cylinder tilted | ON | IntPatch walking |
| cylinder/sphere offset | ON | IntPatch walking |
| plane/torus | ON | torus is never closed form |
| torus/torus, torus/cylinder | ON | IntPatch walking |
| any pair with a freeform (BSpline) face | ON | pre-G4 behavior, unchanged |
| any other surface type | ON | in doubt, probe on |

Coaxiality is conservative by design and uses existing quantities only:
base_tol and the two face tolerances for the linear check,
Precision.Angular for the direction check (modulo sign); the sphere side
contributes its center only (its stored axis direction is arbitrary).
A missed coaxial verdict costs only an extra probe run, while a wrong one
would skip the probe where OCCT walks (I1).

### What was measured

Analytic fuzz, 150 trials, seed 7, kinds box,cyl,sph,cone,torus, 300-point
membership audit per accept, run twice on the same seed: once with the G4
scope (probe ON for the walking group), once with the temporary env switch
BREPKERNEL_G4_MEASURE_LEGACY_PROBE_SCOPE=1 forcing the pre-G4 scope
(probe off for analytic pairs, freeform only). The env switch was removed
right after (6d5ea06); grep confirms no temp switches remain in
src/tools/tests/docs.

Tallies:

| Run | accept | refuse:InsufficientPatchWitnesses | refuse:PatchClassificationInconsistent | refuse:SectionToleranceTooLoose | refuse:SectionCompletenessMismatch | WRONG | CRASH |
|---|---|---|---|---|---|---|---|
| ON (G4 scope) | 142 | 5 | 1 | 1 | 1 | 0 | 0 |
| OFF (legacy scope) | 143 | 5 | 1 | 1 | 0 | 0 | 0 |
| G0 baseline | 143 | 5 | 1 | 1 | 0 | 0 | 0 |

Per-trial comparison: ON differs from OFF and G0 at exactly one trial
(trial 51); OFF matches the G0 baseline per trial everywhere. All 7 G0
baseline refusals reproduce at the same trials in the ON run.

Median per-op time from the fuzz `seconds` field:
ON median 5.774 s (mean 7.055, min 0.280, max 27.719);
OFF median 5.897 s (mean 7.071, min 0.309, max 26.823);
G0 median 5.873 s (mean 7.205). Median delta ON vs OFF: -2.1%, within
noise. The probe adds no measurable median cost on this corpus; its
marginal cost shows on individual walking-group trials (e.g. trial 51
refused in 1.203 s). Medians are dominated by the slow torus accepts
(26-28 s) in both runs.

### Per-refusal investigation: trial 51 (the only new refusal)

- Config: cyl union torus, seed 7 trial 51. ON: refused
  SectionCompletenessMismatch at stage intersection. OFF: accept.
  G0: accept, volume 2.3945942216662788 vs OCCT oracle 2.3945942216662788
  (exact), audit 299 checked / 0 kernel errors / 0 disagreements.
- Logged evidence (refusal message):
  "lower-level intersector exposes 1 trimmed curve component(s) not
  represented by verified Section edges (raw=3, trimmed=1,
  max_distance=0.304838); component 1: nearest_distance=0.2929,
  max_sample_distance=0.304838, match_tolerance=8e-06 (8x band 6.4e-05),
  coverage=0.0000 over 4 in-trim samples, curve_tolerance=1e-06".
- Verdict: TRUE CATCH. The unmatched raw component sits 0.2929-0.3048
  away from every verified Section edge, about 4,700x above the 8x
  separation band (6.4e-05), with coverage 0.0. The G3 false-alarm
  criterion is max distance below 8*tol_i; here it is orders of
  magnitude above, so this is a genuinely unrepresented raw branch, not
  a tolerance artifact. The G0 accept was geometrically correct (exact
  oracle volume, clean audit), so the unrepresented branch did not
  change the final solid: the refusal is conservative, which I1
  explicitly allows, and the probe behaved exactly as G4 intends for a
  torus pair. NO false alarm: G3 stays closed, no tolerance constant
  touched (I3).
- Other new refusal kinds: none.

### Deviations from the plan

- D1 (G4 finishing): the parent's original background fuzz was killed
  mid-run (trial ~125/150 of ON) by a VM/service restart that also wiped
  /tmp; the finishing worker re-ran both commands identically (same seed
  7, same output paths) detached with nohup, and both exited 0. The
  partial first ON run had already shown trial 51 refusing
  SectionCompletenessMismatch, and the re-run reproduced it at the same
  trial, confirming determinism. Completed JSONs were mirrored to
  ~/workspace/g4_fuzz_on.json and ~/workspace/g4_fuzz_off.json because
  /tmp can be wiped.
- D2 (G4): the task counted 16 test files; the suite now has 17 (the new
  tests/test_g4_probe_coverage.py). All 17 exit 0.

### Commands run

```
# ON then OFF, sequential in one detached runner (/tmp/g4_run.sh), 16:30-17:05 UTC
~/workspace/brep-booleans/.venv/bin/python tools/review_probes/fuzz_brep.py --trials 150 --seed 7 --out /tmp/g4_fuzz_on.json > /tmp/g4_fuzz_on.log 2>&1
BREPKERNEL_G4_MEASURE_LEGACY_PROBE_SCOPE=1 ~/workspace/brep-booleans/.venv/bin/python tools/review_probes/fuzz_brep.py --trials 150 --seed 7 --out /tmp/g4_fuzz_off.json > /tmp/g4_fuzz_off.log 2>&1
python3 ~/workspace/g4_analyze.py        # tallies, per-trial diff, medians
for t in tests/test_*.py; do ~/workspace/brep-booleans/.venv/bin/python "$t"; done   # full suite
git grep -P em-dash -- src/brepkernel/intersection.py tests/test_g4_probe_coverage.py docs/REVIEW_LEDGER.md   # 0 matches
grep -rni "MEASURE_LEGACY\|temporary g4\|kill switch" src tools tests docs --include=*.py --include=*.md        # 0 matches
```

### Results (trimmed)

/tmp/g4_fuzz_on.log: TALLY {'accept': 142,
'refuse:InsufficientPatchWitnesses': 5,
'refuse:SectionCompletenessMismatch': 1,
'refuse:PatchClassificationInconsistent': 1,
'refuse:SectionToleranceTooLoose': 1}
/tmp/g4_fuzz_off.log: TALLY {'accept': 143,
'refuse:InsufficientPatchWitnesses': 5,
'refuse:PatchClassificationInconsistent': 1,
'refuse:SectionToleranceTooLoose': 1}
Runner: ON exit: 0; OFF exit: 0 (branch gate/G4-probe-coverage tip 372171e).

- Full suite: 17/17 test files exit 0, 0 [FAIL] lines
  (test_brep_pipeline, test_degenerate, test_freeform_assembly,
  test_freeform_intersection, test_freeform_nurbs, test_freeform_split,
  test_g1_negative_corruption, test_g3_f2_regression,
  test_g4_probe_coverage, test_metamorphic, test_multishell_semantics,
  test_nurbs_adversarial_corpus, test_nurbs_boolean_end_to_end,
  test_regression, test_same_domain, test_step_nurbs_multiface,
  test_stress).
- No em dashes in any touched file, comment, docstring, commit message,
  or ledger entry. No temp switches remain.

### Gate verdicts

| Criterion | Result |
|---|---|
| Analytic fuzz (150): 0 WRONG | PASS (ON: 0 WRONG, 0 CRASH; OFF: 0 WRONG, 0 CRASH) |
| New SectionCompletenessMismatch refusals investigated one by one, true catch vs false alarm | PASS (exactly one new refusal, trial 51 cyl union torus: TRUE CATCH, evidence cited above; no false alarm, G3 stays closed) |
| Median time increase recorded | PASS (ON median 5.774 s vs OFF 5.897 s; no increase, within noise) |
| Full suite green (I5) | PASS (17/17 exit 0) |
| I1-I9 invariants held | PASS (no crash; refusals typed with kind= and stage; no G3 constant touched; env switch removed; local branch only, no push/merge) |

G4 gate: CLOSED.

## G3 REWORK - Per-interval adaptive matching (no blanket acceptance)

- Date: 2026-09-25
- Branch: gate/G3G4-rework (local only; never pushed, never merged)
- Base: 79e3079 "docs: close the G4 gate in the review ledger"
  (gate/G4-probe-coverage tip)
- Environment: ~/workspace/brep-booleans/.venv (CPython 3.12.3),
  numpy==2.5.3, manifold3d==3.5.3, cadquery-ocp==8.0.1.0.0 (OCCT 8.0.1).

### Why (the review finding)

The G3 probe's match rule said: a raw IntTools component is "matched" if
its kept ends match a verified edge within tol_i, OR the whole component
is within the trim. The "kept ends match" clause was a blanket acceptance:
if the two ends were fine, the interior was never checked. A raw curve
could omit (or invent) an interior loop while its ends stayed put, and the
probe would accept it. The I4 partial-omission test (t1) proved this: the
old code accepted a 50%-truncated raw curve via the approximation_gap
blanket.

### What was built

`_match_raw_component` (intersection.py) replaces the blanket rule with
per-interval adaptive matching:

- Each raw component is split into N0=20 initial parameter intervals.
- Each interval is checked with a 3-point stencil (endpoints + midpoint,
  deduplicated). A sample is a strict match iff its distance to the
  verified edges is <= tol_i AND it classifies in-trim at trim_tol.
- An interval where any stencil sample violates the strict match is
  subdivided (binary split) to max_depth=5. Every leaf interval must be
  accounted for: either all its stencil samples strictly match, or (in
  end intervals only) the violating sample satisfies the boundary-tail
  provision.
- The boundary-tail provision (keeps F2 accepted): in an end interval, a
  sample failing the strict match is accounted for iff (a) the nearest
  edge point is an edge END (within 64*eps*scale of a polyline vertex
  that is a topological edge end), (b) the sample is OUT at tight_tol
  (1e-7-ish, the pipeline base tolerance) but in/on at trim_tol, and (c)
  that edge end is within trim_tol of the trim wires. This is the F2
  shape: raw curve end overshoots the trim by 9.4e-05 along the curve
  (transverse 4.6e-09), edge end ON the trim.
- If any leaf violates at max depth, or if there are no verified edges
  at all for a material raw component, the probe raises
  SectionCompletenessMismatch (typed refusal, kind= and stage=).
- Per-t caches avoid redundant BRepExtrema_DistShapeShape and
  BRepClass_FaceClassifier calls across the stencil and subdivision.

`RawIntersectionToleranceTooLoose` (typed refusal): if
`ic.Tolerance()` on any raw IntTools curve exceeds the shared
`_section_tolerance_ceiling()` (max(128*base_tol, 1e-8*scale)), the probe
refuses instead of matching against a sloppy raw curve. The ceiling is
the same one `_verify_section_edge` uses (refactored to call the shared
helper); `_raw_curve_tolerance()` is a seam so tests can inflate it
artificially.

`approximation_gaps` removed from CompletenessProbeReport,
FaceIntersectionResult, ModelIntersectionResult, and the pipeline report
dict. Component records now carry per-interval numbers
(leaf_intervals with per-leaf verdict/t_start/t_end/depth/max_distance,
leaf_interval_count, max_depth_reached, boundary_tail_intervals,
tolerance_ceiling).

### I3 tolerance derivation (no constant tuned to the test)

- tol_i = max(16*base_tol, 4*max_verify_tol, 2*curve_tol): UNCHANGED from
  G3. The per-interval rule uses the same tolerance; only the matching
  logic (blanket -> per-interval) changed.
- trim_tol: the face trim tolerance from the pipeline (4e-5 in F2).
- tight_tol: max(base_tol, face tolerances) (1e-7 in F2). The boundary
  band is (tight_tol, trim_tol]: OUT at tight_tol but in/on at trim_tol.
- 64*eps*scale for edge-end detection: 64 multiples of double epsilon
  times the model scale; a numerical guard, not a geometric tolerance.
- N0=20, max_depth=5: sampling density parameters, not tolerances. Depth
  5 gives 32x refinement (20*32=640 max leaves); adequate to isolate a
  partial omission (the t1 test produces 320 violating leaves at depth
  5, proving the omission is localized, not a sampling artifact).

### I4 failing-tests-first

- tests/test_g3_probe_rework.py::t1 (partial omission): FAILS on pre-fix
  code (old code accepts via approximation_gap blanket), PASSES post-fix
  (SectionCompletenessMismatch, 320 violating leaves at depth 5).
- tests/test_g3_probe_rework.py::t2 (inflated raw tolerance):
  FAILS on pre-fix code (no hook to inflate), PASSES post-fix
  (RawIntersectionToleranceTooLoose when _raw_curve_tolerance is
  monkeypatched to 1.0).
- Both verified by stashing src/ and running (2026-09-25).

### F2 (the false alarm that must stay accepted)

- tests/test_g3_f2_regression.py: NURBS cone INTERSECT NURBS sphere.
  Accepted. Volume 0.540815496 vs OCCT 0.540815493 (rel 4.74e-09, within
  1e-6). Independent membership audit: 299 checked, 0 errors.
- Probe records: 4 components, 1 above strict tol, matched via
  boundary-tail provision (1 tail leaf [0.95, 1.0], maxd=9.41089e-05,
  depth 0, no subdivision needed).

### Commands run

```
~/workspace/brep-booleans/.venv/bin/python tests/test_g3_probe_rework.py
~/workspace/brep-booleans/.venv/bin/python tests/test_g3_f2_regression.py
# pre-fix verification (2026-09-25):
git stash push -- src/ && ~/workspace/brep-booleans/.venv/bin/python tests/test_g3_probe_rework.py  # t1,t2 FAIL
git stash pop
```

### Results (trimmed)

test_g3_probe_rework.py: t1 PASS (SectionCompletenessMismatch, 320
violating leaves), t2 PASS (RawIntersectionToleranceTooLoose).
test_g3_f2_regression.py: ALL PASS (5/5 checks).

### Gate verdicts

| Criterion | Result |
|---|---|
| I4: partial-omission negative test fails pre-fix, passes post-fix | PASS |
| I4: inflated-tolerance negative test fails pre-fix, passes post-fix | PASS |
| F2 accepted, volume within 1e-6 of OCCT, arbiter clean | PASS |
| No tolerance loosened (I3); tol_i formula unchanged | PASS |
| Full suite green (I5) | PARTIAL (4/18 pass; see G4 REWORK entry) |

## G4 REWORK - Probe runs for every pair (classification deleted)

- Date: 2026-09-25
- Branch: gate/G3G4-rework (local only; never pushed, never merged)
- Base: 79e3079 (gate/G4-probe-coverage tip)
- Environment: ~/workspace/brep-booleans/.venv (CPython 3.12.3),
  numpy==2.5.3, manifold3d==3.5.3, cadquery-ocp==8.0.1.0.0 (OCCT 8.0.1).

### Why (the review finding)

G4 classified face pairs and skipped the completeness probe for
"closed-form" pairs (plane/plane, plane/cylinder, coaxial quadrics),
claiming OCCT's IntAna exact intersector needs no verification. That is
using a claim about OCCT's exactness to skip verifying OCCT: the probe
exists because we do not trust the intersector, so exempting the pairs
we trust most is backwards. The classification (quadric axis extraction,
coaxiality checks) was also new code that could itself be wrong, and a
wrong "coaxial" verdict would silently skip the probe where OCCT walks.

### What was built

- Deleted `_pair_needs_completeness_probe` and all classification
  helpers (`_quadric_axis_point`, `_axes_coaxial`, `_quadrics_coaxial`,
  `_QUADRIC_TYPES`, `_PLANE_OR_QUADRIC`): ~120 lines removed.
- The probe now runs unconditionally for every candidate pair that
  produces section curves, analytic and freeform alike.
- `BREPKERNEL_COMPLETENESS_PROBE=0/off/false/no` environment switch
  (`_completeness_probe_enabled()`): honest kill switch for measurement
  and emergency use, honored in both `section_face_pair` and
  `intersect_models`. Not a temporary hack; documented in the test.
- tests/test_g4_probe_coverage.py rewritten: asserts the probe runs
  (raw_curve_count > 0) and accepts (zero unmatched) for plane/plane,
  plane/cylinder, tilted cyl/cyl, and torus/cyl; asserts the env switch
  defaults on, disables with =0 (probe skipped, raw_curve_count == 0),
  and re-enables on restore.

### Commands run

```
~/workspace/brep-booleans/.venv/bin/python tests/test_g4_probe_coverage.py
~/workspace/brep-booleans/.venv/bin/python tools/review_probes/fuzz_brep.py --trials 150 --seed 7 --out /tmp/g4r_fuzz_on.json
BREPKERNEL_COMPLETENESS_PROBE=0 ~/workspace/brep-booleans/.venv/bin/python tools/review_probes/fuzz_brep.py --trials 150 --seed 7 --out /tmp/g4r_fuzz_off.json
for t in tests/test_*.py; do ~/workspace/brep-booleans/.venv/bin/python "$t"; done
git grep -P '\x{2014}' -- src tests docs  # no em dashes
```

### Results (trimmed)

test_g4_probe_coverage.py: 12/12 PASS (probe runs and accepts for all
four analytic pair kinds; env switch on/off/restore).

Fuzz (150 trials, seed 7): NOT COMPLETED. The fuzz was started but trial
0 took 126.78s (baseline median 5.7s) due to severe CPU contention: the
coordinator was running parallel suites in other worktrees (load average
8-9 on 2 CPUs). At 127s/trial, 150 trials would take 5+ hours. The fuzz
was killed after trial 0 (accept, 0 WRONG). The full 150-trial on/off
comparison must be run when the system is quiet; recorded here as an open
item.

Full suite: 4/18 test files verified PASS (the probe-critical ones):
test_g3_probe_rework.py (t1, t2), test_g4_probe_coverage.py (12/12),
test_freeform_intersection.py (t6, t7, t8, t9), test_g3_f2_regression.py
(5/5). The remaining 14 test files were not run to completion due to the
same CPU contention (test_brep_pipeline.py alone exceeded 6 minutes under
load). They must be run when the system is quiet; recorded as an open
item.

### Per-refusal investigation

No new SectionCompletenessMismatch was observed in the completed work:
- t1 (partial omission): 320 violating leaves at depth 5, refused as
  designed (true catch, the I4 negative test).
- t9 (omitted torus loop): 640 violating leaves at max depth 5 (all
  intervals violate), refused as designed (true catch, pre-existing
  negative test still valid under the stricter rule).
- F2: accepted via boundary-tail provision (1 tail leaf), not a refusal.

The G4-baseline trial 51 (cyl union torus, TRUE CATCH under the old rule)
was not re-run in the fuzz; its configuration should be checked manually
when the fuzz runs.

### Gate verdicts

| Criterion | Result |
|---|---|
| Analytic fuzz (150): 0 WRONG | OPEN (not completed; system load) |
| New SectionCompletenessMismatch refusals investigated one by one | PASS (t1, t9 are true catches; F2 accepted via provision) |
| Median time increase recorded | OPEN (not measured; system load) |
| Full suite green (I5) | PARTIAL (4/18 pass; 14 not run due to load) |
| I1-I9 invariants held | PASS (no crash; refusals typed; no G3 constant touched; local branch only) |
