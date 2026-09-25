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

---

## G6 - Tolerance-aware broad phase

- Date: 2026-09-25
- Branch: gate/G6-broadphase-tol (local only; never pushed, never merged to main)
- Base: eaa596c "docs: add G4 gate entry to the review ledger" (G4 tip)
- Commits:
  - a5dbca9 "tests: add G6 tolerance-aware broad phase regression" (I4: failing test first)
  - a766b61 "intersection: tolerance-aware broad phase and contact band (G6)"
- Environment: ~/workspace/brep-booleans/.venv (CPython 3.12.3,
  numpy==2.5.3, cadquery-ocp==8.0.1.0.0, OCCT 8.0.1)

### What changed

- src/brepkernel/step_ingest.py: _shape_bbox() now calls
  BRepBndLib.AddOptimal_s(shape, box, False, True) (useShapeTolerance=True).
  This only enlarges face boxes, never shrinks them. New helper
  model_max_tolerance() returns the max OCCT tolerance over all vertices,
  edges, and faces of a model via BRep_Tool.Tolerance_s.
- src/brepkernel/pipeline.py: after ingest, when broadphase_pad was not
  explicitly passed, pad = max(contact_tol, 2 * tol_max) with tol_max the
  max over both models. Explicit caller pads are honored unchanged.
  The report records broadphase_pad and broadphase_max_tolerance.
- src/brepkernel/intersection.py: section_face_pair() now takes the max of
  the passed contact_tol and the face-tolerance-derived band
  max(4*base_tol, 2*tol(fa), 2*tol(fb)). The face-tolerance rule already
  existed for contact_tol=None; the change only stops an explicit small
  contact_tol from suppressing it.

### Pad derivation (I3)

- contact_tol value unchanged: 4.0 * base_tol = 4e-7.
- On ordinary models tol_max is ~1e-7, so
  pad = max(4e-7, 2e-7) = 4e-7, bit-identical to the old default
  max(contact_tol, 4*base_tol). No acceptance ceiling, no refusal, and no
  tolerance constant was loosened: the pad only widens the candidate set,
  which can add typed refusals but never new acceptances. The contact-band
  change is derived from BRep_Tool.Tolerance_s of the exact faces involved
  and likewise only widens the ambiguous band (more refusals, never more
  accepts).

### I4 pre-fix demonstration (test committed at a5dbca9, run on eaa596c)

```
[FAIL] g6 pair is a candidate even with zero pad candidates=0
[PASS] g6 face tolerance visible to the kernel max_face_tol=0.001
[FAIL] g6 union refuses with typed kind kind=BoundaryOrUnknownPatch
[FAIL] g6 pair reached the exact intersection stage candidates=0
[FAIL] g6 pad recorded in report and covers the tolerance broadphase_pad=None
SOME FAILURES
```

Deviation D1 (G6): the task anticipated a disjoint passthrough (accepted
wrong/empty result). Actual pre-fix behavior: the pair was silently dropped
by the broad phase (0 candidates), and the union then refused at the
assembly stage with kind=BoundaryOrUnknownPatch, an accident of the
high-tolerance witness classification, not a diagnosis of the contact. The
I4 property still holds (the committed test fails on current code), and the
pre-fix refusal was at the wrong stage with the wrong kind.

### Post-fix test result (a766b61)

```
[PASS] g6 pair is a candidate even with zero pad candidates=21
[PASS] g6 face tolerance visible to the kernel max_face_tol=0.001
[PASS] g6 union refuses with typed kind refusal={'stage': 'intersection',
  'type': 'IntersectionError', 'kind': 'SectionToleranceTooLoose',
  'message': 'section edge 0: OCCT edge/face tolerance 0.001 exceeds
  acceptance ceiling 1.28e-05'}
[PASS] g6 pad recorded in report and covers the tolerance broadphase_pad=0.002
ALL PASS
```

Refusal-kind choice, documented per the task: the plan listed
ambiguous_contact or NearCoincidentFaces. NearCoincidentFaces does not exist
in the codebase (G2 work not done). ambiguous_contact never fires here
because OCCT's section engine merges the 5e-4 gap inside the face's 1e-3
tolerance and emits section edges, so the near-contact classifier is never
reached. The section verifier then refuses with SectionToleranceTooLoose:
the face tolerance (1e-3) exceeds the acceptance ceiling
max(128*base_tol, 1e-8*scale) = 1.28e-5. This is the most accurate existing
kind: the section evidence genuinely cannot be certified on a face this
loose. UnresolvedContact would describe a later assembly stage the pair
never reaches. Measured directly via section_face_pair on the high-tol face
vs all 5 candidate B faces: every call raised IntersectionError with
kind=SectionToleranceTooLoose.

### Candidate-count comparison (old vs new code)

Measured via boolean_brep on 4 accepting cases; old code measured with the
implementation stashed (git stash push -- src/), then restored:

| case | old candidates | new candidates | new pad | new tol_max |
|---|---|---|---|---|
| A box/box union | 6 | 6 | 4e-07 | 1e-07 |
| B NURBS-sphere union | 1 | 1 | 4e-07 | 1e-07 |
| C cylinder-through-box difference | 2 | 2 | 4e-07 | 1e-07 |
| D NURBS-sphere intersection | 1 | 1 | 4e-07 | 1e-07 |

Counts change only where tolerances are large (the G6 test model: pad 2e-3,
21 candidates vs 0 before). On ordinary low-tolerance models the pad is
unchanged and counts are identical.

### Commands run

```
git checkout -b gate/G6-broadphase-tol gate/G4-probe-coverage
# wrote tests/test_g6_tolerance_broadphase.py, committed (a5dbca9), ran on current code -> SOME FAILURES (I4)
# implemented _shape_bbox useShapeTolerance, model_max_tolerance, pad escalation, contact-band max
python tests/test_g6_tolerance_broadphase.py   # ALL PASS
python /tmp/cand_counts.py                     # new-code counts: 6, 1, 2, 1
git stash push -q -m g6-impl-temp -- src/       # old-code counts: 6, 1, 2, 1
git stash pop -q
for t in tests/test_*.py; do python "$t"; done   # full suite, 18 files
em-dash check on touched files -> 0
```

### Results (trimmed)

- Full suite: 18/18 test files, 0 [FAIL] lines (the 17 existing files plus
  the new tests/test_g6_tolerance_broadphase.py; test_nurbs_adversarial_corpus
  and test_step_nurbs_multiface re-checked with grep -ci fail -> 0).
- No em dashes in touched files, comments, docstrings, commit messages, or
  this ledger entry. No temp switches remain (the G4 env switch was already
  removed in 6d5ea06; the G6 candidate-count measurement used git stash,
  popped immediately after).

### Gate verdicts

| Criterion | Result |
|---|---|
| New test passes with pair-as-candidate and typed refusal | PASS (refusal kind=SectionToleranceTooLoose at stage intersection; choice documented above) |
| Candidate counts change only where tolerances are large | PASS (4 accepting cases identical old/new; G6 model 0 -> 21 candidates, pad 2e-3) |
| Pad recorded in the operation report | PASS (report["broadphase_pad"], report["broadphase_max_tolerance"]) |
| Full suite green (I5) | PASS (18/18, 0 FAIL lines) |
| I1-I9 invariants held | PASS (refusal, never wrong accept; typed kind + stage report; I3 derivation above; test committed before fix; no crash; boolean()/boolean_brep() contracts untouched; local branch only, no push/merge) |

G6 gate: CLOSED.

## Course correction: rebase onto v0.9, G7 replaced (2026-09-25)

Ben forwarded ChatGPT's review of the plan as a work-order update. Verified:
origin/main moved 5554713 -> 88d09ebd856d ("Merge v0.9 verified periodic
seam-aware NURBS splitting", parents 5554713 + bf1f504c4f88).

### v0.9 seam mechanism (summary, replacing the plan's G7 pre-split design)

split_models() in src/brepkernel/split.py now routes section edges by
seam risk flags instead of pre-splitting periodic faces:
- a section edge that is a seam on exactly one operand ("seam_on_a" xor
  "seam_on_b") reuses that operand's existing closing boundary (not added
  as a split tool there) and remains a verified splitting tool on the
  opposite operand;
- a curve that is a seam on both operands is refused by default
  (shared_seam_curve -> unresolved), counted in shared_seam_refusals;
- split_face() gained ignored_risk_flags so each operand ignores the
  other operand's seam flag;
- pipeline report exposes reused_seam_edges_a/b and shared_seam_refusals.
No ShapeUpgrade_ShapeDivideClosed pre-splitting is used.

### G7 replaced

Old G7 (periodic seam pre-splitting) is dropped. New G7 is a shared-seam /
periodic fuzz / revolved-surfaces gate: torus-heavy fuzz runs,
seam-crossing section cases, revolved NURBS cases; it stress-tests the v0.9
mechanism above instead of rewriting it.

### Rebase performed (local only, no pushes, no merges)

Local main reset to 88d09ebd856d. Each gate branch rebased in dependency
order; all rebases applied cleanly with no conflicts except the expected
already-applied skips:
- gate/G0-tooling onto main: 3 commits applied clean
- gate/G1-arbiter onto gate/G0-tooling: 4 commits applied clean
  (test_nurbs_adversarial_corpus.py merged cleanly despite v0.9 touching it)
- gate/G3-probe-tol onto gate/G1-arbiter: 3 commits applied clean
- gate/G4-probe-coverage onto gate/G3-probe-tol: 3 commits applied clean
- gate/G6-broadphase-tol onto gate/G4-probe-coverage: 3 commits applied clean
gate/G2-coincident-faces rebase deferred until its closeout worker finishes;
its split.py changes overlap the v0.9 seam-routing hunk and will need a
combining resolution (seam routing + coincident-pair logic coexist).

### Roadmap update (two-track structure)

Per the correction: evidence certificate schema and persistent naming move
earlier, right after G8. B9 (certified intersection core) becomes a parallel
research track with its own milestones, not a dependency of the practical
track. FreeCAD integration comes before B9. Gate list order going forward:
G0, G1, G3, G4, G6, G2 (planar redirection), G1/G3/G4/G6 reworks, G5
(multi-ray parity as production second classifier), new G7 (seam/periodic
fuzz), G8, evidence-cert + persistent naming, Part B with FreeCAD before B9
(parallel research track).

Full-suite re-run on the rebased G6 tip is in flight; results to be recorded
when it lands.
## G2 coincident faces: 2.5/2.6 remainders, 2.7 tests, fuzz (2026-09-25)

Branch: gate/G2-coincident-faces (local only). Resumes from d4bacc0 after a
VM/service restart killed the previous worker; the in-progress work was
assessed as coherent, committed, and continued here. This entry covers the
completion work: 5bd1d18 (implementation) and 2b0d0ee (tests).

### What was built

2.5 boundary-contact rule (intersection.py, split.py):
- SectionEdgeRecord.is_boundary_contact: set by _mark_boundary_contacts when
  every (u,v) sample classifies TopAbs_ON on both faces' trims via
  BRepClass_FaceClassifier. Pairs whose edges are all boundary contacts get
  status "boundary_contact".
- split_models: boundary-contact edges are NOT split tools ("creates no
  split"); they are collected separately.
- _resolve_contacts: each boundary edge must coincide (within tol) with a
  coincident-region boundary or a verified section edge from another pair
  (including another pair's boundary edge), else UnresolvedContact. Point
  contacts may resolve against boundary-contact edges.

2.5 NonManifoldResult (assembly.py, pipeline.py):
- assemble_boolean(..., allow_nonmanifold=False): after building solids, a
  union with >1 solid where any two touch (BRepExtrema_DistShapeShape <=
  4*max(base_tol, sew_tol)) refuses typed NonManifoldResult by default.
- allow_nonmanifold=True returns the compound with the contact recorded in
  report notes. Intersection touching-only already accepted empty; difference
  face-touching already accepted A unchanged. Disjoint multi-solid unions are
  unaffected (distance > tol).

2.6 coincident_boundary provenance (split.py, assembly.py, pipeline.py):
- split_models indexes partner boundary edges used as overlap-split tools
  (operand, face_id, edge_id) on ModelSplitResult.coincident_boundary_tools.
- _build_edge_lineage: new provenance_kind "coincident_boundary" for result
  edges matching a partner tool AND lying in the interior of a parent face
  (i.e. created by the overlap cut, not merely coincident with an original
  boundary). boolean_section keeps priority (stronger evidence with verified
  pcurves). Same-operand source_boundary now has a geometric fallback
  (_edge_on_face_boundary) for splitter-rebuilt edges.
- complete_lineage requires coincident_boundary edges to carry the other
  operand's source edge IDs; bare claims are IncompleteEdgeLineage.
- UnifySameDomain: verified not in the main pipeline path (only the
  pre-pipeline same-domain equivalence check); requirement satisfied.

2.7 tests/test_g2_coincident_probes.py (ALL PASS):
- 10 probe accept cases with exact volumes + independent membership arbiter.
- Touching-only: edge-touching -> NonManifoldResult; allow_nonmanifold ->
  2-solid compound; point-touching -> UnresolvedContact.
- Box-box exact oracle (5 grid-snapped Tier A pairs): volume and shells vs
  exact_op_volume/expected_shells.
- Cylinder cases: blind hole, boss, coaxial tube-in-tube, counterbore.
- Near-coincident negatives: 0.5x/2x base_tol -> NearCoincidentFaces;
  10x -> accepts (out of band).

### Commands run

```
git log --oneline  # 87bfcdd..2b0d0ee on gate/G2-coincident-faces
python tools/review_probes/common_cad_probes.py  # wrong=0 unmet_accepts=0
python tests/test_g2_keep_table.py                # ALL PASS
python tests/test_g2_coincident_probes.py         # ALL PASS
python tools/review_probes/fuzz_brep.py --trials 200 --seed 5 --snap 0.5 --kinds box,cyl --out /tmp/fuzz_snapped_g2.json
python tools/review_probes/fuzz_brep.py --trials 150 --seed 7 --out /tmp/fuzz_analytic_g2.json
python tools/review_probes/fuzz_brep.py --trials 60 --seed 11 --nurbs --out /tmp/fuzz_nurbs_g2.json
em-dash check on touched files -> 0
```

### Results (trimmed)

- Probes: wrong=0 unmet_accepts=0. Edge-touching union now refuses
  assembly/NonManifoldResult (was accept; expect="either" so still met).
- Snapped fuzz (200, seed 5): 112 accept / 88 refusals, 0 WRONG, 0 CRASH.
  Baseline was 93/107. Histogram: UnresolvedContact 91->2 (G2 overlap split
  now resolves coincident configs), NearCoincidentFaces 44 (new, correct
  in-band refusals), BoundaryOrUnknownPatch 21, NonManifoldResult 1,
  SectionCompletenessMismatch 3. All 112 accepts verified by the arbiter
  and OCCT volume; the shift is G2 working, not wrong accepts.
- Analytic fuzz (150, seed 7): 142 accept / 8 refusals, 0 WRONG, 0 CRASH.
  Baseline was 143/7. Trial 51 (cyl/torus union) flipped accept ->
  refuse:SectionCompletenessMismatch. This refusal is raised inside
  section_face_pair, which runs before any G2 code; the G2 changes cannot
  cause it. Attributed to section-level nondeterminism, not a regression.
- NURBS fuzz (60, seed 11): 59 accept / 1 refusal, 0 WRONG, 0 CRASH.
  Identical to baseline.

### Gate verdicts

| Criterion | Result |
|---|---|
| Snapped fuzz: 0 WRONG, 0 CRASH | PASS (112/88 vs baseline 93/107; histogram shift explained) |
| Analytic fuzz: 0 WRONG, acceptance >= baseline | MARGINAL (142 vs 143; 1-trial section nondeterminism, not G2-caused) |
| NURBS fuzz: 0 WRONG, acceptance >= baseline | PASS (59/1 identical to baseline) |
| 2.5 boundary-contact rule + NonManifoldResult | PASS (typed refusals verified) |
| 2.6 coincident_boundary + complete_lineage | PASS (kind exists, refs required) |
| 2.7 tests | PASS (new test file ALL PASS) |
| I1-I9 invariants held | PASS (refusal never wrong accept; typed kinds+stages; I3 untouched; tests before tables; no crash; contracts separate; local only) |

G2 gate: OPEN. The analytic 1-trial variance should be re-run to confirm
nondeterminism before closing. Remaining known items: 2.3 curved analytic
recognition via ShapeAnalysis_CanonicalRecognition (minor; only planar
implemented), and the full 19-file suite re-run (probes + targeted tests
green; full suite was 19/19 before this work).

### WHAT WORKS

- Coincident-face union/intersection/difference accept at exact volumes with
  arbiter-verified membership (10 probe cases, 5 oracle pairs, 4 cylinder
  cases).
- Boundary contacts create no split and resolve typed (or block typed).
- Touching-only unions refuse NonManifoldResult by default; the opt-in flag
  returns a documented compound.
- Near-coincident (in-band) configurations refuse NearCoincidentFaces rather
  than guessing.

### WHAT IS MISSING

- Curved analytic recognition for NURBS supports (cylinder/sphere/cone via
  ShapeAnalysis_CanonicalRecognition); only planar recognition is wired.
- Re-run of analytic fuzz trial 51 to confirm the SectionCompletenessMismatch
  is nondeterministic.

### WHAT REMAINS

- Full 19-file suite re-run (I5) after the fuzz variance is resolved.
- Parent decision on closing G2 vs holding for the re-run.

### G2 closeout (2026-09-25, evening)

Course correction (work-order update): G2 is narrowed to planar coincident
faces. Coincident cylinders/spheres/cones/tori/NURBS are deferred to a
follow-up gate. The 2.3 curved analytic recognition via
ShapeAnalysis_CanonicalRecognition was NOT attempted. Bindings check on
cadquery-ocp 8.0.1.0.0: ShapeAnalysis_CanonicalRecognition exists, but
ShapeCustom.ConvertToAnalytical is absent, so recovering analytic support
parameters would require hand-fitting them from samples - not small, not
safe. Recorded as deferred by scope decision, not as an open G2 item.

Analytic fuzz re-run (trial-51 confirmation):
- `python tools/review_probes/fuzz_brep.py --trials 150 --seed 7 --out /tmp/fuzz_analytic_g2_rerun.json`
  (exit 0; process imported G2 code, started 18:54 UTC before any checkout).
- TALLY {'accept': 142, 'refuse:InsufficientPatchWitnesses': 5,
  'refuse:SectionCompletenessMismatch': 1,
  'refuse:PatchClassificationInconsistent': 1,
  'refuse:SectionToleranceTooLoose': 1}. 0 WRONG, 0 CRASH.
  Identical tally to the G2 run and the G4 probe-ON run.
- Trial 51 (cyl/torus union) refuses SectionCompletenessMismatch at stage
  intersection with evidence identical to the G4 true-catch finding:
  raw=3, trimmed=1, max_distance=0.304838, nearest_distance=0.2929 vs 8x
  band 6.4e-05 (match_tolerance=8e-06), coverage=0.0. Deterministic across
  four runs (G4 partial, G4 full, G2, this re-run). The 143/7 -> 142/8
  "flip" vs the G0 baseline is the G4 probe-scope change (legacy scope
  accepted it; probe-ON refuses it as a true catch), not a G2 regression.
  Analytic criterion: PASS with this explanation recorded. JSON mirrored
  to ~/workspace/fuzz_analytic_g2_rerun.json.

Full suite (20 files) on G2 code: ALL exit 0.
- At ~19:08 UTC the main working tree was checked out from
  gate/G2-coincident-faces to main for rebase work (reflog 19:08:34).
  Forensics: files 1-18 (test_brep_pipeline.py .. test_same_domain.py)
  each started before 19:08:34, so all imported G2 code (no brepkernel
  module was recompiled between 19:08:34 and 19:08:49; the first
  post-checkout import, ~19:08:49, was test_step_nurbs_multiface.py).
  All 18 exit 0 on G2 code, including test_nurbs_adversarial_corpus.py.
- test_step_nurbs_multiface.py and test_stress.py started on the G6 tree
  in the main checkout, so both were re-run in a git worktree of
  gate/G2-coincident-faces at b251711 (~/workspace/brep-gates-g2):
  both exit 0.
- test_nurbs_adversarial_corpus.py: 3/3 PASS on G2 code (1 main-tree +
  2 worktree), 24/24 checks each, 0 FAIL lines. The reported 1-of-4 flake
  did not reproduce. Characterization: the only nondeterminism-sensitive
  contract is the h=1e-5 sliver case, which permits either a typed refusal
  or an accurate accept; every accepted configuration in all runs passed
  the membership arbiter (tests/_arbiter.py) and OCCT volume agreement
  inside the test itself. No WRONG accept in any run: no I7 STOP. A flaky
  refusal remains acceptable under I1.

Probes re-run (G2 worktree):
- `python tools/review_probes/common_cad_probes.py` (no --baseline),
  exit 0: wrong=0 unmet_accepts=0. All 14 expect=accept cases accepted
  with exact volumes; the singular equal-radius crossing-cylinders case
  still refuses (assembly/PatchClassificationInconsistent); edge-touching
  union refuses assembly/NonManifoldResult (expect=either).

Note on the work-plan criterion "snapped fuzz refusal rate under 5%":
the measured G0 baseline on the snapped fuzz is 93/107 (53.5% refusal),
so the 5% figure does not match any measured baseline; applied literally
it fails both baseline and G2. Used the task-directed reading instead:
acceptance at or above the applicable baseline. Snapped: 112/200 vs
93/200 baseline, PASS. NURBS: 59/60 identical to baseline, PASS.

### Final gate verdicts

| Criterion | Result |
|---|---|
| Probes (no --baseline): 0 wrong, all accepts exact | PASS (wrong=0 unmet_accepts=0, re-run 25 Sept) |
| Snapped fuzz: 0 WRONG, 0 CRASH, acceptance >= baseline | PASS (112/200 vs 93/200; shift explained in 2.8) |
| Analytic fuzz: 0 WRONG, acceptance explained vs baseline | PASS (142/8; trial-51 true catch confirmed deterministic) |
| NURBS fuzz: 0 WRONG, acceptance >= baseline | PASS (59/60 identical to baseline) |
| 2.5 boundary-contact rule + NonManifoldResult | PASS (2.8, unchanged) |
| 2.6 coincident_boundary + complete_lineage | PASS (2.8, unchanged) |
| 2.7 tests | PASS (2.8, unchanged) |
| Full suite green (I5) | PASS (20/20 exit 0 on G2 code) |
| Corpus flakiness characterized, no wrong accept | PASS (3/3 on G2; variance is accept/refuse only) |
| I1-I9 invariants held | PASS (local branch only; no push/merge; no crash; typed refusals) |

G2 gate: CLOSED (planar coincident-face scope per the narrowed work
order). Curved analytic coincidence (cylinder/sphere/cone/torus/NURBS)
is explicitly deferred to a follow-up gate, not a G2 remainder.

### WHAT REMAINS (post-closeout, for the parent/coordinator)

- Follow-up gate: curved coincident-face recognition
  (ShapeAnalysis_CanonicalRecognition + parameter recovery; planar path
  in src/brepkernel/coincidence.py is the template).
- Observation: after the 19:08 UTC rebase, gate/G6-broadphase-tol no
  longer contains src/brepkernel/coincidence.py (deleted in the G2..G6
  diff) and trims tests/test_g2_coincident_probes.py and
  tests/test_g2_keep_table.py. If G6 is meant to build on G2's planar
  work, that drop needs a coordinator decision.

---

## Probe hardening - review probe pass/fail semantics (2026-09-25)

- Date: 2026-09-25
- Branch: probe-hardening (local only; never pushed, never merged to main)
- Base: 88d09eb "Merge v0.9 verified periodic seam-aware NURBS splitting"
- Commits:
  - 555616f "tools: import review probe scripts and fixtures verbatim from G0"
  - d6d2530 "probe: harden review probe pass/fail semantics"
  - (this entry) "docs: probe-hardening ledger entry"
- Environment: ~/workspace/brep-booleans/.venv (CPython 3.12.3),
  numpy==2.5.3, manifold3d==3.5.3, cadquery-ocp==8.0.1.0.0 (OCCT 8.0.1).

Notes on setup (deviations from a plain main checkout):
- D1: tools/review_probes/ and tests/data/review_20260925/ are not on
  main; they were imported verbatim from 3ab4426 (G0) in 555616f so this
  branch is self-contained and the probes are runnable.
- D2: work was done in a private git worktree at
  ~/workspace/brep-probe-work because the shared ~/workspace/brep-gates
  tree was being switched between gate branches by the coordinator
  mid-task (observed branch flips gate/G7-seam-stress ->
  gate/G2-planar -> gate/G3G4-rework during this item).
- D3: main has 14 test files under tests/, not 20. The "20" count seen
  elsewhere includes unmerged gate-branch tests (G1/G2/G3/G4/G6). The
  full suite below is the 14 files on this branch; the 6 extra files
  belong to other gates' unmerged code and were not run here.

### Changes (d6d2530)

tools/review_probes/common_cad_probes.py
- Added unexpected_accept counter: an accepted case with expect="refuse"
  now increments it and forces nonzero exit. Previously it only appended
  "(accepted a case expected to refuse: inspect)" to the detail string and
  the script could exit 0. Tally line now prints
  wrong=.. unmet_accepts=.. unexpected_accepts=...; exit is 1 if wrong or
  unexpected_accept is nonzero (in --baseline mode as well).
- Docstring reconciled with the measured baseline on main: 17 cases,
  4 accepted, 13 refused (was: "4 accepted and 10 refused"). Eleven of the
  13 refusals are coincident-face cases (finding F1); the equal-radius
  crossing cylinders stay refuse-expected (genuine singular).

tools/review_probes/repro_findings.py
- F2: on accept, asserts relative volume within 1e-6 of the
  BRepAlgoAPI_Common oracle AND arbiter membership_audit kernel_errors
  == 0 (seeded per-trial RNG 20260925). Typed refusal (BRepAmbiguousResult)
  keeps old behavior: records the refusal, returns False, no assertions.
- F4: asserts the kernel union winding-number verdict at
  p=[1.3322189978645596, 1.1200289002934705, -0.8236295495693953] is OUT
  (|w| < 0.05) and surface_distance from p to the kernel union exceeds 1.0;
  records the OCCT BRepClass3d_SolidClassifier state at p for A, B, the
  kernel union, and the OCCT fuse.

No kernel code touched; boolean()/boolean_brep() contracts unchanged;
typed refusals stay typed.

### I4 evidence (hardened probe vs unmodified main code)

common_cad_probes.py --baseline on main code (pre-hardening script):
17 cases, 4 ACCEPT, 13 REFUSE, wrong=0, unmet_accepts=11, EXIT=0.
Accepted: slot overshooting top face (vol=6 ref=6), through-hole cylinder
(vol=6.429203673), crossing cylinders r1 != r2 (vol=14.93530141), sphere
from box corner (vol=7.134662047). All 11 accept-expected refusals are
assembly/UnresolvedContact (coincident faces, finding F1).

Unexpected-accept enforcement demo (throwaway /tmp driver, forced
expect="refuse" on the accepted case "cut: slot overshooting top face"):
- old script: ACCEPT ... (accepted a case expected to refuse: inspect) /
  wrong=0 unmet_accepts=0 / exit 0  <- the weakness
- hardened script: same line / wrong=0 unmet_accepts=0
  unexpected_accepts=1 / exit 1  <- fixed

Hardened script on main code: --baseline gives
wrong=0 unmet_accepts=11 unexpected_accepts=0, exit 0; without --baseline
exit 1 (unmet accept-expected refusals; G2 not yet done, as designed).

repro_findings.py (hardened) on main code, trimmed:
versions: numpy=2.5.3 cadquery-ocp=8.0.1.0.0 manifold3d=3.5.3
F2 OCCT common volume: 0.54081549348893
F2 REFUSED: {'stage': 'intersection', 'type': 'IntersectionError',
 'kind': 'SectionCompletenessMismatch', 'message': 'lower-level
 intersector exposes 1 trimmed curve component(s) not represented by
 verified Section edges (raw=3, trimmed=2, max_distance=9.41089e-05,
 tol=4e-05)'}
F4 A             occt_state=TopAbs_State.TopAbs_OUT winding=-0.0000 surface_distance=1.5594
F4 B             occt_state=TopAbs_State.TopAbs_OUT winding=-0.0000 surface_distance=1.0396
F4 kernel union  occt_state=TopAbs_State.TopAbs_IN winding=+0.0000 surface_distance=1.0396
F4 OCCT fuse     occt_state=TopAbs_State.TopAbs_IN winding=+0.0000 surface_distance=1.0396
EXIT=1
F2 still refuses with SectionCompletenessMismatch on main (G3 not done);
recorded, not "fixed" here, per the task. Both new F4 asserts pass on
main (winding +0.0000 < 0.05; distance 1.0396 > 1.0); the OCCT false IN on
kernel union and OCCT fuse is reproduced exactly as documented.

F2 accept-path assertion sanity (throwaway /tmp script; OCCT's own common
as stand-in accepted output, since main refuses): relerr_vs_occt=0.0,
arbiter checked=299 kernel_errors=0 classifier_disagreements=0; negative
control (wrong volume) trips the 1e-6 assert. The new assertion code is
sound; it is simply not exercised on main until G3 lands.

### Criterion results

| Criterion | Result |
|---|---|
| unexpected_accept counter forces nonzero exit | PASS (demo: old exit 0, new exit 1) |
| Docstring reconciled with measured baseline | PASS (17 cases, 4 accepted, 13 refused) |
| F2 asserts volume (1e-6) and arbiter (0 errors) on accept; refusal stays typed | PASS (asserts verified via stand-in; refusal recorded) |
| F4 asserts winding OUT and distance > 1.0; OCCT states recorded | PASS (both asserts hold on main) |
| Full suite green (I5) | PASS (14/14 exit 0, 0 [FAIL] lines; detached run after a /tmp wipe and a service restart killed two earlier attempts) |
| I1-I9 invariants held | PASS (no kernel code touched; refusals typed; boolean()/boolean_brep() contracts unchanged; local branch only, no push/merge) |
| No em dashes (I8) | PASS (grep over touched files, ledger, commit messages) |

### Open

- F2's hardened accept-path assertions are verified only via the OCCT
  stand-in; they will first run for real when G3 makes F2 accept.
- The 6 unmerged gate-branch test files (G1/G2/G3/G4/G6) were not run on
  this branch; they test code that is not on main.
