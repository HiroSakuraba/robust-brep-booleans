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
grep -c $'\u2014' on touched files -> 0
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

## G6 rework - Per-face broad-phase pad (2026-09-25)

- Branch: gate/G6-broadphase-rework (local only; never pushed, never merged)
- Base: 6bd8031 "docs: record v0.9 rebase, G7 replacement, roadmap update"
  (gate/G6-broadphase-tol tip)
- Commits:
  - a63aec7 "tests: per-face broad-phase pad regression (G6 rework, I4
    failing first)"
  - d98607e "intersection: per-face broad-phase pad (G6 rework)"
- Environment: ~/workspace/brep-booleans/.venv (CPython 3.12.3,
  numpy==2.5.3, cadquery-ocp==8.0.1.0.0, OCCT 8.0.1)
- Worktree note: the shared ~/workspace/brep-gates checkout was moved to
  probe-hardening by the gate coordinator mid-task, so this work ran in a
  separate worktree ~/workspace/brep-gates-g6rework on the same branch.
- Rebase observation (per task): after the v0.9 rebase this branch no
  longer contains src/brepkernel/coincidence.py and the G2 test files are
  gone (19 test files, not 20). This does not matter for the change: the
  scope is pad granularity only (broad-phase candidate generation), which
  touches neither coincidence classification nor the G2 tests.

### Review finding addressed

G6's per-model pad, 2x the largest tolerance anywhere in either model,
is catastrophically pessimistic: one damaged 0.5 mm-tolerance edge makes
every face in a thousand-face model overlap in broad phase. Replaced with
a per-face conservative pad.

### What changed

- src/brepkernel/step_ingest.py:
  - FaceRecord gains tol_face: max OCCT tolerance over the face and its
    incident edges and vertices (BRep_Tool.Tolerance), computed in
    index_shape via the new _face_max_tolerance() helper (both ingest
    branches).
  - New face_broadphase_pads(model, contact_tol): per-face pad array,
    pad_i = contact_tol + tol_face_i.
  - candidate_face_pairs() gains optional pads_a/pads_b; _face_pairs_aabb()
    expands each face box by its own pad (scalar-pad path unchanged).
    The NURBS patch filter uses the pair sum pads_a[i] + pads_b[j].
    _face_pairs_aabb keeps its old return contract (list of pairs) after
    an interim 3-tuple broke test_freeform_nurbs's direct unit test; the
    effective pads are recomputed in candidate_face_pairs.
  - useShapeTolerance=True in _shape_bbox is untouched.
- src/brepkernel/intersection.py: intersect_models() gains optional
  broadphase_face_pads=(pads_a, pads_b), forwarded to
  candidate_face_pairs(). Scalar broadphase_pad path unchanged.
- src/brepkernel/pipeline.py: when broadphase_pad is not explicit,
  per-face pads are computed for both models and passed through; the
  report records broadphase_pad (max per-face pad),
  broadphase_pad_mode ("per_face" | "explicit_scalar"),
  broadphase_contact_tol, broadphase_max_tolerance (unchanged semantics),
  and broadphase_pad_summary {A, B}: faces/min/mean/max plus
  tolerance-driven face ids (faces whose own tolerance exceeds the
  contact band). An explicit broadphase_pad is still honored verbatim as
  a uniform scalar. crosscheck_ops companion runs receive the max
  per-face pad as their explicit scalar (conservative, as before).
- tests/test_g6_tolerance_broadphase.py: one assertion updated to the new
  formula (pad >= FACE_TOL + contact_tol instead of pad >= 2*FACE_TOL);
  the old 2x factor encoded the replaced per-model formula. Functional
  requirements unchanged: the 1e-3-tolerance pair still becomes a
  candidate and still refuses typed.

### Pad derivation (I3)

pad_i = contact_tol + max tolerance over face i and its incident edges
and vertices. Derived from entity tolerances, not tuned: a face whose
boundary entities carry tolerance t is geometrically uncertain over a
band of about t around it (OCCT tolerance semantics), and the exact
contact classifier needs a further contact_tol band so tolerance-near
contacts are not dropped as "disjoint". Additive padding only widens the
candidate set: it can add typed refusals, never new acceptances. The
pair budget for faces i, j is pad_i + pad_j = 2*contact_tol + t_i + t_j,
conservative against the minimal t_i + t_j + contact_tol.

### I4 pre-fix demonstration (test committed at a63aec7, run on 6bd8031)

```
[FAIL] g6r per-face pad API exists face_broadphase_pads importable from step_ingest

SOME FAILURES
```

### Post-fix test result (new tests)

```
[PASS] g6r per-face pad API exists face_broadphase_pads importable from step_ingest
[PASS] g6r dirty edge visible on its incident faces dirty_faces=2
[PASS] g6r dirty faces carry the damaged pad dirty_pads=['0.0005', '0.0005']
[PASS] g6r clean faces keep a small per-face pad max_clean_pad=5e-07 cap=1.4e-06 n_clean=58
[PASS] g6r old per-model pad explodes on the row old_pad=0.001 old_count=21 new_count=0
[PASS] g6r per-face candidate count stays small new_count=0

ALL PASS
```

The row model (10 unit boxes, one edge of box 0 raised to 0.5 mm) shows
the review's failure mode directly: the old per-model pad (1e-3) yields
21 candidates against a neighbor box 5e-4 away, while per-face pads
yield 0. Clean-face pads stay at 5e-7 (contact_tol 4e-7 + OCCT default
face tolerance 1e-7).

Existing G6 test, unchanged behavior:

```
[PASS] g6 pair is a candidate even with zero pad candidates=21
[PASS] g6 face tolerance visible to the kernel max_face_tol=0.001
[PASS] g6 union refuses with typed kind refusal={'stage': 'intersection',
  'type': 'IntersectionError', 'kind': 'SectionToleranceTooLoose',
  'message': 'section edge 0: OCCT edge/face tolerance 0.001 exceeds
  acceptance ceiling 1.28e-05'}
[PASS] g6 pad recorded in report and covers the tolerance broadphase_pad=0.0010004

ALL PASS
```

The G6 dirty model report now shows the per-face distribution: model A
summary faces=6 min=5e-07 mean=8.3375e-04 max=0.0010004
tolerance_driven_faces=5 (the raised +X face plus its 4 neighbors sharing
the raised edges/vertices; the -X face is clean), model B all clean.

### Candidate-count comparison (old vs new code)

Measured via boolean_brep on 4 accepting cases (script
.suite_logs/cand_counts.py, since removed); old code measured with the
implementation stashed (git stash push -- src/), then popped:

| case | old candidates | new candidates | new pad | old pad |
|---|---|---|---|---|
| A box/box union | 6 | 6 | 5e-07 | 4e-07 |
| B NURBS-sphere union | 1 | 1 | 5e-07 | 4e-07 |
| C cylinder-through-box difference | 2 | 2 | 5e-07 | 4e-07 |
| D NURBS-sphere intersection | 1 | 1 | 5e-07 | 4e-07 |

Counts change only where tolerances are large (the row model: 21 -> 0;
the G6 dirty model: 21 candidates, max pad 0.0010004). On ordinary
low-tolerance models the candidate sets are identical; the pad differs
by 1e-7 because the per-face formula now includes each face's own OCCT
tolerance (contact_tol + tol_face) instead of the bare contact_tol
floor. Explicit-pad and crosscheck_ops paths verified separately
(explicit 1e-3 honored verbatim, mode explicit_scalar; crosscheck
accepted with pad 5e-7).

### Full suite (19 test files; G2 files trimmed by the rebase)

- 17 files EXIT0, 0 [FAIL] lines in the sequential run.
- tests/test_brep_pipeline.py: the loop's run died with a 0-byte log
  during the post-reboot memory spike (environment artifact, no output at
  all); dedicated rerun with timeout 1500: ALL PASS, 0 [FAIL].
- tests/test_freeform_nurbs.py: the loop's run caught a real bug I
  introduced (interim _face_pairs_aabb 3-tuple return broke its direct
  unit test); contract restored to list-of-pairs; rerun on final code:
  ALL PASS, 0 [FAIL].
- Result: 19/19 green.

### Commands run

```
git checkout -b gate/G6-broadphase-rework gate/G6-broadphase-tol  # 6bd8031
git worktree add ~/workspace/brep-gates-g6rework gate/G6-broadphase-rework
# wrote tests/test_g6_perface_broadphase.py, committed (a63aec7), ran -> SOME FAILURES (I4)
# implemented per-face pads in step_ingest/intersection/pipeline
python tests/test_g6_perface_broadphase.py        # ALL PASS
python tests/test_g6_tolerance_broadphase.py      # ALL PASS
python .suite_logs/cand_counts.py                 # new-code counts: 6, 1, 2, 1
git stash push -q -m g6r-impl-temp -- src/        # old-code counts: 6, 1, 2, 1
git stash pop -q
# full suite: nohup loop over tests/test_*.py -> 17 EXIT0; dedicated reruns
#   for test_brep_pipeline (ALL PASS) and test_freeform_nurbs (ALL PASS)
grep for U+2014 on touched files -> 0
```

### Gate verdicts

| Criterion | Result |
|---|---|
| Per-face pad: clean faces stay small with one 0.5 mm edge | PASS (max clean pad 5e-7 <= contact_tol + 1e-6; candidate counts 21 -> 0) |
| Existing G6 test keeps passing | PASS (pair is candidate, typed SectionToleranceTooLoose refusal; pad assertion updated to the new formula) |
| Pad formula derived from entity tolerances (I3) | PASS (pad_i = contact_tol + BRep_Tool.Tolerance over face + incident edges/vertices) |
| Per-face pad recorded in the operation report | PASS (broadphase_pad max, broadphase_pad_mode, broadphase_contact_tol, broadphase_max_tolerance, broadphase_pad_summary) |
| Candidate counts change only where tolerances are large | PASS (4 accepting cases identical; row model 21 -> 0) |
| Full suite green (I5) | PASS (19/19, 0 FAIL lines) |
| I1-I9 invariants held | PASS (refusal, never wrong accept; I4 failing test first; no em dashes; local branch only, no push/merge) |

G6 rework: CLOSED.
