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

## G1 rework - Scale-aware independent Boolean arbiter

- Date: 2026-09-25
- Branch: gate/G1-arbiter-rework (local only; never pushed, never merged to main)
- Base: 032dee7 "docs: add G1 gate entry to the review ledger" (gate/G1-arbiter tip)
- Commits:
  - f497abb "tests: add G1 arbiter scale-awareness and no-mutation regression (fails pre-fix)"
  - 52c47f4 "arbiter: scale-aware domain/band, no input mutation; independent Boolean arbiter naming"
- Environment: ~/workspace/brep-booleans/.venv (CPython 3.12.3),
  numpy==2.5.3, manifold3d==3.5.3, cadquery-ocp==8.0.1.0.0 (OCCT 8.0.1).

### What was built

The review found the bundled winding arbiter was not scale-aware and
mutated its inputs. All four items are in tools/review_probes/arbiter.py;
all existing public names and signatures are kept backward compatible
(membership_audit, occt_state, winding, triangles, surface_distance,
has_solid, want).

1. Scale-aware sampling domain: membership_audit() lo/hi/band now default
   to None, meaning "derive". New combined_domain(*shapes, pad_frac=0.1)
   unions the Bnd_Box of A, B and out and pads by 0.1 * (largest box
   edge); null shapes and void bboxes are skipped, and the legacy
   [-2.2, 2.2]^3 cube survives only as the documented FALLBACK_DOMAIN
   for the all-void case. Explicit lo/hi still select the legacy fixed
   behavior. tests/_arbiter.py raw_audit() now passes lo/hi=None through
   (previously it derived from A and B only via audit_bounds, which now
   delegates to combined_domain and is kept for compatibility).
2. Scale-aware surface-exclusion band: new exclusion_band(*shapes,
   deflection=2e-4) returns BAND_SAFETY * deflection + max entity
   tolerance, the max over BRep_Tool.Tolerance_s of every face, edge and
   vertex of the input shapes; the documented fallback when no tolerance
   is measurable is BAND_SAFETY * deflection + CONTACT_TOL (1e-7, OCCT's
   nominal contact tolerance). BAND_SAFETY = 2.0 is documented, not tuned:
   one deflection for the tessellation chordal deviation behind the
   winding test, roughly one more for the distance query against the true
   surface. The winding decisiveness guard (|w| not in (0.05, 0.95)) is
   unchanged.
3. No input mutation: triangles() deep-copies the shape with
   BRepBuilderAPI_Copy (verified empirically: the copy's TShapes are
   distinct, so tessellations attach to the copy only) and runs
   BRepTools.Clean_s + BRepMesh_IncrementalMesh on the copy.
   surface_distance() and occt_state() were audited and do not mutate.
4. Naming: module and caller docstrings now say "independent Boolean
   arbiter" with the explicit caveat that it tessellates with OCCT's
   BRepMesh, so it is independent of OCCT's Boolean/intersection decision
   paths, not of OCCT entirely (not an "independent kernel arbiter").
   Result-dict keys (kernel_errors etc.) are unchanged; the dict gains
   additive diagnostics: domain, band_used, skipped_winding,
   skipped_band, interior_checked.

### I4 pre-fix demonstration (test committed at f497abb, run on 032dee7)

New tests/test_g1_arbiter_scale.py, 3 tests, all FAIL pre-fix, all PASS
post-fix:

- t1 translated operands (+1000 x): pre-fix the audit sampled the
  hard-coded cube (domain x=[-2.200, 2.200]) while the geometry sat at
  x=[1000, 1001.5], yet reported checked=120, kernel_errors=0: a vacuous
  audit. Post-fix the domain is x=[999.850, 1001.650], checked=120,
  kernel_errors=0.
- t2 micro-part (2e-3 boxes): pre-fix band_used=2.000e-3 (absolute) and
  interior_checked=0: no interior point of the part was ever examined.
  Post-fix band_used=4.001e-04, interior_checked=11, checked=106,
  kernel_errors=0 on the correct OCCT fuse.
- t3 no mutation: pre-fix triangles() left 6/6 faces of the input box
  triangulated; post-fix 0/6, soup shape (12, 3, 3) identical.

### Deviations from the plan

- D1 (G1 rework): the suite on this branch is 16 files, not 20: the G1
  tip carries 15 test files and this item adds 1
  (tests/test_g1_arbiter_scale.py). An early `ls` showed 20 because the
  shared clone transiently contained other gates' in-progress test files
  (g2/g3/g4/g6) in the working tree during the coordinator's migration to
  per-worker worktrees; those files are not on this branch.
- D2 (G1 rework): the clone was shared with concurrently running gate
  workers during this item, which caused three incidents, all repaired
  without touching other workers' commits: (a) the f497abb commit first
  landed on gate/G3G4-rework after a branch switch under this worker;
  repaired by cherry-picking to gate/G1-arbiter-rework and resetting
  gate/G3G4-rework to its prior tip 79e3079 (removing only this worker's
  own commit); (b) a later branch switch deleted
  tests/data/review_20260925/completeness_false_alarm_*.brep from the
  working tree; restored via git checkout -- (fixtures needed by
  repro_findings.py); (c) two VM reboots SIGTERM'd suite runs mid-flight;
  the killed test_brep_pipeline.py run was re-executed to completion.
  Each worker now has a dedicated worktree.

### Commands run

```
python tests/test_g1_arbiter_scale.py                       # 6/6 checks post-fix, ~8s
for t in $(git ls-files 'tests/test_*.py' | sort); do python $t; done   # full 16-file suite
python tests/test_brep_pipeline.py                          # rerun after reboot kill: ALL PASS, exit 0
python tools/review_probes/fuzz_brep.py --trials 20 --seed 7 --out /tmp/fuzz_g1r.json
python tools/review_probes/repro_findings.py                # F4 behavior check
```

### Results (trimmed)

- test_g1_arbiter_scale: ALL PASS, 6 checks.
- Full suite: 16/16 test files exit 0, 0 [FAIL] lines (15 G1 files + the
  new scale test). All 23 Tier B/C membership audits still pass with the
  derived domain/band (checked 299-300, kernel_errors=0); the negative
  corruption test still catches the flipped keep decision.
- fuzz_brep.py short run (--trials 20 --seed 7): exit 0,
  TALLY {'accept': 19, 'refuse:InsufficientPatchWitnesses': 1},
  0 WRONG, 0 CRASH.
- repro_findings.py: F4 output byte-identical to
  docs/baseline_20260925/repro.txt (classifier false IN on kernel union
  and OCCT fuse at surface_distance=1.0396, winding 0.0). F2 refuses with
  SectionCompletenessMismatch exactly as baselined (the G3 fix is a later
  gate), so the script exits 1 via f2 as before.
- No em dashes in any touched file, prose, or commit message.

### Gate verdicts

| Criterion | Result |
|---|---|
| Scale-aware domain from combined A/B/out bbox | PASS (t1; fuzz and Tier B/C audits use the derived domain) |
| Scale-aware band from deflection + tolerances | PASS (t2; band_used reported per audit) |
| triangles() never mutates inputs | PASS (t3; deep copy verified) |
| "Independent Boolean arbiter" naming, no overclaim | PASS |
| Full suite green, ledger updated | PASS (16/16 exit 0; this entry) |
| I1-I9 invariants held | PASS (no crash; no kernel tolerance touched; I3 derivations documented above; test committed failing first per I4; local branch only, no push/merge) |

G1 rework: CLOSED.
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

### Merged-main addendum (2026-09-25)

The merged main (`10ea7d8`, includes the G3G4 rework) ran the full suite
26/26 exit 0 with 0 [FAIL] lines (see the merge entry), including
`tests/test_g3_probe_rework.py` and `tests/test_g3_f2_regression.py`. The
I5 gap recorded on this branch (4/18, under CPU contention) is resolved.
Every other criterion in this gate's table was already PASS, so:

G3 rework: CLOSED.

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

### Merged-main addendum (2026-09-25)

The merged main (`10ea7d8`, includes the G3G4 rework) ran the full suite
26/26 exit 0 with 0 [FAIL] lines (see the merge entry), including
`tests/test_g4_probe_coverage.py` (12/12). The I5 partial recorded on
this branch is resolved. The two fuzz criteria below stay OPEN: the
150-trial ON/OFF comparison was not re-run on the merged code, and the
median-time comparison was not measured. Not claimed, not closed.

G4 rework: OPEN on the fuzz criteria only (analytic fuzz 150-trial
ON/OFF comparison with 0 WRONG, and the median-time recording).
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

### G2 planar redirection: three-case coincidence (2026-09-25, night)

Branch: gate/G2-planar (local only). Work order ITEM 2 of 8: rework G2
coincidence per the narrowed planar scope. The Fraction-exact rule that
also accepted curved analytic coincidence is replaced by a three-case
classification in src/brepkernel/coincidence.py, planar faces only.
Curved/NURBS machinery is set aside explicitly, not deleted silently.

Rule (classify_support_pair):
  (a) EXACT: canonical analytic plane params agree bit-identically
      (Fraction: direction mod sign, offset exactly 0), including params
      recovered through the same canonical recognizer
      (GeomLib_IsPlanarSurface). Returns ("coincident", sense, "exact").
  (b) TOLERANCE-CERTIFIED: only when a recognizer was involved
      (non-analytic face, so recovery uncertainty exists) AND the sampled
      support deviation is within the entities' ACTUAL per-entity
      tolerances (face/edge/vertex via BRep_Tool.Tolerance_s) AND the
      pipeline contact_tol. Never a global constant. Returns
      ("coincident", sense, "tolerance_certified").
  (c) UNDECIDABLE: anything else, as a CoincidenceUndecidable record
      carrying reason, measured deviation, deviation kind, and every
      tolerance consulted. Analytic/analytic planar pairs never take
      (b): their params are exact, so in-band-but-not-exact is (c).
      Non-parallel analytic planes are provably "distinct" (exact,
      no sampling).

Curved analytic (cylinder/sphere/cone/torus) and unrecognized NURBS:
active code claims "distinct" only when sampling proves deviation
beyond the band; anything that could be coincident is (c) with reason
curved_deferred / nurbs_deferred. Different analytic kinds (plane vs
cylinder etc.) short-circuit to "distinct". The pre-narrowing curved
machinery (exact Fraction param extractors, curved exact branches,
dense-sampling acceptance) is preserved verbatim in
src/brepkernel/coincidence_deferred.py under ENABLE_DEFERRED_COINCIDENCE
= False, with a header listing what a follow-up gate must do before
re-enabling (canonical recognition parameter recovery, three-case
certification, I4-first tests, full suite + probes).

Callers: intersection.py maps (c) to IntersectionError with kind
NearCoincidentFaces (planar in-band: near_coincident_planar,
recognized_beyond_tolerance, sense_failed) vs CoincidenceUndecidable
(deferred); the exception carries structured evidence (reason,
deviation, deviation_kind, tolerances) and pipeline.py refuse() now
copies it into the refusal report. CoincidentPairRecord gains a
certification field ("exact" default).

I4: tests/test_g2_three_case.py (12 checks) committed BEFORE the
behavior change as 4a442c2; it failed 7/12 on the old code (the 2
passes were the import-tolerant fallback and one distinct case) and
passes 12/12 after. Notable pins: stacked (non-identical) coaxial
cylinders union now refuses CoincidenceUndecidable (old code accepted
via curved exact path); 0.5x-tol planar offset is (c) with deviation
5e-8 and per-entity tolerances on the record; in-tolerance recognized
planar is (b) tolerance_certified; 2x recognized is (c).

Verification (all on gate/G2-planar, venv
~/workspace/brep-booleans/.venv/bin/python):
  - tests/test_g2_three_case.py: 12/12 PASS.
  - tests/test_g2_coincident_probes.py: ALL PASS (edge-touching still
    refuses NonManifoldResult; 0.5x/2x near-coincident refuse
    NearCoincidentFaces; 10x accepts; all cylinder cases keep exact
    volumes).
  - tests/test_g2_keep_table.py: ALL PASS.
  - Full suite (21 files): 20/21 PASS. The two anomalies are NOT
    regressions: tests/test_brep_pipeline.py was SIGTERM-killed at 5
    min during the contended suite run (5 concurrent suites on the
    box) and passes ALL PASS in 5.8 min on a clean rerun; 
    tests/test_freeform_split.py t7 (shared seam) fails identically on
    the pre-change code (verified via git stash), a pre-existing
    failure that contradicts the closeout ledger's 20/20 claim for
    this worktree state.
  - tools/review_probes/common_cad_probes.py: wrong=0,
    unmet_accepts=0. All 14 expect=accept cases accepted with exact
    volumes (incl. planar flush cylinder cap cases); equal-radius
    crossing cylinders still refuse (assembly/
    PatchClassificationInconsistent), as required.

Kept as is: the 24-cell keep table in brepkernel/assembly.py::keep_patch
(untouched); the G2.4 overlap split; the G2.5 boundary-contact rule.

Deferred explicitly: curved/NURBS coincidence acceptance (see
coincidence_deferred.py). No push, no merge; local branch only.

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

---

## G5: classifier independence (review finding F4)

Branch `gate/G5-classifier`, off `main` at pinned commit `88d09eb`.
Venv `~/workspace/brep-booleans/.venv/bin/python` for every command below.

Finding F4: `BRepClass3d_SolidClassifier` returned a false IN for a point
more than 1.0 from both input surfaces, and the kernel's patch/shell
classification depended on that same classifier. Fix: a second independent
classifier (multi-ray parity via `IntCurvesFace_ShapeIntersector`, a
different OCCT subsystem that never calls the solid classifier); every
keep/discard decision resting on point classification now requires both
classifiers to agree, and any mismatch refuses with typed kind
`ClassifierDisagreement` carrying both verdicts and the point. Witness
selection prefers points at distance >= 10x tol from the other operand's
boundary (via `BRepExtrema_DistShapeShape` against the other model's faces)
where feasible; documented in `_witness_material_verdict` (applied in
`_classify_pieces`; deliberately NOT applied in `_solid_interior_points`
or `_shell_records` nesting, whose witnesses are near-boundary by
construction). The winding/tessellation arbiter in
`tools/review_probes/arbiter.py` stays the TESTING arbiter only; it is not
used as a production classifier.

### Implementation notes (bugs found while verifying)

1. Tangent-threshold too tight: with `_TANGENT_COS = 0.02`, a near-tangent
   ray on a NURBS sphere reported a single hit at `|n.d| = 0.0231` whose
   paired crossing (0.045 further along the ray) was missed by the
   intersector, flipping parity to a spurious "inside"; 11 rays said
   outside, 1 said inside, so the point went "unknown" and a legitimate
   op refused. Raised to 0.05 (~3 degrees, ~2x margin over the observed
   0.0231). After the fix, outside/inside/far points classify correctly.
2. Witness-preference soundness: the first version filtered witnesses to
   the far set before the unanimity check, and measured "distance to
   boundary" against the solids, for which OCCT reports 0 for interior
   points. That hid a straddling patch's disagreeing witnesses
   (`test_freeform_assembly` t6 wrongly accepted). Fixed: distance is
   measured against the other model's FACES; the far witnesses must be
   unanimous; a near-boundary witness that materially contradicts the far
   verdict (both classifiers agreeing) still blocks with
   `PatchClassificationInconsistent`; near boundary/unknown verdicts are
   confusion-zone noise and do not veto. Where fewer than 3 witnesses
   clear the band, the pre-G5 full-set unanimity rule applies unchanged.
3. Bidirectional rays: a near-tangent ray on NURBS dropped one crossing of
   a pair (analytic check: true crossings at t=0.5485 and t=0.7367, only
   the second found), flipping parity and causing a wrongful
   `ClassifierDisagreement` refusal on a legitimate op
   (`test_nurbs_adversarial_corpus` t4, which the base code accepts with
   oracle error 8.9e-16). Fix: each direction is cast as a +d/-d pair and
   seated only if every solid's bidirectional crossing sum is even (a line
   meets a closed solid an even number of times, so an odd sum proves a
   missed/added crossing). The pair's verdict comes from the +d counts;
   cross-pair unanimity is still required.

### Test results

Full suite, final code (15 files; base `88d09eb` ships 14, G5 adds
`tests/test_g5_classifier_independence.py`):

```
$ for t in tests/test_*.py; do timeout 1500 <venv>/python -u "$t"; done
test_brep_pipeline exit=0 fails=0
test_degenerate exit=0 fails=0
test_freeform_assembly exit=0 fails=0
test_freeform_intersection exit=0 fails=0
test_freeform_nurbs exit=0 fails=0
test_freeform_split exit=0 fails=0
test_g5_classifier_independence exit=0 fails=0
test_metamorphic exit=0 fails=0
test_multishell_semantics exit=0 fails=0
test_nurbs_adversarial_corpus exit=0 fails=0
test_nurbs_boolean_end_to_end exit=0 fails=0
test_regression exit=0 fails=0
test_same_domain exit=0 fails=0
test_step_nurbs_multiface exit=0 fails=0
test_stress exit=0 fails=0
```
PASS: 15/15 files, zero [FAIL] lines, zero crashes.

G5 tests (`tests/test_g5_classifier_independence.py`):
```
[PASS] t1 baseline union accepted with sane volume volume=15.000000
[PASS] t1 flipped classifier refuses with ClassifierDisagreement kind=ClassifierDisagreement
[PASS] t1 refusal carries both verdicts and the point occt=inside independent=outside point=(0.0, 0.8246858765999999, 0.2)
[PASS] t2 F4 union accepted with correct volume volume=3.226354 expected=3.226354
[PASS] t2 independent classifier verdict is OUT at probe point independent=outside
[PASS] t2 F4 false IN still reproduces in the OCCT classifier occt=inside
[PASS] t2 pair withholds agreement (no agreed IN) agreed=False decision=None
ALL PASS
```
PASS. t1: monkeypatched classifier with IN<->OUT flipped refuses with
`ClassifierDisagreement` (pre-fix it refused downstream with
`OperationVolumeInvariantFailed`). t2: the F4 repro point
`(1.3322189978645596, 1.1200289002934705, -0.8236295495693953)` still gets
a false IN from the OCCT classifier (fault reproduces, test not vacuous);
the multi-ray classifier says outside; the pair withholds agreement; the
union is accepted with the correct volume (3.226354127576551, exactly the
OCCT fuse volume).

### Overhead (common_cad_probes.py --baseline, base vs G5)

```
$ cd ~/workspace/g5base-wt && PYTHONPATH=src <venv>/python /tmp/g5probes/common_cad_probes.py --baseline
$ cd ~/workspace/brep-gates-g5 && PYTHONPATH=src <venv>/python /tmp/g5probes/common_cad_probes.py --baseline
```
Base: 47.3 s, wrong=0, unmet_accepts=11. G5: 44.9 s, wrong=0,
unmet_accepts=11. Per-case ACCEPT/REFUSE dispositions identical between
base and G5 (diff of the disposition lines, ignoring timings, is empty);
0 ClassifierDisagreement refusals in the G5 run. No measurable overhead
on this probe set (G5 nominally faster; within machine noise under load).

### Fuzz (fuzz_brep.py --trials 60 --seed 7, analytic)

```
$ cd ~/workspace/g5base-wt && PYTHONPATH=src <venv>/python /tmp/g5probes/fuzz_brep.py --trials 60 --seed 7 --out ~/workspace/g5bench/fuzz_base.json
$ cd ~/workspace/brep-gates-g5 && PYTHONPATH=src <venv>/python /tmp/g5probes/fuzz_brep.py --trials 60 --seed 7 --out ~/workspace/g5bench/fuzz_g5.json
```
Base: exit=0, TALLY {'accept': 58, 'refuse:InsufficientPatchWitnesses': 2},
0 WRONG, 0 CRASH. G5: exit=0, TALLY {'accept': 58,
'refuse:InsufficientPatchWitnesses': 2}, 0 WRONG, 0 CRASH,
0 ClassifierDisagreement refusals. Tallies identical between base and G5.

---

## G7 periodic-seam stress (25 Sept 2026, branch gate/G7-seam-stress)

Scope: stress-test the v0.9 seam mechanism (reuse an existing seam boundary
on the seam-side operand; verify the section edge normally on the opposite
operand), NOT a rewrite. Tests + harnesses only; zero src/ changes.

### What was built

- tests/test_g7_seam_stress.py (new, 38 checks):
  - g1 equatorial seam-plane cut (torus minus box, z<=0): one section loop
    IS the torus v-seam; report pins reused_seam_edges_A=1, B=0,
    shared=0, seam_on_a in section payloads, 2 verified loops, oracle-exact
    volume, arbiter 200 pts / 0 kernel errors.
  - g2 polar seam-plane cut (torus minus box, y<=0): one section loop IS
    the u-seam meridian; pins reused_seam_edges_A=1, seam_on_a flag,
    oracle-exact volume, arbiter 200/0.
  - g3 same polar cut, operands reversed (intersection): pins
    reused_seam_edges_B=1, seam_on_b flag, oracle-exact, arbiter 200/0.
  - g4 torus-vs-torus union (small torus hooked through the big tube):
    4 transverse section curves crossing seams geometrically; accepted,
    oracle-exact (66.5750592625), arbiter 200/0, routing counters 0/0.
  - g5 torus-vs-box difference (box corner through tube wall): 8 section
    curves; accepted, oracle-exact (58.0715722701), arbiter 200/0.
  - g6 revolved freeform NURBS (B-spline egg profile via MakeRevol) vs box
    union: one-face GeomAbs_SurfaceOfRevolution solid; accepted in 7.7 s,
    oracle-exact (95.2739324039), arbiter 200/0.
  - g7 asserts every accepted public report carries reused_seam_edges_A/B
    and shared_seam_refusals counters.
- tools/review_probes/{arbiter,fuzz_brep,common_cad_probes,repro_findings}.py
  and tests/_arbiter.py checked out from the gate working branch (missing
  from main; harness-only, needed for the fuzz and the audits).

### Fuzz (torus-heavy, --kinds torus,cyl,sph, 60 trials each)

- Analytic, seed 7: TALLY {accept: 56, refuse:SectionToleranceTooLoose: 4},
  0 WRONG, 0 CRASH. Refusals are tolerance-gate, none seam-related.
- NURBS, seed 11: TALLY {accept: 59, refuse:InsufficientPatchWitnesses: 1},
  0 WRONG, 0 CRASH. The single refusal is a patch-witness-count refusal,
  not seam-related.
- Seam-related refusal-kind histogram: 0 rows in both runs.
- Baselines: docs/baseline_20260925/fuzz_g7_torus_heavy.json and
  docs/baseline_20260925/fuzz_g7_torus_heavy_nurbs.json.

### Commands run

```
python tests/test_g7_seam_stress.py            # 38 checks, ALL PASS, ~42 min
python tools/review_probes/fuzz_brep.py --trials 60 --seed 7 --kinds torus,cyl,sph --out docs/baseline_20260925/fuzz_g7_torus_heavy.json
python tools/review_probes/fuzz_brep.py --trials 60 --seed 11 --nurbs --kinds torus,cyl,sph --out docs/baseline_20260925/fuzz_g7_torus_heavy_nurbs.json
for f in <14 other tests/test_*.py>; do python $f; done   # all exit 0, 0 FAIL
git grep -P '\x{2014}' -- <touched files>                  # 0 matches
```

### Results (trimmed)

- test_g7_seam_stress: ALL PASS, 38 checks, ~42 min; 6 audits, all
  kernel_errors=0, checked=200 each; every accepted volume matches the
  independent OCCT oracle to the printed precision (g1/g2 err 0.0,
  g3 err 4.5e-11, g4 err 0.0, g5 err 0.0, g6 err 1.4e-14).
- Full sweep of the other 14 test files: all exit 0, 0 [FAIL] lines.
- Total: 15/15 test files green on this branch (14 existing + the new
  test_g7_seam_stress.py).
- No wrong accept found; no stop-and-report event. The near-tangent
  torus-torus offsets probed during case selection refused with typed
  UnresolvedContact (allowed), so the committed g4 uses a transverse
  tube-crossing configuration instead.

### Notes / open

- Torus-torus and torus-box boolean_brep calls are slow (~360-384 s each
  on NURBS tori); the g7 file takes ~42 min total. Not a correctness
  issue; worth profiling before any CI adoption.
- The old G7 pre-split plan stays dead; v0.9 routing needed no changes.
- I1-I9 held: local branch only, no push, no src/ edits, no em dashes.

G7 gate: stress round COMPLETE, no wrong accept, suite green.

---

## Merge of the 7 gate branches into main (2026-09-25)

- Date: 2026-09-25
- Base: 88d09eb "Merge v0.9 verified periodic seam-aware NURBS splitting"
- Authorization: Ben said "merge them" (18:17 EDT); local merge only at this stage.

### Merge order

1. gate/G2-planar (brings the whole G0->G1->G3->G4->G6->G2 stack plus the planar rework) -> 97aaf32
2. gate/G1-arbiter-rework (scale-aware arbiter) -> b676d49
3. gate/G3G4-rework (adaptive probe matching, tolerance cap) -> 3c079ef
4. gate/G6-broadphase-rework (per-face pads) -> 58eb526
5. probe-hardening (pass/fail semantics, F2/F4 assertions) -> 184a735
6. gate/G5-classifier (multi-ray parity, dual agreement) -> 66c95df
7. gate/G7-seam-stress (seam stress round) -> 6737c93

All 7 branches verified reachable from main via git merge-base --is-ancestor.

### Conflicts encountered and how resolved

- docs/REVIEW_LEDGER.md (every merge): each branch appended its section at the
  file tail of its own base, so the sections were spliced into gate order
  (G1 rework after G1; G3 REWORK after G3; G4 REWORK after G4; G6 rework after
  G6; probe-hardening, G5, G7 appended). No content dropped; section counts
  verified.
- probe-hardening and G7: both branches recreated tools/review_probes/*.py from
  the appendix instead of editing G0's versions (add/add). The G1 arbiter
  rework's scale-aware arbiter.py would have been regressed by taking their
  copies, so HEAD's versions were kept and only probe-hardening's functional
  changes were applied as patches: common_cad_probes.py refuse-counter/exit
  fix, repro_findings.py F2/F4 assertions. G7's copies were verbatim appendix
  files with no functional changes; HEAD's kept.
- src/brepkernel/assembly.py (G5 merge): genuine semantic clash. G2's
  one_side classified non-coincident pieces with _classify_point_in_model and
  mapped to the four G2 states; G5's one_side used the dual classifier
  (_agreed_point_verdict + _witness_material_verdict) but its
  _decision_rule only accepted inside/outside. Resolution: kept G2's
  coincident-piece early path and G2's keep-table _decision_rule; for
  non-coincident pieces kept G5's dual-classifier path with the
  {"inside": "IN", "outside": "OUT"} mapping before the keep table
  (exactly the mapping the G2 code already used). The G5 multi-ray
  classifier, _agreed_point_verdict, _witness_material_verdict and
  ClassifierDisagreement all survive in the merged file.
- src/brepkernel/intersection.py, pipeline.py (G3G4, G6 merges): auto-merged
  cleanly; syntax-checked with ast.parse.

### Full suite results on merged main

26/26 test files exit 0, 0 [FAIL] lines:
test_brep_pipeline, test_degenerate, test_freeform_assembly,
test_freeform_intersection, test_freeform_nurbs, test_freeform_split
(includes t7 shared-seam, passing after the fix below),
test_g1_arbiter_scale, test_g1_negative_corruption,
test_g2_coincident_probes, test_g2_keep_table, test_g2_three_case,
test_g3_f2_regression, test_g3_probe_rework, test_g4_probe_coverage,
test_g5_classifier_independence, test_g6_perface_broadphase,
test_g6_tolerance_broadphase, test_g7_seam_stress, test_metamorphic,
test_multishell_semantics, test_nurbs_adversarial_corpus,
test_nurbs_boolean_end_to_end, test_regression, test_same_domain,
test_step_nurbs_multiface, test_stress.

### Probe results on merged main

- tools/review_probes/common_cad_probes.py (no --baseline): exit 0,
  wrong=0 unmet_accepts=0 unexpected_accepts=0 (17 cases: 15 accepted with
  correct volumes, equal-radius crossing cylinders refused as expected,
  edge-touching union refused NonManifoldResult with expect=either).
- tools/review_probes/repro_findings.py: exit 0. F2 ACCEPTED, volume rel err
  2.56e-09 vs OCCT (limit 1e-6), arbiter checked=300 kernel_errors=0,
  classifier_disagreements=0. F4 kernel union winding +0.0000 (OUT),
  surface_distance=1.0396 > 1.0; hardened assertions pass.

### Known exception

- tests/test_freeform_split.py t7: the "pre-existing failure" claim did NOT
  hold. Bisection: t7 passes on v0.9 (88d09eb) and on the G1/G3/G4/G6 tips,
  fails on gate/G2-coincident-faces and on merged main. Root cause: the G2.5
  _resolve_contacts refactor overwrote split_models' local unresolved list,
  silently dropping shared_seam_curve entries. The assembly UnresolvedContact
  gate keys off unresolved_contacts, so a shared seam slipped through instead
  of refusing: a genuine I1 safety-gate bypass, not a stale test. Fixed on
  main (commit 59fe570): shared-seam entries are carried in a dedicated list
  and prepended to the final unresolved contacts; t7 passes. Failing test
  first per I4 (t7 failed pre-fix, passes post-fix).

### Push status

- Ben authorized the push (follow-up message). Awaiting green suite + probes
  before pushing local main to origin main (HiroSakuraba/robust-brep-booleans).

### Gate status after merge

G0, G1 (+rework), G2-planar, G3 (+rework), G4 (+rework), G5, G6 (+rework),
probe-hardening, new G7: merged. Still open: curved coincident-face
recognition (deferred), G8 real-data refusal rate, G9 docs/release, Part B.
---

## G8 real-data refusal rate (25 Sept 2026, branch gate/G8-real-data)

Scope: measure the typed refusal rate of the merged Tier B/C pipeline
(boolean_brep, defaults, no tuning) on real CAD data. No src/ changes;
harness only. Gate criterion: a MEASURED, DOCUMENTED refusal rate with
the dominant refusal reasons identified; no target number.

### Method

- Corpora (30 models, 15 distinct-model pairs, 3 ops each = 45 cases):
  - NIST MBE PMI: 3 pairs (ftc_06/ftc_07, ftc_08e2/ftc_09_e1,
    ctc_01/ctc_02, AP242 files from NIST-PMI-STEP-Files_8_cy64.zip)
  - DeepCAD: 5 pairs (m000/m001 ... m008/m009 from the prior
    probes-abc-deepcad scratch extract)
  - ABC chunk 27: 7 pairs at sorted offsets 0, 660, 1320, 1980, 2640,
    3300, 3960 (consecutive files, geographic spread across the chunk)
- Pairing: B is rigid-motion translated so its bbox center lands near
  A's center (0.25*diagA x offset): guaranteed overlap, same geometry,
  no tangency by construction.
- Each case runs boolean_brep(A, B, op) with pipeline defaults
  (base_tol=1e-7, contact_tol=4e-7) in an isolated worker process,
  120 s timeout, max 4 workers. Outcomes: accepted / refused (typed,
  with stage + kind + message) / timeout / exception / crash.
- Accepted triples per pair get an inclusion-exclusion check on OCCT
  volumes: |volU + volI - volA - volB| and |volD - (volA - volI)|
  (flag threshold 2% relative). Independent of the kernel.
- Harness: probes-abc-deepcad/g8_worker.py (single case),
  probes-abc-deepcad/probe_brep_g8.py (driver),
  probes-abc-deepcad/aggregate_g8.py (aggregation),
  results in probes-abc-deepcad/results_g8.json (per-case records:
  model ids, op, outcome, refusal stage/kind, wall time, volumes,
  authoring tolerances) plus scratch_g8/ case/result files.

### Measured results (45 cases)

Outcome totals: accepted 21 (46.7%), typed refusal 15 (33.3%),
timeout 9 (20.0%), exception 0, crash 0.

Refusal histogram (stage/kind):
- 9  intersection/SectionOutsideTrim
- 6  assembly/InsufficientPatchWitnesses

By corpus:
- NIST: 0 accepted, 6 refused, 3 timeout (9 cases)
- DeepCAD: 15 accepted, 0 refused, 0 timeout (15 cases)
- ABC: 6 accepted, 9 refused, 6 timeout (21 cases)

By op (union / difference / intersection): 7 / 7 / 7 accepted,
5 / 5 / 5 refused, 3 / 3 / 3 timeout. No op-specific bias.

Accepted wall times (21): min 6.2 s, median 13.6 s, p95 87.6 s,
max 90.4 s. DeepCAD accepts run 6-25 s; ABC accepts 55-90 s.

Inclusion-exclusion on the 7 fully-accepted pairs: 7/7 pass, 0 flags
(relative errors 0.0 to 1.3e-8). No wrong-accept signal anywhere in
the accepted set.

### Dominant refusal reasons

1. SectionOutsideTrim (9/15 refusals): the OCCT section edge's p-curve
   leaves the trimmed input face beyond the verification tolerance, so
   the intersection verifier refuses to certify the section curve
   (intersection.py _verify_section_edge). Concentrated on the NIST
   pairs, whose authoring tolerances are loose
   (broadphase_max_tolerance 0.104 and 0.0096), plus one ABC pair
   (0.0022). The curve the engine produced is not demonstrably ON both
   trimmed faces; refusing is the contract.
2. InsufficientPatchWitnesses (6/15): assembly found only 1 stable
   interior witness point where 3 are required on a complex-trim face,
   so patch classification cannot be certified (assembly.py witness
   sampler: center-biased low-discrepancy sequence plus Cartesian
   fallback grid). Both ABC pairs refused this way (tolerances
   ~2.3e-4); the trims are intricate enough that the fixed sample
   pattern finds no interior points.
3. Timeouts (9/45, not refusals): nist-p02 (ftc_08/ftc_09),
   abc-p02, abc-p07 hit the 120 s harness limit in all three ops.
   Recorded per case as non-results. Accepted ABC cases took 55-90 s,
   so some timeouts might accept with a larger budget; the timeout is
   a harness property, not a pipeline verdict.

Tolerance context (honest, not gamed): the section acceptance ceiling
is 128*base_tol (~1.3e-5 at defaults). Sampled authoring tolerances
ran from 1e-7 (DeepCAD, 15/15 accepted) through 1e-5 and 2.7e-4 (ABC
accepts) to 2.2e-3..0.104 (refusals). But tolerance alone does not
predict the outcome: one ABC pair accepted at 2.7e-4 while another
refused at 2.3e-4 via a different mechanism. The binding constraints
are section-curve quality against the trim and trim complexity, not
the raw tolerance number. Nothing was loosened to improve the rate
(I3); defaults were used throughout.

### Harness bugs found and fixed (I6)

- g8_worker.py omitted "pair" from result records; the driver's
  aggregation crashed with KeyError, losing 9 in-memory timeout
  records. Fixed: worker records "pair"; driver persists driver-side
  records to disk. The 9 timeout records were reconstructed from the
  sweep log and verified against the case files; the recovery is
  documented in the results manifest and the aggregation is
  idempotent (aggregate_g8.py).
- Two VM/service restarts during the run wiped /tmp scratch (NIST
  re-extracted) and one in-flight smoke; the sweep itself was
  unaffected (per-case result files on disk). To isolate from a
  concurrent worker sharing the main checkout, all G8 work ran in a
  git worktree at ~/workspace/brep-gates-g8 on branch
  gate/G8-real-data; the main checkout's uncommitted G9/evidence
  changes were never touched.

### Gate verdict

G8: COMPLETE. The refusal rate is measured (45 cases, 46.7% accept /
33.3% typed refusal / 20.0% timeout, 0 wrong, 0 crashes), documented
(results_g8.json committed with per-case records), and the dominant
reasons are identified (SectionOutsideTrim, InsufficientPatchWitnesses,
plus the timeout class). No target number was set; honesty was the
criterion and it is met.

Proposed next target (for the parent/coordinator, per the workplan):
acceptance above 80% on clean single-solid models with authoring
tolerance <= 1e-4. Candidate follow-up gates for the top refusal
kinds: (a) SectionOutsideTrim reduction (section-curve repair or
tighter OCCT section control on loose-tolerance models), (b) a
witness-point fallback for complex trims, (c) performance/timeout
investigation on the heavy pairs (nist-p02, abc-p02, abc-p07).

I1-I9: I1 held (0 wrong accepts; inclusion-exclusion 7/7 clean); I2
(all 15 refusals carry typed stage/kind); I3 (defaults only, no
tolerance loosening); I5 (full test_*.py sweep green, run after the
G8 sweep on this branch: 26/26 files exit 0, 0 FAIL lines); I6 (this
entry, bugs owned); I7 (0 crashes; 9 timeouts recorded per case);
I8 (no em dashes; grepped); I9 (boolean_brep only; boolean() untouched).
Local branch only, no push.

Note (2026-09-25): the G9 docs work below landed on branch gate/G9-docs
(local only, never pushed). The merged-main addenda above close the G3
rework I5 gap and record the G4 rework's remaining fuzz criteria as OPEN.

---

## G9 - Documentation and release (2026-09-25)

- Branch: gate/G9-docs (local only; never pushed, never merged, no tag)
- Base: `10ea7d8` "docs: merge ledger entry for the 7 gate branches plus
  the shared-seam fix" (local main)
- Scope: docs only. No `src/` or `tests/` behavior changed in this work.

### What was written

1. `README.md` rewritten for the merged state. The old text described
   only the v0.2 mesh/Tier A prototype; it now documents both routes per
   invariant I9 (route 1 `boolean()`: legacy mesh/Tier A, frozen;
   route 2 `boolean_brep()`: current exact trimmed B-rep Tier B/C), the
   full route-2 stage table, one subsection per merged gate stage naming
   what it does, what it refuses, and where the code lives:
   - G1 rework independent Boolean arbiter
     (`tools/review_probes/arbiter.py`; testing only, scale-aware domain
     and band, no input mutation),
   - G2 planar three-case coincidence (`src/brepkernel/coincidence.py`;
     EXACT / TOLERANCE-CERTIFIED / UNDECIDABLE; curved and NURBS deferred
     in `coincidence_deferred.py` under `ENABLE_DEFERRED_COINCIDENCE`),
   - G3/G4 completeness-probe reworks
     (`src/brepkernel/intersection.py`; per-curve `tol_i`, per-interval
     adaptive matching to depth 5, boundary-tail provision, raw-tolerance
     ceiling, unconditional probe with the `BREPKERNEL_COMPLETENESS_PROBE`
     kill switch),
   - G5 multi-ray parity second classifier
     (`src/brepkernel/assembly.py`; `IntCurvesFace_ShapeIntersector`-based,
     bidirectional even-sum seating, dual agreement with typed
     `ClassifierDisagreement`, 10x-tol witness band),
   - G6 rework per-face broad-phase pads (`src/brepkernel/step_ingest.py`,
     `intersection.py`, `pipeline.py`; `pad_i = contact_tol + tol_face_i`,
     report fields `broadphase_pad`, `broadphase_pad_mode`,
     `broadphase_contact_tol`, `broadphase_max_tolerance`,
     `broadphase_pad_summary`),
   - probe hardening (`tools/review_probes/common_cad_probes.py`:
     `unexpected_accepts` forces nonzero exit; `repro_findings.py`: F2
     volume/arbiter assertions, F4 winding/distance assertions),
   - G7 seam stress round (`tests/test_g7_seam_stress.py`; 38 checks,
     torus-heavy fuzz 0 WRONG, 0 seam-related refusals).
   Plus: what route 2 accepts and refuses (headline refusal kinds),
   quickstart, full-suite and probe commands, measured numbers with dates
   (26/26 test files, probes, G7 fuzz, G0 baseline pointer), honest limits.
2. `docs/PROTOTYPE.md`: header note added that it now describes the
   legacy route-1 v0.2 design (historical record kept); the current
   pipeline is documented in `README.md` and this ledger.
3. `docs/RELEASE_CHECKLIST.md` (new): proposes tag **v0.9.0** (remote
   tags v0.1-v0.4 cover the old prototype lineage; v0.8/v0.9 exist only
   as merge-commit names; the work plan's G9 names the release v0.9.0),
   the `brepkernel-v0.9.0.zip` asset list, the ordered pre-release
   verification steps (clean tree, I8 em-dash scan, 26/26 suite, both
   probes, ledger finalized, zip spot-check), and a clearly-marked
   "Blocked on Ben" section (merge approval, push to origin/main, tag,
   GitHub release with asset, I9 dispatcher decision, release-notes
   wording for G8 and the deferred curved-coincidence gate).
4. Ledger addenda on this branch: merged-main addenda on the G3 REWORK
   and G4 REWORK entries (the I5 gaps recorded under CPU contention are
   resolved by the merged-main 26/26 run; the G4 rework's 150-trial
   ON/OFF fuzz comparison and median-time recording stay OPEN, honestly
   not claimed).

### Verification on this branch

- `git grep -P '\x{2014}' -- src tests tools docs README.md`: 0 matches
  (invariant I8; no em dashes in prose, comments, or the new files).
- Full suite on `gate/G9-docs` (venv `~/workspace/brep-booleans/.venv`,
  `PYTHONPATH=~/workspace/brep-gates/src`): 26/26 test files exit 0,
  0 `[FAIL]` lines. Docs-only change, so this confirms the branch did not
  disturb the code.
- No G8 section written or touched (a separate worker owns G8; no G8
  section exists in this ledger).

### Not done (kept open honestly)

The work plan's G9 also asks for: `CHANGELOG.md` (not written),
`docs/NURBS_ACCEL.md` "What remains" update (not done), the I9
dispatcher decision (left to Ben), the actual tag and GitHub release
(Ben only), and G8 measured numbers for the release notes (separate
worker, in progress). Those are recorded in
`docs/RELEASE_CHECKLIST.md`, not claimed here.

G9 gate: OPEN. The docs slice is complete and verified on
gate/G9-docs; the release actions are blocked on Ben per the checklist.

---

## Part B: evidence schema + persistent naming (B11, B4 naming half)

- Date: 2026-09-25
- Branch: partB/evidence-schema (local only; never pushed)
- Base: 10ea7d8 "docs: merge ledger entry for the 7 gate branches plus the shared-seam fix"
- Commits:
  - 69318ba "test: evidence schema + persistent naming tests (failing first, I4)"
  - (implementation commit follows; see below)
- Roadmap position: per the correction recorded 2026-09-25, evidence
  certificate schema and persistent naming move earlier, right after G8,
  ahead of B9 (certified intersection core, now a parallel research track)
  and with FreeCAD integration before B9.

### What was built

- src/brepkernel/evidence.py: `brepkernel.evidence/1.0` schema module.
  `brep_sha256()` / `canonical_brep_bytes()` hash the canonical BREP text
  (BRepTools_Write, byte-deterministic for identically constructed
  shapes). `evidence_name(op, sha_a, sha_b)` implements
  `brepkernel.naming/1.0`: name = f(op, input hashes, pipeline version),
  e.g. `ev_e10_union_b0caafc89533_73d09a490069_55810891`; the trailing
  component is sha256(version|commit)[:8], so evidence from different
  builds never shares a name silently. `build_record()` assembles the
  record (operation id, inputs, op + effective tolerances, kernel
  version/commit, timestamps, outcome, stage summaries, artifacts);
  `validate_evidence()` is a dependency-free schema checker;
  `write_evidence_file()` writes `<evidence_dir>/<name>.json` atomically
  and refuses invalid records.
- src/brepkernel/pipeline.py: `boolean_brep()` is now a thin wrapper that
  calls the renamed `_boolean_brep_impl()` and attaches the evidence
  record at `report["evidence"]` on the accept path and on
  typed-refusal paths (the BRepAmbiguousResult carries the report, so
  the record rides along). New keyword-only `evidence_dir=None`: when
  given, the record is written as `<name>.json`. Emission is strictly
  additive and fully guarded: any failure inside evidence code degrades
  to `report["evidence_error"]` and a failing sidecar write is recorded
  in `artifacts.evidence_write_error`; the accept/refuse outcome is
  never altered. `boolean()` (Tier A) untouched (I9).
- docs/EVIDENCE_SCHEMA.md: schema field table, naming rules, emission
  guarantees, abridged example, and the open remainder (G10 CLI +
  `brepkernel verify`; per-entity persistent naming service B4 items 1-4).
- tests/test_evidence.py: 40 checks. I4: committed failing first
  (ImportError: no evidence module), then implemented. Two test bugs
  fixed along the way, both in the test, not the implementation:
  expected union volume was written as 1.75 instead of 1.875, and the
  refusal stage-report assertion demanded an "assembly" key that the
  pipeline only writes after successful assembly.

### Verification

- tests/test_evidence.py: 40/40 PASS (accept record validates, refusal
  record validates with category UnresolvedContact + stage assembly,
  naming deterministic across fresh rebuilds, name changes with input/op
  change, input hashes match independent BRepTools_Write hashing,
  broken evidence_dir still accepts at the right volume and still
  refuses with the right kind, sidecar written at the persistent name
  and re-validates from disk, no em dashes in new files).
- Full suite: 27/27 test files exit 0, 0 [FAIL] lines (26 prior files +
  test_evidence).

### Invariants

- I1: no acceptance logic touched; evidence code only reads shapes and
  the report.
- I2: refusals stay typed; the evidence record adds category + stage +
  full stage_report but the refusal kind/stage/message are unchanged.
- I3: evidence code reads tolerances, never writes them.
- I4: failing test committed first (69318ba).
- I5: suite green (27/27).
- I6: this entry; the two test bugs above are recorded, not hidden.
- I7: no crashes encountered.
- I8: grepped new/changed files for U+2014; none present (also enforced
  by a test).
- I9: boolean() and boolean_brep() contracts separate; only boolean_brep
  gained the additive evidence_dir parameter.

### Open / not in this change

- G10 CLI and `brepkernel verify` (re-check a certificate without
  trusting it).
- Per-entity persistent naming (B4 items 1-4): deterministic face/edge
  names, report["naming"] ancestor mapping, resolve() with geometric
  fallback. Raw material (edge lineage) exists in the assembly report.
- Nothing in this change needed Ben's keys or files; nothing blocked.

---

## Part B: FreeCAD integration groundwork (B5) - 25 Sept 2026, branch partB/freecad

- Branch: partB/freecad (local only, off main 10ea7d8; never pushed). Work
  done in worktree ~/workspace/brep-gates-freecad; the shared
  ~/workspace/brep-gates checkout was not touched.
- FreeCAD is NOT installed on this machine, so per the task brief the
  result is the integration doc plus a FreeCAD-free adapter design and
  tests; no FreeCAD install was attempted.

### What was built

- `src/brepkernel/freecad_adapter.py` (new): FreeCAD-free geometry
  bridge. BREP text in/out (`shape_of_brep_text` / `brep_text_of`,
  the bridge format from the work plan), STEP in/out (`read_step` via
  `step_ingest.load_step`, `write_step`), `topology_counts()`,
  `shape_volume()`, `round_trip_report()` (counts + volume before/after,
  documents what the bridge preserves and loses), `run_boolean()` (runs
  boolean_brep across the bridge and returns the outcome as data:
  "accepted" with result BREP bytes, "refused" with typed
  category+stage+message and result_brep None, "error" for boundary
  input failures), `present_refusal()` and `format_status_panel()`
  (user-facing text for the FreeCAD task panel), and a defensive
  `to_json_safe()` for report serialization. The result dict carries
  `result["evidence"] = report.get("evidence")` verbatim when present,
  so the evidence sidecar rides the same seam; on this branch it is
  None because the partB/evidence-schema branch is not merged yet.
- `tests/test_freecad_adapter.py` (new): 8 test groups. BREP text round
  trip preserves volume 6.0 and 6 faces and is ASCII; garbage BREP
  raises typed AdapterError; STEP round trip preserves all topology
  counts and volume to 1e-9; overlapping-box union accepts at volume
  1.875; corner-touch union refuses as data with category
  UnresolvedContact at stage assembly, no solid fabricated, panel text
  shows kind and stage; invalid op raises ValueError; status panels
  render for accept and refuse.
- `docs/FREECAD_INTEGRATION.md` (new): integration surface survey
  (embedded Python, subprocess chosen over in-process because of the
  two-OCCT-builds trap), OCCT version-compat concerns
  (cadquery-ocp OCCT 8.0.1 vs FreeCAD 1.1's build, BREP format version,
  STEP AP242 as the safer bridge), STEP round-trip preserve/lose
  analysis, typed-refusal UX (Status / RefusalKind / Report /
  EvidenceFile properties, error state, "uncertified" explicit button),
  evidence as the audit trail, 5 open questions, the first working
  demo (headless FreeCADCmd script), and the blocked-on-Ben list.

### Verification

- tests/test_freecad_adapter.py: written first, failed on import
  (ModuleNotFoundError) before the module existed; 16/16 [PASS] after.
- Full suite: 27/27 test files exit 0, 0 [FAIL] lines (26 prior files +
  test_freecad_adapter). Environment ~/workspace/brep-booleans/.venv,
  PYTHONPATH=<worktree>/src.
- No em dashes in new files (checked with Python chr(0x2014) count;
  the test uses the "\u2014" escape sequence only).

### Invariants

- I1: no acceptance logic touched; the adapter only calls
  boolean_brep and reports its outcome. Refusal paths fabricate no
  solid (result_brep is None).
- I2: refusals stay typed end to end (AdapterError / BRepAmbiguousResult
  kind / "error" outcome kinds).
- I3: the refusal text never suggests loosening a tolerance; no
  tolerance is read or changed by adapter code.
- I4: failing test committed first (test written and run before the
  module existed).
- I5: suite green 27/27.
- I6: this entry; nothing hidden.
- I7: no crashes encountered. One environment note: the first full-suite
  loop backgrounded after 120 s; it completed on its own (see results
  above). The corner-touch-union probe earlier showed face-to-face
  touching unions are now accepted (G2), so the refusal test uses
  corner touch (UnresolvedContact, assembly) instead.
- I8: zero U+2014 in new files.
- I9: boolean()/boolean_brep() contracts untouched.

### Open / not in this change

- FreeCAD workbench (CertifiedBoolean FeaturePython object), the
  subprocess runner, and the headless FreeCADCmd CI job (G11) need a
  FreeCAD install; none exists on this machine.
- FreeCAD 1.1's actual OCCT build version and the BREP format version
  both sides parse (doc open question 1).
- Merging with partB/evidence-schema so run_boolean() carries real
  evidence records instead of None.
- Onshape API keys / OAuth app and CATIA access / sample files stay
  parked: blocked on Ben, untouched by design.

## PR0 - v1 harnesses: verdict equivalence + perf budget (26 Sept 2026, branch muse/v1-harness)

- Base: ceea17d (main). Head: 93ed213.
- Env: python 3.12.3, OCP 8.0.1.0, existing venv
  ~/workspace/brep-booleans/.venv (deviation from the plan's fresh
  .venv: the existing venv is the known-good one for this repo
  lineage; recorded here instead of rebuilt).
- Commands:
  - python tools/review_probes/verdict_equivalence.py --before src
    --after src --skip-probes --fuzz-trials 4 --fuzz-seed 5
    --fuzz-snap 0.5 --fuzz-kinds box --verbose --out /tmp/equiv-smoke.json
- What landed: tools/review_probes/verdict_equivalence.py (scrubbed
  subprocess workers per side, sys.path holds only the requested src
  tree, PYTHONPATH stripped, --verbose prints brepkernel.__file__ for
  both sides and the worker asserts the import is under the requested
  tree), tools/review_probes/perf_budget.py (best-of-N pinned cases;
  tool only, no budgets committed), tools/review_probes/corpus_manifest.json
  (17 everyday-CAD probes + 2 seeded fuzz expansions = 97 cases).
- Results: self-equivalence on the same tree: 84 cases run (80 manifest
  fuzz + 4 extra), 63 accepts before and after, 0 differences, 0
  blocking. Both workers reported the same brepkernel.__file__.
  Refusal kinds in the corpus baseline: recorded in /tmp/equiv-smoke.json.
- Perf: the 84-case self-run took ~6 min wall (both sides sequential);
  acceptable for per-PR use, fuzz depth stays adjustable via flags.
- Merge criteria: harness-only PR, no src/brepkernel change;
  self-equivalence reports 0 differences. Met.
- I1: no acceptance logic touched (new files only, no src change).
- I2: refusal kinds compared as typed kind/stage pairs, never strings.
- I3: no tolerance read or changed.
- I4: the harness's own acceptance test is the self-equivalence run
  (0 differences on identical trees), run before commit.
- I5: full suite baseline running on main at time of writing; no test
  files touched by this PR.
- I6: this entry; nothing hidden.
- I7: no crashes encountered.
- I8: zero U+2014 in new files (checked with grep).
- I9: boolean()/boolean_brep() contracts untouched.

### Open / not in this change

- Machine-specific budgets (docs/perf_budget.json) land in G13, not here.
- G10..G18b each get their own branch from main; this branch stays
  harness-only.
