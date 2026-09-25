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
