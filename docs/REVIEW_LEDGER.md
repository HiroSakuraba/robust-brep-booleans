# Review Ledger

Running record of review-driven gate work: commands run, trimmed output,
and PASS/FAIL per item. Conventions for this file:

- I7 (crash rule): any crash stops the run; the traceback and the inputs
  that triggered it are recorded here before anything else is done.
- I8 (no em dashes): no U+2014 anywhere in changed files or in this file.
- Nothing is pushed to GitHub; all branches and commits stay local.

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
