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
