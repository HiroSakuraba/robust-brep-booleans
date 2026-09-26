# Release checklist: brepkernel v0.9.0

Target state: the merged v0.9 + 7 gate branches on local main (`10ea7d8`
plus the shared-seam fix `59fe570`), with the G9 docs branch merged in.
Nothing here has been pushed: local main is 55 commits ahead of
`origin/main`, and no tag or release exists yet.

## Proposed version tag: v0.9.0

Rationale: the remote tags `v0.1` through `v0.4` cover the old mesh/Tier A
prototype lineage. The names `v0.8` / `v0.9` on this branch exist only as
merge-commit names ("Merge v0.8 exact NURBS Boolean pipeline", "Merge v0.9
verified periodic seam-aware NURBS splitting"). The work plan's G9 names
the release v0.9.0, covering the v0.9 seam-aware NURBS pipeline plus the
seven merged gates: G0 tooling/CI/baseline, G1 independent Boolean arbiter
+ rework, G2 planar coincidence (three-case rule), G3/G4 completeness
probe reworks, G5 multi-ray parity second classifier, G6 per-face
broad-phase pads + rework, probe hardening, and the G7 seam stress round.

## Asset list: brepkernel-v0.9.0.zip

Source zip built from the tagged tree. Contents:

- `src/` (the `brepkernel` package; the release payload)
- `tests/` (26 test files, including `tests/_arbiter.py` and
  `tests/data/review_20260925/`)
- `tools/` (review probes: `arbiter.py`, `fuzz_brep.py`,
  `common_cad_probes.py`, `repro_findings.py`)
- `docs/` (README stays at repo root; `docs/` carries `REVIEW_LEDGER.md`,
  `PROTOTYPE.md`, `NURBS_ACCEL.md`, `RELEASE_CHECKLIST.md`,
  `baseline_20260925/`)
- `README.md`, `LICENSE`, `requirements-freeform.txt`
- `.github/workflows/` (CI definitions; validated locally, remote green
  pending a push)

Excluded: `.venv/`, `__pycache__/`, `*.pyc`, `.git/`,
`probes-abc-deepcad/` (stress-round scratch from the earlier ABC/DeepCAD
round, kept in the repo but not in the release asset).

## Pre-release verification (all local)

Run in order; each must pass before the next step:

1. Working tree clean on the release commit:
   `git status --porcelain` empty; `git log -1 --format=%H` is the tagged
   commit.
2. Em-dash scan (invariant I8): `git grep -P '\x{2014}' -- src tests
   tools docs README.md` returns 0 matches.
3. Full suite: `for t in tests/test_*.py; do PYTHONPATH=src
   .venv/bin/python "$t" || echo "FAILED: $t"; done` with the venv at
   `~/workspace/brep-booleans/.venv`. Required: all 26 files exit 0, 0
   `[FAIL]` lines. (Reference run: 26/26 on merged main, 2026-09-25.)
4. Review probes: `PYTHONPATH=src .venv/bin/python
   tools/review_probes/common_cad_probes.py` must exit 0 with `wrong=0`,
   `unmet_accepts=0`, `unexpected_accepts=0`; `PYTHONPATH=src
   .venv/bin/python tools/review_probes/repro_findings.py` must exit 0
   (F2 accepted with OCCT volume agreement within 1e-6 and arbiter
   `kernel_errors == 0`; F4 assertions pass).
5. Ledger: the G9 entry in `docs/REVIEW_LEDGER.md` records the final
   suite/probe results and the tagged commit.
6. Spot-check the zip: unzip to a temp dir, confirm `src/brepkernel`,
   `tests/`, `tools/`, `docs/`, `README.md`, `LICENSE`,
   `requirements-freeform.txt` are present and the venv is not.

## Blocked on Ben (only he can do these)

- Review and approve the G9 docs branch (`gate/G9-docs`) for merge into
  local main.
- Push local main to `origin/main` (`HiroSakuraba/robust-brep-booleans`,
  public). Local main is 55 commits ahead; the push also unblocks remote
  CI verification, which the G0 entry records as unverified by
  construction so far.
- Create the `v0.9.0` tag on the release commit and push the tag.
  (No tag was created locally for this checklist; tagging is left to Ben.)
- Create the GitHub release `v0.9.0` and attach `brepkernel-v0.9.0.zip`.
- Decide invariant I9: keep `boolean()` (Tier A mesh route) and
  `boolean_brep()` (Tier B/C) contracts separate, or add a dispatcher
  `boolean_any()` that picks the route by input type. The work plan
  defers this decision to G9.
- Decide the release-notes wording for the known-open items: G8
  real-data refusal rate (separate worker, in progress) and the deferred
  curved/NURBS coincidence follow-up gate.

## Out of scope for this release (recorded, not done)

- `CHANGELOG.md`: the work plan's G9 asks for one; not yet written.
- `docs/NURBS_ACCEL.md` "What remains": not yet updated to reflect the
  closed gates.
- G8 measured numbers: not available for the release notes.
