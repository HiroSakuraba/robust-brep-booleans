# brepkernel prototype v0.1 — what was built

A working slice of the "B-reps that don't break" design (v0.1 design doc),
covering analytic solids (box, sphere, cylinder, cone = the design's Tier A).

## Architecture (maps to the design's 7 stages)

| Stage | Module | What it does |
|---|---|---|
| 0 | `ingest.py` | Audits input parameters; `ToleranceLedger` records every epsilon in one place (design G1) |
| 1 | `proxy.py` | `certified_proxy()`: triangle mesh with a *proven* chordal-error bound (asserted, not hoped) |
| 2 | `arrange.py` | Proxy boolean via manifold3d. Exact CGAL/libigl core slots in behind this interface |
| 3 | `classify.py` | **Tier A**: exact degeneracy detection from defining parameters (coincident planes, tangent cylinders) — no tolerance voting |
| 4 | `classify.py` | Winding-style classification with margins using **exact** implicit functions; ON-margin faces are never decided by the proxy |
| 5 | `assemble.py` | Per-face provenance + engine/classifier cross-check |
| 6 | `verify.py` | Independent verification: closure, Euler χ, orientation, Monte-Carlo field check against exact implicits, manifold3d cross-check |

## The core contract

`brepkernel.boolean(A, B, op)` either returns a result that passed **all**
Stage 6 checks, or raises `AmbiguousResult` carrying the partial mesh and
the full ambiguity report. It never silently returns a broken solid.

## Key design decisions (and where the prototype deviates)

1. **Exact implicits as ground truth.** Every analytic solid has a closed-form
   implicit function. Classification and the Monte-Carlo field check use it —
   so topology decisions don't depend on the arrangement engine's numerics.
2. **Margins, not epsilons.** A face whose centroid is within
   `2·(chordal_A + chordal_B)` of a surface is ON — undecided, reported,
   never voted on.
3. **Arrangement engine is swappable.** manifold3d (robust-inexact) stands in
   for the exact libigl/CGAL core: no pip wheel exists for pyigl on this
   machine (conda-only) and building CGAL+swig bindings from source was out
   of scope. `arrange()` is the seam; the J birth-triangle map plugs in there.
4. **Verification is engine-diverse.** The field check samples the *exact*
   implicits (no mesh involved); the cross-check rebuilds through manifold3d
   independently. A bug in the arrangement engine can't confirm itself.

## Test results (23 Sept 2026)

- `tests/test_metamorphic.py`: 13/13 PASS — A∪A=A, A−A=∅, disjoint identities,
  commutativity, translation invariance, inclusion–exclusion, cyl/sphere/cone smoke.
- `tests/test_degenerate.py`: 8/8 PASS, zero silent failures — coincident faces,
  tangent cylinders, point/corner contacts, nested boxes. Degeneracies are
  detected analytically (Tier A) and either resolved or explicitly reported.

## Known limits / next steps

- Freeform (NURBS) faces: not yet — needs OCCT STEP ingest + Tier B/C.
- Exact arrangement core: needs conda or a CGAL build for pyigl.
- Winding-number field (design Section 3): approximated here by exact
  implicits; libigl `fast_winding_number` when the wheel is available.
- `arrange()` currently trusts manifold3d's output topology; Stage 6 is the
  backstop (it caught nothing in tests, but it is the gate that matters).

## Reproduce

```
cd brep-booleans
.venv/bin/python tests/test_metamorphic.py
.venv/bin/python tests/test_degenerate.py
```
