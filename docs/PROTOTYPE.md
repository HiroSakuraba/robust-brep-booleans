# brepkernel prototype v0.2: what was built

NOTE (2026-09-25): this document describes the v0.2 Tier A slice, which is
now the legacy route 1 (`boolean()` in `src/brepkernel/pipeline.py`). The
current pipeline, route 2 (`boolean_brep()`: exact trimmed B-rep, Tier B/C,
gates G0 through G7 merged 2026-09-25), is documented in `../README.md` and
`REVIEW_LEDGER.md`. The text below is kept as the historical record of the
mesh-era prototype.

A working slice of the "B-reps that don't break" design, covering analytic
solids (box, sphere, cylinder, cone: the design's Tier A).

v0.2 reworks the prototype around a review finding: in v0.1 the mesh engine
made every topology decision while the "exact" layer only wrote notes, and
the final gate checked volume rather than shape. Now the engine does the
geometry, the exact layer audits every topology decision per face, and
anything unverifiable blocks the result.

## Architecture (maps to the design's 7 stages)

| Stage | Module | What it does |
|---|---|---|
| 0 | `ingest.py` | Audits input parameters; `ToleranceLedger` records every epsilon in one place (design G1) |
| 1 | `proxy.py` | `certified_proxy()`: triangle mesh with a proven chordal-error bound computed from the actual mesh, not a formula (raises, never asserts) |
| 2 | `arrange.py` | Proxy boolean in float64 (`manifold3d` Mesh64) with per-face origin tags (A: 0..nA-1, B: nA..). Exact CGAL/libigl core slots in behind this interface |
| 3 | `classify.py` | **Tier A**: exact degeneracy detection from defining parameters (no tolerance voting). Exact helpers: `exact_op_volume`, `expected_shells`, `expected_euler` (box-box), `identical_inputs` |
| 4 | `classify.py` | Per-face audit against the *other* solid's exact implicit, with per-operation polarity (union: kept A-face needs B_implicit >= 0, etc.). Verified / ambiguous / violation, per face |
| 5 | `assemble.py` | Attaches provenance and the audit to the arrangement mesh |
| 6 | `verify.py` | 8 checks: V1 directed-edge closure, V2 Euler vs exact prediction, V3 orientation, V4 exact volume (box-box) or labeled-statistical MC, V5 OCCT cross-check (independent engine and kernel), V6 vertex-on-surface, V7 sample membership (own ray caster vs exact implicits), V8 shell count vs expected |

## The core contract

`brepkernel.boolean(A, B, op)` returns a result only if every kept face
verified against the exact implicits and all Stage 6 checks passed.
Otherwise it raises `AmbiguousResult` carrying the partial mesh and the
full report. Ambiguities block: unverifiable faces, engine-decision
violations, and failed checks are never returned as success. Identical
inputs (A op A) resolve exactly via Tier A without touching the engine.

## Key design decisions

1. **Exact implicits audit, not just observe.** Every kept face is checked
   against the other solid's closed-form implicit. The v0.1 audit compared
   against both implicits, which flags every face ON (every result face lies
   on an input surface); v0.2 checks each face against the *other* solid
   only, which is the decision that matters.
2. **Origin tracking.** `arrange()` tags every result face with the input
   face it came from (`face_id` round-trips through the boolean), so the
   audit knows which implicit is "mine" and which is "other".
3. **Margins are principled.** The degeneracy margin is the sum of the two
   certified chordal bounds plus a 1e-9 numeric epsilon, not a tuned
   constant. The sphere/cylinder certificates are worst-case bounds derived
   from the actual tessellation (verified against 300k-sample brute force:
   no understatement).
4. **float64 end to end.** The engine runs Mesh64; Tier A reasons in float64.
   The v0.1 float32 engine fused boxes 1e-9 apart into one solid.
5. **Shape is checked, not just volume.** V6 requires every output vertex
   within the margin of an input surface; V7 requires point-in-mesh to agree
   with exact membership on samples outside the margin band. A
   volume-preserving corner push is rejected.
6. **Cross-check is engine-diverse.** V5 rebuilds the boolean from OCCT
   primitives (independent geometry kernel) and compares volumes. The MC
   field check is kept only where no closed form exists, and is labeled
   statistical.

## Test results (23 Sept 2026)

- `tests/test_metamorphic.py`: 15/15 PASS. Box cases assert exact volumes
  (bitwise); translation uses fractional offsets (0.1, -0.3, 0.7) so float
  rounding is exercised; curved smoke tests must be accepted with sane
  volumes.
- `tests/test_degenerate.py`: 8/8 PASS. Each case declares its expected
  disposition ("accepted" with exact volume, or "ambiguous" raising
  `AmbiguousResult`). Tangent cylinders, point-touching spheres, and
  corner-touching boxes are rejected explicitly; coincident-face union and
  coplanar difference are accepted with exact volumes and all faces
  verified.
- `tests/test_stress.py`: 28/28 PASS (added 23 Sept 2026). Adversarial
  battery: near-degenerate box gaps (accepted at 1e-8, ambiguous at 5e-10
  and 1e-12; the margin is 1e-9 * coord_scale, so 2e-9 sits below it and is
  correctly refused), a 2-micron sliver intersection (exact volume 2e-6),
  nested spheres (union and cavity difference), grazing contacts, cone apex
  on a box face (accepted: 2 clean shells, chi = 4, OCCT agrees), a 3-box
  chain (exact 1.875 / 1.973), 6 randomized box-soup trials (commutativity
  + no-silent-failure invariant), invalid inputs refused loudly, large
  coordinate offsets (vol exactly 15.0), identical spheres, inscribed sphere
  (union accepted with vol exactly 8.0; difference ambiguous), and box with
  cubic cavity (exact 0.875). The headline: zero silent failures found by
  this battery. The tilted/axis-aligned cylinder-through-box cases are ACCEPTED via
  patch-level classification (see below); true tangencies (tangent
  cylinder/sphere probes) are still refused.

## Known limits / next steps

- The audit cannot pin intersection-curve locations tighter than the margin
  band; faces fully inside the band block rather than guess. Patch-level
  classification recovers the important transverse cases: a cylinder through
  a box (tilted or axis-aligned) is now ACCEPTED, because the sliver
  triangles hugging the crossing ellipse belong to patches that touch an
  A/B crossing and contain decisively verified faces (|f| ~ 0.6 * margin at
  tol 1e-3 and 1e-4, so refining does not help; the residual error is
  geometric, bounded by the margin, and cannot change topology). A whole
  patch stuck in the band with no decisive face -- the tangency case -- is
  still refused.
- The audit itself is fast: each face gets a rigorous whole-face lower
  bound from the solid's analytic form (box: max of per-axis pieces;
  sphere/cylinder/cone: point-triangle distance), which verifies most faces
  in microseconds. The audit-only subdivision runs solely as a
  violation-prover on faces whose bound dips below -margin. The tilted
  cylinder-through-box audit runs in 0.03 s (was minutes with subdivision
  alone).
- V5 compares volumes, so it catches gross errors only; V6/V7 cover shape.
- Mixed-kind exact relations (box vs sphere containment, etc.) are not
  decided; those checks stay report-only.
- Freeform (NURBS) faces: not yet, needs OCCT STEP ingest + Tier B/C.
- Exact arrangement core: needs conda or a CGAL build for pyigl.
- Winding-number field (design Section 3): approximated here by exact
  implicits; libigl `fast_winding_number` when the wheel is available.

## Reproduce

```
cd brep-booleans
python3 -m venv .venv
.venv/bin/pip install manifold3d cadquery-ocp
.venv/bin/python tests/test_metamorphic.py
.venv/bin/python tests/test_degenerate.py
.venv/bin/python tests/test_regression.py
.venv/bin/python tests/test_stress.py
```
