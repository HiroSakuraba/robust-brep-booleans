# robust-brep-booleans

A prototype of robust B-rep Boolean operations for CAD. Geometry comes from a
fast mesh engine; every topological verdict is then audited against exact
analytic implicits, and anything that cannot be verified is reported, never
guessed.

This implements the Tier A (analytic solids) slice of the design in `docs/`.

## The idea in one paragraph

Boolean operations on boundary representations fail in practice because
floating-point geometry is used to make topological decisions (is this face
inside or outside?). This prototype splits the job: a mesh engine
(manifold3d, float64) computes the arrangement geometry, then an exact layer
audits every kept face against the closed-form implicits of the input solids.
Faces the exact layer cannot decide are blocked, not guessed. The result is a
pipeline that either returns a verified solid or raises `AmbiguousResult`
with the partial mesh and a full report.

## Pipeline (7 stages)

| Stage | Module | What it does |
|---|---|---|
| 0 | `ingest.py` | Audits input parameters; `ToleranceLedger` records every epsilon in one place |
| 1 | `proxy.py` | `certified_proxy()`: mesh with a *proven* chordal-error bound, computed from the actual mesh (raises, never asserts) |
| 2 | `arrange.py` | Proxy boolean in float64 (`manifold3d` Mesh64) with per-face origin tags. The exact CGAL/libigl core plugs in behind this interface |
| 3 | `classify.py` | **Tier A**: exact degeneracy detection from defining parameters (touching boxes, tangent spheres/cylinders). Exact closed-form volumes, shell counts, and Euler predictions for box-box |
| 4 | `classify.py` | Per-face audit: each kept face is checked against the *other* solid's exact implicit, with per-operation polarity. `s*f > margin` is verified, `< -margin` is a violation, inside the band is ambiguous and blocks |
| 5 | `assemble.py` | Attaches provenance and the audit to the arrangement mesh |
| 6 | `verify.py` | Independent verification, 8 checks: directed-edge closure, Euler vs prediction, orientation, exact/MC volume, OCCT cross-check, vertex-on-surface, sample membership, shell count |

## The contract

```python
from brepkernel import boolean, AmbiguousResult
from brepkernel.solids import Box, Sphere

try:
    mesh, report = boolean(Box([0,0,0],[2,2,2]), Sphere([1,0,0], 1.5), "union")
except AmbiguousResult as e:
    # e.mesh: partial result. e.report: full stage-by-stage report.
    ...
```

`boolean()` returns a mesh only if every kept face verified against the exact
implicits and all Stage 6 checks passed. Otherwise it raises
`AmbiguousResult`. Degenerate inputs (tangent spheres, point contacts) are
rejected explicitly. Identical inputs (A op A) resolve exactly without the
engine.

## Quickstart

```bash
python3 -m venv .venv
.venv/bin/pip install manifold3d cadquery-ocp
.venv/bin/python tests/test_metamorphic.py
.venv/bin/python tests/test_degenerate.py
.venv/bin/python tests/test_regression.py
```

`cadquery-ocp` is needed for the OCCT cross-check (check V5); without it that
check is skipped with a note. `rhino3dm` is needed only for the future
STEP/freeform slice.

## Test results (23 Sept 2026)

- `tests/test_metamorphic.py`: 15/15 pass. Box cases assert exact volumes;
  translation uses fractional offsets so float rounding is exercised.
- `tests/test_degenerate.py`: 8/8 pass. Each case declares its expected
  disposition: certifiable cases must return exact volumes, degenerate cases
  (tangent cylinders/spheres, point contacts) must raise `AmbiguousResult`.
- `tests/test_regression.py`: 5/5 pass, and 5/5 fail on the pre-fix code.
  Covers the 1e-9-apart boxes that used to fuse silently, a volume-preserving
  corner-push shape attack, non-round coordinates, the sphere certificate
  bound, and flipped triangle winding.

## Honest limits

- The mesh engine does the geometric work (intersection curves, face
  splitting); the exact layer audits its topology decisions and blocks what
  it cannot verify. The engine is not exact, and the audit cannot pin cut
  locations tighter than the margin band.
- The OCCT cross-check is engine-diverse (independent kernel, independent
  geometry) but compares volumes, so it catches gross errors, not small
  shape deviations. Small deviations are caught by the vertex-on-surface
  and sample-membership checks instead.
- The Monte-Carlo field check is statistical, used only where no closed-form
  volume exists (non-box pairs). Box pairs are checked against exact volumes.
- Freeform (NURBS) faces, Tier B/C, are not in this slice. Next build: OCCT
  STEP ingest.
- Exact arrangement core (libigl/CGAL): no pyigl wheel on PyPI (conda-only).
  `arrange()` is the seam.

See `docs/PROTOTYPE.md` for the full write-up.

## License

MIT, see `LICENSE`.
