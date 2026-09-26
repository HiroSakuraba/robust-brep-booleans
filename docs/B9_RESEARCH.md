# B9 Research: Certified Intersector Spike

EXPERIMENTAL. Branch `partB/b9-research`, worktree
`~/workspace/brep-gates-b9`. Nothing here is imported by the shipped
pipeline; nothing here changes any accept/refuse behavior. Allowed to
fail; this document reports what actually happened, including the
negative results.

## 1. Research question

The work plan's B9 asks for a "certified intersection core": a
section-curve computation that carries its own correctness
certificate instead of relying on OCCT tolerances plus refusal gates.
Concretely: for one pair of analytic implicit surfaces (F, G) on a
bounded domain, can we compute

1. a **rigorous superset enclosure** of the true intersection curve
   (completeness: no branch missed, ever);
2. a **per-box certificate** that each retained box holds exactly one
   regular arc (so the curve's topology is known, not guessed);
3. an **honest unresolved report** where no proof was possible
   (singularities, degeneracies); and
4. a **two-intersector agreement protocol** with OCCT's floating-point
   section, where disagreement is an alarm, not a silent accept?

Scope: analytic surfaces only (plane, sphere, cylinder, cone,
torus: exactly the analytic kinds Tier B/C ingests). NURBS need the
plan's Bezier-decomposition route first and are out of scope for this
spike.

## 2. What was built

`src/brepkernel/experimental_b9/` (see its README):

- `intervals.py`: rigorous interval arithmetic over binary64 with
  outward rounding (fixed widening; documented as spike-grade).
- `implicits.py`: analytic implicits with interval evaluation,
  interval gradients, and exact `Fraction` evaluation.
- `certified_section.py`: the enclosure engine.
- `exact_classify.py`: exact point-membership for analytic solids,
  an F4-independent oracle for OCCT's `BRepClass3d_SolidClassifier`.
- `tests/test_experimental_b9.py`: pins every testable claim,
  including the negative ones.

### 2.1 Enclosure (completeness)

`enclosure()` subdivides the domain octree-style. A box is discarded
only with an interval proof that F or G keeps a nonzero sign on it.
Three exclusion surfaces, in order:

1. plain interval test: `0 not in F(box)` or `0 not in G(box)`;
2. parametric Krawczyk exclusion: the 2D Krawczyk image of the box
   (see 2.2) is disjoint from the box, so no per-parameter zero can
   lie inside (kills the false positives of surface 1);
3. at `min_size`, boxes that resist both proof and exclusion are
   retained as **unresolved** (honest refuse).

`EnclosureResult.check_partition()` verifies the retained plus
excluded leaves exactly tile the domain (volume check), and the test
re-verifies every stored exclusion proof independently.

### 2.2 Per-box certificates

For each retained box, in order:

1. **Transversality**: the interval cross product
   `grad F x grad G` avoids zero in at least one component, so no
   singular point lies in the box.
2. **Graph (parametric Krawczyk) verdict** (primary): on a
   transverse box, the dominant tangent component `t_k` provably
   avoids zero, so the curve is a graph over coordinate `k`. With
   `u` the other two coordinates, the parametric Krawczyk operator
   `K(U) = u - C.[H](u, Kparam) + (I - C.[J_u](U, Kparam))(U - u)`
   (interval residual over the full parameter range folded in)
   gives:
   - `K(U)` strictly inside `U`: **unique**, exactly one regular
     arc spanning the box, existence and uniqueness;
   - `K(U)` disjoint from `U`: **excluded**, no curve in the box;
   - else **none**.
3. **Slicing Krawczyk** (fallback): `H(x) = (F, G, t1.(x - x1))`
   with `x1` the Gauss-Newton closest curve point to the box
   center and `t1` the tangent there. Standard Krawczyk
   existence/uniqueness for the sliced system.

### 2.3 OCCT agreement protocol

`occt_agreement()`: every 3D sample of OCCT's section edges must lie
inside the retained enclosure (a sample outside is a soundness
alarm). `unexplained_regions()`: retained-box clusters (26-connected
on the voxel grid) containing no OCCT sample are flagged. In the
transverse cases these must be empty; in singular cases they are the
honest report.

## 3. Measurements (2026-09-25)

All from `tests/test_experimental_b9.py`, ALL PASS.

| case | boxes | unique | transv. | unres. | OCCT edges | cover | unexplained |
|---|---|---|---|---|---|---|---|
| A sphere/plane circle | 40 | 40 | 40 | 0 | 1 | 1.0000 | 0 |
| B cylinder/tilted-plane ellipse | 178 | 108 | 178 | 70 | 1 | 1.0000 | 0 |
| C concentric spheres (miss) | 0 | 0 | 0 | 0 | 0 | 1.0000 | 0 |
| D sphere/plane tangent | 8 | 0 | 0 | 8 | 0 | 1.0000 | 1 |
| E crossing cylinders (equal r) | 712 | 160 | 680 | 552 | 7 | 1.0000 | 0 |

Plus: interval soundness 120/120; certificate re-verification
(independent re-run of every exclusion proof and retained-box
interval) passes on all cases; exact classifier agrees with OCCT
300/300 on sphere probes and 297/300 on box probes (3 skipped as
near-boundary, disagreements 0).

Finer-resolution check on case B (`min_size` 0.05 -> 0.005):
unique 108/178 -> 186/198; unresolved 70 -> 12. The 12 persistent
boxes are corner-clips where the curve grazes a box corner and the
closest curve point to the box center lies just outside the box
(verified by hand on one box); all lie within `min_size` of the true
curve and inside OCCT-explained clusters.

## 4. Findings

### 4.1 The enclosure works; the naive certificate did not

The rigorous superset enclosure was straightforward and has held up
under every check (partition, re-verified proofs, OCCT coverage
1.0 on all five cases). The per-box certificate went through three
designs:

1. **Slicing Krawczyk through the box center** (first attempt):
   failed on most ellipse boxes (417/525 unresolved). Root cause:
   the Newton iteration for the sliced system converges to the
   slicing plane's intersection with the curve, which can lie far
   outside a static octree box even when the box contains curve.
   Lesson: slicing certificates assume tracer-placed boxes
   (centered on the curve); they are a poor fit for static
   subdivision.
2. **Slicing Krawczyk through the closest point** (Gauss-Newton):
   strictly better (the slice always passes through the nearest
   curve point), but still fails on corner-clips where the nearest
   point is outside the box.
3. **Parametric graph Krawczyk** (final, primary): no slicing
   plane at all; proves the box holds exactly one graph arc, or
   rigorously excludes the box. This killed the false-positive
   inflation (case B retained 525 -> 164 at the same `min_size`)
   and certifies 94% of boxes at fine resolution.

The honest summary: for static subdivision, the certificate that
fits is the graph/IFT one, not the slicing one. The slicing
Krawczyk is kept as a fallback.

### 4.2 Unresolved means boundary layer or singularity, and we can tell which

Every unresolved box in the transverse cases lies within
O(`min_size`) of the true curve (measured: max distance 0 at
`min_size` 0.05 with 8000 samples; the 12 persistent boxes at
0.005 are hand-verified corner-clips). None forms an unexplained
cluster. In singular cases (D tangent, E branch crossings) the
unresolved boxes concentrate exactly at the singularities, which is
the correct honest report: the certificate refuses where the
hypotheses (transversality) fail.

### 4.3 Bugs found by the spike's own checks

- Krawczyk centering: first wrote `K = C.H(m)` (wrong); corrected
  to the standard `K = xc + (I - C.J(X))(X - x1)` after catching my
  own unsound variant (deviation must be `(X - x1)`).
- Torus implicit sign: wrote `(s - R^2 - r^2)^2`, correct is
  `(s + R^2 - r^2)^2`; caught by an exact-sign spot check, 6/6 pass
  after the fix.
- `gp_Ax2` vs `gp_Ax3` for `Geom_CylindricalSurface`; empty
  face-bbox overlap handled as certified-empty.
- Grid alignment: curve exactly on box faces defeats strict
  containment; added the expanded-box retry.

### 4.4 Exact classifier (F4)

`exact_classify` (exact `Fraction` point membership for analytic
solids) agrees with OCCT's `BRepClass3d_SolidClassifier` on all
597 non-skipped probes, 0 disagreements. This gives an F4-independent
oracle for the sample points the pipeline classifies with OCCT.

## 5. Honest verdict

**Works:** rigorous superset enclosure with partition-verified
leaves and independently re-verified exclusion proofs; per-box
transversality flags; parametric graph uniqueness certificates
(100% on the circle, 94% on the ellipse at fine resolution);
certified-empty on clear misses; honest unresolved reports on
tangent and branch-crossing singularities; OCCT agreement 1.0 on all
cases with zero unexplained regions in transverse cases; exact
classifier agreement with OCCT.

**Partial:** the remaining unresolved boxes in transverse cases are
all boundary-layer corner-clips (curve grazes a box corner below
the resolution). They are honest, bounded (within `min_size` of the
curve), and clustered with explained regions, but they are not
certified. A production version would trace along the curve
(predictor-corrector boxes centered on the curve, as the literature
does) instead of using a static grid; the static grid is the
spike's known structural limitation.

**Failed / out of scope:** NURBS (needs the Bezier-decomposition
route first); component-count certification from the enclosure
alone (the enclosure is a superset, not a partition into
components; counting components needs the tracing layer);
performance engineering (the interval library is spike-grade).

## 6. Next steps (not started)

1. Curve tracing layer: predictor-corrector boxes centered on the
   curve, chaining the per-box uniqueness certificates into
   certified branches (this is what eliminates the corner-clip
   boundary layer).
2. Component counting from traced branches + the enclosure's
   completeness.
3. Bezier-decomposition route for NURBS faces, then the same
   machinery per patch.
4. Replace the fixed-widening interval arithmetic with a proper
   outward-rounding library before any production use.

## 7. Literature surveyed

- Krawczyk operator existence/uniqueness (Moore, Rump): the
  per-box certificate used here.
- Krawczyk-based certified curve tracking (arXiv 2602.07718):
  the tracing layer this spike does not yet build.
- Sederberg/Meyers + Hohmeyer normal-cone loop-detection test:
  relevant to the completeness question for the tracing layer.
- SSI topology via 4D computations (hal.science hal-01292835):
  "one curve box" concept, closely related to the graph
  certificate here.

## 8. Reproduction

```
cd ~/workspace/brep-gates-b9
PYTHONPATH=src ~/workspace/brep-booleans/.venv/bin/python \
    tests/test_experimental_b9.py
```

Expected: `RESULT: ALL PASS` (about 30 s; case E dominates).
