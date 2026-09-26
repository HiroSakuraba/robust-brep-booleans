# experimental_b9: certified intersector research spike (B9)

EXPERIMENTAL. Parallel research track, explicitly allowed to fail.
Nothing here is imported by the shipped pipeline (`brepkernel`
Tier A/B/C); nothing here changes any accept/refuse behavior.

## What it is

A prototype toward the work plan's B9 "certified intersection core":
a section-curve computation that carries its own correctness
certificate instead of relying on OCCT tolerances plus refusal gates.

For one pair of analytic implicit surfaces (F, G) (plane, sphere,
cylinder, cone, torus: exactly the analytic kinds Tier B/C ingests)
on a bounded domain box, `certified_section.enclosure()` computes:

- a **rigorous superset enclosure**: boxes are discarded only with an
  interval-arithmetic proof that F or G keeps a nonzero sign on them,
  so no intersection branch can be missed (this is the certified
  answer to the plan's open completeness problem, for analytic pairs);
- per-box **transversality flags** (interval grad F x grad G avoids
  zero: no singular point inside);
- per-box **uniqueness certificates**: a parametric 2D Krawczyk step
  (interval implicit function theorem) proving the box holds exactly
  one regular arc, a graph over the dominant tangent coordinate;
  plus a **rigorous exclusion** (Krawczyk image disjoint from the
  box) that kills the false positives of the plain interval test,
  and a slicing-plane Krawczyk fallback;
- **unresolved** boxes where the minimum size was reached without a
  proof: the honest "refuse" of this prototype.

`occt_agreement()` runs the plan's "both intersectors" protocol:
every 3D sample of OCCT's section edges must lie inside the retained
enclosure (a sample outside is a soundness alarm), and retained-box
clusters with no OCCT sample are flagged as unexplained regions.

`exact_classify` gives exact (Fraction, no floating point) point
membership for analytic solids: an independent oracle for the sample
points the pipeline classifies with OCCT's BRepClass3d classifier
(review finding F4).

## What it is not

- Not wired into the pipeline. No gate reads these results.
- Not freeform: NURBS need the plan's Bezier-decomposition route
  first (listed as follow-up work in docs/B9_RESEARCH.md).
- Not a production interval library: rounding is handled by a fixed
  outward widening (documented in intervals.py), adequate for a spike.

## Running the spike tests

```
.venv/bin/python tests/test_experimental_b9.py
```

with `PYTHONPATH` including the repo `src` (the test inserts it
itself, following the repo's test convention).
