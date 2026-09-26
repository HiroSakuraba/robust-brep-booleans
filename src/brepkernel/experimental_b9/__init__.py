"""B9 certified intersector research spike (EXPERIMENTAL).

Nothing in this package is imported by the shipped pipeline. It is a
parallel research track toward the plan's B9 "certified intersection
core": a section-curve computation that carries its own correctness
certificate instead of relying on OCCT tolerances plus refusal gates.

Contents:

- intervals.py: rigorous interval arithmetic over binary64 with
  outward rounding.
- implicits.py: analytic implicit surfaces (plane, sphere, cylinder,
  cone, torus) with interval and exact (Fraction) evaluation.
- certified_section.py: exclusion-subdivision enclosure (rigorous
  superset of the true intersection), per-box transversality flags,
  per-box Krawczyk existence/uniqueness certificates, and the OCCT
  agreement protocol.
- exact_classify.py: exact point-membership predicates for analytic
  solids, an independent cross-check of OCCT's BRepClass3d classifier
  (review finding F4).

Status: research prototype. Allowed to fail. See docs/B9_RESEARCH.md
for the question, the measurements, and the honest verdict.
"""

from . import certified_section, exact_classify, implicits, intervals

__all__ = ["certified_section", "exact_classify", "implicits", "intervals"]
