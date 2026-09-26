"""Exact point classification for analytic solids (EXPERIMENTAL, B9).

Review finding F4: BRepClass3d_SolidClassifier once returned IN for a
point more than 1.0 away from both input surfaces, and the kernel uses
that same classifier to decide which patches to keep. For analytic
solids (planes, spheres, cylinders, cones, tori) point membership is
decidable exactly: every bounding surface is a polynomial implicit
with rational coefficients, so the sign at a rational point is exact
integer arithmetic (fractions.Fraction). No tolerances, no floating
point, no classifier heuristics.

This module is the "exact-predicate classification of sample points"
direction from the B9 task: an independent oracle for the sample
points the pipeline classifies with OCCT. Where the two disagree
outside the tolerance band, OCCT is wrong (F4 class); inside the band,
the exact predicate abstains and the pipeline keeps its refusal
behavior unchanged.

Only what the spike needs is implemented: closed solids bounded by
the implicits in implicits.py, with outward-oriented faces.
"""

from __future__ import annotations

from fractions import Fraction

from . import implicits as I


def _fr(x) -> Fraction:
    return x if isinstance(x, Fraction) else Fraction(x)


class ExactSolid:
    """Closed solid as a list of outward-oriented implicit surfaces.

    A rational point p is IN when every surface evaluates <= 0, ON when
    every surface evaluates <= 0 and at least one is == 0, else OUT.
    """

    def __init__(self, name: str, surfaces: list):
        self.name = name
        self.surfaces = list(surfaces)

    def classify(self, p) -> str:
        p = tuple(_fr(v) for v in p)
        on = False
        for s in self.surfaces:
            v = s.exact_value(p)
            if v > 0:
                return "OUT"
            if v == 0:
                on = True
        return "ON" if on else "IN"


def exact_box(lo, hi) -> ExactSolid:
    """Axis-aligned box [lo, hi]."""
    lo = tuple(_fr(v) for v in lo)
    hi = tuple(_fr(v) for v in hi)
    surfs = []
    for i in range(3):
        n = [Fraction(0)] * 3
        n[i] = Fraction(-1)
        surfs.append(I.plane(lo, n))   # x_i >= lo_i  <=> -x_i + lo_i <= 0
        n = [Fraction(0)] * 3
        n[i] = Fraction(1)
        surfs.append(I.plane(hi, n))   # x_i <= hi_i
    return ExactSolid("box", surfs)


def exact_sphere_solid(center, radius) -> ExactSolid:
    return ExactSolid("sphere", [I.sphere(center, radius)])


def exact_cylinder_solid(point, axis, radius, z_lo, z_hi) -> ExactSolid:
    """Finite cylinder: lateral surface plus two cap planes.

    axis must be a unit vector; caps are perpendicular to axis.
    Only axis-aligned (z) cylinders are needed by the spike tests, but
    the construction is general.
    """
    a = tuple(_fr(v) for v in axis)
    p = tuple(_fr(v) for v in point)
    surfs = [I.cylinder(point, axis, radius)]
    # cap planes: outward normals -a at height z_lo, +a at z_hi.
    # plane through q with normal n: n.(x - q) <= 0 inside.
    q_lo = tuple(pi + ai * _fr(z_lo) for pi, ai in zip(p, a))
    q_hi = tuple(pi + ai * _fr(z_hi) for pi, ai in zip(p, a))
    surfs.append(I.plane(q_lo, tuple(-ai for ai in a)))
    surfs.append(I.plane(q_hi, a))
    return ExactSolid("cylinder", surfs)


def compare_with_occt(exact_solid: ExactSolid, occt_shape, points,
                      band: float) -> dict:
    """Classify rational points exactly and with OCCT; compare.

    Points whose exact |signed distance| proxy is inside `band` are
    skipped (the exact predicate abstains there; OCCT legitimately
    reports ON). Every other point must agree, or OCCT misclassified.

    Returns counts and the list of disagreements.
    """
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.gp import gp_Pnt
    from OCP.TopAbs import TopAbs_IN, TopAbs_ON, TopAbs_OUT

    agreed = skipped = 0
    disagreements = []
    for p in points:
        exact = exact_solid.classify(p)
        # Band test: minimum |f| over surfaces, scaled crudely by the
        # gradient norm at the point (float is fine for the band test;
        # the verdict itself stays exact).
        pf = tuple(float(v) for v in p)
        dists = []
        for s in exact_solid.surfaces:
            g = s.float_grad(pf)
            gn = sum(x * x for x in g) ** 0.5
            if gn > 0:
                dists.append(abs(s.float_value(pf)) / gn)
        if dists and min(dists) < band:
            skipped += 1
            continue
        clf = BRepClass3d_SolidClassifier(
            occt_shape, gp_Pnt(*pf), float(band))
        state = clf.State()
        occt = ("IN" if state == TopAbs_IN else
                "ON" if state == TopAbs_ON else
                "OUT" if state == TopAbs_OUT else "UNKNOWN")
        if exact == "ON":
            # Exactly on a face: OCCT may read either side; not evidence.
            agreed += 1
        elif occt == exact:
            agreed += 1
        else:
            disagreements.append((p, exact, occt))
    return {"agreed": agreed, "skipped": skipped,
            "disagreements": disagreements,
            "n": len(points)}
