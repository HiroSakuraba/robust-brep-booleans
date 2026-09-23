"""Stage 4: winding-number classification with margins (Tier A exact).

For analytic solids the implicit function is exact, so inside/outside
tests are exact too. The only inexactness is the *proxy*: a face of the
arrangement mesh can sit up to chordal_error away from the true surface.
Classification therefore carries a margin: a query point whose |implicit|
is below the margin is ON the surface and must not be decided by the
proxy -- it goes to Tier A degeneracy analysis instead.

Provenance labels per face of the result mesh:
  A_ONLY / B_ONLY / INTERIOR (removed by the op) / COINCIDENT
"""

import numpy as np

IN, OUT, ON = 1, -1, 0


def classify_point(p, solid, margin):
    """Exact implicit classification with an ON margin. p: (...,3)."""
    f = solid.implicit(p)
    label = np.where(f < -margin, IN, np.where(f > margin, OUT, ON))
    return label


def op_wants_face(in_a, in_b, op):
    """Does the boolean op keep a face with this (inA, inB) status?

    Faces are kept where exactly one side's membership flips the result.
    Uses the standard classification table for regularized booleans.
    """
    # truth tables: result membership as a function of (a, b)
    if op == "union":
        r = (in_a == IN) | (in_b == IN)
    elif op == "intersection":
        r = (in_a == IN) & (in_b == IN)
    elif op == "difference":  # A - B
        r = (in_a == IN) & (in_b != IN)
    else:
        raise ValueError(f"unknown op {op}")
    return r


def classify_faces(face_points, solidA, solidB, margin, op):
    """Classify each face centroid. Returns (labels, on_mask).

    labels: 'A_ONLY', 'B_ONLY', 'SHARED' per face (which input it came from
    is decided by the arrangement engine; here we record the exact
    (inA, inB) status used for verification).
    on_mask: faces whose centroid is ON either surface within margin --
    these need Tier A degeneracy analysis, not proxy voting.
    """
    in_a = classify_point(face_points, solidA, margin)
    in_b = classify_point(face_points, solidB, margin)
    on_mask = (in_a == ON) | (in_b == ON)
    keep = op_wants_face(in_a, in_b, op)
    return in_a, in_b, keep, on_mask


def tierA_degeneracy(solidA, solidB, margin):
    """Exact degeneracy detection between the *defining parameters*.

    Returns a list of findings like 'coincident_planes', 'tangent_cylinders',
    'coincident_cylindrical_surfaces'. Empty list = general position
    (up to the margin). This is Tier A: analytic, no tolerance voting.
    """
    findings = []
    from .solids import Box
    if isinstance(solidA, Box) and isinstance(solidB, Box):
        # coincident face planes: compared exactly as floats (Tier A)
        for ax in range(3):
            if solidA.lo[ax] == solidB.lo[ax]:
                findings.append(f"coincident_planes: A.lo[{ax}] == B.lo[{ax}]")
            if solidA.hi[ax] == solidB.hi[ax]:
                findings.append(f"coincident_planes: A.hi[{ax}] == B.hi[{ax}]")
            if solidA.lo[ax] == solidB.hi[ax]:
                findings.append(f"touching_planes: A.lo[{ax}] == B.hi[{ax}]")
            if solidA.hi[ax] == solidB.lo[ax]:
                findings.append(f"touching_planes: A.hi[{ax}] == B.lo[{ax}]")
    from .solids import Cylinder
    if isinstance(solidA, Cylinder) and isinstance(solidB, Cylinder):
        # parallel axes?
        if abs(abs(float(solidA.axis @ solidB.axis)) - 1.0) < 1e-12:
            # distance between axis lines
            w = solidB.base - solidA.base
            dist = float(np.linalg.norm(w - (w @ solidA.axis) * solidA.axis))
            rr = solidA.r + solidB.r
            if abs(dist - rr) <= margin:
                findings.append("tangent_cylinders: axis distance == rA + rB")
            if dist <= margin and abs(solidA.r - solidB.r) <= margin:
                findings.append("coincident_cylindrical_surfaces")
    return findings
