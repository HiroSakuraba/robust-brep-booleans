"""Stages 3-4: Tier A exact analysis + per-face audit of the engine.

The old design classified every result face against *both* implicits, but
every result face lies on an input surface, so that always reported ON.
The fix: each result face is classified against the *other* solid only
(origin tracked per triangle from the engine), and the engine's
keep/discard decision for that face is audited against the exact implicit:

  union:        kept A-face must satisfy B_implicit >= 0; kept B-face: A >= 0
  intersection: kept A-face must satisfy B_implicit <= 0; kept B-face: A <= 0
  difference:   kept A-face must satisfy B_implicit >= 0; kept B-face: A <= 0

With the rigorous chordal margin m, evaluated at each kept face's centroid:
  s*f >  m  -> verified (centroid provably on the kept side; the face was
                cut from an input triangle along the mesh crossing, so the
                mesh-world classification transfers to the true solid up to
                the certified resolution band)
  s*f < -m  -> VIOLATION (provably on the wrong side; engine is wrong)
  otherwise -> AMBIGUOUS (centroid inside the margin band; cannot verify,
                blocks the result)

Note the honest boundary: near a surface-surface crossing the engine emits
sliver triangles hugging the crossing curve, and their centroids can land
inside the margin band even for clean transverse intersections. Those faces
block (safe refusal) rather than guess. Pinning crossings tighter than the
band is the deferred exact arrangement core's job (design Section 4).

Tier A (below) detects degeneracies from defining parameters with exact
float predicates -- no tolerance voting.
"""

import numpy as np


def audit_faces(F, V, origins, solidA, solidB, margin, op):
    """Audit the engine's keep decision per face against exact implicits.

    Returns dict with boolean masks: verified, ambiguous, violation, plus
    the raw other-solid implicit values.
    """
    F = np.asarray(F)
    n = len(F)
    if n == 0:
        z = np.zeros(0, dtype=bool)
        return {"verified": z, "ambiguous": z, "violation": z,
                "f_other": np.zeros(0), "origins": origins}
    centroids = V[F].mean(axis=1)
    fromA = (origins == "A")
    f_other = np.empty(n)
    if np.any(fromA):
        f_other[fromA] = solidB.implicit(centroids[fromA])
    if np.any(~fromA):
        f_other[~fromA] = solidA.implicit(centroids[~fromA])
    # required polarity s: keep ⟹ s * f_other >= 0
    if op == "union":
        s = np.ones(n)
    elif op == "intersection":
        s = -np.ones(n)
    elif op == "difference":
        s = np.where(fromA, 1.0, -1.0)
    else:
        raise ValueError(f"unknown op {op!r}")
    g = s * f_other
    verified = g > margin
    violation = g < -margin
    ambiguous = ~(verified | violation)
    return {"verified": verified, "ambiguous": ambiguous,
            "violation": violation, "f_other": f_other, "origins": origins}


def identical_inputs(solidA, solidB):
    """Exact fast path: bitwise-identical defining parameters."""
    return solidA.params_equal(solidB)


def _axis_relations(a_lo, a_hi, b_lo, b_hi):
    """Per-axis interval relation: 'disjoint', 'touch', or 'overlap'."""
    rel = []
    for ax in range(3):
        if a_hi[ax] < b_lo[ax] or b_hi[ax] < a_lo[ax]:
            rel.append("disjoint")
        elif a_hi[ax] == b_lo[ax] or b_hi[ax] == a_lo[ax]:
            rel.append("touch")
        else:
            rel.append("overlap")
    return rel


def tierA_degeneracy(solidA, solidB, margin):
    """Exact degeneracy detection from the *defining parameters*.

    Findings are informational: they explain why faces may be
    unverifiable, but they never block on their own. Blocking comes
    from the per-face audit and Stage 6. No tolerance voting: every
    predicate here is an exact comparison of float values.
    """
    from .solids import Box, Sphere, Cylinder
    findings = []
    if isinstance(solidA, Box) and isinstance(solidB, Box):
        rel = _axis_relations(solidA.lo, solidA.hi, solidB.lo, solidB.hi)
        if "disjoint" in rel:
            findings.append("disjoint_boxes")
        else:
            touches = [ax for ax, r in enumerate(rel) if r == "touch"]
            if touches:
                kinds = {3: "point", 2: "edge", 1: "face"}[len(touches)]
                findings.append(
                    f"touching_boxes: {kinds} contact on axes {touches}")
            # coincident face planes, only meaningful when closures meet
            for ax in range(3):
                if solidA.lo[ax] == solidB.lo[ax]:
                    findings.append(f"coincident_planes: A.lo[{ax}]==B.lo[{ax}]")
                if solidA.hi[ax] == solidB.hi[ax]:
                    findings.append(f"coincident_planes: A.hi[{ax}]==B.hi[{ax}]")
    elif isinstance(solidA, Sphere) and isinstance(solidB, Sphere):
        d2 = float(np.sum((solidA.center - solidB.center) ** 2))
        rs, rd = solidA.r + solidB.r, solidA.r - solidB.r
        if d2 == 0.0:
            findings.append("concentric_spheres")
        elif d2 == rs * rs:
            findings.append("tangent_spheres: externally tangent")
        elif d2 == rd * rd:
            findings.append("tangent_spheres: internally tangent")
    elif isinstance(solidA, Cylinder) and isinstance(solidB, Cylinder):
        cosang = abs(float(solidA.axis @ solidB.axis))
        if cosang == 1.0:  # exactly parallel
            w = solidB.base - solidA.base
            perp = w - (w @ solidA.axis) * solidA.axis
            dist2 = float(perp @ perp)
            rr = solidA.r + solidB.r
            if dist2 == rr * rr:
                findings.append("tangent_cylinders: axis distance == rA+rB")
            elif dist2 == 0.0 and solidA.r == solidB.r:
                findings.append("coincident_cylindrical_surfaces")
    return findings


def exact_separation(solidA, solidB):
    """Exactly decide disjointness where a closed-form predicate exists.

    Returns True (disjoint), False (not disjoint), or None (unknown).
    """
    from .solids import Box, Sphere
    if isinstance(solidA, Box) and isinstance(solidB, Box):
        rel = _axis_relations(solidA.lo, solidA.hi, solidB.lo, solidB.hi)
        if "disjoint" in rel:
            return True
        return False
    if isinstance(solidA, Sphere) and isinstance(solidB, Sphere):
        d2 = float(np.sum((solidA.center - solidB.center) ** 2))
        rs = solidA.r + solidB.r
        if d2 > rs * rs:
            return True
        return False
    # bbox-level: disjoint bboxes are exactly disjoint
    alo, ahi = solidA.bbox()
    blo, bhi = solidB.bbox()
    if np.any(ahi < blo) or np.any(bhi < alo):
        return True
    return None


def expected_shells(solidA, solidB, op):
    """Expected number of connected triangle-shells, or None if unknown."""
    from .solids import Box
    sep = exact_separation(solidA, solidB)
    if op == "union":
        if sep is True:
            return 2
        if sep is False:
            return 1
        return None
    if op == "intersection":
        return None  # decided by emptiness downstream
    if op == "difference":
        if sep is True:
            return 1  # A - B = A
        if isinstance(solidA, Box) and isinstance(solidB, Box):
            if (np.all(solidB.lo > solidA.lo)
                    and np.all(solidB.hi < solidA.hi)):
                return 2  # enclosed cavity: outer shell + inner shell
            return 1  # overlap / touch / B-contains-A(empty, handled upstream)
        return None  # relation not exactly decidable: report only
    raise ValueError(op)


def expected_euler(solidA, solidB, op):
    """Predicted Euler characteristic, or None if not exactly decidable.

    Implemented for box-box, where interval relations decide everything.
    """
    from .solids import Box
    if not (isinstance(solidA, Box) and isinstance(solidB, Box)):
        return None
    rel = _axis_relations(solidA.lo, solidA.hi, solidB.lo, solidB.hi)
    touches = sum(r == "touch" for r in rel)
    if op == "union":
        if "disjoint" in rel:
            return 4
        if touches == 3 or touches == 2:
            return 3  # point / edge contact: wedge of spheres
        return 2  # face contact or overlap
    if op == "intersection":
        if "disjoint" in rel or touches > 0:
            return None  # empty or degenerate; skip
        return 2
    if op == "difference":
        if "disjoint" in rel:
            return 2
        if np.all(solidB.lo > solidA.lo) and np.all(solidB.hi < solidA.hi):
            return 4  # cavity: two shells
        return 2
    raise ValueError(op)


def exact_op_volume(solidA, solidB, op):
    """Closed-form result volume, or None if not exactly decidable.

    Box-box: the intersection is a box (possibly empty/degenerate), so
    union / intersection / difference volumes are exact via inclusion.
    """
    from .solids import Box
    if not (isinstance(solidA, Box) and isinstance(solidB, Box)):
        return None
    lo = np.maximum(solidA.lo, solidB.lo)
    hi = np.minimum(solidA.hi, solidB.hi)
    v_inter = float(np.prod(np.maximum(hi - lo, 0.0)))
    vA, vB = solidA.volume_exact(), solidB.volume_exact()
    if op == "union":
        return vA + vB - v_inter
    if op == "intersection":
        return v_inter
    if op == "difference":
        return vA - v_inter
    raise ValueError(op)
