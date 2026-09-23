"""Metamorphic identity tests.

These hold for ANY correct boolean implementation, so they test the
pipeline without needing ground-truth meshes. Box cases assert EXACT
volumes (boxes are exactly representable); curved cases allow the
certified proxy tolerance. Translation uses fractional offsets so float
rounding is actually exercised.
"""

import sys
import numpy as np

sys.path.insert(0, "src")
from brepkernel import boolean, AmbiguousResult
from brepkernel.solids import Box, Sphere, Cylinder, Cone
from brepkernel.verify import signed_volume


def vol(mesh):
    return abs(signed_volume(mesh["V"], mesh["F"]))


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"[{status}] {name} {detail}")
    return cond


def main():
    ok = True
    A = Box([0, 0, 0], [2, 2, 2])                    # vol 8
    B = Sphere([4, 0, 0], 1.0)                       # disjoint from A
    C = Box([1, 1, 1], [3, 3, 3])                    # overlaps A

    # 1. A u A = A (Tier A exact resolution; volume exact)
    m, r = boolean(A, A, "union")
    ok &= check("A u A = A (volume exact)", vol(m) == 8.0,
                f"vol={vol(m)!r}")
    ok &= check("A u A accepted", r["accepted"])

    # 2. A - A = empty (exact)
    m, r = boolean(A, A, "difference")
    ok &= check("A - A empty", m["empty"] or len(m["F"]) == 0)
    ok &= check("A - A accepted", r["accepted"])

    # 3. disjoint union: volumes add up to the certified proxy tolerance
    # (the sphere proxy is chordal-bounded, not exact)
    m1, _ = boolean(A, B, "union")
    ok &= check("disjoint union volume adds",
                abs(vol(m1) - (A.volume_exact() + B.volume_exact())) < 0.05,
                f"vol={vol(m1):.6f}")

    # 4. commutativity (exact for boxes)
    mu1, _ = boolean(A, C, "union")
    mu2, _ = boolean(C, A, "union")
    ok &= check("union commutes (exact)", vol(mu1) == vol(mu2),
                f"{vol(mu1)!r} vs {vol(mu2)!r}")
    mi1, _ = boolean(A, C, "intersection")
    mi2, _ = boolean(C, A, "intersection")
    ok &= check("intersection commutes (exact)", vol(mi1) == vol(mi2),
                f"{vol(mi1)!r} vs {vol(mi2)!r}")
    ok &= check("union volume exact 15", vol(mu1) == 15.0,
                f"vol={vol(mu1)!r}")

    # 5. disjoint intersection is empty
    m, r = boolean(A, B, "intersection")
    ok &= check("disjoint intersection empty", m["empty"] or len(m["F"]) == 0)
    ok &= check("disjoint intersection accepted", r["accepted"])

    # 6. translation invariance with FRACTIONAL offsets (exercises rounding)
    t = np.array([0.1, -0.3, 0.7])
    At = Box(A.lo + t, A.hi + t)
    Ct = Box(C.lo + t, C.hi + t)
    mt, _ = boolean(At, Ct, "union")
    ok &= check("translation invariance (fractional)",
                abs(vol(mt) - vol(mu1)) < 1e-9,
                f"{vol(mt)!r} vs {vol(mu1)!r}")

    # 7. cylinder/sphere/cone smoke: accepted, sane volumes
    cyl = Cylinder([0, 0, -1], [0, 0, 1], 1.0, 2.0)
    sph = Sphere([0, 0, 0], 1.5)
    m, r = boolean(cyl, sph, "union")
    sane = cyl.volume_exact() < vol(m) < cyl.volume_exact() + sph.volume_exact()
    ok &= check("cyl u sph accepted+sane", r["accepted"] and sane,
                f"vol={vol(m):.4f}")
    m, r = boolean(cyl, Sphere([0, 0, 0], 0.5), "difference")
    ok &= check("cyl - sph accepted", r["accepted"]
                and 0 < vol(m) < cyl.volume_exact(), f"vol={vol(m):.4f}")
    cone = Cone([0, 0, 0], [0, 0, 1], 1.0, 2.0)
    m, r = boolean(cone, Box([-2, -2, -2], [2, 2, 2]), "intersection")
    ok &= check("cone n box accepted", r["accepted"], f"vol={vol(m):.4f}")

    # 8. inclusion-exclusion (exact for boxes)
    lhs = vol(mu1)
    rhs = A.volume_exact() + C.volume_exact() - vol(mi1)
    ok &= check("inclusion-exclusion (exact)", lhs == rhs,
                f"{lhs!r} vs {rhs!r}")

    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
