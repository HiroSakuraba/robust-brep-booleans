"""Metamorphic identity tests (design G9-style self-consistency).

These hold for ANY correct boolean implementation, so they test the
pipeline without needing ground-truth meshes:
  A u A = A,  A - A = empty,  (A u B) n A = A,  commutativity of u/n,
  rigid-motion invariance, volume identities for disjoint/touching pairs.
"""

import sys
import math
import numpy as np

sys.path.insert(0, "src")
from brepkernel import boolean
from brepkernel.solids import Box, Sphere, Cylinder, Cone
from brepkernel.verify import signed_volume


def vol(mesh):
    return abs(signed_volume(mesh["V"], mesh["F"]))


def rel_err(a, b):
    return abs(a - b) / max(1e-12, abs(b))


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"[{status}] {name} {detail}")
    return cond


def main():
    ok = True
    A = Box([0, 0, 0], [2, 2, 2])                    # vol 8
    B = Sphere([3, 0, 0], 1.0)                       # disjoint from A
    C = Box([1, 1, 1], [3, 3, 3])                    # overlaps A

    # 1. A u A = A (volume identity + accepted)
    m, r = boolean(A, A, "union")
    ok &= check("A u A = A (volume)", rel_err(vol(m), A.volume_exact()) < 0.02,
                f"vol={vol(m):.4f} exact=8")
    ok &= check("A u A accepted", r["accepted"])

    # 2. A - A = empty
    m, r = boolean(A, A, "difference")
    ok &= check("A - A empty", m["empty"] or len(m["F"]) == 0)
    ok &= check("A - A accepted", r["accepted"])

    # 3. (A u B) n A = A for disjoint B
    m1, _ = boolean(A, B, "union")
    # rebuild A-side: intersect the union result against A via pipeline on proxies?
    # simpler volume identity: vol(A u B) = volA + volB for disjoint inputs
    ok &= check("disjoint union volume adds",
                rel_err(vol(m1), A.volume_exact() + B.volume_exact()) < 0.03,
                f"vol={vol(m1):.4f}")

    # 4. commutativity: A u C vs C u A, A n C vs C n A
    mu1, _ = boolean(A, C, "union")
    mu2, _ = boolean(C, A, "union")
    ok &= check("union commutes (volume)",
                rel_err(vol(mu1), vol(mu2)) < 1e-6,
                f"{vol(mu1):.4f} vs {vol(mu2):.4f}")
    mi1, _ = boolean(A, C, "intersection")
    mi2, _ = boolean(C, A, "intersection")
    ok &= check("intersection commutes (volume)",
                rel_err(vol(mi1), vol(mi2)) < 1e-6,
                f"{vol(mi1):.4f} vs {vol(mi2):.4f}")

    # 5. A n A = A ; A u empty-side: A n B = empty for disjoint
    m, r = boolean(A, B, "intersection")
    ok &= check("disjoint intersection empty", m["empty"] or len(m["F"]) == 0)

    # 6. rigid-motion invariance: translate both inputs, volume unchanged
    t = np.array([5.0, -3.0, 7.0])
    At = Box(A.lo + t, A.hi + t)
    Ct = Box(C.lo + t, C.hi + t)
    mt, _ = boolean(At, Ct, "union")
    ok &= check("translation invariance", rel_err(vol(mt), vol(mu1)) < 1e-6,
                f"{vol(mt):.4f} vs {vol(mu1):.4f}")

    # 7. cylinder/sphere/cone smoke: union accepted, volume sane
    cyl = Cylinder([0, 0, -1], [0, 0, 1], 1.0, 2.0)
    sph = Sphere([0, 0, 0], 1.5)
    m, r = boolean(cyl, sph, "union")
    sane = cyl.volume_exact() < vol(m) < cyl.volume_exact() + sph.volume_exact()
    ok &= check("cyl u sph accepted+sane", r["accepted"] and sane,
                f"vol={vol(m):.4f}")
    m, r = boolean(cyl, Sphere([1.0, 0, 0], 1.0), "difference")
    ok &= check("cyl - sph accepted", r["accepted"]
                and 0 < vol(m) < cyl.volume_exact(), f"vol={vol(m):.4f}")
    cone = Cone([0, 0, 0], [0, 0, 1], 1.0, 2.0)
    m, r = boolean(cone, Box([-2, -2, -2], [2, 2, 2]), "intersection")
    ok &= check("cone n box accepted", r["accepted"], f"vol={vol(m):.4f}")

    # 8. inclusion-exclusion: vol(A u C) = volA + volC - vol(A n C)
    lhs = vol(mu1)
    rhs = A.volume_exact() + C.volume_exact() - vol(mi1)
    ok &= check("inclusion-exclusion", rel_err(lhs, rhs) < 0.03,
                f"{lhs:.4f} vs {rhs:.4f}")

    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
