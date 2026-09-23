"""Degenerate placement tests (design Stage 3 / Tier A).

Each case declares its EXPECTED disposition:
  - "accepted": must return a certified mesh; exact volumes checked.
  - "ambiguous": must raise AmbiguousResult with a blocking report
    (degenerate geometry the pipeline cannot certify -- reported,
    never silently guessed).

Anything else -- a crash, a silent accept of an ambiguous case, or a
rejection of a certifiable case -- is a failure.
"""

import sys
import numpy as np

sys.path.insert(0, "src")
from brepkernel import boolean, AmbiguousResult
from brepkernel.solids import Box, Sphere, Cylinder
from brepkernel.verify import signed_volume


def vol(mesh):
    return abs(signed_volume(mesh["V"], mesh["F"]))


def run_case(name, A, B, op, expect, exact_vol=None):
    try:
        m, r = boolean(A, B, op, proxy_tol=1e-3)
    except AmbiguousResult as e:
        n_amb = len(e.report.get("ambiguities", []))
        findings = e.report["stages"]["degeneracy"]["findings"]
        if expect == "ambiguous" and n_amb > 0:
            print(f"[PASS] {name}: AmbiguousResult as expected "
                  f"({n_amb} blocking, findings={findings})")
            return True
        print(f"[FAIL] {name}: unexpected AmbiguousResult "
              f"(expected {expect}, findings={findings})")
        return False
    except Exception as e:
        print(f"[FAIL] {name}: unexpected {type(e).__name__}: {e}")
        return False
    if expect != "accepted" or not r["accepted"]:
        print(f"[FAIL] {name}: accepted={r['accepted']} but expected {expect}")
        return False
    detail = f"vol={vol(m):.6f}"
    ok = True
    if exact_vol is not None:
        ok = abs(vol(m) - exact_vol) < 1e-9
        detail += f" exact={exact_vol}"
    n_amb = len(r.get("ambiguities", []))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: accepted, {detail}, "
          f"ambiguities={n_amb}")
    return ok


def main():
    ok = True
    # Face-coincident union: the true answer is a clean box; every kept
    # face verifies, so this is accepted WITH the degeneracy noted.
    ok &= run_case("coincident faces (union)",
                   Box([0, 0, 0], [1, 1, 1]), Box([1, 0, 0], [2, 1, 1]),
                   "union", "accepted", exact_vol=2.0)
    ok &= run_case("coincident faces (intersection)",
                   Box([0, 0, 0], [1, 1, 1]), Box([1, 0, 0], [2, 1, 1]),
                   "intersection", "accepted", exact_vol=0.0)
    ok &= run_case("identical boxes (difference)",
                   Box([0, 0, 0], [1, 1, 1]), Box([0, 0, 0], [1, 1, 1]),
                   "difference", "accepted", exact_vol=0.0)
    ok &= run_case("tangent cylinders (union)",
                   Cylinder([0, 0, 0], [0, 0, 1], 1.0, 2.0),
                   Cylinder([2, 0, 0], [0, 0, 1], 1.0, 2.0),
                   "union", "ambiguous")
    ok &= run_case("point-touching spheres (union)",
                   Sphere([0, 0, 0], 1.0), Sphere([2, 0, 0], 1.0),
                   "union", "ambiguous")
    ok &= run_case("corner-touching boxes (union)",
                   Box([0, 0, 0], [1, 1, 1]), Box([1, 1, 1], [2, 2, 2]),
                   "union", "ambiguous")
    # Coplanar partial overlap: the notch walls all verify exactly,
    # so this is accepted with the coplanar finding noted.
    ok &= run_case("coplanar partial overlap (difference)",
                   Box([0, 0, 0], [2, 2, 2]), Box([1, 0, 0], [3, 1, 1]),
                   "difference", "accepted", exact_vol=7.0)
    ok &= run_case("nested boxes (difference)",
                   Box([0, 0, 0], [4, 4, 4]), Box([1, 1, 1], [2, 2, 2]),
                   "difference", "accepted", exact_vol=63.0)
    print("\nALL PASS (no silent failures)" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
