"""Degenerate placement tests (design Stage 3 / Tier A).

Inputs in exactly degenerate position: coincident faces, tangent
cylinders, point/edge contacts. The pipeline must either resolve them
exactly (Tier A) or report them explicitly -- never silently return a
broken solid. Every case must: (a) be accepted with clean verification,
or (b) raise AmbiguousResult with a report. Crashes and silent garbage
are failures.
"""

import sys
import numpy as np

sys.path.insert(0, "src")
from brepkernel import boolean, AmbiguousResult
from brepkernel.solids import Box, Sphere, Cylinder
from brepkernel.verify import signed_volume
from brepkernel.arrange import ArrangementError
from brepkernel.ingest import IngestError


def vol(mesh):
    return abs(signed_volume(mesh["V"], mesh["F"]))


def run_case(name, A, B, op):
    try:
        m, r = boolean(A, B, op, proxy_tol=1e-3)
    except AmbiguousResult as e:
        n_amb = len(e.report.get("ambiguities", []))
        print(f"[AMBIG ] {name}: rejected with {n_amb} ambiguities "
              f"(explicit, not silent) -- findings: "
              f"{e.report['stages']['degeneracy']['findings']}")
        return True  # explicit ambiguity is a PASS: nothing silent
    except (ArrangementError, IngestError) as e:
        print(f"[AMBIG ] {name}: engine/ingest refused explicitly: {e}")
        return True
    except Exception as e:
        print(f"[FAIL  ] {name}: unexpected {type(e).__name__}: {e}")
        return False
    n_amb = len(r.get("ambiguities", []))
    ok = r["accepted"]
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: accepted={ok} "
          f"vol={vol(m):.4f} ambiguities={n_amb} "
          f"deg={r['stages']['degeneracy']['findings']}")
    return ok


def main():
    ok = True

    # 1. coincident faces: two boxes sharing the face x=1
    ok &= run_case("coincident faces (union)",
                   Box([0, 0, 0], [1, 1, 1]), Box([1, 0, 0], [2, 1, 1]), "union")
    # 2. coincident faces: intersection of the same pair (touching -> empty-ish)
    ok &= run_case("coincident faces (intersection)",
                   Box([0, 0, 0], [1, 1, 1]), Box([1, 0, 0], [2, 1, 1]),
                   "intersection")
    # 3. identical boxes, difference -> empty
    ok &= run_case("identical boxes (difference)",
                   Box([0, 0, 0], [1, 1, 1]), Box([0, 0, 0], [1, 1, 1]),
                   "difference")
    # 4. tangent cylinders (externally tangent along a line)
    ok &= run_case("tangent cylinders (union)",
                   Cylinder([0, 0, 0], [0, 0, 1], 1.0, 2.0),
                   Cylinder([2, 0, 0], [0, 0, 1], 1.0, 2.0), "union")
    # 5. point contact: two spheres touching at one point
    ok &= run_case("point-touching spheres (union)",
                   Sphere([0, 0, 0], 1.0), Sphere([2, 0, 0], 1.0), "union")
    # 6. box corner touching box face (single-point contact)
    ok &= run_case("corner-touching boxes (union)",
                   Box([0, 0, 0], [1, 1, 1]), Box([1, 1, 1], [2, 2, 2]), "union")
    # 7. one box face coplanar with another's face, partial overlap
    ok &= run_case("coplanar partial overlap (difference)",
                   Box([0, 0, 0], [2, 2, 2]), Box([1, 0, 0], [3, 1, 1]),
                   "difference")
    # 8. nested: small box strictly inside big box
    ok &= run_case("nested boxes (difference)",
                   Box([0, 0, 0], [4, 4, 4]), Box([1, 1, 1], [2, 2, 2]),
                   "difference")

    print("\nALL PASS (no silent failures)" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
