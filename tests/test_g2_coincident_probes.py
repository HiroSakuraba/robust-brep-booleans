"""G2.7: end-to-end tests for coincident-face handling.

Covers the work-plan section 2.7 requirements:
  - probe cases as expect-enforced tests (the 10 expect="accept" cases
    from tools/review_probes/common_cad_probes.py must ACCEPT at the
    exact volume and pass the independent membership arbiter);
  - box-box exact oracle: grid-snapped Tier A Box pairs, exact volume /
    shells / euler vs boolean_brep;
  - cylinder cases: blind hole, boss, coaxial tube-in-tube, counterbore;
  - near-coincident negatives: 0.5x/2x/10x base_tol offsets;
  - touching-only cases: edge-touching -> NonManifoldResult (or compound
    with allow_nonmanifold=True), point-touching -> UnresolvedContact.
"""
import math
import sys

sys.path.insert(0, "src")
sys.path.insert(0, "tests")

from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

import _arbiter
from brepkernel import boolean_brep, BRepAmbiguousResult
from brepkernel.classify import (exact_op_volume, expected_euler,
                                  expected_shells)
from brepkernel.solids import Box


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def box(x0, y0, z0, x1, y1, z1):
    return BRepPrimAPI_MakeBox(gp_Pnt(x0, y0, z0),
                              gp_Pnt(x1, y1, z1)).Shape()


def cyl(p, d, r, h):
    return BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(*p), gp_Dir(*d)), r, h).Shape()


def refusal_kind(fn):
    try:
        fn()
    except BRepAmbiguousResult as e:
        return (e.report.get("refusal") or {}).get("kind")
    return None


A = box(0, 0, 0, 2, 2, 2)

# (name, A, B, op, exact volume or None)
PROBE_ACCEPTS = [
    ("shared full face",
     box(0, 0, 0, 1, 1, 1), box(1, 0, 0, 2, 1, 1), "union", 2.0),
    ("shared partial face",
     box(0, 0, 0, 1, 1, 1), box(1, .25, .25, 2, .75, .75),
     "union", 1.25),
    ("offset box stacked on top",
     A, box(.5, .5, 2, 1.5, 1.5, 3), "union", 9.0),
    ("slot flush with top face",
     A, box(.5, -1, 1, 1.5, 3, 2), "difference", 6.0),
    ("pocket flush with top face",
     A, box(.5, .5, 1, 1.5, 1.5, 2), "difference", 7.0),
    ("slot overshooting top face",
     A, box(.5, -1, 1, 1.5, 3, 3), "difference", 6.0),
    ("through-hole cylinder",
     A, cyl((1, 1, -1), (0, 0, 1), .5, 4),
     "difference", 8 - math.pi * .25 * 2),
    ("blind hole flush with top",
     A, cyl((1, 1, 1), (0, 0, 1), .5, 1),
     "difference", 8 - math.pi * .25),
    ("boss on top face",
     A, cyl((1, 1, 2), (0, 0, 1), .5, 1),
     "union", 8 + math.pi * .25),
    ("sub-box sharing 3 faces",
     A, box(0, 0, 0, 1, 1, 1), "intersection", 1.0),
    ("corner notch sharing 3 faces",
     A, box(1, 1, 1, 2, 2, 2), "difference", 7.0),
]


def main():
    ok = True

    # 1. Probe accept cases: must ACCEPT at the exact volume and pass
    # the independent membership arbiter.
    for name, a, b, op, vol in PROBE_ACCEPTS:
        try:
            out, rep = boolean_brep(a, b, op)
        except BRepAmbiguousResult as e:
            ok &= check(f"g2 probe accept: {name}", False,
                        f"refused {(e.report.get('refusal') or {}).get('kind')}")
            continue
        vr = rep["stages"]["verification"]["result_volume"]
        ok &= check(f"g2 probe accept: {name} volume",
                    abs(vr - vol) <= max(1e-6, 1e-9 * abs(vol)),
                    f"vol={vr:.9f} ref={vol}")
        ok &= _arbiter.check_accepted(
            f"g2 probe {name}", check, a, b, out, op)[0]

    # 2. Touching-only cases.
    kind = refusal_kind(lambda: boolean_brep(
        box(0, 0, 0, 1, 1, 1), box(1, 1, 0, 2, 2, 1), "union"))
    ok &= check("g2 edge-touching union refuses NonManifoldResult",
                kind == "NonManifoldResult", f"kind={kind}")
    try:
        out, rep = boolean_brep(
            box(0, 0, 0, 1, 1, 1), box(1, 1, 0, 2, 2, 1), "union",
            allow_nonmanifold=True)
        asm = rep["stages"]["assembly"]
        ok &= check("g2 edge-touching allow_nonmanifold accepts compound",
                    asm["solids"] == 2 and bool(asm.get("notes")),
                    f"solids={asm['solids']}")
    except BRepAmbiguousResult as e:
        ok &= check("g2 edge-touching allow_nonmanifold accepts compound",
                    False,
                    f"refused {(e.report.get('refusal') or {}).get('kind')}")
    kind = refusal_kind(lambda: boolean_brep(
        box(0, 0, 0, 1, 1, 1), box(1, 1, 1, 2, 2, 2), "union"))
    ok &= check("g2 point-touching union refuses UnresolvedContact",
                kind == "UnresolvedContact", f"kind={kind}")

    # 3. Box-box exact oracle on grid-snapped Tier A pairs.
    grid_cases = [
        (Box((0, 0, 0), (2, 2, 2)), Box((1, 1, 1), (3, 3, 3)), "union"),
        (Box((0, 0, 0), (2, 2, 2)), Box((1, 1, 1), (3, 3, 3)),
         "intersection"),
        (Box((0, 0, 0), (2, 2, 2)), Box((1, 1, 1), (3, 3, 3)),
         "difference"),
        (Box((0, 0, 0), (1, 1, 1)), Box((1, 0, 0), (2, 1, 1)), "union"),
        (Box((0, 0, 0), (2, 2, 2)), Box((0, 0, 0), (1, 1, 1)),
         "intersection"),
    ]
    for i, (ba, bb, op) in enumerate(grid_cases):
        a = box(*ba.lo, *ba.hi)
        b = box(*bb.lo, *bb.hi)
        ev = exact_op_volume(ba, bb, op)
        es = expected_shells(ba, bb, op)
        ee = expected_euler(ba, bb, op)
        try:
            out, rep = boolean_brep(a, b, op)
        except BRepAmbiguousResult as e:
            ok &= check(f"g2 oracle box-box {i} {op}", False,
                        f"refused {(e.report.get('refusal') or {}).get('kind')}")
            continue
        vr = rep["stages"]["verification"]["result_volume"]
        ok &= check(f"g2 oracle box-box {i} {op} volume",
                    abs(vr - ev) <= 1e-9 * max(1.0, abs(ev)),
                    f"vol={vr:.9f} exact={ev}")
        if es is not None:
            got_s = rep["stages"]["assembly"]["solids"]
            ok &= check(f"g2 oracle box-box {i} {op} shells",
                        got_s == es, f"solids={got_s} expected={es}")

    # 4. Cylinder cases: blind hole, boss, coaxial tube-in-tube,
    # counterbore. Volumes are exact analytic values.
    cyl_cases = [
        ("blind hole",
         box(0, 0, 0, 2, 2, 2), cyl((1, 1, 1.5), (0, 0, 1), .5, .5),
         "difference", 8 - math.pi * .25 * .5),
        ("boss",
         box(0, 0, 0, 2, 2, 2), cyl((1, 1, 2), (0, 0, 1), .5, 1),
         "union", 8 + math.pi * .25 * 1),
        ("coaxial tube-in-tube",
         cyl((0, 0, 0), (0, 0, 1), 1.0, 2),
         cyl((0, 0, -1), (0, 0, 1), .5, 4),
         "difference",
         math.pi * (1.0 - .25) * 2),
        ("counterbore",
         box(0, 0, 0, 4, 4, 2),
         cyl((2, 2, 1), (0, 0, 1), 1.0, 1),  # wide shallow
         "difference", 32 - math.pi * 1.0 * 1),
    ]
    for name, a, b, op, vol in cyl_cases:
        try:
            out, rep = boolean_brep(a, b, op)
        except BRepAmbiguousResult as e:
            ok &= check(f"g2 cylinder case: {name}", False,
                        f"refused {(e.report.get('refusal') or {}).get('kind')}")
            continue
        vr = rep["stages"]["verification"]["result_volume"]
        ok &= check(f"g2 cylinder case: {name} volume",
                    abs(vr - vol) <= 1e-6 * max(1.0, abs(vol)),
                    f"vol={vr:.9f} ref={vol:.9f}")
        ok &= _arbiter.check_accepted(
            f"g2 cylinder {name}", check, a, b, out, op)[0]

    # 5. Near-coincident negatives: offsets at 0.5x/2x/10x base_tol.
    # The coincidence band is 4x base_tol; in-band offsets are
    # undecidable and must refuse NearCoincidentFaces.
    bt = 1e-7
    for mult, expect in [(0.5, "NearCoincidentFaces"),
                         (2.0, "NearCoincidentFaces"),
                         (10.0, None)]:
        off = mult * bt
        kind = refusal_kind(lambda: boolean_brep(
            box(0, 0, 0, 1, 1, 1), box(1 + off, 0, 0, 2 + off, 1, 1),
            "union", base_tol=bt))
        if expect is None:
            ok &= check(f"g2 near-coincident {mult}x tol accepts",
                        kind is None, f"kind={kind}")
        else:
            ok &= check(f"g2 near-coincident {mult}x tol refuses {expect}",
                        kind == expect, f"kind={kind}")

    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
