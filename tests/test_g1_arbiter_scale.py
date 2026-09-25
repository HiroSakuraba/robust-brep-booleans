"""G1 rework: the independent Boolean arbiter must be scale-aware and must
never mutate its inputs.

Weaknesses demonstrated on the pre-rework arbiter (each test FAILS there):

t1: membership_audit() samples from the hard-coded cube [-2.2, 2.2]^3 no
    matter where the geometry lives. Two overlapping boxes translated
    +1000 in x audit "successfully" (checked == n, zero errors) without a
    single sample point anywhere near the geometry: a vacuous audit.
    The reworked arbiter derives the sample domain from the combined
    bounding box of A, B and out, with a margin.
t2: the absolute surface-exclusion band (2e-3) excludes every interior
    point of a millimetre-scale part, so the audit never examines the
    inside of small parts (interior_checked == 0). The reworked arbiter
    derives the band from the mesh deflection plus the entities'
    tolerances.
t3: triangles() calls BRepTools.Clean_s + BRepMesh_IncrementalMesh on the
    caller's shape, attaching triangulations to the authoritative input
    shape. The reworked triangles() meshes a deep copy and leaves the
    input's triangulation state untouched.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools",
                                "review_probes"))

import numpy as np

import arbiter

from OCP.BRep import BRep_Tool
from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.TopAbs import TopAbs_FACE
from OCP.TopExp import TopExp_Explorer
from OCP.TopLoc import TopLoc_Location
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Pnt


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def box(pmin, pmax):
    return BRepPrimAPI_MakeBox(gp_Pnt(*pmin), gp_Pnt(*pmax)).Shape()


def fuse(a, b):
    return BRepAlgoAPI_Fuse(a, b).Shape()


LEGACY_CUBE = ([-2.2, -2.2, -2.2], [2.2, 2.2, 2.2])


def t1_translated_operands_get_derived_domain():
    """Audit domain must follow the geometry, not the hard-coded cube."""
    a = box((1000.0, 0.0, 0.0), (1001.0, 1.0, 1.0))
    b = box((1000.5, 0.0, 0.0), (1001.5, 1.0, 1.0))
    out = fuse(a, b)
    rng = np.random.default_rng(11)
    res = arbiter.membership_audit(a, b, out, "union", rng, n=120)
    lo, hi = res.get("domain", LEGACY_CUBE)
    covers = lo[0] <= 1000.0 and hi[0] >= 1001.5
    ok = check(
        "g1r translated audit domain covers operand geometry",
        covers,
        f"domain x=[{lo[0]:.3f}, {hi[0]:.3f}] vs geometry x=[1000, 1001.5] "
        f"checked={res['checked']} kernel_errors={res['kernel_errors']}")
    ok &= check(
        "g1r translated audit still decides points and finds no errors",
        res["kernel_errors"] == 0 and res["checked"] >= 100,
        f"checked={res['checked']} kernel_errors={res['kernel_errors']}")
    return ok


def t2_micro_part_gets_derived_band():
    """The exclusion band must scale with deflection/tolerances, and the
    audit must examine interior points of small parts."""
    s = 2e-3
    a = box((0.0, 0.0, 0.0), (s, s, s))
    b = box((s / 2, s / 2, s / 2), (3 * s / 2, 3 * s / 2, 3 * s / 2))
    out = fuse(a, b)
    rng = np.random.default_rng(12)
    res = arbiter.membership_audit(a, b, out, "union", rng, n=300)
    band = res.get("band_used", 2e-3)
    interior = res.get("interior_checked", 0)
    ok = check(
        "g1r micro-part exclusion band is derived, not the absolute 2e-3",
        band < 1e-3,
        f"band_used={band:.3e} (absolute 2e-3 excludes every interior "
        f"point of a {s} part)")
    ok &= check(
        "g1r micro-part audit examines interior points",
        interior > 0,
        f"interior_checked={interior} checked={res['checked']}")
    ok &= check(
        "g1r micro-part audit finds no errors on the correct fuse",
        res["kernel_errors"] == 0,
        f"kernel_errors={res['kernel_errors']} checked={res['checked']}")
    return ok


def triangulated_faces(shape):
    n = 0
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        face = TopoDS.Face(ex.Current())
        loc = TopLoc_Location()
        if BRep_Tool.Triangulation_s(face, loc) is not None:
            n += 1
        ex.Next()
    return n


def t3_triangles_does_not_mutate_input():
    """triangles() must not attach tessellations to the caller's shape."""
    a = box((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    before = triangulated_faces(a)
    soup = arbiter.triangles(a)
    after = triangulated_faces(a)
    ok = check(
        "g1r triangles() leaves input triangulation state unchanged",
        before == 0 and after == 0,
        f"triangulated faces before={before} after={after} "
        f"soup shape={soup.shape}")
    return ok


def main():
    ok = True
    ok &= t1_translated_operands_get_derived_domain()
    ok &= t2_micro_part_gets_derived_band()
    ok &= t3_triangles_does_not_mutate_input()
    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
