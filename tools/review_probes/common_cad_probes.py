"""Everyday CAD operations that exercise coincident (flush) faces.

Review baseline (main at 88d09eb, 2026-09-25): 17 cases, 4 accepted and
13 refused. Eleven refusals are coincident-face cases (finding F1); the
equal-radius crossing cylinders are a genuine singular intersection and
SHOULD keep refusing.

Every "expect" field is the target after Gate G2. Exit status is 1 if any
case is accepted with a wrong volume (must never happen), if any case
marked expect="accept" is refused (G2 not yet met), or if any case marked
expect="refuse" is accepted (needs human inspection: the kernel may have
improved, or it may be wrongly accepting). Use --baseline to only fail on
wrong answers and unexpected accepts, not on refusals of accept-expected
cases.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from brepkernel import boolean_brep, BRepAmbiguousResult  # noqa: E402

from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder,
                             BRepPrimAPI_MakeSphere)
from OCP.GProp import GProp_GProps
from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

PI = math.pi


def vol(s):
    g = GProp_GProps()
    BRepGProp.VolumePropertiesGK_s(s, g, 1e-10, True, True, False, False, False)
    return g.Mass()


def box(x0, y0, z0, x1, y1, z1):
    return BRepPrimAPI_MakeBox(gp_Pnt(x0, y0, z0), gp_Pnt(x1, y1, z1)).Shape()


def cyl(p, d, r, h):
    return BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(*p), gp_Dir(*d)), r, h).Shape()


A = box(0, 0, 0, 2, 2, 2)
CASES = [
    # name, A, B, op, exact volume, expected disposition after G2
    ("union: boxes sharing a full face", box(0, 0, 0, 1, 1, 1), box(1, 0, 0, 2, 1, 1), "union", 2.0, "accept"),
    ("union: boxes sharing a partial face", box(0, 0, 0, 1, 1, 1), box(1, .25, .25, 2, .75, .75), "union", 1.25, "accept"),
    ("union: offset box stacked on top", A, box(.5, .5, 2, 1.5, 1.5, 3), "union", 9.0, "accept"),
    ("cut: slot flush with top face", A, box(.5, -1, 1, 1.5, 3, 2), "difference", 6.0, "accept"),
    ("cut: pocket flush with top face", A, box(.5, .5, 1, 1.5, 1.5, 2), "difference", 7.0, "accept"),
    ("cut: slot overshooting top face", A, box(.5, -1, 1, 1.5, 3, 3), "difference", 6.0, "accept"),
    ("cut: through-hole cylinder", A, cyl((1, 1, -1), (0, 0, 1), .5, 4), "difference", 8 - PI * .25 * 2, "accept"),
    ("cut: blind hole flush with top", A, cyl((1, 1, 1), (0, 0, 1), .5, 1), "difference", 8 - PI * .25, "accept"),
    ("union: boss on top face", A, cyl((1, 1, 2), (0, 0, 1), .5, 1), "union", 8 + PI * .25, "accept"),
    ("common: sub-box sharing 3 faces", A, box(0, 0, 0, 1, 1, 1), "intersection", 1.0, "accept"),
    ("cut: corner notch sharing 3 faces", A, box(1, 1, 1, 2, 2, 2), "difference", 7.0, "accept"),
    ("union: crossing cylinders r1 != r2", cyl((0, 0, -2), (0, 0, 1), 1, 4), cyl((-2, 0, 0), (1, 0, 0), .6, 4), "union", None, "accept"),
    ("union: crossing cylinders equal r (singular)", cyl((0, 0, -2), (0, 0, 1), 1, 4), cyl((-2, 0, 0), (1, 0, 0), 1, 4), "union", None, "refuse"),
    ("cut: sphere from box corner", A, BRepPrimAPI_MakeSphere(gp_Pnt(2.1, 1.3, .9), .8).Shape(), "difference", None, "accept"),
    # Touching-only cases: regularized result is lower-dimensional or
    # non-manifold. Accept only if the kernel supports non-manifold output;
    # otherwise a typed refusal is correct.
    ("union: boxes touching along one edge", box(0, 0, 0, 1, 1, 1), box(1, 1, 0, 2, 2, 1), "union", 2.0, "either"),
    ("common: boxes sharing only a face", box(0, 0, 0, 1, 1, 1), box(1, 0, 0, 2, 1, 1), "intersection", 0.0, "accept"),
    ("cut: tool sharing only a face", box(0, 0, 0, 1, 1, 1), box(1, 0, 0, 2, 1, 1), "difference", 1.0, "accept"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline", action="store_true",
                    help="fail only on wrong answers and unexpected accepts, "
                         "not on refusals of accept-expected cases")
    args = ap.parse_args()
    wrong = unmet = unexpected_accept = 0
    from brepkernel.assembly import _shape_volume  # noqa: F401  (import check)
    from arbiter import version_banner  # noqa: E402
    print(version_banner(), flush=True)
    for name, a, b, op, exact, expect in CASES:
        t = time.perf_counter()
        from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
        occ = {"union": BRepAlgoAPI_Fuse, "difference": BRepAlgoAPI_Cut,
               "intersection": BRepAlgoAPI_Common}[op](a, b).Shape()
        ref = exact if exact is not None else vol(occ)
        try:
            out, rep = boolean_brep(a, b, op)
            v = 0.0 if rep["stages"].get("assembly", {}).get("empty") else vol(out)
            ok = abs(v - ref) <= 1e-7 * max(1.0, abs(ref))
            status = "ACCEPT" if ok else "WRONG"
            detail = f"vol={v:.10g} ref={ref:.10g}"
            wrong += not ok
            if expect == "refuse":
                detail += "  (accepted a case expected to refuse: inspect)"
                unexpected_accept += 1
        except BRepAmbiguousResult as e:
            r = e.report.get("refusal", {})
            status = "REFUSE"
            detail = f"{r.get('stage')}/{r.get('kind')}"
            unmet += expect == "accept"
        print(f"{status:7s} expect={expect:6s} {name:46s} {time.perf_counter()-t:6.2f}s  {detail}", flush=True)
    print(f"wrong={wrong} unmet_accepts={unmet} "
          f"unexpected_accepts={unexpected_accept}")
    if wrong or unexpected_accept:
        return 1
    return 0 if (args.baseline or unmet == 0) else 1


if __name__ == "__main__":
    raise SystemExit(main())
