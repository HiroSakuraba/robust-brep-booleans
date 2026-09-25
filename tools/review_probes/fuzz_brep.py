"""Randomized generic-position fuzz of boolean_brep() with an independent arbiter.

Usage (from repo root):
    python tools/review_probes/fuzz_brep.py --trials 150 --seed 7
    python tools/review_probes/fuzz_brep.py --trials 60 --seed 11 --nurbs
    python tools/review_probes/fuzz_brep.py --trials 60 --seed 3 --snap 0.25
    python tools/review_probes/fuzz_brep.py --trials 60 --seed 5 --snap 0.5 --kinds box

--snap G rounds every translation and dimension to a multiple of G and uses
only axis-aligned rotations by multiples of 90 degrees. That generates the
coincident-face configurations that dominate real CAD (review finding F1).

Outcome classes per trial:
  accept        result accepted, OCCT volume agrees, winding arbiter agrees
  refuse:<kind> typed BRepAmbiguousResult (allowed; counted toward refusal rate)
  WRONG         accepted but the winding-number arbiter found a wrong point
                or the volume disagrees with OCCT. This must never happen.
  CRASH         any other exception. Stop and report (project rule).

Exit status: 1 if any WRONG or CRASH, else 0. A JSON summary is written to
--out (default fuzz_summary.json).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
sys.path.insert(0, os.path.dirname(__file__))

from brepkernel import boolean_brep, BRepAmbiguousResult  # noqa: E402
from arbiter import has_solid, membership_audit, version_banner  # noqa: E402

from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCP.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert, BRepBuilderAPI_Transform
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCone,
                             BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeSphere,
                             BRepPrimAPI_MakeTorus)
from OCP.GProp import GProp_GProps
from OCP.gp import gp_Ax1, gp_Ax2, gp_Dir, gp_Pnt, gp_Trsf, gp_Vec

OCC = {"union": BRepAlgoAPI_Fuse, "difference": BRepAlgoAPI_Cut,
       "intersection": BRepAlgoAPI_Common}
KINDS = ["box", "cyl", "sph", "cone", "torus"]


def vol(shape) -> float:
    if not has_solid(shape):
        return 0.0
    g = GProp_GProps()
    BRepGProp.VolumePropertiesGK_s(shape, g, 1e-10, True, True, False, False, False)
    return float(g.Mass())


def q(x, snap):
    return float(np.round(x / snap) * snap) if snap else float(x)


def place(shape, rng, snap):
    t = gp_Trsf()
    if snap:
        axis = [(1, 0, 0), (0, 1, 0), (0, 0, 1)][rng.integers(3)]
        ang = (math.pi / 2) * int(rng.integers(4))
    else:
        axis = rng.normal(size=3)
        axis = axis / np.linalg.norm(axis)
        ang = rng.uniform(0, 2 * math.pi)
    t.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(*map(float, axis))), ang)
    t2 = gp_Trsf()
    t2.SetTranslation(gp_Vec(*[q(v, snap) for v in rng.uniform(-0.6, 0.6, 3)]))
    s = BRepBuilderAPI_Transform(shape, t, True).Shape()
    return BRepBuilderAPI_Transform(s, t2, True).Shape()


def prim(kind, rng, snap):
    u = lambda lo, hi: max(q(rng.uniform(lo, hi), snap), snap or 0.0) or rng.uniform(lo, hi)
    if kind == "box":
        d = np.array([u(0.6, 1.8) for _ in range(3)])
        s = BRepPrimAPI_MakeBox(gp_Pnt(*(-d / 2)), gp_Pnt(*(d / 2))).Shape()
    elif kind == "cyl":
        r, h = u(0.3, 0.9), u(0.8, 2.0)
        s = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -h / 2), gp_Dir(0, 0, 1)), r, h).Shape()
    elif kind == "sph":
        s = BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), u(0.4, 1.0)).Shape()
    elif kind == "cone":
        h = u(0.8, 1.8)
        s = BRepPrimAPI_MakeCone(gp_Ax2(gp_Pnt(0, 0, -h / 2), gp_Dir(0, 0, 1)),
                                 u(0.4, 0.9), rng.uniform(0.0, 0.3), h).Shape()
    else:
        s = BRepPrimAPI_MakeTorus(u(0.6, 1.0), rng.uniform(0.15, 0.35)).Shape()
    return place(s, rng, snap)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=50)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--nurbs", action="store_true")
    ap.add_argument("--snap", type=float, default=0.0)
    ap.add_argument("--points", type=int, default=300)
    ap.add_argument("--kinds", default=",".join(KINDS),
                    help="comma list from box,cyl,sph,cone,torus")
    ap.add_argument("--out", default="fuzz_summary.json")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    kinds = [k for k in args.kinds.split(",") if k in KINDS]
    print(version_banner(), flush=True)
    rows, tally = [], {}
    for i in range(args.trials):
        ka, kb = str(rng.choice(kinds)), str(rng.choice(kinds))
        op = str(rng.choice(["union", "intersection", "difference"]))
        a, b = prim(ka, rng, args.snap), prim(kb, rng, args.snap)
        if args.nurbs:
            a = BRepBuilderAPI_NurbsConvert(a, True).Shape()
            b = BRepBuilderAPI_NurbsConvert(b, True).Shape()
        # Per-trial arbiter RNG so refusals do not shift later trials.
        arng = np.random.default_rng([args.seed, i])
        t0 = time.perf_counter()
        row = {"trial": i, "A": ka, "B": kb, "op": op}
        try:
            oracle = vol(OCC[op](a, b).Shape())
        except Exception as exc:  # oracle failure is data, not a kernel verdict
            oracle = float("nan")
            row["oracle_error"] = repr(exc)
        try:
            out, rep = boolean_brep(a, b, op)
            v = vol(out)
            audit = membership_audit(a, b, out, op, arng, n=args.points)
            vol_ok = (not math.isfinite(oracle)) or abs(v - oracle) <= 1e-6 * max(1.0, abs(oracle))
            status = "accept" if (audit["kernel_errors"] == 0 and vol_ok) else "WRONG"
            row.update(volume=v, oracle_volume=oracle, audit=audit)
        except BRepAmbiguousResult as exc:
            r = exc.report.get("refusal", {})
            status = f"refuse:{r.get('kind')}"
            row["refusal"] = r
        except Exception as exc:
            status = "CRASH"
            row["error"] = f"{type(exc).__name__}: {exc}"
        row["status"] = status
        row["seconds"] = round(time.perf_counter() - t0, 3)
        rows.append(row)
        tally[status] = tally.get(status, 0) + 1
        print(f"{i:4d} {ka:5s} {op:12s} {kb:5s} {status:36s} {row['seconds']:7.2f}s", flush=True)
        if status == "CRASH":
            print("CRASH: stopping per project rule:", row["error"], flush=True)
            break

    summary = {"args": vars(args), "versions": version_banner(),
               "tally": tally, "rows": rows}
    with open(args.out, "w") as fh:
        json.dump(summary, fh, indent=1, default=str)
    print("TALLY", tally)
    return 1 if any(k in ("WRONG", "CRASH") for k in tally) else 0


if __name__ == "__main__":
    raise SystemExit(main())
