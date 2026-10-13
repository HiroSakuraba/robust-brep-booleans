"""Deterministic fuzz: every accepted brepkernel result is re-checked by SMLib.

Random rigidly placed primitives (box, cylinder, sphere, cone, torus), a
random Boolean, brepkernel first; each accepted result goes through
brepkernel.usd.smlib_backend.cross_check. Needs usd-brep.

    python tools/review_probes/smlib_crosscheck_fuzz.py --trials 30 --seed 7 \
        --out /tmp/smlib_fuzz.json

Exit 1 on any "disagree" or crash. "inconclusive" (SMLib could not
answer) is counted and listed but does not fail the run; pass
--fail-inconclusive to make it fail too.
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
import time
import traceback

import numpy as np


def _second_opinion(sb, a, b, op):
    """SMLib's volume and plain OCCT's for a case brepkernel refused."""
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
    out = {}
    try:
        out["smlib_volume"] = sb._volume(sb.smlib_boolean(a, b, op))
    except Exception as exc:
        out["smlib"] = f"{getattr(exc, 'kind', type(exc).__name__)}: {exc}"[:200]
    try:
        algo = {"union": BRepAlgoAPI_Fuse, "intersection": BRepAlgoAPI_Common,
                "difference": BRepAlgoAPI_Cut}[op](a, b)
        out["occt_bop_volume"] = sb._volume(algo.Shape()) if algo.IsDone() else None
    except Exception as exc:
        out["occt_bop"] = f"{type(exc).__name__}: {exc}"[:200]
    if out.get("smlib_volume") is not None and out.get("occt_bop_volume"):
        v, w = out["smlib_volume"], out["occt_bop_volume"]
        out["relative_difference"] = abs(v - w) / max(abs(v), abs(w), 1e-300)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=30)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", default=None)
    ap.add_argument("--fail-inconclusive", action="store_true")
    ap.add_argument("--second-opinion", action="store_true",
                    help="on a brepkernel refusal, record SMLib's answer and plain "
                         "OCCT's (BRepAlgoAPI) volume; informational only")
    ap.add_argument("--dump-dir", default=None,
                    help="write the operands of every non-agreeing case as .brep files")
    args = ap.parse_args()

    from brepkernel import BRepAmbiguousResult, boolean_brep
    from brepkernel.usd import smlib_backend as sb
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
    from OCP.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCone,
                                 BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeSphere,
                                 BRepPrimAPI_MakeTorus)
    from OCP.gp import gp_Ax1, gp_Ax2, gp_Dir, gp_Pnt, gp_Trsf, gp_Vec

    if not sb.available():
        print("usd-brep is not installed", file=sys.stderr)
        return 2
    rng = np.random.default_rng(args.seed)

    def place(s):
        t = gp_Trsf()
        ax = rng.normal(size=3)
        ax /= np.linalg.norm(ax)
        t.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(*map(float, ax))),
                      float(rng.uniform(0, 6.283)))
        s = BRepBuilderAPI_Transform(s, t, True).Shape()
        t2 = gp_Trsf()
        t2.SetTranslation(gp_Vec(*map(float, rng.uniform(-.6, .6, 3))))
        return BRepBuilderAPI_Transform(s, t2, True).Shape()

    def prim(k):
        if k == 0:
            d = rng.uniform(.6, 1.8, 3)
            return place(BRepPrimAPI_MakeBox(gp_Pnt(*(-d / 2)), gp_Pnt(*(d / 2))).Shape())
        if k == 1:
            r, h = rng.uniform(.3, .9), rng.uniform(.8, 2)
            return place(BRepPrimAPI_MakeCylinder(
                gp_Ax2(gp_Pnt(0, 0, -h / 2), gp_Dir(0, 0, 1)), r, h).Shape())
        if k == 2:
            return place(BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), rng.uniform(.4, 1)).Shape())
        if k == 3:
            h = rng.uniform(.8, 1.8)
            return place(BRepPrimAPI_MakeCone(gp_Ax2(gp_Pnt(0, 0, -h / 2), gp_Dir(0, 0, 1)),
                                              rng.uniform(.4, .9), rng.uniform(0, .3), h).Shape())
        return place(BRepPrimAPI_MakeTorus(rng.uniform(.6, 1.0), rng.uniform(.15, .35)).Shape())

    kinds = ["box", "cyl", "sph", "cone", "torus"]
    ops = ["union", "intersection", "difference"]
    tally = collections.Counter()
    records = []
    for i in range(args.trials):
        ka, kb = int(rng.integers(5)), int(rng.integers(5))
        a, b = prim(ka), prim(kb)
        op = ops[int(rng.integers(3))]
        rec = {"trial": i, "a": kinds[ka], "b": kinds[kb], "op": op}
        try:
            out, _rep = boolean_brep(a, b, op)
        except BRepAmbiguousResult as exc:
            rec["brepkernel"] = "refused"
            rec["kind"] = exc.report.get("refusal", {}).get("kind")
            tally["brepkernel_refused"] += 1
            msg = ""
            if args.second_opinion:
                rec["second_opinion"] = _second_opinion(sb, a, b, op)
                msg = " | " + json.dumps(rec["second_opinion"])
            records.append(rec)
            print(f"{i:03d} {rec['a']}-{rec['b']} {op}: brepkernel refused {rec['kind']}{msg}")
            continue
        except Exception as exc:
            rec["brepkernel"] = "crash"
            rec["error"] = f"{type(exc).__name__}: {exc}"
            rec["traceback"] = traceback.format_exc()
            tally["crash"] += 1
            records.append(rec)
            print(f"{i:03d} CRASH {rec['error']}")
            continue
        t0 = time.perf_counter()
        try:
            cc = sb.cross_check(a, b, op, out)
        except Exception as exc:
            cc = {"status": "crash", "error": f"{type(exc).__name__}: {exc}",
                  "traceback": traceback.format_exc()}
        rec["brepkernel"] = "accepted"
        rec["cross_check"] = cc
        rec["seconds"] = round(time.perf_counter() - t0, 3)
        tally[cc["status"]] += 1
        if args.dump_dir and cc["status"] != "agree":
            import os
            from OCP.BRepTools import BRepTools
            os.makedirs(args.dump_dir, exist_ok=True)
            for tag, shp in (("A", a), ("B", b), ("result", out)):
                path = os.path.join(args.dump_dir, f"trial{i:03d}_{op}_{tag}.brep")
                BRepTools.Write_s(shp, path)
            rec["dumped_to"] = args.dump_dir
        records.append(rec)
        line = f"{i:03d} {rec['a']}-{rec['b']} {op}: {cc['status']} ({rec['seconds']} s)"
        if cc["status"] not in ("agree",):
            line += " " + json.dumps({k: cc.get(k) for k in ("kind", "message", "error")
                                      if cc.get(k)})
            if cc["status"] == "disagree":
                c = cc["comparison"]
                line += (f" dV={c['volume']['difference']:.3g} "
                         f"pts={c['points']['disagreements']} solids={c['solids']}")
        print(line, flush=True)
    print("TALLY", dict(tally))
    if args.out:
        with open(args.out, "w") as fh:
            json.dump({"seed": args.seed, "trials": args.trials, "tally": dict(tally),
                       "smlib_version": sb.smlib_version(), "records": records},
                      fh, indent=1, default=str)
    bad = tally["disagree"] + tally["crash"]
    if args.fail_inconclusive:
        bad += tally["inconclusive"]
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
