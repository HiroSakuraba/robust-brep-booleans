"""Deterministic fuzz: SMLib re-checks every accept; refusals get triaged.

Random rigidly placed primitives (box, cylinder, sphere, cone, torus), a
random Boolean, brepkernel first; each accepted result goes through
brepkernel.usd.smlib_backend.cross_check. Needs usd-brep.

    python tools/review_probes/smlib_crosscheck_fuzz.py --trials 30 --seed 7 \
        --out /tmp/smlib_fuzz.json

--second-opinion: every brepkernel refusal also gets OCCT's own Boolean
and SMLib's answer (brepkernel.advisory) and a triage label:
suspected_over_refusal (the two agree to 1e-6: a regression-queue
candidate), justified (they differ by 1e-4 or more) or unclear. The
METRICS line reports refusal_rate and suspected_over_refusal_rate.
--dump-dir writes the operands of refused and non-agreeing cases as .brep
files so each queue entry can be reproduced.

Exit 1 on any "disagree" or crash. "inconclusive" (SMLib could not
answer) is counted and listed but does not fail the run; pass
--fail-inconclusive to make it fail too. Refusals never fail the run.
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
    """SMLib's volume and plain OCCT's for a case brepkernel refused
    (brepkernel.advisory, flattened for the ledger)."""
    from brepkernel.advisory import second_opinions
    summary, _shapes = second_opinions(a, b, op)
    out = {}
    for name, key in (("smlib", "smlib"), ("occt", "occt_bop")):
        rec = summary["engines"].get(name, {})
        if "volume" in rec:
            out[f"{key}_volume"] = rec["volume"]
        else:
            out[key] = rec.get("error") or rec.get("status")
    out["relative_difference"] = summary["relative_difference"]
    out["triage"] = summary["triage"]
    return out


def _dump(directory, i, op, shapes):
    import os
    from OCP.BRepTools import BRepTools
    os.makedirs(directory, exist_ok=True)
    for tag, shp in shapes:
        BRepTools.Write_s(shp, os.path.join(directory, f"trial{i:03d}_{op}_{tag}.brep"))


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
                    help="write the operands of every refused or non-agreeing case "
                         "as .brep files")
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
                tally["refusal_" + rec["second_opinion"]["triage"]] += 1
                msg = " | " + json.dumps(rec["second_opinion"])
            if args.dump_dir:
                _dump(args.dump_dir, i, op, (("A", a), ("B", b)))
                rec["dumped_to"] = args.dump_dir
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
            _dump(args.dump_dir, i, op, (("A", a), ("B", b), ("result", out)))
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
    n = args.trials
    metrics = {
        "trials": n,
        "refusal_rate": tally["brepkernel_refused"] / n if n else 0.0,
    }
    if args.second_opinion:
        metrics["suspected_over_refusal_rate"] = (
            tally["refusal_suspected_over_refusal"] / n if n else 0.0)
    print("METRICS", json.dumps(metrics))
    if args.out:
        with open(args.out, "w") as fh:
            json.dump({"seed": args.seed, "trials": args.trials, "tally": dict(tally),
                       "metrics": metrics,
                       "smlib_version": sb.smlib_version(), "records": records},
                      fh, indent=1, default=str)
    bad = tally["disagree"] + tally["crash"]
    if args.fail_inconclusive:
        bad += tally["inconclusive"]
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
