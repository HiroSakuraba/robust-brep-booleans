"""Verdict equivalence harness: prove a tree change did not alter semantics.

Compares boolean_brep() verdicts between a --before and an --after source
tree over a deterministic corpus (everyday CAD probes plus fixed-seed fuzz
from tools/review_probes/corpus_manifest.json).

Each side runs in a scrubbed subprocess whose sys.path contains ONLY the
requested src tree (the worker strips PYTHONPATH from its environment and
inserts just the given tree). --verbose prints the imported
brepkernel.__file__ for both sides so a shadowing accident is visible.

Compared per case: verdict (accept/refuse/crash/timeout), refusal kind and
stage, accepted volume (relative tolerance 1e-9), and optionally a cheap
topology summary (solids, shells, faces).

Exit status: 0 when there are no blocking differences, 1 otherwise.
--allow-new-accepts permits refuse->accept transitions (for gates whose
explicit purpose is expanding acceptance). Everything else blocks:
accept->refuse, changed refusal kind, changed accepted volume, any crash
or timeout on either side.

Example:
  python tools/review_probes/verdict_equivalence.py \\
      --before /path/to/base/src --after /path/to/branch/src --verbose
  python tools/review_probes/verdict_equivalence.py \\
      --before ... --after ... --fuzz-trials 200 --fuzz-seed 11 \\
      --allow-new-accepts --out /tmp/equiv.json
"""
from __future__ import annotations

import argparse
import base64
import json
import math
import os
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_MANIFEST = os.path.join(HERE, "corpus_manifest.json")

VOL_REL_TOL = 1e-9
VOL_ABS_FLOOR = 1e-12

# ---------------------------------------------------------------------------
# Worker payload. Runs in a scrubbed subprocess: sys.path gets ONLY the
# requested src tree. Builds every case from parametric specs with OCP
# directly (no imports from any tools/ tree, which could shadow brepkernel).
# ---------------------------------------------------------------------------
WORKER = r'''
import json, math, os, sys, time

SRC = sys.argv[1]
CASES_PATH = sys.argv[2]
VERBOSE = sys.argv[3] == "1"

# Scrub: the parent already removed PYTHONPATH from the environment, but
# drop any leftover path entry that could shadow the requested tree.
sys.path = [p for p in sys.path
            if not (p.endswith("/src") or "brep-gates" in p or "brepkernel" in p)]
sys.path.insert(0, SRC)

import brepkernel
from brepkernel import boolean_brep, BRepAmbiguousResult

if VERBOSE:
    print("WORKER brepkernel.__file__ =", brepkernel.__file__, flush=True)
assert brepkernel.__file__.startswith(SRC), (
    "shadowing accident: %s not under %s" % (brepkernel.__file__, SRC))

from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCone,
                             BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeSphere,
                             BRepPrimAPI_MakeTorus)
from OCP.GProp import GProp_GProps
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS
from OCP.TopAbs import TopAbs_SOLID, TopAbs_SHELL, TopAbs_FACE
from OCP.gp import gp_Ax1, gp_Ax2, gp_Dir, gp_Pnt, gp_Trsf, gp_Vec


def build(spec):
    k = spec["prim"]
    if k == "box":
        s = BRepPrimAPI_MakeBox(gp_Pnt(*spec["p0"]), gp_Pnt(*spec["p1"])).Shape()
    elif k == "cyl":
        s = BRepPrimAPI_MakeCylinder(
            gp_Ax2(gp_Pnt(*spec["p"]), gp_Dir(*spec["d"])),
            spec["r"], spec["h"]).Shape()
    elif k == "sphere":
        s = BRepPrimAPI_MakeSphere(gp_Pnt(*spec["c"]), spec["r"]).Shape()
    elif k == "cone":
        s = BRepPrimAPI_MakeCone(
            gp_Ax2(gp_Pnt(*spec["p"]), gp_Dir(*spec["d"])),
            spec["r1"], spec["r2"], spec["h"]).Shape()
    elif k == "torus":
        s = BRepPrimAPI_MakeTorus(spec["R"], spec["r"]).Shape()
    else:
        raise ValueError("unknown prim %r" % k)
    xf = spec.get("xform")
    if xf:
        t = gp_Trsf()
        ax, ang = xf["rotate"]
        t.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(*ax)), ang)
        s = BRepBuilderAPI_Transform(s, t, True).Shape()
        t2 = gp_Trsf()
        t2.SetTranslation(gp_Vec(*xf["translate"]))
        s = BRepBuilderAPI_Transform(s, t2, True).Shape()
    return s


def volume_of(shape):
    g = GProp_GProps()
    BRepGProp.VolumePropertiesGK_s(shape, g, 1e-10, True, True, False, False, False)
    return float(g.Mass())


def topo_summary(shape):
    def count(typ):
        e = TopExp_Explorer(shape, typ)
        n = 0
        while e.More():
            n += 1
            e.Next()
        return n
    return [count(TopAbs_SOLID), count(TopAbs_SHELL), count(TopAbs_FACE)]


def main():
    with open(CASES_PATH) as fh:
        cases = json.load(fh)
    do_topo = os.environ.get("EQUIV_TOPO") == "1"
    out = []
    for case in cases:
        t0 = time.perf_counter()
        row = {"name": case["name"]}
        try:
            a = build(case["A"])
            b = build(case["B"])
            res, rep = boolean_brep(a, b, case["op"])
            empty = bool(rep.get("stages", {}).get("assembly", {}).get("empty"))
            row["verdict"] = "accept"
            row["volume"] = 0.0 if empty else volume_of(res)
            row["empty"] = empty
            row["topo"] = topo_summary(res) if do_topo else None
        except BRepAmbiguousResult as exc:
            r = (exc.report or {}).get("refusal", {}) or {}
            row["verdict"] = "refuse"
            row["kind"] = r.get("kind")
            row["stage"] = r.get("stage")
        except Exception as exc:  # noqa: BLE001 - crash is data here
            row["verdict"] = "crash"
            row["error"] = "%s: %s" % (type(exc).__name__, str(exc)[:300])
        row["seconds"] = round(time.perf_counter() - t0, 3)
        out.append(row)
    print("EQUIV-RESULTS-BEGIN")
    print(json.dumps(out))
    print("EQUIV-RESULTS-END")


main()
'''


# ---------------------------------------------------------------------------
# Parent: deterministic case expansion (mirrors tools/review_probes/fuzz_brep.py
# geometry, but emits concrete specs so both workers see identical inputs).
# ---------------------------------------------------------------------------
def _q(x, snap):
    return float(np.round(x / snap) * snap) if snap else float(x)


def _prim_spec(kind, rng, snap):
    u = lambda lo, hi: max(_q(rng.uniform(lo, hi), snap), snap or 0.0) or rng.uniform(lo, hi)
    if kind == "box":
        d = [u(0.6, 1.8) for _ in range(3)]
        spec = {"prim": "box", "p0": [-d[0] / 2, -d[1] / 2, -d[2] / 2],
                "p1": [d[0] / 2, d[1] / 2, d[2] / 2]}
    elif kind == "cyl":
        r, h = u(0.3, 0.9), u(0.8, 2.0)
        spec = {"prim": "cyl", "p": [0, 0, -h / 2], "d": [0, 0, 1], "r": r, "h": h}
    elif kind == "sph":
        spec = {"prim": "sphere", "c": [0, 0, 0], "r": u(0.4, 1.0)}
    elif kind == "cone":
        h = u(0.8, 1.8)
        spec = {"prim": "cone", "p": [0, 0, -h / 2], "d": [0, 0, 1],
                "r1": u(0.4, 0.9), "r2": rng.uniform(0.0, 0.3), "h": h}
    else:
        spec = {"prim": "torus", "R": u(0.6, 1.0), "r": rng.uniform(0.15, 0.35)}
    if snap:
        axis = [(1, 0, 0), (0, 1, 0), (0, 0, 1)][int(rng.integers(3))]
        ang = (math.pi / 2) * int(rng.integers(4))
    else:
        axis = rng.normal(size=3)
        axis = axis / np.linalg.norm(axis)
        ang = float(rng.uniform(0, 2 * math.pi))
    tr = [_q(v, snap) for v in rng.uniform(-0.6, 0.6, 3)]
    spec["xform"] = {"rotate": [[float(c) for c in axis], float(ang)],
                     "translate": tr}
    return spec


def _expand_fuzz(cfg):
    rng = np.random.default_rng(cfg["seed"])
    kinds = cfg["kinds"]
    cases = []
    for i in range(cfg["trials"]):
        ka = str(rng.choice(kinds))
        kb = str(rng.choice(kinds))
        op = str(rng.choice(["union", "intersection", "difference"]))
        cases.append({
            "name": "fuzz seed=%s trial=%d (%s %s %s)" % (
                cfg["seed"], i, ka, op, kb),
            "A": _prim_spec(ka, rng, cfg["snap"]),
            "B": _prim_spec(kb, rng, cfg["snap"]),
            "op": op,
        })
    return cases


def run_side(src_dir, cases_path, verbose, timeout, label, do_topo):
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    if do_topo:
        env["EQUIV_TOPO"] = "1"
    cmd = [sys.executable, "-c", WORKER, os.path.abspath(src_dir),
           cases_path, "1" if verbose else "0"]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              timeout=timeout, env=env)
    except subprocess.TimeoutExpired:
        return None, ["%s: worker TIMEOUT after %ss" % (label, timeout)]
    if verbose and proc.stdout:
        for line in proc.stdout.splitlines():
            if line.startswith("WORKER "):
                print("[%s] %s" % (label, line), flush=True)
    if proc.returncode != 0:
        return None, ["%s: worker exited %d\n%s" % (
            label, proc.returncode, proc.stderr[-2000:])]
    try:
        body = proc.stdout.split("EQUIV-RESULTS-BEGIN")[1]
        body = body.split("EQUIV-RESULTS-END")[0]
        return json.loads(body), []
    except (IndexError, ValueError) as exc:
        return None, ["%s: unparseable worker output: %s\n%s" % (
            label, exc, proc.stdout[-2000:])]


def volumes_equal(va, vb):
    if va == vb == 0.0:
        return True
    denom = max(abs(va), abs(vb), VOL_ABS_FLOOR)
    return abs(va - vb) <= VOL_REL_TOL * denom


def compare(before_rows, after_rows, allow_new_accepts):
    diffs = []
    blocking = []
    bmap = {r["name"]: r for r in before_rows}
    amap = {r["name"]: r for r in after_rows}
    for name in bmap:
        b, a = bmap[name], amap.get(name)
        if a is None:
            blocking.append("%s: missing from after-run" % name)
            continue
        bv, av = b["verdict"], a["verdict"]
        if bv == av == "accept":
            if not volumes_equal(b["volume"], a["volume"]):
                blocking.append(
                    "%s: accepted volume changed %r -> %r" %
                    (name, b["volume"], a["volume"]))
            if b.get("topo") and a.get("topo") and b["topo"] != a["topo"]:
                diffs.append("%s: topology summary changed %r -> %r (info)" %
                             (name, b["topo"], a["topo"]))
        elif bv == av == "refuse":
            if (b.get("kind"), b.get("stage")) != (a.get("kind"), a.get("stage")):
                blocking.append(
                    "%s: refusal changed %s/%s -> %s/%s" %
                    (name, b.get("stage"), b.get("kind"),
                     a.get("stage"), a.get("kind")))
        elif bv == "refuse" and av == "accept":
            msg = ("%s: refuse(%s/%s) -> accept vol=%r" %
                   (name, b.get("stage"), b.get("kind"), a.get("volume")))
            if allow_new_accepts:
                diffs.append(msg + " (allowed)")
            else:
                blocking.append(msg)
        elif bv == "accept" and av == "refuse":
            blocking.append(
                "%s: ACCEPT -> REFUSE (%s/%s); accepted vol was %r" %
                (name, a.get("stage"), a.get("kind"), b.get("volume")))
        else:
            blocking.append("%s: verdict %s -> %s" % (name, bv, av))
    return diffs, blocking


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--before", required=True, help="reference src tree")
    ap.add_argument("--after", required=True, help="changed src tree")
    ap.add_argument("--manifest", default=DEFAULT_MANIFEST)
    ap.add_argument("--fuzz-trials", type=int, default=0,
                    help="extra fuzz trials appended (seeded)")
    ap.add_argument("--fuzz-seed", type=int, default=5)
    ap.add_argument("--fuzz-snap", type=float, default=0.5)
    ap.add_argument("--fuzz-kinds", default="box,cyl")
    ap.add_argument("--skip-probes", action="store_true")
    ap.add_argument("--allow-new-accepts", action="store_true")
    ap.add_argument("--topology", action="store_true",
                    help="also compare cheap topology summaries")
    ap.add_argument("--timeout", type=int, default=1800,
                    help="per-side worker timeout seconds")
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--out", default=None, help="JSON summary path")
    args = ap.parse_args()

    with open(args.manifest) as fh:
        manifest = json.load(fh)
    cases = []
    if not args.skip_probes:
        cases.extend(manifest["probes"])
    for cfg in manifest.get("fuzz", []):
        cases.extend(_expand_fuzz(cfg))
    if args.fuzz_trials:
        cases.extend(_expand_fuzz({
            "trials": args.fuzz_trials, "seed": args.fuzz_seed,
            "snap": args.fuzz_snap, "kinds": args.fuzz_kinds.split(",")}))

    with tempfile.TemporaryDirectory(prefix="equiv-") as tmp:
        cases_path = os.path.join(tmp, "cases.json")
        with open(cases_path, "w") as fh:
            json.dump(cases, fh)
        before_rows, errs = run_side(args.before, cases_path, args.verbose,
                                     args.timeout, "before", args.topology)
        problems = list(errs)
        after_rows, errs = run_side(args.after, cases_path, args.verbose,
                                    args.timeout, "after", args.topology)
        problems.extend(errs)

    summary = {"before": args.before, "after": args.after,
               "n_cases": len(cases), "problems": problems,
               "differences": [], "blocking": []}
    rc = 0
    if before_rows is None or after_rows is None:
        summary["blocking"] = problems
        rc = 1
    else:
        diffs, blocking = compare(before_rows, after_rows,
                                  args.allow_new_accepts)
        summary["differences"] = diffs
        summary["blocking"] = blocking
        tb = sum(1 for r in before_rows if r["verdict"] == "accept")
        ta = sum(1 for r in after_rows if r["verdict"] == "accept")
        summary["accepts"] = {"before": tb, "after": ta}
        secs = sum(r.get("seconds", 0) for r in after_rows)
        summary["after_seconds_total"] = round(secs, 1)
        if blocking:
            rc = 1

    print("cases=%d accepts before=%d after=%d" % (
        summary["n_cases"], summary["accepts"]["before"] if before_rows else -1,
        summary["accepts"]["after"] if after_rows else -1))
    for d in summary["differences"]:
        print("  DIFF", d)
    for b in summary["blocking"]:
        print("  BLOCKING", b)
    if not summary["differences"] and not summary["blocking"]:
        print("  no differences")
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(summary, fh, indent=1)
        print("wrote", args.out)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
