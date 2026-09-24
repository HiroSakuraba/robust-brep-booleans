#!/usr/bin/env python3
"""Regression replay: run stage-2 booleans with the NEW arrange() and diff
against the recorded baselines (results_abc.json, results_deepcad.json).

Compares per-op verdicts (arranged/refused) and union volumes.
Writes replay_regression.json.
"""
import importlib.util
import json
import os
import sys
import traceback

import numpy as np

sys.path.insert(0, "/home/hatch/workspace/brep-booleans/src")
from brepkernel.arrange import arrange, ArrangementError
from brepkernel import verify as kv

HERE = os.path.dirname(os.path.abspath(__file__))


def load_probe(name):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


probe_abc = load_probe("probe_abc")
probe_deepcad = load_probe("probe_deepcad")


def run3(V, F, shift):
    Vb = (V + np.asarray(shift, dtype=np.float64)).astype(np.float64)
    out = {}
    for op in ("union", "intersection", "difference"):
        try:
            r = arrange({"V": V, "F": F}, {"V": Vb, "F": F}, op)
            Vr = np.asarray(r["V"], dtype=np.float64)
            Fr = np.asarray(r["F"], dtype=np.int64).reshape(-1, 3)
            out[op] = {"verdict": "arranged",
                       "volume": float(kv.signed_volume(Vr, Fr)),
                       "shells": int(kv.connected_shells(Fr))}
        except ArrangementError as e:
            out[op] = {"verdict": "refused", "detail": str(e)[:120]}
        except Exception as e:
            out[op] = {"verdict": "EXCEPTION",
                       "detail": f"{type(e).__name__}: {e}"[:120]}
    return out


def abc_part():
    recs = json.load(open(os.path.join(HERE, "results_abc.json")))["models"]
    rows = []
    for m in recs:
        p = m.get("step", "")
        if not p or not os.path.exists(p):
            rows.append({"id": m["id"], "skipped": "step missing"})
            continue
        try:
            t = probe_abc.step_to_mesh(p)
        except Exception as e:
            rows.append({"id": m["id"], "skipped": f"tess fail {e}"[:80]})
            continue
        V, F, diag = t["V"], t["F"], t["diag"]
        new = run3(V, F, [0.35 * diag, 0, 0])
        s1 = m["stage1"]
        old = {op: m["stage2"][op].get("verdict") for op in
               ("union", "intersection", "difference")}
        newv = {op: new[op]["verdict"] for op in old}
        row = {"id": m["id"], "shells_in": s1.get("shells"),
               "defective": bool(s1.get("boundary_edges") or
                                 s1.get("nonmanifold_edges") or
                                 not s1.get("directed_closed")),
               "verdict_changed": {op: (old[op], newv[op]) for op in old
                                   if old[op] != newv[op]},
               "new": {op: {"verdict": new[op]["verdict"]} for op in old}}
        if (old["union"] == "arranged" and newv["union"] == "arranged"
                and m["stage2"]["union"].get("volume") is not None):
            ov = m["stage2"]["union"]["volume"]
            nv = new["union"]["volume"]
            row["union_vol_old"] = ov
            row["union_vol_new"] = nv
            row["union_vol_reldiff"] = abs(nv - ov) / max(1e-300, abs(ov))
            row["shells_old"] = m["oracle"].get("kernel_union_shells")
            row["shells_new"] = new["union"]["shells"]
        rows.append(row)
    return rows


def deepcad_part():
    import glob
    from OCP.STEPControl import STEPControl_Reader
    d = json.load(open(os.path.join(HERE, "results_deepcad.json")))
    # match extracts to results entries by (input_tris, diag)
    entries = {}
    for m in d["models"]:
        entries[(m["input_tris"], round(m["diag"], 6))] = m
    rows = []
    for fp in sorted(glob.glob(os.path.join(HERE, "scratch/extract/m*.step"))):
        rdr = STEPControl_Reader()
        if rdr.ReadFile(fp) != 1:
            continue
        rdr.TransferRoots()
        shape = rdr.OneShape()
        try:
            V, F, tinfo = probe_deepcad.tessellate_shape(shape)
        except Exception:
            continue
        key = (len(F), round(tinfo["diag"], 6))
        m = entries.get(key)
        if m is None:
            rows.append({"file": os.path.basename(fp),
                         "skipped": f"no results match for {key}"})
            continue
        ax = m.get("shift_axis", 0)
        shift = [0.0, 0.0, 0.0]
        shift[ax] = 0.35 * m["diag"]
        new = run3(V, F, shift)
        old = {op: m["ops"][op].get("status") for op in
               ("union", "intersection", "difference")}
        newv = {op: new[op]["verdict"] for op in old}
        row = {"name": m["name"], "defective":
               bool(m.get("input_boundary") or m.get("input_nonmanifold")
                    or not m.get("input_closed")),
               "verdict_changed": {op: (old[op], newv[op]) for op in old
                                   if old[op] != newv[op]},
               "new": {op: {"verdict": new[op]["verdict"]} for op in old}}
        if old["union"] == "arranged" and newv["union"] == "arranged":
            ov = m["ops"]["union"]["volume"]
            nv = new["union"]["volume"]
            row["union_vol_reldiff"] = abs(nv - ov) / max(1e-300, abs(ov))
        rows.append(row)
    return rows


def main():
    print("ABC replay...", flush=True)
    abc = abc_part()
    print("DeepCAD replay...", flush=True)
    dc = deepcad_part()
    json.dump({"abc": abc, "deepcad": dc},
              open(os.path.join(HERE, "replay_regression.json"), "w"),
              indent=1)
    # summary
    for tag, rows in (("ABC", abc), ("DeepCAD", dc)):
        n = len(rows)
        changed = [r for r in rows if r.get("verdict_changed")]
        exc = [r for r in rows
               if any(v.get("verdict") == "EXCEPTION" for v in
                      r.get("new", {}).values())]
        voldiff = [(r.get("id") or r.get("name"),
                    r.get("union_vol_reldiff"))
                   for r in rows if r.get("union_vol_reldiff", 0) > 1e-9]
        print(f"{tag}: {n} models replayed, {len(changed)} verdict changes, "
              f"{len(exc)} exceptions, {len(voldiff)} union-vol diffs")
        for r in changed:
            print("   CHANGED:", r.get("id") or r.get("name"),
                  r["verdict_changed"],
                  "defective_input=" + str(r.get("defective")))
        for rid, vd in voldiff[:12]:
            print(f"   VOLDIFF {rid}: rel={vd:.3g}")


if __name__ == "__main__":
    main()
