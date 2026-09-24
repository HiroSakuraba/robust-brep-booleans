#!/usr/bin/env python3
"""Pre-pass sampler for the ABC chunk stress round (24 Sept 2026).

Enumerates STEP files in an extracted chunk dir, draws a uniform random pool,
tessellates each via probe_abc.step_to_mesh (same tessellation pipeline as the
pilot), records defect stats, and selects ~40 worst-case + ~10 clean controls,
copying the selected STEP files into a destination dir for probe_abc.py.

Worst-case ordering: open (non-closed) inputs first, then multi-shell, then
boundary/non-manifold edge counts. Clean controls: closed, 1 shell, 0 defects.
Models above --max-tris are skipped (keeps Stage 2 worker times bounded) and
counted honestly.
"""

import argparse
import json
import os
import random
import shutil
import sys
import time
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import probe_abc  # noqa: E402


def defect_score(r):
    s = r["stats"]
    return (s["boundary_edges"] + s["nonmanifold_edges"]
            + (0 if s["directed_closed"] else 10**9)
            + abs(s["shells"] - 1) * 1000)


def is_clean(r):
    s = r["stats"]
    return (s["directed_closed"] and s["boundary_edges"] == 0
            and s["nonmanifold_edges"] == 0 and s["shells"] == 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chunk-dir", required=True,
                    help="dir of extracted *.step files")
    ap.add_argument("--pool", type=int, default=300,
                    help="uniform random pool size to tessellate")
    ap.add_argument("--seed", type=int, default=20260924)
    ap.add_argument("--n-worst", type=int, default=40)
    ap.add_argument("--n-clean", type=int, default=10)
    ap.add_argument("--max-tris", type=int, default=200000,
                    help="skip models with more tessellated tris (honest count)")
    ap.add_argument("--dest", required=True,
                    help="dir to copy selected STEP files into")
    ap.add_argument("--stats-out", required=True)
    a = ap.parse_args()

    step_files = []
    for root, _dirs, fnames in os.walk(a.chunk_dir):
        for fn in fnames:
            if fn.endswith(".step"):
                step_files.append(os.path.join(root, fn))
    files = sorted(step_files)
    print(f"chunk dir has {len(files)} STEP files", flush=True)
    rng = random.Random(a.seed)
    pool = rng.sample(files, min(a.pool, len(files)))

    recs = []
    t0 = time.time()
    for i, path in enumerate(pool):
        f = os.path.basename(path)
        r = {"id": f[:-5], "file": f, "path": path}
        try:
            m = probe_abc.step_to_mesh(path)
            r["ok"] = True
            r["stats"] = probe_abc.defect_stats(m["V"], m["F"])
        except Exception:
            r["ok"] = False
            r["error"] = traceback.format_exc()[-300:]
        recs.append(r)
        if (i + 1) % 50 == 0:
            print(f"  tessellated {i + 1}/{len(pool)}", flush=True)
    print(f"tessellation pass: {time.time() - t0:.1f}s", flush=True)

    ok = [r for r in recs if r.get("ok") and r["stats"]["n_tris"] <= a.max_tris]
    skipped_big = sum(1 for r in recs
                      if r.get("ok") and r["stats"]["n_tris"] > a.max_tris)
    failed = sum(1 for r in recs if not r.get("ok"))
    print(f"ok={len(ok)} failed_import={failed} skipped_huge={skipped_big}",
          flush=True)

    worst = sorted(ok, key=defect_score, reverse=True)
    cleans = [r for r in ok if is_clean(r)]
    picked, seen = [], set()
    for r in worst[:a.n_worst]:
        picked.append(r)
        seen.add(r["id"])
    for r in cleans[:a.n_clean]:
        if r["id"] not in seen:
            picked.append(r)
            seen.add(r["id"])

    os.makedirs(a.dest, exist_ok=True)
    for r in picked:
        shutil.copy(r["path"], os.path.join(a.dest, r["file"]))
    print(f"selected {len(picked)} "
          f"({sum(1 for r in picked if is_clean(r))} clean) -> {a.dest}",
          flush=True)

    with open(a.stats_out, "w") as fh:
        json.dump({
            "chunk_dir": a.chunk_dir,
            "pool_size": len(pool), "pool_seed": a.seed,
            "ok": len(ok), "failed_import": failed,
            "skipped_huge": skipped_big, "max_tris": a.max_tris,
            "selected_ids": [r["id"] for r in picked],
            "selected": [{"id": r["id"], "stats": r["stats"]} for r in picked],
        }, fh, indent=1, default=str)
    print("wrote", a.stats_out)


if __name__ == "__main__":
    main()
