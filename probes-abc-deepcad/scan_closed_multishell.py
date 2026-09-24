"""Scan an extracted ABC chunk dir for CLEAN CLOSED MULTI-SHELL models.

The chunk-27 scale-up round (sample_chunk.py) selected 40 worst-case + 10
clean controls, but is_clean() requires shells==1, so the pilot's failure
class (clean, closed, multi-shell -- the pilot model was a closed 3-shell
solid whose union came out +18.4% vs the OCCT oracle) was never exercised.
This script fills that gap: tessellate files in the chunk and keep the
first --keep models with directed_closed, 0 boundary edges,
0 non-manifold edges, shells > 1, tris <= --max-tris.

Usage: scan_closed_multishell.py --chunk-dir DIR --dest DIR --stats-out F
                                   [--keep 8] [--max-files 2000] [--max-tris 200000]
"""
import argparse, json, os, random, shutil, sys, time, traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import probe_abc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chunk-dir", required=True)
    ap.add_argument("--dest", required=True)
    ap.add_argument("--stats-out", required=True)
    ap.add_argument("--keep", type=int, default=8)
    ap.add_argument("--max-files", type=int, default=2000,
                    help="cap on files tessellated (honest count)")
    ap.add_argument("--max-tris", type=int, default=200000)
    ap.add_argument("--seed", type=int, default=20260925)
    a = ap.parse_args()

    step_files = []
    for root, _dirs, fnames in os.walk(a.chunk_dir):
        for fn in fnames:
            if fn.endswith(".step"):
                step_files.append(os.path.join(root, fn))
    files = sorted(step_files)
    rng = random.Random(a.seed)
    rng.shuffle(files)  # different order than the scale-up pool (seed 20260924)
    print(f"chunk dir has {len(files)} STEP files; scanning up to {a.max_files}",
          flush=True)

    kept, scanned, failed = [], 0, 0
    t0 = time.time()
    for path in files[:a.max_files]:
        if len(kept) >= a.keep:
            break
        f = os.path.basename(path)
        scanned += 1
        try:
            m = probe_abc.step_to_mesh(path)
            s = probe_abc.defect_stats(m["V"], m["F"])
        except Exception:
            failed += 1
            continue
        if (s["directed_closed"] and s["boundary_edges"] == 0
                and s["nonmanifold_edges"] == 0 and s["shells"] > 1
                and s["n_tris"] <= a.max_tris):
            kept.append({"id": f[:-5], "file": f, "path": path, "stats": s})
            print(f"  KEEP {f[:-5]} shells={s['shells']} tris={s['n_tris']} "
                  f"chi={s['chi']}", flush=True)
        if scanned % 100 == 0:
            print(f"  scanned {scanned}, kept {len(kept)}, "
                  f"{time.time()-t0:.0f}s", flush=True)

    os.makedirs(a.dest, exist_ok=True)
    for r in kept:
        shutil.copy(r["path"], os.path.join(a.dest, r["file"]))
    with open(a.stats_out, "w") as fh:
        json.dump({"chunk_dir": a.chunk_dir, "seed": a.seed,
                   "scanned": scanned, "failed_import": failed,
                   "kept": len(kept), "keep_target": a.keep,
                   "max_tris": a.max_tris,
                   "kept_ids": [r["id"] for r in kept],
                   "kept_stats": [{"id": r["id"], "stats": r["stats"]}
                                  for r in kept]}, fh, indent=1, default=str)
    print(f"scanned {scanned} files, kept {len(kept)} closed multi-shell "
          f"models in {time.time()-t0:.0f}s -> {a.dest}")


if __name__ == "__main__":
    main()
