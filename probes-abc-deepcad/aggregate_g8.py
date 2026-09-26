#!/usr/bin/env python3
"""G8 aggregation: rebuild results_g8.json from per-case result files.

Recovers the driver-side timeout records that were lost when the sweep
driver crashed in its own aggregation step (the worker-written files all
survived; the 9 timeouts are reconstructed from the sweep log and the
case files, and persisted next to the others so the results directory
is complete and re-aggregation is idempotent).
"""

import datetime
import json
import os
import re
import subprocess
import sys

REPO = "/home/hatch/workspace/brep-gates-g8"
sys.path.insert(0, os.path.join(REPO, "probes-abc-deepcad"))
import probe_brep_g8 as driver  # noqa: E402

SCRATCH = os.path.join(REPO, "probes-abc-deepcad", "scratch_g8")
CASE_DIR = os.path.join(SCRATCH, "cases")
RES_DIR = os.path.join(SCRATCH, "results")
OUT_JSON = os.path.join(REPO, "probes-abc-deepcad", "results_g8.json")
LOG = "/tmp/g8sweep.log"


def main():
    # Reconstruct driver-side timeout records from the log.
    timeout_re = re.compile(
        r"\[(\d+)/45\] (\S+): timeout / ([\d.]+)s")
    timeouts = {}
    with open(LOG) as fh:
        for line in fh:
            m = timeout_re.search(line)
            if m:
                timeouts[m.group(2)] = float(m.group(3))

    results = []
    n_reconstructed = 0
    for fname in sorted(os.listdir(CASE_DIR)):
        if not fname.endswith(".json"):
            continue
        case = json.load(open(os.path.join(CASE_DIR, fname)))
        cid = case["case_id"]
        rpath = os.path.join(RES_DIR, cid + ".json")
        if os.path.exists(rpath):
            rec = json.load(open(rpath))
            if "pair" not in rec:
                rec["pair"] = case["pair"]
            results.append(rec)
            continue
        assert cid in timeouts, f"no result and no timeout log for {cid}"
        rec = {"case_id": cid, "corpus": case["corpus"],
               "pair": case["pair"], "model_a": case["model_a"],
               "model_b": case["model_b"], "op": case["op"],
               "outcome": "timeout", "worker_rc": -9,
               "wall_s": timeouts[cid],
               "note": "no result file after 120 s; worker killed "
                       "(record reconstructed from sweep log after "
                       "driver crash; harness bug since fixed)"}
        with open(rpath, "w") as fh:
            json.dump(rec, fh, indent=1)
        results.append(rec)
        n_reconstructed += 1

    results.sort(key=lambda r: r["case_id"])
    summary = driver.summarize(results)
    checks = driver.inclusion_exclusion(results)
    commit = subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"], cwd=REPO,
        capture_output=True, text=True).stdout.strip()
    report = {
        "manifest": {
            "date": datetime.date.today().isoformat(),
            "repo": REPO, "commit": commit,
            "pipeline": "boolean_brep defaults (base_tol=1e-7, "
                        "contact_tol=4e-7, no tuning)",
            "venv": driver.VENV_PY, "src_path": driver.SRC,
            "timeout_s": driver.TIMEOUT,
            "max_workers": driver.MAX_WORKERS,
            "pairing": "distinct models; B rigid-motion translated so "
                       "its bbox center lands near A's center "
                       "(0.25*diagA x offset); guaranteed overlap",
            "recovery_note": (
                f"{n_reconstructed} timeout records reconstructed from "
                "the sweep log after the driver crashed in its own "
                "aggregation step (KeyError on worker records missing "
                "'pair'); worker files were intact. Harness fixed: "
                "worker now records 'pair'; driver persists "
                "driver-side records to disk."),
        },
        "summary": summary,
        "inclusion_exclusion": checks,
        "results": results,
    }
    with open(OUT_JSON, "w") as fh:
        json.dump(report, fh, indent=1)
    print(f"cases: {len(results)}, reconstructed: {n_reconstructed}")
    print(f"wrote {OUT_JSON}")
    print(json.dumps(summary["by_outcome"], indent=1))
    print("refusal histogram:")
    for k, v in sorted(summary["refusal_histogram"].items(),
                       key=lambda kv: -kv[1]):
        print(f"  {v:3d}  {k}")
    print("inclusion-exclusion checks:", len(checks),
          "flags:", sum(1 for c in checks if c["flag"]))


if __name__ == "__main__":
    main()
