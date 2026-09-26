#!/usr/bin/env python3
"""G8 real-data refusal-rate sweep driver.

For each model pair and each of union/difference/intersection, runs one
isolated worker process (probes-abc-deepcad/g8_worker.py) with a 120 s
join so one hang or crash cannot wedge the sweep. At most 4 workers run
at once.

Produces probes-abc-deepcad/results_g8.json:
  manifest, per-case records (model ids, op, outcome, refusal
  stage/kind, wall time, volumes, tolerances), inclusion-exclusion
  sanity on accepted triples, and a summary table with the refusal
  histogram.
"""

import datetime
import json
import os
import subprocess
import sys
import time

REPO = "/home/hatch/workspace/brep-gates-g8"
SRC = os.path.join(REPO, "src")
VENV_PY = "/home/hatch/workspace/brep-booleans/.venv/bin/python"
WORKER = os.path.join(REPO, "probes-abc-deepcad", "g8_worker.py")
SCRATCH = os.path.join(REPO, "probes-abc-deepcad", "scratch_g8")
CASE_DIR = os.path.join(SCRATCH, "cases")
RES_DIR = os.path.join(SCRATCH, "results")
OUT_JSON = os.path.join(REPO, "probes-abc-deepcad", "results_g8.json")

ABC_DIR = "/home/hatch/workspace/abc_chunk_27/step"
DC_DIR = "/home/hatch/workspace/brep-booleans/probes-abc-deepcad/scratch/extract"
NIST_DIR = "/tmp/nist_g8/NIST-PMI-STEP-Files"

TIMEOUT = 120.0
MAX_WORKERS = 4
OPS = ("union", "difference", "intersection")


def model_id(path):
    base = os.path.basename(path)
    for ext in (".stp", ".step"):
        if base.endswith(ext):
            base = base[: -len(ext)]
    return base


def build_cases():
    cases = []

    def add(corpus, pair_id, pa, pb):
        for op in OPS:
            cid = f"{corpus}-{pair_id}-{op}"
            cases.append({
                "case_id": cid, "corpus": corpus, "pair": pair_id,
                "model_a": model_id(pa), "model_b": model_id(pb),
                "path_a": pa, "path_b": pb, "op": op,
                "out": os.path.join(RES_DIR, cid + ".json"),
            })

    nist_pairs = [
        ("nist_ftc_06_asme1_ap242-e2.stp", "nist_ftc_07_asme1_ap242-e2.stp"),
        ("nist_ftc_08_asme1_ap242-e2.stp", "nist_ftc_09_asme1_ap242-e1.stp"),
        ("nist_ctc_01_asme1_ap242-e1.stp", "nist_ctc_02_asme1_ap242-e2.stp"),
    ]
    for i, (fa, fb) in enumerate(nist_pairs, 1):
        pa, pb = os.path.join(NIST_DIR, fa), os.path.join(NIST_DIR, fb)
        assert os.path.exists(pa), pa
        assert os.path.exists(pb), pb
        add("nist", f"p{i:02d}", pa, pb)

    dc_files = [os.path.join(DC_DIR, f"m{i:03d}.step") for i in range(10)]
    for f in dc_files:
        assert os.path.exists(f), f
    for i in range(0, 10, 2):
        add("deepcad", f"p{i // 2 + 1:02d}", dc_files[i], dc_files[i + 1])

    abc_files = []
    for root, _dirs, files in os.walk(ABC_DIR):
        for f in files:
            if f.endswith("_step_000.step"):
                abc_files.append(os.path.join(root, f))
    abc_files.sort()
    assert len(abc_files) >= 4300, len(abc_files)
    for i, off in enumerate([0, 660, 1320, 1980, 2640, 3300, 3960], 1):
        add("abc", f"p{i:02d}", abc_files[off], abc_files[off + 1])
    return cases


def summarize(results):
    hist = {}
    summary = {"total": len(results), "by_outcome": {},
               "refusal_histogram": hist, "by_corpus": {},
               "by_op": {}, "timeouts": 0, "worker_died": 0,
               "exceptions": 0, "accepted": 0, "refused": 0,
               "timeout_s": 0.0, "accept_s": 0.0}
    for r in results:
        oc = r.get("outcome")
        summary["by_outcome"][oc] = summary["by_outcome"].get(oc, 0) + 1
        corp = summary["by_corpus"].setdefault(
            r["corpus"], {"accepted": 0, "refused": 0,
                          "other": 0, "cases": []})
        op = summary["by_op"].setdefault(
            r["op"], {"accepted": 0, "refused": 0, "other": 0})
        if oc == "accepted":
            summary["accepted"] += 1
            corp["accepted"] += 1
            op["accepted"] += 1
            summary["accept_s"] += r.get("wall_s", 0.0)
        elif oc == "refused":
            summary["refused"] += 1
            corp["refused"] += 1
            op["refused"] += 1
            key = f"{r.get('refusal_stage')}/{r.get('refusal_kind')}"
            hist[key] = hist.get(key, 0) + 1
        elif oc == "timeout":
            summary["timeouts"] += 1
            corp["other"] += 1
            op["other"] += 1
            summary["timeout_s"] += r.get("wall_s", 0.0)
        elif oc == "worker_died":
            summary["worker_died"] += 1
            corp["other"] += 1
            op["other"] += 1
        else:
            summary["exceptions"] += 1
            corp["other"] += 1
            op["other"] += 1
        corp["cases"].append(r["case_id"])
    return summary


def inclusion_exclusion(results):
    """On pairs where all three ops accepted, check the volume
    identities: volU + volI = volA + volB and volD = volA - volI.
    Independent of the kernel; uses OCCT volumes only."""
    by_pair = {}
    for r in results:
        if r.get("outcome") != "accepted":
            continue
        by_pair.setdefault((r["corpus"], r["pair"]), {})[r["op"]] = r
    checks = []
    for (corpus, pair), ops in sorted(by_pair.items()):
        if set(ops) != set(OPS):
            continue
        u, d, i_ = ops["union"], ops["difference"], ops["intersection"]
        va, vb = u["volume_a"], u["volume_b"]
        err1 = abs((u["result_volume"] + i_["result_volume"]) - (va + vb))
        err2 = abs(d["result_volume"] - (va - i_["result_volume"]))
        scale = max(va, vb, 1e-12)
        checks.append({"corpus": corpus, "pair": pair,
                       "identity_u_plus_i": {"abs_err": err1,
                                             "rel_err": err1 / scale},
                       "identity_d": {"abs_err": err2,
                                      "rel_err": err2 / scale},
                       "flag": (err1 / scale > 0.02
                                or err2 / scale > 0.02)})
    return checks


def main():
    os.makedirs(CASE_DIR, exist_ok=True)
    os.makedirs(RES_DIR, exist_ok=True)
    cases = build_cases()
    print(f"{len(cases)} cases; workers cap {MAX_WORKERS}; "
          f"timeout {TIMEOUT}s", flush=True)

    running = []  # (case, Popen, t0, state dict with timed_out flag)
    results = []

    def collect():
        for item in list(running):
            case, proc, t0, state = item
            rc = proc.poll()
            if rc is None and time.perf_counter() - t0 > TIMEOUT:
                proc.kill()
                proc.wait()
                state["timed_out"] = True
            if proc.poll() is None:
                continue
            running.remove(item)
            rec = finalize_case(case, proc.returncode, t0,
                                state["timed_out"])
            results.append(rec)
            print(f"[{len(results)}/{len(cases)}] {rec['case_id']}: "
                  f"{rec['outcome']} "
                  f"{rec.get('refusal_stage', '')}/"
                  f"{rec.get('refusal_kind', '')} "
                  f"{rec.get('wall_s', 0):.1f}s", flush=True)

    pending = list(cases)
    while pending or running:
        while pending and len(running) < MAX_WORKERS:
            case = pending.pop(0)
            case_path = os.path.join(CASE_DIR, case["case_id"] + ".json")
            with open(case_path, "w") as fh:
                json.dump(case, fh)
            env = dict(os.environ)
            env["G8_SRC"] = SRC
            env["PYTHONPATH"] = SRC
            proc = subprocess.Popen(
                [VENV_PY, WORKER, case_path], env=env,
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            running.append((case, proc, time.perf_counter(),
                            {"timed_out": False}))
        time.sleep(2.0)
        collect()

    summary = summarize(results)
    checks = inclusion_exclusion(results)
    commit = subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"], cwd=REPO,
        capture_output=True, text=True).stdout.strip()
    report = {
        "manifest": {
            "date": datetime.date.today().isoformat(),
            "repo": REPO, "commit": commit,
            "pipeline": "boolean_brep defaults (base_tol=1e-7, "
                        "contact_tol=4e-7, no tuning)",
            "venv": VENV_PY, "src_path": SRC,
            "timeout_s": TIMEOUT, "max_workers": MAX_WORKERS,
            "pairing": "distinct models; B rigid-motion translated so "
                       "its bbox center lands near A's center "
                       "(0.25*diagA x offset); guaranteed overlap",
        },
        "summary": summary,
        "inclusion_exclusion": checks,
        "results": results,
    }
    with open(OUT_JSON, "w") as fh:
        json.dump(report, fh, indent=1)
    print(f"wrote {OUT_JSON}", flush=True)
    print(json.dumps(summary["by_outcome"], indent=1))
    print("refusal histogram:")
    for k, v in sorted(summary["refusal_histogram"].items(),
                       key=lambda kv: -kv[1]):
        print(f"  {v:3d}  {k}")


def finalize_case(case, rc, t0, timed_out):
    wall = time.perf_counter() - t0
    if os.path.exists(case["out"]):
        rec = json.load(open(case["out"]))
        rec["worker_rc"] = rc
        if rec.get("outcome") is None:
            rec["outcome"] = "exception"
            rec["note"] = "worker wrote a file but no outcome; rc=" + str(rc)
        return rec
    if timed_out:
        rec = {"case_id": case["case_id"], "corpus": case["corpus"],
               "pair": case["pair"], "model_a": case["model_a"],
               "model_b": case["model_b"], "op": case["op"],
               "outcome": "timeout", "worker_rc": rc, "wall_s": wall,
               "note": "no result file after 120 s; worker killed"}
    else:
        rec = {"case_id": case["case_id"], "corpus": case["corpus"],
               "pair": case["pair"], "model_a": case["model_a"],
               "model_b": case["model_b"], "op": case["op"],
               "outcome": "worker_died",
               "worker_rc": rc, "wall_s": wall,
               "note": "worker exited without writing a result file"}
    # Persist driver-side records too: the aggregation step must not
    # depend on driver memory surviving.
    with open(case["out"], "w") as fh:
        json.dump(rec, fh, indent=1)
    return rec


if __name__ == "__main__":
    main()
