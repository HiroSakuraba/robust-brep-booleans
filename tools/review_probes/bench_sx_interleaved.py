#!/usr/bin/env python3
"""Interleaved SX (point-verdict bbox short-circuit) vs S5-base plate benchmark.

Compares wall time plus verdict preservation (volume, region stats) across
alternating runs to average out machine noise.
"""
import subprocess

VENV = "/home/hatch/workspace/brep-booleans/.venv/bin/python"
SX_SRC = "/home/hatch/workspace/brep-gates"
BASE_SRC = "/tmp/wt-s5base"

BENCH = """
import sys, time, json
sys.path.insert(0, "@@SRC@@/src")
sys.path.insert(0, "@@SRC@@/tests")
sys.path.insert(0, "@@SRC@@/tools/review_probes")
import bench_plate256 as bp
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.gp import gp_Pnt
from brepkernel.pipeline import boolean_brep
drilled = bp.build_drilled()
slot = BRepPrimAPI_MakeBox(gp_Pnt(3.0, 3.5, -0.25), 2.0, 1.0, 1.0).Shape()
boolean_brep(drilled, slot, "difference")  # warmup
t0 = time.perf_counter()
_out, _rep = boolean_brep(drilled, slot, "difference")
dt = time.perf_counter() - t0
_asm = _rep["stages"]["assembly"]
print(json.dumps({"t": round(dt, 3), "volume": _asm["volume"],
                  "region_stats": _asm["region_stats"]}))
"""


def one(src):
    code = BENCH.replace("@@SRC@@", src)
    p = subprocess.run([VENV, "-c", code], capture_output=True, text=True,
                       env={"PATH": "/usr/bin:/bin"}, cwd=src)
    if p.returncode != 0:
        raise RuntimeError(f"bench failed for {src}:\n{p.stderr[-2000:]}")
    import json as _json
    return _json.loads(p.stdout.strip().splitlines()[-1])


def median(xs):
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2


sx_ts, base_ts, vols = [], [], []
for i in range(6):
    r = one(SX_SRC)
    assert "point_verdict_bbox_shortcircuit" in open(
        f"{SX_SRC}/src/brepkernel/assembly.py").read(), "SX code not loaded"
    sx_ts.append(r["t"])
    vols.append(("sx", r["volume"], r["region_stats"]))
    print(f"sx[{i}]={r['t']:.3f} vol={r['volume']:.4f}", flush=True)
    r = one(BASE_SRC)
    assert "point_verdict_bbox_shortcircuit" not in open(
        f"{BASE_SRC}/src/brepkernel/assembly.py").read(), "base polluted"
    base_ts.append(r["t"])
    vols.append(("base", r["volume"], r["region_stats"]))
    print(f"base[{i}]={r['t']:.3f} vol={r['volume']:.4f}", flush=True)

print(f"SX   median={median(sx_ts):.3f} samples={[f'{x:.3f}' for x in sx_ts]}")
print(f"base median={median(base_ts):.3f} samples={[f'{x:.3f}' for x in base_ts]}")
v0 = vols[0][1]
assert all(abs(v[1] - v0) < 1e-9 for v in vols), f"VOLUME CHANGED: {vols}"
rs0 = vols[0][2]
assert all(v[2] == rs0 for v in vols), "REGION STATS CHANGED"
print(f"volume identical ({v0:.4f}); region stats identical")
