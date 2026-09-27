#!/usr/bin/env python3
"""Interleaved S5 vs main (S4) plate benchmark, fixed path handling."""
import subprocess
import sys
import os

VENV = "/home/hatch/workspace/brep-booleans/.venv/bin/python"
S5_SRC = "/home/hatch/workspace/brep-gates"
MAIN_SRC = "/tmp/wt-main"

BENCH = """
import sys, time
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
print("{:.3f}".format(time.perf_counter()-t0))
"""


def one(src):
    code = BENCH.replace("@@SRC@@", src)
    # verify the right same_domain is loaded
    check = "import brepkernel.same_domain as m; print('S5' if 'S5' in open(m.__file__).read() else 'main')"
    p = subprocess.run(
        [VENV, "-c", code],
        capture_output=True, text=True,
        env={"PATH": "/usr/bin:/bin"},
        cwd=src,
    )
    out = p.stdout.strip().split()[-1]
    # sanity: confirm which code ran
    q = subprocess.run(
        [VENV, "-c", f"import sys; sys.path.insert(0, '{src}/src'); {check}"],
        capture_output=True, text=True, env={"PATH": "/usr/bin:/bin"},
    )
    tag = q.stdout.strip()
    return float(out), tag


def median(xs):
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2


s5, main = [], []
for i in range(6):
    t, tag = one(S5_SRC)
    assert tag == "S5", f"S5 run loaded {tag}"
    s5.append(t)
    print(f"s5[{i}]={t:.3f} ({tag})", flush=True)
    t, tag = one(MAIN_SRC)
    assert tag == "main", f"main run loaded {tag}"
    main.append(t)
    print(f"main[{i}]={t:.3f} ({tag})", flush=True)

print(f"S5   median={median(s5):.3f} samples={[f'{x:.3f}' for x in s5]}")
print(f"main median={median(main):.3f} samples={[f'{x:.3f}' for x in main]}")
