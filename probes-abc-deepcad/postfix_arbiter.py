#!/usr/bin/env python3
"""Post-fix arbiter check on the multi-shell ABC cases whose volumes changed."""
import importlib.util
import os
import sys

import numpy as np

sys.path.insert(0, "/home/hatch/workspace/brep-booleans/src")
from brepkernel.arrange import arrange
from brepkernel import verify as kv

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "probe_abc", os.path.join(HERE, "probe_abc.py"))
probe_abc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(probe_abc)

CASES = [
    "scratch/pilot/pilot.step",
    "/home/hatch/workspace/abc_chunk_27/closed_multishell/00273581_b5e213618b746c2c6eff9376_step_000.step",
    "/home/hatch/workspace/abc_chunk_27/closed_multishell/00276289_bac679d3b1093419f4d82f65_step_010.step",
    "/home/hatch/workspace/abc_chunk_27/closed_multishell/00277962_57f92b31b54dff10aa5fb4ad_step_000.step",
]

for path in CASES:
    t = probe_abc.step_to_mesh(path)
    V, F, diag = t["V"], t["F"], t["diag"]
    dx = 0.35 * diag
    Vb = V + np.array([dx, 0, 0])
    r = arrange({"V": V, "F": F}, {"V": Vb, "F": F}, "union")
    Vr = np.asarray(r["V"], dtype=np.float64)
    Fr = np.asarray(r["F"], dtype=np.int64).reshape(-1, 3)
    vol = kv.signed_volume(Vr, Fr)
    oc = probe_abc.occt_fuse_oracle(t["shape"], dx)
    arb = probe_abc.arbiter(oc["fused"], oc["bbox"], Vr, Fr,
                            solids6=oc["solids6"])
    tr = arb["truth_or"]
    print(f"{os.path.basename(path)[:30]}: vol err vs OCCT {vol/oc['volume']-1:+.3%}, "
          f"arbiter {arb['n_agree']}/{arb['n_decided']} "
          f"(kin_tout={tr['kernel_in_truth_out']}, kout_tin={tr['kernel_out_truth_in']})")
