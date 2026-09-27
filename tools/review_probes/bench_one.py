#!/usr/bin/env python3
"""Single plate boolean timing for the interleaved S4 benchmark.

Usage: PYTHONPATH=<branch>/src:tests python bench_one.py
Prints one wall-clock number: the boolean_brep(drilled, slot) time.
"""
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "src"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))

import bench_plate256 as bp
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.gp import gp_Pnt
from brepkernel.pipeline import boolean_brep

drilled = bp.build_drilled()
slot = BRepPrimAPI_MakeBox(gp_Pnt(3.0, 3.5, -0.25), 2.0, 1.0, 1.0).Shape()

# Warmup (caches, OCCT init) not timed
boolean_brep(drilled, slot, "difference")

t0 = time.perf_counter()
_out, _rep = boolean_brep(drilled, slot, "difference")
dt = time.perf_counter() - t0
print(f"{dt:.3f}")
