#!/usr/bin/env python3
"""S7 gate benchmark: 64-hole pattern, sequential reference vs optimized.

Plan criterion: optimized at least 2x faster than the sequential
reference before the many-tool path may be labeled optimized.
"""
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
sys.path.insert(0, os.path.dirname(__file__))

from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.gp import gp_Pnt
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_FACE

from brepkernel.many import boolean_brep_many
from brepkernel.assembly import _shape_volume


def count_faces(shape):
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    n = 0
    while ex.More():
        n += 1
        ex.Next()
    return n


def plain_plate():
    return BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 8.0, 8.0, 0.5).Shape()


def hole_tool(i, j, n):
    step = 8.0 / n
    x = step / 2 + i * step
    y = step / 2 + j * step
    return BRepPrimAPI_MakeBox(gp_Pnt(x - 0.15, y - 0.15, -0.25),
                              0.3, 0.3, 1.0).Shape()


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    base = plain_plate()
    tools = [hole_tool(i, j, n) for i in range(n) for j in range(n)]

    t0 = time.perf_counter()
    ref_out, ref_sum = boolean_brep_many(base, tools, "difference")
    t_ref = time.perf_counter() - t0
    assert all(r["accepted"] for r in ref_sum["steps"])
    v_ref = abs(float(_shape_volume(ref_out)))
    f_ref = count_faces(ref_out)

    t0 = time.perf_counter()
    opt_out, opt_sum = boolean_brep_many(base, tools, "difference",
                                         optimized=True)
    t_opt = time.perf_counter() - t0
    assert opt_sum["optimized"] is True
    assert all(r["accepted"] for r in opt_sum["steps"])
    v_opt = abs(float(_shape_volume(opt_out)))
    f_opt = count_faces(opt_out)

    print(f"holes: {n*n}", flush=True)
    print(f"sequential reference: {t_ref:.1f}s vol={v_ref:.6f} faces={f_ref}",
          flush=True)
    print(f"optimized           : {t_opt:.1f}s vol={v_opt:.6f} faces={f_opt}",
          flush=True)
    print(f"speedup: {t_ref/t_opt:.2f}x "
          f"(criterion >= 2.0x: {'HIT' if t_ref/t_opt >= 2.0 else 'MISSED'})",
          flush=True)
    print(f"volume match: {abs(v_ref-v_opt) <= 1e-9*max(1,v_ref)} "
          f"faces match: {f_ref == f_opt}", flush=True)
    print(f"opt detail: skipped={opt_sum['n_skipped']} "
          f"groups={opt_sum['n_groups']} "
          f"combine_ops={opt_sum['n_combine_ops']}", flush=True)


if __name__ == "__main__":
    main()
