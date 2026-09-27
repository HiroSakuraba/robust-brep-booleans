#!/usr/bin/env python3
"""S3 benchmark: two separated 1000+ face operands, best-of-3.

Each operand is a compound of 170 boxes (1020 faces), the compounds
well separated.  Times the old pipeline (fast_paths=False) against the
S3 fast path (fast_paths=True, fast_path_shadow=False); the plan's S3
pass criterion is >= 5x on large completely separated models.
(The old pipeline is ~O(n^2) here, so best_of=2 keeps the probe
practical: ~10 min total.)
"""
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                           "..", "..", "src"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                           "..", "..", "tests"))


def build_compound(n_boxes, x0):
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Compound
    from OCP.gp import gp_Pnt

    builder = BRep_Builder()
    comp = TopoDS_Compound()
    builder.MakeCompound(comp)
    side = int(round(n_boxes ** (1.0 / 3.0))) + 1
    made = 0
    for i in range(side):
        for j in range(side):
            for k in range(side):
                if made >= n_boxes:
                    break
                s = BRepPrimAPI_MakeBox(
                    gp_Pnt(x0 + i * 2.0, j * 2.0, k * 2.0),
                    1.0, 1.0, 1.0).Shape()
                builder.Add(comp, s)
                made += 1
    return comp


def count_faces(shape):
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopAbs import TopAbs_FACE
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    n = 0
    while ex.More():
        n += 1
        ex.Next()
    return n


def time_call(fn, best_of=2):
    ts = []
    for _ in range(best_of):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    return min(ts)


def main():
    from brepkernel.pipeline import boolean_brep

    n_boxes = 170
    a = build_compound(n_boxes, 0.0)
    b = build_compound(n_boxes, 1000.0)
    fa, fb = count_faces(a), count_faces(b)
    print(f"faces: A={fa} B={fb}")
    assert fa >= 1000 and fb >= 1000

    def old():
        return boolean_brep(a, b, "union", fast_paths=False,
                            crosscheck_ops=False)

    def new():
        return boolean_brep(a, b, "union", fast_paths=True,
                            fast_path_shadow=False,
                            crosscheck_ops=False)

    # one warmup each (OCCT caches, first-call effects)
    old()
    new()
    t_old = time_call(old)
    t_new = time_call(new)
    speedup = t_old / t_new if t_new > 0 else float("inf")
    print(f"old pipeline : {t_old:.3f} s")
    print(f"S3 fast path : {t_new:.3f} s")
    print(f"speedup      : {speedup:.2f}x")
    ok = speedup >= 5.0
    print("PASS" if ok else "FAIL (target >= 5x)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
