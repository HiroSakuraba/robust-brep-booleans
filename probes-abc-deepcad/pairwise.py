#!/usr/bin/env python3
"""Pairwise isolation: boolean each shell of A against each shell of B.

Determines whether manifold3d fails on specific shell pairs (geometry issue)
or only on the multi-shell assembly (assembly issue).
"""
import importlib.util
import sys

import numpy as np

sys.path.insert(0, "/home/hatch/workspace/brep-booleans/src")
from brepkernel import verify as kv
from brepkernel.arrange import arrange, ArrangementError

_spec = importlib.util.spec_from_file_location(
    "probe_abc", "/home/hatch/workspace/brep-booleans/probes-abc-deepcad/probe_abc.py")
probe_abc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(probe_abc)


def shell_meshes(V, F):
    labels = kv.shell_face_labels(F)
    out = []
    for k in sorted(np.unique(labels).tolist()):
        m = labels == k
        Fm = F[m]
        used = np.unique(Fm)
        remap = np.full(len(V), -1, dtype=np.int64)
        remap[used] = np.arange(len(used))
        out.append((V[used], remap[Fm]))
    return out


def aabb(V):
    return V.min(axis=0), V.max(axis=0)


def main():
    path = sys.argv[1]
    t = probe_abc.step_to_mesh(path)
    V, F, diag = t["V"], t["F"], t["diag"]
    dx = 0.35 * diag
    A = shell_meshes(V, F)
    B = [(Va + np.array([dx, 0, 0]), Fa) for Va, Fa in A]
    print(f"input shells: {len(A)}, dx={dx:.4g}")
    for i, (Va, Fa) in enumerate(A):
        lo, hi = aabb(Va)
        print(f"A[{i}]: tris={len(Fa)} vol={kv.signed_volume(Va, Fa):.6g} "
              f"bbox=({lo[0]:.1f},{lo[1]:.1f},{lo[2]:.1f})-({hi[0]:.1f},{hi[1]:.1f},{hi[2]:.1f})")
    for i in range(len(A)):
        for j in range(len(B)):
            Va, Fa = A[i]
            Vb, Fb = B[j]
            lo = np.maximum(Va.min(axis=0), Vb.min(axis=0))
            hi = np.minimum(Va.max(axis=0), Vb.max(axis=0))
            overlap = np.prod(np.maximum(hi - lo, 0.0))
            try:
                res = arrange({"V": Va, "F": Fa}, {"V": Vb, "F": Fb}, "union")
            except ArrangementError as e:
                print(f"pair A{i}xB{j}: REFUSED {e}; aabb_overlap_vol={overlap:.6g}")
                continue
            Vr = np.asarray(res["V"], dtype=np.float64)
            Fr = np.asarray(res["F"], dtype=np.int64).reshape(-1, 3)
            nsh = kv.connected_shells(Fr)
            vol = kv.signed_volume(Vr, Fr)
            va = kv.signed_volume(Va, Fa)
            vb = kv.signed_volume(Vb, Fb)
            print(f"pair A{i}xB{j}: shells={nsh} union_vol={vol:.6g} "
                  f"(va+vb={va+vb:.6g}) merged={'YES' if vol < va+vb-1e-6*(va+vb) else 'no'} "
                  f"aabb_overlap_vol={overlap:.6g}")


if __name__ == "__main__":
    main()
