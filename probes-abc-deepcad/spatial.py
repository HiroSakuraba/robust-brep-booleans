#!/usr/bin/env python3
"""Spatial story of the pilot wrong-accept: AABBs, containment, pairwise overlaps."""
import importlib.util
import sys

import numpy as np

sys.path.insert(0, "/home/hatch/workspace/brep-booleans/src")
from brepkernel import verify as kv
from brepkernel.arrange import arrange

_spec = importlib.util.spec_from_file_location(
    "probe_abc", "/home/hatch/workspace/brep-booleans/probes-abc-deepcad/probe_abc.py")
probe_abc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(probe_abc)

from manifold3d import Mesh64, Manifold


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


def M(V, F):
    return Manifold(Mesh64(vert_properties=np.ascontiguousarray(V),
                           tri_verts=np.ascontiguousarray(F.astype(np.uint32))))


def mesh_of(man):
    o = man.to_mesh64()
    return (np.asarray(o.vert_properties, dtype=np.float64),
            np.asarray(o.tri_verts, dtype=np.int64).reshape(-1, 3))


def main():
    t = probe_abc.step_to_mesh(sys.argv[1])
    V, F, diag = t["V"], t["F"], t["diag"]
    dx = 0.35 * diag
    A = shell_meshes(V, F)
    B = [(Va + np.array([dx, 0, 0]), Fa) for Va, Fa in A]
    print(f"dx={dx:.4g}")
    names = []
    for tag, S in (("A", A), ("B", B)):
        for i, (Vs, Fs) in enumerate(S):
            lo, hi = Vs.min(axis=0), Vs.max(axis=0)
            names.append(f"{tag}{i}")
            print(f"{tag}{i}: vol={kv.signed_volume(Vs, Fs):12.2f} "
                  f"bbox=({lo[0]:7.1f},{lo[1]:7.1f},{lo[2]:7.1f})-"
                  f"({hi[0]:7.1f},{hi[1]:7.1f},{hi[2]:7.1f})")
    # containment: is shell i strictly inside shell j? (sample vertices)
    print("\ncontainment (fraction of shell-i verts strictly inside shell-j):")
    Ms = [M(Vs, Fs) for Vs, Fs in A]
    for i in range(len(A)):
        row = []
        for j in range(len(A)):
            if i == j:
                row.append("  --  ")
                continue
            Vi = A[i][0]
            # point-in-mesh via manifold: use to_mesh + ray? use OCP-free approx:
            # winding via ray parity from verify is mesh-based; simpler: use
            # manifold3d's boolean: (Mi - Mj) volume == vol(Mi) -> inside
            diff = Ms[i] - Ms[j]
            Vd, Fd = mesh_of(diff)
            vd = kv.signed_volume(Vd, Fd) if len(Fd) else 0.0
            vi = kv.signed_volume(*A[i])
            row.append("IN " if abs(vd - vi) < 1e-6 * abs(vi) else "out")
        print(f"  A{i}: " + " ".join(row))
    # pairwise intersections within A and across A/B
    print("\npairwise intersection volumes:")
    Mall = [M(Vs, Fs) for Vs, Fs in A + B]
    for i in range(len(Mall)):
        for j in range(i + 1, len(Mall)):
            inter = Mall[i] ^ Mall[j]
            Vi, Fi = mesh_of(inter)
            v = kv.signed_volume(Vi, Fi) if len(Fi) else 0.0
            if abs(v) > 1e-9:
                print(f"  {names[i]} ^ {names[j]} = {v:.4g}")


if __name__ == "__main__":
    main()
