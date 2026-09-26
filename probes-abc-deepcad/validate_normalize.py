#!/usr/bin/env python3
"""Validate the normalize-then-op fix idea on the pilot model.

normalize(M) = fold the mesh's own shells with union (set-union of the input).
Then run the requested op on the two normalized meshes, all-at-once.
Compare union volume to OCCT and run the arbiter (replay of probe_abc.arbiter).
"""
import importlib.util
import sys
import time

import numpy as np

sys.path.insert(0, "/home/hatch/workspace/brep-booleans/src")
from brepkernel import verify as kv
from brepkernel.arrange import arrange, ArrangementError

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


def normalize(V, F):
    """Set-union of the mesh's own shells (fold)."""
    shells = shell_meshes(V, F)
    acc = M(*shells[0])
    for Vs, Fs in shells[1:]:
        acc = acc + M(Vs, Fs)
        if acc.status().name != "NoError":
            raise ArrangementError("normalize fold failed")
    return mesh_of(acc)


def main():
    t = probe_abc.step_to_mesh(sys.argv[1])
    V, F, diag = t["V"], t["F"], t["diag"]
    dx = 0.35 * diag
    Vb = V + np.array([dx, 0, 0])
    t0 = time.time()
    Ua = normalize(V, F)
    Ub = normalize(Vb, F)
    print(f"normalize: A {kv.connected_shells(F)} shells -> "
          f"{kv.connected_shells(Ua[1])} shells, "
          f"vol {kv.signed_volume(*Ua):.6g}; took {time.time()-t0:.2f}s")
    for op in ("union", "intersection", "difference"):
        r0 = time.time()
        res = arrange({"V": Ua[0], "F": Ua[1]}, {"V": Ub[0], "F": Ub[1]}, op)
        Vr = np.asarray(res["V"], dtype=np.float64)
        Fr = np.asarray(res["F"], dtype=np.int64).reshape(-1, 3)
        vol = kv.signed_volume(Vr, Fr)
        print(f"{op}: vol={vol:.6g} shells={kv.connected_shells(Fr)} "
              f"({time.time()-r0:.2f}s)")
        if op == "union":
            # arbiter replay vs OCCT solid-by-solid fused shape
            oc = probe_abc.occt_fuse_oracle(t["shape"], dx)
            arb = probe_abc.arbiter(oc["fused"], oc["bbox"], Vr, Fr,
                                    solids6=oc["solids6"])
            print(f"  OCCT fuse vol={oc['volume']:.6g} "
                  f"err={vol/oc['volume']-1:+.4%}")
            print(f"  arbiter: {arb['n_agree']}/{arb['n_decided']} "
                  f"truth_or={arb['truth_or']}")


if __name__ == "__main__":
    main()
