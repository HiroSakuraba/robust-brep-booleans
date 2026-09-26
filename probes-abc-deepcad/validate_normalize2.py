#!/usr/bin/env python3
"""Validate the normalize-then-op fix on the pilot model (no numpy round-trips).

normalize(M) = fold the input's own shells with union, all as Manifold objects.
Then the requested op runs on the two normalized Manifolds directly.
Compares union volume to OCCT and replays the arbiter.
"""
import importlib.util
import sys
import time

import numpy as np

sys.path.insert(0, "/home/hatch/workspace/brep-booleans/src")
from brepkernel import verify as kv
from brepkernel.arrange import _to_manifold, ArrangementError

_spec = importlib.util.spec_from_file_location(
    "probe_abc", "/home/hatch/workspace/brep-booleans/probes-abc-deepcad/probe_abc.py")
probe_abc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(probe_abc)


def split_shells(V, F, tag0):
    """Per-shell (Manifold, n_faces) with global face ids starting at tag0."""
    from manifold3d import Manifold
    labels = kv.shell_face_labels(F)
    out = []
    k = tag0
    for lab in sorted(np.unique(labels).tolist()):
        m = labels == lab
        Fm = F[m]
        used = np.unique(Fm)
        remap = np.full(len(V), -1, dtype=np.int64)
        remap[used] = np.arange(len(used))
        man, _ = _to_manifold(
            {"V": np.ascontiguousarray(V[used]),
             "F": np.ascontiguousarray(remap[Fm])}, k)
        out.append(man)
        k += int(m.sum())
    return out


def normalize(manifolds):
    acc = manifolds[0]
    for m in manifolds[1:]:
        acc = acc + m
        if acc.status().name != "NoError":
            raise ArrangementError("normalize fold failed")
    return acc


def mesh_of(man):
    o = man.to_mesh64()
    V = np.array(np.asarray(o.vert_properties, dtype=np.float64),
                 dtype=np.float64)
    F = np.array(np.asarray(o.tri_verts, dtype=np.int64).reshape(-1, 3),
                 dtype=np.int64)
    return V, F


def main():
    t = probe_abc.step_to_mesh(sys.argv[1])
    V, F, diag = t["V"], t["F"], t["diag"]
    dx = 0.35 * diag
    Vb = V + np.array([dx, 0, 0])
    nA = len(F)
    t0 = time.time()
    mA = normalize(split_shells(V, F, 0))
    mB = normalize(split_shells(Vb, F, nA))
    Ua, Ub = mesh_of(mA), mesh_of(mB)
    print(f"normalize: A {kv.connected_shells(F)} shells -> "
          f"{kv.connected_shells(Ua[1])} shells vol={kv.signed_volume(*Ua):.6g}; "
          f"took {time.time()-t0:.2f}s")
    for op in ("union", "intersection", "difference"):
        r0 = time.time()
        mR = {"union": mA + mB, "intersection": mA ^ mB,
              "difference": mA - mB}[op]
        assert mR.status().name == "NoError", mR.status()
        Vr, Fr = mesh_of(mR)
        vol = kv.signed_volume(Vr, Fr)
        print(f"{op}: vol={vol:.6g} shells={kv.connected_shells(Fr)} "
              f"({time.time()-r0:.2f}s)")
        if op == "union":
            oc = probe_abc.occt_fuse_oracle(t["shape"], dx)
            arb = probe_abc.arbiter(oc["fused"], oc["bbox"], Vr, Fr,
                                    solids6=oc["solids6"])
            print(f"  OCCT fuse vol={oc['volume']:.6g} "
                  f"err={vol/oc['volume']-1:+.4%}")
            tr = arb["truth_or"]
            print(f"  arbiter: {arb['n_agree']}/{arb['n_decided']} "
                  f"(kin_tout={tr['kernel_in_truth_out']} "
                  f"kout_tin={tr['kernel_out_truth_in']})")


if __name__ == "__main__":
    main()
