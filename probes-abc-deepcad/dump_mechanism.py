#!/usr/bin/env python3
"""Mechanism dump for the ABC/DeepCAD silent wrong-accept class.

For each target STEP: tessellate exactly like probe_abc.py, shift a copy by
0.35*diag along +x, run brepkernel.arrange.arrange('union'), and dump
per-shell signed volumes + per-shell chi for input A and the union result.
Also runs the OCCT BRepAlgoAPI_Fuse oracle for comparison.

Usage: dump_mechanism.py <step1> [<step2> ...] [--op union]
"""
import importlib.util
import json
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


def shell_dump(V, F, tag):
    labels = kv.shell_face_labels(F)
    out = []
    for k in sorted(np.unique(labels).tolist()):
        m = labels == k
        Fm = F[m]
        # shell-local vertices for a meaningful chi
        used = np.unique(Fm)
        remap = np.full(len(V), -1, dtype=np.int64)
        remap[used] = np.arange(len(used))
        Fl = remap[Fm]
        Vl = V[used]
        vol = kv.signed_volume(Vl, Fl)
        chi = kv.euler_chi(Vl, Fl)
        out.append({"shell": int(k), "n_tris": int(m.sum()),
                    "signed_vol": vol, "chi": chi})
    return out


def occt_fuse_solid_by_solid(shapeA, dx):
    """Mirror of probe_abc.occt_fuse_oracle (fuse solid-by-solid: the
    whole-compound Fuse silently drops solids in this OCCT build)."""
    from OCP.gp import gp_Trsf, gp_Vec
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopAbs import TopAbs_SOLID
    from OCP.TopoDS import TopoDS
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    solids = []
    e = TopExp_Explorer(shapeA, TopAbs_SOLID)
    while e.More():
        solids.append(TopoDS.Solid(e.Current()))
        e.Next()
    trsf = gp_Trsf()
    trsf.SetTranslation(gp_Vec(dx, 0.0, 0.0))
    moved = [BRepBuilderAPI_Transform(s, trsf, True).Shape() for s in solids]
    acc = solids[0]
    for s in solids[1:] + moved:
        f = BRepAlgoAPI_Fuse(acc, s)
        if not f.IsDone():
            raise RuntimeError("BRepAlgoAPI_Fuse did not finish")
        acc = f.Shape()
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(acc, props)
    return float(props.Mass())


def main():
    op = "union"
    paths = [p for p in sys.argv[1:] if not p.startswith("--")]
    for path in paths:
        print("=" * 70)
        print("MODEL:", path)
        t = probe_abc.step_to_mesh(path)
        V, F, diag = t["V"], t["F"], t["diag"]
        print(f"input: n_tris={len(F)} diag={diag:.6g}")
        for s in shell_dump(V, F, "A"):
            print(f"  A shell {s['shell']}: tris={s['n_tris']} "
                  f"signed_vol={s['signed_vol']:.6g} chi={s['chi']}")
        print(f"  A total |vol| = {sum(abs(s['signed_vol']) for s in shell_dump(V, F, 'A')):.6g}")
        Vb = V + np.array([0.35 * diag, 0.0, 0.0])
        t0 = time.time()
        try:
            res = arrange({"V": V, "F": F}, {"V": Vb, "F": F}, op)
        except ArrangementError as e:
            print("  arrange REFUSED:", e)
            continue
        Vr = np.asarray(res["V"], dtype=np.float64)
        Fr = np.asarray(res["F"], dtype=np.int64).reshape(-1, 3)
        print(f"  arrange ok in {time.time()-t0:.2f}s: n_tris={len(Fr)}")
        tot = 0.0
        for s in shell_dump(Vr, Fr, "R"):
            print(f"  R shell {s['shell']}: tris={s['n_tris']} "
                  f"signed_vol={s['signed_vol']:.6g} chi={s['chi']}")
            tot += s["signed_vol"]
        print(f"  R total signed vol = {tot:.6g}")
        # OCCT fuse oracle (solid-by-solid, mirroring probe_abc)
        ovol = occt_fuse_solid_by_solid(t["shape"], 0.35 * diag)
        print(f"  OCCT fuse vol = {ovol:.6g}  "
              f"kernel/occt-1 = {tot/ovol-1:+.4%}")


if __name__ == "__main__":
    main()
