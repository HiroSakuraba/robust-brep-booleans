"""Attribution checks: (1) is 00651489 non-manifold at the B-rep level?
(2) are the chi mismatches explained by open fused tessellations?
(3) arbiter inside-fractions (non-triviality)."""
import sys, json
sys.path.insert(0, '/home/hatch/workspace/brep-booleans/probes-abc-deepcad')
import numpy as np
from probe_deepcad import (extract_entry, tessellate_shape, longest_axis,
                           occt_fuse_volume_shells_chi, arbiter, CD_LISTING,
                           edge_stats)
from brepkernel.arrange import arrange
from OCP.STEPControl import STEPControl_Reader
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_ShapeEnum


def load_model(zip_name):
    names = [tuple(n) for n in json.load(open(CD_LISTING))]
    entry = [n for n in names if n[0] == zip_name][0]
    data = extract_entry(entry)
    fp = '/tmp/deepcad_stress/inv.step'
    open(fp, 'wb').write(data)
    rdr = STEPControl_Reader()
    assert rdr.ReadFile(fp) == 1
    rdr.TransferRoots()
    return rdr.OneShape()


def main():
    # (1) B-rep level validity of 00651489
    shape = load_model('cad_step/00651489.step')
    ana = BRepCheck_Analyzer(shape)
    print("00651489 BRepCheck_Analyzer.IsValid:", ana.IsValid(), flush=True)

    def count(kind):
        e = TopExp_Explorer(shape, kind)
        n = 0
        while e.More():
            n += 1
            e.Next()
        return n
    print("solids:", count(TopAbs_ShapeEnum.TopAbs_SOLID),
          "shells:", count(TopAbs_ShapeEnum.TopAbs_SHELL),
          "faces:", count(TopAbs_ShapeEnum.TopAbs_FACE), flush=True)

    # (2) chi-mismatch models: is the OCCT fused tessellation closed?
    for zn in ['cad_step/00645276.step', 'cad_step/00022637.step',
               'cad_step/00892289.step', 'cad_step/00460815.step',
               'cad_step/00669242.step', 'cad_step/00022250.step']:
        shape = load_model(zn)
        V, F, info = tessellate_shape(shape)
        ax, _ = longest_axis(V)
        shift = [0.0, 0.0, 0.0]
        shift[ax] = 0.35 * info['diag']
        oc = occt_fuse_volume_shells_chi(shape, shift)
        # re-tessellate fused shape and check boundary edges
        Vf, Ff, _ = tessellate_shape(oc["fused_shape"])
        nb, nn = edge_stats(Ff)
        print(f"{zn}: fused-tess boundary={nb} nonmanifold={nn} "
              f"fused_tris={len(Ff)}", flush=True)

    # (3) arbiter inside fractions on two models (one clean, one bad)
    rng = np.random.default_rng(7)
    for zn in ['cad_step/00645276.step', 'cad_step/00651489.step']:
        shape = load_model(zn)
        V, F, info = tessellate_shape(shape)
        ax, _ = longest_axis(V)
        shift = [0.0, 0.0, 0.0]
        shift[ax] = 0.35 * info['diag']
        Vb = V + np.array(shift)
        r = arrange({"V": V, "F": F.astype(np.uint32)},
                    {"V": Vb, "F": F.astype(np.uint32)}, "union")
        Vu = np.asarray(r["V"], dtype=np.float64)
        Fu = np.asarray(r["F"], dtype=np.int64).reshape(-1, 3)
        oc = occt_fuse_volume_shells_chi(shape, shift)
        # replicate arbiter internals to get inside fractions
        from OCP.BRepClass3d import BRepClass3d_SolidClassifier
        from OCP.gp import gp_Pnt
        from OCP.TopAbs import TopAbs_State
        from OCP.Bnd import Bnd_Box
        from OCP.BRepBndLib import BRepBndLib
        box = Bnd_Box()
        BRepBndLib.Add_s(oc["fused_shape"], box)
        x0, x1 = box.GetXMin(), box.GetXMax()
        y0, y1 = box.GetYMin(), box.GetYMax()
        z0, z1 = box.GetZMin(), box.GetZMax()
        n = 800
        pts = np.column_stack([rng.uniform(x0, x1, n), rng.uniform(y0, y1, n),
                               rng.uniform(z0, z1, n)])
        diag = info['diag']
        nin = 0
        for x, y, z in pts:
            cl = BRepClass3d_SolidClassifier(oc["fused_shape"], gp_Pnt(x, y, z),
                                             diag * 1e-7)
            if cl.State() == TopAbs_State.TopAbs_IN:
                nin += 1
        print(f"{zn}: OCCT-inside fraction {nin}/{n}", flush=True)


if __name__ == '__main__':
    main()
