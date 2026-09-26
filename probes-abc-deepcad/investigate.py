"""Follow-up investigation: (a) 00651489 input non-manifold: real or weld artifact?
(b) arbiter on the 6 models with kernel-vs-OCCT chi mismatch."""
import sys, json
sys.path.insert(0, '/home/hatch/workspace/brep-booleans/probes-abc-deepcad')
import numpy as np
from probe_deepcad import (extract_entry, tessellate_shape, longest_axis,
                           occt_fuse_volume_shells_chi, arbiter, CD_LISTING,
                           edge_stats)
from brepkernel.arrange import arrange
from brepkernel import verify
from OCP.STEPControl import STEPControl_Reader


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
    # (a) 00651489 input inspection
    shape = load_model('cad_step/00651489.step')
    V, F, info = tessellate_shape(shape)
    print("00651489 tessellated:", V.shape, F.shape, info, flush=True)
    print("edge stats post-weld:", edge_stats(F), flush=True)
    ok, cinfo = verify.check_closure_directed(V, F)
    print("directed closure:", ok, cinfo, flush=True)
    print("chi:", verify.euler_chi(V, F), "shells:", verify.connected_shells(F),
          flush=True)
    sv = verify.shell_signed_volumes(V, F)
    print("signed vols:", sv, flush=True)
    # where are the non-manifold edges?
    from collections import Counter
    c = Counter()
    for a, b, d_ in F:
        for u, v in ((a, b), (b, d_), (d_, a)):
            c[(u, v) if u < v else (v, u)] += 1
    nm = [e for e, k in c.items() if k > 2]
    print("n non-manifold edges:", len(nm), flush=True)
    for e in nm[:4]:
        print("  edge verts:", V[e[0]], V[e[1]], "uses:", c[e], flush=True)

    # (b) arbiter on chi-mismatch models
    targets = ['cad_step/00022250.step', 'cad_step/00645276.step',
               'cad_step/00022637.step', 'cad_step/00892289.step',
               'cad_step/00460815.step', 'cad_step/00669242.step']
    rng = np.random.default_rng(20260924)
    for zn in targets:
        try:
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
            assert oc["ok"], oc.get("error")
            kvol = abs(float(np.sum(np.einsum(
                "ij,ij->i", Vu[Fu[:, 0]],
                np.cross(Vu[Fu[:, 1]], Vu[Fu[:, 2]]))) / 6.0))
            verr = abs(kvol - oc["volume"]) / oc["volume"]
            ab = arbiter(oc["fused_shape"], Vu, Fu, rng, n=800)
            print(f"{zn}: verr={verr*100:.3f}% kchi={verify.euler_chi(Vu,Fu)} "
                  f"ochi={oc['chi']} arbiter_agreement={ab['agreement']} "
                  f"on_surface={ab['n_on_surface']}", flush=True)
        except Exception as e:
            print(f"{zn}: INVESTIGATION FAILED {type(e).__name__}: {e}",
                  flush=True)


if __name__ == '__main__':
    main()
