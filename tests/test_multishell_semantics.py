"""Multi-shell solid-semantics regression tests (24 Sept 2026).

The ABC/DeepCAD independent-oracle rounds found silent wrong-accepts in the
mesh layer: a mesh can be closed, consistently oriented, and have plausible
Euler/volume while denoting the WRONG SET (all four mesh verify functions
passed on the wrong mesh). These tests pin the solid-semantics contract for
anonymous multi-shell meshes:

  1. overlapping positive shells in one operand -> exact union volume
     (the old code double-counted the overlap)
  2. three overlapping shells (ABC pilot structure)
  3. duplicate coincident components -> single volume
  4. disjoint multi-solid assembly -> normalization is a no-op (bit-identical)
  5. outer shell + reversed-orientation inner cavity shell (hollow cube):
     the hole MUST survive union and difference. Cavities are never unioned
     as material.
  6. two tetrahedra welded at exactly one vertex: every edge looks manifold
     (directed closure passes) but the shared vertex link is two disjoint
     loops -> typed refusal.
  7. self-intersecting single shell (closed, edge-manifold) -> typed refusal.
  8. shell-order permutation: permuted input shells must not materially
     change the result (volume / topology / membership).
  9. three-level nesting (material island in void in material, vol 57).
  10. 00277962 failure mode: two shells interpenetrating in a thin slab
      must be FOLDED (union volume, not the sum) -- pins the correct
      expectation an earlier MECHANISM.md entry misread as a regression.
  11. genuinely disjoint shells take the "disjoint passthrough" (status
      pinned, bit-identical).
  12. negative signed result volume (inside-out/corrupted result) is
      refused by the semantic leg for all three ops.

Deterministic, tiny, no corpus needed. Uses arrange() directly on raw
meshes (the layer where the bug lived).
"""

import itertools
import sys

import numpy as np

sys.path.insert(0, "src")
from brepkernel.arrange import arrange, ArrangementError
try:
    from brepkernel.normalize import prepare_operand
except ImportError:  # pre-normalize.py code: fall back to the old _prepare
    from brepkernel.arrange import _prepare as _old_prepare

    def prepare_operand(proxy, tag0, name):
        from manifold3d import Manifold
        man = _old_prepare(proxy, tag0, name)
        out = man.to_mesh64()
        import numpy as np
        return {"manifold": man,
                "V": np.asarray(out.vert_properties, dtype=np.float64),
                "F": np.asarray(out.tri_verts,
                                dtype=np.int64).reshape(-1, 3)}
from brepkernel.verify import signed_volume, connected_shells, _winding_number


# --------------------------------------------------------------------------
# mesh builders
# --------------------------------------------------------------------------

def box_mesh(lo, hi):
    """Outward-oriented unit of a box; positive signed volume."""
    x0, y0, z0 = lo
    x1, y1, z1 = hi
    V = np.array([[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
                  [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]],
                 dtype=float)
    F = np.array([[0, 2, 1], [0, 3, 2],
                  [4, 5, 6], [4, 6, 7],
                  [0, 1, 5], [0, 5, 4],
                  [3, 6, 2], [3, 7, 6],
                  [0, 4, 7], [0, 7, 3],
                  [1, 2, 6], [1, 6, 5]], dtype=int)
    assert signed_volume(V, F) > 0
    return V, F


def tetra_mesh(verts):
    V = np.array(verts, dtype=float)
    F = np.array([[1, 2, 3], [0, 2, 1], [0, 3, 2], [0, 1, 3]], dtype=int)
    if signed_volume(V, F) < 0:
        F = F[:, ::-1].copy()
    assert signed_volume(V, F) > 0
    return V, F


def concat(meshes):
    """Concatenate meshes into one multi-shell mesh."""
    Vs, Fs, off = [], [], 0
    for V, F in meshes:
        Vs.append(V)
        Fs.append(F + off)
        off += len(V)
    return np.vstack(Vs), np.vstack(Fs)


def vol(mesh):
    return abs(signed_volume(mesh["V"], mesh["F"]))


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def expect_refusal(name, kind_fragment, fn):
    try:
        fn()
    except ArrangementError as e:
        ok = kind_fragment.lower() in str(e).lower()
        return check(name, ok, f"refused kind={getattr(e, 'kind', '?')}")
    except Exception as e:  # noqa: BLE001 - wrong exception type is a failure
        return check(name, False, f"wrong exception: {type(e).__name__}: {e}")
    return check(name, False, "no refusal raised")


# --------------------------------------------------------------------------
# 1-3: overlapping / duplicate shells in one operand
# --------------------------------------------------------------------------

def t1_two_overlapping_shells():
    # A1=[0,2]x[0,1]x[0,1] vol 2, A2=[1,3]x[0,1]x[0,1] vol 2, overlap vol 1
    # => A denotes 3. B=[2.5,4]x[0,1]x[0,1] vol 1.5, B cap A = 0.5
    A = concat([box_mesh([0, 0, 0], [2, 1, 1]),
                box_mesh([1, 0, 0], [3, 1, 1])])
    B = box_mesh([2.5, 0, 0], [4, 1, 1])
    dA = {"V": A[0], "F": A[1]}
    dB = {"V": B[0], "F": B[1]}
    ok = True
    m = arrange(dA, dB, "union")
    ok &= check("t1 union vol", abs(vol(m) - 4.0) < 1e-9, f"vol={vol(m)}")
    m = arrange(dA, dB, "intersection")
    ok &= check("t1 intersection vol", abs(vol(m) - 0.5) < 1e-9,
                f"vol={vol(m)}")
    m = arrange(dA, dB, "difference")
    ok &= check("t1 difference vol", abs(vol(m) - 2.5) < 1e-9,
                f"vol={vol(m)}")
    return ok


def t2_three_overlapping_shells():
    # ABC pilot structure: C1,C2 overlap by 1; C2,C3 overlap by 0.5;
    # C1,C3 disjoint. Union of operand = 2+2+1.5-1-0.5 = 4.
    # B=[3.5,5] vol 1.5 meets C3 in 0.5 => total union 5.0.
    A = concat([box_mesh([0, 0, 0], [2, 1, 1]),
                box_mesh([1, 0, 0], [3, 1, 1]),
                box_mesh([2.5, 0, 0], [4, 1, 1])])
    B = box_mesh([3.5, 0, 0], [5, 1, 1])
    m = arrange({"V": A[0], "F": A[1]}, {"V": B[0], "F": B[1]}, "union")
    return check("t2 three-shell union vol", abs(vol(m) - 5.0) < 1e-9,
                 f"vol={vol(m)}")


def t3_duplicate_coincident():
    b = box_mesh([0, 0, 0], [1, 1, 1])
    A = concat([b, b])  # exact duplicate component
    B = box_mesh([2, 0, 0], [3, 1, 1])
    m = arrange({"V": A[0], "F": A[1]}, {"V": B[0], "F": B[1]}, "union")
    ok = check("t3 duplicate union vol", abs(vol(m) - 2.0) < 1e-9,
               f"vol={vol(m)}")
    ok &= check("t3 duplicate shells", connected_shells(m["F"]) == 2,
                f"shells={connected_shells(m['F'])}")
    return ok


# --------------------------------------------------------------------------
# 4: disjoint assembly is untouched
# --------------------------------------------------------------------------

def t4_disjoint_assembly_bit_identical():
    A = concat([box_mesh([0, 0, 0], [1, 1, 1]),
                box_mesh([2, 0, 0], [3, 1, 1])])
    prep = prepare_operand({"V": A[0], "F": A[1]}, 0, "input A")
    same_v = np.array_equal(np.asarray(prep["V"]), A[0])
    same_f = np.array_equal(np.asarray(prep["F"], dtype=np.int64), A[1])
    ok = check("t4 normalization bit-identical", same_v and same_f)
    B = box_mesh([4, 0, 0], [5, 1, 1])
    m = arrange({"V": A[0], "F": A[1]}, {"V": B[0], "F": B[1]}, "union")
    ok &= check("t4 union vol", abs(vol(m) - 3.0) < 1e-9, f"vol={vol(m)}")
    ok &= check("t4 union shells", connected_shells(m["F"]) == 3,
                f"shells={connected_shells(m['F'])}")
    return ok


# --------------------------------------------------------------------------
# 5: hollow cube -- cavities are never unioned as material
# --------------------------------------------------------------------------

def hollow_cube():
    outer = box_mesh([0, 0, 0], [2, 2, 2])          # vol 8, outward
    inner = box_mesh([0.5, 0.5, 0.5], [1.5, 1.5, 1.5])
    inner = (inner[0], inner[1][:, ::-1].copy())      # reversed: cavity
    assert signed_volume(*inner) < 0
    return concat([outer, inner])


def t5_hollow_cube():
    H = hollow_cube()
    B = box_mesh([3, 0, 0], [4, 1, 1])               # disjoint, vol 1
    dH = {"V": H[0], "F": H[1]}
    dB = {"V": B[0], "F": B[1]}
    ok = True
    m = arrange(dH, dB, "union")
    ok &= check("t5 union vol", abs(vol(m) - 8.0) < 1e-9, f"vol={vol(m)}")
    # the hole must survive: cavity point outside, material point inside
    w = _winding_number(np.array([[1.0, 1.0, 1.0], [0.25, 0.25, 0.25]]),
                        m["V"], m["F"])
    ok &= check("t5 cavity survives union", w[0] == 0 and w[1] == 1,
                f"winding={w.tolist()}")
    m = arrange(dH, dB, "difference")
    ok &= check("t5 difference vol", abs(vol(m) - 7.0) < 1e-9,
                f"vol={vol(m)}")
    w = _winding_number(np.array([[1.0, 1.0, 1.0]]), m["V"], m["F"])
    ok &= check("t5 cavity survives difference", w[0] == 0,
                f"winding={w.tolist()}")
    # filling the void partially: union with a box strictly inside the cavity
    F = box_mesh([0.75, 0.75, 0.75], [1.25, 1.25, 1.25])   # vol 0.125
    m = arrange(dH, {"V": F[0], "F": F[1]}, "union")
    ok &= check("t5 partial fill vol", abs(vol(m) - 7.125) < 1e-9,
                f"vol={vol(m)}")
    return ok


# --------------------------------------------------------------------------
# 6-7: validation refusals
# --------------------------------------------------------------------------

def t6_welded_tetras():
    t1 = tetra_mesh([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    t2 = tetra_mesh([[0, 0, 0], [-1, 0, 0], [0, -1, 0], [0, 0, -1]])
    # weld at vertex 0 (shared index 0)
    V = np.vstack([t1[0][:1], t1[0][1:], t2[0][1:]])
    F = np.vstack([t1[1], np.where(t2[1] == 0, 0, t2[1] + 3)])
    B = box_mesh([5, 5, 5], [6, 6, 6])
    return expect_refusal(
        "t6 welded tetras refused", "link",
        lambda: arrange({"V": V, "F": F}, {"V": B[0], "F": B[1]}, "union"))


def t7_self_intersecting_shell():
    V, F = box_mesh([0, 0, 0], [1, 1, 1])
    V = V.copy()
    V[6] = [0.5, 0.5, -1.0]   # pull corner through the bottom face
    B = box_mesh([5, 5, 5], [6, 6, 6])
    return expect_refusal(
        "t7 self-intersecting shell refused", "self-intersect",
        lambda: arrange({"V": V, "F": F}, {"V": B[0], "F": B[1]}, "union"))


# --------------------------------------------------------------------------
# 8: shell-order permutation metamorphic test
# --------------------------------------------------------------------------

def t8_shell_permutation():
    shells = [box_mesh([0, 0, 0], [2, 1, 1]),
              box_mesh([1, 0, 0], [3, 1, 1]),
              box_mesh([2.5, 0, 0], [4, 1, 1])]
    B = box_mesh([3.5, 0, 0], [5, 1, 1])
    dB = {"V": B[0], "F": B[1]}
    vols, nsh, agree = [], [], True
    grids = []
    for perm in [ (0, 1, 2), (2, 1, 0), (1, 2, 0), (2, 0, 1) ]:
        A = concat([shells[i] for i in perm])
        m = arrange({"V": A[0], "F": A[1]}, dB, "union")
        vols.append(vol(m))
        nsh.append(connected_shells(m["F"]))
        # membership on a fixed witness grid over the bbox
        lo = np.array([0.0, 0.0, 0.0])
        g1d = np.linspace(0.05, 4.95, 9)
        gx, gy, gz = np.meshgrid(g1d, np.linspace(0.05, 0.95, 3),
                                 np.linspace(0.05, 0.95, 3), indexing="ij")
        pts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)
        grids.append(_winding_number(pts, m["V"], m["F"]))
    ok = check("t8 permutation volumes",
               max(vols) - min(vols) < 1e-9 * max(vols),
               f"vols={['%.12g' % v for v in vols]}")
    ok &= check("t8 permutation topology", len(set(nsh)) == 1,
                f"shells={nsh}")
    ok &= check("t8 permutation membership",
                all(np.array_equal(grids[0], g) for g in grids[1:]))
    return ok


# --------------------------------------------------------------------------
# 9: three-level nesting -- material island inside a void inside material
# --------------------------------------------------------------------------

def t9_nested_island():
    outer = box_mesh([0, 0, 0], [4, 4, 4])              # vol 64, material
    void = box_mesh([1, 1, 1], [3, 3, 3])
    void = (void[0], void[1][:, ::-1].copy())           # reversed: void
    island = box_mesh([1.5, 1.5, 1.5], [2.5, 2.5, 2.5]) # vol 1, material
    A = concat([outer, void, island])
    prep = prepare_operand({"V": A[0], "F": A[1]}, 0, "input A")
    ok = check("t9 nested island volume",
               abs(prep["volume"] - 57.0) < 1e-9, f"vol={prep['volume']}")
    # island interior must be IN, void ring must be OUT, outer ring IN
    for p, want, tag in [([2.0, 2.0, 2.0], 1, "island"),
                         ([1.2, 1.2, 1.2], 0, "void ring"),
                         ([0.5, 0.5, 0.5], 1, "outer ring")]:
        got = _winding_number(np.array([p]), prep["V"], prep["F"])[0]
        ok &= check(f"t9 nested island {tag}", got == want,
                    f"winding={got}")
    # through arrange: union with a disjoint box
    B = box_mesh([5, 0, 0], [6, 1, 1])
    m = arrange({"V": A[0], "F": A[1]}, {"V": B[0], "F": B[1]}, "union")
    ok &= check("t9 union vol", abs(vol(m) - 58.0) < 1e-9, f"vol={vol(m)}")
    return ok


# --------------------------------------------------------------------------
# 10: 00277962 failure mode -- interpenetrating slab must be FOLDED
# --------------------------------------------------------------------------

def t10_interpenetrating_slab_folded():
    # 00277962 (24 Sept 2026): two shells with the same x/y footprint whose
    # z-ranges overlap in a thin slab genuinely interpenetrate (OCCT
    # pairwise fuse confirms ~overlap volume). The denoted set is the
    # UNION -- overlap counted once -- not the sum of shell volumes.
    # (An earlier MECHANISM.md entry misread the v0.7 result as a
    # regression: it compared against a v0.6 baseline that had never
    # folded the overlap. The folded result matches the f1d2e331 target
    # and is within tessellation bias of OCCT.)
    A = concat([box_mesh([0, 0, 0], [2, 2, 3]),
                box_mesh([0, 0, 2], [2, 2, 5])])
    # analytic: 12 + 12 - overlap(2*2*1 = 4) = 20
    prep = prepare_operand({"V": A[0], "F": A[1]}, 0, "input A")
    ok = check("t10 slab overlap folded",
               prep["report"]["status"] == "normalized",
               f"status={prep['report']['status']}")
    ok &= check("t10 slab union volume", abs(prep["volume"] - 20.0) < 1e-9,
                f"vol={prep['volume']}")
    ok &= check("t10 slab one shell",
                connected_shells(prep["F"]) == 1,
                f"shells={connected_shells(prep['F'])}")
    # ... and NOT the naive sum (24.0) the unfolded code produced
    ok &= check("t10 slab not double-counted",
                abs(prep["volume"] - 24.0) > 1.0,
                f"vol={prep['volume']}")
    # end-to-end through arrange: union with a disjoint box
    B = box_mesh([10, 0, 0], [11, 1, 1])
    m = arrange({"V": A[0], "F": A[1]}, {"V": B[0], "F": B[1]}, "union")
    ok &= check("t10 arrange union vol", abs(vol(m) - 21.0) < 1e-9,
                f"vol={vol(m)}")
    return ok


# --------------------------------------------------------------------------
# 11: genuinely disjoint shells take the passthrough (status pinned)
# --------------------------------------------------------------------------

def t11_disjoint_passthrough_status():
    A = concat([box_mesh([0, 0, 0], [1, 1, 1]),
                box_mesh([5, 0, 0], [6, 1, 1])])
    prep = prepare_operand({"V": A[0], "F": A[1]}, 0, "input A")
    ok = check("t11 passthrough status",
               prep["report"]["status"] == "disjoint passthrough",
               f"status={prep['report']['status']}")
    ok &= check("t11 passthrough bit-identical",
                np.array_equal(np.asarray(prep["V"]), A[0]) and
                np.array_equal(np.asarray(prep["F"]), A[1]))
    ok &= check("t11 passthrough volume", abs(prep["volume"] - 2.0) < 1e-9,
                f"vol={prep['volume']}")
    return ok


# --------------------------------------------------------------------------
# 12: negative signed result volume is refused by the semantic leg
# --------------------------------------------------------------------------

def t12_negative_signed_volume_refused():
    from brepkernel.normalize import check_boolean_semantics
    A = box_mesh([0, 0, 0], [1, 1, 1])
    B = box_mesh([2, 0, 0], [3, 1, 1])
    dA = {"V": A[0], "F": A[1], "volume": 1.0}
    dB = {"V": B[0], "F": B[1], "volume": 1.0}
    bad = box_mesh([0, 0, 0], [1, 1, 1])
    badF = bad[1][:, ::-1].copy()          # inside-out: negative signed vol
    assert signed_volume(bad[0], badF) < 0
    ok = True
    for op in ("union", "intersection", "difference"):
        def run(op=op):
            check_boolean_semantics(bad[0], badF, dA, dB, op, 1.0)
        ok &= expect_refusal(f"t12 negative-volume {op} refused",
                             "negative signed", run)
    return ok


def main():
    ok = True
    ok &= t1_two_overlapping_shells()
    ok &= t2_three_overlapping_shells()
    ok &= t3_duplicate_coincident()
    ok &= t4_disjoint_assembly_bit_identical()
    ok &= t5_hollow_cube()
    ok &= t6_welded_tetras()
    ok &= t7_self_intersecting_shell()
    ok &= t8_shell_permutation()
    ok &= t9_nested_island()
    ok &= t10_interpenetrating_slab_folded()
    ok &= t11_disjoint_passthrough_status()
    ok &= t12_negative_signed_volume_refused()
    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
