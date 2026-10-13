"""T2 and T6: refusals that two independent engines showed to be wrong.

Each case used to refuse; OCCT's own Boolean and NVIDIA SMLib agreed on the
answer, and the cause was found in brepkernel:

* T6 (SectionCompletenessMismatch): the completeness probe paired each raw
  IntTools p-curve with the wrong face (IntTools_FaceFace reorders faces
  internally), so its in-trim test looked at the wrong UV points. The
  pairing is now decided geometrically.
* T2 (InsufficientPatchWitnesses): thin trimmed slivers fill a few percent
  of their UV box and the fixed witness sample missed them; the search
  now escalates to finer grids. The 3-witness rule is unchanged.

Specs come from tools/review_probes/corpus_manifest.json (fuzz seeds 5 and
7); the expected volumes are OCCT BRepAlgoAPI's, which SMLib matched to
better than 2e-6. Every accepted result must also pass an independent
point-membership audit.
"""
import os
import sys

import numpy as np
import pytest

from brepkernel import boolean_brep
from brepkernel.assembly import _all_faces_analytic, _shape_volume
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeSphere
from OCP.gp import gp_Ax1, gp_Ax2, gp_Dir, gp_Pnt, gp_Trsf, gp_Vec

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools", "review_probes"))
import arbiter  # noqa: E402


def build(spec):
    k = spec["prim"]
    if k == "box":
        s = BRepPrimAPI_MakeBox(gp_Pnt(*spec["p0"]), gp_Pnt(*spec["p1"])).Shape()
    elif k == "cyl":
        s = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(*spec["p"]), gp_Dir(*spec["d"])),
                                     spec["r"], spec["h"]).Shape()
    else:
        s = BRepPrimAPI_MakeSphere(gp_Pnt(*spec["c"]), spec["r"]).Shape()
    ax, ang = spec["xform"]["rotate"]
    t = gp_Trsf()
    t.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(*ax)), ang)
    s = BRepBuilderAPI_Transform(s, t, True).Shape()
    t2 = gp_Trsf()
    t2.SetTranslation(gp_Vec(*spec["xform"]["translate"]))
    return BRepBuilderAPI_Transform(s, t2, True).Shape()


CASES = [
    ("T6 seed5 trial3", "union", 2.7548703512592163,
     {"prim": "box", "p0": [-0.5, -0.75, -0.5], "p1": [0.5, 0.75, 0.5],
      "xform": {"rotate": [[1.0, 0.0, 0.0], 3.141592653589793], "translate": [-0.5, -0.5, 0.5]}},
     {"prim": "cyl", "p": [0, 0, -1.0], "d": [0, 0, 1], "r": 0.5, "h": 2.0,
      "xform": {"rotate": [[0.0, 1.0, 0.0], 1.5707963267948966], "translate": [0.5, -0.0, 0.5]}}),
    ("T6 seed5 trial4", "difference", 3.1259987948890036,
     {"prim": "box", "p0": [-0.75, -0.75, -0.75], "p1": [0.75, 0.75, 0.75],
      "xform": {"rotate": [[1.0, 0.0, 0.0], 3.141592653589793], "translate": [-0.5, 0.0, -0.5]}},
     {"prim": "cyl", "p": [0, 0, -0.75], "d": [0, 0, 1], "r": 0.5, "h": 1.5,
      "xform": {"rotate": [[0.0, 1.0, 0.0], 4.71238898038469], "translate": [0.5, -0.5, -0.0]}}),
    ("T6 seed7 trial2", "intersection", 0.7082610948700642,
     {"prim": "cyl", "p": [0, 0, -0.9281992923884972], "d": [0, 0, 1], "r": 0.46055958273827124,
      "h": 1.8563985847769944,
      "xform": {"rotate": [[0.05189409494388281, -0.9967293025347472, 0.061949175769753444],
                           4.6606843177567026],
                "translate": [-0.4902052739243452, 0.04937258565178648, 0.00932668356041999]}},
     {"prim": "sphere", "c": [0, 0, 0], "r": 0.9228036260157284,
      "xform": {"rotate": [[0.815241202649682, 0.3106480625485527, -0.4887530692866599],
                           2.4355624373344904],
                "translate": [-0.21235638449015198, -0.41976032511545774, 0.37960572458289077]}}),
    ("T2 seed7 trial31", "difference", 1.5683203843349014,
     {"prim": "box", "p0": [-0.7163526457167202, -0.5526306847481507, -0.8108886542812301],
      "p1": [0.7163526457167202, 0.5526306847481507, 0.8108886542812301],
      "xform": {"rotate": [[0.08719635092216456, 0.7868692513218176, -0.6109284554758444],
                           0.830781766875637],
                "translate": [-0.3788733247982318, -0.23529836518567376, -0.06298113093528113]}},
     {"prim": "cyl", "p": [0, 0, -0.5670507387277491], "d": [0, 0, 1], "r": 0.8567979469740528,
      "h": 1.1341014774554983,
      "xform": {"rotate": [[0.1852117156482754, -0.5126900731906607, -0.8383588189065451],
                           2.5337158045740846],
                "translate": [-0.5960679097036258, 0.3443608344788055, -0.3054289785370093]}}),
    ("T2 seed7 trial35", "union", 3.701824802377997,
     {"prim": "box", "p0": [-0.6568878819481545, -0.8191871725649438, -0.8202396580415123],
      "p1": [0.6568878819481545, 0.8191871725649438, 0.8202396580415123],
      "xform": {"rotate": [[0.04078148677359839, 0.14117043948419658, 0.9891449728691815],
                           4.7754849100741055],
                "translate": [-0.5121372038246175, -0.2820881253900782, 0.37517134238660155]}},
     {"prim": "cyl", "p": [0, 0, -0.6718901928950246], "d": [0, 0, 1], "r": 0.581733733522104,
      "h": 1.3437803857900492,
      "xform": {"rotate": [[-0.8183706519120203, 0.3790329297010845, -0.431976289038308],
                           2.4642705180243887],
                "translate": [-0.41595172611644315, -0.5017297725627801, 0.08143242813745921]}}),
]


@pytest.mark.parametrize("name,op,expect,sa,sb", CASES, ids=[c[0] for c in CASES])
def test_now_accepted_and_correct(name, op, expect, sa, sb):
    a, b = build(sa), build(sb)
    out, rep = boolean_brep(a, b, op)
    v = abs(_shape_volume(out, 1e-12 if _all_faces_analytic(out) else 1e-10))
    assert abs(v - expect) <= 1e-8 * expect, (v, expect)
    audit = arbiter.membership_audit(a, b, out, op, np.random.default_rng(7), n=200)
    assert audit["kernel_errors"] == 0, audit
    assert rep["stages"]["verification"]["tolerance_ratio"] <= 2.0


def test_pairing_is_geometric():
    """The probe's p-curve pairing must survive IntTools swapping faces."""
    from brepkernel import intersection as ix

    class C2:
        def __init__(self, f):
            self.f = f

        def Value(self, t):
            from OCP.gp import gp_Pnt2d
            return gp_Pnt2d(*self.f(t))

    class C3:
        def Value(self, t):
            return gp_Pnt(t, 0.0, 0.0)

    class Face:
        def __init__(self, shape):
            self.face = shape

    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace
    from OCP.gp import gp_Pln, gp_Ax3
    # plane z=0 (u=x, v=y) and plane y=0 (u=x, v=-z or similar): the curve
    # x-axis lies on both
    fa = Face(BRepBuilderAPI_MakeFace(gp_Pln(gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0))),
                                      -5, 5, -5, 5).Face())
    fb = Face(BRepBuilderAPI_MakeFace(gp_Pln(gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 1, 0), gp_Dir(0, 0, 1))),
                                      -5, 5, -5, 5).Face())
    on_a = C2(lambda t: (t, 0.0))           # (u, v) on fa
    on_b = C2(lambda t: (0.0, t))           # on fb: u along z, v along x
    ca, cb, err = ix._pair_raw_pcurves(C3(), on_b, on_a, fa, fb, 0.0, 1.0)
    assert ca is on_a and cb is on_b and err < 1e-12
