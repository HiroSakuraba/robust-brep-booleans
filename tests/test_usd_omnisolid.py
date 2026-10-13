"""omniSolid BrepArray (.usda) writer and reader: OCCT only, no usd-brep.

Round trips OCCT -> .usda -> OCCT must reproduce the volume and stay
valid, and SMLib-written fixtures (tests/data/usd) must rebuild into valid
OCCT solids with the right volume.
"""
import math
import os

import pytest

import brepkernel._occt_compat  # noqa: F401
from brepkernel.usd import (OmniSolidReadError, OmniSolidWriteError, brep_array_prims,
                            read_omnisolid_usda, write_omnisolid_usda)
from OCP.BRep import BRep_Builder
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCP.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire
from OCP.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCone,
                             BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeSphere,
                             BRepPrimAPI_MakeTorus)
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_SOLID
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS, TopoDS_Compound
from OCP.gp import gp_Ax1, gp_Ax2, gp_Dir, gp_Elips, gp_Pnt, gp_Trsf, gp_Vec

DATA = os.path.join(os.path.dirname(__file__), "data", "usd")


def vol(s):
    # Shapes rebuilt from USD carry fitted UV curves; OCCT's default
    # integration accuracy is not enough for 1e-9 comparisons on them.
    from brepkernel.assembly import _all_faces_analytic, _shape_volume
    return _shape_volume(s, 1e-12 if _all_faces_analytic(s) else 1e-10)


def nsolids(s):
    n, ex = 0, TopExp_Explorer(s, TopAbs_SOLID)
    while ex.More():
        n += 1
        ex.Next()
    return n


def box(x0, y0, z0, x1, y1, z1):
    return BRepPrimAPI_MakeBox(gp_Pnt(x0, y0, z0), gp_Pnt(x1, y1, z1)).Shape()


def moved(s, ax=(1, 1, 1), ang=0.7, off=(0, 0, 0)):
    t = gp_Trsf()
    t.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(*ax)), ang)
    s = BRepBuilderAPI_Transform(s, t, True).Shape()
    t = gp_Trsf()
    t.SetTranslation(gp_Vec(*off))
    return BRepBuilderAPI_Transform(s, t, True).Shape()


def mirrored(s):
    t = gp_Trsf()
    t.SetMirror(gp_Ax2(gp_Pnt(0.1, 0.2, 0.3), gp_Dir(1, 0.2, 0.1)))
    return BRepBuilderAPI_Transform(s, t, True).Shape()


def filleted_box():
    b = box(0, 0, 0, 2, 2, 2)
    f = BRepFilletAPI_MakeFillet(b)
    ex = TopExp_Explorer(b, TopAbs_EDGE)
    while ex.More():
        f.Add(0.3, TopoDS.Edge(ex.Current()))
        ex.Next()
    return f.Shape()


def loft():
    ts = BRepOffsetAPI_ThruSections(True, False)
    for z, a, b in [(0, 2, 1), (1, 1.5, 1.2), (2, 2.2, 0.8)]:
        e = BRepBuilderAPI_MakeEdge(gp_Elips(gp_Ax2(gp_Pnt(0, 0, z), gp_Dir(0, 0, 1)), a, b)).Edge()
        ts.AddWire(BRepBuilderAPI_MakeWire(e).Wire())
    ts.Build()
    return ts.Shape()


def holebox():
    return BRepAlgoAPI_Cut(box(0, 0, 0, 2, 2, 2), BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(1, 1, -1), gp_Dir(0, 0, 1)), 0.5, 4).Shape()).Shape()


def compound(*shapes):
    c = TopoDS_Compound()
    b = BRep_Builder()
    b.MakeCompound(c)
    for s in shapes:
        b.Add(c, s)
    return c


ROUND_TRIP = {
    "box": lambda: box(0, 0, 0, 1, 2, 3),
    "holebox": holebox,
    "cylinder": lambda: BRepPrimAPI_MakeCylinder(1.0, 2.0).Shape(),
    "sphere": lambda: BRepPrimAPI_MakeSphere(1.0).Shape(),
    "torus": lambda: BRepPrimAPI_MakeTorus(3.0, 1.0).Shape(),
    "cone": lambda: BRepPrimAPI_MakeCone(1.0, 0.4, 2.0).Shape(),
    "fillet": filleted_box,
    "loft_bspline": loft,
    "cavity": lambda: BRepAlgoAPI_Cut(box(0, 0, 0, 4, 4, 4), BRepPrimAPI_MakeSphere(
        gp_Pnt(2, 2, 2), 1.0).Shape()).Shape(),
    "far_rotated": lambda: moved(holebox(), off=(1000, -2000, 500)),
    "mirrored_fillet": lambda: mirrored(filleted_box()),
    "mirrored_torus_cone": lambda: mirrored(BRepAlgoAPI_Fuse(
        BRepPrimAPI_MakeTorus(3.0, 1.0).Shape(),
        BRepPrimAPI_MakeCone(gp_Ax2(gp_Pnt(3, 0, -2), gp_Dir(0, 0, 1)), 1.0, 0.3, 4.0).Shape()
    ).Shape()),
    "three_lumps": lambda: compound(
        box(0, 0, 0, 1, 1, 1),
        BRepAlgoAPI_Cut(box(3, 0, 0, 6, 3, 3),
                        BRepPrimAPI_MakeSphere(gp_Pnt(4.5, 1.5, 1.5), 1).Shape()).Shape(),
        BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 5, 0), gp_Dir(0, 0, 1)), 1, 2).Shape()),
}


@pytest.mark.parametrize("name", sorted(ROUND_TRIP))
@pytest.mark.parametrize("per_solid", [False, True])
def test_round_trip(tmp_path, name, per_solid):
    shape = ROUND_TRIP[name]()
    path = str(tmp_path / f"{name}.usda")
    prim = write_omnisolid_usda(shape, path, one_brep_per_solid=per_solid)
    assert prim == "/World/Brep0"
    back = read_omnisolid_usda(path)
    assert BRepCheck_Analyzer(back).IsValid()
    assert nsolids(back) == nsolids(shape)
    v0, v1 = vol(shape), vol(back)
    assert abs(v1 - v0) <= 1e-9 * v0, (v0, v1)


def test_analytic_surfaces_are_written_as_analytic(tmp_path):
    path = str(tmp_path / "c.usda")
    fused = BRepAlgoAPI_Fuse(BRepPrimAPI_MakeCylinder(1.0, 2.0).Shape(),
                             BRepPrimAPI_MakeTorus(3.0, 0.5).Shape()).Shape()
    write_omnisolid_usda(compound(fused, BRepPrimAPI_MakeSphere(gp_Pnt(9, 0, 0), 1).Shape()), path)
    types = set(brep_array_prims(path)["/World/Brep0"]["face:surfaceType"])
    assert {"BrepSurfaceConeAPI", "BrepSurfaceTorusAPI", "BrepSurfaceSphereAPI",
            "BrepSurfacePlaneAPI"} <= types
    assert "BrepSurfaceNurbAPI" not in types


def test_output_is_deterministic(tmp_path):
    a, b = str(tmp_path / "a.usda"), str(tmp_path / "b.usda")
    s = filleted_box()
    write_omnisolid_usda(s, a)
    write_omnisolid_usda(s, b)
    assert open(a).read() == open(b).read()


def test_no_em_dash_or_non_ascii(tmp_path):
    p = str(tmp_path / "x.usda")
    write_omnisolid_usda(holebox(), p)
    open(p, encoding="ascii").read()


def test_writer_refuses_non_solid(tmp_path):
    ex = TopExp_Explorer(box(0, 0, 0, 1, 1, 1), TopAbs_FACE)
    with pytest.raises(OmniSolidWriteError):
        write_omnisolid_usda(ex.Current(), str(tmp_path / "f.usda"))


def test_reader_refuses_usdc(tmp_path):
    p = tmp_path / "x.usda"
    p.write_bytes(b"PXR-USDC" + b"\0" * 64)
    with pytest.raises(OmniSolidReadError, match="usdc"):
        read_omnisolid_usda(str(p))


def test_reader_needs_prim_path_when_ambiguous(tmp_path):
    a = str(tmp_path / "a.usda")
    write_omnisolid_usda(box(0, 0, 0, 1, 1, 1), a, prim_name="A")
    body = open(a).read()
    two = body.replace('def Xform "World"\n{', 'def Xform "World"\n{\n' + body.split(
        'def Xform "World"\n{')[1].rsplit("}", 1)[0].replace('"A"', '"B"'), 1)
    p = str(tmp_path / "two.usda")
    open(p, "w").write(two)
    assert sorted(brep_array_prims(p)) == ["/World/A", "/World/B"]
    with pytest.raises(OmniSolidReadError, match="prim_path"):
        read_omnisolid_usda(p)
    assert abs(vol(read_omnisolid_usda(p, "/World/B")) - 1.0) < 1e-12


SPHERE = 4.0 / 3.0 * math.pi * 1.2 ** 3


def _cylcyl_reference():
    return vol(BRepAlgoAPI_Fuse(
        BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 1.0, 2.0).Shape(),
        BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0.5, 0, 1), gp_Dir(0, 0, 1)), 0.6, 2.0).Shape()
    ).Shape())


@pytest.mark.parametrize("fname,expected,solids", [
    ("sm_boxsph.usda", 8 - SPHERE / 8, 1),
    ("sm_corner_union.usda", 8 + 7 * SPHERE / 8, 1),
    ("sm_cylcyl.usda", None, 1),
    ("sm_split.usda", 12.0, 2),
    ("multi_sm.usda", 1 + (27 - 4.0 / 3.0 * math.pi) + 2 * math.pi, 3),
])
def test_reads_smlib_output(fname, expected, solids):
    if expected is None:
        expected = _cylcyl_reference()
    s = read_omnisolid_usda(os.path.join(DATA, fname))
    assert BRepCheck_Analyzer(s).IsValid()
    assert nsolids(s) == solids
    assert abs(vol(s) - expected) <= 1e-7 * expected, (vol(s), expected)
