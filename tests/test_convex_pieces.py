"""Convex collision pieces for extruded parts (brepkernel.convex_pieces).

Each piece is checked independently of the certificate: OCCT's solid
classifier must put every piece vertex inside or on the part, in the part's
own frame, and the pieces' total volume must match the part within the
polygon error.
"""
import json
import math
import os
import subprocess
import sys

import pytest

pytest.importorskip("scipy")
pytest.importorskip("shapely")

import brepkernel._occt_compat  # noqa: F401
from brepkernel.cli import _volume, write_shape
from brepkernel.convex_pieces import convex_pieces, merge_convex
from brepkernel.mujoco_fit.geom2d import exact_area, frac_points
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCone, BRepPrimAPI_MakeCylinder
from OCP.gp import gp_Ax1, gp_Ax2, gp_Dir, gp_Pnt, gp_Trsf
from OCP.TopAbs import TopAbs_IN, TopAbs_ON

SRC = os.path.join(os.path.dirname(__file__), "..", "src")


def box(*c):
    return BRepPrimAPI_MakeBox(gp_Pnt(*c[:3]), gp_Pnt(*c[3:])).Shape()


def cyl(x, y, z, r, h, d=(0, 0, 1)):
    return BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(x, y, z), gp_Dir(*d)), r, h).Shape()


PARTS = {
    "through_hole": lambda: BRepAlgoAPI_Cut(box(-.025, -.025, 0, .025, .025, .01), cyl(0, 0, -.001, .008, .012)).Shape(),
    "two_holes": lambda: BRepAlgoAPI_Cut(BRepAlgoAPI_Cut(box(-.025, -.025, 0, .025, .025, .01),
                                                         cyl(-.01, 0, -.001, .005, .012)).Shape(),
                                         cyl(.01, 0, -.001, .005, .012)).Shape(),
    "l_shape": lambda: BRepAlgoAPI_Fuse(box(0, 0, 0, .04, .01, .01), box(0, 0, 0, .01, .04, .01)).Shape(),
    "boss": lambda: BRepAlgoAPI_Fuse(box(-.02, -.02, 0, .02, .02, .005), cyl(0, 0, .005, .006, .01)).Shape(),
    "y_axis_ring": lambda: BRepAlgoAPI_Cut(cyl(0, 0, 0, .02, .01, (0, 1, 0)), cyl(0, -.001, 0, .01, .012, (0, 1, 0))).Shape(),
}


def test_merge_convex_square_from_triangles():
    tri = [frac_points([[0, 0], [1, 0], [1, 1]]), frac_points([[0, 0], [1, 1], [0, 1]])]
    out = merge_convex(tri)
    assert len(out) == 1 and len(out[0]) == 4 and exact_area(out[0]) == 1


def _inside(shape, p, tol=1e-9):
    c = BRepClass3d_SolidClassifier(shape, gp_Pnt(*p), tol)
    return c.State() in (TopAbs_IN, TopAbs_ON)


def _prism_volume(piece):
    b, t = piece["bottom"], piece["top"]
    import numpy as np
    b, t = np.array(b), np.array(t)
    d = t[0] - b[0]
    n = np.zeros(3)
    for i in range(len(b)):
        n += np.cross(b[i], b[(i + 1) % len(b)])
    return abs(float(np.dot(n / 2, d)))


@pytest.mark.parametrize("name", sorted(PARTS))
def test_pieces_inside_part_and_volume(name):
    shape = PARTS[name]()
    rep = convex_pieces(shape, max_error=5e-5)
    assert rep["status"] == "certified", rep.get("refusal")
    for piece in rep["pieces"]:
        for p in piece["bottom"] + piece["top"]:
            assert _inside(shape, p), (name, p)
    vol = sum(_prism_volume(p) for p in rep["pieces"])
    exact = _volume(shape)
    assert vol <= exact * (1 + 1e-12)                     # never more material than the part
    worst = max(s["deviation"]["reference_outside_sim"] for s in rep["slabs"])
    assert worst <= 5e-5 + rep["resolution"]
    assert exact - vol < 0.01 * exact
    assert rep["frame_exact"]


def test_slabs_and_topology():
    rep = convex_pieces(PARTS["boss"](), max_error=5e-5)
    assert [s["z"] for s in rep["slabs"]] == [[0.0, 0.005], [0.005, 0.015]]
    rep = convex_pieces(PARTS["two_holes"](), max_error=5e-5)
    assert (rep["slabs"][0]["outer_loops"], rep["slabs"][0]["hole_loops"]) == (1, 2)


def test_tilted_part_and_refusal():
    t = gp_Trsf()
    t.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 1, 0)), math.radians(30))
    shape = BRepBuilderAPI_Transform(PARTS["through_hole"](), t, True).Shape()
    rep = convex_pieces(shape, max_error=5e-5)
    assert rep["status"] == "certified" and not rep["frame_exact"] and rep["warnings"]
    for piece in rep["pieces"]:
        for p in piece["bottom"] + piece["top"]:
            assert _inside(shape, p, 1e-8)
    cone = BRepPrimAPI_MakeCone(gp_Ax2(gp_Pnt(0, 0, -1e-5), gp_Dir(0, 0, 1)), .008, .009, .01002).Shape()
    rep = convex_pieces(BRepAlgoAPI_Cut(box(-.025, -.025, 0, .025, .025, .01), cone).Shape())
    assert rep["status"] == "refused" and rep["refusal"]["kind"] == "NotAnExtrusion" and "pieces" not in rep


def test_cli(tmp_path):
    write_shape(PARTS["through_hole"](), str(tmp_path / "plate.brep"))
    cone = BRepPrimAPI_MakeCone(gp_Ax2(gp_Pnt(0, 0, -1e-5), gp_Dir(0, 0, 1)), .008, .009, .01002).Shape()
    write_shape(BRepAlgoAPI_Cut(box(-.025, -.025, 0, .025, .025, .01), cone).Shape(), str(tmp_path / "cone.brep"))
    env = dict(os.environ, PYTHONPATH=os.path.abspath(SRC))

    def run(*a):
        p = subprocess.run([sys.executable, "-m", "brepkernel", *a], capture_output=True, text=True, env=env,
                           cwd=tmp_path, timeout=600)
        return p.returncode, json.loads(p.stdout)

    code, rep = run("convex-pieces", "plate.brep", "-o", "pieces.json", "--max-error", "5e-5")
    assert code == 0 and rep["status"] == "certified"
    data = json.loads((tmp_path / "pieces.json").read_text())
    assert data["units"] == "meters" and len(data["pieces"]) == rep["output"]["pieces"] > 4
    code, rep = run("convex-pieces", "cone.brep", "-o", "no.json")
    assert code == 2 and rep["certificate"]["refusal"]["kind"] == "NotAnExtrusion"
    assert not (tmp_path / "no.json").exists()
