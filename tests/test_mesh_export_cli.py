"""B-rep -> glTF export and the brepkernel command line.

Exporter: watertight, outward, hard-edged meshes whose volume converges to
the exact volume; correct frame conversion (millimetres +Z up -> meters +Y
up). CLI: stable exit codes (0 done, 2 refused with no output file, 1
error), JSON reports, exact placement of B (rigid, uniform scale, mirror,
general affine).
"""
import json
import math
import os
import subprocess
import sys

import numpy as np
import pytest

import brepkernel._occt_compat  # noqa: F401
from brepkernel.cli import _volume, main, place_shape, read_shape, write_shape
from brepkernel.mesh_export import MeshExportError, read_glb, tessellate, write_glb
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCP.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCone,
                             BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeSphere,
                             BRepPrimAPI_MakeTorus)
from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

SRC = os.path.join(os.path.dirname(__file__), "..", "src")

SHAPES = {
    "box": lambda: BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), gp_Pnt(1, 2, 3)).Shape(),
    "sphere": lambda: BRepPrimAPI_MakeSphere(1.0).Shape(),          # poles
    "cone": lambda: BRepPrimAPI_MakeCone(1.0, 0.0, 2.0).Shape(),     # apex
    "torus": lambda: BRepPrimAPI_MakeTorus(2.0, 0.5).Shape(),        # two seams
    "cavity": lambda: BRepAlgoAPI_Cut(
        BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), gp_Pnt(4, 4, 4)).Shape(),
        BRepPrimAPI_MakeSphere(gp_Pnt(2, 2, 2), 1.0).Shape()).Shape(),
}


def _outward(m):
    P, N, T = m.positions, m.normals, m.indices
    g = np.cross(P[T[:, 1]] - P[T[:, 0]], P[T[:, 2]] - P[T[:, 0]])
    g /= np.linalg.norm(g, axis=1)[:, None]
    vn = N[T].mean(axis=1)
    vn /= np.linalg.norm(vn, axis=1)[:, None]
    return float(np.einsum("ij,ij->i", g, vn).min())


@pytest.mark.parametrize("name", sorted(SHAPES))
def test_watertight_outward_and_convergent(name):
    s = SHAPES[name]()
    exact = _volume(s)
    coarse = tessellate(s, exact_volume=exact)
    fine = tessellate(s, linear_deflection=coarse.linear_deflection / 8,
                      angular_deflection=0.1, exact_volume=exact)
    for m in (coarse, fine):
        assert m.watertight, m.summary()
        assert _outward(m) > 0.5
        assert m.mesh_volume > 0
    err_c = coarse.summary()["volume_relative_error"]
    err_f = fine.summary()["volume_relative_error"]
    assert err_f <= max(err_c, 1e-12)
    assert err_f < 2e-3


def test_hard_edges_between_faces():
    m = tessellate(SHAPES["box"]())
    assert m.vertex_count == 24          # 6 faces x 4 corners: no shared normals
    assert len(m.face_ranges) == 6


def test_mm_z_up_to_meters_y_up():
    box = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), gp_Pnt(200, 100, 50)).Shape()   # mm, z tallest = 50
    m = tessellate(box, scale=0.001, from_z_up=True, exact_volume=_volume(box))
    lo, hi = m.positions.min(axis=0), m.positions.max(axis=0)
    assert np.allclose(hi - lo, [0.2, 0.05, 0.1])       # height (old z) is now y
    assert m.watertight and abs(m.mesh_volume - 0.001) < 1e-12


def test_glb_roundtrip_and_face_ranges(tmp_path):
    p = str(tmp_path / "t.glb")
    m = write_glb(SHAPES["cavity"](), p)
    doc, P, N, T = read_glb(p)
    assert doc["asset"]["version"] == "2.0"
    assert P.shape == m.positions.shape and T.shape == m.indices.shape
    faces = doc["meshes"][0]["extras"]["brepkernel"]["faces"]
    assert sum(c for _f, c in faces) == T.size
    assert np.allclose(P, m.positions, atol=1e-6)


def test_open_shell_is_refused(tmp_path):
    from OCP.TopAbs import TopAbs_FACE
    from OCP.TopExp import TopExp_Explorer
    face = TopExp_Explorer(SHAPES["box"](), TopAbs_FACE).Current()
    with pytest.raises(MeshExportError, match="watertight"):
        write_glb(face, str(tmp_path / "f.glb"))
    assert not (tmp_path / "f.glb").exists()


@pytest.mark.parametrize("m12,kind,vol", [
    ([-1, 0, 0, 5, 0, 1, 0, 0, 0, 0, 1, 0], "rigid_mirror", 6.0),
    ([0, -2, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0], "similarity", 48.0),
    ([2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0.5, 0], "general_affine", 6.0),
])
def test_placement(m12, kind, vol):
    s, k = place_shape(SHAPES["box"](), m12)
    assert k == kind
    assert abs(_volume(s) - vol) < 1e-9 * vol


def _run(*args):
    env = dict(os.environ, PYTHONPATH=os.path.abspath(SRC))
    p = subprocess.run([sys.executable, "-m", "brepkernel", *args], capture_output=True,
                       text=True, env=env, timeout=600)
    return p.returncode, (json.loads(p.stdout) if p.stdout.strip().startswith("{") else None), p


def test_cli_import_boolean_mesh(tmp_path):
    os.chdir(tmp_path)
    write_shape(BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), gp_Pnt(200, 100, 50)).Shape(), "plate.step")
    write_shape(BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 1, 0)), 0.02, 0.2).Shape(),
                "peg.brep")
    code, rep, _ = _run("import", "plate.step", "-o", "plate.brep", "--mesh", "plate.glb")
    assert code == 0 and rep["frame"] == {"scale": 0.001, "up": "z", "target": "meters, +Y up"}
    assert abs(rep["output"]["volume"] - 0.001) < 1e-12 and rep["mesh"]["watertight"]
    code, rep, _ = _run("boolean", "plate.brep", "peg.brep", "--op", "difference",
                        "-o", "hole.brep", "--mesh", "hole.glb",
                        "--b-transform", "1,0,0,0.1,0,1,0,-0.05,0,0,1,-0.05")
    assert code == 0 and rep["status"] == "accepted"
    expect = 0.001 - math.pi * 0.02 ** 2 * 0.05
    assert abs(rep["result"]["volume"] - expect) < 1e-12
    assert rep["mesh"]["watertight"] and os.path.exists("hole.glb")
    shape, _ = read_shape("hole.brep")
    assert abs(_volume(shape) - expect) < 1e-12


def test_cli_refusal_writes_nothing(tmp_path):
    os.chdir(tmp_path)
    write_shape(BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, -2, 0), gp_Dir(0, 1, 0)), 1, 4).Shape(), "a.brep")
    write_shape(BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-2, 0, 0), gp_Dir(1, 0, 0)), 1, 4).Shape(), "b.brep")
    code, rep, _ = _run("boolean", "a.brep", "b.brep", "--op", "intersection", "-o", "r.brep",
                        "--mesh", "r.glb", "--advisory")
    assert code == 2 and rep["status"] == "refused" and rep["refusal"]["kind"]
    assert rep["advisory"]["certified"] is False
    assert not os.path.exists("r.brep") and not os.path.exists("r.glb")


def test_cli_errors_are_structured(tmp_path):
    os.chdir(tmp_path)
    code, rep, _ = _run("boolean", "missing.brep", "missing2.brep", "--op", "union", "-o", "x.brep")
    assert code == 1 and rep["status"] == "error" and "not found" in rep["error"]
    assert main(["frobnicate"]) == 1


def test_cli_report_file_and_clean_stdout(tmp_path):
    os.chdir(tmp_path)
    write_shape(SHAPES["box"](), "b.step")
    code, rep, p = _run("import", "b.step", "-o", "b.brep", "--report", "r.json")
    assert code == 0 and rep is None and p.stdout.strip() == ""
    assert json.load(open("r.json"))["status"] == "ok"


def test_cli_primitive(tmp_path):
    os.chdir(tmp_path)
    code, rep, _ = _run("primitive", "cylinder", "-o", "c.brep", "--radius", "0.25", "--height", "2",
                        "--axis", "x", "--at", "1", "0", "0", "--mesh", "c.glb")
    assert code == 0 and abs(rep["output"]["volume"] - math.pi * 0.25 ** 2 * 2) < 1e-12
    assert rep["mesh"]["watertight"]
    shape, _ = read_shape("c.brep")
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib
    box = Bnd_Box()
    BRepBndLib.Add_s(shape, box)
    lo, hi = box.CornerMin(), box.CornerMax()
    assert abs(lo.X() - 0.0) < 1e-6 and abs(hi.X() - 2.0) < 1e-6       # centred on --at, along x
    code, rep, _ = _run("primitive", "box", "-o", "b.brep", "--size", "1", "0", "1")
    assert code == 1 and "positive" in rep["error"]
