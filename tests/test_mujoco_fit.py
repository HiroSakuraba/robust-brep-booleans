"""mjfit bridge: certify mjfit peg-in-hole models (mode A) and build them from CAD (mode C).

Geometry tests need scipy and shapely; model tests also need mjfit (and MuJoCo).
Every check that the certificate passes has a planted-error twin that it must refuse.
"""
import importlib.util
import math
import xml.etree.ElementTree as ET

import numpy as np
import pytest

pytest.importorskip("scipy")
pytest.importorskip("shapely")

from brepkernel.mujoco_fit.geom2d import (Arc, LoopRegion, OffsetRegion, PieceRegion, deviation,
                                          exact_overlap, exact_union_loops, frac_points, polygon_loop)
from brepkernel.mujoco_fit.reference import polygonize

H = 1e-6


def square(cx, cy, s):
    return np.array([[cx - s, cy - s], [cx + s, cy - s], [cx + s, cy + s], [cx - s, cy + s]])


# --- exact pieces --------------------------------------------------------------

def test_exact_overlap_and_union_loops():
    a = frac_points([[0, 0], [1, 0], [1, 1], [0, 1]])
    b = frac_points([[1, 0], [2, 0], [2, 1], [1, 1]])               # shares an edge
    c = frac_points([[0.5, 0], [1.5, 0], [1.5, 1], [0.5, 1]])
    assert not exact_overlap(a, b) and exact_overlap(a, c)
    assert len(exact_union_loops([a, b])) == 1                      # shared edge cancels
    # four squares around a central square hole: outer loop + hole loop
    ring = [frac_points(p) for p in ([[0, 0], [3, 0], [3, 1], [0, 1]], [[0, 2], [3, 2], [3, 3], [0, 3]],
                                     [[0, 1], [1, 1], [1, 2], [0, 2]], [[2, 1], [3, 1], [3, 2], [2, 2]])]
    assert len(exact_union_loops(ring)) == 2
    # T-junction: one long edge against two short ones still cancels
    t = [frac_points([[0, 0], [2, 0], [2, 1], [0, 1]]), frac_points([[0, 1], [1, 1], [1, 2], [0, 2]]),
         frac_points([[1, 1], [2, 1], [2, 2], [1, 2]])]
    assert len(exact_union_loops(t)) == 1
    # a gap of one part in 2^52 is still a gap
    eps = 2.0 ** -52
    g = [frac_points([[0, 0], [1, 0], [1, 1], [0, 1]]), frac_points([[1 + eps, 0], [2, 0], [2, 1], [1 + eps, 1]])]
    assert len(exact_union_loops(g)) == 2


# --- sampled regions -----------------------------------------------------------

def test_erosion_matches_known_shapes():
    # square of half-side 0.01 eroded by 0.001 is the square of half-side 0.009 (sharp corners)
    sq = LoopRegion([polygon_loop(square(0, 0, 0.01))], H)
    ero = OffsetRegion(sq, 0.001, "erode")
    dv = deviation(PieceRegion([square(0, 0, 0.009)], H, 1e-12), ero)
    assert max(dv.outside, dv.inside) <= dv.resolution and dv.grid_mismatch_far == 0
    # disk r=0.01 eroded by 0.002 is the disk r=0.008
    disk = LoopRegion([[Arc((0, 0), 0.01, 0, 2 * math.pi)]], H)
    ero = OffsetRegion(disk, 0.002, "erode")
    r = np.linalg.norm(ero.boundary(), axis=1)
    assert abs(r.min() - 0.008) < 1e-12 and abs(r.max() - 0.008) < 1e-12
    # reflex corner: L-shaped hole; the eroded boundary has an arc of radius d around it
    L = np.array([[0, 0], [0.02, 0], [0.02, 0.01], [0.01, 0.01], [0.01, 0.02], [0, 0.02]])
    ero = OffsetRegion(LoopRegion([polygon_loop(L)], H), 0.001, "erode")
    d = np.linalg.norm(ero.boundary() - [0.01, 0.01], axis=1)
    assert abs(d.min() - 0.001) < 2 * H


def test_deviation_finds_corner_play_and_covered_holes():
    disk = LoopRegion([[Arc((0, 0), 0.01, 0, 2 * math.pi)]], H)
    n = 32
    ang = np.arange(n) * 2 * math.pi / n
    poly = 0.01 / math.cos(math.pi / n) * np.column_stack([np.cos(ang), np.sin(ang)])
    dv = deviation(PieceRegion([poly], H, 1e-12), disk)
    assert abs(dv.outside - 0.01 * (1 / math.cos(math.pi / n) - 1)) <= dv.resolution
    assert dv.inside <= dv.resolution
    # a plate with a hole against a plate without one: whole-region mismatch
    outer = LoopRegion([polygon_loop(square(0, 0, 0.03))], H)
    holed = LoopRegion([polygon_loop(square(0, 0, 0.03)), [Arc((0, 0), 0.01, 0, 2 * math.pi)]], H)
    full = PieceRegion([square(0, 0, 0.03)], H, 1e-12)
    assert deviation(full, outer).grid_mismatch_far == 0
    assert deviation(full, holed).grid_mismatch_far > 0


def test_polygonize_keeps_free_space():
    # hole loop CCW (free inside): tangent polygon contains the circle
    hole = [Arc((0, 0), 0.008, 0, 2 * math.pi)]
    p = polygonize(hole, 5e-5, H)
    r = np.linalg.norm(p, axis=1)
    edges_mid = np.linalg.norm((p + np.roll(p, -1, axis=0)) / 2, axis=1)
    assert edges_mid.min() >= 0.008 * math.cos(math.pi / len(p)) - 1e-15
    assert r.max() - 0.008 <= 5e-5 + 1e-15
    # island loop CW (free outside): chords inside the island
    island = [Arc((0, 0), 0.001, 0, -2 * math.pi)]
    p = polygonize(island, 5e-5, H)
    assert np.linalg.norm(p, axis=1).max() <= 0.001 + 1e-15
    assert 0.001 - np.linalg.norm((p + np.roll(p, -1, axis=0)) / 2, axis=1).min() <= 5e-5 + 1e-15


# --- mjfit models ----------------------------------------------------------------

needs_mjfit = pytest.mark.skipif(importlib.util.find_spec("mjfit") is None, reason="mjfit not installed")
from brepkernel.mujoco_fit import certify_mjfit_config, certify_mjfit_peg, mjfit_peg_from_cad  # noqa: E402
from brepkernel.mujoco_fit.certify import certify_model  # noqa: E402


def _plant(xml, mesh, fn):
    """Rewrite one inline mesh's vertices with fn(array (n,3)) -> array."""
    root = ET.fromstring(xml)
    for m in root.iter("mesh"):
        if m.get("name") == mesh:
            v = np.array([float(x) for x in m.get("vertex").split()]).reshape(-1, 3)
            m.set("vertex", " ".join(repr(float(x)) for x in fn(v).ravel()))
    return ET.tostring(root, encoding="unicode")


def _move_body(xml, body, dx):
    root = ET.fromstring(xml)
    for b in root.iter("body"):
        if b.get("name") == body:
            p = [float(x) for x in b.get("pos").split()]
            p[0] += dx
            b.set("pos", " ".join(repr(x) for x in p))
    return ET.tostring(root, encoding="unicode")


CIRCLE = dict(hole={"type": "circle", "radius": 0.008, "n": 32}, clearance=1e-4, plate_size=0.05,
              plate_thickness=0.01, peg_length=0.03)


@pytest.fixture(scope="module")
def circle():
    xml, rep = certify_mjfit_peg(**CIRCLE)
    assert rep["status"] == "certified", rep.get("refusal")
    return xml, rep


@needs_mjfit
def test_circle_certificate_numbers(circle):
    _xml, rep = circle
    h = rep["resolution"]
    corner = 0.008 * (1 / math.cos(math.pi / 32) - 1)               # how far a corner sticks out
    plate = rep["slabs"]["plate_hole"]["deviation"]
    assert plate["sim_outside_reference"] <= h                       # hole never smaller than drawn
    assert abs(plate["reference_outside_sim"] - corner) <= h
    assert abs(rep["clearance"]["sim_min"] - 1e-4) <= rep["clearance"]["tolerance"]
    gap = 0.008 - (0.008 - 1e-4) / math.cos(math.pi / 32)            # sim peg corner to the real hole
    assert abs(rep["clearance"]["sim_peg_to_real_plate_min"] - gap) <= h
    assert rep["structure"]["exact"] and rep["slabs"]["plate_hole"]["boundary_loops"] == 2


def _centroid_scale(f):
    def fn(v):                                      # shrink or grow in x and y only
        c = v.mean(axis=0)
        out = v.copy()
        out[:, :2] = c[:2] + (v[:, :2] - c[:2]) * f
        return out
    return fn


def _rotate_hole_vertex(angle, delta):
    def fn(v):
        out = v.copy()
        r = np.linalg.norm(v[:, :2], axis=1)
        a = np.arctan2(v[:, 1], v[:, 0])
        pick = (r < 0.009) & (np.abs(a - angle) < 1e-9)
        out[pick, 0] = r[pick] * np.cos(a[pick] + delta)
        out[pick, 1] = r[pick] * np.sin(a[pick] + delta)
        return out
    return fn


def _pull_back(d):
    def fn(v):
        out = v.copy()
        r = np.linalg.norm(v[:, :2], axis=1)
        pick = r < 0.009
        out[pick, :2] *= ((r[pick] + d) / r[pick])[:, None]
        return out
    return fn


@needs_mjfit
def test_wrong_reference_is_a_shape_mismatch(circle):
    xml, rep = circle
    from brepkernel.mujoco_fit import certify as C
    from brepkernel.mujoco_fit.geom2d import DifferenceRegion
    h = rep["sample_spacing"]
    disk = LoopRegion([[Arc((0, 0), 0.0082, 0, 2 * math.pi)]], h)       # real hole 0.2 mm bigger
    refs = {"plate_hole": DifferenceRegion(LoopRegion([polygon_loop(square(0, 0, 0.025))], h), disk)}
    out = C.certify_model(xml, refs, 1e-4, h, loops={"plate_hole": (1, 1), "peg_insert": (1, 0)},
                          plate_budget=0.008 * (1 / math.cos(math.pi / 32) - 1))
    assert out["refusal"]["kind"] == "ShapeMismatch"


@needs_mjfit
@pytest.mark.parametrize("plant,kind", [
    (lambda x: _plant(x, "plate_3", _centroid_scale(1 - 1e-9)), "PieceGap"),
    (lambda x: _plant(x, "plate_3", _centroid_scale(1 + 1e-9)), "PieceOverlap"),
    (lambda x: _move_body(x, "peg", 1e-6), "ClearanceMismatch"),
    # a slit 8 um wide at the hole, closing toward the plate edge: same topology
    (lambda x: _plant(x, "plate_3", _rotate_hole_vertex(math.radians(3 * 11.25), 1e-3)), "BoundaryMismatch"),
    # a hole vertex moved outward in floating point lands a hair inside the neighbor
    (lambda x: _plant(x, "plate_3", _pull_back(1e-3)), "PieceOverlap"),
    (lambda x: _move_body(x, "peg", 2e-4), "SimInterference"),
    (lambda x: _plant(x, "plate_5", lambda v: v + np.where(v[:, 2:3] > 0, [1e-7, 0, 0], 0)), "NonPrismPiece"),
])
def test_planted_errors_are_refused(circle, plant, kind):
    xml, rep = circle
    from brepkernel.mujoco_fit import certify as C
    from brepkernel.mujoco_fit.geom2d import DifferenceRegion
    planted = plant(xml)
    h = rep["sample_spacing"]
    disk = LoopRegion([[Arc((0, 0), 0.008, 0, 2 * math.pi)]], h)
    refs = {"plate_hole": DifferenceRegion(LoopRegion([polygon_loop(square(0, 0, 0.025))], h), disk),
            "peg_insert": OffsetRegion(disk, 1e-4, "erode")}
    budget = 0.008 * (1 / math.cos(math.pi / 32) - 1)
    out = C.certify_model(planted, refs, 1e-4, h, loops={"plate_hole": (1, 1), "peg_insert": (1, 0)},
                          plate_budget=budget)
    assert out["status"] == "refused" and out["refusal"]["kind"] == kind, out.get("refusal")


@needs_mjfit
@pytest.mark.parametrize("name", ["circle", "db9", "hdmi", "rj45", "star", "usb_a", "usb_c", "xlr3",
                                  "rounded_convex"])
def test_bundled_peg_configs_certify(name):
    xml, rep = certify_mjfit_config(name)
    assert rep["status"] == "certified", rep.get("refusal")
    assert all(v["deviation"] is None or v["deviation"]["grid_mismatch_far"] == 0 for v in rep["slabs"].values())


@needs_mjfit
def test_findings_on_bundled_configs():
    # islands drawn as circumscribed polygons are fatter than the real pins
    _, rep = certify_mjfit_config("db9")
    out = rep["slabs"]["plate_hole"]["deviation"]["sim_outside_reference"]
    assert abs(out - 0.0005 * (1 / math.cos(math.pi / 12) - 1)) <= rep["resolution"]
    # USB-C: clearance (50 um) below the default polygon error (100 um): the sim peg's
    # corners reach past the real receptacle wall
    _, rep = certify_mjfit_config("usb_c")
    r, c = 0.00128, 5e-5
    expected = r - (r - c) / math.cos(math.pi / 10)
    assert abs(rep["clearance"]["sim_peg_to_real_plate_min"] - expected) <= rep["resolution"]
    assert any("overlaps the real plate" in w for w in rep["warnings"])


@needs_mjfit
def test_threads_and_gears_are_not_claimed():
    import mjfit
    for name, kind in mjfit.configs().items():
        if kind in ("thread", "gear"):
            _, rep = certify_mjfit_config(name)
            assert rep["status"] == "refused" and rep["refusal"]["kind"] == "UnsupportedGenerator"
            break
    thread_xml = mjfit.thread(major_diameter=0.010, pitch=0.0015)
    rep = certify_model(thread_xml, {}, 1e-4, 1e-6)
    assert rep["status"] == "refused"


@needs_mjfit
def test_obj_mesh_dir_reads_the_same(tmp_path):
    kw = dict(CIRCLE, mesh_dir=str(tmp_path / "meshes"))
    xml, rep = certify_mjfit_peg(**kw)
    # the model refers to meshes relative to its own folder
    from brepkernel.mujoco_fit.mjcf import read_pieces
    pieces = read_pieces(xml, base_dir=tmp_path)
    assert len(pieces["plate"]) == 32 and len(pieces["peg"]) == 1


# --- mode C: CAD plates --------------------------------------------------------

def _cad():
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCone, BRepPrimAPI_MakeCylinder
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

    def cyl(x, y, z, r, h):
        return BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(x, y, z), gp_Dir(0, 0, 1)), r, h).Shape()

    def box(*c):
        return BRepPrimAPI_MakeBox(gp_Pnt(*c[:3]), gp_Pnt(*c[3:])).Shape()

    cone = BRepPrimAPI_MakeCone(gp_Ax2(gp_Pnt(0, 0, -1e-5), gp_Dir(0, 0, 1)), 0.008, 0.009, 0.01002).Shape()
    return {
        "through": BRepAlgoAPI_Cut(box(-.025, -.025, 0, .025, .025, .01), cyl(0, 0, -.001, .008, .012)).Shape(),
        "blind_pin": BRepAlgoAPI_Fuse(BRepAlgoAPI_Cut(box(-.02, -.02, 0, .02, .02, .01),
                                                      cyl(0, 0, .002, .007, .009)).Shape(),
                                      cyl(0, 0, .001, .001, .009)).Shape(),
        "chamfer": BRepAlgoAPI_Cut(box(-.025, -.025, 0, .025, .025, .01), cone).Shape(),
    }


@needs_mjfit
def test_cad_through_hole(tmp_path):
    from brepkernel.cli import write_shape
    path = str(tmp_path / "plate.brep")                                # BREP files are in meters
    write_shape(_cad()["through"], path)
    xml, rep = mjfit_peg_from_cad(path, clearance=1e-4, max_error=5e-5)
    assert rep["status"] == "certified", rep.get("refusal")
    assert rep["plate_thickness"] == pytest.approx(0.01) and rep["bottom_thickness"] == 0
    plate = rep["slabs"]["plate_hole"]["deviation"]
    assert plate["sim_outside_reference"] <= rep["resolution"]          # hole never smaller than the CAD hole
    assert plate["reference_outside_sim"] <= 5e-5 + rep["resolution"]  # within max_error
    assert rep["clearance"]["sim_peg_to_real_plate_min"] > 0           # the sim peg fits the CAD hole


@needs_mjfit
def test_cad_blind_hole_with_pin_and_step_units(tmp_path):
    from brepkernel.cli import frame_shape, write_shape
    shape = _cad()["blind_pin"]
    path = str(tmp_path / "pin.step")
    write_shape(frame_shape(shape, 1000.0, False), path)                # a STEP in millimeters
    xml, rep = mjfit_peg_from_cad(path, clearance=1e-4, max_error=5e-5)
    assert rep["status"] == "certified", rep.get("refusal")
    assert rep["bottom_thickness"] == pytest.approx(0.002) and rep["plate_thickness"] == pytest.approx(0.008)
    assert rep["slabs"]["plate_hole"]["boundary_loops"] == 3            # outer, hole, pin
    assert rep["slabs"]["plate_hole"]["deviation"]["sim_outside_reference"] <= rep["resolution"]
    assert rep["slabs"]["peg_insert"]["boundary_loops"] == 2            # peg ring around the pin pocket


@needs_mjfit
def test_cad_refusals():
    xml, rep = mjfit_peg_from_cad(_cad()["chamfer"], clearance=1e-4)
    assert xml is None and rep["refusal"]["kind"] == "NotAnExtrusion"
    xml, rep = mjfit_peg_from_cad(_cad()["through"], clearance=0.02)    # bigger than the hole
    assert rep["status"] == "refused"


@needs_mjfit
def test_cli_from_cad_and_refusal(tmp_path):
    import json
    import os
    import subprocess
    import sys
    from brepkernel.cli import write_shape
    src = os.path.join(os.path.dirname(__file__), "..", "src")
    env = dict(os.environ, PYTHONPATH=os.path.abspath(src))
    shapes = _cad()
    write_shape(shapes["through"], str(tmp_path / "ok.brep"))
    write_shape(shapes["chamfer"], str(tmp_path / "no.brep"))

    def run(*a):
        p = subprocess.run([sys.executable, "-m", "brepkernel", *a], capture_output=True, text=True,
                           env=env, cwd=tmp_path, timeout=900)
        return p.returncode, json.loads(p.stdout)

    code, rep = run("mjfit-from-cad", "ok.brep", "--clearance", "1e-4", "-o", "ok.xml")
    assert code == 0 and rep["status"] == "certified" and (tmp_path / "ok.xml").exists()
    import mujoco
    mujoco.MjModel.from_xml_path(str(tmp_path / "ok.xml"))
    code, rep = run("mjfit-from-cad", "no.brep", "--clearance", "1e-4", "-o", "no.xml")
    assert code == 2 and rep["certificate"]["refusal"]["kind"] == "NotAnExtrusion"
    assert not (tmp_path / "no.xml").exists()
