"""brepkernel describe / query / export: the exact answers DVE uses for Measure, exact snapping
and File > Export Exact CAD (STEP).

describe: faces in mesh order with their surface data (plane outward normals, bores flagged as
holes, cones, spheres, tori, B-splines), true corners only (no seam or pole vertices), feature
edges with arc-length midpoints, circles merged across arcs with ``full`` set only when they
close. query: exact closest points, the face in mesh order, outward normals, and the face the
point is in front of at edges and corners. export: placed parts, named, in mm or m and Z or Y
up; part and read-back checks; a failure writes nothing and keeps an existing file.
"""
import json
import math
import os

import pytest

import brepkernel._occt_compat  # noqa: F401
from brepkernel.cad_inspect import describe_shape, query_points
from brepkernel.cli import main, place_shape, write_shape
from brepkernel.mesh_export import _faces
from brepkernel.step_export import _ascii_escape, _read_back
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCone,
                             BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeSphere,
                             BRepPrimAPI_MakeTorus)
from OCP.TopAbs import TopAbs_OUT
from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

R, H = 0.006, 0.02                        # the 12 mm pin DVE's real-brepkernel tests use


def pin():
    return BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, -H / 2, 0), gp_Dir(0, 1, 0)), R, H).Shape()


def plate():
    """50 x 10 x 50 mm plate with a 12 mm bore along Y."""
    box = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), gp_Pnt(0.05, 0.01, 0.05)).Shape()
    bore = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0.025, -0.01, 0.025), gp_Dir(0, 1, 0)), R, 0.03).Shape()
    return BRepAlgoAPI_Cut(box, bore).Shape()


def run(tmp_path, *argv):
    report = tmp_path / "report.json"
    code = main([*map(str, argv), "--report", str(report)])
    return code, json.loads(report.read_text())


# --------------------------------------------------------------------------- describe

def test_pin_circles_faces_and_no_false_corners():
    d = describe_shape(pin())
    assert d["solids"] == 1
    assert d["volume"] == pytest.approx(math.pi * R * R * H, rel=1e-12)
    assert d["area"] == pytest.approx(2 * math.pi * R * (R + H), rel=1e-9)
    full = [c for c in d["circles"] if c["full"] and abs(c["radius"] - R) < 1e-12]
    assert sorted(c["centre"][1] for c in full) == pytest.approx([-H / 2, H / 2])
    side = [f for f in d["faces"] if f["type"] == "cylinder"]
    assert len(side) == 1 and side[0]["radius"] == pytest.approx(R) and side[0]["hole"] is False
    assert side[0]["axis"] == pytest.approx([0, 1, 0])
    assert d["vertices"] == []                            # the seam vertex is not a corner
    assert d["skipped"]["seam_edges"] == 1
    assert all(e["type"] == "circle" for e in d["edges"])  # the seam line is not a feature edge


def test_faces_in_mesh_order_with_outward_plane_normals():
    shape = plate()
    d = describe_shape(shape)
    faces = _faces(shape)
    assert [f["index"] for f in d["faces"]] == list(range(len(faces)))
    for rec, face in zip(d["faces"], faces):
        assert rec["area"] > 0
        if rec["type"] == "plane":
            o, n = rec["origin"], rec["normal"]
            assert math.isclose(math.hypot(*n), 1.0, rel_tol=1e-12)
            probe = gp_Pnt(*(o[i] + 1e-4 * n[i] for i in range(3)))
            cls = BRepClass3d_SolidClassifier(shape, probe, 1e-9)
            assert cls.State() == TopAbs_OUT                # just outside the material
    bore = [f for f in d["faces"] if f["type"] == "cylinder"]
    assert len(bore) == 1 and bore[0]["hole"] is True
    assert len(d["vertices"]) == 8                        # the box corners
    assert {round(e["length"], 9) for e in d["edges"] if e["type"] == "line"} == {0.05, 0.01}
    for e in d["edges"]:
        if e["type"] == "line":
            assert e["mid"] == pytest.approx([(a + b) / 2 for a, b in zip(e["start"], e["end"])])
    assert sum(c["full"] for c in d["circles"]) == 2


def test_partial_arcs_are_not_full_circles():
    box = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), gp_Pnt(1, 1, 1)).Shape()
    corner = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, -1, 0), gp_Dir(0, 1, 0)), 0.3, 3).Shape()
    d = describe_shape(BRepAlgoAPI_Cut(box, corner).Shape())
    assert d["circles"] and not any(c["full"] for c in d["circles"])
    arcs = [e for e in d["edges"] if e["type"] == "circle"]
    for a in arcs:
        assert a["length"] == pytest.approx(math.pi / 2 * 0.3, rel=1e-9)
        # arc-length midpoint: on the circle, halfway round the quarter
        assert math.hypot(a["mid"][0], a["mid"][2]) == pytest.approx(0.3)
        assert a["mid"][0] == pytest.approx(a["mid"][2])


def test_sphere_cone_torus_and_cavity():
    s = describe_shape(BRepPrimAPI_MakeSphere(gp_Pnt(1, 2, 3), 0.01).Shape())
    assert s["faces"][0]["type"] == "sphere" and s["faces"][0]["centre"] == pytest.approx([1, 2, 3])
    assert s["faces"][0]["hole"] is False and s["vertices"] == []
    assert s["skipped"]["degenerate_edges"] == 2
    c = describe_shape(BRepPrimAPI_MakeCone(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 0.01, 0.0, 0.02).Shape())
    cone = [f for f in c["faces"] if f["type"] == "cone"][0]
    assert cone["half_angle"] == pytest.approx(math.atan(0.5))
    assert cone["apex"] == pytest.approx([0, 0, 0.02])
    t = describe_shape(BRepPrimAPI_MakeTorus(2.0, 0.5).Shape())
    assert t["faces"][0]["major_radius"] == 2.0 and t["faces"][0]["minor_radius"] == 0.5
    cavity = BRepAlgoAPI_Cut(BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), gp_Pnt(4, 4, 4)).Shape(),
                             BRepPrimAPI_MakeSphere(gp_Pnt(2, 2, 2), 1.0).Shape()).Shape()
    assert [f["hole"] for f in describe_shape(cavity)["faces"] if f["type"] == "sphere"] == [True]


def test_stretched_part_reports_bspline_faces():
    stretched, kind = place_shape(pin(), [2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0])
    assert kind == "general_affine"
    d = describe_shape(stretched)
    spl = [f for f in d["faces"] if f["type"] == "bspline"]
    assert spl and all(f["u_degree"] >= 1 and f["poles"][0] >= 2 for f in spl)
    assert d["volume"] == pytest.approx(2 * math.pi * R * R * H, rel=1e-6)


# --------------------------------------------------------------------------- query

def test_query_exact_closest_points_and_faces():
    shape = pin()
    d = describe_shape(shape)
    side = next(f["index"] for f in d["faces"] if f["type"] == "cylinder")
    hits = query_points(shape, [(0.01, 0, 0), (0, 0.02, 0), (0, -0.03, 0)])
    assert [h["found"] for h in hits] == [True] * 3
    assert hits[0]["closest"][0] == pytest.approx(R, abs=1e-15)
    assert hits[0]["distance"] == pytest.approx(0.004, abs=1e-15)
    assert hits[0]["face"] == side and hits[0]["normal"] == pytest.approx([1, 0, 0])
    assert hits[1]["closest"] == pytest.approx([0, H / 2, 0]) and hits[1]["normal"] == pytest.approx([0, 1, 0])
    assert hits[2]["closest"] == pytest.approx([0, -H / 2, 0]) and hits[2]["normal"] == pytest.approx([0, -1, 0])
    assert d["faces"][hits[2]["face"]]["type"] == "plane"


def test_query_at_a_corner_picks_the_face_the_point_is_in_front_of():
    shape = plate()
    d = describe_shape(shape)
    above, beside = query_points(shape, [(0.0, 0.03, 0.0), (-0.02, 0.005, -1e-9)])
    assert above["closest"] == pytest.approx([0, 0.01, 0]) and above["normal"] == pytest.approx([0, 1, 0])
    assert d["faces"][above["face"]]["normal"] == pytest.approx([0, 1, 0])
    assert beside["normal"] == pytest.approx([-1, 0, 0])


def test_query_inside_the_bore_lands_on_the_bore():
    shape = plate()
    d = describe_shape(shape)
    (hit,) = query_points(shape, [(0.025, 0.005, 0.0251)])
    assert d["faces"][hit["face"]]["hole"] is True
    assert hit["distance"] == pytest.approx(R - 0.0001, rel=1e-9)


# --------------------------------------------------------------------------- command line

def test_cli_describe_and_query(tmp_path):
    src = tmp_path / "pin.brep"
    write_shape(pin(), str(src))
    code, rep = run(tmp_path, "describe", src)
    assert code == 0 and rep["schema"] == "brepkernel.cli/1" and rep["status"] == "ok"
    assert len(rep["description"]["faces"]) == 3
    code, rep = run(tmp_path, "query", src, "--point", 0.01, 0, 0, "--point", 0, -0.05, 0)
    assert code == 0 and len(rep["results"]) == 2 and rep["results"][1]["point"] == [0, -0.05, 0]
    code, rep = run(tmp_path, "query", src)
    assert code == 1 and "--point" in rep["error"]
    code, rep = run(tmp_path, "describe", tmp_path / "missing.brep")
    assert code == 1 and "not found" in rep["error"]


def _parts(tmp_path, parts):
    path = tmp_path / "parts.json"
    path.write_text(json.dumps({"parts": parts}, ensure_ascii=False))
    return path


def test_cli_export_named_checked_and_converted(tmp_path):
    src = tmp_path / "pin.brep"
    write_shape(pin(), str(src))
    parts = _parts(tmp_path, [
        {"path": str(src), "name": "Pin", "transform": [1, 0, 0, 0.05, 0, 1, 0, 0, 0, 0, 1, 0]},
        {"path": str(src), "name": "Wide pin", "transform": [2, 0, 0, -0.05, 0, 1, 0, 0, 0, 0, 1, 0]},
        {"path": str(src), "name": "Same place", "transform": [1, 0, 0, 0.05, 0, 1, 0, 0, 0, 0, 1, 0]}])
    out = tmp_path / "exports" / "pins.step"
    code, rep = run(tmp_path, "export", "--parts", parts, "-o", out, "--units", "mm", "--up", "z")
    assert code == 0, rep.get("error")
    e = rep["export"]
    pin_mm3 = math.pi * 6.0 * 6.0 * 20.0
    assert [p["placement"] for p in e["parts"]] == ["rigid", "general_affine", "rigid"]
    assert [p["exact"] for p in e["parts"]] == [True, False, True]
    assert e["parts"][0]["volume_mm3"] == pytest.approx(pin_mm3, rel=1e-9)
    assert e["parts"][1]["volume_mm3"] == pytest.approx(2 * pin_mm3, rel=1e-6)
    assert e["readback"]["solids"] == 3 and e["readback"]["relative_error"] < 1e-9
    text = out.read_text()
    assert "ISO-10303-21" in text and "'Pin'" in text and "'Wide pin'" in text and ".MILLI." in text
    assert not [f for f in os.listdir(out.parent) if f != "pins.step"]     # no temporary left


@pytest.mark.parametrize("units,up,axis", [("mm", "z", 2), ("m", "y", 1), ("m", "z", 2)])
def test_cli_export_units_and_up_axis(tmp_path, units, up, axis):
    box = tmp_path / "box.brep"
    write_shape(BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), gp_Pnt(0.01, 0.03, 0.02)).Shape(), str(box))
    parts = _parts(tmp_path, [{"path": str(box), "name": "Box", "transform": [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0]}])
    out = tmp_path / "box.step"
    code, rep = run(tmp_path, "export", "--parts", parts, "-o", out, "--units", units, "--up", up)
    assert code == 0, rep.get("error")
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib
    from OCP.STEPControl import STEPControl_Reader
    r = STEPControl_Reader()
    r.ReadFile(str(out))
    r.TransferRoots()
    b = Bnd_Box()
    BRepBndLib.Add_s(r.OneShape(), b)
    lo, hi = b.CornerMin(), b.CornerMax()
    extent = [hi.X() - lo.X(), hi.Y() - lo.Y(), hi.Z() - lo.Z()]
    assert extent[axis] == pytest.approx(30.0, abs=1e-5)          # the 30 mm height is "up"
    assert ("SI_UNIT(.MILLI.,.METRE.)" in out.read_text()) == (units == "mm")


def test_cli_export_failure_keeps_the_old_file(tmp_path):
    src = tmp_path / "pin.brep"
    write_shape(pin(), str(src))
    out = tmp_path / "pins.step"
    out.write_text("OLD")
    broken = tmp_path / "broken.brep"
    broken.write_text("not a brep")
    for parts, words in [
        ([{"path": str(broken), "name": "Broken"}], "Broken"),
        ([{"path": str(src), "name": "Flat", "transform": [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]}], "singular"),
        ([{"path": str(src), "name": "Bad", "transform": [1, 0, 0]}], "12 finite"),
        ([], "no parts"),
    ]:
        code, rep = run(tmp_path, "export", "--parts", _parts(tmp_path, parts), "-o", out)
        assert code == 1 and words in rep["error"], rep["error"]
        assert out.read_text() == "OLD"
    code, rep = run(tmp_path, "export", "--parts", _parts(tmp_path, [{"path": str(src), "name": "P"}]),
                    "-o", tmp_path / "pins.stl")
    assert code == 1 and ".step" in rep["error"]
    assert sorted(os.listdir(tmp_path)) == ["broken.brep", "parts.json", "pin.brep", "pins.step", "report.json"]


def test_cli_export_names_round_trip(tmp_path):
    src = tmp_path / "pin.brep"
    write_shape(pin(), str(src))
    names = ["Bob's pin", "Zahnrad Ø12 – 歯車", "back\\slash", "rocket \U0001f680"]
    parts = _parts(tmp_path, [{"path": str(src), "name": n, "transform": [1, 0, 0, 0.1 * i, 0, 1, 0, 0, 0, 0, 1, 0]}
                              for i, n in enumerate(names)])
    out = tmp_path / "names.step"
    code, rep = run(tmp_path, "export", "--parts", parts, "-o", out)
    assert code == 0, rep.get("error")
    out.read_bytes().decode("ascii")                       # escaped, plain ASCII
    assert sorted(_read_back(str(out))[2]) == sorted(names)
    assert _read_back(str(out))[3] == ["mm"]


def test_ascii_escape():
    assert _ascii_escape("aØb") == "a\\X2\\00D8\\X0\\b"
    assert _ascii_escape("\U0001f680") == "\\X4\\0001F680\\X0\\"
    assert _ascii_escape("plain") == "plain"


def test_export_restores_occt_settings(tmp_path):
    from OCP.Interface import Interface_Static
    from OCP.STEPCAFControl import STEPCAFControl_Writer
    STEPCAFControl_Writer()                                # registers the settings
    before = {k: Interface_Static.CVal_s(k) for k in ("write.step.unit", "write.step.schema", "xstep.cascade.unit")}
    src = tmp_path / "pin.brep"
    write_shape(pin(), str(src))
    parts = _parts(tmp_path, [{"path": str(src), "name": "P"}])
    assert run(tmp_path, "export", "--parts", parts, "-o", tmp_path / "p.step", "--units", "m")[0] == 0
    assert {k: Interface_Static.CVal_s(k) for k in before} == before
