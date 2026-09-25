"""G7 periodic-seam stress tests (25 Sept 2026).

The old G7 plan (pre-split every periodic face with
ShapeUpgrade_ShapeDivideClosed) is dead. v0.9 solved periodic seams with a
smaller mechanism: a one-sided verified seam edge reuses the existing seam
boundary on the seam-side operand and splits normally on the opposite
operand; a seam on both operands stays unresolved. This file STRESS-TESTS
that mechanism; it does not change it.

Coverage:
  g1 equatorial seam-plane cut, seam on operand A, counters pinned
  g2 polar seam-plane cut, seam on operand A, counters pinned
  g3 polar seam-plane cut with reversed operands, seam on operand B
  g4 torus-vs-torus union, section curves crossing seams
  g5 torus-vs-box difference, seam-crossing section
  g6 revolved freeform NURBS profile vs box, seam-carrying faces
  g7 seam-routing evidence counters present in every public report

Contract per test: accepted cases must agree with the independent OCCT
Boolean oracle volume AND pass the independent generalized-winding-number
membership arbiter with zero kernel errors. A typed refusal is allowed only
where the test says so; any WRONG accept is a stop-and-report event.
"""
import math
import sys

sys.path.insert(0, "src")
sys.path.insert(0, "tests")

import _arbiter

from brepkernel import BRepAmbiguousResult, boolean_brep

from OCP.BRepAlgoAPI import (
    BRepAlgoAPI_Common,
    BRepAlgoAPI_Cut,
    BRepAlgoAPI_Fuse,
)
from OCP.BRepBuilderAPI import (
    BRepBuilderAPI_MakeEdge,
    BRepBuilderAPI_MakeFace,
    BRepBuilderAPI_MakeWire,
    BRepBuilderAPI_NurbsConvert,
)
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import (
    BRepPrimAPI_MakeBox,
    BRepPrimAPI_MakeRevol,
    BRepPrimAPI_MakeTorus,
)
from OCP.GeomAPI import GeomAPI_PointsToBSpline
from OCP.GProp import GProp_GProps
from OCP.collections import Array1_gp_Pnt as TColgp_Array1OfPnt
from OCP.gp import gp_Ax1, gp_Ax2, gp_Dir, gp_Pnt


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}", flush=True)
    return bool(cond)


def to_nurbs(shape):
    conv = BRepBuilderAPI_NurbsConvert(shape, True)
    assert conv.IsDone()
    out = conv.Shape()
    assert BRepCheck_Analyzer(out, True).IsValid()
    return out


def fuse_oracle(a, b):
    x = BRepAlgoAPI_Fuse(a, b)
    x.Build()
    assert x.IsDone()
    out = x.Shape()
    assert BRepCheck_Analyzer(out, True).IsValid()
    return out


def cut_oracle(a, b):
    x = BRepAlgoAPI_Cut(a, b)
    x.Build()
    assert x.IsDone()
    out = x.Shape()
    assert BRepCheck_Analyzer(out, True).IsValid()
    return out


def common_oracle(a, b):
    x = BRepAlgoAPI_Common(a, b)
    x.Build()
    assert x.IsDone()
    out = x.Shape()
    assert BRepCheck_Analyzer(out, True).IsValid()
    return out


def volume(shape):
    g = GProp_GProps()
    BRepGProp.VolumePropertiesGK_s(shape, g, 1e-10, True, True, False, False, False)
    return float(g.Mass())


def assert_oracle_close(name, out, oracle_shape, report, rel=5e-6, abs_tol=1e-9):
    rv = volume(out)
    ov = volume(oracle_shape)
    err = abs(rv - ov)
    lim = max(abs_tol, rel * abs(ov))
    return check(f"{name} volume vs OCCT oracle",
                 report.get("accepted") and err <= lim,
                 f"result={rv:.10f} oracle={ov:.10f} err={err:.3e} lim={lim:.3e}")


def torus():
    return to_nurbs(BRepPrimAPI_MakeTorus(3.0, 1.0).Shape())


def seam_report_checks(name, rep, reused_a, reused_b, shared=0,
                       seam_on_a=False, seam_on_b=False):
    sp = rep["stages"]["split"]
    asm = rep["stages"]["assembly"]
    ver = rep["stages"]["verification"]
    ok = check(
        f"{name} seam-routing evidence pinned",
        not sp["unresolved_contacts"]
        and sp["reused_seam_edges_A"] == reused_a
        and sp["reused_seam_edges_B"] == reused_b
        and sp["shared_seam_refusals"] == shared,
        f"split=A:{sp['reused_seam_edges_A']} B:{sp['reused_seam_edges_B']} "
        f"shared:{sp['shared_seam_refusals']} "
        f"unresolved={sp['unresolved_contacts']}")
    flags = [p["risk_flags"] for p in asm["section_payloads"]]
    ok &= check(
        f"{name} seam risk flags in section payloads",
        (any("seam_on_a" in f for f in flags) if seam_on_a else True)
        and (any("seam_on_b" in f for f in flags) if seam_on_b else True),
        f"flags={flags}")
    ok &= check(
        f"{name} result closed and auditable",
        asm["free_edges"] == 0
        and asm["multiple_edges"] == 0
        and asm["edge_lineage"]["unattributed_edges"] == 0
        and ver["brep_valid"]
        and ver["closed"]
        and ver["manifold_edges"]
        and ver["complete_edge_lineage"]
        and ver["volume_bounds_ok"],
        f"valid={ver.get('brep_valid')} closed={ver.get('closed')}")
    return ok


def g1_equatorial_seam_plane_cut():
    """Torus cut by its equatorial plane: one section loop IS the v-seam.

    The v-seam loop is reused as existing topology on the torus (operand A)
    and remains a split tool on the cutter (operand B).
    """
    tor = torus()
    cutter = BRepPrimAPI_MakeBox(
        gp_Pnt(-5.0, -5.0, -2.0),
        gp_Pnt(5.0, 5.0, 0.0)).Shape()
    try:
        out, rep = boolean_brep(tor, cutter, "difference")
    except BRepAmbiguousResult as exc:
        return check("g1 equatorial cut accepted",
                     False, f"refused: {exc.report.get('refusal')}")
    ok = seam_report_checks("g1", rep, 1, 0, seam_on_a=True)
    ok &= check(
        "g1 two loops from one face pair",
        rep["stages"]["intersection"]["verified_edges"] == 2
        and rep["stages"]["intersection"]["candidate_face_pairs"] == 1,
        f"ix={rep['stages']['intersection']}")
    ok &= assert_oracle_close("g1", out, cut_oracle(tor, cutter), rep)
    ok &= _arbiter.check_accepted(
        "g1", check, tor, cutter, out, "difference", n=200)[0]
    return ok, rep


def g2_polar_seam_plane_cut():
    """Torus cut by the plane y=0: one section loop IS the u-seam meridian.

    The meridian at major angle 0 lies in the cutter plane, so the section
    edge is the torus's existing periodic closing boundary.
    """
    tor = torus()
    cutter = BRepPrimAPI_MakeBox(
        gp_Pnt(-10.0, -10.0, -10.0),
        gp_Pnt(10.0, 0.0, 10.0)).Shape()
    try:
        out, rep = boolean_brep(tor, cutter, "difference")
    except BRepAmbiguousResult as exc:
        return check("g2 polar cut accepted",
                     False, f"refused: {exc.report.get('refusal')}")
    ok = seam_report_checks("g2", rep, 1, 0, seam_on_a=True)
    ok &= assert_oracle_close("g2", out, cut_oracle(tor, cutter), rep)
    ok &= _arbiter.check_accepted(
        "g2", check, tor, cutter, out, "difference", n=200)[0]
    return ok, rep


def g3_polar_seam_plane_cut_reversed():
    """Same polar cut with the torus as operand B: seam routes to B."""
    tor = torus()
    cutter = BRepPrimAPI_MakeBox(
        gp_Pnt(-10.0, -10.0, -10.0),
        gp_Pnt(10.0, 0.0, 10.0)).Shape()
    try:
        out, rep = boolean_brep(cutter, tor, "intersection")
    except BRepAmbiguousResult as exc:
        return check("g3 reversed polar cut accepted",
                     False, f"refused: {exc.report.get('refusal')}")
    ok = seam_report_checks("g3", rep, 0, 1, seam_on_b=True)
    ok &= assert_oracle_close("g3", out, common_oracle(cutter, tor), rep)
    ok &= _arbiter.check_accepted(
        "g3", check, cutter, tor, out, "intersection", n=200)[0]
    return ok, rep


def g4_torus_vs_torus_union():
    """Two NURBS tori, small one hooked through the big one's tube.

    Four transverse section curves cross both tori's seam parameter lines
    geometrically; none sits on a seam, so routing counters stay 0/0.
    """
    tor = torus()
    small = to_nurbs(BRepPrimAPI_MakeTorus(
        gp_Ax2(gp_Pnt(3.0, 0.0, 0.4), gp_Dir(0, 1, 0)), 1.4, 0.55).Shape())
    try:
        out, rep = boolean_brep(tor, small, "union")
    except BRepAmbiguousResult as exc:
        return check("g4 torus-torus union accepted",
                     False, f"refused: {exc.report.get('refusal')}")
    ok = seam_report_checks("g4", rep, 0, 0)
    ok &= assert_oracle_close("g4", out, fuse_oracle(tor, small), rep)
    ok &= _arbiter.check_accepted(
        "g4", check, tor, small, out, "union", n=200)[0]
    return ok, rep


def g5_torus_vs_box_difference():
    """Box corner driven through the torus tube wall: 8 section curves.

    The curves cross the torus seam parameter lines; none coincides with a
    seam, so the v0.9 routing stays inactive but must not misfire.
    """
    tor = torus()
    box = BRepPrimAPI_MakeBox(
        gp_Pnt(2.8, -0.5, -0.5),
        gp_Pnt(5.2, 0.5, 0.5)).Shape()
    try:
        out, rep = boolean_brep(tor, box, "difference")
    except BRepAmbiguousResult as exc:
        return check("g5 torus-box difference accepted",
                     False, f"refused: {exc.report.get('refusal')}")
    ok = seam_report_checks("g5", rep, 0, 0)
    ok &= assert_oracle_close("g5", out, cut_oracle(tor, box), rep)
    ok &= _arbiter.check_accepted(
        "g5", check, tor, box, out, "difference", n=200)[0]
    return ok, rep


def _egg_solid():
    n = 12
    pts = TColgp_Array1OfPnt(1, n + 1)
    for i in range(n + 1):
        th = 2 * math.pi * i / n
        r = 0.55 + 0.18 * math.sin(2 * th)
        pts.SetValue(i + 1, gp_Pnt(2.2 + r * math.cos(th), 0.0, r * math.sin(th)))
    curve = GeomAPI_PointsToBSpline(pts, 3).Curve()
    edge = BRepBuilderAPI_MakeEdge(curve).Edge()
    wire = BRepBuilderAPI_MakeWire(edge).Wire()
    face = BRepBuilderAPI_MakeFace(wire, True).Face()
    egg = BRepPrimAPI_MakeRevol(
        face, gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))).Shape()
    assert BRepCheck_Analyzer(egg, True).IsValid()
    return egg


def g6_revolved_nurbs_vs_box_union():
    """Freeform revolved B-spline profile (surface of revolution, one face).

    The revolve start/end meridian is a periodic seam; the box clip cuts two
    section loops across it.
    """
    egg = _egg_solid()
    box = BRepPrimAPI_MakeBox(
        gp_Pnt(1.0, -3.0, -3.0),
        gp_Pnt(3.4, 3.0, 3.0)).Shape()
    try:
        out, rep = boolean_brep(egg, box, "union")
    except BRepAmbiguousResult as exc:
        return check("g6 revolved NURBS union accepted",
                     False, f"refused: {exc.report.get('refusal')}")
    ok = seam_report_checks("g6", rep, 0, 0)
    ok &= assert_oracle_close("g6", out, fuse_oracle(egg, box), rep)
    ok &= _arbiter.check_accepted(
        "g6", check, egg, box, out, "union", n=200)[0]
    return ok, rep


def g7_seam_evidence_in_every_report(reports):
    """Every accepted public report carries the v0.9 seam-routing counters."""
    ok = True
    for name, rep in reports:
        sp = rep["stages"]["split"]
        ok &= check(
            f"g7 {name} report has seam-routing evidence",
            sp["reused_seam_edges_A"] >= 0
            and sp["reused_seam_edges_B"] >= 0
            and sp["shared_seam_refusals"] >= 0
            and isinstance(sp["unresolved_contacts"], list),
            f"A={sp['reused_seam_edges_A']} B={sp['reused_seam_edges_B']} "
            f"shared={sp['shared_seam_refusals']}")
    return ok


def main():
    ok = True
    reports = []
    for name, fn in (("g1", g1_equatorial_seam_plane_cut),
                     ("g2", g2_polar_seam_plane_cut),
                     ("g3", g3_polar_seam_plane_cut_reversed),
                     ("g4", g4_torus_vs_torus_union),
                     ("g5", g5_torus_vs_box_difference),
                     ("g6", g6_revolved_nurbs_vs_box_union)):
        res = fn()
        if isinstance(res, tuple):
            passed, rep = res
            reports.append((name, rep))
        else:
            passed = res
        ok &= passed
    ok &= g7_seam_evidence_in_every_report(reports)
    print("\nALL PASS" if ok else "\nSOME FAILURES", flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
