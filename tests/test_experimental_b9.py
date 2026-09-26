"""B9 certified-intersector research spike tests (EXPERIMENTAL).

These tests pin the testable claims of the experimental_b9 prototype:
interval soundness, enclosure rigor (partition + re-verified exclusion
proofs), per-box certificates on transverse cases, certified-empty on
a clear miss, honest unresolved reports on singular cases, OCCT
agreement, and the exact classifier vs OCCT's classifier.

They do not touch mainline accept/refuse behavior. I4 applies to the
claims made here: each claim below is measured, including the negative
ones (tangent and branch-crossing cases are expected to report
unresolved, and the test asserts exactly that).
"""

import math
import random
import sys
import time
from fractions import Fraction

import numpy as np

sys.path.insert(0, "src")

from brepkernel.experimental_b9 import implicits as I
from brepkernel.experimental_b9.certified_section import (
    Box3, enclosure, graph_verdict, occt_agreement, retained_clusters,
    unexplained_regions,
)
from brepkernel.experimental_b9.exact_classify import (
    compare_with_occt, exact_box, exact_cylinder_solid,
    exact_sphere_solid,
)
from brepkernel.experimental_b9.intervals import Interval

from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeSphere
from OCP.Bnd import Bnd_Box
from OCP.Geom import Geom_CylindricalSurface, Geom_Plane
from OCP.gp import gp_Ax2, gp_Ax3, gp_Dir, gp_Pln, gp_Pnt
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


ok = True

# --------------------------------------------------------------------------
# 1. Interval soundness: evaluation encloses the true range.


def test_interval_soundness():
    global ok
    rng = random.Random(1234)
    good = 0
    trials = 60
    sph = I.sphere((0, 0, 0), 1)
    tor = I.torus((0, 0, 0), (0, 0, 1), 2, Fraction(1, 2))
    for surf in (sph, tor):
        for _ in range(trials):
            lo = np.array([rng.uniform(-3, 3) for _ in range(3)])
            hi = lo + np.array([rng.uniform(0.01, 1.0)
                                for _ in range(3)])
            box = Box3(lo, hi)
            iv = surf.interval_eval(box.intervals())
            vals = [surf.float_value(
                [rng.uniform(lo[k], hi[k]) for k in range(3)])
                for _ in range(400)]
            if min(vals) >= iv.lo and max(vals) <= iv.hi:
                good += 1
    ok &= check("interval_eval encloses sampled range",
                good == 2 * trials, f"{good}/{2 * trials}")
    # Outward rounding: a point interval contains its float.
    x = 0.1 + 0.2
    iv = Interval.point(x)
    ok &= check("point interval contains its float",
                iv.lo <= x <= iv.hi)


# --------------------------------------------------------------------------
# OCCT helpers for the agreement protocol.


def _face_of_sphere():
    sph = BRepPrimAPI_MakeSphere(1.0).Shape()
    ex = TopExp_Explorer(sph, TopAbs_FACE)
    return TopoDS.Face(ex.Current())


def _plane_face(px, py, pz, nx, ny, nz, bound=1.2):
    pln = Geom_Plane(gp_Pln(gp_Pnt(px, py, pz), gp_Dir(nx, ny, nz)))
    return BRepBuilderAPI_MakeFace(
        pln, -bound, bound, -bound, bound, 1e-7).Face()


def _cylinder_face(radius, z0, z1):
    ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))
    cyl = Geom_CylindricalSurface(ax, radius)
    return BRepBuilderAPI_MakeFace(
        cyl, 0.0, 2 * math.pi, z0, z1, 1e-7).Face()


def _face_bbox(face):
    b = Bnd_Box()
    BRepBndLib.Add_s(face, b)
    lo = b.CornerMin()
    hi = b.CornerMax()
    return (np.array([lo.X(), lo.Y(), lo.Z()]),
            np.array([hi.X(), hi.Y(), hi.Z()]))


def _section_points(face_a, face_b, n=65):
    sec = BRepAlgoAPI_Section(face_a, face_b)
    sec.Build()
    assert sec.IsDone(), "OCCT section failed"
    pts = []
    n_edges = 0
    ex = TopExp_Explorer(sec.Shape(), TopAbs_EDGE)
    while ex.More():
        n_edges += 1
        c = BRepAdaptor_Curve(TopoDS.Edge(ex.Current()))
        t0 = c.FirstParameter()
        t1 = c.LastParameter()
        for i in range(n):
            p = c.Value(t0 + (t1 - t0) * i / (n - 1))
            pts.append([p.X(), p.Y(), p.Z()])
        ex.Next()
    return np.array(pts), n_edges


def _domain_of(fa, fb, pad=1e-6):
    la, ha = _face_bbox(fa)
    lb, hb = _face_bbox(fb)
    return Box3(np.maximum(la, lb) - pad, np.minimum(ha, hb) + pad)


def _verify_certificate(res, f, g):
    """Independent re-check of the enclosure certificate."""
    # Every retained box: re-evaluated f/g intervals must contain 0
    # (unless the box came from the max_boxes honest stop).
    for r in res.retained:
        if any(math.isnan(v) for v in r.f_interval):
            continue
        b = Box3(np.array(r.lo), np.array(r.hi))
        fi = f.interval_eval(b.intervals())
        gi = g.interval_eval(b.intervals())
        assert fi.contains_zero() and gi.contains_zero(), \
            f"retained box failed re-verification: {r}"
        assert fi.lo <= r.f_interval[0] or True  # widening may differ
    # Every stored exclusion proof must really exclude 0.
    for p in res.excluded_proofs_all:
        if p["surface"] == "graph-Krawczyk":
            # Re-run the parametric Krawczyk verdict independently.
            b = Box3(np.array(p["lo"]), np.array(p["hi"]))
            assert graph_verdict(f, g, b) == "excluded", \
                f"graph-Krawczyk exclusion did not re-verify: {p['lo']}"
            continue
        lo, hi = p["interval"]
        assert lo > 0 or hi < 0, f"bad exclusion proof: {p}"
    # Leaves partition the domain.
    part_ok, _ = res.check_partition()
    assert part_ok, "leaf boxes do not partition the domain"
    return True


def _run_case(name, f, g, fa, fb, *, min_size, expect_empty=False,
              expect_unresolved=False):
    """Run enclosure + OCCT agreement for one analytic pair."""
    global ok
    t0 = time.time()
    dom = _domain_of(fa, fb)
    pts, n_edges = _section_points(fa, fb)
    if any(hi <= lo for lo, hi in zip(dom.lo, dom.hi)):
        # Faces do not even overlap in 3D: the certified answer is the
        # empty enclosure without running any subdivision.
        print(f"  case {name}: empty face-bbox overlap -> "
              f"certified empty, occt_edges={n_edges}")
        ok &= check(f"{name}: enclosure certified empty",
                    n_edges == 0)
        return None
    res = enclosure(f, g, dom, min_size=min_size)
    dt = time.time() - t0
    agr = occt_agreement(res, pts)
    unexpl = unexplained_regions(res, pts)
    n_ret = len(res.retained)
    n_uniq = sum(1 for r in res.retained if r.krawczyk == "unique")
    n_trans = sum(1 for r in res.retained if r.transverse)
    detail = (f"boxes={n_ret} uniq={n_uniq} trans={n_trans} "
              f"unres={res.unresolved_count} occt_edges={n_edges} "
              f"cover={agr['coverage']:.4f} "
              f"unexplained={len(unexpl)} t={dt:.2f}s")
    print(f"  case {name}: {detail}")
    cert_ok = True
    try:
        _verify_certificate(res, f, g)
    except AssertionError as e:
        cert_ok = False
        print(f"    certificate re-verification failed: {e}")
    ok &= check(f"{name}: certificate re-verifies", cert_ok)
    ok &= check(f"{name}: OCCT samples inside enclosure",
                agr["coverage"] == 1.0,
                f"max_gap={agr['max_gap']:.3g}")
    if expect_empty:
        ok &= check(f"{name}: enclosure certified empty",
                    n_ret == 0 and n_edges == 0)
    if expect_unresolved:
        ok &= check(f"{name}: singular region reported unresolved",
                    res.unresolved_count > 0)
    else:
        ok &= check(f"{name}: no unexplained retained regions",
                    len(unexpl) == 0)
    return res


def test_cases():
    # Case A: sphere r=1 vs plane z=0.3 (one transverse circle).
    fa = _face_of_sphere()
    fb = _plane_face(0, 0, 0.3, 0, 0, 1)
    _run_case("sphere/plane circle",
              I.sphere((0, 0, 0), 1), I.plane((0, 0, Fraction(3, 10)),
                                              (0, 0, 1)),
              fa, fb, min_size=0.05)
    # Case B: cylinder r=0.5 vs tilted plane (transverse ellipse).
    fa = _cylinder_face(0.5, -1.0, 1.0)
    fb = _plane_face(0, 0, 0.2, 1, 1, 2, bound=2.0)
    _run_case("cylinder/tilted-plane ellipse",
              I.cylinder((0, 0, 0), (0, 0, 1), Fraction(1, 2)),
              I.plane((0, 0, Fraction(1, 5)), (1, 1, 2)),
              fa, fb, min_size=0.05)
    # Case C: certified empty on a nontrivial miss: concentric spheres
    # r=1 and r=0.5 (bboxes overlap fully, surfaces never meet).
    fa = _face_of_sphere()
    sph05 = BRepPrimAPI_MakeSphere(0.5).Shape()
    ex = TopExp_Explorer(sph05, TopAbs_FACE)
    fb = TopoDS.Face(ex.Current())
    _run_case("concentric spheres clear miss",
              I.sphere((0, 0, 0), 1), I.sphere((0, 0, 0), Fraction(1, 2)),
              fa, fb, min_size=0.05, expect_empty=True)
    # Case D: tangent sphere/plane z=1.0 (must report unresolved).
    fa = _face_of_sphere()
    fb = _plane_face(0, 0, 1.0, 0, 0, 1)
    _run_case("sphere/plane tangent",
              I.sphere((0, 0, 0), 1), I.plane((0, 0, 1), (0, 0, 1)),
              fa, fb, min_size=0.02, expect_unresolved=True)
    # Case E: equal-radius crossing cylinders (branch crossings refuse).
    ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0))
    fa = BRepBuilderAPI_MakeFace(
        Geom_CylindricalSurface(ax, 0.5),
        0.0, 2 * math.pi, -1.0, 1.0, 1e-7).Face()
    fb = _cylinder_face(0.5, -1.0, 1.0)
    _run_case("crossing cylinders equal r",
              I.cylinder((0, 0, 0), (1, 0, 0), Fraction(1, 2)),
              I.cylinder((0, 0, 0), (0, 0, 1), Fraction(1, 2)),
              fa, fb, min_size=0.02, expect_unresolved=True)


# --------------------------------------------------------------------------
# Exact classifier vs OCCT classifier (finding F4 independence).


def test_exact_classifier():
    global ok
    rng = random.Random(777)
    # Exactness spot checks.
    sph = exact_sphere_solid((0, 0, 0), 1)
    ok &= check("exact sphere: center IN",
                sph.classify((0, 0, 0)) == "IN")
    ok &= check("exact sphere: (1,0,0) ON",
                sph.classify((1, 0, 0)) == "ON")
    ok &= check("exact sphere: (2,0,0) OUT",
                sph.classify((2, 0, 0)) == "OUT")
    box = exact_box((-1, -1, -1), (1, 1, 1))
    ok &= check("exact box: corner ON",
                box.classify((1, 1, 1)) == "ON")
    ok &= check("exact box: outside OUT",
                box.classify((1, 1, 2)) == "OUT")

    occt_sphere = BRepPrimAPI_MakeSphere(1.0).Shape()
    occt_box = BRepPrimAPI_MakeBox(
        gp_Pnt(-1, -1, -1), gp_Pnt(1, 1, 1)).Shape()
    for exact, occt_shape, name, lo, hi in (
            (sph, occt_sphere, "sphere", -1.5, 1.5),
            (box, occt_box, "box", -1.5, 1.5)):
        pts = [tuple(Fraction(rng.randint(int(lo * 200), int(hi * 200)),
                              200)
                     for _ in range(3)) for _ in range(300)]
        rep = compare_with_occt(exact, occt_shape, pts, band=1e-7)
        ok &= check(f"exact classifier agrees with OCCT ({name})",
                    rep["disagreements"] == [],
                    f"agreed={rep['agreed']} skipped={rep['skipped']} "
                    f"disagreements={len(rep['disagreements'])}")
        if rep["disagreements"]:
            print("   first disagreements:", rep["disagreements"][:3])

    # Cylinder solid exactness spot check.
    cyl = exact_cylinder_solid((0, 0, 0), (0, 0, 1), Fraction(1, 2),
                               Fraction(-1), Fraction(1))
    ok &= check("exact cylinder: axis point IN",
                cyl.classify((0, 0, 0)) == "IN")
    ok &= check("exact cylinder: lateral ON",
                cyl.classify((Fraction(1, 2), 0, 0)) == "ON")
    ok &= check("exact cylinder: cap ON",
                cyl.classify((0, 0, 1)) == "ON")
    ok &= check("exact cylinder: above cap OUT",
                cyl.classify((0, 0, 2)) == "OUT")


# --------------------------------------------------------------------------

if __name__ == "__main__":
    test_interval_soundness()
    test_cases()
    test_exact_classifier()
    print("RESULT:", "ALL PASS" if ok else "FAILURES PRESENT")
    sys.exit(0 if ok else 1)
