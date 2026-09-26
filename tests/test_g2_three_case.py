"""G2-planar: unit tests for the three-case coincidence classification.

Replaces the "exact coincidence via fractions.Fraction on float analytic
parameters" rule with:
  (a) EXACT: bit-identical or canonical analytic parameters;
  (b) TOLERANCE-CERTIFIED: deviation within the entities' ACTUAL
      tolerances (each face's own OCCT tolerance, plus contact_tol);
  (c) UNDECIDABLE: anything else -> typed refusal carrying the measured
      deviation and the tolerances consulted.

Scope is PLANAR only. Curved analytic (cylinder/sphere/cone/torus) and
NURBS coincidence machinery is deferred (see
src/brepkernel/coincidence_deferred.py): any curved pair that is not
provably distinct is case (c).

Written before the rework (I4); the marked tests fail until it lands.
"""
import sys

sys.path.insert(0, "src")

from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
from OCP.BRepTools import BRepTools
from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCP.Geom import Geom_BSplineSurface
try:  # OCP 8.x
    from OCP.collections import Array2_gp_Pnt, Array1_double, Array1_int
except ImportError:  # OCP 7.8
    from OCP.TColgp import TColgp_Array2OfPnt as Array2_gp_Pnt
    from OCP.TColStd import TColStd_Array1OfReal as Array1_double
    from OCP.TColStd import TColStd_Array1OfInteger as Array1_int
from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt
from OCP.TopoDS import TopoDS

from brepkernel import boolean_brep, BRepAmbiguousResult
from brepkernel.coincidence import classify_support_pair
try:
    from brepkernel.coincidence import CoincidenceUndecidable
except ImportError:  # pre-rework module: the new typed record is absent
    CoincidenceUndecidable = None
from brepkernel.step_ingest import FaceRecord, index_shape


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def box(x0, y0, z0, x1, y1, z1):
    return BRepPrimAPI_MakeBox(gp_Pnt(x0, y0, z0),
                              gp_Pnt(x1, y1, z1)).Shape()


def cyl(p, d, r, h):
    return BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(*p), gp_Dir(*d)), r, h).Shape()


def bspline_planar_face(x0, tol=1e-7):
    """Exactly-planar degree-1x1 BSpline face at x=x0 (non-analytic)."""
    P = Array2_gp_Pnt(1, 2, 1, 2)
    P.SetValue(1, 1, gp_Pnt(x0, -1.0, -1.0))
    P.SetValue(2, 1, gp_Pnt(x0, 1.0, -1.0))
    P.SetValue(1, 2, gp_Pnt(x0, -1.0, 1.0))
    P.SetValue(2, 2, gp_Pnt(x0, 1.0, 1.0))
    U = Array1_double(1, 2); U.SetValue(1, 0.0); U.SetValue(2, 1.0)
    V = Array1_double(1, 2); V.SetValue(1, 0.0); V.SetValue(2, 1.0)
    MU = Array1_int(1, 2); MU.SetValue(1, 2); MU.SetValue(2, 2)
    MV = Array1_int(1, 2); MV.SetValue(1, 2); MV.SetValue(2, 2)
    surf = Geom_BSplineSurface(P, U, V, MU, MV, 1, 1, False, False)
    return TopoDS.Face(BRepBuilderAPI_MakeFace(surf, tol).Face())


def face_record(face, fid=0):
    st = getattr(BRepAdaptor_Surface(face).GetType(), "name")
    return FaceRecord(
        fid, -1, -1, face,
        getattr(face.Orientation(), "name", str(face.Orientation())),
        st, tuple(float(x) for x in BRepTools.UVBounds_s(face)),
        (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), None)


def planar_face_of(shape, axis, sign):
    """First planar face of `shape` whose outward normal is sign*axis."""
    m = index_shape(shape)
    for fr in m.faces:
        if fr.surface_type != "GeomAbs_Plane":
            continue
        n = BRepAdaptor_Surface(
            TopoDS.Face(fr.face)).Plane().Axis().Direction()
        v = (n.X(), n.Y(), n.Z())
        if fr.orientation == "TopAbs_REVERSED":
            v = (-v[0], -v[1], -v[2])
        want = [0.0, 0.0, 0.0]
        want[axis] = sign
        if all(abs(v[i] - want[i]) < 1e-12 for i in range(3)):
            return fr
    raise AssertionError("planar face not found")


def refusal_kind(fn):
    try:
        fn()
    except BRepAmbiguousResult as e:
        return (e.report.get("refusal") or {}).get("kind")
    return None


def is_undecidable(r):
    return (CoincidenceUndecidable is not None
            and isinstance(r, CoincidenceUndecidable))


BT = 1e-7


def main():
    ok = True

    # (a) EXACT: shared-face boxes classify coincident via canonical
    # analytic parameters.
    fa = planar_face_of(box(0, 0, 0, 1, 1, 1), 0, +1.0)
    fb = planar_face_of(box(1, 0, 0, 2, 1, 1), 0, -1.0)
    r = classify_support_pair(fa, fb, BT)
    ok &= check("g2t exact shared face is case (a)",
                isinstance(r, tuple) and len(r) == 3 and r[0] == "coincident"
                and r[2] == "exact",
                f"-> {r if not is_undecidable(r) else r.reason}")

    # (a) is symmetric: reversed operand order still classifies.
    r = classify_support_pair(fb, fa, BT)
    ok &= check("g2t exact shared face symmetric",
                isinstance(r, tuple) and len(r) == 3 and r[0] == "coincident"
                and r[2] == "exact")

    # (c): the old Fraction rule ACCEPTED exactly-coincident curved
    # supports (cylinder/cylinder here). Curved coincidence is deferred:
    # the new rule must refuse instead of accepting.
    ma = index_shape(cyl((0, 0, 0), (0, 0, 1), 1.0, 2))
    mb = index_shape(cyl((0, 0, 0), (0, 0, 1), 1.0, 2))
    la = next(f for f in ma.faces if f.surface_type == "GeomAbs_Cylinder")
    lb = next(f for f in mb.faces if f.surface_type == "GeomAbs_Cylinder")
    r = classify_support_pair(la, lb, BT)
    ok &= check("g2t coincident cylinders are case (c) [deferred]",
                is_undecidable(r)
                and r.reason == "curved_deferred",
                f"-> {r if not is_undecidable(r) else r.reason}")

    # (c) end-to-end: union of coaxial equal-radius cylinders (stacked;
    # NOT identical inputs, so the Tier A identical fast path does not
    # apply) refuses with kind CoincidenceUndecidable. The old code
    # accepted this via the curved exact-Fraction path.
    kind = refusal_kind(lambda: boolean_brep(
        cyl((0, 0, 0), (0, 0, 1), 1.0, 2),
        cyl((0, 0, 1), (0, 0, 1), 1.0, 1), "union"))
    ok &= check("g2t stacked-cylinder union refuses CoincidenceUndecidable",
                kind == "CoincidenceUndecidable", f"kind={kind}")

    # (c): near-coincident planar pair inside the band carries the
    # measured deviation and every tolerance consulted.
    fa = planar_face_of(box(0, 0, 0, 1, 1, 1), 0, +1.0)
    fb = planar_face_of(box(1 + 0.5 * BT, 0, 0, 2 + 0.5 * BT, 1, 1),
                        0, -1.0)
    r = classify_support_pair(fa, fb, BT)
    is_c = is_undecidable(r)
    ok &= check("g2t 0.5x-tol planar offset is case (c)", is_c,
                f"-> {type(r).__name__}")
    if is_c:
        ok &= check("g2t case (c) carries measured deviation",
                    abs(r.deviation - 0.5 * BT) <= 1e-12,
                    f"deviation={r.deviation:.3g}")
        names = {t["name"] for t in r.tolerances}
        ok &= check("g2t case (c) records per-entity tolerances",
                    {"face_tolerance_a", "face_tolerance_b",
                     "contact_tol", "coincidence_band"} <= names,
                    f"tolerances={sorted(names)}")
        vals = {t["name"]: t["value"] for t in r.tolerances}
        ok &= check("g2t case (c) uses actual OCCT tolerances",
                    abs(vals["face_tolerance_a"] - 1e-7) < 1e-18
                    and abs(vals["face_tolerance_b"] - 1e-7) < 1e-18,
                    f"face tols={vals['face_tolerance_a']:.3g},"
                    f"{vals['face_tolerance_b']:.3g}")

    # (b) TOLERANCE-CERTIFIED: recognized-planar BSpline face within the
    # entities' actual tolerances of an analytic plane. The fitted plane
    # is not bit-identical, so case (a) cannot apply; the deviation
    # (5e-8) is within the faces' own 1e-7 tolerances.
    fa = planar_face_of(box(0, 0, 0, 1, 1, 1), 0, +1.0)
    fb = face_record(bspline_planar_face(1.0 + 0.5 * BT), fid=1)
    r = classify_support_pair(fa, fb, BT)
    ok &= check("g2t in-tolerance recognized planar is case (b)",
                isinstance(r, tuple) and len(r) == 3 and r[0] == "coincident"
                and r[2] == "tolerance_certified",
                f"-> {r if not is_undecidable(r) else r.reason}")

    # (c): recognized-planar pair beyond the entities' tolerances but
    # inside the band refuses (the old sampling path accepted this).
    fa = planar_face_of(box(0, 0, 0, 1, 1, 1), 0, +1.0)
    fb = face_record(bspline_planar_face(1.0 + 2.0 * BT), fid=1)
    r = classify_support_pair(fa, fb, BT)
    ok &= check("g2t 2x-tol recognized planar is case (c)",
                is_undecidable(r),
                f"-> {type(r).__name__}")

    # Distinct planar pair stays distinct (no behavior change).
    fa = planar_face_of(box(0, 0, 0, 1, 1, 1), 0, +1.0)
    fb = planar_face_of(box(1 + 10.0 * BT, 0, 0, 2 + 10.0 * BT, 1, 1),
                        0, -1.0)
    r = classify_support_pair(fa, fb, BT)
    ok &= check("g2t 10x-tol planar offset is distinct", r == "distinct",
                f"-> {type(r).__name__ if is_undecidable(r) else r}")

    # Curved pairs that are provably different stay "distinct" (the
    # section still handles them; only coincidence acceptance is
    # deferred, not distinctness).
    ma = index_shape(cyl((0, 0, 0), (0, 0, 1), 1.0, 2))
    mb = index_shape(cyl((0, 0, -1), (0, 0, 1), 0.5, 4))
    la = next(f for f in ma.faces if f.surface_type == "GeomAbs_Cylinder")
    lb = next(f for f in mb.faces if f.surface_type == "GeomAbs_Cylinder")
    r = classify_support_pair(la, lb, BT)
    ok &= check("g2t different-radii cylinders stay distinct",
                r == "distinct", f"-> {type(r).__name__ if is_undecidable(r) else r}")

    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
