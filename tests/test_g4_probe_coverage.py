"""G4 regression: completeness-probe coverage for analytic pairs.

The probe must run wherever OCCT uses a numerical walking intersector,
not only for NURBS pairs. Classification rule under test
(intersection._pair_needs_completeness_probe):

- probe ON: freeform pairs (pre-G4 behavior, unchanged), non-coaxial
  quadric/quadric pairs, anything involving a torus or another
  non-quadric surface;
- probe OFF: plane/plane, plane/cylinder, plane/sphere, plane/cone, and
  coaxial quadric pairs (closed form via IntAna).

Coaxiality is conservative on purpose: a missed coaxial verdict only
costs an extra probe run, while a wrong one would skip the probe where
OCCT walks (I1).

The test also checks the behavior end to end through section_face_pair:
a plane/cylinder pair must report raw_curve_count == 0 (probe off) while
a tilted cylinder/cylinder pair must run the probe (raw_curve_count > 0)
and accept with zero unmatched components.
"""
import sys

sys.path.insert(0, "src")

from brepkernel.intersection import (
    _pair_needs_completeness_probe,
    section_face_pair,
)
from brepkernel.step_ingest import index_shape

from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCP.BRepPrimAPI import (
    BRepPrimAPI_MakeBox,
    BRepPrimAPI_MakeCone,
    BRepPrimAPI_MakeCylinder,
    BRepPrimAPI_MakeSphere,
    BRepPrimAPI_MakeTorus,
)
from OCP.gp import gp_Ax1, gp_Ax2, gp_Dir, gp_Pnt, gp_Trsf, gp_Vec

failures = []


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    if not cond:
        failures.append(name)


def pick(shape, want):
    for f in index_shape(shape).faces:
        if want in f.surface_type:
            return f
    raise AssertionError(f"no {want} face found")


def moved(shape, dx=0.0, dy=0.0, dz=0.0):
    t = gp_Trsf()
    t.SetTranslation(gp_Vec(dx, dy, dz))
    return BRepBuilderAPI_Transform(shape, t, True).Shape()


box = BRepPrimAPI_MakeBox(4, 4, 4).Shape()
cyl = BRepPrimAPI_MakeCylinder(
    gp_Ax2(gp_Pnt(0, 0, -2), gp_Dir(0, 0, 1)), 1.0, 6.0).Shape()
sph = BRepPrimAPI_MakeSphere(
    gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 2.0).Shape()
cone = BRepPrimAPI_MakeCone(
    gp_Ax2(gp_Pnt(0, 0, -2), gp_Dir(0, 0, 1)), 1.0, 2.0, 6.0).Shape()
tor = BRepPrimAPI_MakeTorus(
    gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 2.0, 0.7).Shape()

fp = pick(box, "Plane")
fc = pick(cyl, "Cylinder")
fs = pick(sph, "Sphere")
fk = pick(cone, "Cone")
ft = pick(tor, "Torus")

fc_off = pick(moved(cyl, 1.2, 0, 0), "Cylinder")
fs_off = pick(moved(sph, 1.5, 0, 0), "Sphere")

tilt = gp_Trsf()
tilt.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 1, 0)), 0.35)
fc_tilt = pick(BRepBuilderAPI_Transform(cyl, tilt, True).Shape(), "Cylinder")
fk_tilt = pick(BRepBuilderAPI_Transform(cone, tilt, True).Shape(), "Cone")

CASES = [
    # (name, face_a, face_b, probe_expected)
    ("plane/plane", fp, fp, False),
    ("plane/cylinder", fp, fc, False),
    ("plane/sphere", fp, fs, False),
    ("plane/cone", fp, fk, False),
    ("plane/torus", fp, ft, True),
    ("cyl/cyl coaxial", fc, fc, False),
    ("cyl/cone coaxial", fc, fk, False),
    ("cyl/sph coaxial", fc, fs, False),
    ("sph/sph", fs, fs, False),
    ("sph/sph offset", fs, fs_off, False),
    ("cyl/cyl offset", fc, fc_off, True),
    ("cyl/cyl tilted", fc, fc_tilt, True),
    ("cone/cone tilted", fk, fk_tilt, True),
    ("cone/cyl tilted", fk_tilt, fc, True),
    ("cyl/sph offset", fc_off, fs, True),
    ("torus/torus", ft, ft, True),
    ("torus/cyl", ft, fc, True),
]

for name, a, b, want in CASES:
    got = _pair_needs_completeness_probe(a, b, base_tol=1e-7)
    check(f"classify {name}", got == want, f"want={want} got={got}")
    # Classification is symmetric in the pair order.
    got_rev = _pair_needs_completeness_probe(b, a, base_tol=1e-7)
    check(f"classify {name} reversed", got_rev == want,
          f"want={want} got={got_rev}")

# End to end through section_face_pair.
r_off = section_face_pair(fp, fc)
check("plane/cyl probe off end-to-end",
      r_off.raw_curve_count == 0 and r_off.status == "curve",
      f"status={r_off.status} raw={r_off.raw_curve_count}")

r_on = section_face_pair(fc, fc_tilt)
check("tilted cyl/cyl probe on end-to-end",
      r_on.raw_curve_count > 0,
      f"raw={r_on.raw_curve_count}")
check("tilted cyl/cyl probe accepts",
      r_on.raw_unmatched_components == 0,
      f"unmatched={r_on.raw_unmatched_components} "
      f"status={r_on.status}")

print(f"{len(failures)} failures")
sys.exit(1 if failures else 0)
