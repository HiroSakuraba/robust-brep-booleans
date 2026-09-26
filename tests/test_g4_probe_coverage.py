"""G4 rework: the completeness probe runs for EVERY face pair.

The G4 classification (IntAna closed-form pairs skip the probe) is gone:
it used a claim about OCCT's exactness to skip verifying OCCT. The probe
now runs unconditionally for every candidate pair that produces section
curves, for analytic and freeform pairs alike.

This test checks end to end through section_face_pair that the probe runs
(raw_curve_count > 0) and accepts (zero unmatched components) for
representative analytic pair kinds: plane/plane, plane/cylinder,
tilted cylinder/cylinder, and torus/cylinder. It also checks the
BREPKERNEL_COMPLETENESS_PROBE environment kill switch: with it set to 0
the probe is skipped everywhere (raw_curve_count == 0), and by default
it is on.
"""
import os
import sys

sys.path.insert(0, "src")

from brepkernel.intersection import (
    _completeness_probe_enabled,
    section_face_pair,
)
from brepkernel.step_ingest import index_shape

from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCP.BRepPrimAPI import (
    BRepPrimAPI_MakeBox,
    BRepPrimAPI_MakeCylinder,
    BRepPrimAPI_MakeTorus,
)
from OCP.gp import gp_Ax1, gp_Ax2, gp_Dir, gp_Pnt, gp_Trsf, gp_Vec

failures = []


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    if not cond:
        failures.append(name)


def pick(shape, want, skip=0):
    seen = 0
    for f in index_shape(shape).faces:
        if want in f.surface_type:
            if seen == skip:
                return f
            seen += 1
    raise AssertionError(f"no {want} face found")


box = BRepPrimAPI_MakeBox(4, 4, 4).Shape()
cyl = BRepPrimAPI_MakeCylinder(
    gp_Ax2(gp_Pnt(0, 0, -2), gp_Dir(0, 0, 1)), 1.0, 6.0).Shape()
tor = BRepPrimAPI_MakeTorus(
    gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 2.0, 0.7).Shape()
# Shift the cylinder through the torus tube so the pair intersects.
t = gp_Trsf()
t.SetTranslation(gp_Vec(2.0, 0, 0))
cyl_in_tor = BRepBuilderAPI_Transform(cyl, t, True).Shape()

# Two intersecting (non-parallel, adjacent) box planes: search the
# deterministic face order for a pair whose section is a curve.
planes = [f for f in index_shape(box).faces if "Plane" in f.surface_type]
fp0, fp1 = None, None
for i in range(len(planes)):
    for j in range(i + 1, len(planes)):
        r = section_face_pair(planes[i], planes[j])
        if r.status == "curve" and r.raw_curve_count > 0:
            fp0, fp1 = planes[i], planes[j]
            break
    if fp0 is not None:
        break
assert fp0 is not None, "no intersecting plane pair on the box"
fc = pick(cyl, "Cylinder")
ft = pick(tor, "Torus")
fc_in_tor = pick(cyl_in_tor, "Cylinder")

tilt = gp_Trsf()
tilt.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 1, 0)), 0.35)
fc_tilt = pick(BRepBuilderAPI_Transform(cyl, tilt, True).Shape(), "Cylinder")

PAIRS = [
    ("plane/plane", fp0, fp1),
    ("plane/cylinder", fp0, fc),
    ("cyl/cyl tilted", fc, fc_tilt),
    ("torus/cyl", ft, fc_in_tor),
]

for name, a, b in PAIRS:
    r = section_face_pair(a, b)
    check(f"{name}: probe ran",
          r.raw_curve_count > 0,
          f"status={r.status} raw={r.raw_curve_count}")
    check(f"{name}: probe accepts",
          r.raw_unmatched_components == 0,
          f"unmatched={r.raw_unmatched_components} "
          f"status={r.status}")

# The kill switch defaults to on ...
check("probe enabled by default", _completeness_probe_enabled(),
      f"env={os.environ.get('BREPKERNEL_COMPLETENESS_PROBE')!r}")
# ... and BREPKERNEL_COMPLETENESS_PROBE=0 disables it everywhere.
os.environ["BREPKERNEL_COMPLETENESS_PROBE"] = "0"
try:
    check("probe disabled by env switch",
          not _completeness_probe_enabled(), "")
    r = section_face_pair(fc, fc_tilt)
    check("env switch skips probe end-to-end",
          r.raw_curve_count == 0 and r.status.startswith("curve"),
          f"status={r.status} raw={r.raw_curve_count}")
finally:
    del os.environ["BREPKERNEL_COMPLETENESS_PROBE"]
check("probe re-enabled after env restore",
      _completeness_probe_enabled(), "")

print(f"{len(failures)} failures")
sys.exit(1 if failures else 0)
