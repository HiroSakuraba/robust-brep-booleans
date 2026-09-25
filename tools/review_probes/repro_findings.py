"""Reproducers for review findings F2 (completeness-probe false alarm) and
F4 (OCCT solid classifier false IN). Run from repo root.

F2: NURBS cone INTERSECT NURBS sphere. Baseline behavior: refused with
    SectionCompletenessMismatch, max_distance ~9.41e-5 > tol 4e-5, although the
    raw component is the same branch and OCCT's result volume (~0.54082)
    matches an independent winding-number Monte Carlo estimate (~0.542 +/- 0.03).
    Target after G3: accepted, volume within 1e-6 relative of OCCT, and the
    arbiter finds zero errors. The probe must still catch the omitted torus
    loop in tests/test_nurbs_adversarial_corpus.py.

F4: NURBS sphere UNION NURBS sphere. Baseline: kernel result is correct
    (winding number 0 at the probe point; volume equals vA + vB - vCommon),
    but BRepClass3d_SolidClassifier reports IN for a point more than 1.0 away
    from both input surfaces, on both the kernel result and OCCT's own fuse.
    Target after G4: the kernel's own patch/shell classification never relies
    on this classifier alone.
"""
import os
import sys

HERE = os.path.dirname(__file__)
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, HERE)

import numpy as np  # noqa: E402
from brepkernel import boolean_brep, BRepAmbiguousResult  # noqa: E402
from arbiter import occt_state, surface_distance, triangles, winding  # noqa: E402

from OCP.BRep import BRep_Builder
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Fuse
from OCP.BRepGProp import BRepGProp
from OCP.BRepTools import BRepTools
from OCP.GProp import GProp_GProps
from OCP.TopoDS import TopoDS_Shape

DATA = os.path.join(ROOT, "tests", "data", "review_20260925")


def load(name):
    s = TopoDS_Shape()
    if not BRepTools.Read_s(s, os.path.join(DATA, name), BRep_Builder()):
        raise RuntimeError(f"could not read {name}")
    return s


def vol(s):
    g = GProp_GProps()
    BRepGProp.VolumePropertiesGK_s(s, g, 1e-10, True, True, False, False, False)
    return g.Mass()


def f2():
    a, b = load("completeness_false_alarm_A.brep"), load("completeness_false_alarm_B.brep")
    print("F2 OCCT common volume:", vol(BRepAlgoAPI_Common(a, b).Shape()))
    try:
        out, rep = boolean_brep(a, b, "intersection")
        print("F2 ACCEPTED volume:", vol(out))
        return True
    except BRepAmbiguousResult as e:
        print("F2 REFUSED:", e.report.get("refusal"))
        return False


def f4():
    a, b = load("classifier_false_in_A.brep"), load("classifier_false_in_B.brep")
    p = np.array([1.3322189978645596, 1.1200289002934705, -0.8236295495693953])
    out, _ = boolean_brep(a, b, "union")
    for name, s in (("A", a), ("B", b), ("kernel union", out),
                    ("OCCT fuse", BRepAlgoAPI_Fuse(a, b).Shape())):
        print(f"F4 {name:13s} occt_state={occt_state(s, p)!s:22s} "
              f"winding={winding(triangles(s), p):+.4f} "
              f"surface_distance={surface_distance(s, p):.4f}")


if __name__ == "__main__":
    from arbiter import version_banner  # noqa: E402
    print(version_banner(), flush=True)
    ok = f2()
    f4()
    raise SystemExit(0 if ok else 1)
