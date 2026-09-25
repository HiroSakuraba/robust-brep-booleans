"""Reproducers for review findings F2 (completeness-probe false alarm) and
F4 (OCCT solid classifier false IN). Run from repo root.

F2: NURBS cone INTERSECT NURBS sphere. Baseline behavior: refused with
    SectionCompletenessMismatch, max_distance ~9.41e-5 > tol 4e-5, although the
    raw component is the same branch and OCCT's result volume (~0.54082)
    matches an independent winding-number Monte Carlo estimate (~0.542 +/- 0.03).
    Target after G3: accepted, volume within 1e-6 relative of OCCT, and the
    arbiter finds zero errors; on accept the probe asserts both, using a
    per-trial seeded RNG for the arbiter. On typed refusal the probe records
    the refusal and returns False (no assertions run). The probe must still
    catch the omitted torus loop in tests/test_nurbs_adversarial_corpus.py.

F4: NURBS sphere UNION NURBS sphere. Baseline: kernel result is correct
    (winding number 0 at the probe point; volume equals vA + vB - vCommon),
    but BRepClass3d_SolidClassifier reports IN for a point more than 1.0 away
    from both input surfaces, on both the kernel result and OCCT's own fuse.
    The probe asserts the kernel's own verdict (winding number OUT at the
    probe point, surface distance above 1.0) and records the OCCT classifier
    state at the probe point for A, B, the kernel union, and the OCCT fuse.
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
from arbiter import (occt_state, surface_distance, triangles, winding,  # noqa: E402
                    membership_audit)

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
    oracle_vol = vol(BRepAlgoAPI_Common(a, b).Shape())
    print("F2 OCCT common volume:", oracle_vol)
    try:
        out, rep = boolean_brep(a, b, "intersection")
    except BRepAmbiguousResult as e:
        print("F2 REFUSED:", e.report.get("refusal"))
        return False
    kv = vol(out)
    rel = abs(kv - oracle_vol) / max(1.0, abs(oracle_vol))
    print(f"F2 ACCEPTED volume: {kv}  relerr_vs_occt={rel:.3g}")
    assert rel <= 1e-6, f"F2 volume {kv} not within 1e-6 of OCCT {oracle_vol}"
    # Seeded per-trial RNG so the audit is reproducible.
    aud = membership_audit(a, b, out, "intersection",
                           np.random.default_rng(20260925))
    print(f"F2 arbiter: checked={aud['checked']} "
          f"kernel_errors={aud['kernel_errors']} "
          f"classifier_disagreements={aud['classifier_disagreements']}")
    if aud["error_points"]:
        print("F2 arbiter error_points:", aud["error_points"])
    assert aud["kernel_errors"] == 0, \
        f"F2 arbiter kernel_errors={aud['kernel_errors']}"
    return True


def f4():
    a, b = load("classifier_false_in_A.brep"), load("classifier_false_in_B.brep")
    p = np.array([1.3322189978645596, 1.1200289002934705, -0.8236295495693953])
    out, _ = boolean_brep(a, b, "union")
    rows = []
    for name, s in (("A", a), ("B", b), ("kernel union", out),
                    ("OCCT fuse", BRepAlgoAPI_Fuse(a, b).Shape())):
        rows.append((name, str(occt_state(s, p)), winding(triangles(s), p),
                     surface_distance(s, p)))
    for name, st, w, d in rows:
        print(f"F4 {name:13s} occt_state={st:22s} "
              f"winding={w:+.4f} surface_distance={d:.4f}")
    # rows[2] is the kernel union row: (name, occt_state, winding, distance).
    w, d = rows[2][2], rows[2][3]
    assert abs(w) < 0.05, f"F4 kernel union winding verdict at p is not OUT: w={w}"
    assert d > 1.0, f"F4 kernel union surface_distance at p not > 1.0: d={d}"


if __name__ == "__main__":
    from arbiter import version_banner  # noqa: E402
    print(version_banner(), flush=True)
    ok = f2()
    f4()
    raise SystemExit(0 if ok else 1)
