"""G5 regression: classifier independence (review finding F4).

F4: BRepClass3d_SolidClassifier returned a false IN for a point more than
1.0 from both input surfaces, and the kernel's patch/shell classification
depended on that same classifier. G5 adds a second, independent
multi-ray-parity classifier and requires both classifiers to agree for
every keep/discard decision that rests on point classification.

t1 (fault injection, committed failing pre-G5): monkeypatch the OCCT
solid classifier to flip every IN<->OUT verdict. A single trusted
classifier would then silently keep the wrong patches and accept a wrong
solid. The test requires the operation to REFUSE with the typed kind
ClassifierDisagreement, carrying both verdicts and the probe point.
"""
import sys
import os

sys.path.insert(0, "src")

import numpy as np

from brepkernel import boolean_brep, BRepAmbiguousResult
from brepkernel.assembly import classify_point_two_classifier
from brepkernel.step_ingest import index_shape

import OCP.BRepClass3d as _bcm
from OCP.BRep import BRep_Builder
from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.BRepTools import BRepTools
from OCP.GProp import GProp_GProps
from OCP.TopAbs import TopAbs_IN, TopAbs_OUT
from OCP.TopoDS import TopoDS_Shape
from OCP.gp import gp_Pnt


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def volume(shape):
    g = GProp_GProps()
    BRepGProp.VolumePropertiesGK_s(shape, g, 1e-10, True, True, False,
                                   False, False)
    return float(g.Mass())


def overlapping_boxes():
    a = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 2.0, 2.0, 2.0).Shape()
    b = BRepPrimAPI_MakeBox(gp_Pnt(1, 1, 1), 2.0, 2.0, 2.0).Shape()
    return a, b


class _FlippedClassifier(_bcm.BRepClass3d_SolidClassifier):
    """Fault injection: every material verdict is inverted.

    Subclassing keeps the real OCCT geometry engine underneath and only
    flips the reported TopAbs state, which is exactly the F4 failure mode
    (a wrong verdict from an otherwise working classifier).
    """

    def State(self):
        st = super().State()
        if st == TopAbs_IN:
            return TopAbs_OUT
        if st == TopAbs_OUT:
            return TopAbs_IN
        return st


def _install_flip():
    # assembly.py imports the classifier from the OCP.BRepClass3d module at
    # call time, so replacing the module attribute redirects every
    # classification site (patch witnesses, shell witnesses, nesting).
    orig = _bcm.BRepClass3d_SolidClassifier
    _bcm.BRepClass3d_SolidClassifier = _FlippedClassifier
    return orig


def _restore(orig):
    _bcm.BRepClass3d_SolidClassifier = orig


def t1_flipped_classifier_refuses_typed():
    """A lying OCCT classifier must cause refusal, never a wrong solid."""
    a, b = overlapping_boxes()
    out, _ = boolean_brep(a, b, "union")
    v_expected = volume(out)
    ok = check("t1 baseline union accepted with sane volume",
               abs(v_expected - 15.0) < 1e-6,
               f"volume={v_expected:.6f}")

    orig = _install_flip()
    try:
        out2, _ = boolean_brep(a, b, "union")
    except BRepAmbiguousResult as exc:
        refusal = exc.report.get("refusal", {})
        kind = refusal.get("kind")
        ok &= check("t1 flipped classifier refuses with ClassifierDisagreement",
                    kind == "ClassifierDisagreement",
                    f"kind={kind}")
        cause = exc.cause
        has_verdicts = (getattr(cause, "occt_verdict", None) in
                        ("inside", "outside") and
                        getattr(cause, "independent_verdict", None) in
                        ("inside", "outside"))
        pt = getattr(cause, "point", None)
        has_point = (pt is not None and len(tuple(pt)) == 3 and
                     all(np.isfinite(tuple(pt))))
        ok &= check("t1 refusal carries both verdicts and the point",
                    has_verdicts and has_point,
                    f"occt={getattr(cause, 'occt_verdict', None)} "
                    f"independent={getattr(cause, 'independent_verdict', None)} "
                    f"point={pt}")
        return ok
    finally:
        _restore(orig)

    v_wrong = volume(out2)
    return check("t1 flipped classifier must not accept a wrong solid",
                 False,
                 f"silently accepted volume={v_wrong:.6f} "
                 f"(expected {v_expected:.6f})")


def t2_f4_probe_point_dual_verdict():
    """F4 repro: the kernel's independent classifier says OUT at the probe point.

    Documented F4 probe point (review_20260925 fixtures): more than 1.0
    from both input surfaces, true winding number 0, yet
    BRepClass3d_SolidClassifier reports IN on the union result.  The test
    asserts the kernel's own two-classifier machinery:
    - the independent multi-ray parity verdict at the probe point is OUT;
    - the F4 false IN still reproduces in the OCCT classifier (so the test
      is not vacuous);
    - the pair withholds agreement, so a production decision resting on
      this point would refuse with ClassifierDisagreement instead of
      trusting the false IN;
    - the F4 union itself is still accepted with the correct volume (the
      probe point is not a decision witness; the kernel never classifies
      result points in production).
    """
    data = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "data", "review_20260925")

    def load(name):
        s = TopoDS_Shape()
        if not BRepTools.Read_s(s, os.path.join(data, name),
                                BRep_Builder()):
            raise RuntimeError(f"could not read fixture {name}")
        return s

    a = load("classifier_false_in_A.brep")
    b = load("classifier_false_in_B.brep")
    p = np.array([1.3322189978645596, 1.1200289002934705,
                  -0.8236295495693953])
    out, _ = boolean_brep(a, b, "union")
    v_out = volume(out)
    fuse = BRepAlgoAPI_Fuse(a, b)
    fuse.Build()
    assert fuse.IsDone(), "OCCT fuse failed on F4 fixtures"
    v_expected = volume(fuse.Shape())
    ok = check("t2 F4 union accepted with correct volume",
               abs(v_out - v_expected) <= 1e-6 * max(v_expected, 1.0),
               f"volume={v_out:.6f} expected={v_expected:.6f}")

    model = index_shape(out)
    r = classify_point_two_classifier(p, model, 1e-7)
    ok &= check("t2 independent classifier verdict is OUT at probe point",
                r["independent"] == "outside",
                f"independent={r['independent']}")
    ok &= check("t2 F4 false IN still reproduces in the OCCT classifier",
                r["occt"] == "inside",
                f"occt={r['occt']}")
    ok &= check("t2 pair withholds agreement (no agreed IN)",
                r["agreed"] is False and r["decision"] is None,
                f"agreed={r['agreed']} decision={r['decision']}")
    return ok


def main():
    ok = True
    ok &= t1_flipped_classifier_refuses_typed()
    ok &= t2_f4_probe_point_dual_verdict()
    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
