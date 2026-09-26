"""G10 Problem B: the optional center-of-mass witness must be cheap and safe.

Finding interior witness points for shell-nesting decisions uses a
high-accuracy adaptive volume integration solely to obtain a
center-of-mass candidate. That integration must (a) not run at all when
the ordinary face-based witness strategy already produced enough
points, (b) run at reduced accuracy when it does run, and (c) never
turn an otherwise certifiable operation into a VolumeIntegrationFailed
refusal when the integration itself fails: a failed optional witness
is skipped, not fatal.
"""
import os
import sys
from unittest import mock

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import OCP.BRepClass  # noqa: E402
import OCP.BRepGProp  # noqa: E402
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox  # noqa: E402
from OCP.TopAbs import TopAbs_OUT  # noqa: E402
from OCP.gp import gp_Pnt  # noqa: E402

from brepkernel.assembly import AssemblyError, _solid_interior_points  # noqa: E402


def _box():
    return BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), gp_Pnt(1, 1, 1)).Shape()


class _AlwaysOutFaceClassifier:
    """Starve the face-based witness loop: no UV sample reads as inside."""

    def __init__(self, *args, **kwargs):
        pass

    def State(self):
        return TopAbs_OUT


def _wrap_gk(impl):
    real = OCP.BRepGProp.BRepGProp.VolumePropertiesGK_s
    calls = []

    def wrapper(*args, **kwargs):
        calls.append(args)
        return impl(real, args, kwargs)

    return wrapper, calls


def _t1_skip_when_face_witnesses_suffice(fails):
    def impl(real, args, kwargs):
        return real(*args, **kwargs)

    wrapper, calls = _wrap_gk(impl)
    with mock.patch.object(
            OCP.BRepGProp.BRepGProp, "VolumePropertiesGK_s", wrapper):
        pts = _solid_interior_points(_box(), 1e-7)
    if len(pts) < 3:
        fails.append("t1: only %d witness points on a plain box" % len(pts))
    if calls:
        fails.append("t1: expensive COM integration ran %d time(s) although "
                     "face witnesses already sufficed" % len(calls))
    print("t1 skip-when-sufficient: %s" %
          ("FAIL" if any(f.startswith("t1") for f in fails) else "ok"),
          flush=True)


def _t2_failure_is_skipped_not_fatal(fails):
    def impl(real, args, kwargs):
        return -1.0  # integration reports failure

    wrapper, calls = _wrap_gk(impl)
    try:
        with mock.patch.object(
                OCP.BRepGProp.BRepGProp, "VolumePropertiesGK_s", wrapper), \
             mock.patch.object(
                OCP.BRepClass, "BRepClass_FaceClassifier",
                _AlwaysOutFaceClassifier):
            # Face loop starved -> COM path executes -> integration fails ->
            # must NOT raise VolumeIntegrationFailed; the frac-grid fallback
            # still certifies this ordinary box.
            pts = _solid_interior_points(_box(), 1e-7)
    except AssemblyError as exc:
        fails.append("t2: raised %s (%s) although the box is certifiable" %
                     (exc.kind, exc))
        print("t2 failure-skipped: FAIL", flush=True)
        return
    if len(calls) != 1:
        fails.append("t2: COM fallback did not execute (calls=%d)" %
                     len(calls))
    if len(pts) < 3:
        fails.append("t2: only %d witness points after skipped COM" %
                     len(pts))
    print("t2 failure-skipped: %s" %
          ("FAIL" if any(f.startswith("t2") for f in fails) else "ok"),
          flush=True)


def _t3_reduced_accuracy(fails):
    def impl(real, args, kwargs):
        return real(*args, **kwargs)

    wrapper, calls = _wrap_gk(impl)
    with mock.patch.object(
            OCP.BRepGProp.BRepGProp, "VolumePropertiesGK_s", wrapper), \
         mock.patch.object(
            OCP.BRepClass, "BRepClass_FaceClassifier",
            _AlwaysOutFaceClassifier):
        pts = _solid_interior_points(_box(), 1e-7)
    if len(calls) != 1:
        fails.append("t3: COM fallback did not execute (calls=%d)" %
                     len(calls))
    else:
        eps = calls[0][2]
        if abs(eps - 1e-4) > 1e-12:
            fails.append("t3: COM witness integration ran at eps=%r, want "
                         "1e-4" % (eps,))
    if len(pts) < 3:
        fails.append("t3: only %d witness points" % len(pts))
    print("t3 reduced-accuracy: %s" %
          ("FAIL" if any(f.startswith("t3") for f in fails) else "ok"),
          flush=True)


def _t4_no_volume_integration_failed_kind(fails):
    def impl(real, args, kwargs):
        return -1.0

    wrapper, _ = _wrap_gk(impl)
    try:
        with mock.patch.object(
                OCP.BRepGProp.BRepGProp, "VolumePropertiesGK_s", wrapper), \
             mock.patch.object(
                OCP.BRepClass, "BRepClass_FaceClassifier",
                _AlwaysOutFaceClassifier):
            _solid_interior_points(_box(), 1e-7)
    except AssemblyError as exc:
        if exc.kind == "VolumeIntegrationFailed":
            fails.append("t4: optional COM witness failure raised "
                         "VolumeIntegrationFailed")
    print("t4 no-fatal-kind: %s" %
          ("FAIL" if any(f.startswith("t4") for f in fails) else "ok"),
          flush=True)


def main():
    fails = []
    _t1_skip_when_face_witnesses_suffice(fails)
    _t2_failure_is_skipped_not_fatal(fails)
    _t3_reduced_accuracy(fails)
    _t4_no_volume_integration_failed_kind(fails)
    print("failures=%d" % len(fails))
    for f in fails:
        print("FAIL " + f)
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
