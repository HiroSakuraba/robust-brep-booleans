"""G10: rigid-motion invariance of boolean_brep() verdicts.

A rigid motion (rotation + translation) applied identically to both
operands must not change the verdict: accept stays accept with the same
volume (1e-9 rel). A refusal may resolve into an accept under rotation
(the plan explicitly allows new accepts: only accept -> refuse flips are
forbidden); every such new accept is verified against the independent
OCCT boolean oracle volume. Refusals that stay refusals must keep their
kind and stage.

Background: plane coincidence was decided by bit-identical canonical
plane parameters. A rigid motion re-rounds those doubles, so two faces
coincident before the motion could read as non-coincident after it and
flip the verdict. The fix makes the offset test representation-relative
(64 ULPs of the location magnitude: what a rigid motion can introduce
through float64 rounding) instead of bit-exact.

25 base cases x 3 rigid motions = 75 checks.
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from brepkernel import boolean_brep, BRepAmbiguousResult  # noqa: E402

from OCP.BRepAlgoAPI import (BRepAlgoAPI_Common, BRepAlgoAPI_Cut,  # noqa: E402
                             BRepAlgoAPI_Fuse)
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform  # noqa: E402
from OCP.BRepGProp import BRepGProp  # noqa: E402
from OCP.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder,  # noqa: E402
                             BRepPrimAPI_MakeSphere)
from OCP.GProp import GProp_GProps  # noqa: E402
from OCP.TopExp import TopExp_Explorer  # noqa: E402
from OCP.TopAbs import TopAbs_SOLID  # noqa: E402
from OCP.gp import gp_Ax1, gp_Dir, gp_Pnt, gp_Trsf, gp_Vec  # noqa: E402


def _box(x0, y0, z0, x1, y1, z1):
    return BRepPrimAPI_MakeBox(gp_Pnt(x0, y0, z0), gp_Pnt(x1, y1, z1)).Shape()


def _cyl(px, py, pz, dx, dy, dz, r, h):
    from OCP.gp import gp_Ax2
    return BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(px, py, pz), gp_Dir(dx, dy, dz)), r, h).Shape()


def _sph(cx, cy, cz, r):
    return BRepPrimAPI_MakeSphere(gp_Pnt(cx, cy, cz), r).Shape()


def _has_solid(shape):
    ex = TopExp_Explorer(shape, TopAbs_SOLID)
    return ex.More()


def _volume(shape):
    if not _has_solid(shape):
        return 0.0
    g = GProp_GProps()
    BRepGProp.VolumePropertiesGK_s(shape, g, 1e-10, True, True, False, False, False)
    return float(g.Mass())


def _rigid(axis, angle_deg, tx, ty, tz):
    t = gp_Trsf()
    n = math.sqrt(sum(c * c for c in axis))
    t.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0),
                         gp_Dir(axis[0] / n, axis[1] / n, axis[2] / n)),
                  math.radians(angle_deg))
    t2 = gp_Trsf()
    t2.SetTranslation(gp_Vec(tx, ty, tz))
    t.Multiply(t2)
    return t


# Three fixed rigid motions: a generic rotation, a second generic
# rotation, and a rotation plus a large translation (stresses
# representation rounding at scale).
MOTIONS = [
    ("r37", _rigid((1, 2, 3), 37, 10, -7, 3)),
    ("r113", _rigid((-2, 1, 0.5), 113, -50, 21, -13)),
    ("r251_big", _rigid((0.3, -0.8, 0.5), 251, 1000, 2000, -1500)),
]


def _move(shape, trsf):
    return BRepBuilderAPI_Transform(shape, trsf, True).Shape()


def _verdict(a, b, op):
    try:
        out, _rep = boolean_brep(a, b, op)
        return ("accept", _volume(out))
    except BRepAmbiguousResult as exc:
        r = (exc.report or {}).get("refusal", {}) or {}
        return ("refuse", r.get("kind"), r.get("stage"))
    except Exception as exc:  # noqa: BLE001
        return ("crash", type(exc).__name__, str(exc)[:120])


# 25 base cases. The first ~18 are coincidence-sensitive (shared faces);
# the rest are controls that must not change either.
BASE_CASES = [
    ("union full shared face",
     lambda: _box(0, 0, 0, 1, 1, 1), lambda: _box(1, 0, 0, 2, 1, 1), "union"),
    ("union partial shared face",
     lambda: _box(0, 0, 0, 1, 1, 1), lambda: _box(1, 0.25, 0.25, 2, 0.75, 0.75), "union"),
    ("union stacked",
     lambda: _box(0, 0, 0, 2, 2, 2), lambda: _box(0.5, 0.5, 2, 1.5, 1.5, 3), "union"),
    ("diff slot flush top",
     lambda: _box(0, 0, 0, 2, 2, 2), lambda: _box(0.5, -1, 1, 1.5, 3, 2), "difference"),
    ("diff pocket flush top",
     lambda: _box(0, 0, 0, 2, 2, 2), lambda: _box(0.5, 0.5, 1, 1.5, 1.5, 2), "difference"),
    ("diff slot overshoot top",
     lambda: _box(0, 0, 0, 2, 2, 2), lambda: _box(0.5, -1, 1, 1.5, 3, 3), "difference"),
    ("diff corner notch 3 faces",
     lambda: _box(0, 0, 0, 2, 2, 2), lambda: _box(1, 1, 1, 2, 2, 2), "difference"),
    ("isect sub-box 3 faces",
     lambda: _box(0, 0, 0, 2, 2, 2), lambda: _box(0, 0, 0, 1, 1, 1), "intersection"),
    ("isect boxes one shared face",
     lambda: _box(0, 0, 0, 1, 1, 1), lambda: _box(1, 0, 0, 2, 1, 1), "intersection"),
    ("diff tool one shared face",
     lambda: _box(0, 0, 0, 1, 1, 1), lambda: _box(1, 0, 0, 2, 1, 1), "difference"),
    ("union edge touch",
     lambda: _box(0, 0, 0, 1, 1, 1), lambda: _box(1, 1, 0, 2, 2, 1), "union"),
    ("union boss on face",
     lambda: _box(0, 0, 0, 2, 2, 2),
     lambda: _cyl(1, 1, 2, 0, 0, 1, 0.5, 1), "union"),
    ("diff blind hole flush",
     lambda: _box(0, 0, 0, 2, 2, 2),
     lambda: _cyl(1, 1, 1, 0, 0, 1, 0.5, 1), "difference"),
    ("diff through hole",
     lambda: _box(0, 0, 0, 2, 2, 2),
     lambda: _cyl(1, 1, -1, 0, 0, 1, 0.5, 4), "difference"),
    ("union cyl side flush",
     lambda: _box(0, 0, 0, 2, 2, 2),
     lambda: _cyl(3, 1, 1, -1, 0, 0, 0.5, 1), "union"),
    ("diff box flush side",
     lambda: _box(0, 0, 0, 2, 2, 2), lambda: _box(2, 0.5, 0.5, 3, 1.5, 1.5), "difference"),
    ("union three boxes chain",
     lambda: _box(0, 0, 0, 3, 1, 1), lambda: _box(1, 0, 0, 2, 1, 1), "union"),
    ("isect identical boxes",
     lambda: _box(0, 0, 0, 1, 1, 1), lambda: _box(0, 0, 0, 1, 1, 1), "intersection"),
    # Controls: no coincidence involved; verdicts must still not move.
    ("union disjoint",
     lambda: _box(0, 0, 0, 1, 1, 1), lambda: _box(5, 5, 5, 6, 6, 6), "union"),
    ("union overlapping",
     lambda: _box(0, 0, 0, 2, 2, 2), lambda: _box(1, 1, 1, 3, 3, 3), "union"),
    ("diff nested",
     lambda: _box(0, 0, 0, 3, 3, 3), lambda: _box(1, 1, 1, 2, 2, 2), "difference"),
    ("isect overlapping",
     lambda: _box(0, 0, 0, 2, 2, 2), lambda: _box(1, 1, 1, 3, 3, 3), "intersection"),
    ("union box sphere overlap",
     lambda: _box(0, 0, 0, 2, 2, 2), lambda: _sph(1, 1, 1, 1.2), "union"),
    ("diff sphere bite",
     lambda: _box(0, 0, 0, 2, 2, 2), lambda: _sph(2.1, 1.3, 0.9, 0.8), "difference"),
    ("union crossing cyls",
     lambda: _cyl(0, 0, -2, 0, 0, 1, 1.0, 4),
     lambda: _cyl(-2, 0, 0, 1, 0, 0, 0.6, 4), "union"),
]


def _occt_oracle_volume(fa, fb, op):
    """Independent oracle volume via OCCT's own booleans on base inputs.

    Returns None when OCCT cannot build the result.
    """
    try:
        a, b = fa(), fb()
        if op == "union":
            mk = BRepAlgoAPI_Fuse(a, b)
        elif op == "difference":
            mk = BRepAlgoAPI_Cut(a, b)
        else:
            mk = BRepAlgoAPI_Common(a, b)
        mk.Build()
        if not mk.IsDone():
            return None
        return _volume(mk.Shape())
    except Exception:  # noqa: BLE001
        return None


def _check_equal(v0, v1, fa, fb, op):
    """Return (err, note); err is None on pass.

    accept -> accept: volumes must agree to 1e-9 relative.
    accept -> refuse: hard failure (verdict regression).
    refuse -> refuse: kind and stage must be stable.
    refuse -> accept: allowed (new accepts), but the accepted volume must
      agree with the independent OCCT oracle to 1e-9 relative.
    """
    if v0[0] == "accept" and v1[0] == "accept":
        x, y = v0[1], v1[1]
        denom = max(abs(x), abs(y), 1e-12)
        if abs(x - y) > 1e-9 * denom:
            return ("volume %r -> %r" % (x, y), None)
        return (None, None)
    if v0[0] == "accept":
        return ("verdict accept -> %s (REGRESSION)" % (v1[0],), None)
    if v1[0] != "accept":
        if v0[1:] != v1[1:]:
            return ("refusal %r -> %r" % (v0[1:], v1[1:]), None)
        return (None, None)
    oracle = _occt_oracle_volume(fa, fb, op)
    if oracle is None:
        return ("new accept but OCCT oracle failed; manual review needed",
                None)
    x, y = v1[1], oracle
    denom = max(abs(x), abs(y), 1e-12)
    if abs(x - y) > 1e-9 * denom:
        return ("new accept vol %r disagrees with OCCT oracle %r" % (x, y),
                None)
    return (None, "new-accept, OCCT oracle vol %.15g" % y)


def main():
    assert len(BASE_CASES) == 25, len(BASE_CASES)
    assert len(MOTIONS) == 3, len(MOTIONS)
    fails = []
    notes = []
    n = 0
    for name, fa, fb, op in BASE_CASES:
        a, b = fa(), fb()
        v0 = _verdict(a, b, op)
        for mname, trsf in MOTIONS:
            n += 1
            v1 = _verdict(_move(a, trsf), _move(b, trsf), op)
            err, note = _check_equal(v0, v1, fa, fb, op)
            status = "ok" if err is None else "FAIL"
            extra = " [%s]" % note if note else ""
            print("%-28s %-8s %-7s -> %-7s %s%s" %
                  (name, mname, v0[0], v1[0], status, extra), flush=True)
            if err is not None:
                fails.append((name, mname, v0, v1, err))
            elif note:
                notes.append((name, mname, note))
    print("checks=%d failures=%d new_accepts=%d" %
          (n, len(fails), len(notes)))
    for name, mname, v0, v1, err in fails:
        print("FAIL %s %s: %s (before=%r after=%r)" % (name, mname, err, v0, v1))
    for name, mname, note in notes:
        print("note %s %s: %s" % (name, mname, note))
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
