"""Planted-error audit for a possible 'consensus' certification mode.

    python tools/review_probes/planted_error_audit.py [N_CASES] [N_POINTS]

Prints one RESULT line (JSON): per planted error type, how many of the
N_CASES were caught by each check and by any check. Measured 2026-10-13
(25 cases, 300 points): 274/275 caught, 0 false alarms; see the T1-T6
entry in docs/REVIEW_LEDGER.md.

Question: if an outside engine hands us a result, would generic,
construction-independent checks catch it when it is wrong?

Checks (none uses brepkernel's own construction evidence):
  valid      BRepCheck valid, closed manifold solids
  volume     brepkernel's operation volume bounds (_operation_volume_bounds)
  membership arbiter.membership_audit, n points, 0 errors required
  boundary   points on every result face lie on A's or B's surface

Planted errors on correct results R (brepkernel-accepted):
  wrong_op         the result of a different operation
  bump_s / dent_s  R union / minus a sphere of radius s*L at a boundary point
  shift_d          R translated by d*L
"""
import json
import sys

import numpy as np

_HERE = __import__("os").path.dirname(__import__("os").path.abspath(__file__))
sys.path.insert(0, __import__("os").path.join(_HERE, "..", "..", "src"))
sys.path.insert(0, _HERE)
import arbiter  # noqa: E402
from brepkernel import BRepAmbiguousResult, boolean_brep  # noqa: E402
from brepkernel.assembly import _face_points, _shape_volume, _all_faces_analytic  # noqa: E402
from brepkernel.pipeline import _operation_volume_bounds  # noqa: E402
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse  # noqa: E402
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform  # noqa: E402
from OCP.BRepCheck import BRepCheck_Analyzer  # noqa: E402
from OCP.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCone,  # noqa: E402
                             BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeSphere,
                             BRepPrimAPI_MakeTorus)
from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID  # noqa: E402
from OCP.TopExp import TopExp_Explorer  # noqa: E402
from OCP.TopoDS import TopoDS  # noqa: E402
from OCP.gp import gp_Ax1, gp_Ax2, gp_Dir, gp_Pnt, gp_Trsf, gp_Vec  # noqa: E402

BOP = {"union": BRepAlgoAPI_Fuse, "intersection": BRepAlgoAPI_Common,
       "difference": BRepAlgoAPI_Cut}


def vol(s):
    if not TopExp_Explorer(s, TopAbs_SOLID).More():
        return 0.0
    return abs(_shape_volume(s, 1e-12 if _all_faces_analytic(s) else 1e-10))


def size(*shapes):
    lo, hi = arbiter.combined_domain(*shapes, pad_frac=0.0)
    return float(np.linalg.norm(np.asarray(hi) - np.asarray(lo)))


def check_valid(r):
    return (TopExp_Explorer(r, TopAbs_SOLID).More()
            and BRepCheck_Analyzer(r).IsValid())


def check_volume(a, b, r, op, L):
    ok, _bounds, _tol = _operation_volume_bounds(op, vol(a), vol(b), vol(r), max(L, 1.0), 1e-7)
    return bool(ok)


def check_membership(a, b, r, op, seed, n):
    res = arbiter.membership_audit(a, b, r, op, np.random.default_rng(seed), n=n)
    return res["kernel_errors"] == 0


def check_boundary(a, b, r, L):
    bound = 1e-5 * L
    ex = TopExp_Explorer(r, TopAbs_FACE)
    while ex.More():
        f = TopoDS.Face(ex.Current())
        ex.Next()
        try:
            pts = _face_points(f, 1e-7, max_points=3, min_points=1)
        except Exception:
            continue
        for p in pts:
            q = np.asarray(p)
            if min(arbiter.surface_distance(a, q), arbiter.surface_distance(b, q)) > bound:
                return False
    return True


def boundary_point(r, rng):
    faces = []
    ex = TopExp_Explorer(r, TopAbs_FACE)
    while ex.More():
        faces.append(TopoDS.Face(ex.Current()))
        ex.Next()
    f = faces[int(rng.integers(len(faces)))]
    pts = _face_points(f, 1e-7, max_points=3, min_points=1)
    return pts[0]


def plant(a, b, r, op, kind, mag, rng, L):
    if kind == "wrong_op":
        other = [o for o in BOP if o != op][int(rng.integers(2))]
        return BOP[other](a, b).Shape()
    if kind in ("bump", "dent"):
        c = boundary_point(r, rng)
        s = BRepPrimAPI_MakeSphere(gp_Pnt(*map(float, c)), mag * L).Shape()
        return (BRepAlgoAPI_Fuse if kind == "bump" else BRepAlgoAPI_Cut)(r, s).Shape()
    if kind == "shift":
        d = rng.normal(size=3)
        d = d / np.linalg.norm(d) * mag * L
        t = gp_Trsf()
        t.SetTranslation(gp_Vec(*map(float, d)))
        return BRepBuilderAPI_Transform(r, t, True).Shape()
    raise ValueError(kind)


def cases(rng, n):
    def place(s):
        t = gp_Trsf()
        ax = rng.normal(size=3)
        ax /= np.linalg.norm(ax)
        t.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(*map(float, ax))), float(rng.uniform(0, 6.283)))
        s = BRepBuilderAPI_Transform(s, t, True).Shape()
        t2 = gp_Trsf()
        t2.SetTranslation(gp_Vec(*map(float, rng.uniform(-.6, .6, 3))))
        return BRepBuilderAPI_Transform(s, t2, True).Shape()

    def prim(k):
        if k == 0:
            d = rng.uniform(.6, 1.8, 3)
            return place(BRepPrimAPI_MakeBox(gp_Pnt(*(-d / 2)), gp_Pnt(*(d / 2))).Shape())
        if k == 1:
            r, h = rng.uniform(.3, .9), rng.uniform(.8, 2)
            return place(BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -h / 2), gp_Dir(0, 0, 1)), r, h).Shape())
        if k == 2:
            return place(BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), rng.uniform(.4, 1)).Shape())
        if k == 3:
            h = rng.uniform(.8, 1.8)
            return place(BRepPrimAPI_MakeCone(gp_Ax2(gp_Pnt(0, 0, -h / 2), gp_Dir(0, 0, 1)),
                                              rng.uniform(.4, .9), rng.uniform(0, .3), h).Shape())
        return place(BRepPrimAPI_MakeTorus(rng.uniform(.6, 1.0), rng.uniform(.15, .35)).Shape())
    out = []
    while len(out) < n:
        a, b = prim(int(rng.integers(5))), prim(int(rng.integers(5)))
        op = ["union", "intersection", "difference"][int(rng.integers(3))]
        try:
            r, _ = boolean_brep(a, b, op)
        except BRepAmbiguousResult:
            continue
        if vol(r) < 1e-3:
            continue
        out.append((a, b, op, r))
    return out


def main():
    n_cases = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    n_pts = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    rng = np.random.default_rng(20261013)
    plants = [("none", 0), ("wrong_op", 0), ("bump", 0.2), ("bump", 0.05), ("bump", 0.01),
              ("dent", 0.2), ("dent", 0.05), ("dent", 0.01),
              ("shift", 1e-2), ("shift", 1e-3), ("shift", 1e-4)]
    table = {f"{k}{'_' + str(m) if m else ''}": {"n": 0, "valid": 0, "volume": 0, "membership": 0,
                                                 "boundary": 0, "any": 0} for k, m in plants}
    for ci, (a, b, op, r) in enumerate(cases(rng, n_cases)):
        L = size(a, b)
        for kind, mag in plants:
            key = f"{kind}{'_' + str(mag) if mag else ''}"
            try:
                bad = r if kind == "none" else plant(a, b, r, op, kind, mag, rng, L)
            except Exception as exc:
                print("plant failed", key, exc, flush=True)
                continue
            caught = {}
            caught["valid"] = not check_valid(bad)
            caught["volume"] = (not caught["valid"]) and not check_volume(a, b, bad, op, L)
            caught["membership"] = (not caught["valid"]) and not check_membership(a, b, bad, op, ci, n_pts)
            caught["boundary"] = (not caught["valid"]) and not check_boundary(a, b, bad, L)
            if caught["valid"]:
                caught["volume"] = caught["membership"] = caught["boundary"] = True
            t = table[key]
            t["n"] += 1
            for k in ("valid", "volume", "membership", "boundary"):
                t[k] += int(caught[k])
            t["any"] += int(any(caught.values()))
        print(f"case {ci} done", flush=True)
    print("RESULT", json.dumps(table))


if __name__ == "__main__":
    main()
