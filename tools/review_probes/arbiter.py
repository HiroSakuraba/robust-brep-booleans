"""Engine-independent point-membership arbiter for Boolean results.

Two independent inside/outside tests are provided:

  * occt_state(): OCCT BRepClass3d_SolidClassifier (fast, but the kernel
    itself also depends on it, so it is NOT independent evidence).
  * winding(): generalized winding number of a fine OCCT tessellation.
    Independent of every Boolean/intersection algorithm. Only trustworthy
    for points farther from the surface than the tessellation deflection;
    callers must skip points inside that band.

A Boolean result is judged wrong at a point only when the winding-number
verdict on the RESULT disagrees with the set-operation of the winding-number
verdicts on the INPUTS. OCCT-classifier disagreements are reported
separately as "classifier_disagreements" and are not counted as kernel
errors (see review finding F4: the OCCT classifier produced a false IN).
"""
from __future__ import annotations

import numpy as np


def version_banner() -> str:
    """One-line version summary for embedding at the top of probe reports."""
    try:
        import importlib.metadata as md
        ocp = md.version("cadquery-ocp")
        mfd = md.version("manifold3d")
    except Exception:
        ocp, mfd = "unknown", "unknown"
    return f"versions: numpy={np.__version__} cadquery-ocp={ocp} manifold3d={mfd}"

from OCP.BRep import BRep_Tool
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.BRepExtrema import BRepExtrema_DistShapeShape
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.BRepTools import BRepTools
from OCP.TopAbs import (TopAbs_FACE, TopAbs_IN, TopAbs_ON, TopAbs_OUT,
                        TopAbs_REVERSED, TopAbs_SOLID)
from OCP.TopExp import TopExp_Explorer
from OCP.TopLoc import TopLoc_Location
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Pnt


def has_solid(shape) -> bool:
    return (not shape.IsNull()) and TopExp_Explorer(shape, TopAbs_SOLID).More()


def triangles(shape, deflection: float = 2e-4) -> np.ndarray:
    """Oriented triangle soup (n, 3, 3) of a shape's faces."""
    BRepTools.Clean_s(shape)
    BRepMesh_IncrementalMesh(shape, deflection, False, 0.1, True)
    out = []
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        face = TopoDS.Face(ex.Current())
        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation_s(face, loc)
        if tri is None:
            raise RuntimeError("face has no triangulation")
        trsf = loc.Transformation()
        rev = face.Orientation() == TopAbs_REVERSED
        for k in range(1, tri.NbTriangles() + 1):
            i1, i2, i3 = tri.Triangle(k).Get()
            p = [np.array(tri.Node(j).Transformed(trsf).Coord()) for j in (i1, i2, i3)]
            out.append([p[0], p[2], p[1]] if rev else p)
        ex.Next()
    return np.asarray(out, dtype=np.float64).reshape(-1, 3, 3)


def winding(tris: np.ndarray, q: np.ndarray) -> float:
    """Generalized winding number (Van Oosterom-Strackee solid angles)."""
    if len(tris) == 0:
        return 0.0
    a = tris[:, 0] - q
    b = tris[:, 1] - q
    c = tris[:, 2] - q
    la, lb, lc = (np.linalg.norm(x, axis=1) for x in (a, b, c))
    num = np.einsum("ij,ij->i", a, np.cross(b, c))
    den = (la * lb * lc + np.einsum("ij,ij->i", a, b) * lc
           + np.einsum("ij,ij->i", b, c) * la + np.einsum("ij,ij->i", c, a) * lb)
    return float(np.sum(2.0 * np.arctan2(num, den)) / (4.0 * np.pi))


def occt_state(shape, p, tol: float = 1e-7):
    c = BRepClass3d_SolidClassifier(shape)
    c.Perform(gp_Pnt(*map(float, p)), tol)
    return c.State()


def surface_distance(shape, p) -> float:
    """Distance from p to the shape's faces (not to the solid volume)."""
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Compound
    comp = TopoDS_Compound()
    bld = BRep_Builder()
    bld.MakeCompound(comp)
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        bld.Add(comp, ex.Current())
        ex.Next()
    v = BRepBuilderAPI_MakeVertex(gp_Pnt(*map(float, p))).Vertex()
    d = BRepExtrema_DistShapeShape(v, comp)
    return float(d.Value())


def want(op: str, ia: bool, ib: bool) -> bool:
    return {"union": ia or ib, "intersection": ia and ib,
            "difference": ia and not ib}[op]


def membership_audit(a, b, out, op, rng, n=300, lo=-2.2, hi=2.2,
                     band=2e-3, deflection=2e-4) -> dict:
    """Return counts of kernel errors and OCCT-classifier disagreements."""
    ta, tb = triangles(a, deflection), triangles(b, deflection)
    to = triangles(out, deflection) if has_solid(out) else np.zeros((0, 3, 3))
    res = {"checked": 0, "kernel_errors": 0, "classifier_disagreements": 0,
           "error_points": [], "classifier_points": []}
    for q in rng.uniform(lo, hi, size=(n, 3)):
        wa, wb, wo = winding(ta, q), winding(tb, q), winding(to, q)
        # Skip points whose winding number is not decisively 0 or 1
        # (inside the tessellation band) and points near any input surface.
        if any(0.05 < abs(w) < 0.95 for w in (wa, wb, wo)):
            continue
        if surface_distance(a, q) < band or surface_distance(b, q) < band:
            continue
        ia, ib, io = abs(wa) > 0.5, abs(wb) > 0.5, abs(wo) > 0.5
        res["checked"] += 1
        if io != want(op, ia, ib):
            res["kernel_errors"] += 1
            res["error_points"].append(q.tolist())
        if has_solid(out):
            st = occt_state(out, q)
            if st in (TopAbs_IN, TopAbs_OUT) and (st == TopAbs_IN) != io:
                res["classifier_disagreements"] += 1
                res["classifier_points"].append(q.tolist())
    return res
