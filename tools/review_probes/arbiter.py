"""Independent Boolean arbiter: point-membership verdicts for Boolean results.

Two inside/outside tests are provided:

  * occt_state(): OCCT BRepClass3d_SolidClassifier (fast, but the kernel
    itself also depends on it, so it is NOT independent evidence).
  * winding(): generalized winding number of a fine OCCT tessellation.
    Independent of OCCT's Boolean/intersection decision paths (it never
    calls the intersector or the Boolean operators). It still tessellates
    with OCCT's BRepMesh, so it is not independent of OCCT entirely; call
    this an independent Boolean arbiter, not an independent kernel
    arbiter. Only trustworthy for points farther from the surface than
    the tessellation deflection; callers must skip points inside that
    band.

A Boolean result is judged wrong at a point only when the winding-number
verdict on the RESULT disagrees with the set-operation of the winding-number
verdicts on the INPUTS. OCCT-classifier disagreements are reported
separately as "classifier_disagreements" and are not counted as kernel
errors (see review finding F4: the OCCT classifier produced a false IN).

Scale awareness: the sample domain is derived from the combined bounding
box of the two operands and the result (never a hard-coded cube), and the
surface-exclusion band is derived from the mesh deflection plus the
entities' tolerances (never an absolute constant). The arbiter never
mutates its input shapes: tessellation runs on a deep copy.
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

# OCP 7.8 spells TopoDS casts with a trailing _s; reuse the kernel's shim
# when it is importable so this testing tool runs on both OCP versions.
try:
    import brepkernel._occt_compat  # noqa: F401
except ImportError:
    from OCP.TopoDS import TopoDS as _T
    for _n in ("Vertex", "Edge", "Wire", "Face", "Shell", "Solid", "Compound"):
        if not hasattr(_T, _n) and hasattr(_T, _n + "_s"):
            setattr(_T, _n, getattr(_T, _n + "_s"))
from OCP.Bnd import Bnd_Box
from OCP.BRep import BRep_Tool
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepBuilderAPI import BRepBuilderAPI_Copy, BRepBuilderAPI_MakeVertex
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.BRepExtrema import BRepExtrema_DistShapeShape
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.BRepTools import BRepTools
from OCP.TopAbs import (TopAbs_EDGE, TopAbs_FACE, TopAbs_IN, TopAbs_ON,
                        TopAbs_OUT, TopAbs_REVERSED, TopAbs_SOLID,
                        TopAbs_VERTEX)
from OCP.TopExp import TopExp_Explorer
from OCP.TopLoc import TopLoc_Location
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Pnt

# Fraction of the largest combined-bbox edge used as the sampling-domain
# margin on every side. Documented, not tuned to any test: it only needs to
# be wide enough that near-surface sample points exist around the parts.
DOMAIN_PAD_FRAC = 0.1

# Legacy hard-coded cube, kept ONLY as the documented fallback when every
# input shape is null or has a void bounding box (nothing to audit then).
FALLBACK_DOMAIN = ((-2.2, -2.2, -2.2), (2.2, 2.2, 2.2))

# Safety factor on the mesh deflection in the exclusion band: one
# deflection covers the chordal deviation of the tessellation used by the
# winding test, roughly one more covers the distance query against the
# true surface. Documented, not tuned to any test.
BAND_SAFETY = 2.0

# Documented fallback entity tolerance when no shape carries a measurable
# face/edge/vertex tolerance: OCCT's nominal contact tolerance.
CONTACT_TOL = 1e-7


def has_solid(shape) -> bool:
    return (not shape.IsNull()) and TopExp_Explorer(shape, TopAbs_SOLID).More()


def combined_domain(*shapes, pad_frac: float = DOMAIN_PAD_FRAC):
    """Sample domain from the combined bbox of the given shapes.

    Returns (lo, hi) float arrays: the smallest axis-aligned box covering
    every non-void input bbox, padded by pad_frac * (largest box edge) on
    every side. Null shapes and void bboxes are skipped; if every input is
    void the documented FALLBACK_DOMAIN cube is returned.
    """
    box = Bnd_Box()
    for s in shapes:
        if s is None or s.IsNull():
            continue
        BRepBndLib.Add_s(s, box)
    if box.IsVoid():
        return (np.array(FALLBACK_DOMAIN[0], dtype=float),
                np.array(FALLBACK_DOMAIN[1], dtype=float))
    lo = np.array(box.CornerMin().Coord(), dtype=float)
    hi = np.array(box.CornerMax().Coord(), dtype=float)
    pad = float(pad_frac) * float(np.max(hi - lo))
    return lo - pad, hi + pad


def _max_entity_tolerance(shape) -> float:
    """Largest BRep_Tool tolerance over the shape's faces, edges, vertices."""
    tmax = 0.0
    for ttype, cast in ((TopAbs_FACE, TopoDS.Face),
                        (TopAbs_EDGE, TopoDS.Edge),
                        (TopAbs_VERTEX, TopoDS.Vertex)):
        ex = TopExp_Explorer(shape, ttype)
        while ex.More():
            tmax = max(tmax, float(BRep_Tool.Tolerance_s(cast(ex.Current()))))
            ex.Next()
    return tmax


def exclusion_band(*shapes, deflection: float = 2e-4) -> float:
    """Scale-aware surface-exclusion band for the winding audit.

    band = BAND_SAFETY * deflection + max entity tolerance, where the max
    runs over the face/edge/vertex tolerances (BRep_Tool.Tolerance_s) of
    every non-null shape passed. The tessellation behind the winding test
    can deviate from the true surface by up to the mesh deflection, and
    the B-rep itself is only accurate to its entity tolerances, so points
    closer than this band to an input surface carry no verdict. If no
    shape yields a measurable tolerance the documented fallback is
    BAND_SAFETY * deflection + CONTACT_TOL.
    """
    tmax = 0.0
    for s in shapes:
        if s is None or s.IsNull():
            continue
        tmax = max(tmax, _max_entity_tolerance(s))
    if tmax <= 0.0:
        tmax = CONTACT_TOL
    return BAND_SAFETY * float(deflection) + tmax


# G13: per-audit tessellation cache.  Keyed by (id(shape), deflection)
# with a strong reference to the shape so id() reuse after GC cannot
# alias a dead entry.  Cleared at the start of every membership_audit;
# entries never outlive one audit.
_TRIANGLE_CACHE: dict = {}


def clear_triangle_cache() -> None:
    """Drop cached tessellations.  Called once per membership audit."""
    _TRIANGLE_CACHE.clear()


def triangles(shape, deflection: float = 2e-4) -> np.ndarray:
    """Oriented triangle soup (n, 3, 3) of a shape's faces.

    The input shape is never modified: it is deep-copied
    (BRepBuilderAPI_Copy, which duplicates the TShapes, so tessellations
    attach to the copy only) and the copy is cleaned and tessellated.

    G13: tessellations are cached per audit (see clear_triangle_cache,
    called at the start of membership_audit) so a shape appearing more
    than once in one audit is tessellated once.
    """
    key = (id(shape), float(deflection))
    hit = _TRIANGLE_CACHE.get(key)
    if hit is not None and hit[0] is shape:
        return hit[1]
    work = BRepBuilderAPI_Copy(shape).Shape()
    BRepTools.Clean_s(work)
    BRepMesh_IncrementalMesh(work, deflection, False, 0.1, True)
    out = []
    ex = TopExp_Explorer(work, TopAbs_FACE)
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
    tris = np.asarray(out, dtype=np.float64).reshape(-1, 3, 3)
    _TRIANGLE_CACHE[key] = (shape, tris)
    return tris


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


def membership_audit(a, b, out, op, rng, n=300, lo=None, hi=None,
                     band=None, deflection=2e-4) -> dict:
    """Return counts of kernel errors and OCCT-classifier disagreements.

    lo/hi default to None, meaning the sample domain is derived from the
    combined bounding box of a, b and out via combined_domain() (with the
    documented margin); pass explicit bounds for the legacy behavior.
    band defaults to None, meaning the surface-exclusion band is derived
    via exclusion_band(a, b, deflection=deflection); pass an explicit
    value for the legacy absolute behavior.

    G13: the per-audit tessellation cache is cleared on entry, so each
    distinct shape in this audit is tessellated at most once.
    """
    clear_triangle_cache()
    if lo is None or hi is None:
        dlo, dhi = combined_domain(a, b, out)
        if lo is None:
            lo = dlo
        if hi is None:
            hi = dhi
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    if band is None:
        band = exclusion_band(a, b, deflection=deflection)
    ta, tb = triangles(a, deflection), triangles(b, deflection)
    to = triangles(out, deflection) if has_solid(out) else np.zeros((0, 3, 3))
    res = {"checked": 0, "kernel_errors": 0, "classifier_disagreements": 0,
           "error_points": [], "classifier_points": [],
           "domain": [lo.tolist(), hi.tolist()],
           "band_used": float(band),
           "skipped_winding": 0, "skipped_band": 0, "interior_checked": 0}
    for q in rng.uniform(lo, hi, size=(n, 3)):
        wa, wb, wo = winding(ta, q), winding(tb, q), winding(to, q)
        # Skip points whose winding number is not decisively 0 or 1
        # (inside the tessellation band) and points near any input surface.
        if any(0.05 < abs(w) < 0.95 for w in (wa, wb, wo)):
            res["skipped_winding"] += 1
            continue
        if surface_distance(a, q) < band or surface_distance(b, q) < band:
            res["skipped_band"] += 1
            continue
        ia, ib, io = abs(wa) > 0.5, abs(wb) > 0.5, abs(wo) > 0.5
        res["checked"] += 1
        if ia or ib or io:
            res["interior_checked"] += 1
        if io != want(op, ia, ib):
            res["kernel_errors"] += 1
            res["error_points"].append(q.tolist())
        if has_solid(out):
            st = occt_state(out, q)
            if st in (TopAbs_IN, TopAbs_OUT) and (st == TopAbs_IN) != io:
                res["classifier_disagreements"] += 1
                res["classifier_points"].append(q.tolist())
    return res
