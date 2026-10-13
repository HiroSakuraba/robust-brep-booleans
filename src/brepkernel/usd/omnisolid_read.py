"""Read omniSolid BrepArray prims (.usda text) into OCCT solids.

This is the reverse of omnisolid_write. It needs neither pxr nor usd-brep:
the arrays are parsed from the .usda text that usd-brep's
``usd.export_brep`` writes (or that omnisolid_write produced).

What is rebuilt exactly:
  * every 3D edge curve (NURBS, trimmed to edge:range) and vertex;
  * every face surface (plane, cylinder, cone, sphere, torus, NURBS);
  * loops, edge sharing, shells, cavities and region membership.

What is computed rather than read: OCCT pcurves (UV curves). They are
built from the exact 3D edges by _face_pcurves (exact on planes; sampled,
unwrapped, fitted per knot span and verified in 3D elsewhere), so the UV
conventions of the source surfaces never enter the OCCT model. Authored
UV curves, when present, are ignored.

Every rebuilt solid is checked with BRepCheck_Analyzer, must have a
positive volume, and must keep its tolerances within 10x the file
tolerance; anything else raises OmniSolidReadError rather than returning
a model that only looks right.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass

from .. import _occt_compat  # noqa: F401  (TopoDS cast spellings on OCP 7.8)

UNUSED_VERTEX = 9999999


class OmniSolidReadError(ValueError):
    """The BrepArray could not be rebuilt as a valid OCCT solid.

    When a shape was assembled but rejected, it is kept on ``.shape`` for
    diagnosis; it must not be used as a result.
    """

    def __init__(self, message: str, shape=None):
        super().__init__(message)
        self.shape = shape


# ---------------------------------------------------------------------------
# .usda text parsing
# ---------------------------------------------------------------------------

_DEF = re.compile(r'def\s+(\w+)\s+"([^"]+)"')
_ATTR = re.compile(r'uniform\s+(\S+)\s+(\S+)\s*=\s*\[(.*?)\]', re.S)


def _values(typ: str, body: str) -> list:
    body = body.strip()
    if not body:
        return []
    if typ.startswith("token"):
        return re.findall(r'"([^"]*)"', body)
    if "(" in body:
        conv = int if typ.startswith(("int", "uint")) else float
        return [tuple(conv(float(x)) if conv is int else float(x) for x in g.split(","))
                for g in re.findall(r"\(([^)]*)\)", body)]
    parts = [p.strip() for p in body.split(",") if p.strip()]
    if typ.startswith(("double", "float")):
        return [float(p) for p in parts]
    return [int(float(p)) for p in parts]


def _block_end(text: str, start: int) -> int:
    """Index just past the brace block whose '{' is the first at or after start."""
    i = text.index("{", start)
    depth = 0
    while i < len(text):
        c = text[i]
        if c == '"':
            i = text.index('"', i + 1)
        elif c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return i + 1
        i += 1
    raise OmniSolidReadError("unbalanced braces in .usda")


def brep_array_prims(path: str) -> dict:
    """{prim path: {attribute name: values}} for every BrepArray in the file.

    Handles the flat layout usd-brep and omnisolid_write produce
    (prims nested in Xforms, one level of metadata parentheses). Binary
    .usdc is not supported; export with a .usda filename.
    """
    with open(path, "rb") as fh:
        head = fh.read(8)
    if head.startswith(b"PXR-USDC"):
        raise OmniSolidReadError(f"{path} is binary usdc; export to .usda instead")
    text = open(path, encoding="utf-8").read()
    out = {}

    def walk(lo: int, hi: int, prefix: str) -> None:
        pos = lo
        while True:
            m = _DEF.search(text, pos, hi)
            if not m:
                return
            kind, name = m.group(1), m.group(2)
            # skip the optional "( ... )" metadata block before the body
            j = m.end()
            while text[j].isspace():
                j += 1
            if text[j] == "(":
                depth = 0
                while True:
                    if text[j] == "(":
                        depth += 1
                    elif text[j] == ")":
                        depth -= 1
                        if depth == 0:
                            break
                    elif text[j] == '"':
                        j = text.index('"', j + 1)
                    j += 1
            body_start = text.index("{", j)
            end = _block_end(text, body_start)
            prim = f"{prefix}/{name}"
            if kind == "BrepArray":
                body = text[body_start:end]
                out[prim] = {n: _values(t, b) for t, n, b in _ATTR.findall(body)}
            else:
                walk(body_start + 1, end - 1, prim)
            pos = end

    walk(0, len(text), "")
    return out


# ---------------------------------------------------------------------------
# OCCT geometry builders
# ---------------------------------------------------------------------------

def _arrays():
    try:
        from OCP.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
        from OCP.TColStd import (TColStd_Array1OfInteger, TColStd_Array1OfReal,
                                 TColStd_Array2OfReal)
        return (TColgp_Array1OfPnt, TColgp_Array2OfPnt, TColStd_Array1OfReal,
                TColStd_Array2OfReal, TColStd_Array1OfInteger)
    except ImportError:                       # OCP 8 moved collections
        from OCP.collections import (Array1_double, Array1_gp_Pnt, Array1_int,
                                     Array2_double, Array2_gp_Pnt)
        return Array1_gp_Pnt, Array2_gp_Pnt, Array1_double, Array2_double, Array1_int


def _array1_pnt2d():
    try:
        from OCP.TColgp import TColgp_Array1OfPnt2d
        return TColgp_Array1OfPnt2d
    except ImportError:
        from OCP.collections import Array1_gp_Pnt2d
        return Array1_gp_Pnt2d


def _knots(full: list):
    """Full knot vector -> (distinct knots, multiplicities)."""
    ks, ms = [], []
    for k in full:
        if ks and abs(k - ks[-1]) <= 1e-12 * max(1.0, abs(k)):
            ms[-1] += 1
        else:
            ks.append(k)
            ms.append(1)
    return ks, ms


def _nurbs_curve(cvs, weights, knots, order):
    from OCP.Geom import Geom_BSplineCurve
    from OCP.gp import gp_Pnt
    A1P, _A2P, A1R, _A2R, A1I = _arrays()
    n = len(cvs)
    if len(knots) != n + order:
        raise OmniSolidReadError(f"curve knot count {len(knots)} != {n} + {order}")
    ks, ms = _knots(knots)
    # schema knots list every knot once per multiplicity, with the end
    # knots repeated `order` times; OCCT wants degree + 1 at the ends
    poles, w = A1P(1, n), A1R(1, n)
    for i, (p, wi) in enumerate(zip(cvs, weights), 1):
        poles.SetValue(i, gp_Pnt(*p))
        w.SetValue(i, float(wi))
    kk, mm = A1R(1, len(ks)), A1I(1, len(ms))
    for i, (k, m) in enumerate(zip(ks, ms), 1):
        kk.SetValue(i, float(k))
        mm.SetValue(i, int(m))
    rational = any(abs(x - 1.0) > 1e-15 for x in weights)
    if rational:
        return Geom_BSplineCurve(poles, w, kk, mm, order - 1)
    return Geom_BSplineCurve(poles, kk, mm, order - 1)


def _nurbs_surface(cvs, weights, uk, vk, uo, vo, nu, nv):
    from OCP.Geom import Geom_BSplineSurface
    from OCP.gp import gp_Pnt
    _A1P, A2P, A1R, A2R, A1I = _arrays()
    if len(cvs) != nu * nv or len(uk) != nu + uo or len(vk) != nv + vo:
        raise OmniSolidReadError("inconsistent NURBS surface array sizes")
    poles, w = A2P(1, nu, 1, nv), A2R(1, nu, 1, nv)
    for i in range(nu):
        for j in range(nv):
            poles.SetValue(i + 1, j + 1, gp_Pnt(*cvs[i * nv + j]))   # U rows
            w.SetValue(i + 1, j + 1, float(weights[i * nv + j]))
    def kv(full):
        ks, ms = _knots(full)
        kk, mm = A1R(1, len(ks)), A1I(1, len(ms))
        for i, (k, m) in enumerate(zip(ks, ms), 1):
            kk.SetValue(i, float(k))
            mm.SetValue(i, int(m))
        return kk, mm
    ukk, umm = kv(uk)
    vkk, vmm = kv(vk)
    if any(abs(x - 1.0) > 1e-15 for x in weights):
        return Geom_BSplineSurface(poles, w, ukk, vkk, umm, vmm, uo - 1, vo - 1)
    return Geom_BSplineSurface(poles, ukk, vkk, umm, vmm, uo - 1, vo - 1)


def _ax3(origin, axis, ref):
    from OCP.gp import gp_Ax3, gp_Dir, gp_Pnt
    return gp_Ax3(gp_Pnt(*origin), gp_Dir(*axis), gp_Dir(*ref))


@dataclass
class _Cursor:
    """Running index into one packed per-type array family."""
    i: int = 0

    def take(self) -> int:
        k = self.i
        self.i += 1
        return k


def _surfaces(A: dict) -> list:
    """One Geom_Surface per face, in face order."""
    from OCP.Geom import (Geom_ConicalSurface, Geom_CylindricalSurface, Geom_Plane,
                          Geom_SphericalSurface, Geom_ToroidalSurface)
    cur = {k: _Cursor() for k in ("plane", "cone", "cylinder", "sphere", "torus", "nurb")}
    nurb_cv = nurb_uk = nurb_vk = 0
    out = []
    g = lambda key: A.get(key, [])
    for st in A["face:surfaceType"]:
        if st == "BrepSurfacePlaneAPI":
            k = cur["plane"].take()
            out.append(Geom_Plane(_ax3(g("brep:surface:plane:origin")[k],
                                       g("brep:surface:plane:axis")[k],
                                       g("brep:surface:plane:refDirection")[k])))
        elif st == "BrepSurfaceConeAPI":
            k = cur["cone"].take()
            ax = _ax3(g("brep:surface:cone:origin")[k], g("brep:surface:cone:axis")[k],
                      g("brep:surface:cone:refDirection")[k])
            r, a = g("brep:surface:cone:radius")[k], g("brep:surface:cone:semiAngle")[k]
            if abs(a) < 1e-12:
                out.append(Geom_CylindricalSurface(ax, r))
            else:
                out.append(Geom_ConicalSurface(ax, a, r))
        elif st == "BrepSurfaceCylinderAPI":
            k = cur["cylinder"].take()
            out.append(Geom_CylindricalSurface(
                _ax3(g("brep:surface:cylinder:origin")[k], g("brep:surface:cylinder:axis")[k],
                     g("brep:surface:cylinder:refDirection")[k]),
                g("brep:surface:cylinder:radius")[k]))
        elif st == "BrepSurfaceSphereAPI":
            k = cur["sphere"].take()
            out.append(Geom_SphericalSurface(
                _ax3(g("brep:surface:sphere:center")[k], g("brep:surface:sphere:axis")[k],
                     g("brep:surface:sphere:refDirection")[k]),
                g("brep:surface:sphere:radius")[k]))
        elif st == "BrepSurfaceTorusAPI":
            k = cur["torus"].take()
            out.append(Geom_ToroidalSurface(
                _ax3(g("brep:surface:torus:origin")[k], g("brep:surface:torus:axis")[k],
                     g("brep:surface:torus:refDirection")[k]),
                g("brep:surface:torus:majorRadius")[k], g("brep:surface:torus:minorRadius")[k]))
        elif st == "BrepSurfaceNurbAPI":
            k = cur["nurb"].take()
            nu = g("brep:surface:nurb:uVertexCount")[k]
            nv = g("brep:surface:nurb:vVertexCount")[k]
            uo = g("brep:surface:nurb:uOrder")[k]
            vo = g("brep:surface:nurb:vOrder")[k]
            n = nu * nv
            out.append(_nurbs_surface(
                g("brep:surface:nurb:controlVertices")[nurb_cv:nurb_cv + n],
                g("brep:surface:nurb:weights")[nurb_cv:nurb_cv + n],
                g("brep:surface:nurb:uKnots")[nurb_uk:nurb_uk + nu + uo],
                g("brep:surface:nurb:vKnots")[nurb_vk:nurb_vk + nv + vo],
                uo, vo, nu, nv))
            nurb_cv += n
            nurb_uk += nu + uo
            nurb_vk += nv + vo
        else:
            raise OmniSolidReadError(f"unsupported surface type {st}")
    return out


def _curves(A: dict) -> list:
    """(Geom_BSplineCurve, t0, t1) per edge, in edge order."""
    out = []
    cv = kn = 0
    k = 0
    rng = A["edge:range"]
    for ei, ct in enumerate(A["edge:curveType"]):
        if ct != "BrepCurve3dNurbAPI":
            raise OmniSolidReadError(f"unsupported edge curve type {ct}")
        n = A["brep:edge3dNurb:curve3d:nurb:vertexCount"][k]
        o = A["brep:edge3dNurb:curve3d:nurb:order"][k]
        c = _nurbs_curve(A["brep:edge3dNurb:curve3d:nurb:controlVertices"][cv:cv + n],
                         A["brep:edge3dNurb:curve3d:nurb:weights"][cv:cv + n],
                         A["brep:edge3dNurb:curve3d:nurb:knots"][kn:kn + n + o], o)
        cv += n
        kn += n + o
        k += 1
        out.append((c, float(rng[2 * ei]), float(rng[2 * ei + 1])))
    return out


# ---------------------------------------------------------------------------
# topology
# ---------------------------------------------------------------------------

class _UseUV:
    """Unwrapped UV samples of one edge use, in traversal order."""

    def __init__(self, ei, same, ts, uv, sing):
        self.ei, self.same = ei, same
        self.ts, self.uv, self.sing = ts, uv, sing      # lists, traversal order
        self.placed = False

    def start(self):
        return self.uv[0], self.sing[0]

    def end(self):
        return self.uv[-1], self.sing[-1]

    def shift(self, du, dv):
        self.uv = [(u + du, v + dv) for u, v in self.uv]


def _face_pcurves(B, face, surface, edges, loops, tol, fi, size=1.0):
    """Exact-enough UV curves for every edge use of `face`.

    OCCT's projection (ShapeFix) treats each edge on its own, so an edge
    that starts on a periodic seam can land on the wrong copy of the
    period and cut back across the face; the same-parameter step then
    hides the error by inflating tolerances (seen: 0.068). Instead:

    * planes: the 3D B-spline's poles are mapped affinely into the plane
      frame, which is exact;
    * other surfaces: the exact 3D edge is sampled, each sample is
      inverted onto the surface, periodic parameters are unwrapped
      continuously, samples at singular points (poles, apexes) take their
      singular parameter from their neighbours, and a 2D B-spline is fitted
      per knot span of the edge, on the edge's own parameter, refined
      towards 1e-10 of the part size and never accepted above tol / 4;
    * seams (an edge used twice by the face) get two curves one period
      apart, each use on the side its loop direction requires (the face
      lies on the left of every loop in UV);
    * every other use is shifted by whole periods to meet its neighbours
      at shared non-singular vertices, starting from the seams.
    A use whose fit cannot be verified raises OmniSolidReadError.
    """
    from OCP.BRep import BRep_Tool
    from OCP.Geom import Geom_Plane
    from OCP.ShapeAnalysis import ShapeAnalysis_Surface
    from OCP.gp import gp_Pnt, gp_Vec

    if isinstance(surface, Geom_Plane):
        _plane_pcurves(B, face, surface, edges, loops, tol)
        return

    sas = ShapeAnalysis_Surface(surface)
    u0, u1, v0, v1 = surface.Bounds()
    per_u = surface.IsUPeriodic() or surface.IsUClosed()
    per_v = surface.IsVPeriodic() or surface.IsVClosed()
    pu = surface.UPeriod() if surface.IsUPeriodic() else (u1 - u0)
    pv = surface.VPeriod() if surface.IsVPeriodic() else (v1 - v0)

    def singular(u, v):
        p, du, dv = gp_Pnt(), gp_Vec(), gp_Vec()
        surface.D1(u, v, p, du, dv)
        scale = max(1.0, abs(p.X()), abs(p.Y()), abs(p.Z()))
        return du.Magnitude() < 1e-9 * scale, dv.Magnitude() < 1e-9 * scale

    def near(x, ref, period):
        return x + period * round((ref - x) / period)

    counts: dict = {}
    for lp in loops:
        for ei, same in lp:
            counts.setdefault(ei, []).append(same)
    for ei, ss in counts.items():
        if len(ss) > 2 or (len(ss) == 2 and ss[0] == ss[1]):
            raise OmniSolidReadError(f"face {fi}: edge {ei} used {len(ss)} times "
                                     f"with orientations {ss}")

    def breaks(ei):
        """Knots of the 3D curve inside the edge range: the UV curve may
        only lose smoothness there."""
        e = edges[ei]
        c3, (a, b) = BRep_Tool.Curve_s(e, 0.0, 0.0), BRep_Tool.Range_s(e)
        ks = [a]
        if hasattr(c3, "NbKnots"):
            for i in range(1, c3.NbKnots() + 1):
                k = c3.Knot(i)
                if a + 1e-12 * (b - a) < k < b - 1e-12 * (b - a):
                    ks.append(k)
        ks.append(b)
        return ks

    def sample(ei, same, m):
        """m intervals per knot span, ascending, then put in traversal order."""
        e = edges[ei]
        c3 = BRep_Tool.Curve_s(e, 0.0, 0.0)
        ks = breaks(ei)
        ts = [ks[0]]
        for lo, hi in zip(ks[:-1], ks[1:]):
            ts += [lo + (hi - lo) * k / m for k in range(1, m + 1)]
        if not same:
            ts = ts[::-1]
        uv, sing = [], []
        prev = None
        for t in ts:
            q = sas.ValueOfUV(c3.Value(t), tol)
            u, v = q.X(), q.Y()
            if prev is not None:
                if per_u:
                    u = near(u, prev[0], pu)
                if per_v:
                    v = near(v, prev[1], pv)
            su, sv = singular(u, v)
            uv.append([u, v])
            sing.append((su, sv))
            if not su and not sv:
                prev = (u, v)
        # singular samples take the singular coordinate from a neighbour
        for axis in (0, 1):
            good = [k for k, sg in enumerate(sing) if not sg[axis]]
            if not good:
                continue
            for k, sg in enumerate(sing):
                if sg[axis]:
                    nb = min(good, key=lambda g: abs(g - k))
                    uv[k][axis] = uv[nb][axis]
        return _UseUV(ei, same, ts, [tuple(x) for x in uv],
                      [bool(sg[0] or sg[1]) for sg in sing])

    def fit(use):
        """Piecewise 2D B-spline: one fit per knot span of the 3D curve,
        joined C0 at those knots, on the edge's own parameter."""
        ts, uv = use.ts, use.uv
        if ts[0] > ts[-1]:
            ts, uv = ts[::-1], uv[::-1]
        ks = breaks(use.ei)
        pieces = []
        for lo, hi in zip(ks[:-1], ks[1:]):
            idx = [k for k, t in enumerate(ts) if lo - 1e-15 <= t <= hi + 1e-15]
            if len(idx) < 2:
                return None, float("inf")
            pieces.append(_fit_piece([ts[k] for k in idx], [uv[k] for k in idx]))
            if pieces[-1] is None:
                return None, float("inf")
        c2 = _join_c0(pieces)
        c3 = BRep_Tool.Curve_s(edges[use.ei], 0.0, 0.0)
        dev = 0.0
        for k in range(len(ts) - 1):
            for f in (0.25, 0.5, 0.75):
                t = ts[k] + f * (ts[k + 1] - ts[k])
                q = c2.Value(t)
                dev = max(dev, surface.Value(q.X(), q.Y()).Distance(c3.Value(t)))
        return c2, dev

    uses = [[(ei, same, 16) for ei, same in lp] for lp in loops]

    built = [[sample(ei, same, n) for ei, same, n in row] for row in uses]

    # Seams fix the parameter window of the face. A u-seam is an iso-u
    # line used twice; the face lies left of every loop in UV, so the use
    # travelling +v is the high copy and the other sits one period lower.
    # The face then occupies exactly [low copy, low copy + period] in u.
    window = {0: None, 1: None}
    seam_axis = {}
    pairs: dict = {}
    for row in built:
        for use in row:
            if len(counts[use.ei]) == 2:
                pairs.setdefault(use.ei, []).append(use)
    for ei, pair in pairs.items():
        a = pair[0]
        du = a.uv[-1][0] - a.uv[0][0]
        dv = a.uv[-1][1] - a.uv[0][1]
        if per_u and abs(du) <= abs(dv):
            axis, period, along = 0, pu, 1
        elif per_v:
            axis, period, along = 1, pv, 0
        else:
            raise OmniSolidReadError(f"face {fi}: edge {ei} is used twice but the "
                                     f"surface is not closed")
        seam_axis[ei] = axis
        c = sum(p[axis] for p in a.uv) / len(a.uv)
        a_forward = (a.uv[-1][along] > a.uv[0][along])
        # u-seam: +v traveller is high; v-seam: +u traveller is low
        a_high = a_forward if axis == 0 else not a_forward
        lo = c - period if a_high else c
        if not surface.IsUPeriodic() and axis == 0:
            lo = u0
        if not surface.IsVPeriodic() and axis == 1:
            lo = v0
        if window[axis] is not None:
            lo = window[axis]       # every piece of a seam shares one window
        window[axis] = lo
        for use in pair:
            fwd = use.uv[-1][along] > use.uv[0][along]
            high = fwd if axis == 0 else not fwd
            target = lo + (period if high else 0.0)
            use.uv = [((target, q[1]) if axis == 0 else (q[0], target)) for q in use.uv]

    def into_window(use, axis, period):
        lo = window[axis]
        vals = [q[axis] for q in use.uv]
        mid = 0.5 * (min(vals) + max(vals))
        k = math.floor((mid - lo) / period)
        if k:
            use.shift(-k * period if axis == 0 else 0.0, -k * period if axis == 1 else 0.0)

    def mid(use, axis):
        vals = [q[axis] for q in use.uv]
        return 0.5 * (min(vals) + max(vals))

    for li, row in enumerate(built):
        prev = None
        for use in row:
            for axis, per, period in ((0, per_u, pu), (1, per_v, pv)):
                if not per:
                    continue
                if seam_axis.get(use.ei) == axis:
                    continue                      # constant coordinate already set
                if window[axis] is not None:
                    into_window(use, axis, period)
                elif prev is not None:
                    (pe, psing) = prev.end()
                    if not psing:
                        # continuity at a regular vertex
                        x0 = use.uv[0][axis]
                        d = near(x0, pe[axis], period) - x0
                    else:
                        # through a pole the coordinate jumps; take the copy
                        # nearest the previous use (a face spans < 1 period
                        # in a direction without a seam)
                        m = mid(use, axis)
                        d = near(m, mid(prev, axis), period) - m
                    use.shift(d if axis == 0 else 0.0, d if axis == 1 else 0.0)
            prev = use

    # Repair: a use lying exactly on a window boundary (for example an edge
    # running along the seam line that this face uses only once) has two
    # admissible copies; pick the one that connects to its neighbours.
    eps = 1e-9 * max(1.0, abs(pu) if per_u else 1.0, abs(pv) if per_v else 1.0)
    for row in built:
        for _sweep in range(2):
            for k, use in enumerate(row):
                for axis, per, period in ((0, per_u, pu), (1, per_v, pv)):
                    if not per or window[axis] is None or seam_axis.get(use.ei) == axis:
                        continue
                    lo, hi = window[axis] - 1e-7, window[axis] + period + 1e-7
                    prev_use, next_use = row[k - 1], row[(k + 1) % len(row)]
                    for ref, mine in ((prev_use.end(), use.start()),
                                      (next_use.start(), use.end())):
                        (rp, rsing), (mp, msing) = ref, mine
                        if rsing or msing:
                            continue
                        gap = rp[axis] - mp[axis]
                        if abs(gap) < eps:
                            break                         # already connected
                        d = period * round(gap / period)
                        if d and abs(gap - d) < 1e-6 * period:
                            vals = [q[axis] + d for q in use.uv]
                            if lo <= min(vals) and max(vals) <= hi:
                                use.shift(d if axis == 0 else 0.0, d if axis == 1 else 0.0)
                                break

    # fit, verify, attach
    fitted: dict = {}
    for row, spec in zip(built, uses):
        for use, (_ei, _same, n0) in zip(row, spec):
            c2, dev = fit(use)
            n = n0
            goal = max(1e-12, 1e-10 * size)       # aim: no visible volume change
            while dev > goal and n < 256:
                n *= 2
                shift_uv = use.uv[0]
                fresh = sample(use.ei, use.same, n)
                du = (near(fresh.uv[0][0], shift_uv[0], pu) - fresh.uv[0][0]) if per_u else 0.0
                dv = (near(fresh.uv[0][1], shift_uv[1], pv) - fresh.uv[0][1]) if per_v else 0.0
                fresh.shift(du, dv)
                if len(counts[use.ei]) == 2:     # seams keep their exact constant
                    if per_u and abs(use.uv[0][0] - use.uv[-1][0]) < 1e-12:
                        fresh.uv = [(use.uv[0][0], p[1]) for p in fresh.uv]
                    else:
                        fresh.uv = [(p[0], use.uv[0][1]) for p in fresh.uv]
                use.ts, use.uv, use.sing = fresh.ts, fresh.uv, fresh.sing
                c2, dev = fit(use)
            if c2 is None or dev > 0.25 * tol:
                raise OmniSolidReadError(
                    f"face {fi}: UV curve for edge {use.ei} deviates {dev:.3g} "
                    f"from the 3D edge (limit {0.25 * tol:.3g})")
            fitted.setdefault(use.ei, {})[use.same] = c2
    for ei, by in fitted.items():
        e = edges[ei]
        if len(by) == 2:
            B.UpdateEdge(e, by[True], by[False], face, BRep_Tool.Tolerance_s(e))
        else:
            B.UpdateEdge(e, next(iter(by.values())), face, BRep_Tool.Tolerance_s(e))


def _fit_piece(ts, uv):
    """2D B-spline through (ts, uv) on parameter [ts[0], ts[-1]], ends exact."""
    from OCP.Geom2dAPI import Geom2dAPI_PointsToBSpline
    from OCP.GeomAbs import GeomAbs_C2
    from OCP.gp import gp_Pnt2d
    A2 = _array1_pnt2d()
    _A1P, _A2P, A1R, _A2R, _A1I = _arrays()
    n = len(ts)
    pts, prm = A2(1, n), A1R(1, n)
    for k, (t, (u, v)) in enumerate(zip(ts, uv), 1):
        pts.SetValue(k, gp_Pnt2d(u, v))
        prm.SetValue(k, t)
    ap = Geom2dAPI_PointsToBSpline(pts, prm, 1 if n < 4 else 3, min(8, n - 1), GeomAbs_C2, 1e-10)
    if not ap.IsDone():
        return None
    c = _remap_knots(ap.Curve(), ts[0], ts[-1])
    c.SetPole(1, gp_Pnt2d(*uv[0]))
    c.SetPole(c.NbPoles(), gp_Pnt2d(*uv[-1]))
    return c


def _join_c0(pieces):
    """Concatenate consecutive non-rational 2D B-splines into one, C0 at joins."""
    from OCP.Geom2d import Geom2d_BSplineCurve
    if len(pieces) == 1:
        return pieces[0]
    A2 = _array1_pnt2d()
    _A1P, _A2P, A1R, _A2R, A1I = _arrays()
    d = max(p.Degree() for p in pieces)
    for p in pieces:
        if p.Degree() < d:
            p.IncreaseDegree(d)
    poles, knots, mults = [], [], []
    for i, p in enumerate(pieces):
        pp = [p.Pole(k) for k in range(1, p.NbPoles() + 1)]
        kk = [p.Knot(k) for k in range(1, p.NbKnots() + 1)]
        mm = [p.Multiplicity(k) for k in range(1, p.NbKnots() + 1)]
        if i == 0:
            poles += pp
            knots += kk
            mults += mm
        else:
            poles += pp[1:]               # shared end point
            mults[-1] = d                 # C0 join
            knots += kk[1:]
            mults += mm[1:]
    mults[-1] = d + 1
    P, K, M = A2(1, len(poles)), A1R(1, len(knots)), A1I(1, len(mults))
    for i, q in enumerate(poles, 1):
        P.SetValue(i, q)
    for i, (k, m) in enumerate(zip(knots, mults), 1):
        K.SetValue(i, k)
        M.SetValue(i, m)
    return Geom2d_BSplineCurve(P, K, M, d)


def _remap_knots(c, t0, t1):
    """The same 2D B-spline with its knots mapped affinely onto [t0, t1].

    Geom2dAPI_PointsToBSpline normalizes the given parameters to [0, 1];
    the UV curve must share the 3D edge's parameter.
    """
    from OCP.Geom2d import Geom2d_BSplineCurve
    A2 = _array1_pnt2d()
    _A1P, _A2P, A1R, _A2R, A1I = _arrays()
    a, b = c.FirstParameter(), c.LastParameter()
    n, m = c.NbPoles(), c.NbKnots()
    poles, w = A2(1, n), A1R(1, n)
    for i in range(1, n + 1):
        poles.SetValue(i, c.Pole(i))
        w.SetValue(i, c.Weight(i))
    kk, mm = A1R(1, m), A1I(1, m)
    for i in range(1, m + 1):
        kk.SetValue(i, t0 + (c.Knot(i) - a) * (t1 - t0) / (b - a))
        mm.SetValue(i, c.Multiplicity(i))
    if c.IsRational():
        return Geom2d_BSplineCurve(poles, w, kk, mm, c.Degree())
    return Geom2d_BSplineCurve(poles, kk, mm, c.Degree())


def _plane_pcurves(B, face, plane, edges, loops, tol):
    """Exact UV curves on a plane: map the 3D B-spline poles affinely."""
    from OCP.BRep import BRep_Tool
    from OCP.Geom2d import Geom2d_BSplineCurve
    from OCP.gp import gp_Pnt2d
    A2 = _array1_pnt2d()
    _A1P, _A2P, A1R, _A2R, A1I = _arrays()
    ax = plane.Position()
    o, x, y = ax.Location(), ax.XDirection(), ax.YDirection()
    done = set()
    for lp in loops:
        for ei, _same in lp:
            if ei in done:
                continue
            done.add(ei)
            c = BRep_Tool.Curve_s(edges[ei], 0.0, 0.0)
            n = c.NbPoles()
            poles = A2(1, n)
            for i in range(1, n + 1):
                p = c.Pole(i)
                d = (p.X() - o.X(), p.Y() - o.Y(), p.Z() - o.Z())
                poles.SetValue(i, gp_Pnt2d(d[0] * x.X() + d[1] * x.Y() + d[2] * x.Z(),
                                           d[0] * y.X() + d[1] * y.Y() + d[2] * y.Z()))
            kk, mm = A1R(1, c.NbKnots()), A1I(1, c.NbKnots())
            for i in range(1, c.NbKnots() + 1):
                kk.SetValue(i, c.Knot(i))
                mm.SetValue(i, c.Multiplicity(i))
            if c.IsRational():
                w = A1R(1, n)
                for i in range(1, n + 1):
                    w.SetValue(i, c.Weight(i))
                c2 = Geom2d_BSplineCurve(poles, w, kk, mm, c.Degree())
            else:
                c2 = Geom2d_BSplineCurve(poles, kk, mm, c.Degree())
            B.UpdateEdge(edges[ei], c2, face, BRep_Tool.Tolerance_s(edges[ei]))


def _build(A: dict) -> list:
    """List of OCCT shapes (one per Brep in the prim)."""
    from OCP.BRep import BRep_Builder, BRep_Tool
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.BRepLib import BRepLib
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    from OCP.ShapeFix import ShapeFix_Wire
    from OCP.TopAbs import TopAbs_FORWARD, TopAbs_REVERSED
    from OCP.TopoDS import TopoDS_Compound, TopoDS_Face, TopoDS_Shell, \
        TopoDS_Solid, TopoDS_Vertex, TopoDS_Wire
    from OCP.gp import gp_Pnt

    need = ("brep:regionCount", "region:shellCount", "region:type", "shell:faceuseCount",
            "faceuse:faceIndex", "faceuse:orientationType", "face:loopCount",
            "face:surfaceType", "loop:edgeuseCount", "edgeuse:edgeIndex",
            "edgeuse:orientationType", "edge:curveType", "edge:range",
            "edge:vertexIndices", "brep:vertexPoint:point:position")
    missing = [k for k in need if k not in A]
    if missing:
        raise OmniSolidReadError(f"BrepArray lacks {missing}")
    if any(c for c in A.get("shell:wireEdgeCount", [])):
        raise OmniSolidReadError("wire edges are not supported")
    tols = A.get("brep:intersectTol3d") or [2e-5]
    tol = float(max(tols))
    B = BRep_Builder()

    verts = []
    for p in A["brep:vertexPoint:point:position"]:
        v = TopoDS_Vertex()
        B.MakeVertex(v, gp_Pnt(*p), tol)
        verts.append(v)

    edges = []
    for ei, (curve, t0, t1) in enumerate(_curves(A)):
        vs, ve = A["edge:vertexIndices"][ei]
        pts_ = A["brep:vertexPoint:point:position"]
        gap = max(curve.Value(t0).Distance(gp_Pnt(*pts_[vs])),
                  curve.Value(t1).Distance(gp_Pnt(*pts_[ve])))
        if gap > 10.0 * tol:
            raise OmniSolidReadError(
                f"edge {ei} ends {gap:.3g} from its vertex (file tolerance {tol:.3g})")
        for vi in (vs, ve):
            if gap > 0.5 * BRep_Tool.Tolerance_s(verts[vi]):
                B.UpdateVertex(verts[vi], 2.0 * gap)
        mk = BRepBuilderAPI_MakeEdge(curve, verts[vs], verts[ve], t0, t1)
        if not mk.IsDone():
            raise OmniSolidReadError(f"edge {ei}: MakeEdge failed ({int(mk.Error())})")
        e = mk.Edge()
        B.UpdateEdge(e, max(tol, 1e-7))
        edges.append(e)

    pts = A["brep:vertexPoint:point:position"]
    size = max(1.0, max((max(p[i] for p in pts) - min(p[i] for p in pts)) for i in range(3))) \
        if pts else 1.0
    surfaces = _surfaces(A)
    faces = []
    eu = 0
    lp = 0
    for fi, nloops in enumerate(A["face:loopCount"]):
        f = TopoDS_Face()
        B.MakeFace(f, surfaces[fi], tol)
        wires = []
        loops = []
        for _ in range(nloops):
            n = A["loop:edgeuseCount"][lp]
            lp += 1
            if n == 0:
                raise OmniSolidReadError(f"face {fi}: vertex-only loops are not supported")
            w = TopoDS_Wire()
            B.MakeWire(w)
            loop = []
            for _k in range(n):
                ei = A["edgeuse:edgeIndex"][eu]
                same = A["edgeuse:orientationType"][eu] == "same"
                B.Add(w, edges[ei].Oriented(TopAbs_FORWARD if same else TopAbs_REVERSED))
                loop.append((ei, same))
                eu += 1
            loops.append(loop)
            w.Closed(True)
            B.Add(f, w)
            wires.append(w)
        _face_pcurves(B, f, surfaces[fi], edges, loops, tol, fi, size)
        # degenerate edges where a loop passes through a pole
        f2 = TopoDS_Face()
        B.MakeFace(f2, surfaces[fi], tol)
        for w in wires:
            sfw = ShapeFix_Wire(w, f, tol)
            sfw.FixReorderMode = 0
            sfw.FixConnectedMode = 0
            sfw.FixDegenerated()
            B.Add(f2, sfw.Wire())
        faces.append(f2)
    if eu != len(A["edgeuse:edgeIndex"]):
        raise OmniSolidReadError("edgeuse count does not match the loops")

    shapes = []
    r = sh = fu = 0
    for nregions in A["brep:regionCount"]:
        solids = []
        for _ in range(nregions):
            rtype = A["region:type"][r]
            nshells = A["region:shellCount"][r]
            r += 1
            shells = []
            for _s in range(nshells):
                n = A["shell:faceuseCount"][sh]
                sh += 1
                uses = [(A["faceuse:faceIndex"][fu + k], A["faceuse:orientationType"][fu + k])
                        for k in range(n)]
                fu += n
                shells.append(uses)
            if rtype != "solidRegion":
                continue
            so = TopoDS_Solid()
            B.MakeSolid(so)
            for uses in shells:
                shell = TopoDS_Shell()
                B.MakeShell(shell)
                for fi, side in uses:
                    # material on the surface-normal side -> OCCT face reversed
                    ori = TopAbs_REVERSED if side == "same" else TopAbs_FORWARD
                    B.Add(shell, faces[fi].Oriented(ori))
                shell.Closed(True)
                B.Add(so, shell)
            solids.append(so)
        if not solids:
            raise OmniSolidReadError("Brep has no solidRegion")
        if len(solids) == 1:
            shape = solids[0]
        else:
            shape = TopoDS_Compound()
            B.MakeCompound(shape)
            for so in solids:
                B.Add(shape, so)
        BRepLib.SameParameter_s(shape, tol, True)
        grown = _max_tolerance(shape)
        if grown > 10.0 * tol:
            raise OmniSolidReadError(
                f"edge or vertex tolerance grew to {grown:.3g} (file tolerance "
                f"{tol:.3g}); the rebuilt geometry is not faithful", shape)
        chk = BRepCheck_Analyzer(shape)
        if not chk.IsValid():
            raise OmniSolidReadError("rebuilt shape fails BRepCheck_Analyzer", shape)
        g = GProp_GProps()
        BRepGProp.VolumeProperties_s(shape, g, 1e-9, True)
        if not g.Mass() > 0:
            raise OmniSolidReadError(f"rebuilt shape has non-positive volume {g.Mass():.6g}",
                                     shape)
        shapes.append(shape)
    return shapes


def _max_tolerance(shape) -> float:
    from OCP.BRep import BRep_Tool
    from OCP.TopAbs import TopAbs_EDGE, TopAbs_VERTEX
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS
    mx = 0.0
    for kind, cast, get in ((TopAbs_EDGE, TopoDS.Edge, BRep_Tool.Tolerance_s),
                            (TopAbs_VERTEX, TopoDS.Vertex, BRep_Tool.Tolerance_s)):
        ex = TopExp_Explorer(shape, kind)
        while ex.More():
            mx = max(mx, float(get(cast(ex.Current()))))
            ex.Next()
    return mx


def read_omnisolid_usda(path: str, prim_path: str | None = None):
    """OCCT shape(s) from a BrepArray prim in a .usda file.

    With prim_path None the file must hold exactly one BrepArray prim.
    Returns a single shape when the prim holds one Brep, else a compound.
    """
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Compound
    prims = brep_array_prims(path)
    if not prims:
        raise OmniSolidReadError(f"no BrepArray prim in {path}")
    if prim_path is None:
        if len(prims) != 1:
            raise OmniSolidReadError(f"{len(prims)} BrepArray prims; pass prim_path "
                                     f"(one of {sorted(prims)})")
        prim_path = next(iter(prims))
    if prim_path not in prims:
        raise OmniSolidReadError(f"{prim_path} not found (have {sorted(prims)})")
    shapes = _build(prims[prim_path])
    if len(shapes) == 1:
        return shapes[0]
    comp = TopoDS_Compound()
    B = BRep_Builder()
    B.MakeCompound(comp)
    for s in shapes:
        B.Add(comp, s)
    return comp
