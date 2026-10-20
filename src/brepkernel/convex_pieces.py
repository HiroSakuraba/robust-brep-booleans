"""Convex collision pieces for extruded CAD parts (needs scipy and shapely).

Physics engines collide convex shapes. A part with a hole needs several
convex pieces, or the hole fills in. For a part that is a straight extrusion
(every face is a flat cap across one axis or a wall parallel to it), this
module cuts each constant cross-section into convex polygons and extrudes
them, so the union of the pieces is the part up to the polygon error of its
arcs.

Arcs are polygonized on the material side (see
:func:`brepkernel.mujoco_fit.reference.polygonize`): the pieces never have
material where the part does not, so holes are never smaller than drawn.
The result is certified like the mjfit models: exact convexity, no overlaps
and the right boundary topology (outer loops, holes) in rational arithmetic;
the shape against the part's true cross-sections with a stated bound.

Parts that are not extrusions are refused; a caller can fall back to a
single convex hull.
"""
from __future__ import annotations

import numpy as np

from .mujoco_fit.geom2d import (LoopRegion, PieceRegion, deviation, exact_area, exact_hull, exact_overlap,
                                exact_union_loops, frac_points, loop_signature, sample_loop)
from .mujoco_fit.reference import (Refusal, cap_levels, extrusion_axis, orient, polygonize, section_loops,
                                   to_axis_frame)

SCHEMA = "brepkernel.convex_pieces/1"
MAX_PIECES = 4096


def _cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def _convex(poly) -> bool:
    n = len(poly)
    return all(_cross(poly[i - 1], poly[i], poly[(i + 1) % n]) >= 0 for i in range(n))


def _drop_collinear(poly):
    out = []
    n = len(poly)
    for i in range(n):
        if _cross(poly[i - 1], poly[i], poly[(i + 1) % n]) != 0:
            out.append(poly[i])
    return out


def merge_convex(triangles) -> list:
    """Hertel-Mehlhorn: drop shared diagonals while both sides stay convex (exact)."""
    polys = {}
    for k, t in enumerate(triangles):
        t = list(t)
        if exact_area(t) < 0:
            t.reverse()
        polys[k] = t
    edge_owner = {}
    for k, p in polys.items():
        for i in range(len(p)):
            edge_owner[(p[i], p[(i + 1) % len(p)])] = k
    changed = True
    while changed:
        changed = False
        for (u, v), a in list(edge_owner.items()):
            b = edge_owner.get((v, u))
            if b is None or a == b or a not in polys or b not in polys:
                continue
            pa, pb = polys[a], polys[b]
            ia = pa.index(v)
            ra = pa[ia:] + pa[:ia]                       # v ... u
            ib = pb.index(u)
            rb = pb[ib:] + pb[:ib]                       # u ... v
            if ra[-1] != u or rb[-1] != v:
                continue
            merged = ra + rb[1:-1]
            if not _convex(merged):
                continue
            del polys[b]
            polys[a] = merged
            for i in range(len(merged)):
                edge_owner[(merged[i], merged[(i + 1) % len(merged)])] = a
            edge_owner.pop((u, v), None)
            edge_owner.pop((v, u), None)
            changed = True
    return [_drop_collinear(p) for p in polys.values()]


def _nesting(loops, h):
    import shapely
    polys = [sample_loop(lp, h).points for lp in loops]
    geoms = [shapely.Polygon(p) for p in polys]
    depth, parent = [], []
    for i, p in enumerate(polys):
        inside = [j for j, g in enumerate(geoms) if j != i and g.contains(shapely.Point(*p[0]))]
        depth.append(len(inside))
        parent.append(min(inside, key=lambda j: geoms[j].area) if inside else None)
    return depth, parent


def _slab_pieces(loops, max_error, hc):
    """Exact convex polygons (Fractions) tiling one cross-section."""
    import shapely
    depth, parent = _nesting(loops, hc)
    polys = {}
    for i, lp in enumerate(loops):
        outer = depth[i] % 2 == 0                    # material inside an outer loop
        p = polygonize(orient(lp, not outer, hc), max_error, hc)
        polys[i] = p[::-1] if outer else p           # shells CCW, holes CW
    tris = []
    for i in range(len(loops)):
        if depth[i] % 2:
            continue
        holes = [polys[j] for j in range(len(loops)) if parent[j] == i and depth[j] == depth[i] + 1]
        region = shapely.Polygon(polys[i], holes)
        if not region.is_valid:
            raise Refusal("InvalidSection", f"a cross-section outline is not a valid region "
                          f"({shapely.is_valid_reason(region)})")
        for t in shapely.constrained_delaunay_triangles(region).geoms:
            tris.append(frac_points(np.asarray(t.exterior.coords)[:-1]))
    signature = (sum(1 for d in depth if d % 2 == 0), sum(1 for d in depth if d % 2))
    return merge_convex(tris), signature


def convex_pieces(shape, max_error: float = 1e-4, sample_spacing: float | None = None) -> dict:
    """Report with ``pieces`` (prisms as bottom/top rings in the shape's frame) or a refusal."""
    rep = {"schema": SCHEMA, "status": "certified", "max_error": max_error, "slabs": [], "warnings": []}
    try:
        from OCP.Bnd import Bnd_Box
        from OCP.BRepBndLib import BRepBndLib
        box = Bnd_Box()
        BRepBndLib.Add_s(shape, box)
        size = box.CornerMin().Distance(box.CornerMax())
        axis = extrusion_axis(shape)
        sz, m = to_axis_frame(shape, axis, shift_to_zero=False)
        levels = cap_levels(sz, 1e-9 * size)
        if len(levels) < 2:
            raise Refusal("NotAnExtrusion", "no pair of caps across the extrusion axis")
        h = sample_spacing or max(min(max_error, size) / 50, size * 2e-7)
        hc = size / 4000
        rot, trans = m[:, :3], m[:, 3]
        pieces = []
        for z0, z1 in zip(levels, levels[1:]):
            loops = section_loops(sz, (z0 + z1) / 2, 1e-7 * size)
            if not loops:
                continue                               # a gap between two separate parts
            polys, want = _slab_pieces(loops, max_error, hc)
            if len(pieces) + len(polys) > MAX_PIECES:
                raise Refusal("TooManyPieces", f"more than {MAX_PIECES} convex pieces")
            # --- certificate for this slab (exact, then measured) ---
            hulls = []
            for p in polys:
                hl = exact_hull(p)
                if len(hl) != len(p):
                    raise Refusal("NonConvexPiece", "a merged piece is not strictly convex")
                hulls.append(hl)
            for i in range(len(hulls)):
                for j in range(i + 1, len(hulls)):
                    if exact_overlap(hulls[i], hulls[j]):
                        raise Refusal("PieceOverlap", f"pieces {i} and {j} overlap at z = {z0:g}..{z1:g}")
            sig = loop_signature(exact_union_loops(hulls))
            if tuple(sig) != tuple(want):
                raise Refusal("PieceGap", f"pieces at z = {z0:g}..{z1:g} have {sig[0]} outer and {sig[1]} inner "
                              f"loops, the part has {want[0]} and {want[1]}")
            fl = [np.array([[float(x), float(y)] for x, y in hl]) for hl in hulls]
            dev = deviation(PieceRegion(fl, h, 1e-12), LoopRegion(loops, h))
            res = dev.resolution
            slab = {"z": [z0, z1], "pieces": len(hulls), "outer_loops": sig[0], "hole_loops": sig[1],
                    "deviation": dev.as_dict()}
            rep["slabs"].append(slab)
            if dev.grid_mismatch_far or dev.extra["sim_edge_far_at"] is not None:
                raise Refusal("RegionMismatch", f"pieces at z = {z0:g}..{z1:g} differ from the part's section")
            if dev.outside > res:
                raise Refusal("PiecesOutsidePart", f"pieces reach {dev.outside:.3g} m past the part at z = "
                              f"{z0:g}..{z1:g}; collision must stay inside the material")
            if dev.inside > max_error + res:
                raise Refusal("ShapeMismatch", f"the part reaches {dev.inside:.3g} m past the pieces "
                              f"(allowed {max_error:.3g} m)")
            for p in fl:
                ring0 = np.column_stack([p, np.full(len(p), z0)])
                ring1 = np.column_stack([p, np.full(len(p), z1)])
                # back to the shape's frame: m maps shape -> axis frame (rotation + translation)
                pieces.append({"bottom": ((ring0 - trans) @ rot).tolist(), "top": ((ring1 - trans) @ rot).tolist()})
        rep["pieces"] = pieces
        rep["piece_count"] = len(pieces)
        rep["axis"] = list(map(float, axis))
        rep["resolution"] = 2 * h
        rep["sample_spacing"] = h
        exact_frame = bool(np.all(np.isin(rot, (-1.0, 0.0, 1.0))) and np.all(trans == 0))
        rep["frame_exact"] = exact_frame
        if not exact_frame:
            rep["warnings"].append("the extrusion axis is not a coordinate axis: piece coordinates were rotated "
                                   "back after the exact checks (float rounding, about 1e-16 relative)")
    except Refusal as exc:
        rep["status"] = "refused"
        rep["refusal"] = {"kind": exc.kind, "message": str(exc)}
        rep.pop("pieces", None)
    return rep


__all__ = ["convex_pieces", "merge_convex", "SCHEMA"]
