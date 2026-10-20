"""Certificate for an mjfit peg-in-hole model against the true parts."""
from __future__ import annotations

import re

import numpy as np

from .geom2d import (PieceRegion, Region, deviation, exact_area, exact_hull, exact_overlap, exact_union_loops,
                     loop_signature, min_distance_between)
from .mjcf import read_pieces

SCHEMA = "brepkernel.mjfit/1"

ROLES = [
    ("plate", re.compile(r"^(plate|island)_\d+$"), "plate_hole"),
    ("plate", re.compile(r"^bottom$"), "plate_bottom"),
    ("peg", re.compile(r"^(peg|peg_ring)_\d+$"), "peg_insert"),
    ("peg", re.compile(r"^peg_back_\d+$"), "peg_back"),
    ("peg", re.compile(r"^handle_\d+$"), "handle"),
]

ROLE_TEXT = {
    "plate_hole": "plate around the hole",
    "plate_bottom": "plate bottom",
    "peg_insert": "peg (the part that enters the hole)",
    "peg_back": "peg above the island pockets",
    "handle": "handle",
}


class _Refuse(Exception):
    def __init__(self, kind, message, **evidence):
        super().__init__(message)
        self.kind, self.evidence = kind, evidence


def _prism(piece):
    """(exact CCW hull, z0, z1) of a piece that is a z-prism, else _Refuse."""
    zs = sorted({v[2] for v in piece.verts})
    if len(zs) != 2:
        raise _Refuse("NonPrismPiece", f"piece {piece.name} is not a straight prism along z "
                      f"({len(zs)} vertex heights); this check covers extruded parts only")
    rings = [{(v[0], v[1]) for v in piece.verts if v[2] == z} for z in zs]
    if rings[0] != rings[1]:
        raise _Refuse("NonPrismPiece", f"piece {piece.name}: top and bottom outlines differ")
    hull = exact_hull(list(rings[0]))
    if len(hull) < 3:
        raise _Refuse("DegeneratePiece", f"piece {piece.name} has no area")
    return hull, zs[0], zs[1]


def _floats(hull):
    return np.array([[float(x), float(y)] for x, y in hull])


def certify_model(xml: str, refs: dict[str, Region | None], clearance: float, h: float,
                  base_dir=None, loops: dict[str, tuple[int, int]] | None = None,
                  plate_budget: float | None = None) -> dict:
    """Check the model in `xml` against reference regions per role.

    Structure (exact): every collision piece is a convex prism along z, and
    no two pieces of one body share interior points. Clearance (exact up to
    float distance evaluation): the peg never overlaps the plate, and the
    smallest gap equals `clearance`. Shape (measured, within 2 h): each
    slab's union of pieces against its reference region. Topology (exact):
    each slab's pieces have the expected numbers of outer and inner boundary
    loops (`loops`). Plate slabs must stay within `plate_budget` of the real
    part (the polygon error the outlines allow); peg slabs are reported.
    """
    res = 2 * h                     # every measured number is within this of the true value
    rep: dict = {"schema": SCHEMA, "status": "certified", "resolution": res, "sample_spacing": h,
                 "clearance_requested": clearance, "warnings": [], "slabs": {}}
    try:
        bodies = read_pieces(xml, base_dir)
        groups: dict[str, list] = {}
        unknown = []
        for body, pieces in bodies.items():
            for p in pieces:
                for b, rx, role in ROLES:
                    if b == body and rx.match(p.name):
                        groups.setdefault(role, []).append(p)
                        break
                else:
                    unknown.append(f"{body}/{p.name}")
        if unknown:
            raise _Refuse("UnrecognizedModel", "pieces outside the mjfit peg layout: " + ", ".join(unknown[:5]))
        if "plate_hole" not in groups or "peg_insert" not in groups:
            raise _Refuse("UnrecognizedModel", "need mjfit's plate and peg bodies (plate_*, peg_*)")

        # --- structure (exact) -------------------------------------------
        prisms = {role: [(p, *_prism(p)) for p in ps] for role, ps in groups.items()}
        float32 = 0.0
        for role, items in prisms.items():
            for p, hull, z0, z1 in items:
                f = _floats(hull)
                diam = float(np.linalg.norm(f.max(axis=0) - f.min(axis=0)) + abs(float(z1 - z0)))
                float32 = max(float32, diam * 2.0 ** -24)
        for body in ("plate", "peg"):
            items = [it for role, its in prisms.items() for it in its if it[0].body == body]
            for i in range(len(items)):
                pi, hi, a0, a1 = items[i]
                for j in range(i + 1, len(items)):
                    pj, hj, b0, b1 = items[j]
                    if a0 < b1 and b0 < a1 and exact_overlap(hi, hj):
                        raise _Refuse("PieceOverlap", f"{body} pieces {pi.name} and {pj.name} overlap")
        rep["structure"] = {"pieces": {r: len(v) for r, v in prisms.items()},
                            "convex_prisms": True, "overlapping_pairs": 0, "exact": True}
        rep["mujoco_float32_bound"] = float32

        # --- slabs against references -------------------------------------
        regions = {}
        for role, items in prisms.items():
            zr = sorted({(float(z0), float(z1)) for _, _, z0, z1 in items})
            if len(zr) != 1:
                raise _Refuse("UnrecognizedModel", f"{ROLE_TEXT[role]} pieces have different heights")
            polys = [_floats(hl) for _, hl, _, _ in items]
            union = exact_union_loops([hl for _, hl, _, _ in items])
            regions[role] = PieceRegion(polys, h, eps=1e-12, loops=union)
            sig = loop_signature(union)
            entry = {"pieces": len(items), "z": list(zr[0]),
                     "area_sim": float(sum(exact_area(hl) for _, hl, _, _ in items)),
                     "boundary_loops": len(union), "outer_loops": sig[0], "hole_loops": sig[1]}
            want = (loops or {}).get(role)
            if want is not None and tuple(sig) != tuple(want):
                raise _Refuse("PieceGap", f"the {ROLE_TEXT[role]} pieces have {sig[0]} outer and {sig[1]} "
                              f"inner boundary loops, the part has {want[0]} and {want[1]}: pieces leave "
                              "a gap or are disconnected", slab=role)
            ref = refs.get(role)
            if ref is not None:
                dev = deviation(regions[role], ref)
                entry["deviation"] = dev.as_dict()
                if dev.extra["sim_edge_far_at"] is not None:
                    x, y = dev.extra["sim_edge_far_at"]
                    raise _Refuse("BoundaryMismatch", f"the {ROLE_TEXT[role]} has an edge at ({x:.6g}, {y:.6g}) "
                                  f"{dev.extra['sim_edge_to_reference_max']:.3g} m from any real edge (allowed "
                                  f"{dev.extra['sim_edge_budget']:.3g} m): a notch or slit between pieces",
                                  slab=role)
                if plate_budget is not None and role.startswith("plate"):
                    worst = max(dev.outside, dev.inside)
                    if worst > plate_budget + res:
                        raise _Refuse("ShapeMismatch", f"the {ROLE_TEXT[role]} is {worst:.3g} m from the real "
                                      f"part; its outlines allow {plate_budget:.3g} m", slab=role)
                if dev.grid_mismatch_far:
                    raise _Refuse("RegionMismatch", f"the {ROLE_TEXT[role]} differs from the reference "
                                  "by a whole region, not just along its edges",
                                  slab=role, points=dev.grid_mismatch_far)
            else:
                entry["deviation"] = None
            rep["slabs"][role] = entry

        # --- clearance -----------------------------------------------------
        peg = [it for it in prisms["peg_insert"]]
        plate = [it for it in prisms["plate_hole"]]
        for pi, hi, _, _ in peg:
            for pj, hj, _, _ in plate:
                if exact_overlap(hi, hj):
                    raise _Refuse("SimInterference", f"peg piece {pi.name} overlaps plate piece {pj.name}")
        dmin, at = min_distance_between([_floats(h_) for _, h_, _, _ in peg],
                                        [_floats(h_) for _, h_, _, _ in plate])
        sim_peg, sim_plate = regions["peg_insert"], regions["plate_hole"]
        if sim_plate.contains(sim_peg.boundary()[:1])[0] and dmin > 0:
            raise _Refuse("SimInterference", "the peg lies inside the plate material")
        bpeg = sim_peg.boundary()
        play = float(sim_plate.dist_to_boundary(bpeg).max())
        tol = 1e-12 + 1e-9 * clearance
        rep["clearance"] = {"sim_min": dmin, "sim_min_at": at, "sim_max_play": play,
                            "tolerance": tol}
        if abs(dmin - clearance) > tol:
            raise _Refuse("ClearanceMismatch", f"the smallest peg to plate gap is {dmin:.6g} m, "
                          f"not the requested {clearance:.6g} m", sim_min=dmin)
        ref_plate = refs.get("plate_hole")
        if ref_plate is not None:
            inside = ref_plate.contains(bpeg)
            d = ref_plate.dist_to_boundary(bpeg)
            signed = np.where(inside, -d, d)
            k = int(np.argmin(signed))
            rep["clearance"]["sim_peg_to_real_plate_min"] = float(signed[k])
            rep["clearance"]["sim_peg_to_real_plate_min_at"] = bpeg[k].tolist()
            if signed[k] < -res:
                rep["warnings"].append(
                    f"the simulated peg overlaps the real plate by up to {-signed[k]:.3g} m at "
                    f"({bpeg[k][0]:.6g}, {bpeg[k][1]:.6g}): it is fatter than the real hole allows there")

        # --- findings ------------------------------------------------------
        dv = (rep["slabs"].get("plate_hole") or {}).get("deviation")
        if dv and dv["sim_outside_reference"] > res:
            x, y = dv["sim_outside_reference_at"]
            rep["warnings"].append(
                f"the simulated plate has material up to {dv['sim_outside_reference']:.3g} m beyond the "
                f"real part at ({x:.6g}, {y:.6g}): the hole, or the gap around an island, is smaller "
                "in simulation than drawn there")
        dv = (rep["slabs"].get("peg_insert") or {}).get("deviation")
        if dv and dv["sim_outside_reference"] > res:
            x, y = dv["sim_outside_reference_at"]
            rep["warnings"].append(
                f"the simulated peg is up to {dv['sim_outside_reference']:.3g} m fatter than the ideal "
                f"peg (the real hole shrunk by the clearance) at ({x:.6g}, {y:.6g})")
    except _Refuse as exc:
        rep["status"] = "refused"
        rep["refusal"] = {"kind": exc.kind, "message": str(exc), **exc.evidence}
    return rep


def default_resolution(clearance: float, max_error: float, size: float) -> float:
    """Sample spacing: a fiftieth of the smaller of clearance and max_error."""
    h = min(x for x in (clearance, max_error) if x > 0) / 50 if (clearance > 0 or max_error > 0) else size / 1e4
    return max(h, size * 2e-7)
