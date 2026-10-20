"""Read the convex collision pieces out of an MJCF model.

Only what mjfit writes is supported: mesh geoms whose vertices are given
inline (``vertex=``) or in OBJ files, bodies placed with ``pos`` (and an
identity ``quat``). Coordinates are kept as exact rationals: MJCF numbers
written with ``repr`` are exact binary fractions, and body positions are
added exactly, so the structural checks see precisely what MuJoCo reads.
"""
from __future__ import annotations

import xml.etree.ElementTree as ET
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

VISUAL_CLASSES = {"mjfit_visual"}


class MjcfError(ValueError):
    pass


@dataclass
class Piece:
    body: str
    name: str
    verts: list                       # exact (x, y, z) Fractions in world frame

    def as_float(self) -> np.ndarray:
        return np.array([[float(c) for c in v] for v in self.verts])


def _nums(text: str) -> list[Fraction]:
    return [Fraction(float(t)) for t in text.split()]


def _obj_vertices(path: Path) -> list[Fraction]:
    out = []
    for line in path.read_text().splitlines():
        if line.startswith("v "):
            out += _nums(line[2:])
    return out


def read_pieces(xml: str, base_dir: str | Path | None = None) -> dict[str, list[Piece]]:
    """{body name: [Piece, ...]} for every collision mesh geom."""
    root = ET.fromstring(xml)
    meshdir = Path(base_dir or ".")
    comp = root.find("compiler")
    if comp is not None and comp.get("meshdir"):
        meshdir = meshdir / comp.get("meshdir")
    meshes = {}
    asset = root.find("asset")
    for m in (asset.findall("mesh") if asset is not None else []):
        if m.get("vertex"):
            v = _nums(m.get("vertex"))
        elif m.get("file"):
            v = _obj_vertices(meshdir / m.get("file"))
        else:
            continue
        if len(v) % 3:
            raise MjcfError(f"mesh {m.get('name')}: vertex count is not a multiple of 3")
        meshes[m.get("name")] = [tuple(v[i:i + 3]) for i in range(0, len(v), 3)]

    bodies: dict[str, list[Piece]] = {}

    def place(el, origin):
        for key in ("euler", "axisangle", "xyaxes", "zaxis"):
            if el.get(key):
                raise MjcfError(f"{el.tag} {el.get('name')}: rotated frames are not supported")
        quat = el.get("quat")
        if quat and [float(q) for q in quat.split()] != [1.0, 0.0, 0.0, 0.0]:
            raise MjcfError(f"{el.tag} {el.get('name')}: rotated frames are not supported")
        pos = _nums(el.get("pos")) if el.get("pos") else [Fraction(0)] * 3
        return tuple(o + p for o, p in zip(origin, pos))

    def walk(body, origin, name):
        here = place(body, origin) if body.tag == "body" else origin
        for g in body.findall("geom"):
            if g.get("class") in VISUAL_CLASSES:
                continue
            if g.get("contype") == "0" and g.get("conaffinity") == "0":
                continue
            if g.get("type") != "mesh" or g.get("mesh") not in meshes:
                continue
            at = place(g, here)
            verts = [tuple(a + c for a, c in zip(at, v)) for v in meshes[g.get("mesh")]]
            bodies.setdefault(name, []).append(Piece(name, g.get("name") or g.get("mesh"), verts))
        for child in body.findall("body"):
            walk(child, here, child.get("name") or name)

    wb = root.find("worldbody")
    if wb is None:
        raise MjcfError("no <worldbody>")
    walk(wb, (Fraction(0),) * 3, "world")
    return bodies
