"""Uncertified second opinions for operations brepkernel refuses.

A refusal says "brepkernel could not certify this", not "this has no
answer". When asked, brepkernel computes what two outside engines return
for the same operation and how far apart they are:

* "occt": OCCT's own Boolean (BRepAlgoAPI), always available;
* "smlib": NVIDIA SMLib through brepkernel.usd.smlib_backend, when
  usd-brep is installed.

Nothing here is certified and nothing here changes a verdict. The two
engines share no code, so their agreement is evidence that the case had a
definite answer (and that the refusal may be over-cautious); their
disagreement is evidence the case is genuinely contested.

    triage: "suspected_over_refusal"  both answered, agree to AGREE_REL
            "justified"               both answered, differ by CONTESTED_REL
            "unclear"                 anything else
"""
from __future__ import annotations

from . import _occt_compat  # noqa: F401

AGREE_REL = 1e-6
CONTESTED_REL = 1e-4

_OPS = ("union", "intersection", "difference")


def triage(relative_difference):
    if relative_difference is None:
        return "unclear"
    if relative_difference <= AGREE_REL:
        return "suspected_over_refusal"
    if relative_difference >= CONTESTED_REL:
        return "justified"
    return "unclear"


def _volume(shape) -> float:
    from .usd.smlib_backend import _volume as v
    return v(shape)


def _occt(shapeA, shapeB, op):
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
    from OCP.BRepCheck import BRepCheck_Analyzer
    algo = {"union": BRepAlgoAPI_Fuse, "intersection": BRepAlgoAPI_Common,
            "difference": BRepAlgoAPI_Cut}[op](shapeA, shapeB)
    if not algo.IsDone():
        raise RuntimeError("BRepAlgoAPI reported an error")
    shape = algo.Shape()
    if shape.IsNull() or not BRepCheck_Analyzer(shape).IsValid():
        raise RuntimeError("BRepAlgoAPI result is not a valid B-rep")
    return shape


def second_opinions(shapeA, shapeB, op: str, *, engines=("occt", "smlib")):
    """(summary dict, {engine: shape}) for one operation. Never raises for
    an engine failure; failures are recorded in the summary."""
    if op not in _OPS:
        raise ValueError(f"op must be one of {_OPS}")
    summary = {"certified": False, "engines": {}}
    shapes = {}
    for name in engines:
        rec = {}
        try:
            if name == "occt":
                shp = _occt(shapeA, shapeB, op)
            elif name == "smlib":
                from .usd import smlib_backend as sb
                if not sb.available():
                    rec.update(status="unavailable")
                    summary["engines"][name] = rec
                    continue
                shp = sb.smlib_boolean(shapeA, shapeB, op)
                rec["engine_version"] = sb.smlib_version()
            else:
                raise ValueError(f"unknown engine {name!r}")
            rec.update(status="answered", volume=_volume(shp))
            shapes[name] = shp
        except ValueError:
            raise
        except Exception as exc:
            rec.update(status="failed",
                       error=f"{getattr(exc, 'kind', type(exc).__name__)}: {exc}"[:300])
        summary["engines"][name] = rec
    vols = [r["volume"] for r in summary["engines"].values() if "volume" in r]
    rd = None
    if len(vols) == 2:
        a, b = vols
        rd = abs(a - b) / max(abs(a), abs(b), 1e-300)
    summary["relative_difference"] = rd
    summary["triage"] = triage(rd)
    return summary, shapes
