"""OCP version compatibility shim (G18a).

OCP 7.8 and 8.0 spell some static casts differently.  This module
centralizes the differences so algorithm files do not scatter
try/except ImportError blocks.  Import it before any module that uses
the cast spellings (brepkernel/__init__.py does this).

Compatibility policy: a numerical tolerance difference between OCCT 7.8
and 8.0 may legitimately change an acceptance to a typed refusal; that
is recorded, not tuned away.  A wrong accept or crash on either version
is a blocker.
"""

from OCP.TopoDS import TopoDS

# OCP 8.0 exposes TopoDS.Vertex/Edge/... directly; OCP 7.8 only has the
# trailing-_s static-method spellings.  Alias the new names when missing.
for _name in (
    "Vertex", "Edge", "Wire", "Face", "Shell", "Solid", "CompSolid",
    "Compound",
):
    if not hasattr(TopoDS, _name) and hasattr(TopoDS, _name + "_s"):
        setattr(TopoDS, _name, getattr(TopoDS, _name + "_s"))
del _name


def occt_major_minor():
    """Return (major, minor) of the installed cadquery-ocp, e.g. (8, 0)."""
    try:
        from importlib.metadata import version as _pkg_version
        parts = str(_pkg_version("cadquery-ocp")).split(".")
        return (int(parts[0]), int(parts[1]) if len(parts) > 1 else 0)
    except Exception:
        return (0, 0)
