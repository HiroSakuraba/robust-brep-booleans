"""brepkernel: robust B-rep booleans prototype.

Topology is decided only by exact tools (analytic implicits here);
geometry is a bounded approximation (certified proxy meshes).
Never silently returns a broken solid.
"""

# G18a: apply the OCP compatibility shim before any module that uses
# the TopoDS cast spellings.
from . import _occt_compat  # noqa: F401

from .pipeline import (boolean, AmbiguousResult, boolean_brep,
                       BRepAmbiguousResult)
from . import solids
from . import evidence

__all__ = ["boolean", "AmbiguousResult", "boolean_brep",
           "BRepAmbiguousResult", "solids", "evidence"]

# G18a: version is single-sourced from the installed package metadata.
# Do not hard-code a version string here.
try:
    from importlib.metadata import PackageNotFoundError, version
    __version__ = version("brepkernel")
except PackageNotFoundError:
    __version__ = "0+unknown"
