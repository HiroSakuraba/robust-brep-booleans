"""brepkernel: robust B-rep booleans prototype.

Topology is decided only by exact tools (analytic implicits here);
geometry is a bounded approximation (certified proxy meshes).
Never silently returns a broken solid.
"""

from .pipeline import (boolean, AmbiguousResult, boolean_brep,
                       BRepAmbiguousResult)
from . import solids

__all__ = ["boolean", "AmbiguousResult", "boolean_brep",
           "BRepAmbiguousResult", "solids"]
__version__ = "0.1.0"
