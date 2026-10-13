"""OpenUSD B-rep exchange (NVIDIA omniSolid BrepArray) and the optional SMLib engine.

Pure brepkernel, no extra dependency:
    write_omnisolid_usda(shape, path)   OCCT solid(s) -> BrepArray .usda
    read_omnisolid_usda(path)           BrepArray .usda -> OCCT solid(s)
    export_usd / import_usd             the same two, under short names

With ``pip install usd-brep`` (brepkernel[usd]):
    brepkernel.usd.smlib_backend        SMLib conversion, Booleans, cross_check
    boolean_brep(..., independent_engine="smlib")
"""
from .omnisolid_read import OmniSolidReadError, brep_array_prims, read_omnisolid_usda
from .omnisolid_write import OmniSolidWriteError, write_omnisolid_usda

export_usd = write_omnisolid_usda
import_usd = read_omnisolid_usda

__all__ = [
    "OmniSolidReadError", "OmniSolidWriteError", "brep_array_prims",
    "export_usd", "import_usd", "read_omnisolid_usda", "write_omnisolid_usda",
]
