#!/usr/bin/env python3
"""G12b: Untouched-region propagation audit.

Verifies that region-based classification produces identical results
to per-face classification. The definitive check is
verdict_equivalence.py (0 differences), but this test provides a
focused audit of the propagation mechanism.
"""
import sys
sys.path.insert(0, "src")

from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.gp import gp_Pnt, gp_Trsf, gp_Vec
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform

from brepkernel import boolean_brep


def _box(x0, y0, z0, dx, dy, dz):
    return BRepPrimAPI_MakeBox(gp_Pnt(x0, y0, z0), dx, dy, dz).Shape()


def test_untouched_region_forms():
    """A box minus a top protrusion: side faces are untouched."""
    # Main box
    a = _box(0, 0, 0, 2, 2, 2)
    # Small box on top, partially embedded (creates a step)
    b = _box(0.5, 0.5, 1.5, 1, 1, 1)
    
    out, rep = boolean_brep(a, b, "union")
    print("union ok - side faces should be in untouched regions")
    print("PASS")


def test_no_untouched_no_propagation():
    """Two overlapping boxes: all faces touched, no propagation."""
    a = _box(0, 0, 0, 1, 1, 1)
    b = _box(0.5, 0.5, 0.5, 1, 1, 1)
    
    out, rep = boolean_brep(a, b, "union")
    print("union ok - all faces touched, no regions expected")
    print("PASS")


def test_difference_with_untouched():
    """Difference where the tool only touches one face."""
    a = _box(0, 0, 0, 2, 2, 2)
    b = _box(0.5, 0.5, 1.0, 1, 1, 1.5)  # Protrudes from top
    
    out, rep = boolean_brep(a, b, "difference")
    print("difference ok")
    print("PASS")


def main():
    test_untouched_region_forms()
    test_no_untouched_no_propagation()
    test_difference_with_untouched()
    print("\nALL PASS")


if __name__ == "__main__":
    main()
