# omniSolid BrepArray fixtures written by SMLib

These `.usda` files were written by NVIDIA's SMLib through
`usd_brep.usd.export_brep` / `export_breps` (pip package `usd-brep` 1.0.0,
Apache-2.0). They let `tests/test_usd_omnisolid.py` check brepkernel's
reader against real SMLib output without usd-brep installed.

| file | made with | exact volume |
|---|---|---|
| `sm_boxsph.usda` | box (0,0,0) 2x2x2 minus sphere at (2,2,2) r 1.2 | 8 - (1/8)(4/3)pi 1.2^3 |
| `sm_corner_union.usda` | the same box union the same sphere (boundary runs through a sphere pole) | 8 + (7/8)(4/3)pi 1.2^3 |
| `sm_cylcyl.usda` | cylinder r 1 h 2 union cylinder at (0.5,0,1) r 0.6 h 2 | (checked against OCCT) |
| `sm_split.usda` | box 4x2x2 minus a slab that cuts it in two: one Brep, two solid regions | 12 |
| `multi_sm.usda` | three solids packed as three Breps in one prim (`export_breps`) | 1 + (27 - (4/3)pi) + 2pi |

Regenerate with the snippets in `docs/usd_brep_backend.md`.
