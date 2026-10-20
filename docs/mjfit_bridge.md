# mjfit bridge: certified MuJoCo collision geometry

[mjfit](https://github.com/younghyopark/mjfit) (MIT, Younghyo Park) builds
MuJoCo peg-in-hole models out of convex prisms. MuJoCo only collides convex
shapes, so mjfit cuts the plate around the hole into convex wedges and makes
the peg the hole shrunk by the clearance. Its claim is that the pieces add up
to exactly the part and the peg sits exactly one clearance from the hole.
Curves become polygons drawn around them, so a round hole is never smaller
than drawn.

`brepkernel.mujoco_fit` does two things with it:

- **Mode A, certify mjfit models.** Run `mjfit.peg` and check the model
  against the true parts: real circles, stadiums and Bezier corners rebuilt
  from the same shape specs.
- **Mode C, CAD to MuJoCo.** Read a CAD plate (STEP, BREP or omniSolid USD),
  confirm it is an extrusion with one hole, turn its outline into polygons
  that never make the hole smaller, run `mjfit.peg`, and check the result
  against the CAD part.

Install with `pip install -e .[mjfit]` (mjfit pinned to a tested commit, plus
scipy and shapely).

## Use

```python
from brepkernel.mujoco_fit import certify_mjfit_config, certify_mjfit_peg, mjfit_peg_from_cad

xml, rep = certify_mjfit_config("usb_c")                 # a bundled config, YAML path or dict
xml, rep = certify_mjfit_peg(hole={"type": "circle", "radius": 0.008}, clearance=1e-4)
xml, rep = mjfit_peg_from_cad("plate.step", clearance=1e-4, max_error=5e-5)
rep["status"]        # "certified" or "refused" (then rep["refusal"] = {kind, message})
```

```bash
brepkernel mjfit-certify usb_c -o usb_c.xml --set clearance=8e-5
brepkernel mjfit-from-cad plate.step --clearance 1e-4 --max-error 5e-5 -o plate.xml
```

The MJCF file is written only when the model is certified. The CLI exit codes
are 0 (certified), 2 (refused, nothing written) and 1 (error), and the JSON
report has the certificate under `certificate`.

Give outlines as spec dicts (`{"type": "circle", ...}`) when you want the
true curves checked. A vertex array carries no curve, so it is its own truth.

## What the certificate checks

**Exact checks.** These run in rational arithmetic on the numbers in the MJCF
file. Those are `repr` floats, so they read back exactly as the binary
values MuJoCo gets.

- Every collision piece is a convex prism along z. Threads are helical and
  are refused (`NonPrismPiece`).
- No two pieces of one body overlap (`PieceOverlap`).
- Each slab's pieces form the right topology: the plate around the hole has
  one outer boundary per island plus the plate itself, and exactly one hole.
  Shared edges are cancelled exactly, so a gap of one part in 2^52 between
  two wedges is caught (`PieceGap`).
- The peg never overlaps the plate (`SimInterference`).

**Clearance.** The smallest gap between the peg and the plate in the model
must equal the requested clearance to within 1e-12 m (`ClearanceMismatch`).
Distances come from float segment arithmetic, which is good to about 1e-18 m
at these sizes.

**Measured against the real part.** Both boundaries are sampled at a spacing
h (by default a fiftieth of the smaller of the clearance and the polygon
error). Every measured number is within `resolution` = 2h of the true value.
Per slab:

- `sim_outside_reference`: how far the model's material reaches past the
  real part. For the plate, this is where the hole, or the gap around an
  island, is smaller in simulation than in reality.
- `reference_outside_sim`: how far the real part reaches past the model. For
  the plate, this is the extra play at the polygon corners.
- `sim_edge_to_reference_max`: no model edge may sit far from every real edge.
  This catches a thin notch cut between two pieces (`BoundaryMismatch`).
- A grid check catches whole-region differences, such as a hole that is
  missing entirely (`RegionMismatch`).
- When the outlines come from a hole spec or a CAD plate, the plate must stay
  within the polygon error they allow (`ShapeMismatch`). With a peg spec, the
  real hole is a design choice, so the plate is reported but not held to a
  bound.

The peg is compared with the ideal peg: the real hole shrunk by the
clearance, with round corners where the hole has inside corners. The
certificate also reports `sim_peg_to_real_plate_min`, the smallest gap
between the simulated peg and the real hole. A negative value means the
simulated peg would not fit the real part.

`mujoco_float32_bound` is how far MuJoCo's float32 storage of mesh vertices
can move a vertex (piece size times 2^-24). This is a bound, not a check.

## Mode C: which CAD plates work

- **Shape:** a straight extrusion along some axis (found automatically, so a
  tilted plate works). Every face must be a flat cap across the axis or a
  wall parallel to it (`NotAnExtrusion`).
- **Hole:** one through hole, or one blind hole with a flat bottom. Islands
  inside the hole are allowed, such as a pin standing on the floor of a
  blind hole (`UnsupportedPlate`).
- **Outline:** the outer outline must be convex after polygonizing, as mjfit
  requires (`MjfitRejected`). Outline edges must be lines or circle arcs
  (`UnsupportedOutlineCurve`).
- **Units:** STEP is read as millimeters, BREP and USD as meters.

Polygons always stay inside the material. Arcs around free space get a polygon
drawn outside them; arcs around material get chords. So the hole and the gap
around each island are never smaller than in the CAD part, and the
certificate confirms it (`sim_outside_reference` stays at rounding level).

The peg is generated from the hole and the clearance. The default peg length
is the larger of 30 mm and twice the plate thickness. Other `mjfit.peg`
options pass through (`peg_length`, `handle_height`, `spawn_height`, ...).

## Findings on mjfit's bundled configs

All nine bundled peg configs certify. The measurements found four things:

| config | finding |
|---|---|
| `usb_c` | Clearance 50 um is below the default polygon error (100 um). The simulated peg's corners reach 13.3 um into the real receptacle wall (`sim_peg_to_real_plate_min` = -1.33e-5 m). The island is 54 um fatter than the real tongue. With `clearance=8e-5` the peg clears the real hole by 18 um. |
| `db9`, `xlr3` | The pin islands are drawn as polygons around the circles, so each pin is fatter in simulation: 17.6 um (db9, 12 sides on r = 0.5 mm) and 19.6 um (xlr3). The hole is never smaller, but the gap around a pin is. |
| `rounded_convex` | Bezier outlines are sampled on the curve (mjfit documents this), so the polygonal peg is inside the true peg and the hole made from it is up to 98 um smaller than the true peg grown by the clearance. The clearance is 1 mm, so this is about a tenth of it. |
| `rj45`, `star`, `usb_a` | At inside corners of the hole, mjfit's sharp-cornered peg is thinner than the largest peg that keeps the clearance (by 0.29 x clearance at a 90 degree corner). This is extra play, not interference. |

Each number was checked by hand against the polygon geometry (for example
r (1/cos(pi/n) - 1) for corner reach), and the tests pin them.

## Limits

- Peg-in-hole only. Threads (helical pieces) and gears (involute reference)
  are not covered.
- Mode C handles one hole per plate, lines and circle arcs in outlines, and
  extrusions only.
- The peg in mode C is mjfit's offset of the polygonized hole, not a CAD peg.
  A CAD peg file is not compared yet.
- The measured checks are sampling with a stated bound, not exact. The
  topology, overlap and clearance checks are exact.
