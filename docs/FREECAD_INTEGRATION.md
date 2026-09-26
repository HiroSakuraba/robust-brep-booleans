# FreeCAD Integration (Part B, B5): groundwork

FreeCAD is the cheapest and most complete integration point in Part B:
it runs on OCCT (the same kernel family as brepkernel), it is extended
in Python, and its Part module already speaks the same exchange
formats brepkernel uses. This doc records the survey of the integration
surface, the seam design, the refusal UX, and what a first working
demo looks like. FreeCAD is NOT installed on this machine, so there is
no FreeCAD binary to test against; the adapter module and its tests are
deliberately FreeCAD-free (see below).

Status: groundwork only. No workbench, no FreeCADCmd CI job yet.

## 1. How FreeCAD's embedded Python would call this kernel

FreeCAD exposes a full Python API: `Part.Shape` wraps a `TopoDS_Shape`
from FreeCAD's own OCCT build. A `Part::FeaturePython` object with an
`execute()` method is the standard pattern for a parametric feature.

Two integration shapes are possible:

- **In-process:** a workbench imports brepkernel inside FreeCAD's
  interpreter. This is a trap and is rejected: brepkernel uses
  `cadquery-ocp` pinned to OCCT 8.0.1, while FreeCAD ships its own OCCT
  build. Objects from two different OCCT builds cannot be passed to
  each other in one process, and FreeCAD's bundled Python will not
  match the brepkernel venv's Python or numpy builds either.
- **Subprocess (chosen):** FreeCAD runs brepkernel as a child process
  (its own venv, its own Python, its own OCCT) and exchanges data as
  files plus a JSON report. This is the robust option and the one the
  adapter simulates: `run_boolean()` in
  `src/brepkernel/freecad_adapter.py` executes the pipeline in-process
  today, but its signature is exactly the subprocess contract (BREP
  bytes in, outcome dict out), so the production runner only has to
  move the same bytes across a process boundary.

### Version-compat concerns (OCCT)

- brepkernel side: OCCT 8.0.1 via cadquery-ocp==8.0.1.0.0.
- FreeCAD 1.1 side: OCCT version must be confirmed from FreeCAD's own
  release notes or `FreeCAD.Version()` at integration time (FreeCAD 1.0
  shipped OCCT 7.8.x; 1.1 may have moved). Record both versions in
  every evidence record's `kernel` block, since the evidence schema
  already carries `occt`.
- OCCT BREP text is versioned. `BRepTools::Write` accepts a format
  version argument; if FreeCAD's OCCT predates 8.0, write the result
  with an older format version and test the round trip both ways
  before trusting it. The adapter's `round_trip_report()` is the tool
  for that check (counts + volume before/after).
- STEP (AP242) is the safer bridge when OCCT versions differ: it is a
  neutral format both sides read through the same `STEPControl` code
  path. BREP text stays the primary bridge only while both builds
  parse it identically.

## 2. The seam

### Data in

| Route | FreeCAD produces | Adapter consumes |
|---|---|---|
| Primary | `brep_a = obj_a.Shape.exportBrepToString()` | `shape_of_brep_text(bytes)` |
| CAD exchange | Export selection as STEP AP242 | `read_step(path)` (uses `step_ingest.load_step`) |

Units: FreeCAD models in mm. The bridge converts once at ingest and
records the conversion; the tolerance ledger (B3) must carry the host's
declared resolution (OCCT `Precision::Confusion()` = 1e-7 in model
units, usually mm). Per-entity tolerances survive the BREP bridge
verbatim; the STEP bridge carries them through OCCT's uncertainty
mapping, which the adapter must verify with `round_trip_report()`.

### Data out

| Artifact | Produced by | Consumed by |
|---|---|---|
| Result BREP text | `brep_text_of(shape)` | `Part.Shape(); shape.importBrepFromString(...)` |
| Outcome JSON | `report_json(result)` | `CertifiedBoolean.Report` (JSON string property) |
| Evidence sidecar | `result["evidence"]` when the pipeline attached one (schema `brepkernel.evidence/1.0`) | written next to the FreeCAD document as `part.evidence.json`; the `EvidenceFile` property holds its path |

### Evidence as the audit trail

The evidence record is the artifact that makes the FreeCAD result
auditable after the fact: input SHA-256s (over the canonical BREP
text, byte-identical on both sides of the bridge), the op, tolerances,
the accept certificate (volume, solids, shells, validity) or the typed
refusal (category, stage, full stage report). The audit command
re-runs a FreeCAD Boolean through brepkernel, compares with the
independent arbiter, and reports agreement or disagreement points;
the evidence name (`ev_e10_<op>_...`, deterministic and
content-derived) is the ledger key that ties the two runs together.

On this branch the evidence schema branch is not merged, so
`run_boolean()` returns `evidence: None` with the seam documented;
once the branches merge, `report["evidence"]` rides through verbatim.

## 3. STEP round trip: what a round trip preserves and loses

Tested here with OCCT 8.0.1 both ways (`STEPControl_Writer` out,
`step_ingest.load_step` in) on analytic boxes: solids/shells/faces/
edges/vertices counts identical, volume preserved to 1e-9.

- **Preserved:** geometry, topology (solid/shell/face/edge/vertex
  counts), per-entity OCCT tolerances (subject to STEP uncertainty
  quantization; verify per corpus), orientation, units.
- **Lost:** parametric history (the FreeCAD feature tree stays on the
  FreeCAD side), topological naming references (FreeCAD's toponaming
  ids do not survive; brepkernel's deterministic evidence naming is a
  separate scheme, see B4), colors and appearance unless routed through
  `STEPCAFControl`, and assembly/product structure unless routed
  through STEP AP242 PDM entities.
- **Rule:** every import runs validity check (`BRepCheck_Analyzer`),
  recorded healing (`ShapeFix` with before/after diff), and unit
  normalization; healing is never silent and lands in the evidence
  certificate (B8, B11).

## 4. Where typed refusals surface in the FreeCAD UX

Workbench object `CertifiedBoolean` (`Part::FeaturePython`):

- Properties: `Base`, `Tool`, `Operation`, `Status`
  (Accepted / Refused / Error), `RefusalKind`, `Report` (JSON string),
  `EvidenceFile` (path).
- `execute()` calls brepkernel via the subprocess bridge. On refusal,
  the feature goes into the error state and the task panel shows the
  rendered refusal: `present_refusal()` states the category
  (e.g. `UnresolvedContact`), the pipeline stage (e.g. `assembly`),
  and the plain message. It states plainly that no solid was produced
  and that the user may run FreeCAD's own Boolean via an explicit
  button labeled "uncertified".
- Invariants that bind the UX:
  - I1: the feature never shows a solid on a refused operation. There
    is no fallback to `Part::Cut` / `Part::Fuse` on refusal, silent
    or otherwise.
  - I2: the refusal category is typed end to end: kernel refusal kind
    -> `RefusalKind` property -> task panel text. The same string
    appears in the evidence record and in the report JSON.
  - I3: the panel never suggests re-running with a looser tolerance.
- The "Audit" command: select a FreeCAD Boolean feature, re-run it
  with brepkernel, compare with the arbiter, report agreement or the
  disagreement points. This is the "kernel that checks kernels" use
  case from B6, running against FreeCAD's own OCCT Booleans.

## 5. Open questions

1. Which OCCT build does FreeCAD 1.1 actually ship, and does its
   `Part.Shape.importBrepFromString` parse BREP text written by OCCT
   8.0.1? If not, which format version do both sides read?
2. Does FreeCAD's STEP exporter write AP242 with units and the
   uncertainty measure consistent with the tolerance ledger, or must
   the adapter normalize units explicitly?
3. Headless CI: is `FreeCADCmd` from the AppImage or conda-forge
   reliable enough to rebuild all `common_cad_probes.py` cases as
   features and check status + volume (the G11 pass criterion)?
4. Toponaming interaction: when brepkernel returns a result solid into
   a parametric FreeCAD tree, downstream features reference
   brepkernel-named entities; how do FreeCAD's toponaming ids and
   brepkernel's deterministic evidence naming coexist? (Deferred to
   G16; the groundwork only needs the result import to be a static
   imported body.)
5. Refusal wording: the current `present_refusal()` text is plain
   English. Does the workbench need localized or structured
   (machine-readable) refusal fields for downstream features to react
   to? The typed `RefusalKind` property already covers the
   machine-readable half.

## 6. First working demo

A headless script (no GUI FreeCAD needed), runnable in CI:

1. `FreeCADCmd` builds two overlapping boxes and a corner-touching
   pair as `Part` features.
2. Export each operand with `exportBrepToString()`; feed the bytes to
   the subprocess kernel runner (`run_boolean` contract).
3. Accepted case: import the result BREP text, assert volume 1.875
   and `Status == Accepted`.
4. Refused case: assert `Status == Refused`,
   `RefusalKind == UnresolvedContact`, and that no shape was created.
5. Audit command: re-run a native FreeCAD `Part::Fuse` through the
   kernel and report agreement.

Pass criterion (G11): headless FreeCAD CI rebuilds all probes;
refusals appear as feature errors; the audit command works on a
FreeCAD Boolean.

## 7. Blocked on Ben

Nothing in this groundwork needs Ben's keys or files. Onshape API
keys / OAuth app (B6, G13) and CATIA access or sample files (B7, G14)
are separate Part B items and stay parked until Ben provides them.
The FreeCAD track needs only a FreeCAD install for the demo and CI,
which is downloadable, not gated on Ben.

## Files

- `src/brepkernel/freecad_adapter.py`: FreeCAD-free bridge module
  (BREP text in/out, STEP in/out, `run_boolean()`, refusal rendering,
  `round_trip_report()`).
- `tests/test_freecad_adapter.py`: 8 test groups, all passing without
  a FreeCAD install.
