# Evidence schema and persistent naming

Part B work-plan items B11 (evidence certificate schema) and the
persistent-naming half of B4. Implemented in
`src/brepkernel/evidence.py`, wired into `boolean_brep()` in
`src/brepkernel/pipeline.py`, tested in `tests/test_evidence.py`.

## Schema: `brepkernel.evidence/1.1`

Every `boolean_brep()` call attaches one evidence record at
`report["evidence"]`, on the accept path and on typed-refusal paths
(the refusal exception carries the report, so the record rides along).
The record is a plain JSON-serializable dict:

| Field | Content |
|---|---|
| `schema` | `"brepkernel.evidence/1.1"` |
| `name` | Persistent deterministic name (see below) |
| `naming` | Naming scheme id, formula, op, input hashes, pipeline version and commit |
| `operation_id` | Unique per-call id (uuid4 hex); distinguishes calls that share a name |
| `kernel` | `{name, version, commit, occt, python}` identifying the pipeline build |
| `certification` | `{mode, completeness_probe, allow_nonmanifold}` plus `"unprobed": true` when the completeness probe did not run (see below) |
| `inputs` | Two operand records: `{role, sha256, brep_bytes, format, units, solids, shells, faces}`; `sha256` is the SHA-256 of the canonical BREP text |
| `operation` | `{op, params}` with the effective tolerances and options (base_tol, contact_tol, fuzzy, broadphase_pad + mode, chord_tol, tangent_sin_tol, max_section_tol, area_rel_tol, sew_tol, allow_nonmanifold, crosscheck flags, parallel, use_obb) |
| `timestamps` | `{started_utc, finished_utc, duration_ms}` (ISO 8601 UTC) |
| `outcome` | Accept: `{verdict: "accepted", certificate: {...}}`. Refusal: `{verdict: "refused", category, stage, message, stage_report}` |
| `stages` | Compact per-stage summary on the accept path |
| `artifacts` | `{evidence_file, result_file, log_file}`; paths when written, else null |

The accept certificate summarizes the result without re-running anything:
`result_sha256` (content hash of the result B-rep), `result_brep_bytes`,
`volume`, `solids`, `shells`, `brep_valid`, `closed`, `manifold_edges`,
`complete_edge_lineage`, `volume_bounds_ok`, and `resolution` (one of
`exact_topological_identity`, `strict_same_domain`, `full_pipeline`).

A refusal record carries the typed `category` (the refusal kind, e.g.
`UnresolvedContact`), the pipeline `stage` where it refused, the message,
and the full `stage_report` (the pipeline stages dict, defensively
converted to JSON-safe data).

`evidence.validate_evidence(record)` checks a record against the schema
and returns a list of problem strings (empty means valid). The sidecar
writer refuses to write an invalid record.

## Persistent naming: `brepkernel.naming/1.0`

The evidence name is deterministic and content-derived:

```
name = ev_e10_<op>_<sha256(canonical_brep(A))[:12]>_<sha256(canonical_brep(B))[:12]>_<sha256(version|commit)[:8]>
```

Example: `ev_e10_union_b0caafc89533_73d09a490069_55810891`.

Rules:

- The canonical BREP text is what `BRepTools_Write` emits for the shape
  with triangulations and normals suppressed and the format version
  pinned to 1; the mutable 7-bit TShape flags are normalized to zeros
  before hashing (G17).  It is byte-deterministic for identically
  constructed shapes, so the input hashes are stable across meshing,
  validity checks, deep copies, processes, and machines, and they change
  if and only if the geometry/topology input changes.
- `version` is the pipeline version (`0.9.0`); `commit` is the git HEAD
  SHA at record build time (`"unknown"` outside a git checkout). Any
  code change that moves HEAD changes the trailing code component, so
  evidence from different builds never shares a name silently.
- The name covers the op, the two input hashes (A/B order matters:
  `difference` is not commutative), and the pipeline identity. It
  deliberately excludes timestamps, tolerances, options, and the
  outcome: the name identifies the *operation*, and the record carries
  the rest. Two calls with the same inputs under the same code always
  produce the same evidence name; their `operation_id` and `timestamps`
  tell the calls apart.
- Names are lowercase hex plus underscores: safe as filenames and as
  ledger keys.

When `evidence_dir` is passed to `boolean_brep()`, the record is written
to `<evidence_dir>/<name>.json` (atomic write via temp file + rename).
Re-running the same operation under the same code overwrites the same
file with a fresh record of the latest run; the in-record
`operation_id` and timestamps distinguish the calls.

## Emission guarantees

- Evidence emission never alters the accept/refuse outcome. The wrapper
  builds the record after the pipeline returns or raises, inside a
  guard: any failure in evidence code (hashing, building, writing)
  degrades to `report["evidence_error"]` and the original outcome is
  preserved verbatim. A failing sidecar write (bad `evidence_dir`) is
  recorded in `record["artifacts"]["evidence_write_error"]`, never raised.
- `boolean()` (Tier A mesh API) is untouched; the evidence path applies
  to `boolean_brep()` only (I9).
- No tolerance is read or changed by evidence code (I3).

## Certification block (G17)

`record["certification"]` is copied verbatim from the pipeline's
certification dict:

| Key | Content |
|---|---|
| `mode` | `"strict"` (this pipeline); `"imported_tolerant"` reserved for future import paths |
| `completeness_probe` | `true` when the intersection stage ran its completeness probe |
| `allow_nonmanifold` | The `boolean_brep()` option, as a bool |
| `unprobed` | Present and `true` only when the completeness probe did not run |

An accepted record generated with the completeness probe disabled
carries the conspicuous `"unprobed": true` marker; validators require
it, and CLI tooling must display it.  Evidence validation is additive:
a missing or malformed certification block is reported by
`validate_evidence()` but can never change the Boolean verdict.

## Example (abridged)

```json
{
  "schema": "brepkernel.evidence/1.0",
  "name": "ev_e10_union_b0caafc89533_73d09a490069_55810891",
  "naming": {
    "scheme": "brepkernel.naming/1.0",
    "formula": "ev_e10_<op>_<sha256(canonical_brep(A))[:12]>_<sha256(canonical_brep(B))[:12]>_<sha256(version|commit)[:8]>",
    "op": "union",
    "input_hashes": ["b0caafc8...", "73d09a49..."],
    "pipeline_version": "0.9.0",
    "pipeline_commit": "69318baa..."
  },
  "operation_id": "9f2ac41bd07e44aa8c1e90aa3f1d55c1",
  "kernel": {"name": "brepkernel", "version": "0.9.0",
             "commit": "69318baa...", "occt": "8.0.1", "python": "3.12.3"},
  "inputs": [
    {"role": "A", "sha256": "b0caafc8...",
     "brep_bytes": 2490, "format": "OCCT BREP text via BRepTools_Write",
     "units": "mm", "solids": 1, "shells": 1, "faces": 6},
    {"role": "B", "sha256": "73d09a49...", "brep_bytes": 2490,
     "format": "OCCT BREP text via BRepTools_Write",
     "units": "mm", "solids": 1, "shells": 1, "faces": 6}
  ],
  "operation": {"op": "union", "params": {"base_tol": 1e-07, "...": "..."}},
  "timestamps": {"started_utc": "2026-09-25T23:04:17.556849+00:00",
                 "finished_utc": "2026-09-25T23:04:18.597634+00:00",
                 "duration_ms": 1040.8},
  "outcome": {
    "verdict": "accepted",
    "certificate": {
      "result_sha256": "3f2502bd...",
      "result_brep_bytes": 4821,
      "volume": 1.875,
      "solids": 1, "shells": 1,
      "brep_valid": true, "closed": true, "manifold_edges": true,
      "complete_edge_lineage": true, "volume_bounds_ok": true,
      "resolution": "full_pipeline"
    }
  },
  "stages": {"ingest": {"...": "..."}, "verification": {"...": "..."}},
  "artifacts": {"evidence_file": null, "result_file": null, "log_file": null}
}
```

## Future work (not this change)

- G10 CLI (`brepkernel boolean A.step B.step --op difference -o out.step
  --evidence out.json`) and `brepkernel verify part.step
  part.evidence.json` (re-check the certificate without trusting it).
- Per-entity persistent naming (faces/edges, B4 items 1-4 in the plan):
  deterministic entity names from operation id, operand, parent entity
  names, section references, and piece index; `report["naming"]` mapping
  output entities to input ancestors; `resolve(name, new_result)` with
  geometric fallback. The edge-lineage records in the assembly report
  are the raw material; the naming service is still open.
