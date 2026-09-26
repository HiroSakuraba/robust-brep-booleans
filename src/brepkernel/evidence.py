"""Machine-readable evidence records for brepkernel boolean outcomes.

Implements the work-plan Part B items B11 (evidence certificate schema)
and the persistent-naming half of B4: every boolean_brep() call attaches
an evidence record to its report, on the accept path and on typed-refusal
paths, and each record carries a deterministic, content-derived name so
the same inputs under the same code always produce the same evidence name.

Schema: ``brepkernel.evidence/1.0``
Naming: ``brepkernel.naming/1.0``

Evidence emission never changes the accept/refuse outcome. The pipeline
wrapper calls build_record() inside a guard: any failure in evidence
code degrades to report["evidence_error"] and the original outcome
(result or refusal) is preserved verbatim.
"""

import hashlib
import math
import os
import platform
import re
import subprocess
import tempfile
import uuid
from datetime import datetime, timezone

import numpy as np
from OCP.BRepTools import BRepTools

SCHEMA_ID = "brepkernel.evidence/1.0"
NAMING_SCHEME_ID = "brepkernel.naming/1.0"

# Pipeline generation this schema describes. The gates branch builds on
# the v0.9 exact trimmed-B-rep pipeline.
KERNEL_VERSION = "0.9.0"

VALID_OPS = ("union", "intersection", "difference")

_NAME_RE = re.compile(
    r"^ev_e10_(union|intersection|difference)_"
    r"[0-9a-f]{12}_[0-9a-f]{12}_[0-9a-f]{8}$")
_SHA_RE = re.compile(r"^[0-9a-f]{64}$")

NAMING_FORMULA = (
    "ev_e10_<op>_<sha256(canonical_brep(A))[:12]>_"
    "<sha256(canonical_brep(B))[:12]>_<sha256(version|commit)[:8]>")


def _repo_root():
    return os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def git_commit():
    """Best-effort pipeline commit hash; 'unknown' outside a git checkout."""
    try:
        out = subprocess.run(
            ["git", "-C", _repo_root(), "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=10)
        if out.returncode == 0 and re.fullmatch(r"[0-9a-f]{40}",
                                                out.stdout.strip()):
            return out.stdout.strip()
    except Exception:
        pass
    return "unknown"


def occt_version():
    """OCCT version from the pinned cadquery-ocp distribution metadata."""
    try:
        from importlib.metadata import version as _pkg_version
        raw = _pkg_version("cadquery-ocp")
        parts = str(raw).split(".")
        if len(parts) >= 3:
            return ".".join(parts[:3])
        return str(raw)
    except Exception:
        return "unknown"


def kernel_info():
    """Version/commit block identifying the pipeline that produced a record."""
    return {
        "name": "brepkernel",
        "version": KERNEL_VERSION,
        "commit": git_commit(),
        "occt": occt_version(),
        "python": platform.python_version(),
    }


def canonical_brep_bytes(shape):
    """Canonical BREP text of a shape, as bytes.

    BRepTools_Write emits entities in construction order, which is
    deterministic for identically constructed shapes, so the bytes are a
    stable content fingerprint of the operand.
    """
    fd, path = tempfile.mkstemp(suffix=".brep")
    os.close(fd)
    try:
        BRepTools.Write_s(shape, path)
        with open(path, "rb") as f:
            return f.read()
    finally:
        try:
            os.unlink(path)
        except OSError:
            pass


def brep_sha256(shape):
    """Hex SHA-256 of the canonical BREP text of a shape."""
    return hashlib.sha256(canonical_brep_bytes(shape)).hexdigest()


def evidence_name(op, sha_a, sha_b, version=None, commit=None):
    """Persistent, deterministic, content-derived evidence name.

    name = f(op, input hashes, pipeline version): the same op on the
    same operand B-reps under the same pipeline version and commit
    always yields the same name. Timestamps, tolerances and the outcome
    are recorded inside the record, never in the name.
    """
    if op not in VALID_OPS:
        raise ValueError(f"unknown op {op!r}")
    if not (_SHA_RE.match(sha_a or "") and _SHA_RE.match(sha_b or "")):
        raise ValueError("input hashes must be 64-char lowercase hex digests")
    version = KERNEL_VERSION if version is None else str(version)
    commit = git_commit() if commit is None else str(commit)
    code = hashlib.sha256(
        f"{version}|{commit}".encode("utf-8")).hexdigest()[:8]
    return f"ev_e10_{op}_{sha_a[:12]}_{sha_b[:12]}_{code}"


def new_operation_id():
    """Unique per-call operation id (distinguishes calls sharing a name)."""
    return uuid.uuid4().hex


def utcnow_iso():
    return datetime.now(timezone.utc).isoformat()


def _jsonable(x, _depth=0):
    """Defensively convert arbitrary report values to JSON-safe data.

    Non-finite floats become None; numpy scalars/arrays become plain
    Python values; unknown objects fall back to repr. Never raises.
    """
    try:
        if _depth > 64:
            return "<max-depth>"
        if x is None or isinstance(x, (bool, int, str)):
            return x
        if isinstance(x, float):
            return x if math.isfinite(x) else None
        if isinstance(x, dict):
            return {str(k): _jsonable(v, _depth + 1)
                    for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [_jsonable(v, _depth + 1) for v in x]
        if isinstance(x, (set, frozenset)):
            return sorted((_jsonable(v, _depth + 1) for v in x), key=repr)
        if isinstance(x, np.ndarray):
            return _jsonable(x.tolist(), _depth + 1)
        if isinstance(x, np.generic):
            return _jsonable(x.item(), _depth + 1)
        return repr(x)
    except Exception:
        try:
            return repr(x)
        except Exception:
            return "<unrepresentable>"


def _ingest_counts(report, role):
    ing = report.get("stages", {}).get("ingest", {}).get(role, {})
    if not isinstance(ing, dict):
        return {"solids": None, "shells": None, "faces": None}
    return {k: ing.get(k) for k in ("solids", "shells", "faces")}


def _resolution(stages):
    if stages.get("identity", {}).get("exact_topological_identity"):
        return "exact_topological_identity"
    if stages.get("same_domain", {}).get("equivalent"):
        return "strict_same_domain"
    return "full_pipeline"


def _stage_summary(stages):
    """Compact per-stage summary for the accept path (certificate view)."""
    summary = {}
    ing = stages.get("ingest")
    if isinstance(ing, dict):
        summary["ingest"] = {
            r: _ingest_counts({"stages": {"ingest": ing}}, r)
            for r in ("A", "B")}
    for key in ("intersection", "topology_stability", "crosscheck"):
        if isinstance(stages.get(key), dict):
            summary[key] = _jsonable(stages[key])
    sp = stages.get("split")
    if isinstance(sp, dict):
        sp = dict(sp)
        uc = sp.get("unresolved_contacts", [])
        sp["unresolved_contacts"] = list(uc)[:64]
        sp["unresolved_contact_count"] = len(uc) if isinstance(uc, list) \
            else None
        summary["split"] = _jsonable(sp)
    asm = stages.get("assembly")
    if isinstance(asm, dict):
        el = asm.get("edge_lineage", {})
        el = el if isinstance(el, dict) else {}
        summary["assembly"] = {
            "selected_faces": asm.get("selected_faces"),
            "shells": asm.get("shells"),
            "solids": asm.get("solids"),
            "free_edges": asm.get("free_edges"),
            "multiple_edges": asm.get("multiple_edges"),
            "volume": asm.get("volume"),
            "empty": asm.get("empty"),
            "notes": list(asm.get("notes", []))[:16],
            "section_sampling": _jsonable(asm.get("section_sampling")),
            "edge_lineage": {
                k: el.get(k) for k in (
                    "result_edges", "boolean_section_edges",
                    "source_boundary_edges", "unattributed_edges",
                    "coincident_boundary_edges")},
            "decisions": len(asm.get("decisions", []))
            if isinstance(asm.get("decisions"), list) else None,
        }
    ver = stages.get("verification")
    if isinstance(ver, dict):
        summary["verification"] = _jsonable(ver)
    for key in ("identity", "same_domain"):
        if isinstance(stages.get(key), dict):
            summary[key] = _jsonable(stages[key])
    return summary


def build_record(*, op, input_a, input_b, params, report,
                 result_shape=None, started_utc=None, finished_utc=None,
                 duration_ms=None, operation_id=None):
    """Build the evidence record dict for one boolean_brep() outcome.

    input_a/input_b: {"sha256": hex, "brep_bytes": int}.
    params: effective tolerance/option dict recorded under operation.params.
    report: the pipeline report dict (stages, timings, refusal).
    result_shape: accepted result shape, or None on the refusal path.
    A refusal is detected from report.get("refusal").
    """
    if op not in VALID_OPS:
        raise ValueError(f"unknown op {op!r}")
    version = KERNEL_VERSION
    commit = git_commit()
    info = kernel_info()
    info["commit"] = commit
    name = evidence_name(op, input_a["sha256"], input_b["sha256"],
                         version=version, commit=commit)
    stages = report.get("stages", {}) if isinstance(report, dict) else {}
    refusal = report.get("refusal") if isinstance(report, dict) else None

    inputs = []
    for role, blk, counts_role in (("A", input_a, "A"), ("B", input_b, "B")):
        inputs.append({
            "role": role,
            "sha256": blk["sha256"],
            "brep_bytes": int(blk["brep_bytes"]),
            "format": "OCCT BREP text via BRepTools_Write",
            "units": "mm",
            **_ingest_counts(report if isinstance(report, dict) else {},
                             counts_role),
        })

    if refusal:
        outcome = {
            "verdict": "refused",
            "category": refusal.get("kind"),
            "stage": refusal.get("stage"),
            "message": refusal.get("message"),
            "stage_report": _jsonable(stages),
        }
        stage_block = {}
    else:
        ver = stages.get("verification", {})
        asm = stages.get("assembly", {})
        result_sha = brep_sha256(result_shape) if result_shape is not None \
            else None
        result_bytes = (len(canonical_brep_bytes(result_shape))
                        if result_shape is not None else None)
        outcome = {
            "verdict": "accepted",
            "certificate": {
                "result_sha256": result_sha,
                "result_brep_bytes": result_bytes,
                "volume": ver.get("result_volume"),
                "solids": asm.get("solids"),
                "shells": asm.get("shells"),
                "brep_valid": ver.get("brep_valid"),
                "closed": ver.get("closed"),
                "manifold_edges": ver.get("manifold_edges"),
                "complete_edge_lineage": ver.get("complete_edge_lineage"),
                "volume_bounds_ok": ver.get("volume_bounds_ok"),
                "resolution": _resolution(stages),
            },
        }
        stage_block = _stage_summary(stages)

    record = {
        "schema": SCHEMA_ID,
        "name": name,
        "naming": {
            "scheme": NAMING_SCHEME_ID,
            "formula": NAMING_FORMULA,
            "op": op,
            "input_hashes": [input_a["sha256"], input_b["sha256"]],
            "pipeline_version": version,
            "pipeline_commit": commit,
        },
        "operation_id": operation_id or new_operation_id(),
        "kernel": info,
        "inputs": inputs,
        "operation": {"op": op, "params": _jsonable(params)},
        "timestamps": {
            "started_utc": started_utc,
            "finished_utc": finished_utc,
            "duration_ms": duration_ms,
        },
        "outcome": outcome,
        "stages": stage_block,
        "artifacts": {
            "evidence_file": None,
            "result_file": None,
            "log_file": None,
        },
    }
    return _jsonable(record)


def validate_evidence(record):
    """Check a record against brepkernel.evidence/1.0.

    Returns a list of problem strings; the empty list means valid.
    """
    problems = []

    def _err(msg):
        problems.append(msg)

    if not isinstance(record, dict):
        return ["record is not a dict"]
    if record.get("schema") != SCHEMA_ID:
        _err(f"schema must be {SCHEMA_ID!r}")
    name = record.get("name")
    if not isinstance(name, str) or not _NAME_RE.match(name):
        _err("name must match the persistent naming format")
    naming = record.get("naming")
    if not isinstance(naming, dict):
        _err("naming block missing")
    else:
        if naming.get("scheme") != NAMING_SCHEME_ID:
            _err(f"naming.scheme must be {NAMING_SCHEME_ID!r}")
        if naming.get("op") not in VALID_OPS:
            _err("naming.op must be a known op")
        if not naming.get("pipeline_version"):
            _err("naming.pipeline_version missing")
        if not naming.get("pipeline_commit"):
            _err("naming.pipeline_commit missing")
    if not isinstance(record.get("operation_id"), str) \
            or not record.get("operation_id"):
        _err("operation_id must be a non-empty string")
    kernel = record.get("kernel")
    if not isinstance(kernel, dict):
        _err("kernel block missing")
    else:
        for k in ("version", "commit", "occt", "python"):
            if not kernel.get(k):
                _err(f"kernel.{k} missing")
    inputs = record.get("inputs")
    if not isinstance(inputs, list) or len(inputs) != 2:
        _err("inputs must be a list of two operand records")
    else:
        for i, blk in enumerate(inputs):
            if not isinstance(blk, dict):
                _err(f"inputs[{i}] is not a dict")
                continue
            if blk.get("role") not in ("A", "B"):
                _err(f"inputs[{i}].role must be 'A' or 'B'")
            if not isinstance(blk.get("sha256"), str) \
                    or not _SHA_RE.match(blk["sha256"]):
                _err(f"inputs[{i}].sha256 must be 64-char hex")
        if isinstance(naming, dict):
            hashes = naming.get("input_hashes")
            got = [b.get("sha256") for b in inputs
                   if isinstance(b, dict)]
            if hashes != got:
                _err("naming.input_hashes must match inputs[*].sha256")
    operation = record.get("operation")
    if not isinstance(operation, dict):
        _err("operation block missing")
    else:
        if operation.get("op") not in VALID_OPS:
            _err("operation.op must be a known op")
        if isinstance(naming, dict) and naming.get("op") != operation.get("op"):
            _err("naming.op and operation.op must agree")
        params = operation.get("params")
        if not isinstance(params, dict) or "base_tol" not in params:
            _err("operation.params must be a dict containing base_tol")
    ts = record.get("timestamps")
    if not isinstance(ts, dict) or not ts.get("started_utc") \
            or not ts.get("finished_utc"):
        _err("timestamps must carry started_utc and finished_utc")
    outcome = record.get("outcome")
    if not isinstance(outcome, dict):
        _err("outcome block missing")
    else:
        verdict = outcome.get("verdict")
        if verdict == "accepted":
            cert = outcome.get("certificate")
            if not isinstance(cert, dict):
                _err("accepted outcome needs a certificate block")
            else:
                if not isinstance(cert.get("result_sha256"), str) \
                        or not _SHA_RE.match(cert["result_sha256"]):
                    _err("certificate.result_sha256 must be 64-char hex")
                if not isinstance(cert.get("brep_valid"), bool):
                    _err("certificate.brep_valid must be a bool")
        elif verdict == "refused":
            if not outcome.get("category") \
                    or not isinstance(outcome.get("category"), str):
                _err("refused outcome needs a category string")
            if not outcome.get("stage") \
                    or not isinstance(outcome.get("stage"), str):
                _err("refused outcome needs a stage string")
            if not isinstance(outcome.get("stage_report"), dict):
                _err("refused outcome needs a stage_report dict")
        else:
            _err("outcome.verdict must be 'accepted' or 'refused'")
    if not isinstance(record.get("stages"), dict):
        _err("stages block must be a dict")
    if not isinstance(record.get("artifacts"), dict):
        _err("artifacts block must be a dict")
    return problems


def write_evidence_file(record, evidence_dir):
    """Write the record as <evidence_dir>/<name>.json; return the path.

    The filename is the persistent evidence name, so re-running the same
    operation under the same code overwrites the same file with a fresh
    record of the latest run (the in-record operation_id and timestamps
    distinguish the calls).
    """
    import json

    problems = validate_evidence(record)
    if problems:
        raise ValueError(
            "refusing to write an invalid evidence record: "
            + "; ".join(problems))
    os.makedirs(evidence_dir, exist_ok=True)
    path = os.path.join(evidence_dir, record["name"] + ".json")
    tmp_path = path + ".tmp"
    with open(tmp_path, "w", encoding="utf-8") as f:
        json.dump(record, f, indent=2, sort_keys=True)
        f.write("\n")
    os.replace(tmp_path, path)
    return path
