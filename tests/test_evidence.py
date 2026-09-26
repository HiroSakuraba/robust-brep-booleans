"""Part B evidence schema (B11) + persistent naming tests.

Contract under test:
- every boolean_brep() call emits a machine-readable evidence record at
  report["evidence"], on the accept path AND on typed-refusal paths
  (the refusal exception carries the report, so the record rides along);
- the record validates against brepkernel.evidence/1.1;
- the record name is deterministic and content-derived:
  name = f(op, input hashes, pipeline version), so the same inputs under
  the same code always produce the same evidence name;
- evidence emission never alters the accept/refuse outcome.
"""
import hashlib
import json
import os
import re
import sys
import tempfile

sys.path.insert(0, "src")
sys.path.insert(0, "tests")

from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.BRepTools import BRepTools
from OCP.gp import gp_Pnt

from brepkernel import boolean_brep, BRepAmbiguousResult
from brepkernel import evidence


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def box(x0, y0, z0, x1, y1, z1):
    return BRepPrimAPI_MakeBox(gp_Pnt(x0, y0, z0),
                              gp_Pnt(x1, y1, z1)).Shape()


def canonical_sha256(shape):
    # G17: use the pipeline's canonical bytes (pinned format, no
    # triangulations/normals, normalized TShape flags).
    return evidence.brep_sha256(shape)


NAME_RE = re.compile(
    r"^ev_e10_(union|intersection|difference)_"
    r"[0-9a-f]{12}_[0-9a-f]{12}_[0-9a-f]{8}$")


def main():
    ok = True

    # 1. Schema module exists with the expected identifiers.
    ok &= check("evidence schema id",
                evidence.SCHEMA_ID == "brepkernel.evidence/1.1",
                evidence.SCHEMA_ID)
    ok &= check("naming scheme id",
                evidence.NAMING_SCHEME_ID == "brepkernel.naming/1.0",
                evidence.NAMING_SCHEME_ID)

    # 2. Accept path emits a valid evidence record.
    a = box(0, 0, 0, 1, 1, 1)
    b = box(0.5, 0.5, 0.5, 1.5, 1.5, 1.5)
    out, rep = boolean_brep(a, b, "union")
    rec = rep.get("evidence")
    ok &= check("accept path attaches evidence record",
                isinstance(rec, dict), type(rec).__name__)
    problems = evidence.validate_evidence(rec)
    ok &= check("accept record validates against schema",
                problems == [], str(problems)[:200])
    ok &= check("accept record JSON serializable",
                bool(json.dumps(rec)), "")
    outcome = rec.get("outcome", {})
    ok &= check("accept verdict", outcome.get("verdict") == "accepted",
                str(outcome.get("verdict")))
    cert = outcome.get("certificate", {})
    ok &= check("certificate has result sha256",
                re.fullmatch(r"[0-9a-f]{64}",
                             str(cert.get("result_sha256") or "")) is not None,
                str(cert.get("result_sha256"))[:16])
    ok &= check("certificate brep_valid", cert.get("brep_valid") is True,
                str(cert.get("brep_valid")))
    ok &= check("certificate volume matches report",
                abs(float(cert.get("volume", -1))
                    - float(rep["stages"]["verification"]["result_volume"]))
                < 1e-12, str(cert.get("volume")))
    ok &= check("record op", rec.get("operation", {}).get("op") == "union",
                str(rec.get("operation", {}).get("op")))
    ok &= check("record has operation_id",
                bool(rec.get("operation_id")), "")
    ok &= check("kernel has version and commit",
                bool(rec.get("kernel", {}).get("version"))
                and bool(rec.get("kernel", {}).get("commit")),
                str(rec.get("kernel")))
    ok &= check("tolerances recorded",
                rec.get("operation", {}).get("params", {}).get("base_tol")
                == 1e-7, str(rec.get("operation", {}).get("params")))
    ts = rec.get("timestamps", {})
    ok &= check("timestamps present",
                bool(ts.get("started_utc")) and bool(ts.get("finished_utc")),
                str(ts))
    ok &= check("artifacts block present",
                isinstance(rec.get("artifacts"), dict), "")

    # 3. Input hashes are content hashes of the operand B-reps.
    inputs = rec.get("inputs", [])
    ok &= check("two inputs recorded", len(inputs) == 2, str(len(inputs)))
    ok &= check("input A hash matches canonical BREP",
                inputs[0].get("sha256") == canonical_sha256(a),
                (inputs[0].get("sha256") or "")[:16])
    ok &= check("input B hash matches canonical BREP",
                inputs[1].get("sha256") == canonical_sha256(b),
                (inputs[1].get("sha256") or "")[:16])
    ok &= check("input roles", [i.get("role") for i in inputs] == ["A", "B"],
                str([i.get("role") for i in inputs]))

    # 4. Persistent naming: deterministic, content-derived.
    name1 = rec.get("name", "")
    ok &= check("name format", NAME_RE.match(name1) is not None, name1)
    out2, rep2 = boolean_brep(box(0, 0, 0, 1, 1, 1),
                              box(0.5, 0.5, 0.5, 1.5, 1.5, 1.5), "union")
    name2 = rep2.get("evidence", {}).get("name", "")
    ok &= check("same inputs same code -> same name", name1 == name2,
                f"{name1} vs {name2}")
    ok &= check("operation_id differs per call",
                rep.get("evidence", {}).get("operation_id")
                != rep2.get("evidence", {}).get("operation_id"), "")
    out3, rep3 = boolean_brep(box(0, 0, 0, 1, 1, 1),
                              box(0.5, 0.5, 0.5, 1.5, 1.5, 2.5), "union")
    name3 = rep3.get("evidence", {}).get("name", "")
    ok &= check("changed input -> changed name", name3 != name1,
                f"{name3} vs {name1}")
    out4, rep4 = boolean_brep(box(0, 0, 0, 1, 1, 1),
                              box(0.5, 0.5, 0.5, 1.5, 1.5, 1.5),
                              "intersection")
    name4 = rep4.get("evidence", {}).get("name", "")
    ok &= check("changed op -> changed name", name4 != name1,
                f"{name4} vs {name1}")
    ok &= check("naming block documents formula",
                "formula" in rec.get("naming", {}),
                str(sorted(rec.get("naming", {}).keys())))

    # 5. Refusal path emits a typed evidence record; outcome unchanged.
    try:
        boolean_brep(box(0, 0, 0, 1, 1, 1), box(1, 1, 1, 2, 2, 2), "union")
        ok &= check("point-touching union still refuses", False,
                    "accepted unexpectedly")
        rec_r = None
        kind = None
    except BRepAmbiguousResult as e:
        kind = (e.report.get("refusal") or {}).get("kind")
        rec_r = e.report.get("evidence")
    ok &= check("point-touching union refuses UnresolvedContact",
                kind == "UnresolvedContact", str(kind))
    ok &= check("refusal path attaches evidence record",
                isinstance(rec_r, dict), str(type(rec_r)))
    problems = evidence.validate_evidence(rec_r)
    ok &= check("refusal record validates against schema",
                problems == [], str(problems)[:200])
    outc = rec_r.get("outcome", {})
    ok &= check("refusal verdict", outc.get("verdict") == "refused",
                str(outc.get("verdict")))
    ok &= check("refusal category typed",
                outc.get("category") == "UnresolvedContact",
                str(outc.get("category")))
    ok &= check("refusal stage recorded", bool(outc.get("stage")),
                str(outc.get("stage")))
    ok &= check("refusal carries stage report",
                isinstance(outc.get("stage_report"), dict)
                and len(outc.get("stage_report", {})) > 0,
                str(sorted(outc.get("stage_report", {}).keys())))
    ok &= check("refusal name is deterministic",
                rec_r.get("name") == evidence.evidence_name(
                    "union", canonical_sha256(box(0, 0, 0, 1, 1, 1)),
                    canonical_sha256(box(1, 1, 1, 2, 2, 2))),
                str(rec_r.get("name")))

    # 6. Emission never alters the outcome, even when the sidecar
    # write fails.
    try:
        out5, rep5 = boolean_brep(box(0, 0, 0, 1, 1, 1),
                                  box(0.5, 0.5, 0.5, 1.5, 1.5, 1.5),
                                  "union",
                                  evidence_dir="/nonexistent-dir-xyz-42")
        vol5 = rep5["stages"]["verification"]["result_volume"]
        ok &= check("broken evidence_dir still accepts at right volume",
                    abs(vol5 - 1.875) < 1e-9, f"vol={vol5}")
        ok &= check("evidence still attached on write failure",
                    isinstance(rep5.get("evidence"), dict), "")
    except Exception as e:
        ok &= check("broken evidence_dir never raises", False,
                    f"{type(e).__name__}: {e}")
    try:
        boolean_brep(box(0, 0, 0, 1, 1, 1), box(1, 1, 1, 2, 2, 2), "union",
                     evidence_dir="/nonexistent-dir-xyz-42")
        ok &= check("broken evidence_dir refusal path still refuses",
                    False, "accepted unexpectedly")
    except BRepAmbiguousResult as e:
        ok &= check("broken evidence_dir refusal path still refuses",
                    (e.report.get("refusal") or {}).get("kind")
                    == "UnresolvedContact", "")

    # 7. evidence_dir writes a sidecar named by the persistent name.
    with tempfile.TemporaryDirectory() as tmp:
        out6, rep6 = boolean_brep(box(0, 0, 0, 1, 1, 1),
                                  box(0.5, 0.5, 0.5, 1.5, 1.5, 1.5),
                                  "union", evidence_dir=tmp)
        rec6 = rep6.get("evidence", {})
        expected = os.path.join(tmp, rec6.get("name", "") + ".json")
        ok &= check("sidecar written at persistent name",
                    os.path.isfile(expected), expected)
        with open(expected) as f:
            disk = json.load(f)
        ok &= check("sidecar content validates",
                    evidence.validate_evidence(disk) == [], "")
        ok &= check("sidecar name matches record",
                    disk.get("name") == rec6.get("name"), "")
        ok &= check("record references sidecar path",
                    rec6.get("artifacts", {}).get("evidence_file")
                    == expected, str(rec6.get("artifacts")))

    # 8. No em dashes in the new evidence sources (I8).
    for rel in ("src/brepkernel/evidence.py", "docs/EVIDENCE_SCHEMA.md",
                "tests/test_evidence.py"):
        p = os.path.join(os.getcwd(), rel)
        with open(p, encoding="utf-8") as f:
            body = f.read()
        ok &= check(f"no em dash in {rel}", "\u2014" not in body, "")

    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
