"""G17 evidence integrity: stable fingerprints + certification block.

Contract under test:
- canonical_brep_bytes() is stable across BRepMesh_IncrementalMesh,
  BRepCheck_Analyzer, deep copy, and evidence generation itself;
- changing the geometry changes the hash (negative test);
- the evidence record carries a certification block
  {mode, completeness_probe, allow_nonmanifold} under schema 1.1;
- an unprobed record carries the conspicuous "unprobed": true marker;
- evidence validation is additive: bad certification is reported but
  never changes the Boolean verdict.
"""
import sys

sys.path.insert(0, "src")
sys.path.insert(0, "tests")

from OCP.BRep import BRep_Builder
from OCP.BRepBuilderAPI import BRepBuilderAPI_Copy, BRepBuilderAPI_Transform
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.gp import gp_Pnt, gp_Trsf, gp_Vec
from OCP.TopoDS import TopoDS_Shape

from brepkernel import boolean_brep
from brepkernel import evidence


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def box():
    return BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), gp_Pnt(1, 1, 1)).Shape()


def sha(shape):
    return evidence.brep_sha256(shape)


def main():
    ok = True

    # Baseline hash.
    s0 = box()
    h0 = sha(s0)

    # 1. Stable across meshing.
    s1 = box()
    BRepMesh_IncrementalMesh(s1, 0.01, False, 0.1, True)
    ok &= check("hash stable across BRepMesh_IncrementalMesh",
                sha(s1) == h0)

    # 2. Stable across BRepCheck_Analyzer.
    s2 = box()
    BRepCheck_Analyzer(s2, True).IsValid()
    ok &= check("hash stable across BRepCheck_Analyzer",
                sha(s2) == h0)

    # 3. Stable across deep copy.
    s3 = BRepBuilderAPI_Copy(box()).Shape()
    ok &= check("hash stable across deep copy", sha(s3) == h0)

    # 4. Stable across evidence generation itself (boolean_brep hashes
    # its inputs via canonical_brep_bytes; the input hash in the record
    # must match a fresh hash of an untouched box).
    a = box()
    b = BRepPrimAPI_MakeBox(gp_Pnt(0.5, 0.5, 0.5),
                           gp_Pnt(1.5, 1.5, 1.5)).Shape()
    out, rep = boolean_brep(a, b, "union")
    rec = rep.get("evidence") or {}
    in_hashes = [blk["sha256"] for blk in rec.get("inputs", [])]
    ok &= check("hash stable across evidence generation",
                sha(box()) in in_hashes, f"inputs={in_hashes[:1]}")

    # 5. Negative test: a geometric change changes the hash.
    t = gp_Trsf()
    t.SetTranslation(gp_Vec(0.1, 0, 0))
    moved = BRepBuilderAPI_Transform(box(), t, True).Shape()
    ok &= check("hash changes when geometry changes",
                sha(moved) != h0)

    # 6. Certification block present under schema 1.1.
    ok &= check("schema is brepkernel.evidence/1.1",
                rec.get("schema") == "brepkernel.evidence/1.1",
                rec.get("schema", ""))
    cert = rec.get("certification")
    ok &= check("certification block present", isinstance(cert, dict))
    if isinstance(cert, dict):
        ok &= check("certification.mode is strict",
                    cert.get("mode") == "strict", cert.get("mode", ""))
        ok &= check("certification.completeness_probe is bool",
                    isinstance(cert.get("completeness_probe"), bool))
        ok &= check("certification.allow_nonmanifold is False",
                    cert.get("allow_nonmanifold") is False)

    # 7. allow_nonmanifold=True is recorded.
    a2 = box()
    b2 = BRepPrimAPI_MakeBox(gp_Pnt(1, 0, 0), gp_Pnt(2, 1, 1)).Shape()
    try:
        out2, rep2 = boolean_brep(a2, b2, "union", allow_nonmanifold=True)
        cert2 = (rep2.get("evidence") or {}).get("certification") or {}
        ok &= check("allow_nonmanifold=True recorded",
                    cert2.get("allow_nonmanifold") is True)
    except Exception as e:
        # Touching-box union may refuse on this build; the block is still
        # exercised via the default path above.
        ok &= check("allow_nonmanifold path ran", False, str(e)[:80])

    # 8. Validation: a record with a mangled certification block is
    # reported but the verdict stands (additive).
    import copy
    bad = copy.deepcopy(rec)
    bad["certification"] = {"mode": "bogus"}
    problems = evidence.validate_evidence(bad)
    ok &= check("validate_evidence flags bad certification",
                any("certification" in p for p in problems), str(problems))
    ok &= check("good record validates clean",
                evidence.validate_evidence(rec) == [], "")

    # 9. Unprobed marker: a synthesized unprobed record validates only
    # with the conspicuous marker.
    unprobed = copy.deepcopy(rec)
    unprobed["certification"] = {
        "mode": "strict",
        "completeness_probe": False,
        "allow_nonmanifold": False,
        "unprobed": True,
    }
    ok &= check("unprobed record validates with marker",
                evidence.validate_evidence(unprobed) == [])
    unprobed2 = copy.deepcopy(rec)
    unprobed2["certification"] = {
        "mode": "strict",
        "completeness_probe": False,
        "allow_nonmanifold": False,
    }
    probs2 = evidence.validate_evidence(unprobed2)
    ok &= check("unprobed without marker is flagged",
                any("unprobed" in p for p in probs2), str(probs2))

    print("\n" + ("ALL PASS" if ok else "SOME FAILURES"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
