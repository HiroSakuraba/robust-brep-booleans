"""FreeCAD integration adapter tests (no FreeCAD needed).

Tests the BREP-text bridge format, STEP round-trip fidelity, the
run_boolean() subprocess-style boundary (accept/refuse as data, never a
silent wrong solid), and the user-facing refusal rendering.
"""
import json
import os
import sys
import tempfile

sys.path.insert(0, "src")

from brepkernel.freecad_adapter import (
    AdapterError,
    brep_text_of,
    format_status_panel,
    present_refusal,
    report_json,
    round_trip_report,
    run_boolean,
    shape_of_brep_text,
    shape_volume,
    topology_counts,
    write_step,
    read_step,
)

from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.gp import gp_Pnt


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def box(lo, hi):
    return BRepPrimAPI_MakeBox(gp_Pnt(*lo), gp_Pnt(*hi)).Shape()


def t1_brep_text_round_trip():
    a = box((0, 0, 0), (1, 2, 3))
    data = brep_text_of(a)
    b = shape_of_brep_text(data)
    ok = check("fc1 brep text round trip",
               abs(shape_volume(b) - 6.0) < 1e-12
               and topology_counts(b)["faces"] == 6,
               f"vol={shape_volume(b)} faces={topology_counts(b)['faces']}")
    ok &= check("fc1b brep text is ascii",
                all(c < 128 for c in data),
                f"len={len(data)}")
    return ok


def t2_bad_brep_text_is_typed():
    try:
        shape_of_brep_text(b"this is not a brep file")
        return check("fc2 bad brep raises AdapterError", False)
    except AdapterError as e:
        return check("fc2 bad brep raises AdapterError", True,
                     f"kind={e.kind}")
    except Exception as e:  # noqa: BLE001
        return check("fc2 bad brep raises AdapterError", False,
                     f"wrong exception {type(e).__name__}")


def t3_step_round_trip_fidelity():
    a = box((0, 0, 0), (1, 2, 3))
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "a.step")
        write_step(a, p)
        model = read_step(p)
    before = topology_counts(a)
    after = topology_counts(model.shape)
    ok = check("fc3 step round trip preserves topology and volume",
               before == after
               and abs(shape_volume(model.shape) - 6.0) < 1e-9,
               f"before={before} after={after}")
    return ok


def t4_round_trip_report_matched():
    a = box((0, 0, 0), (1, 2, 3))
    rep = round_trip_report(a, fmt="brep")
    ok = check("fc4 round trip report matched",
               rep["matched"] and rep["volume_ok"] and rep["counts_equal"],
               f"matched={rep['matched']}")
    ok &= check("fc4b report json serializes",
                json.dumps(rep) is not None)
    return ok


def t5_run_boolean_accept():
    a = brep_text_of(box((0, 0, 0), (1, 1, 1)))
    b = brep_text_of(box((0.5, 0.5, 0.5), (1.5, 1.5, 1.5)))
    res = run_boolean(a, b, "union")
    ok = check("fc5 accepted union",
               res["outcome"] == "accepted" and res["result_brep"] is not None,
               f"outcome={res['outcome']}")
    out = shape_of_brep_text(res["result_brep"])
    ok &= check("fc5b union volume 1.875",
                abs(shape_volume(out) - 1.875) < 1e-9,
                f"vol={shape_volume(out)}")
    ok &= check("fc5c no refusal block",
                res["refusal"] is None)
    ok &= check("fc5d report json serializes",
                report_json(res) is not None)
    return ok


def t6_run_boolean_refusal_is_data():
    # Two boxes touching at one corner: the kernel must refuse, and the
    # refusal must arrive as data, never as a crash or a wrong solid.
    a = brep_text_of(box((0, 0, 0), (1, 1, 1)))
    b = brep_text_of(box((1, 1, 1), (2, 2, 2)))
    res = run_boolean(a, b, "union")
    ok = check("fc6 corner touch refused",
               res["outcome"] == "refused"
               and res["refusal"]["category"] == "UnresolvedContact"
               and res["refusal"]["stage"] == "assembly"
               and res["result_brep"] is None,
               f"refusal={res['refusal']}")
    text = present_refusal(res["refusal"])
    ok &= check("fc6b panel shows kind and stage",
                "UnresolvedContact" in text and "assembly" in text,
                repr(text[:80]))
    ok &= check("fc6c no em dash anywhere",
                "\u2014" not in text and "\u2014" not in format_status_panel(res))
    return ok


def t7_invalid_op_raises_valueerror():
    a = brep_text_of(box((0, 0, 0), (1, 1, 1)))
    try:
        run_boolean(a, a, "fuse")
        return check("fc7 invalid op raises ValueError", False)
    except ValueError:
        return check("fc7 invalid op raises ValueError", True)
    except Exception as e:  # noqa: BLE001
        return check("fc7 invalid op raises ValueError", False,
                     f"wrong exception {type(e).__name__}")


def t8_status_panel_shape():
    a = brep_text_of(box((0, 0, 0), (1, 1, 1)))
    b = brep_text_of(box((0.5, 0.5, 0.5), (1.5, 1.5, 1.5)))
    panel = format_status_panel(run_boolean(a, b, "union"))
    ok = check("fc8 accept panel has status",
               "Status: Accepted" in panel and "Volume" in panel)
    c = brep_text_of(box((1, 1, 1), (2, 2, 2)))
    panel2 = format_status_panel(run_boolean(a, c, "union"))
    ok &= check("fc8b refuse panel has kind",
                "Status: Refused" in panel2 and "RefusalKind:" in panel2)
    return ok


def main():
    ok = True
    ok &= t1_brep_text_round_trip()
    ok &= t2_bad_brep_text_is_typed()
    ok &= t3_step_round_trip_fidelity()
    ok &= t4_round_trip_report_matched()
    ok &= t5_run_boolean_accept()
    ok &= t6_run_boolean_refusal_is_data()
    ok &= t7_invalid_op_raises_valueerror()
    ok &= t8_status_panel_shape()
    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
