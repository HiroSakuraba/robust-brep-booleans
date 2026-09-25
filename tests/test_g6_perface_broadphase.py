"""G6 rework: the broad-phase pad is per-face, not per-model.

The G6 per-model pad (2x the largest tolerance anywhere in either model)
is catastrophically pessimistic: one damaged 0.5 mm-tolerance edge makes
every face in the model overlap in broad phase. The rework replaces it
with a per-face conservative pad:

    pad_i = contact_tol + max tolerance over face i and its incident
            edges and vertices   (BRep_Tool.Tolerance on face/edge/vertex)

This test builds a row of 10 boxes where exactly ONE edge of box 0 is
raised to a horrible 0.5 mm tolerance and asserts:

  t1: the per-face pad API exists;
  t2: clean faces (not incident to the damaged edge) keep a small pad
      (<= contact_tol + 1e-6) while the damaged faces carry ~0.5 mm;
  t3: candidate counts with per-face pads do not explode the way the old
      per-model pad did (new count <= old count, and small in absolute
      terms). The inclusion is a theorem, not tuning: every per-face pad
      is pointwise <= the old per-model pad max(contact_tol, 2*tol_max).
"""
import sys

sys.path.insert(0, "src")

from brepkernel.step_ingest import (
    candidate_face_pairs,
    index_shape,
    model_max_tolerance,
)

try:
    from brepkernel.step_ingest import face_broadphase_pads
    HAS_PERFACE = True
except ImportError:
    HAS_PERFACE = False

from OCP.BRep import BRep_Builder
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.TopAbs import TopAbs_EDGE
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS, TopoDS_Compound
from OCP.gp import gp_Pnt

DIRTY_TOL = 5e-4          # 0.5 mm horrible edge tolerance
CONTACT_TOL = 4e-7        # 4 * base_tol, the pipeline default
CLEAN_PAD_CAP = CONTACT_TOL + 1e-6   # OCCT default face tol is 1e-7


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


def make_row_with_dirty_edge(n_boxes=10):
    """Row of unit boxes along x; one edge of box 0 gets DIRTY_TOL."""
    builder = BRep_Builder()
    comp = TopoDS_Compound()
    builder.MakeCompound(comp)
    boxes = []
    for k in range(n_boxes):
        bx = BRepPrimAPI_MakeBox(
            gp_Pnt(2.0 * k, 0.0, 0.0), 1.0, 1.0, 1.0).Shape()
        boxes.append(bx)
        builder.Add(comp, bx)
    # Damage exactly one edge of box 0.
    ex = TopExp_Explorer(boxes[0], TopAbs_EDGE)
    assert ex.More(), "box 0 has no edges"
    dirty = TopoDS.Edge(ex.Current())
    builder.UpdateEdge(dirty, DIRTY_TOL)
    return comp


def make_neighbor():
    """Box hugging the +X end of the row with a 5e-4 gap."""
    return BRepPrimAPI_MakeBox(
        gp_Pnt(2.0 * 10 - 1.0 + 5e-4, 0.0, 0.0), 1.0, 1.0, 1.0).Shape()


def t1_perface_api_exists():
    return check("g6r per-face pad API exists", HAS_PERFACE,
                 "face_broadphase_pads importable from step_ingest")


def t2_clean_face_pads_stay_small():
    m = index_shape(make_row_with_dirty_edge())
    pads = face_broadphase_pads(m, CONTACT_TOL)
    tols = [f.tol_face for f in m.faces]
    dirty = [i for i, t in enumerate(tols) if t >= DIRTY_TOL]
    clean = [i for i, t in enumerate(tols) if t < DIRTY_TOL]
    ok = check("g6r dirty edge visible on its incident faces",
               len(dirty) == 2,
               f"dirty_faces={len(dirty)}")
    ok &= check("g6r dirty faces carry the damaged pad",
                all(pads[i] >= CONTACT_TOL + DIRTY_TOL for i in dirty),
                f"dirty_pads={[f'{pads[i]:.3g}' for i in dirty]}")
    ok &= check("g6r clean faces keep a small per-face pad",
                all(pads[i] <= CLEAN_PAD_CAP for i in clean),
                f"max_clean_pad={max(pads[i] for i in clean):.3g} "
                f"cap={CLEAN_PAD_CAP:.3g} n_clean={len(clean)}")
    return ok


def t3_candidate_counts_do_not_explode():
    ma = index_shape(make_row_with_dirty_edge())
    mb = index_shape(make_neighbor())
    tol_max = max(model_max_tolerance(ma), model_max_tolerance(mb))
    old_pad = max(CONTACT_TOL, 2.0 * tol_max)
    old_count = len(candidate_face_pairs(ma, mb, pad=old_pad))
    pads_a = face_broadphase_pads(ma, CONTACT_TOL)
    pads_b = face_broadphase_pads(mb, CONTACT_TOL)
    new_count = len(candidate_face_pairs(ma, mb, pads_a=pads_a,
                                         pads_b=pads_b))
    ok = check("g6r old per-model pad explodes on the row",
               old_pad >= 2.0 * DIRTY_TOL and old_count > new_count,
               f"old_pad={old_pad:.3g} old_count={old_count} "
               f"new_count={new_count}")
    ok &= check("g6r per-face candidate count stays small",
                new_count <= 6,
                f"new_count={new_count}")
    return ok


def main():
    ok = True
    ok &= t1_perface_api_exists()
    if not ok:
        print("\nSOME FAILURES")
        return 1
    ok &= t2_clean_face_pads_stay_small()
    ok &= t3_candidate_counts_do_not_explode()
    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
