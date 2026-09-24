#!/usr/bin/env python3
"""ABC-dataset independent-oracle stress round for the v0.6 B-rep booleans kernel.

Replicates the method of the NIST clean-geometry and Thingi10K garbage rounds:
  STEP import (OCP) -> tessellate (deflection = diag*1e-3) -> drop zero-area
  tris -> wind per face.Orientation() -> weld ->
  Stage 1: 4 mesh verify fns + defect stats (boundary/non-manifold edges, chi,
           shells, genus estimate) ->
  Stage 2: union/intersection/difference vs a translated copy of itself
           (0.35*diag x-shift, guaranteed overlap), each op in an ISOLATED
           worker process with a 120 s timeout ->
  Oracle (clean inputs only): OCCT BRepAlgoAPI_Fuse of the same pair ->
           union volume % error, solid/shell counts, chi from OCCT tessellation,
           inclusion-exclusion on arranged results ->
  Arbiter (on any union volume mismatch >1% or shell-count disagreement):
           800 random bbox points classified by exact OCCT B-rep
           point-in-solid (BRepClass3d_SolidClassifier, no booleans) vs the
           kernel union mesh via ray-cast parity; report agreement fraction.

Honesty rules (same as the earlier rounds):
  * certified/arranged / refused / unresolved / exceptions / timeouts / crashes
    are separate counts. Zero silent wrong accepts is the bar.
  * Garbage (defective Stage 1) input has no ground truth: arranged results are
    screened for structural soundness (closed, non-empty, finite positive
    volume) and NOT oracle-checked.

Dataset access (documented honestly, 24 Sept 2026):
  * ABC is distributed as 100 chunk 7z files (~1 GB each) on archive.nyu.edu,
    which does NOT honor HTTP Range requests (returns 200 + full body), so a
    subset cannot be extracted without downloading a whole chunk.
  * The documented per-file URL pattern
    https://deep-geometry.github.io/abc-dataset/data/<id>_step_000.step
    returns 404 for every tested model id (0/80 hit rate on a uniform random
    pool drawn from meta-chunk ids of chunks 0/33/66/99).
  * The gh-pages site hosts exactly ONE sample STEP model (plus its features/
    stats yml): 00000050_80d90bfdd2e74e709956122a.
  => This run is a PILOT on the single available ABC model. The script takes a
  directory of *.step files and scales to the full protocol unchanged.

Usage:
  probe_abc.py --steps /tmp/abc_stress/step --out results_abc.json [--seed N]
"""

import argparse
import json
import math
import os
import sys
import time
import traceback
import multiprocessing as mp

import numpy as np

KERNEL_SRC = "/home/hatch/workspace/brep-booleans/src"
sys.path.insert(0, KERNEL_SRC)
from brepkernel import verify  # noqa: E402  (mesh verify fns under test)

WORKER_TIMEOUT = 120.0
ARBITER_N = 800
MISMATCH_PCT = 1.0


# ---------------------------------------------------------------- STEP -> mesh

def step_to_mesh(path):
    """Import STEP via OCP, tessellate, clean, weld. Returns dict or raises."""
    from OCP.STEPControl import STEPControl_Reader
    from OCP.BRepMesh import BRepMesh_IncrementalMesh
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopAbs import TopAbs_FACE, TopAbs_REVERSED
    from OCP.TopoDS import TopoDS
    from OCP.BRep import BRep_Tool
    from OCP.TopLoc import TopLoc_Location
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib

    t0 = time.time()
    r = STEPControl_Reader()
    st = r.ReadFile(path)
    if "RetDone" not in str(st):
        raise RuntimeError(f"STEP read failed: {st}")
    r.TransferRoots()
    shape = r.OneShape()

    box = Bnd_Box()
    BRepBndLib.Add_s(shape, box)
    lo = (box.GetXMin(), box.GetYMin(), box.GetZMin())
    hi = (box.GetXMax(), box.GetYMax(), box.GetZMax())
    diag = math.dist(lo, hi)
    if not (diag > 0 and math.isfinite(diag)):
        raise RuntimeError("degenerate bbox")

    deflection = diag * 1e-3
    BRepMesh_IncrementalMesh(shape, deflection)

    raw = []          # list of (p0,p1,p2) float triples
    n_brep_faces = 0
    n_raw_tris = 0
    ef = TopExp_Explorer(shape, TopAbs_FACE)
    while ef.More():
        face = TopoDS.Face(ef.Current())
        n_brep_faces += 1
        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation_s(face, loc)
        if tri is not None:
            trsf = loc.Transformation()  # OCP 8.0: no args
            rev = face.Orientation() == TopAbs_REVERSED
            for ti in range(1, tri.NbTriangles() + 1):
                n1, n2, n3 = tri.Triangle(ti).Get()
                if len({n1, n2, n3}) < 3:
                    continue  # degenerate node triple (pole)
                pts = []
                for nn in (n1, n2, n3):
                    p = tri.Node(nn)
                    p.Transform(trsf)
                    pts.append((p.X(), p.Y(), p.Z()))
                if rev:
                    pts = [pts[0], pts[2], pts[1]]
                raw.append(pts)
                n_raw_tris += 1
        ef.Next()

    # drop zero-area triangles BEFORE welding (BRepMesh pole/collapsed-edge tris)
    area_tol = 1e-14 * diag * diag
    kept = []
    n_zero_area = 0
    for pts in raw:
        (x0, y0, z0), (x1, y1, z1), (x2, y2, z2) = pts
        ax, ay, az = x1 - x0, y1 - y0, z1 - z0
        bx, by, bz = x2 - x0, y2 - y0, z2 - z0
        cx, cy, cz = ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx
        if (cx * cx + cy * cy + cz * cz) <= area_tol:
            n_zero_area += 1
            continue
        kept.append(pts)

    # weld at scale-aware quantization
    scale = max(1.0, max(abs(c) for t in kept for p in t for c in p)) if kept else 1.0
    q = scale * 1e-9
    verts, V, F = {}, [], []
    n_weld_degen, n_dup = 0, 0
    seen = set()

    def vid(p):
        key = (round(p[0] / q), round(p[1] / q), round(p[2] / q))
        i = verts.get(key)
        if i is None:
            i = len(V)
            verts[key] = i
            V.append(p)
        return i

    for pts in kept:
        t = (vid(pts[0]), vid(pts[1]), vid(pts[2]))
        if len(set(t)) < 3:
            n_weld_degen += 1
            continue
        key = tuple(sorted(t))
        if key in seen:
            n_dup += 1
            continue
        seen.add(key)
        F.append(t)

    return {
        "V": np.array(V, dtype=np.float64),
        "F": np.array(F, dtype=np.int64).reshape(-1, 3),
        "diag": diag,
        "shape": shape,  # kept alive for the OCCT oracle in the main process
        "n_brep_faces": n_brep_faces,
        "n_raw_tris": n_raw_tris,
        "n_zero_area_dropped": n_zero_area,
        "n_weld_degenerate_dropped": n_weld_degen,
        "n_duplicate_dropped": n_dup,
        "tess_time_s": time.time() - t0,
    }


def defect_stats(V, F):
    """Boundary/non-manifold edge counts, chi, shells, per-shell volumes."""
    if len(F):
        e = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
        e = np.sort(e, axis=1)
        _, counts = np.unique(e, axis=0, return_counts=True)
        n_boundary = int(np.sum(counts == 1))
        n_nonmanifold = int(np.sum(counts >= 3))
    else:
        n_boundary = n_nonmanifold = 0
    closed, cinfo = verify.check_closure_directed(V, F)
    chi = verify.euler_chi(V, F)
    shells = verify.connected_shells(F)
    vols = {str(k): float(v) for k, v in verify.shell_signed_volumes(V, F).items()}
    # genus estimate per shell for the closed part (orientable assumption)
    return {
        "n_verts": len(V), "n_tris": len(F),
        "boundary_edges": n_boundary,
        "nonmanifold_edges": n_nonmanifold,
        "directed_closed": bool(closed), "closure_info": cinfo,
        "chi": chi, "shells": shells,
        "shell_signed_volumes": vols,
        "total_abs_volume": float(sum(abs(v) for v in vols.values())),
    }


# ------------------------------------------------------- Stage 2: worker ops

def _worker_main(payload, conn):
    try:
        sys.path.insert(0, payload["src"])
        import numpy as np
        from brepkernel.arrange import arrange, ArrangementError
        from brepkernel import verify as kv

        V = payload["V"]
        F = payload["F"]
        diag = payload["diag"]
        op = payload["op"]
        Vb = V + np.array([0.35 * diag, 0.0, 0.0])
        try:
            res = arrange({"V": V, "F": F}, {"V": Vb, "F": F}, op)
        except ArrangementError as e:
            conn.send({"verdict": "refused", "detail": str(e)[:300]})
            return
        Vr = np.asarray(res["V"], dtype=np.float64)
        Fr = np.asarray(res["F"], dtype=np.int64).reshape(-1, 3)
        closed, cinfo = kv.check_closure_directed(Vr, Fr)
        vol = kv.signed_volume(Vr, Fr)
        finite = bool(np.all(np.isfinite(Vr))) and bool(np.isfinite(vol))
        nonempty = len(Fr) > 0
        conn.send({
            "verdict": "arranged",
            "n_tris": int(len(Fr)),
            "closed": bool(closed), "closure_info": cinfo,
            "signed_volume": float(vol), "volume": float(abs(vol)),
            "finite": finite, "nonempty": nonempty,
            "sane": bool(closed and finite and nonempty and vol != 0.0),
            "Vr": Vr, "Fr": Fr,  # returned for oracle comparison of the union
        })
    except Exception:
        conn.send({"verdict": "exception", "detail": traceback.format_exc()[-800:]})


def run_op_in_worker(V, F, diag, op):
    ctx = mp.get_context("spawn")
    parent_conn, child_conn = ctx.Pipe(duplex=False)
    payload = {"V": np.asarray(V), "F": np.asarray(F),
               "diag": float(diag), "op": op, "src": KERNEL_SRC}
    p = ctx.Process(target=_worker_main, args=(payload, child_conn))
    t0 = time.time()
    p.start()
    child_conn.close()
    # IMPORTANT: never join() before draining the pipe. join-then-recv
    # deadlocks when the result message exceeds the 64KB pipe buffer: the
    # child blocks in send() while the parent blocks in join().
    msg = None
    got_eof = False
    deadline = t0 + WORKER_TIMEOUT
    while True:
        if parent_conn.poll(1.0):
            try:
                msg = parent_conn.recv()
            except EOFError:
                got_eof = True
            break
        if not p.is_alive():
            break
        if time.time() >= deadline:
            break
    dt = time.time() - t0
    if msg is None:
        if p.is_alive():
            p.terminate()
            p.join(10)
            return {"verdict": "timeout",
                    "detail": f"no worker message within {WORKER_TIMEOUT}s",
                    "time_s": dt}
        p.join(5)
        if got_eof or (p.exitcode == 0):
            return {"verdict": "exception",
                    "detail": "worker exited without sending a message",
                    "time_s": dt}
        return {"verdict": "worker_death",
                "detail": f"exitcode={p.exitcode}", "time_s": dt}
    p.join(10)
    msg["time_s"] = dt
    return msg


# ------------------------------------------------------------------ OCCT oracle

def occt_fuse_oracle(shapeA, dx):
    """Fuse shapeA with its x-translated copy. Returns dict or raises.

    NOTE (24 Sept 2026): BRepAlgoAPI_Fuse(compoundA, compoundB) SILENTLY
    drops solids in this OCCT build (3-solid compound in -> 1 solid out,
    IsDone()=True). The oracle therefore explodes the input to solids and
    fuses solid-by-solid (union is associative); every pairwise fuse is
    checked with IsDone().
    """
    from OCP.gp import gp_Trsf, gp_Vec
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopAbs import TopAbs_SOLID, TopAbs_SHELL
    from OCP.TopoDS import TopoDS
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib

    t0 = time.time()
    solids = []
    e = TopExp_Explorer(shapeA, TopAbs_SOLID)
    while e.More():
        solids.append(TopoDS.Solid(e.Current()))
        e.Next()
    if not solids:
        raise RuntimeError("oracle: no solids in input shape")
    trsf = gp_Trsf()
    trsf.SetTranslation(gp_Vec(dx, 0.0, 0.0))
    moved = [BRepBuilderAPI_Transform(s, trsf, True).Shape() for s in solids]
    acc = solids[0]
    for s in solids[1:] + moved:
        f = BRepAlgoAPI_Fuse(acc, s)
        if not f.IsDone():
            raise RuntimeError("BRepAlgoAPI_Fuse did not finish")
        acc = f.Shape()
    fused = acc

    def count(t):
        e = TopExp_Explorer(fused, t)
        n = 0
        while e.More():
            n += 1
            e.Next()
        return n

    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(fused, props)
    vol = props.Mass()
    box = Bnd_Box()
    BRepBndLib.Add_s(fused, box)
    bbox = [(box.GetXMin(), box.GetYMin(), box.GetZMin()),
            (box.GetXMax(), box.GetYMax(), box.GetZMax())]
    n_solids = count(TopAbs_SOLID)
    n_shells = count(TopAbs_SHELL)
    # chi from an OCCT tessellation of the fused shape (independent geometry)
    topo = verify._occt_topology(fused, math.dist(*bbox) * 1e-3)
    return {
        "volume": float(vol), "n_solids": n_solids, "n_shells": n_shells,
        "chi": topo[2], "bbox": bbox, "fused": fused,
        "solids6": solids + moved,  # kept in-process only (not serialized)
        "n_input_solids": len(solids),
        "time_s": time.time() - t0,
    }


# ------------------------------------------------------------------ arbiter

def _ray_parity(pts, V, F):
    """Majority vote over 3 ray directions; returns (inside_mask, undecided)."""
    tv = V[F]  # (T,3,3)
    v0, v1, v2 = tv[:, 0], tv[:, 1], tv[:, 2]
    e1 = v1 - v0
    e2 = v2 - v0
    dirs = [np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([0.31, 0.47, 0.83])]
    dirs = [d / np.linalg.norm(d) for d in dirs]
    votes = np.zeros(len(pts), dtype=int)
    for d in dirs:
        pvec = np.cross(d, e2)          # (T,3)
        det = np.einsum("ij,ij->i", e1, pvec)
        ok = np.abs(det) > 1e-30
        inv = np.zeros_like(det)
        inv[ok] = 1.0 / det[ok]
        for i, o in enumerate(pts):
            tvec = o - v0
            u = np.einsum("ij,ij->i", tvec, pvec) * inv
            qvec = np.cross(tvec, e1)
            vv = np.dot(qvec, d) * inv
            t = np.einsum("ij,ij->i", e2, qvec) * inv
            hit = ok & (u >= 0.0) & (u <= 1.0) & (vv >= 0.0) & (u + vv <= 1.0) & (t > 1e-9)
            if np.count_nonzero(hit) % 2 == 1:
                votes[i] += 1
    inside = votes >= 2
    outside = votes <= 1
    undecided = ~(inside | outside)
    return inside, undecided


def arbiter(fused_shape, fused_bbox, Vr, Fr, solids6=None, n=ARBITER_N, seed=7):
    """800-point independent classification: exact OCCT vs kernel ray parity.

    Two references: (a) the OCCT fused shape via BRepClass3d_SolidClassifier
    (no booleans in the classification itself); (b) when solids6 is given,
    a boolean-free truth: point is IN iff inside ANY input solid.
    Returns agreement fractions for both.
    """
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.gp import gp_Pnt
    from OCP.TopAbs import TopAbs_IN, TopAbs_ON

    rng = np.random.default_rng(seed)
    lo = np.array(fused_bbox[0])
    hi = np.array(fused_bbox[1])
    pts = rng.uniform(lo, hi, size=(n, 3))

    def classify(pt):
        gp = gp_Pnt(float(pt[0]), float(pt[1]), float(pt[2]))
        st = BRepClass3d_SolidClassifier(fused_shape, gp, 1e-9).State()
        return st

    kept, exact_labels = [], []
    n_on = 0
    for p in pts:
        st = classify(p)
        if st == TopAbs_ON:
            n_on += 1
            continue
        kept.append(p)
        exact_labels.append(st == TopAbs_IN)
    kept = np.array(kept)
    exact_labels = np.array(exact_labels, dtype=bool)
    inside, ray_und = _ray_parity(kept, Vr, Fr)
    decided = ~ray_und
    n_decided = int(np.sum(decided))
    agree = int(np.sum(inside[decided] == exact_labels[decided]))
    out = {
        "n_points": n,
        "n_exact_on_surface_undecided": n_on,
        "n_ray_undecided": int(np.sum(ray_und)),
        "n_decided": n_decided,
        "n_agree": agree,
        "agreement_vs_fused": float(agree / n_decided) if n_decided else None,
    }
    if solids6:
        truth, kept2 = [], []
        for p in pts:
            gp = gp_Pnt(float(p[0]), float(p[1]), float(p[2]))
            states = [BRepClass3d_SolidClassifier(s, gp, 1e-9).State()
                      for s in solids6]
            if any(s == TopAbs_ON for s in states):
                continue
            kept2.append(p)
            truth.append(any(s == TopAbs_IN for s in states))
        kept2 = np.array(kept2)
        truth = np.array(truth, dtype=bool)
        inside2, und2 = _ray_parity(kept2, Vr, Fr)
        dec2 = ~und2
        n2 = int(np.sum(dec2))
        ag2 = int(np.sum(inside2[dec2] == truth[dec2]))
        dis = np.where((inside2 != truth) & dec2)[0]
        kin = int(sum(1 for i in dis if inside2[i] and not truth[i]))
        kout = int(sum(1 for i in dis if not inside2[i] and truth[i]))
        out["truth_or"] = {
            "n_decided": n2, "n_agree": ag2,
            "agreement": float(ag2 / n2) if n2 else None,
            "kernel_in_truth_out": kin, "kernel_out_truth_in": kout,
        }
    return out


# ------------------------------------------------------------------ driver

def process_model(mid, step_path):
    rec = {"id": mid, "step": step_path}
    try:
        m = step_to_mesh(step_path)
    except Exception:
        rec["stage1"] = {"verdict": "import_or_tessellation_failed",
                         "detail": traceback.format_exc()[-600:]}
        return rec
    V, F, diag = m["V"], m["F"], m["diag"]
    stats = defect_stats(V, F)
    s1 = {
        "verdict": "ok",
        "n_brep_faces": m["n_brep_faces"],
        "n_raw_tris": m["n_raw_tris"],
        "n_zero_area_dropped": m["n_zero_area_dropped"],
        "n_weld_degenerate_dropped": m["n_weld_degenerate_dropped"],
        "n_duplicate_dropped": m["n_duplicate_dropped"],
        "tess_time_s": round(m["tess_time_s"], 2),
        "diag": diag,
    }
    s1.update(stats)
    try:
        verify.check_closure_directed(V, F)
        verify.euler_chi(V, F)
        verify.connected_shells(F)
        verify.shell_signed_volumes(V, F)
        s1["verify_exceptions"] = 0
    except Exception:
        s1["verify_exceptions"] = 1
        s1["verify_traceback"] = traceback.format_exc()[-600:]
    rec["stage1"] = s1
    # oracle-eligible: input is a valid closed solid (any number of shells);
    # defective (open / non-manifold) inputs have no ground truth.
    rec["clean"] = bool(s1["directed_closed"] and s1["boundary_edges"] == 0
                        and s1["nonmanifold_edges"] == 0 and s1["shells"] == 1
                        and len(F) > 0)
    rec["oracle_eligible"] = bool(s1["directed_closed"]
                                  and s1["boundary_edges"] == 0
                                  and s1["nonmanifold_edges"] == 0
                                  and len(F) > 0)
    rec["_V"], rec["_F"], rec["_shape"] = V, F, m["shape"]
    return rec


def run_stage2(rec):
    V, F = rec["_V"], rec["_F"]
    diag = rec["stage1"]["diag"]
    ops = {}
    for op in ("union", "intersection", "difference"):
        out = run_op_in_worker(V, F, diag, op)
        # don't serialize raw meshes into JSON by default; keep union mesh
        # in memory for the oracle/arbiter comparison
        if out.get("verdict") == "arranged":
            rec[f"_mesh_{op}"] = (out.pop("Vr"), out.pop("Fr"))
        ops[op] = out
    rec["stage2"] = ops


def run_oracle(rec):
    """OCCT fuse oracle + inclusion-exclusion. Only for non-defective inputs."""
    out = {"eligible": bool(rec.get("oracle_eligible"))}
    if not out["eligible"]:
        out["reason"] = ("input defective (Stage 1) -> no ground truth; "
                         "arranged results screened structurally only")
        rec["oracle"] = out
        return
    shapeA = rec["_shape"]
    dx = 0.35 * rec["stage1"]["diag"]
    try:
        oc = occt_fuse_oracle(shapeA, dx)
    except Exception:
        out["oracle_error"] = traceback.format_exc()[-600:]
        rec["oracle"] = out
        return
    out["occt_volume"] = oc["volume"]
    out["occt_n_solids"] = oc["n_solids"]
    out["occt_n_shells"] = oc["n_shells"]
    out["occt_chi"] = oc["chi"]
    out["occt_time_s"] = round(oc["time_s"], 2)
    u = rec["stage2"].get("union", {})
    if u.get("verdict") == "arranged":
        kv = u["volume"]
        ov = oc["volume"]
        out["kernel_union_volume"] = kv
        out["union_vol_err_pct"] = (abs(kv - ov) / ov * 100.0) if ov else None
        ku_shells = verify.connected_shells(rec["_mesh_union"][1])
        ku_chi = verify.euler_chi(rec["_mesh_union"][0], rec["_mesh_union"][1])
        out["kernel_union_shells"] = ku_shells
        out["kernel_union_chi"] = ku_chi
        out["shell_count_agrees"] = (ku_shells == oc["n_shells"])
        out["chi_agrees"] = (oc["chi"] is not None and ku_chi == oc["chi"])
        mismatch = ((out["union_vol_err_pct"] is not None
                     and out["union_vol_err_pct"] > MISMATCH_PCT)
                    or not out["shell_count_agrees"])
        out["arbiter_triggered"] = bool(mismatch)
        if mismatch:
            out["arbiter"] = arbiter(oc["fused"], oc["bbox"],
                                     rec["_mesh_union"][0], rec["_mesh_union"][1],
                                     solids6=oc["solids6"])
    # inclusion-exclusion on arranged results (self-consistency, all inputs)
    ie = {}
    s2 = rec["stage2"]
    vols = {op: s2[op].get("volume") for op in ("union", "intersection", "difference")}
    va = rec["stage1"]["total_abs_volume"]
    if all(s2[op].get("verdict") == "arranged" for op in vols):
        # B is a pure translation of A -> volB == volA
        ie["U_resid"] = abs(vols["union"] - (2 * va - vols["intersection"]))
        ie["U_resid_rel"] = ie["U_resid"] / vols["union"] if vols["union"] else None
        ie["D_resid"] = abs(vols["difference"] - (va - vols["intersection"]))
        ie["D_resid_rel"] = (ie["D_resid"] / vols["difference"]
                             if vols["difference"] else None)
    out["inclusion_exclusion"] = ie
    rec["oracle"] = out


def select_models(recs, n_worst=40, n_clean=10):
    """Worst-case-first selection + clean controls (degrades gracefully)."""
    ok = [r for r in recs if r["stage1"].get("verdict") == "ok"]
    if not ok:
        return []
    def score(r):
        s = r["stage1"]
        return (s["boundary_edges"] + s["nonmanifold_edges"]
                + (0 if s["directed_closed"] else 10**9)
                + abs(s["shells"] - 1) * 1000)
    worst = sorted(ok, key=score, reverse=True)
    cleans = [r for r in worst if r["clean"]]
    picked, seen = [], set()
    for r in worst[:n_worst] + cleans[:n_clean]:
        if r["id"] not in seen:
            picked.append(r)
            seen.add(r["id"])
    return picked


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", required=True,
                    help="directory of *.step files (model id = basename)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=20260924)
    ap.add_argument("--n-worst", type=int, default=40)
    ap.add_argument("--n-clean", type=int, default=10)
    ap.add_argument("--chunk", default=None,
                    help="chunk filename, recorded in meta (added 24 Sept 2026)")
    a = ap.parse_args()

    files = sorted(f for f in os.listdir(a.steps) if f.endswith(".step"))
    print(f"found {len(files)} STEP files in {a.steps}", flush=True)
    t_all = time.time()
    recs = []
    for f in files:
        mid = f[:-5]
        print(f"--- {mid}", flush=True)
        rec = process_model(mid, os.path.join(a.steps, f))
        recs.append(rec)
        s1 = rec["stage1"]
        print(f"    stage1: {s1.get('verdict')} "
              f"tris={s1.get('n_tris')} closed={s1.get('directed_closed')} "
              f"bnd={s1.get('boundary_edges')} nm={s1.get('nonmanifold_edges')} "
              f"shells={s1.get('shells')} chi={s1.get('chi')}", flush=True)

    picked = select_models(recs, a.n_worst, a.n_clean)
    print(f"selected {len(picked)}/{len(recs)} models for Stage 2", flush=True)
    for r in picked:
        r["selected"] = True
    for r in recs:
        r.setdefault("selected", False)

    for r in picked:
        print(f"--- stage2 {r['id']}", flush=True)
        run_stage2(r)
        for op, o in r["stage2"].items():
            print(f"    {op}: {o['verdict']} "
                  f"{('vol=%.6g closed=%s sane=%s' % (o.get('volume'), o.get('closed'), o.get('sane'))) if o['verdict']=='arranged' else o.get('detail','')[:100]} "
                  f"({o.get('time_s',0):.1f}s)", flush=True)
        run_oracle(r)
        oc = r.get("oracle", {})
        if oc.get("eligible") and "union_vol_err_pct" in oc:
            print(f"    oracle: vol_err={oc['union_vol_err_pct']:.4f}% "
                  f"shells k={oc.get('kernel_union_shells')} occt={oc.get('occt_n_shells')} "
                  f"arbiter={oc.get('arbiter_triggered')}", flush=True)

    # strip non-serializable fields
    for r in recs:
        for k in list(r.keys()):
            if k.startswith("_"):
                del r[k]
        r.pop("_V", None)

    summary = {
        "n_step_files": len(files),
        "n_stage1_ok": sum(1 for r in recs if r["stage1"].get("verdict") == "ok"),
        "n_stage1_failed": sum(1 for r in recs if r["stage1"].get("verdict") != "ok"),
        "n_verify_exceptions": sum(1 for r in recs
                                   if r["stage1"].get("verify_exceptions")),
        "n_selected": len(picked),
        "stage2": {},
        "oracle_compared": 0,
        "arbiters_run": 0,
    }
    for op in ("union", "intersection", "difference"):
        c = {}
        for r in picked:
            v = r["stage2"][op]["verdict"]
            c[v] = c.get(v, 0) + 1
        summary["stage2"][op] = c
    for r in picked:
        oc = r.get("oracle", {})
        if oc.get("eligible") and "union_vol_err_pct" in oc:
            summary["oracle_compared"] += 1
        if oc.get("arbiter"):
            summary["arbiters_run"] += 1
    summary["wall_time_s"] = round(time.time() - t_all, 1)

    result = {
        "meta": {
            "date": "2026-09-24",
            "kernel": "/home/hatch/workspace/brep-booleans (v0.6 mesh machinery)",
            "chunk": a.chunk,  # abc chunk filename, None for the pilot
            "protocol": "NIST/Thingi10K replication on ABC dataset (pilot)",
            "sampling": ("uniform random pool of 80 ids from meta chunks 0/33/66/99; "
                         "per-file STEP URLs 404 for all 80 (0% hit rate); "
                         "gh-pages hosts exactly 1 sample STEP; archive.nyu.edu "
                         "chunk 7z files (~1GB) do not support HTTP ranges; "
                         "whole-chunk download declined per task constraints. "
                         "Pilot runs on the single available ABC STEP model."),
            "hit_rate": "0/80 per-file URLs; 1/1 gh-pages sample",
        },
        "summary": summary,
        "models": recs,
    }
    with open(a.out, "w") as fh:
        json.dump(result, fh, indent=1, default=str)
    print("SUMMARY:", json.dumps(summary, indent=1))
    print("wrote", a.out)


if __name__ == "__main__":
    main()
