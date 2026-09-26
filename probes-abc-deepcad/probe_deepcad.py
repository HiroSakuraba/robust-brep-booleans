#!/usr/bin/env python3
"""DeepCAD corpus stress round — v0.6 B-rep booleans mesh machinery.

Replicates the NIST (clean) / Thingi10K (garbage) protocol on a sample of the
DeepCAD corpus (neka-nat/DeepCAD-STEP, cad_step.zip) as STEP files:

  Stage 1: STEP import -> OCCT tessellation (deflection = diag*1e-3, drop
           zero-area triangles, wind per face.Orientation(), weld) ->
           run the 4 kernel mesh verify functions; record defect stats.
           Zero exceptions/hangs allowed.
  Stage 2: union / intersection / difference of each mesh against a translated
           copy of itself (0.35*diag shift along +x, guaranteed overlap), each
           op in an ISOLATED worker process with a 120 s timeout.
           Outcome per op: arranged (closed? directed-closure? non-empty?
           finite positive volume?) | refused (typed ArrangementError) |
           exception | timeout | worker death.
  Oracle:  OCCT BRepAlgoAPI_Fuse of the same pair -> union volume % error,
           shell/solid counts, chi from an OCCT tessellation of the fuse.
           Inclusion-exclusion D = A-I, U = A+B-I on arranged results.
  Arbiter: for ANY volume mismatch > 1% or shell-count disagreement, an
           independent 800-point classifier: exact OCCT B-rep point-in-solid
           (BRepClass3d_SolidClassifier, no booleans) vs the kernel union mesh
           via ray-cast parity. Agreement fraction decides kernel vs oracle.

Honesty: certified/arranged, refused, unresolved, exceptions, timeouts and
crashes are separate counts. Refusals are recorded, never dressed up.

Outputs: results_deepcad.json next to this script.
"""

import sys, os, json, math, time, random, traceback
import zipfile

SCRATCH = "/home/hatch/workspace/brep-booleans/probes-abc-deepcad/scratch"
ZIP_PATH = os.path.join(SCRATCH, "cad_step.zip")  # not downloaded; see below
CD_LISTING = os.path.join(SCRATCH, "cd_listing.json")
ZIP_URL = "https://huggingface.co/datasets/neka-nat/DeepCAD-STEP/resolve/main/cad_step.zip"
OUT_JSON = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results_deepcad.json")

sys.path.insert(0, "/home/hatch/workspace/brep-booleans/src")

import numpy as np

SEED = 20260924
POOL_N = 80          # uniform random pool extracted from the zip listing
SELECT_WORST = 40    # worst-case by defect stats
SELECT_CLEAN = 10    # clean controls
OP_TIMEOUT = 120     # seconds per boolean op, isolated worker
ARBITER_N = 800
SHIFT_FRac = 0.35    # self-pair translation = 0.35 * bbox diag

# ---------------------------------------------------------------- tessellation

def tessellate_shape(shape):
    """Tessellate every face of an OCCT shape.

    Returns (V, F, info): V float64 (n,3), F int64 (m,3) welded, info dict with
    raw triangle count, zero-area dropped, face count. Raises on failure.
    """
    from OCP.BRepMesh import BRepMesh_IncrementalMesh
    from OCP.BRep import BRep_Tool
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopAbs import TopAbs_ShapeEnum, TopAbs_Orientation
    from OCP.TopoDS import TopoDS
    from OCP.TopLoc import TopLoc_Location
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib

    box = Bnd_Box()
    BRepBndLib.Add_s(shape, box)
    x0, x1 = box.GetXMin(), box.GetXMax()
    y0, y1 = box.GetYMin(), box.GetYMax()
    z0, z1 = box.GetZMin(), box.GetZMax()
    diag = math.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2 + (z1 - z0) ** 2)
    if not (diag > 0) or not math.isfinite(diag):
        raise ValueError("degenerate bbox")

    deflection = diag * 1e-3
    BRepMesh_IncrementalMesh(shape, deflection)

    tris = []
    n_faces = 0
    n_zero = 0
    exp = TopExp_Explorer(shape, TopAbs_ShapeEnum.TopAbs_FACE)
    while exp.More():
        face = TopoDS.Face(exp.Current())
        n_faces += 1
        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation_s(face, loc)
        if tri is not None and tri.NbTriangles() > 0:
            nn = tri.NbNodes()
            pts = np.empty((nn + 1, 3))
            for i in range(1, nn + 1):
                p = tri.Node(i)
                pts[i] = (p.X(), p.Y(), p.Z())
            trsf = loc.Transformation()  # OCP 8.0: no args
            if trsf.Form() != 0:  # not gp_Identity: apply location
                for i in range(1, nn + 1):
                    p = tri.Node(i).Transformed(trsf)
                    pts[i] = (p.X(), p.Y(), p.Z())
            rev = (face.Orientation() == TopAbs_Orientation.TopAbs_REVERSED)
            for k in range(1, tri.NbTriangles() + 1):
                n1, n2, n3 = tri.Triangle(k).Get()
                a, b, c = pts[n1], pts[n2], pts[n3]
                cr = np.cross(b - a, c - a)
                if float(np.dot(cr, cr)) == 0.0:
                    n_zero += 1
                    continue  # zero-area pole/collapsed triangle: drop
                if rev:
                    tris.append((a, c, b))
                else:
                    tris.append((a, b, c))
        exp.Next()

    if not tris:
        raise ValueError("no triangles tessellated")
    T = np.asarray(tris, dtype=np.float64).reshape(-1, 3)

    # weld vertices (1e-9 absolute in model units; DeepCAD is mm-scale)
    key = np.round(T, 9)
    ukey, inv = np.unique(key.reshape(-1, 3), axis=0, return_inverse=True)
    F = inv.reshape(-1, 3)
    # drop degenerate faces created by welding (duplicate vertex indices)
    keep = (F[:, 0] != F[:, 1]) & (F[:, 1] != F[:, 2]) & (F[:, 0] != F[:, 2])
    F = F[keep]
    info = {"diag": diag, "raw_tris": len(tris) + n_zero, "zero_dropped": n_zero,
            "weld_degenerate_dropped": int(np.sum(~keep)),
            "n_faces_occt": n_faces}
    return ukey.astype(np.float64), F.astype(np.int64), info


def edge_stats(F):
    """Boundary edges (used once) and non-manifold edges (used 3+ times)."""
    from collections import Counter
    c = Counter()
    for a, b, d in F:
        c[(a, b) if a < b else (b, a)] += 1
        c[(b, d) if b < d else (d, b)] += 1
        c[(d, a) if d < a else (a, d)] += 1
    n_bnd = sum(1 for v in c.values() if v == 1)
    n_nm = sum(1 for v in c.values() if v > 2)
    return n_bnd, n_nm


def mesh_volume(V, F):
    if len(F) == 0:
        return 0.0
    v0, v1, v2 = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    return float(np.sum(np.einsum("ij,ij->i", v0, np.cross(v1, v2))) / 6.0)


def stage1_record(name, V, F, tinfo):
    from brepkernel import verify
    rec = {"name": name, "n_verts": int(len(V)), "n_tris": int(len(F))}
    rec.update({k: (float(v) if isinstance(v, float) else v)
                for k, v in tinfo.items()})
    n_bnd, n_nm = edge_stats(F)
    rec["boundary_edges"] = n_bnd
    rec["nonmanifold_edges"] = n_nm
    t0 = time.time()
    try:
        ok, cinfo = verify.check_closure_directed(V, F)
        rec["directed_closure"] = bool(ok)
        rec["closure_info"] = str(cinfo)[:300]
    except Exception as e:
        rec["directed_closure"] = None
        rec["closure_info"] = f"EXCEPTION: {type(e).__name__}: {e}"
    try:
        rec["chi"] = int(verify.euler_chi(V, F))
    except Exception as e:
        rec["chi"] = None
        rec["chi_error"] = f"{type(e).__name__}: {e}"
    try:
        shells = verify.connected_shells(F)
        rec["n_shells"] = len(shells) if hasattr(shells, "__len__") else int(shells)
    except Exception as e:
        rec["n_shells"] = None
        rec["shells_error"] = f"{type(e).__name__}: {e}"
    try:
        sv = verify.shell_signed_volumes(V, F)
        rec["signed_volume"] = float(sum(sv.values()))
        rec["n_shell_vols"] = len(sv)
    except Exception as e:
        rec["signed_volume"] = None
        rec["svol_error"] = f"{type(e).__name__}: {e}"
    # genus estimate for closed orientable: chi = 2 - 2g (per shell); report total
    if rec.get("chi") is not None and rec.get("n_shells"):
        rec["genus_est"] = rec["n_shells"] - rec["chi"] // 2
    rec["stage1_time_s"] = round(time.time() - t0, 2)
    rec["verify_exception"] = any(k.endswith("_error") for k in rec)
    return rec

# ------------------------------------------------------------- boolean worker

def _worker_bool(Va, Fa, Vb, Fb, op):
    """Runs in an isolated process. Returns a plain-dict result."""
    import sys as _s
    _s.path.insert(0, "/home/hatch/workspace/brep-booleans/src")
    import numpy as _np
    from brepkernel.arrange import arrange, ArrangementError
    from brepkernel import verify as _v
    try:
        r = arrange({"V": _np.asarray(Va, dtype=_np.float64),
                     "F": _np.asarray(Fa, dtype=_np.uint32)},
                    {"V": _np.asarray(Vb, dtype=_np.float64),
                     "F": _np.asarray(Fb, dtype=_np.uint32)}, op)
    except ArrangementError as e:
        return {"status": "refused", "error": str(e)[:300]}
    except Exception as e:
        return {"status": "exception",
                "error": f"{type(e).__name__}: {str(e)[:300]}"}
    V = _np.asarray(r["V"], dtype=_np.float64)
    F = _np.asarray(r["F"], dtype=_np.int64).reshape(-1, 3)
    out = {"status": "arranged", "n_verts": int(len(V)),
           "n_tris": int(len(F)), "empty": bool(len(F) == 0)}
    if len(F) == 0:
        out["volume"] = 0.0
        out["finite_positive"] = True
    else:
        try:
            ok, _ = _v.check_closure_directed(V, F)
            out["closed"] = bool(ok)
        except Exception as e:
            out["closed"] = None
            out["closed_error"] = f"{type(e).__name__}: {e}"
        try:
            out["chi"] = int(_v.euler_chi(V, F))
            _sh = _v.connected_shells(F)
            out["n_shells"] = len(_sh) if hasattr(_sh, "__len__") else int(_sh)
        except Exception as e:
            out["topo_error"] = f"{type(e).__name__}: {e}"
        v0, v1, v2 = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
        vol = float(_np.sum(_np.einsum("ij,ij->i", v0, _np.cross(v1, v2))) / 6.0)
        out["volume"] = vol
        out["finite_positive"] = bool(math.isfinite(vol) and vol > 0)
        out["V"] = V.tolist()
        out["F"] = F.tolist()
    return out


def longest_axis(V):
    ext = V.max(axis=0) - V.min(axis=0)
    ax = int(np.argmax(ext))
    return ax, float(ext[ax])


def run_op_isolated(V, F, shift, op):
    """Boolean mesh vs its shifted copy in an isolated worker, 120 s cap."""
    import multiprocessing as mp
    Vb = (V + np.asarray(shift, dtype=np.float64)).astype(np.float64)
    ctx = mp.get_context("spawn")
    pool = ctx.Pool(1)
    try:
        ar = pool.apply_async(_worker_bool, (V, F, Vb, F, op))
        try:
            return ar.get(OP_TIMEOUT)
        except mp.TimeoutError:
            return {"status": "timeout"}
    except Exception as e:
        return {"status": "worker_death",
                "error": f"{type(e).__name__}: {str(e)[:200]}"}
    finally:
        pool.terminate()
        pool.join()

# ------------------------------------------------------------------ OCCT oracle

def occt_fuse_volume_shells_chi(shape_a, shift):
    """OCCT BRepAlgoAPI_Fuse of shape vs its translated copy.

    Returns dict with union volume, solid/shell counts, chi of the fused shape
    (from its own tessellation), and the fused shape itself.
    """
    from OCP.gp import gp_Trsf, gp_Vec
    from OCP.TopLoc import TopLoc_Location
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopAbs import TopAbs_ShapeEnum
    from OCP.GProp import GProp_GProps
    from OCP.BRepGProp import BRepGProp

    trsf = gp_Trsf()
    trsf.SetTranslation(gp_Vec(*shift))
    shape_b = shape_a.Moved(TopLoc_Location(trsf))

    fuse = BRepAlgoAPI_Fuse(shape_a, shape_b)
    fuse.Build()
    if not fuse.IsDone():
        return {"ok": False, "error": "BRepAlgoAPI_Fuse not done"}
    fr = fuse.Shape()
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(fr, props)
    vol = props.Mass()

    def count(kind):
        e = TopExp_Explorer(fr, kind)
        n = 0
        while e.More():
            n += 1
            e.Next()
        return n

    out = {"ok": True, "volume": float(vol),
           "n_solids": count(TopAbs_ShapeEnum.TopAbs_SOLID),
           "n_shells": count(TopAbs_ShapeEnum.TopAbs_SHELL),
           "fused_shape": fr, "shape_b": shape_b}
    try:
        Vf, Ff, _ = tessellate_shape(fr)
        from brepkernel import verify
        out["chi"] = int(verify.euler_chi(Vf, Ff))
        out["fused_tris"] = int(len(Ff))
    except Exception as e:
        out["chi"] = None
        out["chi_error"] = f"{type(e).__name__}: {e}"
    return out


def arbiter(fused_shape, V, F, rng, n=ARBITER_N):
    """800-point independent arbiter.

    OCCT B-rep point-in-solid (BRepClass3d_SolidClassifier, no booleans) vs
    kernel union mesh ray-cast parity. Returns agreement fraction.
    """
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.gp import gp_Pnt
    from OCP.TopAbs import TopAbs_State
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib

    box = Bnd_Box()
    BRepBndLib.Add_s(fused_shape, box)
    x0, x1 = box.GetXMin(), box.GetXMax()
    y0, y1 = box.GetYMin(), box.GetYMax()
    z0, z1 = box.GetZMin(), box.GetZMax()
    pts = np.column_stack([
        rng.uniform(x0, x1, n), rng.uniform(y0, y1, n),
        rng.uniform(z0, z1, n)])
    diag = math.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2 + (z1 - z0) ** 2)
    tol = diag * 1e-7

    occt = np.empty(n, dtype=np.int8)
    for i, (x, y, z) in enumerate(pts):
        cl = BRepClass3d_SolidClassifier(fused_shape, gp_Pnt(x, y, z), tol)
        s = cl.State()
        occt[i] = 1 if s == TopAbs_State.TopAbs_IN else (
            -1 if s == TopAbs_State.TopAbs_OUT else 0)

    # kernel: ray-cast parity (+x ray), vectorized Moller-Trumbore per point
    v0 = V[F[:, 0]]
    e1 = V[F[:, 1]] - v0
    e2 = V[F[:, 2]] - v0
    kern = np.empty(n, dtype=np.int8)
    d = np.array([1.0, 0.0, 0.0])
    for i, p in enumerate(pts):
        pp = p + np.array([1e-9 * diag, 2e-9 * diag, 3e-9 * diag])  # off surface
        h = np.cross(d, e2)
        a = np.einsum("ij,ij->i", e1, h)
        ok = np.abs(a) > 1e-30
        f = 1.0 / np.where(ok, a, 1.0)
        s_ = pp - v0
        u = f * np.einsum("ij,ij->i", s_, h)
        q = np.cross(s_, e1)
        v_ = f * np.dot(q, d)
        t = f * np.einsum("ij,ij->i", e2, q)
        hit = ok & (u >= 0) & (v_ >= 0) & (u + v_ <= 1) & (t > 1e-9)
        kern[i] = 1 if (np.count_nonzero(hit) % 2) else -1

    decided = occt != 0
    agree = float(np.mean(occt[decided] == kern[decided])) if decided.any() else None
    return {"n": n, "n_on_surface": int(n - decided.sum()),
            "agreement": agree}

# ------------------------------------------------------- ranged zip access
# The full zip (3.5 GB compressed / 19.3 GB uncompressed, 212,550 entries) is
# impractical to pull through this proxy, so entries are fetched by HTTP byte
# range using the central directory (downloaded once as cd.bin, parsed into
# cd_listing.json). Sampling is uniform-random over the FULL entry listing.

def fetch_range(url, start, end, tries=4):
    import urllib.request
    last = None
    for _ in range(tries):
        try:
            req = urllib.request.Request(
                url, headers={"Range": f"bytes={start}-{end}"})
            with urllib.request.urlopen(req, timeout=180) as r:
                if r.status != 206:
                    raise RuntimeError(f"range request -> HTTP {r.status}")
                return r.read()
        except Exception as e:
            last = e
    raise RuntimeError(f"range fetch failed: {last}")


def extract_entry(entry):
    """Download one zip entry by byte range; return its raw file bytes."""
    import struct, zlib
    name, lho, comp_size, ucomp_size = entry
    fname_len = len(name.encode())
    # local header + slack for the local extra field, then compressed data
    raw = fetch_range(ZIP_URL, lho, lho + 30 + fname_len + 65536 + comp_size)
    if raw[:4] != b'PK\x03\x04':
        raise ValueError("bad local header signature")
    method, = struct.unpack('<H', raw[8:10])
    fnl, efl = struct.unpack('<2H', raw[26:30])
    data_start = 30 + fnl + efl
    comp = raw[data_start:data_start + comp_size]
    if len(comp) != comp_size:
        raise ValueError("short entry data")
    if method == 8:
        d = zlib.decompressobj(-15)
        data = d.decompress(comp) + d.flush()
    elif method == 0:
        data = comp
    else:
        raise ValueError(f"unsupported zip method {method}")
    if len(data) != ucomp_size:
        raise ValueError("uncompressed size mismatch")
    return data


# ------------------------------------------------------------------------ main

def main():
    t_start = time.time()
    rng = random.Random(SEED)
    nprng = np.random.default_rng(SEED)

    z_names = json.load(open(CD_LISTING))
    names = [tuple(n) for n in z_names
             if n[0].lower().endswith((".step", ".stp"))]
    print(f"zip entries: {len(z_names)}, STEP files: {len(names)}", flush=True)
    listing = {"zip_url": ZIP_URL, "n_entries": len(z_names),
               "n_step": len(names),
               "access": ("ranged HTTP fetch of sampled entries; central "
                          "directory parsed from cd.bin (see cd_listing.json)")}

    pool_names = rng.sample(names, min(POOL_N, len(names)))
    extract_dir = os.path.join(SCRATCH, "extract")
    os.makedirs(extract_dir, exist_ok=True)

    pool = []
    for i, entry in enumerate(pool_names):
        nm = entry[0]
        rec = {"zip_name": nm}
        try:
            data = extract_entry(entry)
            head = data[:2000].decode("ascii", errors="ignore")
            rec["step_header_ok"] = ("ISO-10303-21" in head)
            if not rec["step_header_ok"]:
                rec["stage"] = "bad_header"
                pool.append(rec)
                continue
            fp = os.path.join(extract_dir, f"m{i:03d}.step")
            with open(fp, "wb") as f:
                f.write(data)
            rec["bytes"] = len(data)

            from OCP.STEPControl import STEPControl_Reader
            rdr = STEPControl_Reader()
            if rdr.ReadFile(fp) != 1:
                rec["stage"] = "step_read_fail"
                pool.append(rec)
                continue
            rdr.TransferRoots()
            shape = rdr.OneShape()
            V, F, tinfo = tessellate_shape(shape)
            if len(F) == 0:
                rec["stage"] = "no_tris"
                pool.append(rec)
                continue
            rec.update(stage1_record(nm, V, F, tinfo))
            rec["stage"] = "ok"
            rec["_V"] = V
            rec["_F"] = F
            rec["_shape"] = shape
        except Exception as e:
            rec["stage"] = "exception"
            rec["error"] = f"{type(e).__name__}: {str(e)[:300]}"
        pool.append(rec)
        if (i + 1) % 10 == 0:
            print(f"  pool {i + 1}/{len(pool_names)}", flush=True)

    ok_pool = [r for r in pool if r["stage"] == "ok"]
    print(f"pool ok: {len(ok_pool)}/{len(pool)}", flush=True)

    # ---- selection: ~40 worst-case by defects, ~10 clean controls
    def defect_key(r):
        return (r["boundary_edges"] + r["nonmanifold_edges"],
                r["n_tris"], abs(r.get("chi") or 0))

    defective = sorted([r for r in ok_pool
                        if r["boundary_edges"] + r["nonmanifold_edges"] > 0
                        or not r.get("directed_closure")],
                       key=defect_key, reverse=True)
    clean = [r for r in ok_pool
             if r["boundary_edges"] == 0 and r["nonmanifold_edges"] == 0
             and r.get("directed_closure") and r.get("n_shells") == 1]
    clean_sorted = sorted(clean, key=lambda r: r["n_tris"], reverse=True)

    selected = defective[:SELECT_WORST]
    # fill remainder with most complex if corpus is clean
    if len(selected) < SELECT_WORST:
        rest = sorted([r for r in ok_pool if r not in selected],
                      key=defect_key, reverse=True)
        selected += rest[:SELECT_WORST - len(selected)]
    n_ctrl = min(SELECT_CLEAN, len(clean_sorted))
    # spread controls across complexity: take evenly spaced
    if n_ctrl:
        idx = [round(k * (len(clean_sorted) - 1) / max(1, n_ctrl - 1))
               for k in range(n_ctrl)]
        controls = [clean_sorted[j] for j in idx if clean_sorted[j] not in selected]
    else:
        controls = []
    selected += controls

    selection = {
        "seed": SEED, "pool_n": len(pool_names),
        "n_worst": len([r for r in selected if r in defective[:SELECT_WORST]]),
        "n_controls": len(controls),
        "n_selected": len(selected),
        "note": ("DeepCAD pool had little garbage; worst-case slots filled by "
                 "most complex (face count / |chi|) models" if len(defective) < SELECT_WORST
                 else "worst-case slots filled by defect stats"),
        "names": [r["zip_name"] for r in selected],
    }
    print(f"selected {len(selected)} "
          f"({selection['n_worst']} worst + {selection['n_controls']} controls)", flush=True)

    # ---- stage 2
    results = []
    counts = {"arranged": 0, "arranged_not_closed": 0, "arranged_bad_volume": 0,
              "refused": 0, "exception": 0, "timeout": 0, "worker_death": 0,
              "oracle_fail": 0, "arbiter_runs": 0}
    for si, r in enumerate(selected):
        m = {"name": r["zip_name"], "diag": r["diag"],
             "input_tris": r["n_tris"], "input_shells": r.get("n_shells"),
             "input_chi": r.get("chi"),
             "input_boundary": r["boundary_edges"],
             "input_nonmanifold": r["nonmanifold_edges"],
             "input_closed": r.get("directed_closure")}
        V, F, shape = r["_V"], r["_F"], r["_shape"]
        # shift along the longest bbox axis: 0.35*diag still overlaps because
        # the longest extent is >= diag/sqrt(3) > 0.35*diag
        ax, _ext = longest_axis(V)
        shift = [0.0, 0.0, 0.0]
        shift[ax] = SHIFT_FRac * r["diag"]
        m["shift_axis"] = ax
        vol_a = abs(r["signed_volume"]) if r.get("signed_volume") else None
        m["vol_a"] = vol_a
        m["ops"] = {}
        t0 = time.time()
        for op in ("union", "intersection", "difference"):
            or_ = run_op_isolated(V, F, shift, op)
            o = {"status": or_["status"]}
            st = or_["status"]
            if st == "arranged":
                counts["arranged"] += 1
                o.update({k: or_.get(k) for k in
                          ("n_verts", "n_tris", "empty", "closed", "chi",
                           "n_shells", "volume", "finite_positive")})
                if or_.get("empty"):
                    counts["arranged_empty"] = counts.get("arranged_empty", 0) + 1
                else:
                    if or_.get("closed") is not True:
                        counts["arranged_not_closed"] += 1
                    if not or_.get("finite_positive"):
                        counts["arranged_bad_volume"] += 1
                if st == "arranged" and "V" in or_ and op == "union":
                    m["_union_mesh"] = (np.array(or_["V"]), np.array(or_["F"]))
            elif st == "refused":
                counts["refused"] += 1
                o["error"] = or_.get("error")
            elif st == "exception":
                counts["exception"] += 1
                o["error"] = or_.get("error")
            elif st == "timeout":
                counts["timeout"] += 1
            else:
                counts["worker_death"] += 1
                o["error"] = or_.get("error")
            m["ops"][op] = o

        # oracle (main process)
        try:
            oc = occt_fuse_volume_shells_chi(shape, shift)
        except Exception as e:
            oc = {"ok": False, "error": f"{type(e).__name__}: {str(e)[:200]}"}
        m["oracle"] = {k: v for k, v in oc.items()
                       if k not in ("fused_shape", "shape_b")}
        if not oc.get("ok"):
            counts["oracle_fail"] += 1

        # inclusion-exclusion on arranged results
        vols = {op: m["ops"][op].get("volume")
                for op in ("union", "intersection", "difference")}
        m["inclusion_exclusion"] = None
        if all(isinstance(vols[op], float) for op in vols) and vol_a:
            # U should equal A+B-I ; D should equal A-I ; B is congruent to A
            ie_u = abs(vols["union"] - (2 * vol_a - vols["intersection"]))
            ie_d = abs(vols["difference"] - (vol_a - vols["intersection"]))
            scale = max(vol_a, 1e-30)
            m["inclusion_exclusion"] = {
                "resid_union": ie_u / scale, "resid_diff": ie_d / scale}

        # kernel vs oracle comparison on union
        m["oracle_compare"] = None
        need_arbiter = False
        if oc.get("ok") and isinstance(vols["union"], float) and oc["volume"] > 0:
            verr = abs(vols["union"] - oc["volume"]) / oc["volume"]
            ku_shells = m["ops"]["union"].get("n_shells")
            cmp_ = {"vol_rel_err": verr,
                    "kernel_shells": ku_shells,
                    "occt_solids": oc.get("n_solids"),
                    "occt_shells": oc.get("n_shells"),
                    "kernel_chi": m["ops"]["union"].get("chi"),
                    "occt_chi": oc.get("chi")}
            m["oracle_compare"] = cmp_
            if verr > 0.01:
                need_arbiter = True
                cmp_["flag"] = "VOLUME_MISMATCH"
            if ku_shells is not None and oc.get("n_solids") is not None \
                    and ku_shells != oc["n_solids"]:
                need_arbiter = True
                cmp_["flag"] = (cmp_.get("flag", "") + "+SHELL_COUNT_MISMATCH").strip("+")

        if need_arbiter and "_union_mesh" in m and oc.get("ok"):
            counts["arbiter_runs"] += 1
            try:
                Vu, Fu = m["_union_mesh"]
                ab = arbiter(oc["fused_shape"], Vu, Fu, nprng)
                m["arbiter"] = ab
            except Exception as e:
                m["arbiter"] = {"error": f"{type(e).__name__}: {str(e)[:200]}"}
        m.pop("_union_mesh", None)

        m["model_time_s"] = round(time.time() - t0, 1)
        results.append(m)
        # free heavy refs
        r.pop("_V", None); r.pop("_F", None); r.pop("_shape", None)
        print(f"  model {si + 1}/{len(selected)} {m['name'].split('/')[-1][:40]} "
              f"ops={[m['ops'][o]['status'] for o in ('union','intersection','difference')]} "
              f"{m['model_time_s']}s", flush=True)

    out = {
        "meta": {
            "date": "2026-09-24", "corpus": "neka-nat/DeepCAD-STEP cad_step.zip",
            "kernel": "/home/hatch/workspace/brep-booleans (v0.6 mesh machinery)",
            "seed": SEED, "op_timeout_s": OP_TIMEOUT,
            "shift_frac": SHIFT_FRac, "arbiter_n": ARBITER_N,
            "wall_time_s": round(time.time() - t_start, 1),
        },
        "listing": listing,
        "pool_summary": {
            "n": len(pool),
            "stages": {s: sum(1 for r in pool if r["stage"] == s)
                       for s in set(r["stage"] for r in pool)},
            "n_verify_exceptions": sum(1 for r in pool
                                       if r.get("verify_exception")),
        },
        "selection": selection,
        "counts": counts,
        "models": results,
    }
    # strip any leftover heavy keys from pool records (not saved anyway)
    with open(OUT_JSON, "w") as f:
        json.dump(out, f)
    print("wrote", OUT_JSON, flush=True)
    print(json.dumps(counts, indent=1), flush=True)


if __name__ == "__main__":
    main()
