"""Parity of DV2_Vertex (the WCT entry point) with the offline predictions of a SPACE-POINT model (format
dv2max-bundle-v2, e.g. aurora_space), on the DL-vertex dump the model was evaluated on.  validate_r7.py's three
steps; only extract differs (it also stores each call's space points), run passes them, report is validate_r7's
with the bundle's x shift.

  extract  (needs h5py) read the export shards <export>/<prefix>_<sample>.s????.h5 and their .space.h5 files
           (sbnd-vertexer aurora/convert_aurora.py), keep the calls in --pred (key, x, y, z[, score, score_raw]),
           store per call the exact C++ arguments: cloud x, y, z, q, n_vertex_rows, trad (NaN unless trad_valid),
           q_scale, q_offset, and the space points (x, y, z, q float32, cluster_id int32: T_dlvtx_space rows of
           the call's space_cluster_id, as the chain collects them); plus label and offline prediction -> calls.npz
           An event id repeated in a later shard is skipped (the cache kept the first copy).
  run      DV2_Vertex.DV2_Vertex(bundle, cloud bytes ..., space bytes ...) per call -> results.jsonl
  report   validate_r7.report (distance to the offline pick in the label frame, score agreement, < 1 cm vs label
           for both, flips, CPU time) -> report[_<name>].json / .md

usage:
  python validate_space.py extract --out DIR --export <dlvtx>/aurora_export/nocut
         --pred <sbnd-vertexer>/dv2max/runs/aurora_space/ranker_refit/refit_test_fold0.csv [--max-calls N]
  python validate_space.py run --out DIR --bundle <wire-cell-data>/dv2max/sbnd/aurora_space/bundle.json [--workers N]
  python validate_space.py report --out DIR --bundle <bundle.json> [--results NAME.jsonl]
"""
import argparse
import csv
import glob
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import validate_r7 as V  # noqa: E402

SAMPLES = {"nue": "nuecc", "bnb": "bnb"}          # aurora_<dir> shard names, nue first (as the cache)


# ---------------------------------------------------------------------------- extract
def extract(a):
    import h5py
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    pred = {}
    with open(a.pred) as f:
        for r in csv.DictReader(f):
            pred[r["key"]] = r
    cols = {k: [] for k in ("x", "y", "z", "q", "sx", "sy", "sz", "sq", "scid")}
    off, soff = [0], [0]
    meta = {k: [] for k in ("key", "nv", "trad", "qsc", "qof", "lab", "ps", "opred", "oscore", "oraw")}
    seen = set()
    for sample, sd in SAMPLES.items():
        for p in sorted(glob.glob(str(Path(a.export) / f"{a.prefix}_{sd}.s????.h5"))):
            with h5py.File(p, "r") as h, h5py.File(p[:-3] + ".space.h5", "r") as hs:
                c = h["calls"]
                co = c["cloud_off"][()]
                rows = hs["calls/export_row"][()]
                so = hs["calls/space_off"][()]
                skeys = [k.decode() if isinstance(k, bytes) else k for k in hs["calls/key"][()]]
                r_, s_, e_, nu_, pa_ = (c[k][()] for k in ("runNo", "subRunNo", "eventNo", "nu_index", "pass"))
                lab_ok = c["label_ok"][()]
                want = [(j, int(rows[j])) for j, k in enumerate(skeys) if k in pred and k not in seen]
                if not want:
                    continue
                cloud = {k: h["cloud"][k] for k in "xyzq"}
                space = {k: hs["space"][k][()] for k in ("x", "y", "z", "q", "cluster_id")}
                nvr, tv, tx = c["n_vertex_rows"][()], c["trad_valid"][()], c["trad_xyz"][()]
                qs_, qo_, lx = c["q_scale"][()], c["q_offset"][()], c["label_xyz"][()]
                for j, i in want:
                    key = skeys[j]
                    if key != f"{sample}:{r_[i]}:{s_[i]}:{e_[i]}:{nu_[i]}:p{pa_[i]}" or not lab_ok[i]:
                        raise RuntimeError(f"{p}: space call {j} ({key}) does not match export call {i}")
                    seen.add(key)
                    s, e = int(co[i]), int(co[i + 1])
                    for k in "xyzq":
                        cols[k].append(cloud[k][s:e].astype(np.float32))
                    off.append(off[-1] + e - s)
                    a0, b0 = int(so[j]), int(so[j + 1])
                    for k, src in (("sx", "x"), ("sy", "y"), ("sz", "z"), ("sq", "q")):
                        cols[k].append(space[src][a0:b0].astype(np.float32))
                    cols["scid"].append(space["cluster_id"][a0:b0].astype(np.int32))
                    soff.append(soff[-1] + b0 - a0)
                    meta["key"].append(key)
                    meta["nv"].append(int(nvr[i]))
                    meta["trad"].append(tx[i].astype(np.float32) if tv[i] == 1 else np.full(3, np.nan, np.float32))
                    meta["qsc"].append(float(qs_[i]))
                    meta["qof"].append(float(qo_[i]))
                    meta["lab"].append(lx[i].astype(np.float64))
                    meta["ps"].append(int(pa_[i]))
                    pr = pred[key]
                    meta["opred"].append([float(pr["x"]), float(pr["y"]), float(pr["z"])])
                    meta["oscore"].append(float(pr["score"]) if pr.get("score", "") != "" else np.nan)
                    meta["oraw"].append(float(pr["score_raw"]) if pr.get("score_raw", "") != "" else np.nan)
            print(f"[extract] {Path(p).name}: {len(meta['key'])} calls so far", flush=True)
    n = len(meta["key"])
    missing = set(pred) - seen
    if missing:
        raise SystemExit(f"{len(missing)} of {len(pred)} --pred calls not found in the export, e.g. {sorted(missing)[:3]}")
    order = np.arange(n)
    if a.max_calls and a.max_calls < n:
        order = np.sort(np.random.default_rng(a.seed).choice(n, a.max_calls, replace=False))
    lens, slens = np.diff(off)[order], np.diff(soff)[order]

    def cat(k, dt):
        return np.concatenate([cols[k][i] for i in order]).astype(dt) if len(order) else np.zeros(0, dt)
    np.savez(out / "calls.npz", key=np.asarray(meta["key"])[order],
             x=cat("x", np.float32), y=cat("y", np.float32), z=cat("z", np.float32), q=cat("q", np.float32),
             off=np.concatenate([[0], np.cumsum(lens)]).astype(np.int64),
             sx=cat("sx", np.float32), sy=cat("sy", np.float32), sz=cat("sz", np.float32), sq=cat("sq", np.float32),
             scid=cat("scid", np.int32), soff=np.concatenate([[0], np.cumsum(slens)]).astype(np.int64),
             n_vertex_rows=np.asarray(meta["nv"])[order], trad=np.asarray(meta["trad"], np.float32).reshape(-1, 3)[order],
             q_scale=np.asarray(meta["qsc"])[order], q_offset=np.asarray(meta["qof"])[order],
             label=np.asarray(meta["lab"]).reshape(-1, 3)[order], pass_=np.asarray(meta["ps"])[order],
             off_pred=np.asarray(meta["opred"]).reshape(-1, 3)[order], off_score=np.asarray(meta["oscore"])[order],
             off_score_raw=np.asarray(meta["oraw"])[order], off_source=np.full(len(order), ""),
             pred_file=np.asarray(str(Path(a.pred).resolve())))
    print(f"[extract] {n} calls of {a.pred}, kept {len(order)}, {int(slens.sum())} space points -> "
          f"{out / 'calls.npz'}", flush=True)


# ---------------------------------------------------------------------------- run
_W = {}


def _init(bundle):
    os.environ["DV2_THREADS"] = "1"
    sys.path.insert(0, str(V.PYDIR))
    import DV2_Vertex as dv
    _W["dv"], _W["bundle"] = dv, bundle
    t0, c0 = time.perf_counter(), time.process_time()
    dv._vertexer(bundle)
    _W["load"] = (time.perf_counter() - t0, time.process_time() - c0)


def _call(item):
    i, key, x, y, z, q, nv, trad, qs, qo, sx, sy, sz, sq, scid = item
    dv = _W["dv"]
    t0, c0 = time.perf_counter(), time.process_time()
    out = dv.DV2_Vertex(_W["bundle"], x.tobytes(), y.tobytes(), z.tobytes(), q.tobytes(), int(nv), trad.tobytes(),
                        float(qs), float(qo), "float32", sx.tobytes(), sy.tobytes(), sz.tobytes(), sq.tobytes(),
                        scid.tobytes())
    wall, cpu = time.perf_counter() - t0, time.process_time() - c0
    v = np.frombuffer(out, np.float32)
    last = getattr(dv, "_LAST", {}) or {}
    return {"i": int(i), "key": str(key), "out": [float(c) for c in v], "wall_s": wall, "cpu_s": cpu,
            "n": int(len(x)), "n_space": int(len(sq)), "source": last.get("source"),
            "score_raw": last.get("score_raw"), "stage_ms": last.get("stage_ms"), "reason": last.get("reason"),
            "input_sha1": last.get("input_sha1"), "pid": os.getpid(), "load": _W["load"]}


def run(a):
    out = Path(a.out)
    z = V._load(out / "calls.npz")
    n = len(z["key"]) if a.first is None else min(a.first, len(z["key"]))
    pick = np.arange(n)
    if a.sample:
        pick = np.sort(np.random.default_rng(a.seed).choice(n, min(a.sample, n), replace=False))
    res_path = out / (a.results or "results.jsonl")
    done = set()
    if res_path.exists():
        for line in open(res_path):
            done.add(json.loads(line)["key"])
    o, so = z["off"], z["soff"]
    items = [(i, z["key"][i], z["x"][o[i]:o[i + 1]], z["y"][o[i]:o[i + 1]], z["z"][o[i]:o[i + 1]],
              z["q"][o[i]:o[i + 1]], z["n_vertex_rows"][i], z["trad"][i], z["q_scale"][i], z["q_offset"][i],
              z["sx"][so[i]:so[i + 1]], z["sy"][so[i]:so[i + 1]], z["sz"][so[i]:so[i + 1]], z["sq"][so[i]:so[i + 1]],
              z["scid"][so[i]:so[i + 1]])
             for i in pick if z["key"][i] not in done]
    # longest calls first, so the pool does not end on one long straggler
    items.sort(key=lambda t: -(len(t[2]) + len(t[13])))
    print(f"[run] {len(items)} calls to do ({len(done)} done) with {a.workers} workers, bundle {a.bundle}", flush=True)
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(k, "1")
    t0 = time.time()
    with open(res_path, "a") as f:
        if a.workers <= 1:
            _init(a.bundle)
            it, pool = map(_call, items), None
        else:
            import multiprocessing as mp
            pool = mp.get_context("spawn").Pool(a.workers, initializer=_init, initargs=(a.bundle,))
            it = pool.imap_unordered(_call, items, chunksize=1)
        for k, r in enumerate(it):
            f.write(json.dumps(r) + "\n")
            f.flush()
            if (k + 1) % 200 == 0 or k + 1 == len(items):
                print(f"[run] {k + 1}/{len(items)} ({time.time() - t0:.0f} s)", flush=True)
        if pool is not None:
            pool.close()
            pool.join()


# ---------------------------------------------------------------------------- report
def report(a):
    man = json.loads(Path(a.bundle).read_text())
    V.X_SHIFT_CM = float(man["x_shift_cm"])                  # validate_r7.to_label: the bundle's label frame
    V.MODEL_NAME = f"{man['model']} (space points)"
    V.report(a)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=("extract", "run", "report"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--export", default=None)
    ap.add_argument("--prefix", default="aurora")
    ap.add_argument("--pred", default=None)
    ap.add_argument("--max-calls", type=int, default=None)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--bundle", default=None)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--first", type=int, default=None)
    ap.add_argument("--sample", type=int, default=None)
    ap.add_argument("--results", default=None)
    a = ap.parse_args()
    {"extract": extract, "run": run, "report": report}[a.step](a)


if __name__ == "__main__":
    main()
