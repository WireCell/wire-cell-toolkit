"""Parity of DV2_Vertex (the WCT entry point) with the offline r7_cloud predictions, on the r7 DL-vertex dump.

Three steps (the run step needs only the in-chain runtime: numpy + the DV2_Vertex dependencies):
  extract  (needs h5py) read the r7 export, select the fold-k calls that are in pred_dv2max.csv (passes 0 and 1),
           and store per call the exact C++ arguments (raw cloud floats x, y, z, q, n_vertex_rows, trad = NaN unless
           trad_valid, q_scale, q_offset), the label and the offline prediction -> calls.npz
  run      call DV2_Vertex.DV2_Vertex(bundle, bytes...) for every call (workers in parallel, DV2_THREADS=1 each),
           record the returned vertex, the CPU / wall time of the call and the stage times -> results.jsonl
  report   compare with the offline picks: distance to the offline pick (in the label frame: x -= sign(x) * 0.31),
           score agreement, < 1 cm efficiency vs label for both (same calls, per pass), top-1 flips, CPU time
           -> report[_<name>].json / .md (named after the results file)
A call without a DV2_Vertex prediction (empty bytes) is scored as the chain does it: the traditional vertex (moved to
the label frame), i.e. the offline tag_nu fallback.

usage:
  python validate_r7.py extract --out DIR [--fold 0] [--max-calls N] [--seed 0]
  python validate_r7.py run --out DIR --bundle <W>/cfg/dv2max/r7_cloud/bundle.json [--workers 16] [--first N]
                           [--sample N --seed S] [--results NAME.jsonl]
  python validate_r7.py report --out DIR [--results NAME.jsonl]
"""
import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PYDIR = HERE.parent                                  # pyutil/python: DV2_Vertex.py and dv2max_vertexer/
EXPORT = "/lus/eagle/projects/neutrinoGPU/snehadri/dlvtx/r7_export"
PRED = "/lus/eagle/projects/neutrinoGPU/snehadri/sbnd-vertexer/dv2max/runs/r7_cloud/pred_dv2max.csv"
INPUTS = {"nue": ["r7_mc-nuecc.h5"], "bnb": ["r7_mc-cv.h5", "r7_mc-cv-fix.h5"]}
X_SHIFT_CM = 0.31


def to_label(p):
    p = np.array(p, np.float64, copy=True).reshape(-1, 3)
    p[:, 0] -= np.sign(p[:, 0]) * X_SHIFT_CM
    return p


# ---------------------------------------------------------------------------- extract
def extract(a):
    import h5py
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    pred = {}
    with open(a.pred) as f:
        for r in csv.DictReader(f):
            if int(r["fold"]) == a.fold:
                pred[r["key"]] = r
    keys, xs, ys, zs, qs, off = [], [], [], [], [], [0]
    nv, trad, qsc, qof, lab, ps, opred, oscore, oraw, osrc, ousrc = [], [], [], [], [], [], [], [], [], [], []
    for sample, files in INPUTS.items():
        for fn in files:
            with h5py.File(Path(a.export) / fn, "r") as h:
                c = h["calls"]
                ok = c["label_ok"][()] & np.isin(c["pass"][()], (0, 1))
                r_, s_, e_, nu_, pa_ = (c[k][()] for k in ("runNo", "subRunNo", "eventNo", "nu_index", "pass"))
                co = c["cloud_off"][()]
                for i in np.flatnonzero(ok):
                    key = f"{sample}:{r_[i]}:{s_[i]}:{e_[i]}:{nu_[i]}:p{pa_[i]}"
                    if key not in pred:
                        continue
                    s, e = int(co[i]), int(co[i + 1])
                    for lst, k in ((xs, "x"), (ys, "y"), (zs, "z"), (qs, "q")):
                        lst.append(h["cloud"][k][s:e].astype(np.float32))
                    off.append(off[-1] + e - s)
                    keys.append(key)
                    nv.append(int(c["n_vertex_rows"][i]))
                    trad.append(c["trad_xyz"][i].astype(np.float32) if c["trad_valid"][i] == 1
                                else np.full(3, np.nan, np.float32))
                    qsc.append(float(c["q_scale"][i]))
                    qof.append(float(c["q_offset"][i]))
                    lab.append(c["label_xyz"][i].astype(np.float64))
                    ps.append(int(pa_[i]))
                    p = pred[key]
                    opred.append([float(p["x"]), float(p["y"]), float(p["z"])])
                    oscore.append(float(p["score"]) if p["score"] != "" else np.nan)
                    oraw.append(float(p["score_raw"]) if p.get("score_raw", "") != "" else np.nan)
                    osrc.append(p.get("pred_source", ""))
                    ousrc.append(p.get("source", ""))
    n = len(keys)
    order = np.arange(n)
    if a.max_calls and a.max_calls < n:
        order = np.sort(np.random.default_rng(a.seed).choice(n, a.max_calls, replace=False))
    cat = lambda L: np.concatenate([L[i] for i in order]) if len(order) else np.zeros(0, np.float32)  # noqa: E731
    lens = np.diff(off)[order]
    np.savez(out / "calls.npz", key=np.asarray(keys)[order], x=cat(xs), y=cat(ys), z=cat(zs), q=cat(qs),
             off=np.concatenate([[0], np.cumsum(lens)]).astype(np.int64), n_vertex_rows=np.asarray(nv)[order],
             trad=np.asarray(trad, np.float32).reshape(-1, 3)[order], q_scale=np.asarray(qsc)[order],
             q_offset=np.asarray(qof)[order], label=np.asarray(lab).reshape(-1, 3)[order],
             pass_=np.asarray(ps)[order], off_pred=np.asarray(opred).reshape(-1, 3)[order],
             off_score=np.asarray(oscore)[order], off_score_raw=np.asarray(oraw)[order],
             off_source=np.asarray(osrc)[order], off_union_source=np.asarray(ousrc)[order])
    print(f"[extract] fold {a.fold}: {n} calls in {a.pred}, kept {len(order)} -> {out / 'calls.npz'}", flush=True)


# ---------------------------------------------------------------------------- run
def _load(path):
    """calls.npz as a dict of arrays (an NpzFile re-reads an array on every item access)."""
    with np.load(path) as f:
        return {k: f[k] for k in f.files}


_W = {}


def _init(bundle):
    os.environ["DV2_THREADS"] = "1"
    sys.path.insert(0, str(PYDIR))
    import DV2_Vertex as dv
    _W["dv"], _W["bundle"] = dv, bundle
    t0, c0 = time.perf_counter(), time.process_time()
    dv._vertexer(bundle)                             # load once, outside the per-call timing
    _W["load"] = (time.perf_counter() - t0, time.process_time() - c0)


def _call(item):
    i, key, x, y, z, q, nv, trad, qs, qo = item
    dv = _W["dv"]
    t0, c0 = time.perf_counter(), time.process_time()
    out = dv.DV2_Vertex(_W["bundle"], x.tobytes(), y.tobytes(), z.tobytes(), q.tobytes(), int(nv), trad.tobytes(),
                        float(qs), float(qo), "float32")
    wall, cpu = time.perf_counter() - t0, time.process_time() - c0
    v = np.frombuffer(out, np.float32)
    last = getattr(dv, "_LAST", {}) or {}
    return {"i": int(i), "key": str(key), "out": [float(c) for c in v], "wall_s": wall, "cpu_s": cpu,
            "n": int(len(x)), "source": last.get("source"), "score_raw": last.get("score_raw"),
            "stage_ms": last.get("stage_ms"), "reason": last.get("reason"), "pid": os.getpid(), "load": _W["load"]}


def run(a):
    out = Path(a.out)
    z = _load(out / "calls.npz")
    n = len(z["key"]) if a.first is None else min(a.first, len(z["key"]))
    pick = np.arange(n)
    if a.sample:                                     # a random subset (e.g. for a clean 1-worker timing run)
        pick = np.sort(np.random.default_rng(a.seed).choice(n, min(a.sample, n), replace=False))
    res_path = out / (a.results or "results.jsonl")
    done = set()
    if res_path.exists():
        for line in open(res_path):
            done.add(json.loads(line)["key"])
    off = z["off"]
    items = [(i, z["key"][i], z["x"][off[i]:off[i + 1]], z["y"][off[i]:off[i + 1]], z["z"][off[i]:off[i + 1]],
              z["q"][off[i]:off[i + 1]], z["n_vertex_rows"][i], z["trad"][i], z["q_scale"][i], z["q_offset"][i])
             for i in pick if z["key"][i] not in done]
    print(f"[run] {len(items)} calls to do ({len(done)} done) with {a.workers} workers, bundle {a.bundle}", flush=True)
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    t0 = time.time()
    with open(res_path, "a") as f:
        if a.workers <= 1:
            _init(a.bundle)
            it = map(_call, items)
            pool = None
        else:
            import multiprocessing as mp
            pool = mp.get_context("spawn").Pool(a.workers, initializer=_init, initargs=(a.bundle,))
            it = pool.imap_unordered(_call, items, chunksize=1)
        for k, r in enumerate(it):
            f.write(json.dumps(r) + "\n")
            f.flush()
            if (k + 1) % 50 == 0 or k + 1 == len(items):
                print(f"[run] {k + 1}/{len(items)} ({time.time() - t0:.0f} s)", flush=True)
        if pool is not None:
            pool.close()
            pool.join()


# ---------------------------------------------------------------------------- report
def _q(a, qs=(0.5, 0.9, 0.99, 1.0)):
    a = np.asarray(a, np.float64)
    a = a[np.isfinite(a)]
    return {f"p{int(100 * q)}" if q < 1 else "max": (round(float(np.quantile(a, q)), 6) if len(a) else None) for q in qs}


def report(a):
    out = Path(a.out)
    z = _load(out / "calls.npz")
    pos = {k: i for i, k in enumerate(z["key"])}
    R = [json.loads(line) for line in open(out / (a.results or "results.jsonl"))]
    R = {r["key"]: r for r in R}
    idx = np.array(sorted(pos[k] for k in R))
    keys = z["key"][idx]
    mine = np.full((len(idx), 3), np.nan)
    score = np.full(len(idx), np.nan)
    empty = np.zeros(len(idx), bool)
    for j, k in enumerate(keys):
        o = R[k]["out"]
        if len(o) == 4:
            mine[j] = to_label(o[:3])[0]
            score[j] = o[3]
        else:
            empty[j] = True
            mine[j] = to_label(z["trad"][idx[j]].astype(np.float64))[0]     # the chain keeps the traditional vertex
    offp = z["off_pred"][idx]
    lab = z["label"][idx]
    ps = z["pass_"][idx]
    d = np.linalg.norm(mine - offp, axis=1)
    d = np.where(np.isfinite(d), d, np.where(np.isnan(mine).any(1) & np.isnan(offp).any(1), 0.0, np.inf))
    dl_m = np.linalg.norm(mine - lab, axis=1)
    dl_o = np.linalg.norm(offp - lab, axis=1)
    dl_m = np.where(np.isfinite(dl_m), dl_m, np.inf)
    dl_o = np.where(np.isfinite(dl_o), dl_o, np.inf)
    osrc = z["off_source"][idx]
    rep = {"n_calls": int(len(idx)), "n_empty": int(empty.sum()),
           "n_offline_fallback_tag_nu": int((osrc == "tag_nu").sum()),
           "empty_vs_offline_fallback_same": int((empty & (osrc == "tag_nu")).sum())}
    rep["dist_to_offline_pick"] = {f"frac_lt_{r:g}": round(float(np.mean(d < r)), 5) for r in (0.01, 0.05, 0.1, 1.0)}
    rep["dist_to_offline_pick"].update(_q(d))
    rep["n_dist_gt_0.1"] = int((d > 0.1).sum())
    ok = ~empty & np.isfinite(z["off_score"][idx])
    ds = np.abs(score - z["off_score"][idx])[ok]
    rep["score_abs_diff"] = _q(ds)
    raw = np.array([R[k].get("score_raw") if R[k].get("score_raw") is not None else np.nan for k in keys], float)
    rep["score_raw_abs_diff"] = _q(np.abs(raw - z["off_score_raw"][idx])[ok])
    src = np.array([str(R[k].get("source")) for k in keys])
    if "off_union_source" in z:                       # candidate source of the pick (cand / tag / clusend / img / cnn)
        lm = osrc == "dv2max_lambdamart"
        rep["pick_source_same_frac"] = round(float(np.mean(src[lm] == z["off_union_source"][idx][lm])), 5)
    eff = {}
    for name, m in (("all", np.ones(len(idx), bool)), ("pass0", ps == 0), ("pass1", ps == 1)):
        hm, ho = dl_m[m] < 1.0, dl_o[m] < 1.0
        eff[name] = {"n": int(m.sum()), "dv2_vertex": round(float(hm.mean()), 5) if m.any() else None,
                     "offline": round(float(ho.mean()), 5) if m.any() else None,
                     "only_dv2_vertex": int((hm & ~ho).sum()), "only_offline": int((~hm & ho).sum())}
        for r in (0.5, 2.0, 5.0):
            eff[name][f"lt{r:g}_dv2_vertex"] = round(float(np.mean(dl_m[m] < r)), 5) if m.any() else None
            eff[name][f"lt{r:g}_offline"] = round(float(np.mean(dl_o[m] < r)), 5) if m.any() else None
    rep["eff_vs_label"] = eff
    big = d > 0.1
    rep["flips_gt_0.1cm"] = {"n": int(big.sum()), "frac": round(float(big.mean()), 5),
                             "hit_changes": int(((dl_m < 1) != (dl_o < 1))[big].sum()),
                             "dist_quantiles": _q(d[big]) if big.any() else None}
    cpu = np.array([R[k]["cpu_s"] for k in keys])
    wall = np.array([R[k]["wall_s"] for k in keys])
    rep["cpu_s_per_call"] = {"median": round(float(np.median(cpu)), 3), "p90": round(float(np.quantile(cpu, 0.9)), 3),
                             "max": round(float(cpu.max()), 3), "mean": round(float(cpu.mean()), 3)}
    rep["wall_s_per_call"] = {"median": round(float(np.median(wall)), 3), "p90": round(float(np.quantile(wall, 0.9)), 3),
                              "max": round(float(wall.max()), 3)}
    st = [R[k]["stage_ms"] for k in keys if R[k].get("stage_ms")]
    if st:
        rep["stage_ms_median"] = {s: round(float(np.median([x[s] for x in st])), 1) for s in st[0]}
    loads = {R[k]["pid"]: R[k]["load"] for k in keys}
    rep["model_load_s"] = {"n_workers": len(loads), "wall_median": round(float(np.median([v[0] for v in loads.values()])), 2)}
    npts = np.array([R[k]["n"] for k in keys])
    rep["n_points"] = _q(npts)
    worst = np.argsort(-np.where(np.isfinite(d), d, 1e9))[:15]
    rep["largest_disagreements"] = [{"key": str(keys[j]), "d_cm": round(float(d[j]), 4),
                                     "d_label_dv2": round(float(dl_m[j]), 3), "d_label_off": round(float(dl_o[j]), 3),
                                     "score": None if not np.isfinite(score[j]) else round(float(score[j]), 4),
                                     "off_score": round(float(z["off_score"][idx[j]]), 4), "source": src[j],
                                     "off_source": str(z["off_union_source"][idx[j]] if "off_union_source" in z
                                                       else z["off_source"][idx[j]])} for j in worst if d[j] > 0.01]
    stem = Path(a.results or "results.jsonl").stem.replace("results", "report") or "report"
    if not stem.startswith("report"):
        stem = "report_" + stem
    (out / f"{stem}.json").write_text(json.dumps(rep, indent=1))
    e = rep["eff_vs_label"]
    lines = [f"# DV2_Vertex vs offline r7_cloud ({rep['n_calls']} calls, {rep['n_empty']} without prediction)", "",
             "distance to the offline pick: " + ", ".join(f"{k} {v}" for k, v in rep["dist_to_offline_pick"].items()),
             f"score |diff|: {rep['score_abs_diff']}", "",
             "| calls | n | eff<1cm DV2_Vertex | eff<1cm offline | only DV2_Vertex | only offline |", "|---|---|---|---|---|---|"]
    for name in ("pass0", "pass1", "all"):
        lines.append(f"| {name} | {e[name]['n']} | {e[name]['dv2_vertex']} | {e[name]['offline']} | "
                     f"{e[name]['only_dv2_vertex']} | {e[name]['only_offline']} |")
    lines += ["", f"top-1 flips (> 0.1 cm): {rep['flips_gt_0.1cm']}",
              f"CPU s per call (1 thread): {rep['cpu_s_per_call']}; wall {rep['wall_s_per_call']}",
              f"stage ms (median): {rep.get('stage_ms_median')}"]
    (out / f"{stem}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(json.dumps(rep, indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=("extract", "run", "report"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--export", default=EXPORT)
    ap.add_argument("--pred", default=PRED)
    ap.add_argument("--fold", type=int, default=0)
    ap.add_argument("--max-calls", type=int, default=None)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--bundle", default=None)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--first", type=int, default=None, help="run only the first N calls of calls.npz")
    ap.add_argument("--sample", type=int, default=None, help="run a random subset of N calls (--seed)")
    ap.add_argument("--results", default=None, help="results file name inside --out (default results.jsonl)")
    a = ap.parse_args()
    {"extract": extract, "run": run, "report": report}[a.step](a)


if __name__ == "__main__":
    main()
