#!/usr/bin/env python
"""Write data/dv2max_golden.json: frozen held-out DV2-MAX calls and what the current code returns for them.

The golden tests (test_dv2max.py, doctest_dv2_vertex.cxx) replay these calls and require the same answer.  Run in
the DV2-MAX runtime env (numpy, scipy, torch, lightgbm) with the bundle:

  python make_dv2max_golden.py --calls <DIR>/calls.npz --bundle <wire-cell-data>/dv2max/sbnd/r7_cloud/bundle.json \\
         --out data/dv2max_golden.json [--workers 16]

calls.npz is the fold-0 extract of dv2max_vertexer/validate_r7.py (calls the r7_cloud bundle was NOT trained on).

Selection (deterministic, seed 0), so that the golden answers are robust on other machines:
  - real calls with 20-1500 points, stratified: sample (bnb, nue) x pass (0, 1) x cloud-size tercile, plus calls
    whose traditional vertex is missing (NaN);
  - a call qualifies only if its top-1 ranker score beats every candidate more than 0.1 cm away by at least
    --margin (raw LambdaMART score): floating-point differences between machines / library versions then cannot
    flip the pick;
  - three synthetic calls with fewer than data.min_points (5) points: no prediction.
Per call it stores the exact DV2_Vertex arguments (float32 values; a missing traditional vertex is null) and the returned [x, y, z, score] (raw reco frame,
empty = no prediction), with score_raw, n_cand, the margin and, for information, the label and the offline pick.
"""
import argparse
import datetime
import hashlib
import json
import os
import platform
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PYDIR = HERE.parent / "python"
SAMPLES, PASSES, PER_CELL, N_NANTRAD = ("bnb", "nue"), (0, 1), 6, 2

_W = {}


def _init(bundle):
    os.environ["DV2_THREADS"] = "1"
    sys.path.insert(0, str(PYDIR))
    import DV2_Vertex as dv
    _W["dv"], _W["bundle"], _W["v"] = dv, bundle, dv._vertexer(bundle)


def _eval(item):
    """Run one call through Vertexer.run (records) and DV2_Vertex (the public answer)."""
    i, x, y, z, q, nv, trad, qs, qo = item
    dv, v = _W["dv"], _W["v"]
    from dv2max_vertexer.vertexer import event_from_call
    ev = event_from_call(x, y, z, q, nv, trad, qs, qo, v.x_shift_cm, v.link_cm)
    r = v.run(ev, with_records=True)
    margin = None
    if r is not None:
        s, P = np.asarray(r["s"], np.float64), np.asarray(r["rows"], np.float64)
        top = int(np.argmax(s))
        far = np.linalg.norm(P - P[top], axis=1) > 0.1
        margin = float(s[top] - s[far].max()) if far.any() else float("inf")
    out = dv.DV2_Vertex(_W["bundle"], x.tobytes(), y.tobytes(), z.tobytes(), q.tobytes(), int(nv), trad.tobytes(),
                        float(qs), float(qo), "float32")
    return {"i": int(i), "margin": margin, "expect": [float(c) for c in np.frombuffer(out, np.float32)],
            "score_raw": None if r is None else float(r["score_raw"]), "n_cand": None if r is None else int(r["n_cand"]),
            "source": None if r is None else str(r["source"])}


def f32num(c):
    """A float32 value as the shortest decimal that reads back (via float64) as the same float32 -- half the JSON
    size of repr(float(c))."""
    c = np.float32(c)
    s = np.format_float_positional(c, unique=True, trim="-")
    return float(s) if np.float32(float(s)) == c else float(c)


def f32list(a):
    return [f32num(c) for c in np.asarray(a, np.float32).ravel()]


def f32list_nan(a):
    """As f32list, NaN -> null (strict JSON; the readers turn null back into NaN)."""
    return [None if np.isnan(c) else f32num(c) for c in np.asarray(a, np.float32).ravel()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True)
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--out", default=str(HERE / "data" / "dv2max_golden.json"))
    ap.add_argument("--margin", type=float, default=0.5)
    ap.add_argument("--pool", type=int, default=24, help="candidates evaluated per stratum")
    ap.add_argument("--workers", type=int, default=16)
    a = ap.parse_args()

    with np.load(a.calls) as f:
        z = {k: f[k] for k in f.files}
    off, n = z["off"], np.diff(z["off"])
    sample = np.array([k.split(":")[0] for k in z["key"]])
    nantrad = np.isnan(z["trad"]).any(1)
    ok = (n >= 20) & (n <= 1500)
    edges = np.quantile(n[ok], [1 / 3, 2 / 3])
    size = np.digitize(n, edges)                     # 0 / 1 / 2 = cloud-size tercile
    rng = np.random.default_rng(0)
    strata = {}
    for s in SAMPLES:
        for p in PASSES:
            for t in range(3):
                idx = np.flatnonzero(ok & ~nantrad & (sample == s) & (z["pass_"] == p) & (size == t))
                strata[(s, p, t)] = list(rng.permutation(idx)[:a.pool])
    strata[("nantrad",)] = list(rng.permutation(np.flatnonzero(ok & nantrad)))
    pool = sorted({int(i) for v in strata.values() for i in v})
    items = [(i, z["x"][off[i]:off[i + 1]], z["y"][off[i]:off[i + 1]], z["z"][off[i]:off[i + 1]],
              z["q"][off[i]:off[i + 1]], z["n_vertex_rows"][i], z["trad"][i], z["q_scale"][i], z["q_offset"][i])
             for i in pool]
    print(f"[golden] evaluating {len(items)} candidate calls with {a.workers} workers", flush=True)
    import multiprocessing as mp
    with mp.get_context("spawn").Pool(a.workers, initializer=_init, initargs=(a.bundle,)) as pl:
        res = {r["i"]: r for r in pl.imap_unordered(_eval, items, chunksize=1)}

    def good(i):
        r = res[i]
        return r["margin"] is not None and r["margin"] >= a.margin and len(r["expect"]) == 4

    picked = []
    for key, idx in strata.items():
        want = N_NANTRAD if key == ("nantrad",) else PER_CELL // 3
        picked += [int(i) for i in idx if good(i)][:want]
    margins = np.array([res[i]["margin"] for i in pool if res[i]["margin"] is not None])
    print(f"[golden] margin over the pool: median {np.median(margins):.2f}, >= {a.margin}: "
          f"{(margins >= a.margin).mean():.2f}; picked {len(picked)}", flush=True)

    calls = []
    for i in sorted(picked):
        r = res[i]
        calls.append({"key": str(z["key"][i]), "x": f32list(z["x"][off[i]:off[i + 1]]),
                      "y": f32list(z["y"][off[i]:off[i + 1]]), "z": f32list(z["z"][off[i]:off[i + 1]]),
                      "q": f32list(z["q"][off[i]:off[i + 1]]), "n_vertex_rows": int(z["n_vertex_rows"][i]),
                      "trad": f32list_nan(z["trad"][i]), "q_scale": float(z["q_scale"][i]),
                      "q_offset": float(z["q_offset"][i]), "expect": f32list(r["expect"]), "score_raw": r["score_raw"],
                      "n_cand": r["n_cand"], "source": r["source"], "margin": r["margin"],
                      "label": f32list(z["label"][i]), "offline_pick": f32list(z["off_pred"][i])})
    # synthetic: fewer than data.min_points points -> no prediction
    for npt in (0, 3, 4):
        g = np.random.default_rng(100 + npt)
        p = g.uniform([-150, -150, 50], [150, 150, 450], (npt, 3)).astype(np.float32)
        calls.append({"key": f"synthetic:{npt}_points", "x": f32list(p[:, 0]), "y": f32list(p[:, 1]),
                      "z": f32list(p[:, 2]), "q": f32list(np.full(npt, 2000.0)), "n_vertex_rows": min(npt, 1),
                      "trad": f32list_nan([np.nan] * 3), "q_scale": 0.1, "q_offset": -1000.0, "expect": [],
                      "score_raw": None, "n_cand": None, "source": None, "margin": None})

    import lightgbm
    import scipy
    import torch
    meta = {"generated": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
            "generator": "pyutil/test/make_dv2max_golden.py", "calls_file": str(a.calls),
            "bundle_sha256": hashlib.sha256(Path(a.bundle).read_bytes()).hexdigest(),
            "selection": {"points": [20, 1500], "per_stratum": PER_CELL // 3, "nan_trad": N_NANTRAD,
                          "margin_min": a.margin, "seed": 0, "synthetic_no_prediction": [0, 3, 4]},
            "env": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
                    "torch": torch.__version__, "lightgbm": lightgbm.__version__, "machine": platform.machine()},
            "tolerance": {"position_cm": 0.01, "score": 1e-3},
            "note": "expect = DV2_Vertex output [x, y, z, P(<1 cm)] in the raw reco frame (cm); empty = no prediction"}
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps({"meta": meta, "calls": calls}, allow_nan=False) + "\n")
    print(f"[golden] wrote {a.out}: {len(calls)} calls, {Path(a.out).stat().st_size / 1e3:.0f} kB", flush=True)


if __name__ == "__main__":
    main()
