#!/usr/bin/env python
"""Write data/dv2max_space_golden.json: frozen hold-out calls of a SPACE-POINT bundle (aurora_space) and what the
current code returns for them.  test_dv2max_space.py (TestSpaceGolden) replays them.  Run in the DV2-MAX runtime env
(numpy, scipy, torch, lightgbm) with the bundle:

  python make_dv2max_space_golden.py --calls <DIR>/calls.npz --bundle <wire-cell-data>/dv2max/sbnd/aurora_space/bundle.json \\
         --out data/dv2max_space_golden.json [--workers 16]

calls.npz is dv2max_vertexer/validate_space.py's extract of the aurora_space validation hold-out (calls the bundle was
NOT trained on).  Selection (deterministic, seed 0), small calls so the fixture stays small:
  - cloud + space points <= --max-points, stratified sample (bnb, nue) x pass (0, 1) x (with / without space points),
    plus calls whose traditional vertex is missing;
  - a call qualifies only if its top-1 ranker score beats every candidate more than 0.1 cm away by at least --margin
    (raw LambdaMART score), so machine / library differences cannot flip the pick.
Arrays are stored as base64 of their little-endian bytes (float32; cluster ids int32), so the replay is bit-exact.
"""
import argparse
import base64
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
_W = {}


def _init(bundle):
    os.environ["DV2_THREADS"] = "1"
    sys.path.insert(0, str(PYDIR))
    import DV2_Vertex as dv
    _W["dv"], _W["bundle"], _W["v"] = dv, bundle, dv._vertexer(bundle)


def _eval(item):
    i, x, y, z, q, nv, trad, qs, qo, sx, sy, sz, sq, sc = item
    dv, v = _W["dv"], _W["v"]
    from dv2max_vertexer.vertexer import event_from_call
    ev = event_from_call(x, y, z, q, nv, trad, qs, qo, v.x_shift_cm, v.link_cm, space=(sx, sy, sz, sq, sc),
                         space_cfg=v.space)
    r = v.run(ev, with_records=True)
    margin = None
    if r is not None:
        s, P = np.asarray(r["s"], np.float64), np.asarray(r["rows"], np.float64)
        top = int(np.argmax(s))
        far = np.linalg.norm(P - P[top], axis=1) > 0.1
        margin = float(s[top] - s[far].max()) if far.any() else float("inf")
    out = dv.DV2_Vertex(_W["bundle"], x.tobytes(), y.tobytes(), z.tobytes(), q.tobytes(), int(nv), trad.tobytes(),
                        float(qs), float(qo), "float32", sx.tobytes(), sy.tobytes(), sz.tobytes(), sq.tobytes(),
                        sc.tobytes())
    last = dv._LAST
    return {"i": int(i), "out": [float(c) for c in np.frombuffer(out, np.float32)], "margin": margin,
            "score_raw": last.get("score_raw"), "n_cand": last.get("n_cand"), "source": last.get("source")}


def b64(a):
    return base64.b64encode(np.ascontiguousarray(a).tobytes()).decode()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True)
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--out", default=str(HERE / "data" / "dv2max_space_golden.json"))
    ap.add_argument("--per-cell", type=int, default=2)
    ap.add_argument("--n-nantrad", type=int, default=2)
    ap.add_argument("--max-points", type=int, default=2500)
    ap.add_argument("--margin", type=float, default=0.5)
    ap.add_argument("--workers", type=int, default=1)
    a = ap.parse_args()
    with np.load(a.calls) as f:
        z = {k: f[k] for k in f.files}
    o, so = z["off"], z["soff"]
    n = len(z["key"])
    npts, nsp = np.diff(o), np.diff(so)
    small = (npts + nsp <= a.max_points) & (npts >= 20)
    rng = np.random.default_rng(0)
    samp = np.array([k.split(":")[0] for k in z["key"]])
    nantrad = np.isnan(z["trad"]).any(1)
    cells = []
    for s in ("bnb", "nue"):
        for p in (0, 1):
            for has in (True, False):
                m = small & (samp == s) & (z["pass_"] == p) & ((nsp > 0) == has) & ~nantrad
                cells.append(((s, p, "space" if has else "no_space"), np.flatnonzero(m)))
    cells.append((("any", "any", "nan_trad"), np.flatnonzero(small & nantrad)))
    cand = []
    for name, idx in cells:
        cand.append((name, rng.permutation(idx)[: 8 * max(a.per_cell, a.n_nantrad)]))
    items = []
    for _, idx in cand:
        for i in idx:
            items.append((int(i), z["x"][o[i]:o[i + 1]], z["y"][o[i]:o[i + 1]], z["z"][o[i]:o[i + 1]],
                          z["q"][o[i]:o[i + 1]], int(z["n_vertex_rows"][i]), z["trad"][i], float(z["q_scale"][i]),
                          float(z["q_offset"][i]), z["sx"][so[i]:so[i + 1]], z["sy"][so[i]:so[i + 1]],
                          z["sz"][so[i]:so[i + 1]], z["sq"][so[i]:so[i + 1]], z["scid"][so[i]:so[i + 1]]))
    if a.workers <= 1:
        _init(a.bundle)
        res = list(map(_eval, items))
    else:
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(a.workers, initializer=_init, initargs=(a.bundle,)) as pool:
            res = pool.map(_eval, items, chunksize=1)
    R = {r["i"]: r for r in res}
    calls = []
    for name, idx in cand:
        want = a.n_nantrad if name[2] == "nan_trad" else a.per_cell
        got = 0
        for i in idx:
            r = R[int(i)]
            if r["margin"] is None or r["margin"] < a.margin or len(r["out"]) != 4:
                continue
            calls.append({"key": str(z["key"][i]), "cell": list(map(str, name)), "n_points": int(npts[i]),
                          "n_space": int(nsp[i]), "n_vertex_rows": int(z["n_vertex_rows"][i]),
                          "q_scale": float(z["q_scale"][i]), "q_offset": float(z["q_offset"][i]),
                          "trad": b64(z["trad"][i].astype(np.float32)),
                          **{k: b64(z[k][o[i]:o[i + 1]].astype(np.float32)) for k in ("x", "y", "z", "q")},
                          **{k: b64(z[k][so[i]:so[i + 1]].astype(np.float32)) for k in ("sx", "sy", "sz", "sq")},
                          "scid": b64(z["scid"][so[i]:so[i + 1]].astype(np.int32)),
                          "out": r["out"], "score_raw": r["score_raw"], "n_cand": r["n_cand"],
                          "source": r["source"], "margin": r["margin"],
                          "label": [float(c) for c in z["label"][i]], "offline_pick": [float(c) for c in z["off_pred"][i]]})
            got += 1
            if got == want:
                break
    sha = hashlib.sha256(Path(a.bundle).read_bytes()).hexdigest()
    import lightgbm
    import scipy
    import torch
    G = {"bundle": "dv2max/sbnd/aurora_space/bundle.json", "bundle_sha256": sha,
         "made": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
         "made_by": "pyutil/test/make_dv2max_space_golden.py", "source_calls": str(Path(a.calls).resolve()),
         "selection": vars(a), "tolerance_cm": 1e-3, "score_tolerance": 1e-3,
         "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
                      "torch": torch.__version__, "lightgbm": lightgbm.__version__},
         "encoding": "base64 of little-endian bytes: float32 for x y z q sx sy sz sq trad, int32 for scid",
         "calls": calls}
    Path(a.out).write_text(json.dumps(G, indent=0) + "\n")
    print(f"{len(calls)} calls ({sum(c['n_space'] > 0 for c in calls)} with space points) -> {a.out} "
          f"({Path(a.out).stat().st_size / 1e3:.0f} kB)")


if __name__ == "__main__":
    main()
