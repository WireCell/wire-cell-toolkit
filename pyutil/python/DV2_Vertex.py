'''
DV2-MAX neutrino vertex for one DL-vertex call (called by WCPPyUtil::DV2_Vertex through the embedded interpreter).

DV2_Vertex(weights, x, y, z, q, n_vertex_rows, trad, q_scale, q_offset, dtype="float32") -> bytes
  weights        absolute path of the model bundle manifest (e.g. <WIRECELL_PATH>/dv2max/r7_cloud/bundle.json)
  x, y, z, q     bytes of `dtype` arrays: the call's point cloud (cm, raw reco frame; q = dQ * q_scale + q_offset,
                 i.e. the floats the DL-vertex dump writes); the first n_vertex_rows points are the PR-graph
                 vertex rows
  n_vertex_rows  int
  trad           bytes of 3 `dtype` values: the traditional main vertex (cm, raw reco frame), NaN if none
  q_scale, q_offset  float
Returns bytes of 4 `dtype` values [x, y, z, score]: the DV2-MAX top-1 vertex in the RAW reco frame and its calibrated
P(< 1 cm) in [0, 1]; empty bytes when no prediction is made (fewer than 5 points, no candidate, or an internal error,
which is reported on stderr; set DV2_STRICT=1 to raise instead).  The caller then keeps the traditional vertex.

Environment: DV2_THREADS (torch / LightGBM threads, default 1; interop threads are set to 1), DV2_HOOK_LOG=<path>
(append one JSON line per call: n, n_vertex_rows, pred [raw frame], score, score_raw, source, n_cand, stage_ms, ms; the
call that loads the models also carries first_call_ms = import / load breakdown), DV2_DEBUG (print a summary per call to stdout), DV2_VERIFY=0 (skip the sha256 check
of the bundle files at load).  Models are loaded once per process and manifest path.
'''
import time

_T0 = time.perf_counter()
import json  # noqa: E402
import os  # noqa: E402
import sys  # noqa: E402

import numpy as np  # noqa: E402

_CACHE = {}
_THREADS = None
_LAST = {}                 # the last call's record (what DV2_HOOK_LOG writes), for tests


_FIRST = {"import_numpy": round(1e3 * (time.perf_counter() - _T0), 1)}   # load-time breakdown (ms), reported once


def _threads():
    global _THREADS
    if _THREADS is None:
        t0 = time.perf_counter()
        import torch
        _FIRST["import_torch"] = round(1e3 * (time.perf_counter() - t0), 1)
        try:
            n = max(1, int(os.environ.get("DV2_THREADS", "1") or 1))
        except ValueError:
            n = 1
        torch.set_num_threads(n)
        try:
            torch.set_num_interop_threads(1)
        except RuntimeError:            # already set, or inter-op work has started in this process
            pass
        _THREADS = n
    return _THREADS


def _vertexer(weights):
    v = _CACHE.get(weights)
    if v is None:
        n = _threads()
        t0 = time.perf_counter()
        from dv2max_vertexer.vertexer import Vertexer
        t1 = time.perf_counter()
        v = Vertexer(weights, verify=os.environ.get("DV2_VERIFY", "1") != "0", threads=n)
        _FIRST.update(import_dv2max_vertexer=round(1e3 * (t1 - t0), 1),
                      load=round(1e3 * (time.perf_counter() - t1), 1), load_steps=v.load_ms)
        _CACHE[weights] = v
    return v


def _hook(rec):
    path = os.environ.get("DV2_HOOK_LOG")
    if path:
        with open(path, "a") as f:
            f.write(json.dumps(rec) + "\n")


def DV2_Vertex(weights, x, y, z, q, n_vertex_rows, trad, q_scale, q_offset, dtype="float32"):
    global _LAST
    t0 = time.perf_counter()
    dt = np.dtype(dtype)
    X = np.frombuffer(x, dtype=dt)
    Y = np.frombuffer(y, dtype=dt)
    Z = np.frombuffer(z, dtype=dt)
    Q = np.frombuffer(q, dtype=dt)
    T = np.frombuffer(trad, dtype=dt)
    nv = int(n_vertex_rows)
    rec = {"n": int(len(Q)), "n_vertex_rows": nv}
    try:
        if not (len(X) == len(Y) == len(Z) == len(Q)) or len(T) != 3 or not 0 <= nv <= len(Q):
            raise ValueError(f"bad input sizes x {len(X)} y {len(Y)} z {len(Z)} q {len(Q)} trad {len(T)} "
                             f"n_vertex_rows {nv}")
        v = _vertexer(weights)
        from dv2max_vertexer.vertexer import event_from_call, unshift_x
        ev = event_from_call(X, Y, Z, Q, nv, T, float(q_scale), float(q_offset), v.x_shift_cm, v.link_cm)
        r = v.run(ev)
        if r is None:
            out = b""
            rec.update(pred=None, reason="no prediction (too few points or no candidate)")
        else:
            p = unshift_x(r["xyz"], v.x_shift_cm)
            out = np.array([p[0], p[1], p[2], r["score"]], dtype=dt).tobytes()
            rec.update(pred=[float(c) for c in p], score=r["score"], score_raw=r["score_raw"], source=r["source"],
                       n_cand=r["n_cand"], stage_ms=r["ms"])
    except Exception as e:                 # noqa: BLE001 - the chain keeps the traditional vertex
        if os.environ.get("DV2_STRICT", "0") not in ("", "0"):
            raise
        print(f"DV2_Vertex: {type(e).__name__}: {e} (no DV2-MAX vertex for this call)", file=sys.stderr, flush=True)
        out = b""
        rec.update(pred=None, reason=f"error: {type(e).__name__}: {e}")
    rec["ms"] = round(1e3 * (time.perf_counter() - t0), 2)
    if "load" in _FIRST and not _FIRST.get("reported"):
        rec["first_call_ms"] = {k: v for k, v in _FIRST.items() if k != "reported"}
        _FIRST["reported"] = True
    _LAST = rec
    _hook(rec)
    if os.environ.get("DV2_DEBUG"):
        print("DV2_Vertex:", json.dumps(rec), flush=True)
    return out
