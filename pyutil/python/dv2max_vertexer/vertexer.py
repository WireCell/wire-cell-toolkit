"""DV2-MAX r7_cloud inference for ONE DL-vertex call (one neutrino candidate's cloud), CPU, torch backend.

Pipeline (= the offline r7_cloud run for a test event of the deployed family; see README.md):
  1. event_from_call   the call's cloud -> Event exactly as r7/prep_r7_cache.py builds a cache row: q decoded to
                       dQ = (q - q_offset) / q_scale clipped at 0, 3 cm connected components as cluster ids, the
                       first n_vertex_rows points as PR candidates, the traditional vertex as tag_nu, every reco
                       input moved x_shift_cm toward the cathode (x -= sign(x) * x_shift_cm)
  2. features          features.compute_features, then the cache's storage precision (point-level features float16,
                       the others float32), because the models and the ranker were fed cached features
  3. pass 1            whole-event ensemble (views x shifts) -> pooled pass-1 record (modes, heat, votes)
  4. pass 2            crops framed at the pass-1 modes (views x shifts) -> pass-2 record
  5. ranker            candidates + features of the deployed family's variant (ranker.event_union), LightGBM
                       LambdaMART score, top-1 after NMS, isotonic map of the top-1 score -> P(< 1 cm)
The returned position is in the LABEL frame (= the shifted reco frame the models work in); DV2_Vertex.py moves it
back to the raw reco frame.
"""
import contextlib
import hashlib
import io
import json
import time
from pathlib import Path

import numpy as np
import torch

from . import baselines, infer, inputs, ranker
from .event import Event
from .features import FEATURE_SPEC, compute_features
from .model import DV2Net

FORMAT = "dv2max-bundle-v1"
F16_MAX = float(np.finfo(np.float16).max)


# ---------------------------------------------------------------------------- event (r7/prep_r7_cache.py)
def xshift(a, x_shift_cm):
    """r7/prep_r7_cache.xshift: move x toward the cathode by x_shift_cm (float32 arithmetic)."""
    a = np.array(a, np.float32, copy=True).reshape(-1, 3)
    a[:, 0] -= np.sign(a[:, 0]) * np.float32(x_shift_cm)
    return a


def unshift_x(xyz, x_shift_cm):
    """Inverse of xshift for a label-frame position (x += sign(x) * x_shift_cm)."""
    p = np.asarray(xyz, np.float64).reshape(3).copy()
    p[0] += np.sign(p[0]) * x_shift_cm
    return p


def cluster_ids(xyz, link_cm):
    """r7/prep_r7_cache.cluster_ids: connected components at link_cm, renumbered by size (largest = 0)."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree
    n = len(xyz)
    if n == 0:
        return np.zeros(0, np.int32)
    pairs = cKDTree(xyz).query_pairs(link_cm, output_type="ndarray")
    g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    _, lab = connected_components(g, directed=False)
    order = np.argsort(-np.bincount(lab), kind="stable")
    remap = np.empty_like(order)
    remap[order] = np.arange(len(order))
    return remap[lab].astype(np.int32)


def event_from_call(x, y, z, q, n_vertex_rows, trad, q_scale, q_offset, x_shift_cm=0.31, link_cm=3.0, key="call"):
    """Truth-free Event of one DL-vertex call, as r7/prep_r7_cache.build_event (reco inputs x-shifted).
    x, y, z, q: float32 arrays of the cloud (q encoded, dQ * q_scale + q_offset); trad: (3,) traditional vertex
    (NaN if none)."""
    xyz = np.stack([np.asarray(x, np.float32), np.asarray(y, np.float32), np.asarray(z, np.float32)], 1)
    qq = np.asarray(q, np.float32).astype(np.float64)
    dq = np.clip((qq - float(q_offset)) / float(q_scale), 0, None).astype(np.float32)
    nv = int(n_vertex_rows)
    xyz_s = xshift(xyz, x_shift_cm)
    t = np.asarray(trad, np.float32).reshape(3)
    tag = xshift(t, x_shift_cm).reshape(3) if np.all(np.isfinite(t)) else np.full(3, np.nan, np.float32)
    return Event(key=key, sample="", index=0, rse=np.zeros(3, np.int64), fold=-1, points_xyz=xyz_s, points_q=dq,
                 points_cluster_id=cluster_ids(xyz_s, link_cm), cand_xyz=xyz_s[:nv].copy(), cand_q=dq[:nv].copy(),
                 tag_nu=tag, wc_reco=[], wc_status="missing", truth=None)


def cached_features(ev):
    """compute_features at the precision the sbndvtx cache stores (sbndvtx.cache.ShardWriter / CacheDataset):
    'point' features float16 (fp32_features is empty), all others float32; event-level arrays shape (dim,)."""
    feats = compute_features(ev.strip_truth())
    out = {}
    for name, dim, level in FEATURE_SPEC:
        a = np.asarray(feats[name])
        if a.dtype.kind == "f":
            if level == "point":
                a = np.clip(a, -F16_MAX, F16_MAX).astype(np.float16)
            else:
                a = a.astype(np.float32)
        out[name] = a.reshape(dim) if level == "event" else a.reshape(-1, dim)
    return out


# ---------------------------------------------------------------------------- bundle
def sha256(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


_INIT_FNS = ("uniform_", "normal_", "constant_", "ones_", "zeros_", "kaiming_uniform_", "kaiming_normal_",
             "xavier_uniform_", "xavier_normal_", "trunc_normal_", "orthogonal_")


@contextlib.contextmanager
def skip_param_init():
    """torch.nn.init functions as no-ops while a model is built.  Every parameter and persistent buffer is then
    replaced by load_state_dict(strict=True, assign=True), so the random initialisation (~0.2 s per 20 M-parameter
    model) is wasted work; outputs are unchanged (checked bit for bit on fold-0 calls)."""
    saved = {n: getattr(torch.nn.init, n) for n in _INIT_FNS if hasattr(torch.nn.init, n)}
    try:
        for n in saved:
            setattr(torch.nn.init, n, lambda tensor, *a, **k: tensor)
        yield
    finally:
        for n, f in saved.items():
            setattr(torch.nn.init, n, f)


def load_state(path, want_sha=None):
    """State dict of a bundle .pt file.  With want_sha the file is read ONCE: hashed, checked, then deserialised from
    memory; without it torch.load maps the file (mmap=True)."""
    if want_sha is None:
        return torch.load(path, map_location="cpu", weights_only=True, mmap=True)
    data = Path(path).read_bytes()
    got = hashlib.sha256(data).hexdigest()
    if got != want_sha:
        raise ValueError(f"{path}: sha256 {got} != manifest {want_sha}")
    return torch.load(io.BytesIO(data), map_location="cpu", weights_only=True)


class Vertexer:
    """Models, ranker and settings of one bundle (bundle.json + sibling files), loaded once."""

    def __init__(self, manifest, verify=True, threads=None):
        t0 = time.perf_counter()
        ms = {}

        def mark(k):
            nonlocal t0
            now = time.perf_counter()
            ms[k] = round(1e3 * (now - t0), 1)
            t0 = now

        self.path = Path(manifest).resolve()
        self.man = json.loads(self.path.read_text())
        if self.man.get("format") != FORMAT:
            raise ValueError(f"{self.path}: format {self.man.get('format')!r}, expected {FORMAT!r}")
        d = self.path.parent
        files = self.man["files"]

        def read_text(name):
            data = (d / files[name]["path"]).read_bytes()
            if verify:
                got = hashlib.sha256(data).hexdigest()
                if got != files[name]["sha256"]:
                    raise ValueError(f"{d / files[name]['path']}: sha256 {got} != manifest {files[name]['sha256']}")
            return data.decode()

        self.cfg = self.man["config"]
        self.x_shift_cm = float(self.man["x_shift_cm"])
        self.link_cm = float(self.man["link_cm"])
        self.min_points = int(self.cfg["data"]["min_points"])
        self.family = tuple(self.man["family_folds"])
        self.threads = threads
        self.models = {}
        mark("manifest")
        for which in ("pass1", "pass2"):
            pc = self.cfg[which]
            with skip_param_init():
                m = DV2Net(pc, inputs.n_voxel_features(self.cfg["inputs"], pc),
                           inputs.n_point_features(self.cfg["inputs"]), "torch")
            mark(f"build_{which}")
            sd = load_state(d / files[which]["path"], files[which]["sha256"] if verify else None)
            mark(f"read_{which}")
            # assign=True: the module takes the loaded tensors (no copy into its own, never-initialised ones)
            m.load_state_dict(sd, strict=True, assign=True)
            m.set_compute_dtype(torch.float32)
            self.models[which] = m.eval()
            del sd
            mark(f"load_state_dict_{which}")
        import lightgbm as lgb
        mark("import_lightgbm")
        self.booster = lgb.Booster(model_str=read_text("ranker"))
        mark("ranker_booster")
        iso = json.loads(read_text("isotonic"))
        self.iso_x = np.asarray(iso["X_thresholds"], np.float64)
        self.iso_y = np.asarray(iso["y_thresholds"], np.float64)
        self.iso_lo, self.iso_hi = float(iso["X_min"]), float(iso["X_max"])
        self.feature_names = list(self.man["ranker"]["feature_names"])
        if self.booster.num_feature() != len(self.feature_names):
            raise ValueError(f"ranker has {self.booster.num_feature()} features, manifest {len(self.feature_names)}")
        self.feature_fns = tuple(f.with_spec(FEATURE_SPEC) if isinstance(f, baselines.CacheFeatureAggregates)
                                 and f.spec is None else f for f in baselines.DEFAULT_FEATURE_FNS)
        rc = self.cfg["ranker"]
        if rc.get("feature_set", "all") != "all" or rc.get("model", "lambdamart") != "lambdamart":
            raise ValueError("this build deploys ranker.feature_set 'all' with the LambdaMART model only")
        mark("isotonic_and_checks")
        self.load_ms = ms

    def calibrate(self, s):
        """IsotonicRegression(out_of_bounds='clip').predict: clip to [X_min, X_max], linear interpolation."""
        return float(np.interp(np.clip(s, self.iso_lo, self.iso_hi), self.iso_x, self.iso_y))

    # ------------------------------------------------------------------ stages
    def pass1(self, ev):
        pc = self.cfg["pass1"]
        prep = infer.pass1_prepare(ev, self.cfg, pc)
        outs = infer.forward_members([self.models["pass1"]], prep["batch"], "cpu", None)
        return infer.pass1_finish(outs, prep, pc["infer"])

    def pass2(self, ev, rec1):
        pc = self.cfg["pass2"]
        prep = infer.pass2_prepare(ev, rec1, self.cfg, pc)
        if prep is None:
            return infer.empty_pass2()
        outs = infer.forward_members([self.models["pass2"]], prep["batch"], "cpu", None)
        return infer.pass2_finish(outs, prep, pc["infer"])

    def ranker_rows(self, ev, rec1, rec2):
        """(positions (n,3), sources (n,), X (n, F) float32) of the deployed family's variant, rows in the order
        of the offline union variant (base pool, then the expanded rows)."""
        fam = self.family
        rc = self.cfg["ranker"]
        edges = (self.cfg["pass1"]["target"]["shell_edges"], self.cfg["pass2"]["target"]["shell_edges"])
        df, X, names, rows_v, block_v = ranker.event_union(ev, {fam: (rec1, rec2)}, [fam], rc, edges,
                                                           self.feature_fns)
        full = list(names) + ["v_" + c for c in ranker.VARIANT_FEATURES]
        if full != self.feature_names:
            raise ValueError("ranker feature names differ from the bundle's")
        rv = rows_v[0]
        U = np.stack([df["x"], df["y"], df["z"]], 1)
        Xr = np.concatenate([np.asarray(X, np.float32)[rv], block_v[0]], 1).astype(np.float32)
        return U[rv], np.asarray(df["source"])[rv], Xr

    def run(self, ev, with_records=False):
        """Final DV2-MAX vertex of a truth-free Event.  Returns None (fewer than data.min_points points or no
        candidate) or a dict: xyz (label frame), score (P(< 1 cm)), score_raw, source, n_cand, ms (stage times)."""
        t = {}
        t0 = time.perf_counter()
        if ev.n_points < self.min_points:
            return None
        ev.features = cached_features(ev)
        t["features"] = time.perf_counter()
        with torch.no_grad():
            rec1 = self.pass1(ev)
            t["pass1"] = time.perf_counter()
            rec2 = self.pass2(ev, rec1)
            t["pass2"] = time.perf_counter()
        P, src, X = self.ranker_rows(ev, rec1, rec2)
        t["ranker_rows"] = time.perf_counter()
        if len(P) == 0:
            return None
        s = self.booster.predict(X, num_threads=int(self.threads or 1))
        rc = self.cfg["ranker"]
        from .decode import nms
        keep = nms(P, s, rc["nms_cm"], rc["top_k_out"])
        i = int(keep[0])
        t["ranker"] = time.perf_counter()
        out = {"xyz": P[i].astype(np.float64), "score": self.calibrate(s[i]), "score_raw": float(s[i]),
               "source": str(src[i]), "n_cand": int(len(P))}
        prev, ms = t0, {}
        for k, v in t.items():
            ms[k] = round(1e3 * (v - prev), 2)
            prev = v
        ms["total"] = round(1e3 * (prev - t0), 2)
        out["ms"] = ms
        if with_records:
            out.update(rec1=rec1, rec2=rec2, rows=P, sources=src, X=X, s=s)
        return out
