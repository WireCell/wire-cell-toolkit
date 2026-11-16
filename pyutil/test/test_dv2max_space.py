#!/usr/bin/env python
"""Tests of the DV2-MAX space-point input (bundle format dv2max-bundle-v2, e.g. aurora_space): the call's space
points (Bee clustering-layer points of the imaging clusters touching its cloud) added to the model input.

Needs numpy and scipy; the bundle groups also torch and lightgbm and the bundles (wire-cell-data
dv2max/sbnd/aurora_space/bundle.json and dv2max/sbnd/r7_cloud/bundle.json on WIRECELL_PATH, or DV2MAX_SPACE_BUNDLE /
DV2MAX_BUNDLE).  Missing requirements skip the group with the reason.  By hand:

  python pyutil/test/test_dv2max_space.py -v

Groups:
  TestAddSpace      event_from_call with space points == the training cache's rule (a verbatim copy of
                    sbnd-vertexer aurora/prep_aurora_cache.add_space and the fg1 helpers it calls, below), point by
                    point; no space point; unmatched cloud points; cluster-size ties         (no model)
  TestManifest      format / "space" section checks at load                                (torch, no model file)
  TestDigest        input_digest: cloud-only digests unchanged, space points change it     (no model)
  TestSpaceAPI      DV2_Vertex / DV2_Load with the space bundle: the load report's space section, a call
                    without space points / with partial ones / a cloud-only bundle given space points raise,
                    an empty space set runs, the hook record counts the space points
  TestSpaceGolden   replay of data/dv2max_space_golden.json (make_dv2max_space_golden.py): hold-out calls the
                    bundle was not trained on must give the recorded vertex and score
"""
import base64
import hashlib
import importlib
import json
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
PYDIR = HERE.parent / "python"
GOLDEN = HERE / "data" / "dv2max_space_golden.json"
SPACE_REL = "dv2max/sbnd/aurora_space/bundle.json"
CLOUD_REL = "dv2max/sbnd/r7_cloud/bundle.json"
sys.path.insert(0, str(PYDIR))


def _missing(mods):
    miss = []
    for m in mods:
        try:
            importlib.import_module(m)
        except Exception as e:          # noqa: BLE001 - any import failure means "not available here"
            miss.append(f"{m} ({type(e).__name__})")
    return miss


MISS_NP = _missing(("numpy", "scipy"))
MISS_ALL = MISS_NP + _missing(("torch", "lightgbm"))
if not MISS_NP:
    import numpy as np


def _find(env, rel):
    p = os.environ.get(env)
    if p:
        return Path(p) if Path(p).is_file() else None
    for d in os.environ.get("WIRECELL_PATH", "").split(":"):
        if d and (Path(d) / rel).is_file():
            return Path(d) / rel
    return None


SPACE_BUNDLE = None if MISS_ALL else _find("DV2MAX_SPACE_BUNDLE", SPACE_REL)
CLOUD_BUNDLE = None if MISS_ALL else _find("DV2MAX_BUNDLE", CLOUD_REL)
needs_np = unittest.skipIf(bool(MISS_NP), f"missing Python packages: {', '.join(MISS_NP)}")
needs_all = unittest.skipIf(bool(MISS_ALL), f"missing Python packages: {', '.join(MISS_ALL)}")
needs_space = unittest.skipIf(bool(MISS_ALL) or SPACE_BUNDLE is None,
                              f"missing {', '.join(MISS_ALL)}" if MISS_ALL else
                              f"no {SPACE_REL} on WIRECELL_PATH (set DV2MAX_SPACE_BUNDLE)")
needs_both = unittest.skipIf(bool(MISS_ALL) or SPACE_BUNDLE is None or CLOUD_BUNDLE is None,
                             f"missing {', '.join(MISS_ALL)}" if MISS_ALL else
                             f"needs both {SPACE_REL} and {CLOUD_REL} on WIRECELL_PATH")
_STATE = {}


def dv():
    if "dv" not in _STATE:
        os.environ.setdefault("DV2_THREADS", "1")
        _STATE["dv"] = importlib.import_module("DV2_Vertex")
    return _STATE["dv"]


# ---------------------------------------------------------------------------- the training rule, verbatim
# sbnd-vertexer fg1/prep_fg1_cache.py xshift / cluster_ids (LINK_CM 3.0) and aurora/prep_aurora_cache.py add_space,
# with the h5 reads (hs["space"][k][a:b], hs["space/q"], hs["space/cluster_id"]) replaced by the arrays S_raw, q_raw,
# cid_raw.  Kept as a copy on purpose: the port in dv2max_vertexer.vertexer must reproduce it, not share it.
REF_LINK_CM, REF_MATCH_CM = 3.0, 1.0


def ref_xshift(a, shift):
    a = np.array(a, np.float32, copy=True).reshape(-1, 3)
    a[:, 0] -= np.sign(a[:, 0]) * np.float32(shift)
    return a


def ref_cluster_ids(xyz):
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree
    n = len(xyz)
    if n == 0:
        return np.zeros(0, np.int32)
    pairs = cKDTree(xyz).query_pairs(REF_LINK_CM, output_type="ndarray")
    g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    _, lab = connected_components(g, directed=False)
    order = np.argsort(-np.bincount(lab), kind="stable")
    remap = np.empty_like(order)
    remap[order] = np.arange(len(order))
    return remap[lab].astype(np.int32)


def ref_add_space(ev, S_raw, q_raw, cid_raw, q_ratio, shift):
    from scipy.spatial import cKDTree
    if len(q_raw) == 0:                          # b <= a
        return ev
    S = ref_xshift(S_raw, shift)
    qs = np.asarray(q_raw).astype(np.float32)
    cid = np.asarray(cid_raw)
    u, inv, cnt = np.unique(cid, return_inverse=True, return_counts=True)
    order = np.argsort(-cnt, kind="stable")
    rank = np.empty_like(order)
    rank[order] = np.arange(len(order))
    cs = rank[inv].astype(np.int32)
    C = ev.points_xyz
    d, nn = cKDTree(S).query(C)
    cc = np.where(d <= REF_MATCH_CM, cs[nn], -1).astype(np.int32)
    un = cc < 0
    if un.any():
        cc[un] = ref_cluster_ids(C[un]) + len(u)
    ev.points_xyz = np.concatenate([S, C]).astype(np.float32)
    ev.points_q = np.concatenate([qs, (ev.points_q * q_ratio).astype(np.float32)])
    ev.points_cluster_id = np.concatenate([cs, cc])
    return ev


CFG = {"touch_cm": 1.0, "match_cm": 1.0, "q_ratio": 0.2474433061809911, "link_cm": 3.0}
SHIFT = 0.29


def synthetic_call(seed, n_cloud=400, n_clusters=5, near=True):
    """A cloud of tracks and, near it (near=True) or far from it, n_clusters space-point clusters (float32)."""
    rng = np.random.default_rng(seed)
    t = rng.uniform(0, 60, n_cloud)
    d = rng.normal(size=3)
    d /= np.linalg.norm(d)
    cloud = (np.array([-50.0, 20.0, 200.0]) + t[:, None] * d + rng.normal(0, 0.3, (n_cloud, 3))).astype(np.float32)
    cloud[: n_cloud // 3, 0] *= -1                       # points on both sides of the cathode (x sign matters)
    q = (rng.uniform(0, 5e4, n_cloud) * 0.1 - 1000).astype(np.float32)
    S, Q, CID = [], [], []
    ids = rng.choice(np.arange(3, 400), n_clusters, replace=False)
    for k, cid in enumerate(sorted(ids)):
        m = int(rng.integers(20, 300))
        base = cloud[rng.integers(0, n_cloud)] + (0 if near else 500.0)
        S.append(base + rng.normal(0, 1.5, (m, 3)))
        Q.append(rng.uniform(0, 4e3, m))
        CID.append(np.full(m, cid))
    S = np.concatenate(S).astype(np.float32)
    return (cloud, q, int(rng.integers(0, 20)), np.array([cloud[0, 0], cloud[0, 1], cloud[0, 2]], np.float32),
            S, np.concatenate(Q).astype(np.float32), np.concatenate(CID).astype(np.int32))


# ---------------------------------------------------------------------------- no model
@needs_np
class TestAddSpace(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from dv2max_vertexer import vertexer as V
        cls.V = V

    def both(self, cloud, q, nv, trad, S, Q, CID, qs=0.1, qo=-1000.0):
        V = self.V
        sp = (S[:, 0], S[:, 1], S[:, 2], Q, CID)
        got = V.event_from_call(cloud[:, 0], cloud[:, 1], cloud[:, 2], q, nv, trad, qs, qo, SHIFT, 3.0,
                                space=sp, space_cfg=CFG)
        base = V.event_from_call(cloud[:, 0], cloud[:, 1], cloud[:, 2], q, nv, trad, qs, qo, SHIFT, 3.0)
        want = ref_add_space(base, S, Q, CID, CFG["q_ratio"], SHIFT)
        return got, want

    def same(self, got, want):
        for f in ("points_xyz", "points_q", "points_cluster_id", "cand_xyz", "cand_q", "tag_nu"):
            a, b = getattr(got, f), getattr(want, f)
            self.assertEqual(a.dtype, b.dtype, f)
            np.testing.assert_array_equal(a, b, err_msg=f)

    def test_matches_the_training_rule_point_by_point(self):
        for seed in range(6):
            got, want = self.both(*synthetic_call(seed))
            self.same(got, want)

    def test_space_points_come_first_then_the_cloud(self):
        cloud, q, nv, trad, S, Q, CID = synthetic_call(1)
        got, _ = self.both(cloud, q, nv, trad, S, Q, CID)
        self.assertEqual(got.n_points, len(S) + len(cloud))
        np.testing.assert_array_equal(got.points_xyz[: len(S)], self.V.xshift(S, SHIFT))
        np.testing.assert_array_equal(got.points_q[: len(S)], Q)
        self.assertEqual(got.n_cand, nv)                  # candidates stay the cloud's vertex rows

    def test_no_space_point_keeps_the_cloud_only_event_unscaled(self):
        cloud, q, nv, trad, S, Q, CID = synthetic_call(2)
        e = np.zeros(0, np.float32)
        got, want = self.both(cloud, q, nv, trad, e.reshape(0, 3), e, np.zeros(0, np.int32))
        self.same(got, want)
        base = self.V.event_from_call(cloud[:, 0], cloud[:, 1], cloud[:, 2], q, nv, trad, 0.1, -1000.0, SHIFT, 3.0)
        np.testing.assert_array_equal(got.points_q, base.points_q)   # dQ NOT multiplied by q_ratio

    def test_unmatched_cloud_points_get_components_after_the_space_clusters(self):
        cloud, q, nv, trad, S, Q, CID = synthetic_call(3, near=False)
        got, want = self.both(cloud, q, nv, trad, S, Q, CID)
        self.same(got, want)
        n_sc = len(np.unique(CID))
        self.assertTrue((got.points_cluster_id[len(S):] >= n_sc).all())

    def test_cluster_ids_are_ranked_by_size_ties_by_increasing_id(self):
        cloud, q, nv, trad, S, Q, CID = synthetic_call(4)
        m = len(S) // 2                                   # two equal-size clusters: lower id ranks first
        CID = np.full(len(S), 9999, np.int32)             # (an odd point out: a 1-point cluster, ranked last)
        CID[:m], CID[m:2 * m] = 50, 7
        got, want = self.both(cloud, q, nv, trad, S, Q, CID)
        self.same(got, want)
        self.assertEqual(int(got.points_cluster_id[m]), 0)       # id 7
        self.assertEqual(int(got.points_cluster_id[0]), 1)       # id 50

    def test_missing_space_points_raise_for_a_space_bundle(self):
        cloud, q, nv, trad, *_ = synthetic_call(5)
        with self.assertRaises(ValueError):
            self.V.event_from_call(cloud[:, 0], cloud[:, 1], cloud[:, 2], q, nv, trad, 0.1, -1000.0, SHIFT, 3.0,
                                   space=None, space_cfg=CFG)


@needs_np
class TestDigest(unittest.TestCase):
    def test_cloud_only_digest_is_unchanged(self):
        D = importlib.import_module("DV2_Vertex")         # imports numpy only (torch on the first call)
        x = np.arange(6, dtype=np.float32)
        t = np.full(3, np.nan, np.float32)
        h = hashlib.sha1()
        for b in (x.tobytes(),) * 4 + (t.tobytes(),):
            h.update(b)
        h.update(f"{2}:{0.1!r}:{-1000.0!r}:float32".encode())
        self.assertEqual(D.input_digest(*(x.tobytes(),) * 4, 2, t.tobytes(), 0.1, -1000.0), h.hexdigest())
        sp = (x.tobytes(),) * 4 + (np.arange(6, dtype=np.int32).tobytes(),)
        with_space = D.input_digest(*(x.tobytes(),) * 4, 2, t.tobytes(), 0.1, -1000.0, "float32", sp)
        self.assertNotEqual(with_space, h.hexdigest())
        h.update(b":space:")
        for b in sp:
            h.update(b)
        self.assertEqual(with_space, h.hexdigest())


@needs_all
class TestManifest(unittest.TestCase):
    """The format checks run before any model file is read, so a manifest alone is enough."""

    def load(self, man):
        from dv2max_vertexer.vertexer import Vertexer
        d = Path(tempfile.mkdtemp(prefix="dv2space_"))
        try:
            (d / "bundle.json").write_text(json.dumps(man))
            return Vertexer(d / "bundle.json")
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def test_v2_needs_a_space_section(self):
        with self.assertRaisesRegex(ValueError, "space"):
            self.load({"format": "dv2max-bundle-v2"})

    def test_v2_space_section_needs_every_key(self):
        sp = dict(CFG)
        sp.pop("q_ratio")
        with self.assertRaisesRegex(ValueError, "space"):
            self.load({"format": "dv2max-bundle-v2", "space": sp})

    def test_touch_must_be_positive(self):
        with self.assertRaisesRegex(ValueError, "touch_cm"):
            self.load({"format": "dv2max-bundle-v2", "space": dict(CFG, touch_cm=0)})

    def test_v1_with_a_space_section_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "v2"):
            self.load({"format": "dv2max-bundle-v1", "space": CFG})

    def test_unknown_format_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "format"):
            self.load({"format": "dv2max-bundle-v3"})


# ---------------------------------------------------------------------------- with the bundles
def _args(cloud, q, nv, trad, qs=0.1, qo=-1000.0):
    return (cloud[:, 0].tobytes(), cloud[:, 1].tobytes(), cloud[:, 2].tobytes(), q.tobytes(), int(nv), trad.tobytes(),
            float(qs), float(qo), "float32")


def _sp(S, Q, CID):
    return (np.ascontiguousarray(S[:, 0]).tobytes(), np.ascontiguousarray(S[:, 1]).tobytes(),
            np.ascontiguousarray(S[:, 2]).tobytes(), Q.tobytes(), CID.astype(np.int32).tobytes())


@needs_space
class TestSpaceAPI(unittest.TestCase):
    def setUp(self):
        self.saved = {k: os.environ.get(k) for k in ("DV2_STRICT", "DV2_HOOK_LOG")}
        os.environ.pop("DV2_STRICT", None)
        os.environ.pop("DV2_HOOK_LOG", None)

    def tearDown(self):
        for k, v in self.saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v

    def test_load_report_has_the_space_section(self):
        rep = json.loads(dv().DV2_Load(str(SPACE_BUNDLE)))
        man = json.loads(Path(SPACE_BUNDLE).read_text())
        self.assertEqual(man["format"], "dv2max-bundle-v2")
        self.assertEqual(rep["space"], {k: float(man["space"][k]) for k in ("touch_cm", "match_cm", "q_ratio", "link_cm")})
        self.assertGreater(rep["space"]["touch_cm"], 0)

    def test_call_without_space_points_raises(self):
        cloud, q, nv, trad, *_ = synthetic_call(10)
        with self.assertRaisesRegex(ValueError, "space points"):
            dv().DV2_Vertex(str(SPACE_BUNDLE), *_args(cloud, q, nv, trad))

    def test_partial_space_arguments_raise(self):
        cloud, q, nv, trad, S, Q, CID = synthetic_call(11)
        sp = _sp(S, Q, CID)
        with self.assertRaisesRegex(ValueError, "all five"):
            dv().DV2_Vertex(str(SPACE_BUNDLE), *_args(cloud, q, nv, trad), *sp[:3])

    def test_space_size_mismatch_raises(self):
        cloud, q, nv, trad, S, Q, CID = synthetic_call(12)
        sp = list(_sp(S, Q, CID))
        sp[3] = Q[:-1].tobytes()
        with self.assertRaisesRegex(ValueError, "space sizes"):
            dv().DV2_Vertex(str(SPACE_BUNDLE), *_args(cloud, q, nv, trad), *sp)

    def test_empty_space_set_runs_and_hook_counts_the_points(self):
        cloud, q, nv, trad, S, Q, CID = synthetic_call(13)
        e = np.zeros(0, np.float32)
        d = Path(tempfile.mkdtemp(prefix="dv2space_"))
        try:
            os.environ["DV2_HOOK_LOG"] = str(d / "hook.jsonl")
            out0 = dv().DV2_Vertex(str(SPACE_BUNDLE), *_args(cloud, q, nv, trad), e.tobytes(), e.tobytes(),
                                   e.tobytes(), e.tobytes(), np.zeros(0, np.int32).tobytes())
            out1 = dv().DV2_Vertex(str(SPACE_BUNDLE), *_args(cloud, q, nv, trad), *_sp(S, Q, CID))
            recs = [json.loads(line) for line in open(d / "hook.jsonl")]
        finally:
            shutil.rmtree(d, ignore_errors=True)
        self.assertEqual([r["n_space"] for r in recs], [0, len(Q)])
        for out in (out0, out1):
            self.assertIn(len(out), (0, 16))
        self.assertNotEqual(recs[0]["input_sha1"], recs[1]["input_sha1"])

    @needs_both
    def test_cloud_only_bundle_given_space_points_raises(self):
        cloud, q, nv, trad, S, Q, CID = synthetic_call(14)
        with self.assertRaisesRegex(ValueError, "takes no space points"):
            dv().DV2_Vertex(str(CLOUD_BUNDLE), *_args(cloud, q, nv, trad), *_sp(S, Q, CID))


def b64f(s, dt):
    return np.frombuffer(base64.b64decode(s), dt)


@needs_space
class TestSpaceGolden(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not GOLDEN.is_file():
            raise unittest.SkipTest(f"no {GOLDEN.name}")
        cls.G = json.loads(GOLDEN.read_text())

    def test_fixture_matches_this_bundle(self):
        sha = hashlib.sha256(Path(SPACE_BUNDLE).read_bytes()).hexdigest()
        self.assertEqual(self.G["bundle_sha256"], sha, "fixture made with another bundle: rerun make_dv2max_space_golden.py")

    def test_held_out_calls_reproduce(self):
        tol = float(self.G["tolerance_cm"])
        for c in self.G["calls"]:
            with self.subTest(key=c["key"]):
                a = {k: b64f(c[k], np.float32) for k in ("x", "y", "z", "q", "sx", "sy", "sz", "sq", "trad")}
                cid = b64f(c["scid"], np.int32)
                out = np.frombuffer(dv().DV2_Vertex(
                    str(SPACE_BUNDLE), a["x"].tobytes(), a["y"].tobytes(), a["z"].tobytes(), a["q"].tobytes(),
                    int(c["n_vertex_rows"]), a["trad"].tobytes(), float(c["q_scale"]), float(c["q_offset"]),
                    "float32", a["sx"].tobytes(), a["sy"].tobytes(), a["sz"].tobytes(), a["sq"].tobytes(),
                    cid.tobytes()), np.float32)
                want = np.asarray(c["out"], np.float32)
                self.assertEqual(len(out), len(want))
                if len(want):
                    self.assertLess(float(np.linalg.norm(out[:3].astype(float) - want[:3])), tol)
                    self.assertLess(abs(float(out[3]) - float(want[3])), float(self.G["score_tolerance"]))


if __name__ == "__main__":
    unittest.main()
