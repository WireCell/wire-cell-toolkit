#!/usr/bin/env python
"""Tests of the DV2-MAX vertexer: pyutil/python/DV2_Vertex.py and pyutil/python/dv2max_vertexer/.

Needs numpy, scipy, torch (CPU) and lightgbm, and for most groups the model bundle (wire-cell-data
dv2max/sbnd/r7_cloud/bundle.json, found on WIRECELL_PATH, or DV2MAX_BUNDLE=<path of bundle.json>).  Tests whose
requirements are missing are skipped with the reason.  WCT runs test/test_*.py as atomic tests; by hand:

  python pyutil/test/test_dv2max.py -v

Groups:
  TestFrameAndEvent  the call -> Event step: x shift toward the cathode and back, charge decoding, PR-vertex
                     candidate rows, traditional vertex, 3 cm connected components   (no model)
  TestNMS            greedy NMS of the ranker candidates                               (no model)
  TestBundle         bundle loading: contents, sha256 verification, format and feature-name checks
  TestCalibration    the isotonic map raw score -> P(< 1 cm)
  TestAPI            DV2_Vertex / DV2_Load: input checks, DV2_STRICT, no prediction, output frame, repeatability,
                     DV2_HOOK_LOG, DV2_THREADS
  TestPipeline       per-stage invariants on real held-out calls (pass 1, ranker rows, top-1 pick, calibration)
  TestGolden         replay of data/dv2max_golden.json (make_dv2max_golden.py): held-out calls the bundle was not
                     trained on must give the recorded vertex and score
"""
import hashlib
import importlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
PYDIR = HERE.parent / "python"
GOLDEN = HERE / "data" / "dv2max_golden.json"
BUNDLE_REL = "dv2max/sbnd/r7_cloud/bundle.json"
sys.path.insert(0, str(PYDIR))


def _missing_deps():
    miss = []
    for m in ("numpy", "scipy", "torch", "lightgbm"):
        try:
            importlib.import_module(m)
        except Exception as e:          # noqa: BLE001 - any import failure means "not available here"
            miss.append(f"{m} ({type(e).__name__})")
    return miss


MISSING = _missing_deps()
if not MISSING:
    import numpy as np


def _find_bundle():
    p = os.environ.get("DV2MAX_BUNDLE")
    if p:
        return Path(p) if Path(p).is_file() else None
    for d in os.environ.get("WIRECELL_PATH", "").split(":"):
        if d and (Path(d) / BUNDLE_REL).is_file():
            return Path(d) / BUNDLE_REL
    return None


BUNDLE = None if MISSING else _find_bundle()
needs_deps = unittest.skipIf(bool(MISSING), f"missing Python packages: {', '.join(MISSING)}")
needs_bundle = unittest.skipIf(bool(MISSING) or BUNDLE is None,
                               f"missing {', '.join(MISSING)}" if MISSING else
                               f"no {BUNDLE_REL} on WIRECELL_PATH (set DV2MAX_BUNDLE)")

_STATE = {}


def dv():
    """The DV2_Vertex module (imported once)."""
    if "dv" not in _STATE:
        os.environ.setdefault("DV2_THREADS", "1")
        _STATE["dv"] = importlib.import_module("DV2_Vertex")
    return _STATE["dv"]


def vertexer():
    """The Vertexer of the bundle, as DV2_Vertex caches it."""
    return dv()._vertexer(str(BUNDLE))


def golden():
    if "golden" not in _STATE:
        _STATE["golden"] = json.loads(GOLDEN.read_text())
    return _STATE["golden"]


def f32(a):
    return np.asarray([np.nan if v is None else v for v in a], np.float32)


def call(c, **env):
    """DV2_Vertex on a golden call (dict); returns the float32 output array.  env: temporary environment."""
    saved = {k: os.environ.get(k) for k in env}
    os.environ.update({k: str(v) for k, v in env.items()})
    try:
        out = dv().DV2_Vertex(str(BUNDLE), f32(c["x"]).tobytes(), f32(c["y"]).tobytes(), f32(c["z"]).tobytes(),
                              f32(c["q"]).tobytes(), int(c["n_vertex_rows"]), f32(c["trad"]).tobytes(),
                              float(c["q_scale"]), float(c["q_offset"]), "float32")
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
    return np.frombuffer(out, np.float32)


def event_of(c):
    from dv2max_vertexer.vertexer import event_from_call
    v = vertexer()
    return event_from_call(f32(c["x"]), f32(c["y"]), f32(c["z"]), f32(c["q"]), c["n_vertex_rows"], f32(c["trad"]),
                           c["q_scale"], c["q_offset"], v.x_shift_cm, v.link_cm)


def real_calls():
    return [c for c in golden()["calls"] if not c["key"].startswith("synthetic")]


# ---------------------------------------------------------------------------- no model
@needs_deps
class TestFrameAndEvent(unittest.TestCase):
    def setUp(self):
        from dv2max_vertexer import vertexer as vx
        self.vx = vx

    def test_xshift_moves_x_toward_the_cathode(self):
        a = np.array([[10, 1, 2], [-10, 3, 4], [0, 5, 6]], np.float32)
        b = self.vx.xshift(a, 0.31)
        s = np.float32(0.31)
        np.testing.assert_array_equal(b[:, 0], [np.float32(10) - s, np.float32(-10) + s, 0])
        np.testing.assert_array_equal(b[:, 1:], a[:, 1:])
        self.assertEqual(b.dtype, np.float32)

    def test_unshift_inverts_xshift(self):
        g = np.random.default_rng(1)
        for _ in range(200):
            p = g.uniform([-200, -200, 0], [200, 200, 500]).astype(np.float32)
            if abs(p[0]) < 1:
                continue
            back = self.vx.unshift_x(self.vx.xshift(p, 0.31)[0], 0.31)
            np.testing.assert_allclose(back, p, atol=2e-5)

    def _event(self, xyz, q, nv=0, trad=(np.nan,) * 3, qs=0.1, qo=-1000.0):
        xyz = np.asarray(xyz, np.float32)
        return self.vx.event_from_call(xyz[:, 0], xyz[:, 1], xyz[:, 2], np.asarray(q, np.float32), nv,
                                       np.asarray(trad, np.float32), qs, qo, 0.31, 3.0)

    def test_charge_is_decoded_and_clipped_at_zero(self):
        dq = np.array([0, 500, 12000, 3.5e4], np.float32)
        q = dq * np.float32(0.1) - np.float32(1000)
        q = np.append(q, np.float32(-2000))                       # below the offset -> dQ 0
        ev = self._event(np.zeros((5, 3)) + [[50, 0, 100]], q)
        np.testing.assert_allclose(ev.points_q, np.append(dq, 0), rtol=1e-5, atol=1e-2)
        self.assertTrue(np.all(ev.points_q >= 0))

    def test_vertex_rows_are_candidates_and_trad_is_shifted(self):
        xyz = np.array([[10, 0, 100], [-20, 5, 200], [30, 1, 300], [40, 2, 310]], np.float32)
        ev = self._event(xyz, np.full(4, 1000), nv=2, trad=(12, 3, 150))
        np.testing.assert_array_equal(ev.cand_xyz, self.vx.xshift(xyz[:2], 0.31))
        np.testing.assert_array_equal(ev.points_xyz, self.vx.xshift(xyz, 0.31))
        np.testing.assert_array_equal(ev.tag_nu, self.vx.xshift([12, 3, 150], 0.31)[0])
        self.assertTrue(ev.has_tag_nu)
        self.assertEqual(ev.n_points, 4)

    def test_missing_traditional_vertex(self):
        ev = self._event(np.array([[10, 0, 100]] * 3, np.float32), np.full(3, 1000), nv=1)
        self.assertTrue(np.all(np.isnan(ev.tag_nu)))
        self.assertFalse(ev.has_tag_nu)

    def test_clusters_are_3cm_connected_components_largest_first(self):
        a = [[0, 0, 100 + 2.9 * i] for i in range(6)]            # linked chain (2.9 cm < 3 cm)
        b = [[50, 0, 100 + 2.9 * i] for i in range(3)]           # a smaller chain far away
        c = [[0, 30, 100], [0, 33.1, 100]]                       # 3.1 cm apart: not linked
        ev = self._event(a + b + c, np.full(11, 1000))
        ids = ev.points_cluster_id
        self.assertEqual(set(ids[:6]), {0})                     # largest component = 0
        self.assertEqual(set(ids[6:9]), {1})
        self.assertEqual(len({ids[9], ids[10]}), 2)
        self.assertEqual(len(set(ids)), 4)


@needs_deps
class TestNMS(unittest.TestCase):
    def setUp(self):
        from dv2max_vertexer.decode import nms
        self.nms = nms

    def test_highest_score_first_and_suppression_radius(self):
        P = [[0, 0, 0], [0.5, 0, 0], [5, 0, 0], [5.2, 0, 0], [20, 0, 0]]
        s = [1.0, 3.0, 2.0, 0.5, -1.0]
        self.assertEqual(list(self.nms(P, s, 1.0)), [1, 2, 4])  # 0 and 3 are within 1 cm of a higher one

    def test_top_k(self):
        P = [[0, 0, 10 * i] for i in range(6)]
        self.assertEqual(list(self.nms(P, [6, 5, 4, 3, 2, 1], 1.0, 3)), [0, 1, 2])

    def test_ties_keep_input_order(self):
        self.assertEqual(list(self.nms([[0, 0, 0], [10, 0, 0]], [1.0, 1.0], 1.0)), [0, 1])


# ---------------------------------------------------------------------------- bundle
@needs_bundle
class TestBundle(unittest.TestCase):
    def test_contents(self):
        v = vertexer()
        man = json.loads(BUNDLE.read_text())
        self.assertEqual(man["format"], "dv2max-bundle-v1")
        self.assertEqual(len(v.feature_names), v.booster.num_feature())
        self.assertEqual(len(v.feature_names), man["ranker"]["n_features"])
        self.assertAlmostEqual(v.x_shift_cm, 0.31)
        self.assertAlmostEqual(v.link_cm, 3.0)
        self.assertEqual(v.min_points, 5)
        self.assertEqual(v.family, (0,))
        self.assertEqual(set(v.models), {"pass1", "pass2"})

    def _copy(self, edit_manifest=None, edit_isotonic=False):
        """A bundle in a temp dir: big files symlinked, manifest and isotonic map copied (and edited)."""
        d = Path(tempfile.mkdtemp(prefix="dv2max_bundle_"))
        self.addCleanup(shutil.rmtree, d, True)
        man = json.loads(BUNDLE.read_text())
        for name, f in man["files"].items():
            src = BUNDLE.parent / f["path"]
            if name == "isotonic":
                t = src.read_text()
                if edit_isotonic:
                    iso = json.loads(t)
                    iso["y_thresholds"][-1] = 0.5
                    t = json.dumps(iso)
                (d / f["path"]).write_text(t)
            else:
                (d / f["path"]).symlink_to(src)
        if edit_manifest:
            edit_manifest(man)
        (d / "bundle.json").write_text(json.dumps(man))
        return d / "bundle.json"

    def _load(self, path, verify=True):
        from dv2max_vertexer.vertexer import Vertexer
        return Vertexer(path, verify=verify, threads=1)

    def test_an_intact_copy_loads(self):
        v = self._load(self._copy())
        self.assertEqual(v.feature_names, vertexer().feature_names)

    def test_edited_file_fails_the_sha256_check(self):
        with self.assertRaisesRegex(ValueError, "sha256"):
            self._load(self._copy(edit_isotonic=True))

    def test_wrong_weights_hash_fails(self):
        def bad(m):
            m["files"]["pass1"]["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "sha256"):
            self._load(self._copy(edit_manifest=bad))

    def test_verify_false_skips_the_check(self):
        v = self._load(self._copy(edit_isotonic=True), verify=False)
        self.assertEqual(v.iso_y[-1], 0.5)

    def test_wrong_format_is_rejected(self):
        def bad(m):
            m["format"] = "dv2max-bundle-v0"
        with self.assertRaisesRegex(ValueError, "format"):
            self._load(self._copy(edit_manifest=bad))

    def test_feature_name_mismatch_is_rejected(self):
        def bad(m):
            m["ranker"]["feature_names"] = m["ranker"]["feature_names"][:-1]
        with self.assertRaisesRegex(ValueError, "features"):
            self._load(self._copy(edit_manifest=bad))


@needs_bundle
class TestCalibration(unittest.TestCase):
    def test_monotone_clipped_probability(self):
        v = vertexer()
        grid = np.linspace(v.iso_lo - 5, v.iso_hi + 5, 2001)
        p = np.array([v.calibrate(s) for s in grid])
        self.assertTrue(np.all(np.diff(p) >= 0))
        self.assertTrue(np.all((p >= 0) & (p <= 1)))
        self.assertEqual(v.calibrate(v.iso_lo - 100), v.calibrate(v.iso_lo))
        self.assertEqual(v.calibrate(v.iso_hi + 100), v.calibrate(v.iso_hi))
        np.testing.assert_allclose(p, np.interp(np.clip(grid, v.iso_lo, v.iso_hi), v.iso_x, v.iso_y))


# ---------------------------------------------------------------------------- DV2_Vertex API
@needs_bundle
class TestAPI(unittest.TestCase):
    def setUp(self):
        self.c = real_calls()[0]
        self._strict = os.environ.pop("DV2_STRICT", None)

    def tearDown(self):
        if self._strict is not None:
            os.environ["DV2_STRICT"] = self._strict

    def _bad_sizes(self):
        x = np.zeros(5, np.float32).tobytes()
        return dv().DV2_Vertex(str(BUNDLE), x, np.zeros(4, np.float32).tobytes(), x, x, 0,
                               np.zeros(3, np.float32).tobytes(), 0.1, -1000.0, "float32")

    def test_bad_input_raises_by_default(self):
        with self.assertRaisesRegex(ValueError, "bad input sizes"):
            self._bad_sizes()

    def test_dv2_strict_0_returns_no_prediction_instead(self):
        os.environ["DV2_STRICT"] = "0"
        self.assertEqual(self._bad_sizes(), b"")
        self.assertTrue(dv()._LAST["reason"].startswith("error: ValueError"))

    def test_fewer_than_min_points_is_no_prediction(self):
        for c in golden()["calls"]:
            if c["key"].startswith("synthetic"):
                with self.subTest(c["key"]):
                    self.assertEqual(len(call(c)), 0)
                    self.assertTrue(dv()._LAST["reason"].startswith("no prediction"))

    def test_output_is_the_raw_frame_vertex_and_calibrated_score(self):
        from dv2max_vertexer.vertexer import unshift_x
        out = call(self.c)
        self.assertEqual(out.dtype, np.float32)
        self.assertEqual(len(out), 4)
        v = vertexer()
        r = v.run(event_of(self.c))
        np.testing.assert_allclose(out[:3], unshift_x(r["xyz"], v.x_shift_cm), atol=1e-4)
        self.assertAlmostEqual(float(out[3]), r["score"], places=6)
        self.assertTrue(0 <= out[3] <= 1)

    def test_same_call_twice_gives_the_same_bytes(self):
        self.assertEqual(call(self.c).tobytes(), call(self.c).tobytes())

    def test_hook_log_writes_one_record_per_call(self):
        with tempfile.TemporaryDirectory() as d:
            log = Path(d) / "hook.jsonl"
            call(self.c, DV2_HOOK_LOG=log)
            call(golden()["calls"][-1], DV2_HOOK_LOG=log)          # a synthetic no-prediction call
            recs = [json.loads(line) for line in log.read_text().splitlines()]
        self.assertEqual(len(recs), 2)
        for k in ("n", "n_vertex_rows", "pred", "score", "score_raw", "source", "n_cand", "stage_ms", "ms"):
            self.assertIn(k, recs[0])
        self.assertIsNone(recs[1]["pred"])

    def test_dv2_load_reports_the_load_breakdown(self):
        rep = json.loads(dv().DV2_Load(str(BUNDLE)))
        for k in ("import_numpy", "import_torch", "import_dv2max_vertexer", "load", "load_steps"):
            self.assertIn(k, rep)

    def _threads_in_subprocess(self, dv2_threads, with_dv2=True):
        """(intra, interop) torch threads of a fresh process after DV2_Vertex set them up (or plain torch)."""
        code = ("import os, sys; sys.path.insert(0, %r); " % str(PYDIR)
                + ("import DV2_Vertex as d; d._threads(); " if with_dv2 else "")
                + "import torch; print(torch.get_num_threads(), torch.get_num_interop_threads())")
        env = {k: v for k, v in os.environ.items() if k not in ("DV2_THREADS", "OMP_NUM_THREADS")}
        if dv2_threads is not None:
            env["DV2_THREADS"] = dv2_threads
        out = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True, check=True)
        return tuple(int(x) for x in out.stdout.split())

    def test_dv2_threads_sets_the_torch_thread_count(self):
        self.assertEqual(self._threads_in_subprocess("3"), (3, 1))
        self.assertEqual(self._threads_in_subprocess(None), (1, 1))      # default

    def test_dv2_threads_0_leaves_torch_alone(self):
        self.assertEqual(self._threads_in_subprocess("0"), self._threads_in_subprocess(None, with_dv2=False))

    def test_load_report_records_the_thread_change(self):
        rep = json.loads(dv().DV2_Load(str(BUNDLE)))
        tt = rep["torch_threads"]
        self.assertEqual(set(tt), {"DV2_THREADS", "before", "after"})
        n = int(os.environ.get("DV2_THREADS", "1") or 1)
        self.assertEqual(tt["after"]["intra"], n if n > 0 else tt["before"]["intra"])


# ---------------------------------------------------------------------------- pipeline + golden
@needs_bundle
class TestPipeline(unittest.TestCase):
    def test_stage_invariants_on_held_out_calls(self):
        v = vertexer()
        for c in real_calls():
            with self.subTest(c["key"]):
                ev = event_of(c)
                r = v.run(ev, with_records=True)
                self.assertIsNotNone(r)
                modes = np.asarray(r["rec1"]["modes"])
                self.assertGreater(len(modes), 0)
                self.assertTrue(np.all(np.isfinite(modes[:, :3])))
                X, P, s = r["X"], np.asarray(r["rows"]), np.asarray(r["s"])
                self.assertEqual(X.dtype, np.float32)
                self.assertEqual(X.shape, (r["n_cand"], len(v.feature_names)))
                self.assertEqual(len(P), r["n_cand"])
                self.assertEqual(len(r["sources"]), r["n_cand"])
                self.assertTrue(np.all(np.isfinite(P)))
                # the pick is the best-scored candidate, and its score is calibrated
                np.testing.assert_array_equal(r["xyz"], P[int(np.argmax(s))].astype(np.float64))
                self.assertEqual(r["score_raw"], float(s.max()))
                self.assertEqual(r["score"], v.calibrate(r["score_raw"]))
                # every PR-graph vertex row is a candidate; the pick lies near the cloud
                self.assertGreaterEqual(r["n_cand"], min(c["n_vertex_rows"], 1))
                lo, hi = ev.points_xyz.min(0) - 10, ev.points_xyz.max(0) + 10
                self.assertTrue(np.all((r["xyz"] >= lo) & (r["xyz"] <= hi)))


@needs_bundle
class TestGolden(unittest.TestCase):
    def test_fixture_matches_this_bundle(self):
        sha = hashlib.sha256(BUNDLE.read_bytes()).hexdigest()
        self.assertEqual(golden()["meta"]["bundle_sha256"], sha,
                         "the golden answers were made with another bundle: rerun make_dv2max_golden.py")

    def test_fixture_covers_the_strata(self):
        calls = real_calls()
        keys = [c["key"] for c in calls]
        self.assertGreaterEqual(len(calls), 20)
        self.assertTrue(any(k.startswith("bnb:") for k in keys) and any(k.startswith("nue:") for k in keys))
        self.assertTrue(any(k.endswith(":p0") for k in keys) and any(k.endswith(":p1") for k in keys))
        self.assertTrue(any(c["trad"][0] is None for c in calls))
        self.assertTrue(any(c["key"].startswith("synthetic") for c in golden()["calls"]))

    def test_held_out_calls_reproduce(self):
        tol = golden()["meta"]["tolerance"]
        worst_d, worst_s = 0.0, 0.0
        for c in golden()["calls"]:
            with self.subTest(c["key"]):
                out = call(c)
                if not c["expect"]:
                    self.assertEqual(len(out), 0)
                    continue
                self.assertEqual(len(out), 4)
                e = np.asarray(c["expect"], np.float64)
                d = float(np.linalg.norm(out[:3].astype(np.float64) - e[:3]))
                ds = abs(float(out[3]) - e[3])
                worst_d, worst_s = max(worst_d, d), max(worst_s, ds)
                self.assertLess(d, tol["position_cm"], f"vertex moved {d:.4g} cm")
                self.assertLess(ds, tol["score"], f"score changed by {ds:.3g}")
        print(f"\n[golden] {len(golden()['calls'])} calls: max vertex deviation {worst_d:.3g} cm, "
              f"max score deviation {worst_s:.3g}", file=sys.stderr)


if __name__ == "__main__":
    unittest.main()
