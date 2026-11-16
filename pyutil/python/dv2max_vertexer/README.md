# dv2max_vertexer: DV2-MAX inference for the in-chain DL-vertex step

Inference-only copy of the DV2-MAX vertexer, called by WCT through `pyutil/python/DV2_Vertex.py`
(`WCPPyUtil::DV2_Vertex`). One call = one neutrino candidate's DL-vertex cloud. The result is the DV2-MAX
top-1 vertex in the raw reco frame plus a calibrated P(< 1 cm).

- Runtime dependencies: numpy, scipy, torch (CPU), lightgbm. No spconv, torch_scatter, pandas, h5py or sklearn.
  Tested versions:

  | environment | Python | numpy | scipy | torch | lightgbm | checked with |
  |---|---|---|---|---|---|---|
  | pyenv-dl (Ubuntu 22.04 container, glibc 2.35) | 3.14.5 | 2.5.3 | 1.18.1 | 2.9.1+cpu | 4.7.0 | validation on 5,038 held-out calls, the golden fixture, the tests |
  | SBND UPS stack (SL7, glibc 2.17) | 3.9 | 2.0.2 | 1.13.1 | 2.6.0 | 4.6.0 | review of PR #543: identical results in-chain and standalone, no warnings |

  lightgbm 4.7 has no wheel for glibc 2.17; 4.6.0 reads this bundle's ranker (`ranker_fold0.txt`) and gives the
  same results.
- Nothing is imported from sbnd-vertexer.
- Model bundles: `dv2max/sbnd/aurora_space/` (SBND default, space points, format v2) and `dv2max/sbnd/r7_cloud/`
  (cloud only) in [wire-cell-data](https://github.com/WireCell/wire-cell-data) (`bundle.json` and the files it
  lists, each with its sha256), not in this repository. A `dv2_*` method needs the selected wire-cell-data
  directory on `WIRECELL_PATH`; without it the job fails at configure ("not found on WIRECELL_PATH").

## Choosing the vertexing method

One option, `vertex_method`, on every SBND PR job: `wct-pr-perevt.jsonnet`, `wct-pr.jsonnet`, the LArSoft
library `wcls-img-clus-matching-pr-lib.jsonnet`, or `pr()` in `clus.jsonnet`. The SBND jobs default to
`dv2_split_move`; pass another value to change it:

```
wire-cell ... -c pgrapher/experiment/sbnd/wct-pr-perevt.jsonnet                               # dv2_split_move
wire-cell ... --tla-str vertex_method=dv2_snap -c pgrapher/experiment/sbnd/wct-pr-perevt.jsonnet  # snap
wire-cell ... --tla-str vertex_method=dv2_split -c pgrapher/experiment/sbnd/wct-pr-perevt.jsonnet # split
wire-cell ... --tla-str vertex_method=original -c pgrapher/experiment/sbnd/wct-pr-perevt.jsonnet  # old chain
```

| `vertex_method` | main neutrino vertex |
|---|---|
| `original` | the pre-DV2-MAX production chain: SCN net + re-rank, dual chain, kink/junction snaps, `improve_vertex`, main-vertex graph audit, conn3 stitch |
| `dv2_snap` | the DV2-MAX point, snapped to the nearest PR-graph vertex |
| `dv2_split` | the nearest PR-graph vertex or trajectory-fit point, whichever is nearer; a fit point splits its segment there |
| `dv2_split_move` (default) | the nearest PR-graph vertex if within 1.5 cm, else as `dv2_split`; then the vertex is moved onto the DV2-MAX point |

- The `dv2_*` methods run the DV2-MAX model and that placement only. None of the original vertex steps run:
  no SCN, re-rank, dual chain, kink/junction snap, `improve_vertex`, graph audit or conn3 stitch.
- Everything before the vertex step is unchanged, because it builds the model's input cloud. Everything after
  it is unchanged too: clustering, directions, showers, taggers and the output trees.
- The model is imported and loaded at configure. A missing package, a bad bundle or a sha256 mismatch fails the
  job at startup.
- A Python error during an event fails the job. `DV2_STRICT=0` instead prints the error on stderr and keeps the
  traditional vertex for that call.
- If DV2-MAX makes no prediction (fewer than 5 points or no candidate), the event keeps the traditional vertex, with
  the original refinement, and the log has a WARN line.
- Each placement is logged at INFO: `vertex_method <m>: main vertex = split|graph vertex at (x, y, z) cm, d cm
  from the DV2-MAX point ...`, plus one line per move for `dv2_split_move`.
- Known behaviour, by design:
  - candidates are every PR-graph vertex and fit point of the event, so the main cluster can switch, as in the
    legacy SCN path;
  - `dv2_snap` and `dv2_split` do not refit trajectories; the fits are those the vertex step saw;
  - `dv2_split` can split a track one fit spacing (about 0.6 cm) from an existing vertex, which leaves a short
    stub segment between the two;
  - a split long-muon segment stays the long muon: its two halves replace it.
- `dv2_bundle` (SBND jsonnet default `dv2max/sbnd/aurora_space/bundle.json`, resolved through `WIRECELL_PATH`)
  picks the model; `--tla-str dv2_bundle=dv2max/sbnd/r7_cloud/bundle.json` selects the cloud-only model. The C++
  default (configs that set no `dv2_bundle`) stays `r7_cloud`.
  `clus.jsonnet` holds the sha256 of each known bundle manifest (`bundle.json`, which carries the sha256 of every
  other bundle file). `TaggerCheckNeutrino` fails at configure if the manifest it loads has another sha256
  (`dv2_bundle_sha256`), and `Trun/dl_weights` records `<vertex_method>:<dv2_bundle>@sha256:<that hash>`.
  `dl_weights` stays the SCN net of `original`.
- `original` compiles to the same configuration as master (b7bd2a1a), byte for byte. The C++ default of
  `TaggerCheckNeutrino` stays `original`, so detectors other than SBND are unchanged.
- The dv2 methods need this directory and `DV2_Vertex.py` on `PYTHONPATH`, a Python with numpy, scipy,
  torch and lightgbm, and libpython preloaded, as the SCN net does.

## Entry point

```python
DV2_Vertex(weights, x, y, z, q, n_vertex_rows, trad, q_scale, q_offset, dtype="float32",
           space_x=None, space_y=None, space_z=None, space_q=None, space_cluster_id=None) -> bytes
```

| argument | meaning |
|---|---|
| `weights` | absolute path of `bundle.json` |
| `x, y, z, q` | bytes of the cloud floats: cm, raw reco frame; `q = dQ * q_scale + q_offset`; the first `n_vertex_rows` points are the PR-graph vertex rows |
| `trad` | bytes of 3 floats: the traditional vertex, NaN if none |
| `space_*` | only for a space-point bundle (below): bytes of the call's space points, x / y / z / q floats and int32 cluster ids; all five or none |

It returns bytes of `[x, y, z, score]` in the raw reco frame. It returns empty bytes when no prediction is
made: fewer than 5 points or no candidate. The caller then keeps the traditional vertex. Any other error raises;
`DV2_STRICT=0` instead reports it on stderr and returns empty bytes.

`DV2_Load(weights)` imports the model code and loads the bundle. `TaggerCheckNeutrino::configure` calls it, so a
broken environment fails at startup. It always raises on error.

The C++ wrappers (`WCPPyUtil::DV2_Load`, `DV2_Vertex`, and `SCN_Vertex` of the original method) start the
interpreter once per process and then release the GIL (`WCPPyUtil/Interpreter.h`); each call takes the GIL for
its own duration, so calls are safe from any thread (they run one at a time).

Models are loaded once per process and manifest path. The sha256 of the bundle files is checked at load
(`DV2_VERIFY=0` skips the check). Other environment variables:

| variable | effect |
|---|---|
| `DV2_THREADS` | torch / LightGBM threads, default 1, interop threads 1. torch's thread setting belongs to the whole process (SCN_Vertex of the original method also sets 1); `0` leaves torch's settings as the process chose them. The load report records the values before and after |
| `DV2_HOOK_LOG=<path>` | append one JSON line per call: n points, n_vertex_rows, `input_sha1`, prediction, score, source, stage and total ms. `input_sha1` (`DV2_Vertex.input_digest`) is a sha1 of the call's exact arguments, so an in-chain call can be matched with its standalone replay; `validate_r7.py` records it too |
| `DV2_DEBUG` | print that line to stdout as well |

## Pipeline

This is what `vertexer.Vertexer.run` does. It is the offline r7_cloud path for a test event of family f0.

1. **Event.** Built from the call as `r7/prep_r7_cache.py build_event` does:
   - `dQ = clip((q - q_offset) / q_scale, 0)`;
   - 3 cm connected components, renumbered by size, are the cluster ids;
   - the first `n_vertex_rows` points are the PR candidates, and `trad` is tag_nu;
   - every reco input is moved 0.31 cm toward the cathode: `x -= sign(x) * 0.31`.
2. **Features.** `sbndvtx.features.compute_features`, rounded to the cache's storage precision (point-level
   float16, the rest float32). The models and the ranker were trained on cached features.
3. **Pass 1.** 1 cm, view `id`, 4 grid shifts, pooled vote decoding.
4. **Pass 2.** 0.3 cm crops at up to 4 pass-1 modes, views `id` and `fy`, 4 shifts.
5. **Ranker.** The f0 variant's candidates (`cand`, `tag`, `clusend`, plus imaging points near the CNN anchors,
   plus the CNN modes) and the 382 union features. Then the LightGBM LambdaMART score, top-1 after 2 cm NMS, and
   the isotonic map of the top-1 score to P(< 1 cm).
6. **Output frame.** The pick is in the label frame. `DV2_Vertex` returns `x + sign(x) * 0.31`.

## Vendored files

All files below were copied from `R = /lus/eagle/projects/neutrinoGPU/snehadri/sbnd-vertexer`. R is not a git
repository, so the sources are identified by sha256 prefix and mtime (UTC). Every source predates the
r7_cloud run of 2026-10-11/12, except `prep_r7_cache.py`.

| file here | source in R | sha256[:16], mtime | changes |
|---|---|---|---|
| `features.py` | `common/sbndvtx/features/__init__.py` | 484147cf3fa1eede, 2026-09-24T22:45Z | none (flat package: relative imports resolve here) |
| `geometry.py` | `common/sbndvtx/features/geometry.py` | a9c66c6ebade3f9e, 2026-09-24T21:41Z | none |
| `local.py` | `common/sbndvtx/features/local.py` | a1d88aaee854a36f, 2026-09-24T22:08Z | none |
| `fitgraph.py` | `common/sbndvtx/features/fitgraph.py` | da2866ae286a425d, 2026-09-24T23:21Z | none |
| `tokens.py` | `common/sbndvtx/features/tokens.py` | 400513167140ed8a, 2026-09-24T22:21Z | none |
| `event.py` | `common/sbndvtx/event.py` | e2efb3d66e21e1e5, 2026-09-24T23:21Z | Event dataclass only; truth / cache I/O dropped |
| `mirror.py` | `common/sbndvtx/augment.py` | f107522b755f5aaf, 2026-09-27T02:10Z | `odd_columns`, `mirror_features`, `mirror_event` only; truth must be None |
| `baselines.py` | `common/sbndvtx/eval/baselines.py` | ce8903e98e1875da, 2026-09-24T23:26Z | pool and pool-feature functions only; bodies unchanged |
| `pool_ext.py` | `common/sbndvtx/eval/pool_ext.py` | f29db0624e95b987, 2026-09-25T14:47Z | `ExtPool`, `extension_points`, `expand_pool` only |
| `inputs.py` | `dv2max/dv2max/inputs.py` | 90bf0752130dd1e4, 2026-10-06T19:30Z | imports only |
| `decode.py` | `dv2max/dv2max/decode.py` | 88061943b6afb64d, 2026-09-27T06:39Z | none |
| `crop.py` | `dv2max/dv2max/crop.py` | bf601814116e676d, 2026-09-27T03:22Z | none |
| `sparse.py` | `dv2max/dv2max/sparse.py` | cf7b62287aeea027, 2026-09-27T22:16Z | spconv adapters and probe removed; torch backend only |
| `model.py` | `dv2max/dv2max/model.py` | 222f34111d8bcce5, 2026-09-27T05:15Z | `torch_scatter` replaced by `scatter.py` |
| `scatter.py` | new | | pure-torch `scatter_mean` / `scatter_max` with torch_scatter semantics (empty groups give 0) |
| `infer.py` | `dv2max/dv2max/infer.py` | b4c2ad8ef5939d8f, 2026-09-27T06:39Z | DataLoader dataset dropped; collate from `batch.py`, mirror from `mirror.py` |
| `batch.py` | `dv2max/dv2max/data.py` | e81bf481c3fa67cb, 2026-09-30T16:09Z | `collate`, `to_device` only |
| `ranker.py` | `dv2max/dv2max/ranker.py` | 1a6b4c8cc2987f18, 2026-10-06T19:30Z | CNN candidates/features, `variant_pool`, `event_union` only; `event_union` returns a dict, no labels |
| `vertexer.py` | `r7/prep_r7_cache.py` (`xshift`, `cluster_ids`, `build_event`) | 7f1a087855e2da49, 2026-10-19T15:48Z | new glue: bundle loading, cache precision, LightGBM and isotonic at inference |

Notes on the table:

- **`prep_r7_cache.py` mtime.** Its `--passes` option was added after the cache was built. The validation below
  shows that the event logic is unchanged.
- **`threadpoolctl`.** It is optional here. The original used it only to limit BLAS threads.

`make_bundle.py` is a build tool (it needs pyyaml + torch, e.g. the sbnd-vertexer venv). It writes the bundle from
`R/dv2max/runs/r7_cloud` and the ranker re-fit. `validate_r7.py` is the parity test.

## Model bundle (wire-cell-data `dv2max/sbnd/r7_cloud/`)

| file | contents |
|---|---|
| `pass1.pt`, `pass2.pt` | state dicts of `runs/r7_cloud/models/pass{1,2}/f0/seed0/final.pt` (SWA weights, fp32, all keys; ~80 MB each) |
| `ranker_fold0.txt` | LightGBM text model of the fold-0 LambdaMART |
| `isotonic_fold0.json` | isotonic thresholds of the fold-0 calibrator |
| `bundle.json` | format, model name, family f0, trained-on folds, x shift, every inference setting (resolved config sections `inputs`, `pass1`, `pass2`, `ranker`, `data.min_points`), the 382 ranker feature names in order, the re-fit report, sha256 of every file |

### Ranker re-fit

The offline run saved only the predictions of the ranker and the isotonic map, not the models. Both were re-fit
from the saved union table (`runs/r7_cloud/ranker/union`) with the repo's own code and settings:
`dv2max.ranker.load_outer(cfg, 0)`, then `_fit`, `calibrator` and `configs/r7_cloud.yaml`. The ranker settings
hash is 7e4ef1400541, the same as offline. Training: 20,320 events / 1,649,599 rows, 382 features.

**Why the re-fit is not bit-exact.**
- LightGBM 4.7 picks col-wise or row-wise histograms by timing both at the start of every fit.
- Its multithreaded training is also not run-to-run reproducible. On a busy node (sophia-gpu-07), two identical
  fits with the mode pinned differed by up to 1.1 in raw score.

**What was kept.** The fit was run on an idle node (sophia-gpu-01) with the mode pinned per fit:
- the early-stopping fit in row-wise mode reproduces the offline NDCG@1 curve exactly (600 values, best
  iteration 560);
- the refit with 560 trees in col-wise mode is the saved model.

**Check of the saved model.** On the offline fold-0 test rows (5,030 events, 415,757 rows):

| quantity | result |
|---|---|
| top-1 picks | 5030 / 5030 identical (distance 0) |
| order of the offline top-10 list | identical in 100 % of events |
| raw score vs offline top-10 hypotheses, \|diff\| | median 8e-7, p99 1.1e-5, max 0.12 |
| calibrated score, \|diff\| | median 0, max 3.3e-5 |

The isotonic map was re-fit the same way. Both calibration halves stop at the same iteration in either mode
(126 and 419). Every mode combination gives the same map to within 3.3e-5.

## Space-point bundles (format `dv2max-bundle-v2`, e.g. `aurora_space`)

A v2 bundle's model also takes the call's **space points**: the Bee clustering-layer points of every in-scope
imaging cluster that has a point within `space.touch_cm` (1 cm) of a cloud point. This is the r7_space input
(sbnd-vertexer `r7/export_r7_space.py`), recorded in production as `T_dlvtx_space` (HaiwangYu's
SbndPrMagnifyTrackingVisitor `dl_vtx_dump_space`), on which `aurora_space` was trained.

- **Collected in C++** by `PatternAlgorithms::dv2_collect_space` at the DL call, from the call's grouping, with
  `BeeClusterPoints` (the same functions as the `T_dlvtx_space` writer): clusters in increasing cluster id, each
  cluster's points in Bee order, Bee charge. Scope = the PR Bee clustering layer's (`TaggerCheckNeutrino`
  `dv2_space_pcname` / `dv2_space_coords` / `dv2_space_filter`, set by `clus.jsonnet`). `TaggerCheckNeutrino`
  learns that the bundle takes space points, and the touch distance, from the `space` section of the `DV2_Load`
  report; a cloud-only bundle (v1, r7_cloud) gets no space points and the C++ collects none.
- **Added in Python** by `vertexer.add_space`, exactly as the training cache was built
  (`r7/prep_r7_space_cache.add_space`; `aurora/prep_aurora_cache.add_space`): space points first (x-shifted, raw
  q, ids renumbered by cluster size), then the cloud (`dQ * space.q_ratio`, the id of the nearest space point within
  `space.match_cm`, else 3 cm components numbered after the space clusters). A call with no touching cluster is
  the cloud-only event with dQ unscaled (the training rule).
- **Checks.** A v2 bundle called without space points, a v1 bundle called with them, or partial / ragged space
  arguments raise. `input_digest` covers the space bytes; cloud-only digests are unchanged.
- **The bundle** adds a `space` section (`touch_cm`, `match_cm`, `q_ratio`, `link_cm`, the cache it was measured
  on, the rule) and uses that cache's x shift (`aurora_space`: 0.29 cm, `q_ratio` 0.2474433061809911).
  `make_bundle.py --model-dir <saved ranker dir> --space-cache <cache_info.json>` builds it.
- **Validation:** `validate_space.py` (validate_r7's steps; extract also reads the `.space.h5` shards).

### Space-point validation results (`aurora_space`, 2026-11-14)

**Hold-out replay.** All 30,000 calls of the `aurora_space` hold-out (15,000 labeled events never trained on:
7,397 BNB + 7,603 nue CC; pass 0 and pass 1) were run through `DV2_Vertex` with their space-point bytes, in
pyenv-dl at 1 thread per worker on Polaris. The offline predictions came from the training pipeline (GPU).

| quantity | value |
|---|---|
| picks < 0.01 / 0.05 / 0.1 / 1 cm from the offline pick | 0.927 / 0.961 / 0.987 / 0.999 |
| median / p99 distance | 3.4e-4 / 0.12 cm |
| calibrated score \|diff\|: median / p90 / p99 / max | 0 / 0.0011 / 0.035 / 0.20 |

| calls | n | < 1 cm DV2_Vertex | < 1 cm offline | only DV2_Vertex | only offline |
|---|---|---|---|---|---|
| pass 0 | 15000 | 0.9293 | 0.9288 | 13 | 5 |
| pass 1 | 15000 | 0.9287 | 0.9284 | 11 | 6 |
| all | 30000 | 0.9290 | 0.9286 | 24 | 11 |

The 385 calls (1.3 %) more than 0.1 cm from the offline pick are near-tie top-1 flips (median 0.17 cm). Every
call gave a prediction.

**CPU time per call** (1 thread, 1,500 random hold-out calls, 32 workers on 32 cores): median 1.8 s, mean 3.2 s,
p90 7.4 s, max 87 s; stage medians: pass 1 0.34 s, pass 2 1.2 s, features + ranker 0.16 s. In the chain run
below, on the same 339 calls and the same node: median 1.7 s vs 0.53 s with `r7_cloud` (3.2x per call), mean
2.3 s vs 0.56 s.

**In the chain** (the 372 `r7/inchain` fold-0 events, build `@vaurora`):
- `original` and `dv2_split_move` with `r7_cloud` give the same `tracking-pr.root` as the PR #543 build on every
  finished event (336 and 337; the same 36 and 35 failed, the known `ParticleInfo` crash).
- With `aurora_space` the call collects the space points; 339 calls in 337 events, all with space points.
- **Same points as the writer.** A test-only build that also writes `T_dlvtx_space` (HaiwangYu's
  sbnd-dlvtx-space-44 merged in) was run on the same events. For all 339 calls the `DV2_Vertex` input digest
  (cloud + space bytes) equals the digest of the dumped `T_dlvtx_call` / `T_dlvtx_cloud` / `T_dlvtx_space` rows
  (1,331,448 space points), and every listed cluster passes the 1 cm touch test. Its DL inputs and outputs are
  identical to `@vaurora`'s on all 337 events.
- The 339 dumped calls replayed outside the chain give the in-chain float32 output bit for bit in 334 calls; the
  other 5 agree to 1.6e-5 cm and 3e-8 in score (the call-order effect, see Repeatability).

## Tests

Automated tests (run with the WCT test suite; each skips with the reason when what it needs is missing):

| file | what it tests | needs |
|---|---|---|
| `pyutil/test/test_dv2max.py` | the algorithm, 33 tests: call -> Event (x shift and back, charge decoding, PR-vertex candidate rows, traditional vertex, 3 cm clusters), NMS, bundle loading (contents, sha256 check, format, feature names), the isotonic calibration, the `DV2_Vertex` / `DV2_Load` API (input checks, `DV2_STRICT`, no prediction, output frame, repeatability, `DV2_HOOK_LOG`, `DV2_THREADS`), per-stage invariants on held-out calls, the golden replay, and run-to-run repeatability in fresh processes | numpy, scipy, torch, lightgbm; the bundle on `WIRECELL_PATH` (or `DV2MAX_BUNDLE`) for 24 of them |
| `pyutil/test/doctest_dv2_vertex.cxx` (`wcdoctest-pyutil`) | the C++ wrapper against a fake module (argument marshalling, empty = no prediction, every Python error reaches C++ with its type and message), and 5 golden calls through the wrapper with the real model | the golden case: the bundle and the Python packages |
| `clus/test/doctest_vertex_method.cxx` | method names, the payload check behind the fallback (empty, wrong size, NaN), the configure-time failures (bad method, missing bundle, failed import) | -- |
| `pyutil/test/test_dv2max_space.py` | space-point bundles, 20 tests: `add_space` against a verbatim copy of the training rule (order, charge scale, cluster ids by size, unmatched-point components, the no-space case), the input digest (cloud-only digests unchanged), the v2 manifest checks, the `DV2_Vertex` space API (load report, missing / partial / ragged space arguments, the empty set, a v1 bundle given space points), and the golden replay | numpy, scipy, torch, lightgbm; the `aurora_space` bundle for the API and golden tests |
| `pyutil/test/doctest_dv2_vertex_space.cxx` (`wcdoctest-pyutil`) | the space overload of the C++ wrapper against a fake module: the five arrays in order (ids as int32), the empty set, no space arguments from the cloud-only overload, size checks before Python is called, a Python error's message | -- |
| `clus/test/doctest_dv2_space.cxx` (`wcdoctest-clus`) | reading the touch distance from the `DV2_Load` report (and its errors), and the touch test: a cluster is taken whole if any point is within `touch_cm` (inclusive), across cell edges | -- |
| `clus/test/doctest_pattern_recognition.cxx`, suite "dv2 placements [B]" | `dv2_place_point` and `dv2_place_main_vertex` on the fixture [B] PR graph: snap, split at a fit point (incl. the long-muon sets), split_move reuse within 1.5 cm, the move onto P and when it declines | fixture [B] |

The golden fixture `pyutil/test/data/dv2max_golden.json` holds 26 held-out fold-0 calls (BNB and nue, both
passes, three cloud sizes, two without a traditional vertex) and 3 synthetic no-prediction calls, with the answer
of the code at the time; tolerances 0.01 cm and 0.001 in score.  Only calls whose top-1 candidate wins by a
ranker-score margin of at least 0.5 are used, so the answers do not depend on floating-point details.  Regenerate
it with `pyutil/test/make_dv2max_golden.py` (see its header) when the model or the inference changes on purpose.
The space-point fixture `pyutil/test/data/dv2max_space_golden.json` holds 10 `aurora_space` hold-out calls (BNB
and nue, both passes, two without a traditional vertex; 50-329 cloud and 330-1,804 space points), selected the
same way (margin >= 0.5), tolerances 0.001 cm and 0.001 in score; made by `pyutil/test/make_dv2max_space_golden.py`.
Every hold-out call has space points, so the no-space case is covered by the unit tests only.
Accuracy itself (< 1 cm efficiency on all 5,038 held-out calls) is too slow for the test suite: that is
`validate_r7.py` below.

## Repeatability

- Within one process the same call gives the same bytes (`test_dv2max.py`).
- Across processes, which is how wire-cell runs it (one process per job): the same calls in the same order give
  the same bytes; in another order they agree to about 1e-5 cm and 1e-7 in score (`TestRepeat`; the golden
  replay saw at most 1.4e-5 cm and 3e-8).  CPU torch at 1 thread, fp32.
- In the chain (2026-11-06): `dv2_snap`, `dv2_split` and `dv2_split_move` rerun on all 372 fold-0 events of
  `r7/inchain` with later builds gave the same `tracking-pr.root`, every tree and branch, as the 2026-10-30 runs
  (337 events finished; the same 35 failed the same way, a known `ParticleInfo` energy-momentum crash also seen with
  `original`).
- `validate_r7.py` on all 5,038 held-out calls (2026-11-06) reproduced the 2026-10-29 results: every accuracy and
  agreement number identical.

## How to validate (`validate_r7.py`)

```bash
python validate_r7.py extract --out DIR                 # needs h5py: the r7 export -> DIR/calls.npz (fold 0)
python validate_r7.py run --out DIR --bundle <wire-cell-data>/dv2max/sbnd/r7_cloud/bundle.json --workers 16 --results results_x.jsonl
python validate_r7.py report --out DIR --results results_x.jsonl      # -> DIR/report_x.json / .md
```

- `run` needs only the runtime environment. It calls `DV2_Vertex.DV2_Vertex` with the exact C++ argument bytes
  built from the raw dumped cloud floats.
- `report` moves each result back to the label frame. It then compares with `pred_dv2max.csv`: distance to the
  offline pick, score agreement, < 1 cm efficiency vs `label_xyz` for both (same calls, per pass), top-1 flips,
  and CPU time per call.

## Validation results (2026-10-29)

All 5,038 fold-0 calls of `pred_dv2max.csv` (2,519 pass 0 + 2,519 pass 1) were run through `DV2_Vertex` on a
compute node, 16 workers at 1 thread each. The offline run was on GPU, with spconv and fp16 autocast. This run is
on CPU, with the torch backend in fp32.

**Distance to the offline pick**

| quantity | value |
|---|---|
| fraction < 0.01 / 0.05 / 0.1 / 1 cm | 0.930 / 0.949 / 0.968 / 0.995 |
| median | 2.5e-4 cm |
| p99 | 0.30 cm |

The 8 calls without a prediction (fewer than 5 points) are exactly the 8 offline tag_nu fallbacks.

**Score agreement**

| score | median \|diff\| | p90 | p99 | max |
|---|---|---|---|---|
| calibrated | 0 | 0.0012 | 0.054 | 0.15 |
| raw | 1e-6 | 0.055 | 0.24 | 0.73 |

**< 1 cm efficiency vs `label_xyz` (same calls)**

| calls | n | DV2_Vertex | offline | only DV2_Vertex | only offline |
|---|---|---|---|---|---|
| pass 0 | 2519 | 0.8495 | 0.8484 | 8 | 5 |
| pass 1 | 2519 | 0.8547 | 0.8535 | 4 | 1 |
| all | 5038 | 0.8521 | 0.8509 | 12 | 6 |

**The 161 calls (3.2 %) more than 0.1 cm from the offline pick**
- 137 of them are within 1 cm; 24 (0.5 %) are 1 cm or more.
- They are near-ties. The offline pick is always among the deployed candidates: within 0.05 cm in 128 calls,
  moved 0.05-1 cm in 33 (the CNN mode moved).
- In the deployed scores, the offline pick is a median 0.07 below the top-1 (10-90 %: 0.006-0.24).
- 120 of the 161 pick a different CNN hypothesis.
- Only 18 change hit status (12 gained, 6 lost), so the efficiency does not change beyond noise.

**Intermediate parity**
- The event and features built from the raw call are bit-identical to the cache rows.
- Given the stored CNN records, the ranker candidates and all 382 features are bit-identical to the offline union
  rows: 400 / 400 events in the sbnd-vertexer venv.
- CNN records vs the stores (CPU fp32 vs GPU fp16): top pass-1 mode |diff| median 3e-4 cm, best pass-2 mode
  3e-4 cm, same number of crops.

**In-chain environment (pyenv-dl: Python 3.14, numpy 2.5.3, scipy 1.18.1, torch 2.9.1+cpu, lightgbm 4.7.0)**
- The picks are identical to the sbnd-vertexer venv run on all 5,038 calls. Calibrated score max \|diff\| is
  0.054, so the numbers above hold unchanged.
- The only intermediate difference: `clus_dend` / `clus_frac` are 0 vs ~1e-16 for candidates at a cluster PCA
  extreme. This is BLAS noise (OpenBLAS 0.3.30 vs 0.3.34). It changed no pick.

**CPU time per call** (pyenv-dl, 1 thread, 16 workers): median 0.63 s, p90 0.93 s, max 1.97 s. Stage medians:

| stage | ms |
|---|---|
| features | 6 |
| pass 1 | 226 |
| pass 2 | 374 |
| ranker rows | 17 |
| ranker | 5 |

The model load is a one-time cost per process. Imports dominate it (see below).

## First-call cost (one wire-cell process per event)

Measured in pyenv-dl at 1 thread, fresh processes on a quiet compute node with a warm file cache (mean of 3). The
`DV2_HOOK_LOG` record of the call that loads the models carries this breakdown as `first_call_ms`.

| step | before | after |
|---|---|---|
| `import numpy` (DV2_Vertex import) | 116 ms | 116 ms |
| `import torch` + thread setup | 1.56 s | 1.56 s |
| `import dv2max_vertexer` (scipy.spatial / sparse, package) | 360 ms | 360 ms |
| model load | 859 ms | 640 ms |
| first call | 670 ms | 670 ms |

Model load, in detail:

| part | before | after |
|---|---|---|
| sha256 check | separate pass, 158 ms | done on the bytes `torch.load` reads (read + hash + deserialise: 230 ms per model) |
| DV2Net build | 315 ms of random init | 50 ms (`skip_param_init`) |
| `torch.load` | 263 ms | (in the read above) |
| `load_state_dict` | 38 ms | 37 ms (`assign=True`) |
| lightgbm import + Booster | 87 ms | 92 ms |

The first call takes 670 ms; a second call takes 615 ms.

Outputs are bit-identical before and after on 100 fold-0 calls.

- **Cold file cache** (first process on a node): the same steps take ~18 s, mostly imports read from Lustre:
  numpy 2.3 s, torch 12 s, scipy 1.7 s.
- **Login node**: `import torch` alone ranged from 3 to 64 s.

So the first-call time is set by where the Python environment lives, not by DV2-MAX. Options:
- stage pyenv-dl on node-local storage;
- use a squashfs image in the container;
- run several events per process.

