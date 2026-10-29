# dv2max_vertexer: DV2-MAX (r7_cloud) inference for the in-chain DL-vertex step

Inference-only copy of the DV2-MAX vertexer, called by WCT through `pyutil/python/DV2_Vertex.py`
(`WCPPyUtil::DV2_Vertex`). One call = one neutrino candidate's DL-vertex cloud. The result is the DV2-MAX
top-1 vertex in the raw reco frame plus a calibrated P(< 1 cm).

- Runtime dependencies: numpy, scipy, torch (CPU), lightgbm. No spconv, torch_scatter, pandas, h5py or sklearn.
- Nothing is imported from sbnd-vertexer.
- Model bundle: `cfg/dv2max/r7_cloud/` (`bundle.json` and the files it lists, each with its sha256).

## Choosing the vertexing method

One option, `vertex_method`, on every SBND PR job: `wct-pr-perevt.jsonnet`, `wct-pr.jsonnet`, the LArSoft
library `wcls-img-clus-matching-pr-lib.jsonnet`, or `pr()` in `clus.jsonnet`:

```
wire-cell ... --tla-str vertex_method=dv2_split_move -c pgrapher/experiment/sbnd/wct-pr-perevt.jsonnet
```

| `vertex_method` | main neutrino vertex |
|---|---|
| `original` (default) | the production chain: SCN net + re-rank, dual chain, kink/junction snaps, `improve_vertex`, main-vertex graph audit, conn3 stitch |
| `dv2_snap` | the DV2-MAX point, snapped to the nearest PR-graph vertex |
| `dv2_split` | the nearest PR-graph vertex or trajectory-fit point, whichever is nearer; a fit point splits its segment there |
| `dv2_split_move` | the nearest PR-graph vertex if within 1.5 cm, else as `dv2_split`; then the vertex is moved onto the DV2-MAX point |

- The `dv2_*` methods run the DV2-MAX model and that placement only. None of the original vertex steps run:
  no SCN, re-rank, dual chain, kink/junction snap, `improve_vertex`, graph audit or conn3 stitch.
- Everything before the vertex step is unchanged, because it builds the model's input cloud. Everything after
  it is unchanged too: clustering, directions, showers, taggers and the output trees.
- If DV2-MAX makes no prediction, the event keeps the traditional vertex, with the original refinement, and the
  log has a WARN line. Set `DV2_STRICT=1` in validation runs to make a Python error fail the job instead.
- Each placement is logged at INFO: `vertex_method <m>: main vertex = split|graph vertex at (x, y, z) cm, d cm
  from the DV2-MAX point ...`, plus one line per move for `dv2_split_move`.
- Known behaviour, by design:
  - candidates are every PR-graph vertex and fit point of the event, so the main cluster can switch, as in the
    legacy SCN path;
  - `dv2_snap` and `dv2_split` do not refit trajectories; the fits are those the vertex step saw;
  - `dv2_split` can split a track one fit spacing (about 0.6 cm) from an existing vertex, which leaves a short
    stub segment between the two;
  - a split long-muon segment stays the long muon: its two halves replace it.
- `dv2_bundle` (default `dv2max/r7_cloud/bundle.json`, resolved through `WIRECELL_PATH`) picks the model.
  `dl_weights` stays the SCN net of `original`.
- `original` compiles to the same configuration as master, byte for byte.
- The dv2 methods need this directory and `DV2_Vertex.py` on `PYTHONPATH`, a Python with numpy, scipy,
  torch and lightgbm, and libpython preloaded, as the SCN net does.

## Entry point

```python
DV2_Vertex(weights, x, y, z, q, n_vertex_rows, trad, q_scale, q_offset, dtype="float32") -> bytes
```

| argument | meaning |
|---|---|
| `weights` | absolute path of `bundle.json` |
| `x, y, z, q` | bytes of the cloud floats: cm, raw reco frame; `q = dQ * q_scale + q_offset`; the first `n_vertex_rows` points are the PR-graph vertex rows |
| `trad` | bytes of 3 floats: the traditional vertex, NaN if none |

It returns bytes of `[x, y, z, score]` in the raw reco frame. It returns empty bytes when no prediction is
made: fewer than 5 points, no candidate, or an internal error (reported on stderr; set `DV2_STRICT=1` to raise
instead). The caller then keeps the traditional vertex.

Models are loaded once per process and manifest path. The sha256 of the bundle files is checked at load
(`DV2_VERIFY=0` skips the check). Other environment variables:

| variable | effect |
|---|---|
| `DV2_THREADS` | torch / LightGBM threads, default 1; interop threads are always 1 |
| `DV2_HOOK_LOG=<path>` | append one JSON line per call: n points, n_vertex_rows, prediction, score, source, stage and total ms |
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

## Model bundle (`cfg/dv2max/r7_cloud/`)

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

## How to validate (`validate_r7.py`)

```bash
python validate_r7.py extract --out DIR                 # needs h5py: the r7 export -> DIR/calls.npz (fold 0)
python validate_r7.py run --out DIR --bundle <W>/cfg/dv2max/r7_cloud/bundle.json --workers 16 --results results_x.jsonl
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

