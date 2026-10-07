# Evaluating DNN-ROI models with `evaluate.snakefile`

`evaluate.snakefile` runs trained models over existing h5 samples and makes
metrics and plots. It does not make the h5 files: use `Snakefile` in this
directory for that (LArSoft, depo generation, wire-cell sigproc, tar to h5).

How you run it depends on the model:

- **xvunet** takes the rec/tru h5 files for all three views (`-0`, `-1`, `-2`)
  at once and outputs one image covering every channel of the APA.
- **dnnroi_custom** runs on a single view. Each view needs its own model, and
  the output covers only that view's channels.

To compare the two, the three per-view dnnroi outputs are concatenated into the
xvunet channel layout (the **merged** style below).

## Setup

You need a Python environment with:

- `snakemake` 9.x (the file uses `pathvars`, e.g. `<results>`);
- `torch`, `numpy`, `matplotlib`, `h5py`;
- `wcpy` (wire-cell-python) on `PATH`, since the rules call `wcpy dnn run_one` / `run_n`.

A wire-cell-python checkout's `.venv` has all of these:

```
export PATH=/path/to/wire-cell-python/.venv/bin:$PATH
```

Each run looks like this:

```
snakemake -s /path/to/evaluate.snakefile \
    --configfile my-eval.yaml --directory WORKDIR \
    --cores 4 [--resources gpu=1] -- TARGET ...
```

- Run `snakemake -s evaluate.snakefile help` for a summary of the config keys
  and targets.
- Put the workdir on a disk with room. A model output `.pt` is about 46 MB per
  sample, so a full cosmics comparison is tens of GB per model. Home areas with
  quotas fill up.
- Every inference job uses several GB of host RAM.
  - On CPU, keep `--cores` low; 4 parallel jobs exhausted 30 GB.
  - On GPU, set `device: gpu`, pick a card with `CUDA_VISIBLE_DEVICES`, and pass
    `--resources gpu=1` so jobs take turns. xvunet needs more than 8 GB of GPU
    memory on cosmics.
- Each wcpy call spends most of its time starting up, so cosmics are run one
  directory (all events) per process.

## Samples

| Sample | `{sample}` name | Inputs |
|---|---|---|
| Line depos | `line_Pplane_t1-A-t2-B` | `input_dir/linedepos-pdhd-Pplane-t1_A-t2_B-g4-{rec,tru}-{0,1,2}.h5` and `...-g4-trio.h5` |
| Cosmics | `cosmics_N_M_evE` | `cosmics_dir/N_M/cosmics_*-g4-{rec,tru}-{0,1,2}.h5`, event group `E`, and `...-g4-trio.h5` |

- The line angles per plane are fixed in `ANGLES` near the top of the snakefile.
- Cosmics use every `N_M` directory under `cosmics_dir`, or the list in
  `cosmics_dirs`, with `cosmics_events` events each (default 10).
- **Always set `cosmics_dirs` to directories held out of training.** The
  campaign cfgs train on `0_0`–`1999_0`, so evaluate on e.g. `2000_0`–`2009_0`.
  Left unset, every directory is used: thousands of them, training ones
  included.
- Quote the names, e.g. `cosmics_dirs: ["2000_0", "2001_0"]`. YAML reads a bare
  `2000_0` as the integer 20000.

## Config shared by both

```yaml
device: cpu                 # or gpu
input_dir: /path/to/linedepos/test        # line-sample h5 files
cosmics_dir: /path/to/cosmics_trios       # holds N_M/ subdirectories
cosmics_events: 10
pathvars:
  results: results_MODELTAG               # keeps each model's outputs apart
```

- `<results>` in a target means the `pathvars.results` value. It defaults to `results`.
- The model `cfg` files are the training cfgs, e.g.
  `wire-cell-python/wirecell/dnn/cfg/spng_campaign/`. Only their `[model]` and
  `[run_one_dataset]` sections matter here, because the input files come from
  the rules, not from `[train] files`.

## xvunet

```yaml
app: xvunet
model: /path/to/training_results_0_NNN.pt
cfg: /path/to/jcal-uniter-....cfg
```

Each target below is listed with what it produces.

**Line samples**

| Target | Output |
|---|---|
| `all_eff_pur_xvu` | pixel efficiency/purity vs angle, every plane (`<results>/all_eff_pur_Pplane-xvu.{png,npz}`) |
| `all_plots_P_xvu` | the above for plane `P`, plus feature/prediction/label displays (`xvu-line-displays_Pplane.tar`) |
| `<results>/threshold-T-xvu-all_{roi,pixel}_eff_pur_Pplane.png` | ROI or pixel inefficiency/impurity vs angle at threshold `T` |
| `all_roi_tables`, `<results>/xvu-{true,reco}_roi_table_info_Pplane-thresh-0.5.npz` | ROI tables and ROI lengths |
| `<results>/threshold-T-xvu-trio-merged-crossplane_Pplane_high_end.png` | trio cross-plane rates per angle |

**Cosmics**

| Target | Output |
|---|---|
| `<results>/threshold-T-xvu-trio-merged-crossplane_cosmics.{npz,png}` | trio rates for every event |

**Threshold scans and F-beta**

These use `run_n` over a list of files, so they need two extra config keys:

- `paths`: a shell glob of the rec/tru files. Exclude the trio file, e.g.
  `.../2000_0/cosmics_1_1-g4-*-?.h5`.
- `nentries`: how many entries to run.

| Target | Output |
|---|---|
| `aggregated_scan_{roi,pixel}_plane_P_fbetas.png` | F-beta vs threshold, summed over the entries |

The scan runs from `threshold_start` (0.1) to 1 in steps of `threshold_step`
(0.1). Threshold 0 is left out on purpose: every pixel passes there, so each
channel becomes one ROI that contains some truth, and the ROI scan reports a
meaningless R = P = 1.

```
snakemake ... --config 'paths=/path/2000_0/cosmics_1_1-g4-*-?.h5' nentries=10 \
    -- aggregated_scan_roi_plane_u_fbetas.png
```

The xvunet default `[run_one_dataset]` takes the sample ID from the file name.
Every cosmics directory holds the same file names, so `paths` should cover one
directory, unless the cfg's `[run_one_dataset]` sets `rec_file_res`/`tru_file_res`
that take the ID from the directory.

**Epochs**

Set `cpt_dir`, `proc` and `cluster` to point at `checkpoint_{proc}_{cluster}_{epoch}.pt`.

| Target | Output |
|---|---|
| `all_epochs_xvu` | per-plane eff/pur for every checkpoint |
| `<results>/eff_pur_epochs_Pplane-xvu.pdf` | eff/pur vs epoch |

## dnnroi_custom

### One per-view model

Set `app`, `model` and `cfg` to the model of the view you want:

```yaml
app: dnnroi_custom
model: /path/to/training_results_0_NNNN.pt
cfg: /path/to/jcal-dnnroi-...-Pplane-pdhd.cfg
test_tpcs: [0]
```

For line samples, `run_one` feeds only the sample plane's own view (`P` →
view 0/1/2). The output covers that view's channels only (800, or 480 for w),
starting at row 0. So:

- run only the plane that matches the model;
- keep `test_tpcs: [0]`, which scores rows 0–800 of the output.

| Target | Output |
|---|---|
| `all_plots_P` | `<results>/all_eff_pur_Pplane.{png,npz}` and the displays in `line-displays_Pplane.tar` |
| `<results>/eff_pur_epochs_Pplane.pdf` | with `cpt_dir`/`proc`/`cluster` set |

For threshold scans, point `paths` at that view's files only (e.g.
`.../2000_0/cosmics_1_1-g4-*-0.h5`) and use plane `all_planes` in the target,
e.g. `aggregated_scan_roi_plane_all_planes_fbetas.png`. `all_planes` covers
every row, so it works whatever the width of the single-view output. Plane
`u`/`v`/`w` would cut at the xvunet offsets, which do not apply here.

Not valid for single-view outputs:

- the `xvu-...` targets;
- `noxvu` trio targets, because trio channel numbers are APA-global and the
  output is not.

### All three views, merged (to compare with xvunet)

List the three per-view models. `run_view` runs each on its own view, and
`merge_views` concatenates them into the xvunet layout (`dnnroi_dir/merged-test_{sample}.pt`):

```yaml
view_app: dnnroi_custom
view_models:                       # views 0 (u), 1 (v), 2 (w)
  - {model: /path/to/u.pt, cfg: /path/to/jcal-dnnroi-...-uplane-pdhd.cfg}
  - {model: /path/to/v.pt, cfg: /path/to/jcal-dnnroi-...-vplane-pdhd.cfg}
  - {model: /path/to/w.pt, cfg: /path/to/jcal-dnnroi-...-wplane-pdhd.cfg}
view_widths: [800, 800, 960]
dnnroi_dir: dnnroi
```

- The PDHD w model covers the first 480 of the 960 w channels. `merge_views`
  zero-pads it to `view_widths`; the second w face has no simulated signal.
- The merged outputs do not depend on the xvunet model, so they are reused when
  you compare several xvunet models in the same workdir.

| Target | Output |
|---|---|
| `<results>/threshold-T-merged-trio-merged-crossplane_{Pplane_high_end,cosmics}.{npz,png}` | trio rates |
| `<results>/threshold-T-merged-trio-crossplane-{sample}.npz` | one sample |

## Comparing xvunet with the merged dnnroi models

Put both the xvunet keys (`app`, `model`, `cfg`) and `view_models` in the same
config, then run:

```
snakemake ... -- all_trio_compare
```

This writes `<results>/threshold-0.5-trio-compare_{u,v,w}plane_high_end.{png,npz}`,
plus `_cosmics` when `cosmics_dir` is set. Each one shows, for every
`compare_styles` entry (default `[xvu, merged]`):

- the fraction of true trios predicted in each combination of views, pooled
  over the set and weighted by the number of valid trios;
- the per-sample U & V & W fraction.

A true trio is a (u, v, w, tick) with truth in all three views. The
combinations are:

| Combo | Predicted in |
|---|---|
| `uvw` | all three views |
| `uv` | U and V, not W (and likewise `vw`, `wu`) |
| `just_u` | only U (and likewise `just_v`, `just_w`) |

Start with a small target to check the setup before the full set, e.g. one
sample:

```
snakemake ... -- results_TAG/threshold-0.5-{xvu,merged}-trio-crossplane-cosmics_2000_0_ev0.npz
```

### Several xvunet models

List them under `xvu_models` with a tag each (letters, digits, `_`). Every tag
becomes a style next to `merged`, with outputs in `<results>/TAG-test_{sample}.pt`:

```yaml
xvu_models:
  trio_c4:   {model: /path/to/training_results_0_634.pt, cfg: /path/to/jcal-uniter-trio-chunk4-band0.cfg}
  uniter_c4: {model: /path/to/training_results_0_636.pt, cfg: /path/to/jcal-uniter-...-legacy-attn-chunk4-band0.cfg}
compare_styles: [trio_c4, uniter_c4, merged]
op_compare:
  - [trio_c4, perplane]
  - [trio_c4, single]
  - [uniter_c4, perplane]
  - [uniter_c4, single]
  - [merged, perplane]
```

- The top-level `model`/`cfg` still define the model tagged `xvu`, which the
  `xvu-...` line-sample targets use.
- The merged dnnroi outputs are computed once and shared by all of them.

### At each model's best-F1 threshold

A fixed threshold is not a fair comparison: the models are calibrated
differently, and every trio rate rises towards 1 as the threshold drops. Instead,
compare each model at its own operating point:

```
snakemake ... -- all_trio_compare_op
```

How it works:

1. `cosmics_dirs` is split into a **tune** set (`cosmics_tune_dirs`, default the
   first half) and a held-out **eval** set (`cosmics_eval_dirs`, default the rest).
2. For each style, every tune sample is scanned over thresholds (ROI and pixel)
   in the ranges `u`, `v`, `w0` (the w face with signal) and `uvw0` (all three
   together). The scans are summed over the tune set, and the threshold with
   the best F1 is taken.
3. The trio rates on the eval set and on the line samples are measured at those
   thresholds. There are two kinds of operating point:
   - `perplane`: each view uses its own best threshold (u, v, w0);
   - `single`: one threshold, the best over `uvw0`, for all three views.

`op_compare` lists the (style, mode) pairs to put side by side. The default is
xvunet `perplane`, xvunet `single`, and merged dnnroi `perplane`.

| Output | Contents |
|---|---|
| `<results>/op-trio-report.pdf` | everything below in one PDF (see next) |
| `<results>/op-{roi,pixel}F1-trio-compare_{cosmics_eval,Pplane_high_end}.{png,npz}` | the comparison for one set; the legend shows each entry's thresholds |
| `<results>/op-{roi,pixel}F1_{perplane,single}-{style}.json` | the chosen thresholds and their F1 |

`op-trio-report.pdf` covers the sets in `report_sets`. By default that is
`cosmics_eval` plus the high-end-angle line samples of each plane
(`uplane_high_end`, `vplane_high_end`, `wplane_high_end`). Its pages:

1. Setup: tune and eval sets, threshold range, and each model's checkpoint and cfg.
2. Summary: the pooled U & V & W rate of every `op_compare` entry, for each
   set and F1 kind.
3. For each set:
   - for ROI F1 and for pixel F1, a table with one row per entry: its thresholds
     (u/v/w), their F1 on the tune set, the pooled trio rates, and the number of
     valid trios;
   - the two comparison plots. For line sets, the per-sample panel is per angle.
   - For line sets only: efficiency and purity (recall and precision) on the
     sample's own plane (u, v or w0). ROI metrics are at each entry's best ROI
     F1 threshold and pixel metrics at its best pixel F1 threshold. There is a
     summary table: pooled over the angles (summed counts, e.g.
     `0.938 (181/193)`), plus the mean and minimum. There are tables of the
     counts per angle (found/true, real/reco), and plots of both against angle
     with the number of true ROIs or pixels under each angle
     (`<results>/op-{roi,pixel}F1-effpur-compare_Pplane_high_end.{png,npz}`).
     A line track at these angles crosses only about 12–60 channels of its own
     plane, with one ROI per channel. So per-angle ROI rates move in large steps,
     and the models often tie exactly.
4. For each model and each F1 kind, F-beta vs threshold on the tune set for u,
   v, w0 and uvw0.

For a quick check, evaluate cosmics on one directory and tune on another:

```
snakemake ... --config 'cosmics_tune_dirs=["2000_0"]' 'cosmics_eval_dirs=["2001_0"]' \
    -- results/op-trio-report.pdf
```

To leave the line samples out, add `'report_sets=["cosmics_eval"]'`.

On a large GPU, several inference jobs fit at once (xvunet takes about 7 GB):
`--resources gpu=4` lets four run together. Raise `--cores` too, because the
scan and trio jobs run on the CPU.
