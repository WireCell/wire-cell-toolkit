# Deploying master of 2026-10-01: what changes in production and how to validate it

On 2026-10-01 `master` was fast-forwarded from `51b5a1fc` to the `apply-pointcloud` branch head,
the commit that adds this file. This note is for whoever deploys the new master into a production
chain. It lists what changes, what has already been checked, and what still has to be checked in
your environment before the new build replaces the old one.

```bash
git log --oneline --no-merges 51b5a1fc..<new master>   # the 22 commits this note covers, plus this file
```

## 1. What changes in production output

Relative to `51b5a1fc`, exactly **three** behaviour changes reach a production job. All three are SBND only.

| # | Where | Change | Commit | Output effect |
|---|---|---|---|---|
| 1 | SBND step 1 (LArSoft, `flash_source='hits'`) and the SBND standalone hit-flash chain | `SBNDOpFlashFinder` picks the prompt-time bin only among bins with ≥ `prompt_min_pe` (6) PE and ≥ `prompt_min_hits` (3) hits (SimpleFlashAlgo `MinPECoinc`/`MinMultCoinc`), falling back to the brightest bin. **The default is on, there is no config key in the SBND jobs, and the compiled config does not change.** | `4afd5cac` | Flash times move on some flashes; Q/L matching downstream may follow. |
| 2 | SBND PR: step 2 `wct-pr.jsonnet`, `wct-pr-perevt.jsonnet`, the obsolete 1-step | `main_vertex_swap_apply: true` on `TaggerCheckNeutrino` (apply the traditional main-vertex path's cluster swap) | `2bb3b96b` | Main cluster / vertex changes on the events where a swap occurs. |
| 3 | Same SBND PR jobs | `nu_particle_links: true` on `UbooneTaggerOutputVisitor` (`root_particle_links` in `pr()`) | `19e6d39a` (C++ in `2c528dc4`) | `tracking-pr.root` gains the `T_segment` tree, per-row identity columns in `T_kine` (`kine_particle_id`, …) and `act_n_seg*` in `T_tagger`. Additive: existing branches are unchanged. |

Everything else in the range is default-off, behaviour-neutral, or for other detectors:
- **Default-off knobs:**
  - `CascadeDeghosting` and `TorchTensorSetService` (img, pytorch);
  - the `ClusteringIsolated` radius knobs;
  - the `OpFlashFinder` late-light knobs, whose defaults equal the former literals;
  - the ICARUS `QLMatching`, `TrackFitting` and `dead_region` knobs;
  - `main_vertex_swap_discard_clean`;
  - the log-only `kine_overlap_probe`.
- **Behaviour-neutral fixes:**
  - `ClusterArrays::to_arrays` now edits the cluster graph in place;
  - `sep_family` is zero-filled before serialisation (a crash fix);
  - sigproc test fixes.

Already in `51b5a1fc` and therefore **not new** with this deploy:
- PR #535: the 2-step chain, `reset_shower_ids_per_event`, and step 1's `wclsTruthInformationAttacher`, which needs a larwirecell that provides it;
- `r_erase_stale_contained` ON.

## 2. What has already been checked (2026-10-01)

| Check | Result |
|---|---|
| Compiled configs, old master vs new master: the 28 production consumers (PDHD, PDVD, uBooNE and SBND jobs) plus the two 2-step jobs | Only the SBND PR jobs move: step 2, `wct-pr-perevt`, and the obsolete 1-step job. Each moves by exactly the two keys of changes 2–3. **Step 1 is byte-identical.** Every PDHD, PDVD and uBooNE artifact is byte-identical. |
| Unit tests | All 16 per-package `wcdoctest-<pkg>` (waf) pass; `clus` is 471 cases with 0 failed. |
| Strict release build (`WCT_BUILD_MODE=release`, `-Werror -Wpedantic` confirmed in `flags.make`, ROOT on) | 0 compile errors. All 18 libraries build. The only failures are link failures in 8 pytorch test programs and the aggregate `wcdoctest`: unresolved symbols in the system `libblas`/`liblapack`/`libhdf5`/`libcurl` pulled in by the system libtorch 1.13. That is an environment problem, the same as in the 2026-09-25 build. |
| Output A/B of the merge (branch before vs after taking `51b5a1fc`), at the content of archive members and ROOT trees | SBND PR: 67 events, identical. SBND img→clus→Q/L: 19 events, 202/202 archives identical. uBooNE: 35 events, Bee zips identical, tagger output 34/35 (see §4). PDHD/PDVD clustering: 120/120 archives identical. PDHD/PDVD PR: 4 jobs, identical. |

Not checked here:
- LArSoft (`lar`) jobs;
- the effect of change 1 on Q/L matching rates;
- the step-2 multi-event process against old master.

## 3. Validation checklist for the deployment

Run every comparison by **content**, never by file bytes. Tar and zip archives embed timestamps, so
compare the member names plus payload hashes. Compare ROOT files tree by tree and branch by branch.

**1. Build and test** in your environment: cmake with `WCT_WITH_TESTS=ON`, then `build/wcdoctest`.
Every failure should be one you can attribute to the environment, such as the libtorch link above.

**2. Config proof.** Compile your production entry points with old and new master and diff the JSON,
node by node. Expect exactly the diffs in the §2 table. Step 2, for example:
```bash
wcsonnet -P <cfg> -P <wire-cell-data> -A input=qlpctree.tar.gz -A reality=sim \
    <cfg>/pgrapher/experiment/sbnd/wct-pr.jsonnet > step2.json
```
Step 1 needs the external variables your `.fcl` supplies (`reality`, `recobwire_tags`, `trace_tags`,
`summary_tags`, `input_mask_tags`, `output_mask_tags`, `opflash{0,1}_input_label`,
`ophit{0,1}_input_label`, `enable_tracking_root`, …), passed with `-V`.

**3. SBND step 1, old vs new**, on the same reco1 events (≥ 50, MC and data):
- Expect differences only in the flash content (times) of `qlpctree.tar.gz` and `mabc.zip`, and in whatever Q/L matching makes of them.
- Report the fraction of flashes whose time moved, and the change in matched bundles.
- **Attribution check:** run new master with the old prompt-time rule. Every bin then qualifies, so it reduces to the previous "brightest bin, earliest on ties" rule. This is by reading the code; this run confirms it. The output should then equal old master:
  ```jsonnet
  // step1_oldflash.jsonnet
  (import 'pgrapher/experiment/sbnd/wcls-img-clus-matching-pr-lib.jsonnet')(
      flash_source='hits', xtpc_sc1_light_gate=true, xtpc_sc1_overpred_max=2.9, stage='ql',
      ff={prompt_min_pe: 0, prompt_min_hits: 0, prompt_min_opdets: 0})
  ```

**4. SBND step 2, old vs new**, on the **same** `qlpctree.tar.gz`:
- Expect `T_segment` and the new `T_kine`/`T_tagger` branches to be added.
- Every branch both versions have should be identical, except on events where the main-vertex swap applies (change 2).
- **Attribution check:** put a copy of `pgrapher/experiment/sbnd/clus.jsonnet` first on the config path (`-P overlay -P <cfg>`) with `main_vertex_swap_apply: false` and `root_particle_links=false`.
  - The compiled step 2 then equals old master's, except for an explicit `main_vertex_swap_apply: false` where old master had no key. That is the same C++ default.
  - `tracking-pr.root` should then equal old master's.
  - For the per-event job the switches are TLAs: `--tla-code main_vertex_swap_apply=false -A root_particle_links=false`.

**5. PDHD, PDVD, uBooNE.** Their compiled configs do not change. A spot check of a few events per
detector should give identical archive and ROOT content.

## 4. Known run-to-run differences: not caused by this deploy

Before treating any difference as real, run the **old** build twice on the same input.

| Detector | Known difference |
|---|---|
| uBooNE run 5384 event 6805 | `kine_pio_*_2` and `kine_pio_angle` alternate between two states, 14.81° and 109.51°. The pre-merge binary alone gave 14.81 twice and 109.51 four times in six runs. |
| DL/SCN vertex | Not bit-stable; keep it out of byte comparisons. |
| SP, run to run | Use the same binary, threads and input; compare runs pairwise before comparing builds. |
