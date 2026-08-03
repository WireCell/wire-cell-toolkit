# NF/SP/DNNROI performance profile — round 0 (PDHD)

**Status: measurement only, no code change.** This is the baseline profiling
round for the full standalone NF+SP+DNNROI+L1SP chain
(`pdhd/run_nf_sp_dnnroi_evt.sh` production settings: fp32 6-ch DNN model,
L1SP on in `dnn` mode, all 4 APAs sequential in one `wire-cell` process),
in the spirit of `clus/docs/imgclus-optimization-log.md` and
`flash/docs/light-perf-round1.md`.

Motivation: per `pdhd/docs/pdhd-pipeline-resource-profile.md`, NF/SP/DNNROI is
the **largest wall-clock stage of the whole pipeline for a typical event**
(~92 s median / ~1.9–2.4 GB peak, flat across events) yet — unlike imaging,
clustering, QLMatching and the light chain — it never had a dedicated
optimization round. The only prior profile (`docs/l1sp/L1SPFilterPD.md`) was
L1SP-focused and single-APA.

## 0. Method

- Driver: `pdhd/profile_nfsp.sh <run> <evt>` (new; mirrors `profile_light.sh`).
  Pre-compiles the cfg with `wcsonnet` (SIGPROF kills gojsonnet at config
  compile otherwise), writes outputs to scratch, and samples
  `VmRSS`/`VmHWM` at 2 Hz alongside.
- CPU: gperftools `LD_PRELOAD=libtcmalloc_and_profiler.so.4`,
  `CPUPROFILE_FREQUENCY=1000`, analyzed with `google-pprof --text [--cum|--focus]`.
- Heap: jemalloc sampling (`MALLOC_CONF=prof:true,lg_prof_sample:19,`
  `lg_prof_interval:33,prof_final:true`, env `JEHEAP=<dir>`), peak = interval
  dump with max in-use, analyzed with jeprof.
- Events: run 029107 charge idx 0 (DAQ evt 983, typical) and idx 4
  (DAQ evt 1015, the bright outlier). NF/SP wall is known population-flat
  (89–103 s over 30 events), so two events suffice for attribution.
- Baseline (un-profiled, from existing `work/029107_*/` logs + `time_*.txt`):
  wall 82–102 s single-job; `VmHWM` 2.11–2.63 GB over 30 events (median ~2.27).

## 1. CPU — where the ~90 s goes (evt idx 0, 19701 samples)

| Component | share | notes |
|---|---|---|
| `OmnibusSigProc` (SP, 4 APAs) | **41.7%** | FFT-dominated; `decon_2D_init` = 35.6% of SP, then tightROI/hits/charge/tighterROI/looseROI ~5–8% each; `restore_baseline`+`Waveform::percentile` ~9% of SP each |
| `FrameFileSink` output | **20.9%** | **19.3% = bzip2 compress** (`BZ2_bzCompress`) of the output tar.bz2 |
| OpenMP worker-thread spin | **11.6%** | `omp_get_num_procs` self time on `clone3/start_thread` threads = libtorch intraop pool spin-waiting between the (tiny) DNN calls. Wasted core-seconds, not wall — but real on shared batch nodes |
| `OmnibusNoiseFilter` (NF) | **8.7%** | 63.5% of it `PDHD::OneChannelNoise` (per-channel FFT filters), 35.6% `CoherentNoiseSub` (`Waveform::percentile_binned` = 32.5% of NF) |
| Input bzip2 decompress | ~4.2% | `FrameFileSource` reading orig tar.bz2 |
| config/init | ~4% | field response + jsonnet-compiled JSON load |
| `L1SPFilterPD` | 1.9% | matches the L1SP doc's earlier finding (1.6%) |
| DNNROI torch inference | ~1.1% | the DNN itself is CPU-negligible |

**FFT total at the `FftwDFT` level: ~29.3%** of the whole job
(`inv1b` 14.3% + `fwd1b` 9.1% + `inv1d` 3.1% + `fwd1d` 2.8%) — and all of it
runs through the *widening* `DftTools::fwd_r2c`/`inv_c2r` paths (complex c2c
transforms on real data). The native r2c/c2r plan infrastructure added for the
light chain (`IDFT::fwd_r2c_1d/inv_c2r_1d`, commits 935a515b/b4894c8e) covers
only the **1d** entry points; SP's hot paths are the **batched 1b** variants.

## 2. Memory — where the ~2.3 GB peak comes from (evt idx 0)

RSS-vs-time (2 Hz sampler): ramps to **~1.3 GB** through NF+SP of the first
APA, then a **single +1.05 GB step at the first DNN-ROI inference**, then
**dead flat ~2.36 GB** to the end (allocator retains the transient; later DNN
calls reuse it). Peak = first-DNN-call transient sitting on the SP working set.

jemalloc peak dump (1.65 GB in-use, scaled):

| Owner | MB | breakdown |
|---|---|---|
| `DNNROIFinding` | **811 (49%)** | 736 MB `c10::alloc_cpu`, of which **645 MB inside one `at::native::_convolution` (mkldnn transposed-conv decoder)** — activation/workspace transient of a full-plane (800×1500×6ch) forward pass |
| `OmnibusSigProc` | **673 (41%)** | 660 MB `SimpleTrace` buffers: `save_data` 300 + `save_roi` 185 + `save_mproi` 124 + `save_ext_roi` 54 — the many tagged trace sets (gauss/wiener/tight_lf/mp2/mp3/…) held concurrently to feed DNNROI+L1SP |
| everything else | ~165 | Eigen scratch 31, bzip2 codecs ~8, torch model weights, config JSON |

## 3. Ranked levers (round 1 candidates)

1. **Batched real-FFT path in SP** (~29% job CPU ceiling, realistic −10–15%):
   extend the existing `IDFT` r2c interface with batched `fwd_r2c_1b/inv_c2r_1b`
   (FftwDFT: cached native `fftwf_plan_many_dft_r2c/c2r` plans) and route
   `OmnibusSigProc`/NF `DftTools` calls through it behind a `use_real_dft`-style
   knob (default OFF = bit-identical, per the light-chain precedent).
2. **Output codec** (19.3% job CPU): `FrameFileSink` bzip2 → zstd (or lower
   bzip2 block size / level knob). Needs coordination: everything downstream
   reads `…-frames-anodeN.tar.bz2` (imaging `FrameFileSource`, python
   `magnify`-style tooling). A `.tar.zst` sibling path keeps A/B trivial.
3. **Torch thread-pool spin** (11.6% core-seconds, ~0 wall single-job):
   `OMP_NUM_THREADS`/`OMP_WAIT_POLICY=passive` (or
   `at::set_num_threads(1)` in `TorchService`) — matters for batch throughput
   (`-P 6` style parallel sweeps share cores). CAUTION: thread-count changes
   mkldnn reduction order → must A/B the DNN output for drift before adopting.
4. **DNN peak transient** (−0.5–1.0 GB RSS ceiling): tile the DNNROI forward
   pass (chunk the 800×1500 image) or evaluate the existing `-P int8` preset's
   footprint. Result-changing at tile seams → knob + physics validation.
   Alternatively cheap and exact: `c10::emptyCache()`-style release after each
   call does nothing for VmHWM (peak already paid) — only tiling shrinks it.
5. **NF `percentile_binned`** (2.8% job): coherent-group median via
   `nth_element` on the raw vector instead of the binned histogram — small,
   bit-exactness needs care; low priority.
6. NOT worth it at this scale: L1SP internals (1.9%), DNNROI inference math
   (1.1%), input decompress (4%, would ride along with lever 2 if inputs are
   ever re-encoded).

## 4. Busy-event check (evt idx 4 / DAQ 1015)

Same shape within noise — the chain is activity-independent, so the typical-event
attribution carries over to the whole population:

| Component | idx 0 (typical) | idx 4 (bright outlier) |
|---|---|---|
| OmnibusSigProc | 41.7% | 43.6% |
| FrameFileSink (bzip2) | 20.9% (19.3%) | 19.3% (17.7%) |
| OMP spin | 11.6% | 12.3% |
| OmnibusNoiseFilter | 8.7% | 8.5% |
| L1SP / DNNROI | 1.9% / 1.1% | 1.9% / 1.1% |
| heap peak (jemalloc) | 1649 MB: DNNROI 811 / SP 673 | 1653 MB: DNNROI 799 / SP 679 |
| VmHWM (profiled run) | 2.37 GB | 2.44 GB |

## 5. Artifacts

- Profiles: `/home/xqian/tmp/nfsp_pdhd_029107_{0,4}.prof` (CPU),
  `/home/xqian/tmp/nfsp_heap_evt{0,4}/jeprof.*.heap` (jemalloc; peak = i17/i8).
- RSS traces: `/home/xqian/tmp/prof_nfsp_pdhd_029107_*/rss_*.csv`.
- Driver: `pdhd/profile_nfsp.sh` (wcp-porting-img repo).
