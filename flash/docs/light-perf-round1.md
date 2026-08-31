# Light-chain performance round 1 (2026-07-08)

Profiling-driven CPU/RSS optimization of the optical reconstruction chain
(OpWaveformSource → OpDecon → [OpRoi] → OpHitFinder → OpHitMerge →
OpFlashFinder) for PDHD and PDVD.  Method: gperftools CPU profiles +
tcmalloc heap profiles + a /proc RSS sampler ported into the run drivers
(`PDVD_RESMON`/`PDHD_RESMON`, `pdvd/run_light_evt.sh`,
`pdhd/run_light_allpd_evt.sh`); byte-identity A/B via
`abtest/hash_archive.py` on the opflash tar.gz inner npy payloads (the
chain is re-run deterministic, verified before any change).

## Result summary

| event | wall before | wall after | peak RSS before | after |
|---|---|---|---|---|
| PDVD 39252/298567 | 2.31–2.39 s | 1.63–1.69 s | 0.63–0.64 GB | 0.35 GB |
| PDVD 39349/22877 (busiest) | 2.53 s | 1.78 s | 0.62 GB | 0.33 GB |
| PDHD 29107/983 | 3.08–3.23 s | 1.56–1.62 s | 0.43 GB | 0.39 GB |
| PDHD 29107/1015 (saturation evt) | 4.01 s | 1.88 s | 0.43 GB | 0.39 GB |

Wall is the wire-cell process (RESMON), not the driver (python decoana /
jsonnet compile excluded).

## What the baseline profiles showed

- PDVD spent **39% of CPU and 70% of peak heap in `OpDecon::configure`**:
  each of the 3 branch instances FFT'd and stored ALL 51 SPE templates at
  its own `samples` (cathode: 51 × 468864-sample spectra ≈ 190 MB + the
  padded waves, for 16 used channels).
- All flash FFTs ran as **full complex transforms on real data** through
  `DftTools::fwd_r2c`/`inv_c2r` (2× flops, fresh buffers per call; ~2.3 GB
  allocation churn per PDVD event).  PDHD's default path did 5 full-length
  transforms per record (fwd, inv, auto-scale inv, post-filter fwd+inv).
- Lesser: OpRoi made 4 full-length copies per trace (dead input copy + 3
  by-value medians); OpFlashFinder used `vector::erase` per dropped flash
  and allocated a hits-sized `used` vector per candidate flash; the PDVD
  source deserialized the unused x/y/z branches 3× per event.

## Changes

**Byte-identical (no knobs; all baseline archive hashes reproduced):**

1. `OpDecon`: lazy per-channel template spectra — `configure()` keeps only
   the unpadded wave; the padded FFT (and `wi_eps`) is built on a
   template's first use (commit `935a515b`).
2. `OpDecon`: with `fixed_snr > 0` and flat noise the Wiener filter is
   record-independent, so the AutoScale normalization (one extra inverse
   FFT per record) is cached per template (`b636c59b`).
3. `OpRoi`: pass the decon straight to the DFT (dead copy), reuse median
   scratch buffers, in-place MAD median; `Aux::SimpleTrace` gained a move
   ctor used by OpDecon/OpRoi/PDVDOpWaveformSource (`b636c59b`, `9054bb65`).
4. `OpFlashFinder`: mark+compact instead of per-flash `vector::erase` in
   `remove_late_light` and the quality cut; `refine_hits_in_flash` reuses
   a caller scratch `used` vector (`9054bb65`).
5. `PDVDOpWaveformSource`: pass 1 enables only the 6 selection/timing
   branches (`9054bb65`).

**Knob-gated, default OFF → bit-identical; enabled in the pdhd +
protodunevd light cfgs after validation (`b4894c8e` code, cfg commit
after):**

6. **Real FFTs** — `IDFT::fwd_r2c_1d`/`inv_c2r_1d` (default impls widen to
   the complex path; `FftwDFT` overrides with cached native r2c/c2r
   plans), `DftTools::fwd_r2c_real`/`inv_c2r_real`, and OpDecon/OpRoi knob
   `use_real_dft`.
7. **Folded post-filter** — OpDecon knob `fold_postfilter` multiplies the
   Gauss post-filter into the decon spectrum (saves fwd+inv per record).
   The post baseline correction keeps its exact pre-filter semantics: the
   head-pedestal mean is a linear functional of the spectrum, evaluated
   with precomputed geometric weights (`m_ped_w`) at O(N) — measuring the
   pedestal on the *filtered* wave instead was tried first and moved
   marginal flashes (221→206 on PDHD 983), so it was rejected.

## Knob-ON validation

- PDVD decon vs the independent float64 python mirror
  (`pdvd/pd_plot/wct_wi_validate.py`): max rel diff 9.9e-7 (PMT), 3.4e-7
  (membrane), 2.6e-7 (cathode) — the same float32-template level as the
  complex path.
- Flash level vs the pre-round baselines: PDVD 39252/298567 identical 398
  flashes, |Δt| ≤ 0.08 ns, rel ΔPE ≤ 8e-4.  PDHD 29107/983 identical 221
  flashes, |Δt| ≤ 3e-3 ns, rel ΔPE ≤ 4e-5.  PDHD 29107/1015 (railed
  channels) identical 768 flashes; 4 dim flashes sit at hit-split
  boundaries and moved (worst 92 ns / 4% PE), the other 764 at ≤0.9 ns /
  ≤3.4e-4 — round-off amplified by a discrete threshold on clipped
  flat-tops, accepted and documented.
- Knobs OFF reproduce every baseline hash bit-identically (re-verified
  with the final code).

## What remains (deliberately not done)

- FFT work still dominates both chains (~50% + OpRoi ~25%); further gains
  would need half-spectrum storage/processing (r2c layout end-to-end) or
  smooth-length padding — measured FFTW at 468864/343808 was not the
  bottleneck vs. the transform count, so `samples` was left alone.
- The PDVD 3-instance source still scans the tree metadata 3×/event
  (~2% CPU after the branch trim — not worth restructuring).
- PDHD `read_snippets` cost is ROOT `TKey::ReadObj` deserialization, not
  the bin copy loop.
- `OpHitFinder` is ~2% — reserve/cast-buffer micro-opts skipped.

## Reproduce

```bash
# profiles (write to scratch, never clobber work/)
pdvd/profile_light.sh 39252 298567 [out.prof]
pdhd/profile_light.sh 29107 983    [out.prof]   # needs prior allpd run
HEAPOUT=/path/heap pdvd/profile_light.sh ...    # tcmalloc heap mode
google-pprof --text --cum $(which wire-cell) out.prof | head -40
# RSS trace per run: work/<dir>/light_rss_*.csv + light_resource_*.txt
# A/B: python3 abtest/hash_archive.py <opflash tar.gz> ...
```
