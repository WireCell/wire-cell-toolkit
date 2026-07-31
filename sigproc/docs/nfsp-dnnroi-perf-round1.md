# NF/SP/DNNROI performance — round 1 (PDHD): r2c knob + zstd sink codec

Round 0 (`nfsp-dnnroi-perf-round0.md`) attributed the ~78 s / ~2.3 GB
NF+SP+DNNROI+L1SP job and ranked levers.  This round implements and
gates the top levers.  All changes are **knob-gated, default OFF; the
OFF path is byte-identical** to the pre-change baseline (verified, see
§4).

## 1. Lever 1 — native batched real FFTs (`use_real_dft`)

New interface + implementation:

- `IDFT::fwd_r2c_1b / inv_c2r_1b` — batched real<->complex transforms
  following the `fwd1b/inv1b` array/axis conventions (full-size
  Hermitian spectrum out / Hermitian-assumed spectrum in, `1/n`
  normalization).  Default implementations widen to the complex 1b
  transforms; `FftwDFT` overrides with cached native
  `fftwf_plan_many_dft_r2c/c2r` plans (`FFTW_ESTIMATE|PRESERVE_INPUT|
  UNALIGNED`, deterministic).
- `DftTools::fwd_r2c_real / inv_c2r_real` array-with-axis overloads.
- `OmnibusSigProc` config knob `use_real_dft` (default false) routes
  every real<->complex time-axis call (`decon_2D_init`, the 7 ROI
  filter stages, response/electronics spectra) through the native
  path.  Wire-axis complex<->complex transforms are untouched.
- jsonnet: `pdhd/wct-nf-sp-dnnroi.jsonnet` TLA `sp_use_real_dft`
  (default false; key omitted when off so the compiled config is
  byte-identical).

**Numerics**: native vs widened agree to ~2.5e-7 relative (float
round-off), both axes, odd/even lengths and the production shapes
(standalone check + job A/B).

### What the round-0 estimate got wrong

Round 0 put "FFT = 29.3% of job CPU" from gperftools samples and
projected −10–15% wall.  Two corrections from this round:

1. **The sampling was biased.** gperftools' ITIMER_PROF signals
   coalesce badly in this job (~19 k samples over ~270 core-seconds)
   because the idle libtorch OMP pool inflates core-seconds; sample
   *shares* over-weight the single-threaded stretches.  By the WCT
   node timers the honest budget is: SP ≈ 34 s node-CPU of which
   time-axis FFT ≈ 12–14 s, sink ≈ 17 s, NF ≈ 7 s (per event, 4 APAs).
2. **The strided time-axis geometry is memory-bound, not flop-bound.**
   At the per-plane shapes (~1029–1200 wires × 6174 ticks, transform
   along the strided axis) FFTW's strided c2c and r2c cost nearly the
   same (fwd −6%…0%, inv −16…19% measured through DftTools).  Halving
   flops does not halve time.  A transpose-to-contiguous strategy that
   wins 4–8x on *large* batch counts (2560+) was tried and **rejected**:
   at production shapes the two extra full-array passes plus ~75 MB
   scratch churn per call made SP *slower* (36.5 vs 34.2 s node CPU).
   Keep the direct strided plans.

**Measured effect (run 029107 evt idx 0):** SP node CPU 34.7 → 33.2 s;
job wall ~ −2 s (~−3%), within run-to-run noise of ±2 s.  Kept because
it is free when off, helps every OmnibusSigProc consumer, and the same
interface serves future half-spectrum work — but it is NOT the
headline lever round 0 hoped for.

**Output effect when ON:** `raw` frames identical; `gauss`/`wiener`
frames differ on ~1.2% of samples (max ~0.7% of frame scale), with
gauss and wiener flipping identically → ROI-boundary shifts from
round-off, same character as the light-chain `use_real_dft`
validation.  Production adoption requires the downstream physics A/B
(imaging + clustering population) — until then the knob stays OFF.

## 2. Lever 2 — output archive codec (`sp_sink_ext`)

bzip2 compression in `FrameFileSink` is ~17 s node-CPU per event
(≈20% of wall; the sink runs on the main thread).  The codec is chosen
by *filename extension* (custard), so the knob is the outname:

- `util/.../custard/custard_boost.hpp`: added `zst|tzst` branches
  (boost::iostreams zstd, available since boost 1.71; this build's
  boost 1.85 already links libzstd).  Both `input_filters` and
  `output_filters`, so `FrameFileSource` reads `.tar.zst` transparently.
- `pdhd/wct-nf-sp-dnnroi.jsonnet` TLA `sp_sink_ext` (default
  `tar.bz2` = byte-identical config).

Codec measurements on a real SP frame archive (184 MB tar payload,
one anode):

| codec | size | compress time |
|---|---|---|
| bzip2 -1 (production) | 61.6 MB | 3.87 s |
| gzip -1 | 61.3 MB | 1.85 s |
| **zstd -1..-3** | **60.2 MB** | **0.14 s (27x)** |

Payload bytes are identical — only the container changes.  Rollout
coordination when flipping production: downstream drivers glob
`*-frames-anode*.tar.bz2` explicitly (`run_img_evt.sh`,
`run_clus_evt.sh`, `run_sp_to_magnify_evt.sh`) and python tooling
needs the `zstandard` module (python's `tarfile` has no zst support).

**Measured effect (full-ON job = r2c + zstd):** see §4.

## 3. Checked and dropped

- **NF `PDHDOneChannelNoise` r2c knob — tried and REVERTED.**  The
  per-channel filter FFTs are ~4 s node-CPU (≤2 s gain), but the FFT
  round-off is *amplified* by the post-filter baseline estimator:
  `Waveform::median_binned` quantizes to bins of `(vmax-vmin)/nticks`
  anchored at the waveform extremes, so an O(1e-4)-ADC round-off can
  move the whole-channel baseline by up to a bin (tens of ADC on a
  channel with a bright pulse).  Measured with the knob ON: `raw`
  frames shift on ~100% of samples (max 22 ADC / 0.8% of scale) and
  gauss/wiener grow 28%-of-scale outliers through the SP ROI
  thresholds.  Not a round-off-class knob; the gain does not justify
  the physics-validation burden.  (Contrast SP's knob: `raw` frames
  identical, gauss/wiener bounded at 0.7%.)
- **EFF-CORE-4 `c_data_afterfilter` workspace reuse**: the fill loops
  and allocations are ≤0.1% self-time each in the round-1 profile —
  the round-0 "8 allocs/plane" concern is not measurable.  Dropped.
- **Transpose-to-contiguous strided FFT strategy**: see §1; rejected
  at production shapes.
- **FFTW_MEASURE plans**: 1.8x better than ESTIMATE on strided r2c but
  still worse than the levers above, and planner output is
  run-to-run non-deterministic (violates the determinism policy).

## 4. Gates

- **Knob-OFF byte-identity**: with all knobs at defaults the compiled
  config is byte-identical (`cmp` on wcsonnet output) and the four
  output archives hash identical (`hash_archive.py`) to the
  pre-change baseline — which itself reproduces the production
  `work/029107_0/` archives exactly.
- **Numerics**: standalone `fwd_r2c_1b/inv_c2r_1b` vs widened path
  ≤2.6e-7 relative on production shapes, both axes.
- **Wall/RSS + full-ON A/B**: (numbers in §5).

## 5. Results (run 029107 evt idx 0, single job, quiet node)

| config | wall | sink node-CPU | SP node-CPU | notes |
|---|---|---|---|---|
| baseline (pre-change binary) | 77.5 s | 17.0 s | 34.7 s | == production `work/029107_0/` archives |
| knob-OFF (new binary) | 77.3 s | 17.0 s | 34.2 s | outputs byte-identical to baseline (x3 gates) |
| SP r2c ON | 75.3 s | — | 33.2 s | `raw` identical; gauss/wiener 1.2% samples, ≤0.7% scale |
| **SP r2c + zstd sink ON** | **60.6 s (−22%)** | **1.9 s** | 33.1 s | archive 61.6→60.1 MB/anode; payload verified vs bz2 |

Wall noise is ±2 s run-to-run; the −17 s decomposes as zstd −15 s +
r2c −2 s.  RSS is unchanged (the DNNROI transient still sets the
peak).  The `.tar.zst` payload was decompressed and compared
member-wise against the bz2 output: identical modulo the r2c
round-off pattern above.

## 6. Artifacts

- A/B work dirs: `/home/xqian/tmp/nfsp_r2c_ab/{base,off,on,fullon}_029107_0`
  (cfg.json + outputs + time.txt), `cmp_frames.py` member-wise
  comparator.
- Standalone checks: session scratchpad `check_r2c_1b.cxx` (numerics +
  production-shape timing), `bench_r2c_1b.cxx` (strided/contiguous
  FFTW variants).
