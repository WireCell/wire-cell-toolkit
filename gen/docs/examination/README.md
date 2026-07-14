# gen/ Code Examination — Summary Index

Date: 2026-04-12

This directory contains detailed examination documents for all ~72 source files
in the `gen/` directory (signal generation and simulation code). The examination
covers potential bugs, efficiency concerns, and algorithm documentation.

## Documents

| Document | Files Covered | Bugs | Efficiency Issues |
|----------|--------------|------|-------------------|
| [response-and-convolution.md](response-and-convolution.md) | 8 file pairs (PlaneImpactResponse, ImpactTransform, ImpactData, ColdElecResponse, WarmElecResponse, RCResponse, JsonElecResponse, ResponseSys) | 4 | 6 |
| [diffusion-and-deposition.md](diffusion-and-deposition.md) | 10 file pairs (GaussianDiffusion, BinnedDiffusion, BinnedDiffusion_transform, Diffuser, Diffusion, DepoTransform, DepoSplat, DepoFluxSplat, DepoPlaneX, DepoFramer) | 5+ | multiple |
| [noise-modeling.md](noise-modeling.md) | 8 file pairs (EmpiricalNoiseModel, GroupNoiseModel, AddNoise, AddGroupNoise, NoiseSource, Noise, SilentNoise, PerChannelVariation) | 6 medium, 6 minor | multiple |
| [drifting-and-transport.md](drifting-and-transport.md) | 6 files (Drifter, DepoSetDrifter, TrackDepos, BlipSource, RecombinationModels, TransportedDepo) | 5 | minor |
| [wire-and-anode-geometry.md](wire-and-anode-geometry.md) | 11 file pairs (AnodePlane, AnodeFace, MegaAnodePlane, WirePlane, WireParams, WireGenerator, WireSummary, WireSummarizer, WireSchemaFile, WireBoundedDepos, WireSource) | 8 | 3 |
| [digitization-and-pipeline.md](digitization-and-pipeline.md) | 29 file pairs (Digitizer, Reframer, FrameFanin/Fanout, FrameSummer, MultiDuctor, Fourdee, DuctorFramer, Scaler, TruthSmearer, TruthTraceID, Retagger, Misconfigure, DepoFanout, DepoMerger, DepoChunker, etc.) | 13 | 5 |

## Top-Priority Bugs (by likely impact)

### Correctness — wrong physics output

1. **Diffuser.cxx:124** — Transverse diffusion sigma uses `centimeter2` instead of `centimeter`, producing wrong diffusion widths.
2. **Drifter.cxx:26,29** — `DepoTimeCompare` compares `lhs->pos().x()` to itself (copy-paste error), breaking spatial tie-breaking in the depo set ordering.
3. **DepoFluxSplat.cxx:352** — `break` instead of `continue` skips all remaining planes when a depo is outside one plane's pitch range.
4. **Diffusion.cxx:16-37** — Copy constructor and `operator=` fail to copy `lbin`/`tbin`, making copied objects return wrong positions.
5. **ResponseSys.cxx:12** — Constructor writes `tick` (double) into `"nticks"` config key instead of the `nticks` int.
6. **TruthTraceID.cxx:153-218** — Wire iteration loop incorrectly nested inside depo loop, processing partial data.

### Crashes / undefined behavior

7. **Scaler.cxx:111** — No null/EOS check on depo pointer causes null dereference at end-of-stream.
8. **Scaler.cxx:152** — Unconditional dereference of `depo->prior()` which can be null for root depositions.
9. **Digitizer.cxx:136-151** — Invalid WPID channels leave nullptr in output trace vector, crashing downstream.
10. **RecombinationModels.cxx:54,97** — Division by zero when `dX=0` (the default parameter value).
11. **DepoSetRotate.cxx:56-66** — When `m_rotate` is false, all depos are silently dropped.

### Silent wrong behavior

12. **FrameSummer.cxx:54-57** — Adjusted `newtwo` frame constructed but never used; `align`/`offset` config options are non-functional.
13. **WireBoundedDepos.cxx:39 vs 67** — Config defines `"wires"` but code reads `"regions"`, causing silent no-op filtering.
14. **NoiseSource.cxx:154** — Config key `"m_nsamples"` vs `"nsamples"` means nsamples can never be configured.
15. **BinnedDiffusion_transform.cxx:320** — Hash `channel * 100000 + abs_tbin` collides for readouts > 100k ticks.

### Memory leaks

16. **WarmElecResponse / ResponseSys** — Raw `new` with no `delete` on reconfigure.
17. **WireSummary.cxx:189** — `m_cache` never deleted in destructor.
18. **WirePlane.h:33** — Raw `Pimpos*` never deleted.
19. **BlipSource.cxx** — Raw `new`/`delete` for distributions; calling `configure()` twice leaks.

## Top Efficiency Concerns

1. **ImpactTransform.cxx** — ~90 lines of duplicated code between paired-groups loop and central-group block; response spectra repeatedly inverse-FFT'd, truncated, and re-FFT'd per group instead of precomputed once.
2. **Scaler.cxx** — O(N^2 log N) sorting pattern.
3. **MultiDuctor.cxx** — Json::Value access in hot loops.
4. **TruthTraceID.cxx** — Redundant filter regeneration.
5. **MegaAnodePlane::channels()** — Copies/concatenates all sub-anode channels on every call with no caching or deduplication.
6. **AddGroupNoise.cxx:122** — Full complex vector copy + redundant inverse FFT per channel in same group instead of computing once per group.
