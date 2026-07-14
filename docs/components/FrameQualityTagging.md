# FrameQualityTagging

Assesses frame quality by detecting noisy or overly busy time regions using Wiener and Gaussian filtered traces; optionally updates the channel mask map and returns the frame or an empty frame if quality is too poor.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `FrameQualityTagging` |
| Concrete class | `WireCell::Img::FrameQualityTagging` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | type/name of the IAnodePlane component (default: "AnodePlane") |
| `gauss_trace_tag` | tag for Gaussian-filtered traces |
| `wiener_trace_tag` | tag for Wiener-filtered traces |
| `tick0` | first tick index to consider |
| `nticks` | total number of ticks to analyze |
| `nrebin` | number of ticks to rebin when scanning for noise |
| `length_cut` | minimum channel-span for an ROI to count |
| `time_cut` | minimum gap between separate noisy regions |
| `ch_threshold` | minimum charge-weighted channel count per rebinned time bin |
| `n_cover_cut1` | first cover-count threshold for declaring a noisy region |
| `n_fire_cut1` | first fire-count threshold for declaring a noisy region |
| `n_cover_cut2` | second (alternative) cover-count threshold |
| `n_fire_cut2` | second (alternative) fire-count threshold |
| `flag_corr` | if 1, propagate detected noisy time ranges back into the channel mask map |
| `global_threshold` | sum-of-plane fired-fraction threshold above which the frame is declared too busy |
