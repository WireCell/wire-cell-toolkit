# ChargeErrorFrameEstimator

Estimates per-ROI charge error traces from an input tagged trace set by looking up pre-computed error waveforms indexed by ROI length, then appends those error traces to the output frame under a new tag.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `ChargeErrorFrameEstimator` |
| Concrete class | `WireCell::Img::ChargeErrorFrameEstimator` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | type/name of the IAnodePlane component (default: "AnodePlane") |
| `errors` | type/name of the IWaveformMap component providing per-plane ROI error waveforms (default: "WaveformMap") |
| `planes` | ordered list of plane names used to look up waveforms (default: ["u","v","w"]) |
| `tick` | waveform sampling period (default: 0.5 us) |
| `nbins` | number of bins in the ROI-length binning (default: 1001) |
| `intag` | trace tag selecting input signal traces (default: "gauss") |
| `outtag` | trace tag applied to the output error traces (default: "gauss_error") |
| `rebin` | rebinning factor that must match the error waveform's period/tick ratio (default: 4) |
| `fudge_factors` | per-plane multiplicative scale factors applied to the error estimates (default: [1,1,1]) |
| `time_limits` | [min, max] ROI length in ticks used to clamp the error lookup index (default: [12,800]) |
