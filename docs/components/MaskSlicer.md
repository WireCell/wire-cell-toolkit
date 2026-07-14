# MaskSlicer

Slices a frame into time slices using Wiener-filter-based activity thresholding per channel, supporting active, dummy, and masked plane configurations; outputs all slices from one frame as a single ISliceFrame.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `MaskSlicer` |
| Concrete class | `WireCell::Img::MaskSlicer` |
| Node category | queuedout |
| Primary interface | `IFrameSlicer` |
| Input type(s) | `IFrame` |
| Output type(s) | `ISlice` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | type/name of the IAnodePlane component |
| `tick_span` | number of ticks summed into each output slice |
| `wiener_tag` | trace tag for Wiener-filtered waveforms used for activity thresholding |
| `summary_tag` | trace tag whose summary (RMS) provides per-channel thresholds |
| `charge_tag` | trace tag for charge waveforms assigned to slice activity |
| `error_tag` | trace tag for charge error waveforms |
| `active_planes` | array of plane indices (0=U,1=V,2=W) that contribute real activity |
| `dummy_planes` | array of plane indices filled with dummy_charge/dummy_error |
| `masked_planes` | array of plane indices where "bad" channel mask ranges are filled with masked_charge/masked_error |
| `dummy_charge` | charge value assigned to dummy-plane channels |
| `dummy_error` | error value assigned to dummy-plane channels |
| `masked_charge` | charge value assigned to masked channel time ranges |
| `masked_error` | error value assigned to masked channel time ranges |
| `nthreshold` | per-plane multiplier applied to per-channel RMS to form the activity threshold |
| `default_threshold` | per-plane fallback threshold when nthreshold*RMS evaluates to zero |
| `min_tbin` | first tick bin to process (0 means determine from input frame) |
| `max_tbin` | last tick bin to process (0 means determine from input frame) |
