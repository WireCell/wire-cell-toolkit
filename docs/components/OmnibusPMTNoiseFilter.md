# OmnibusPMTNoiseFilter

Identifies and removes PMT-flash-induced noise from collection and induction plane waveforms by detecting characteristic negative pulses in collection traces and correlating them with induction-plane responses.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `OmnibusPMTNoiseFilter` |
| Concrete class | `WireCell::SigProc::OmnibusPMTNoiseFilter` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | type:name of the IAnodePlane component (required) |
| `pad_window` | number of bins to search on each side of a PMT ROI for the adaptive baseline |
| `min_window_length` | minimum bin length of a PMT ROI to be processed |
| `threshold` | threshold in units of per-channel RMS for PMT signal detection |
| `intraces` | input trace tag (default "quiet") |
| `outtraces` | output trace tag (default "raw") |
