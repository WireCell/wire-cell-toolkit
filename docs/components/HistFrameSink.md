# HistFrameSink

Saves frame traces as 2D ROOT TH2F histograms (one per wire plane) in a ROOT file, with time in microseconds on the X axis and channel number on the Y axis, applying a configurable unit scale.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `HistFrameSink` |
| Concrete class | `WireCell::Root::HistFrameSink` |
| Node category | sink |
| Primary interface | `IFrameSink` |
| Input type(s) | `IFrame` |
| Output type(s) | `(none)` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `filename` | printf-style pattern for the output ROOT file name, with %d replaced by frame ident (default "histframe-%02d.root") |
| `anode` | type:name of the IAnodePlane component (default "AnodePlane") |
| `units` | voltage/charge unit scale factor applied to waveform samples before filling histogram (default mV) |
