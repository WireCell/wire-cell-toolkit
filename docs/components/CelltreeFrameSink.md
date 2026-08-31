# CelltreeFrameSink

Saves frame waveform data to a ROOT celltree file (Event/Sim TTree) as used by the MicroBooNE event display tools, supporting raw/gauss/wiener tagged traces, trace summaries, and channel mask maps.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `CelltreeFrameSink` |
| Concrete class | `WireCell::Root::CelltreeFrameSink` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | type:name of the IAnodePlane component (default "AnodePlane") |
| `output_filename` | path of the output ROOT file (required) |
| `frames` | array of trace tags to save as waveform branches |
| `cmmtree` | array of channel mask map keys to save as bad-channel branches |
| `root_file_mode` | ROOT TFile open mode, e.g. "RECREATE" or "UPDATE" (default "RECREATE") |
| `nsamples` | number of time samples per waveform (required) |
| `nrebin` | rebinning factor along time axis (default 1) |
| `summaries` | array of trace tags for which to save per-channel trace summary vectors |
| `summary_operator` | object mapping tag to aggregation operator, "sum" or "set" (default "sum") |
