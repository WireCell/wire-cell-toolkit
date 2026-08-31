# MagnifySink

Saves frame waveform data as per-plane 2D ROOT TH2F histograms in the Magnify format, optionally copying existing objects from an input ROOT file, saving run-info metadata, trace summaries, and channel mask maps as TTrees.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `MagnifySink` |
| Concrete class | `WireCell::Root::MagnifySink` |
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
| `input_filename` | optional path of an input ROOT file for shunting objects to the output |
| `shunt` | array of TObject names to copy from input file to output file |
| `frames` | array of trace tags defining which waveforms are saved |
| `trace_has_tag` | if false, use untagged traces even when a tag name is given (default true) |
| `cmmtree` | array of [cmmkey, treename] pairs mapping channel mask keys to TTree names |
| `root_file_mode` | ROOT TFile open mode (default "RECREATE") |
| `runinfo` | JSON object with run-info key/value pairs to write into the Trun TTree, or null |
| `nrebin` | rebinning factor along time axis (default 1) |
| `summaries` | array of trace tags for which to save per-channel trace summary 1D histograms |
| `summary_operator` | object mapping tag to aggregation operator (default "set") |
