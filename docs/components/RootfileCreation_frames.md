# RootfileCreation_frames

Creates (or recreates) an empty ROOT output file at configuration time, then passes all input frames through unchanged; used to pre-create a ROOT file before downstream components write into it. Note: despite its name this component operates on IDepo, not IFrame.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `RootfileCreation_frames` |
| Concrete class | `WireCell::Root::RootfileCreation_frames` |
| Node category | function |
| Primary interface | `IDepoFilter` |
| Input type(s) | `IDepo` |
| Output type(s) | `IDepo` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `output_filename` | path of the ROOT file to create (required) |
| `root_file_mode` | ROOT TFile open mode (default "RECREATE") |
