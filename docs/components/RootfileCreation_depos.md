# RootfileCreation_depos

Creates (or recreates) an empty ROOT output file at configuration time, then passes all input depositions through unchanged; used to pre-create a ROOT file before downstream components write into it.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `RootfileCreation_depos` |
| Concrete class | `WireCell::Root::RootfileCreation_depos` |
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
