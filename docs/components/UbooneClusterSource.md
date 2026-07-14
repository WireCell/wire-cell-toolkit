# UbooneClusterSource

Reads MicroBooNE-format WCP ROOT celltree files and, in combination with upstream IBlobSets, assembles a full point-cloud tree (clusters, blobs, slices, optical/flash data) and serializes it as an ITensorSet.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `UbooneClusterSource` |
| Concrete class | `WireCell::Root::UbooneClusterSource` |
| Node category | queuedout |
| Primary interface | `IBlobTensoring` |
| Input type(s) | `IBlobSet` |
| Output type(s) | `ITensorSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `input` | path or array of paths to input WCP ROOT files (required) |
| `kind` | blob kind, "live" or "dead" (default "live") |
| `sampler` | optional type:name of an IBlobSampler for populating point-cloud nodes |
| `datapath` | path template for the output tensor data (supports %d for frame ident) |
| `anode` | type:name of the IAnodePlane component (default "AnodePlane") |
| `time_offset` | time offset to apply (default 0) |
| `drift_speed` | drift speed value |
