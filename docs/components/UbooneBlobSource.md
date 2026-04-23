# UbooneBlobSource

Reads MicroBooNE-format WCP ROOT celltree files and produces IBlobSet objects for live or dead blobs, reconstructing ray-grid blobs, slice activity maps, and optional 2-view channel bodging from bad-channel masks.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `UbooneBlobSource` |
| Concrete class | `WireCell::Root::UbooneBlobSource` |
| Node category | source |
| Primary interface | `IBlobSetSource` |
| Input type(s) | `(none)` |
| Output type(s) | `IBlobSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `input` | path or array of paths to input WCP ROOT files (required) |
| `kind` | blob kind to load, "live" or "dead" (default "live") |
| `views` | string of plane letters to include, e.g. "uvw" for 3-view |
| `anode` | type:name of the IAnodePlane component (default "AnodePlane") |
| `frame_eos` | if true, append a null (EOS) IBlobSet marker after each frame's worth of blob sets |
