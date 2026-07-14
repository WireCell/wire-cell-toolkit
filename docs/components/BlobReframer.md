# BlobReframer

Converts a cluster of blobs back into a frame by projecting each blob's slice activity onto channel waveforms, producing one output frame covering only the channel/time regions spanned by the blobs.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `BlobReframer` |
| Concrete class | `WireCell::Img::BlobReframer` |
| Node category | function |
| Primary interface | `IClusterFramer` |
| Input type(s) | `ICluster` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `frame_tag` | tag string applied to the output frame (default: "reframe") |
