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

## Output time and length

The output frame takes its ident, time and tick from the frame reached via
`blob->slice()->frame()`.  Output traces are dense and start at tbin 0 (the
frame time).  Their length runs to the end of that frame's traces, if it has
any, and at least to the end of the last slice.  The slice frame may have no
traces, for example MaskSlicer's `TracelessFrame` or a frame rebuilt from a
cluster file.  A cluster with no blobs gives an empty frame.
