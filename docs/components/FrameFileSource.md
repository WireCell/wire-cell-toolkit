# FrameFileSource

Reads IFrame objects from a streaming archive file containing NumPy arrays, reconstructing tagged frames including channel mask maps.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `FrameFileSource` |
| Concrete class | `WireCell::Sio::FrameFileSource` |
| Node category | source |
| Primary interface | `IFrameSource` |
| Input type(s) | `(none)` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `inname` | input archive file path (required) |
| `tags` | array of trace tags to load; if empty, all tagged and untagged traces are loaded |
| `frame_tags` | array of frame-level tags to apply to every output frame |
