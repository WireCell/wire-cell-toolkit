# BlobSetReframer

Projects each blob in a BlobSet into frame traces, distributing blob charge uniformly over its wire and time extent, with optional per-anode trace tagging and selectable charge source.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `BlobSetReframer` |
| Concrete class | `WireCell::Img::BlobSetReframer` |
| Node category | function |
| Primary interface | `IBlobSetFramer` |
| Input type(s) | `IBlobSet` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `frame_tag` | tag applied to the output frame (default: "") |
| `tick` | output frame sampling period in time units (default: 0.5 us) |
| `source` | selects charge source; "blob" uses blob value/error, "activity" uses slice channel activity (default: "blob") |
| `measure` | selects which scalar to use as charge; "value" or "error" (default: "value") |
