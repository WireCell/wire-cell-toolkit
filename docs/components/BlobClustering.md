# BlobClustering

Accumulates a stream of IBlobSets and emits one ICluster per frame (on EOS or frame boundary), building a graph with slice, blob, wire, channel, and blob-blob edges using a configurable geometric proximity policy.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `BlobClustering` |
| Concrete class | `WireCell::Img::BlobClustering` |
| Node category | queuedout |
| Primary interface | `IClustering` |
| Input type(s) | `IBlobSet` |
| Output type(s) | `ICluster` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `policy` | clustering policy to use for blob-blob edge formation; one of "simple", "uboone", or "uboone_local" (default: "uboone") |
