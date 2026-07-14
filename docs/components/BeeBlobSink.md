# BeeBlobSink

Samples blobs from each IBlobSet using one or more IBlobSampler components and writes the resulting 3D point data to a Bee JSON zip file, grouping points by frame ident to form Bee "events".

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `BeeBlobSink` |
| Concrete class | `WireCell::Sio::BeeBlobSink` |
| Node category | sink |
| Primary interface | `IBlobSetSink` |
| Input type(s) | `IBlobSet` |
| Output type(s) | `(none)` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `outname` | output file path for the Bee zip/JSON file (required) |
| `geom` | canonical Bee detector name string (required) |
| `type` | algorithm label written into the Bee file (default "unknown") |
| `run` | run number for Bee RSE metadata (default 0) |
| `sub` | subrun number for Bee RSE metadata (default 0) |
| `evt` | event number offset added to frame ident for Bee RSE metadata (default 0) |
| `samplers` | type:name string or array of type:name strings for IBlobSampler components (required) |
