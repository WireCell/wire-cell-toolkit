# JsonBlobSetSink

Writes each blob set to a JSON file containing blob corner coordinates (converted to drift-x, y, z positions) and associated charge/uncertainty values.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `JsonBlobSetSink` |
| Concrete class | `WireCell::Img::JsonBlobSetSink` |
| Node category | sink |
| Primary interface | `IBlobSetSink` |
| Input type(s) | `IBlobSet` |
| Output type(s) | `(none)` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `filename` | output file name pattern; a "%" will be substituted with the blob set ident number (default: "blobs-%02d.json") |
| `face` | face ident filter; set to -1 to accept all faces (default: 0) |
| `drift_speed` | drift velocity used to convert slice start time to x-coordinate (default: 1.6 mm/us) |
