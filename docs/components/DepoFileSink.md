# DepoFileSink

Writes IDepoSet contents to a streaming archive file as pairs of NumPy arrays (depo_data and depo_info) encoding deposition kinematics and metadata, one pair per IDepoSet received.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `DepoFileSink` |
| Concrete class | `WireCell::Sio::DepoFileSink` |
| Node category | sink |
| Primary interface | `IDepoSetSink` |
| Input type(s) | `IDepoSet` |
| Output type(s) | `(none)` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `outname` | output archive file path; extension determines compression format (required) |
