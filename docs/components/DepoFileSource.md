# DepoFileSource

Reads IDepoSet objects from a streaming archive file containing paired NumPy arrays (depo_data and depo_info), reconstructing depositions with prior-depo relationships and applying an optional charge scale factor.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `DepoFileSource` |
| Concrete class | `WireCell::Sio::DepoFileSource` |
| Node category | source |
| Primary interface | `IDepoSetSource` |
| Input type(s) | `(none)` |
| Output type(s) | `IDepoSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `inname` | input archive file path (required) |
| `scale` | multiplicative factor applied to each deposition charge (default 1.0) |
