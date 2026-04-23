# NumpyDepoSetLoader

Loads energy depositions from a NumPy zip (.npz) file in WCT depo format, emitting one IDepoSet per stored event set and sending EOS when all sets have been read.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `NumpyDepoSetLoader` |
| Concrete class | `WireCell::Sio::NumpyDepoSetLoader` |
| Node category | source |
| Primary interface | `IDepoSetSource` |
| Input type(s) | `(none)` |
| Output type(s) | `IDepoSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `filename` | path to the input NumPy zip file (default "wct-depos.npz") |
