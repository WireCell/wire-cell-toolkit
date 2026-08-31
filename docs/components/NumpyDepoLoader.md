# NumpyDepoLoader

Loads energy depositions from a NumPy zip (.npz) file in WCT depo format, streaming depositions event-by-event and sending EOS after each set.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `NumpyDepoLoader` |
| Concrete class | `WireCell::Sio::NumpyDepoLoader` |
| Node category | source |
| Primary interface | `IDepoSource` |
| Input type(s) | `(none)` |
| Output type(s) | `IDepo` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `filename` | path to the input NumPy zip file (default "wct-frame.npz") |
