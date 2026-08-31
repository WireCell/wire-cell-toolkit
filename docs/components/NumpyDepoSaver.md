# NumpyDepoSaver

Accumulates input depositions and on EOS saves them as paired NumPy 2D arrays (depo_data and depo_info) appended to a NumPy zip file, then passes depositions through unchanged.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `NumpyDepoSaver` |
| Concrete class | `WireCell::Sio::NumpyDepoSaver` |
| Node category | function |
| Primary interface | `IDepoFilter` |
| Input type(s) | `IDepo` |
| Output type(s) | `IDepo` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `filename` | path to the output NumPy zip file; writing is always in append mode |
