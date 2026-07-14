# DepoChunker

Collects depositions into fixed-duration time windows that advance continuously; emits accumulated depos as an IDepoSet when an incoming depo falls outside the current window.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `DepoChunker` |
| Concrete class | `WireCell::Gen::DepoChunker` |
| Node category | queuedout |
| Primary interface | `IDepoCollector` |
| Input type(s) | `IDepo` |
| Output type(s) | `IDepoSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `gate` | two-element array [t_start, t_stop] defining the initial window; the duration is reused for all subsequent windows |
