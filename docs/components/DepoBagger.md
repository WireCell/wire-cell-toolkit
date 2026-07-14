# DepoBagger

Collects all input depositions within a fixed time gate and emits them as a single IDepoSet when EOS is received; if gate is [0,0], all depositions are accepted.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `DepoBagger` |
| Concrete class | `WireCell::Gen::DepoBagger` |
| Node category | queuedout |
| Primary interface | `IDepoCollector` |
| Input type(s) | `IDepo` |
| Output type(s) | `IDepoSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `gate` | two-element array [t_start, t_stop] in WCT time units defining the acceptance window; [0,0] accepts all |
