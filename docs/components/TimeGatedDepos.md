# TimeGatedDepos

Passes or rejects depositions based on whether their time falls within a configurable time gate; the gate can advance by a period at each EOS, allowing periodic gating across multiple readout windows.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `TimeGatedDepos` |
| Concrete class | `WireCell::Gen::TimeGatedDepos` |
| Node category | queuedout |
| Primary interface | `IDrifter` |
| Input type(s) | `IDepo` |
| Output type(s) | `IDepo` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `accept` | "accept" to pass depositions inside the gate, any other value to reject them (default: "accept") |
| `start` | gate start time (default: 0.0) |
| `duration` | gate duration (default: 0.0) |
| `period` | amount added to start at each EOS, enabling periodic gates (default: 0.0) |
