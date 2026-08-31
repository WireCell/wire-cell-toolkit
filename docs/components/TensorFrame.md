# TensorFrame

Converts an ITensorSet to an IFrame by locating the frame tensor at a configurable datapath regex and applying an optional linear amplitude decoding transform.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `TensorFrame` |
| Concrete class | `WireCell::Aux::TensorFrame` |
| Node category | function |
| Primary interface | `ITensorSetFrame` |
| Input type(s) | `ITensorSet` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `dpre` | Regular expression matched against tensor datapaths to find the frame tensor; default "frames/[[:digit:]]+$". |
| `baseline` | Baseline subtracted after decoding; default 0. |
| `scale` | Divisor applied during decoding; default 1. |
| `offset` | Value subtracted from encoded sample before dividing by scale; default 0. |
