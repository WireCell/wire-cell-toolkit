# DuctorFramer

Wraps an IDuctor and an IFrameFanin to produce a single IFrame from a set of depositions, feeding each depo through the ductor and merging resulting frames via the fanin component.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `DuctorFramer` |
| Concrete class | `WireCell::Gen::DuctorFramer` |
| Node category | function |
| Primary interface | `IDepoFramer` |
| Input type(s) | `IDepoSet` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `ductor` | type:name of the IDuctor component to use for simulation |
| `fanin` | type:name of the IFrameFanin component used to merge multiple frames produced by the ductor |
