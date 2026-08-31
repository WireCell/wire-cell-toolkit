# JsonDepoSource

Loads energy depositions from a JSON file, converting charge using a configurable recombination model, and emits them as IDepo objects in time order.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `JsonDepoSource` |
| Concrete class | `WireCell::Sio::JsonDepoSource` |
| Node category | source |
| Primary interface | `IDepoSource` |
| Input type(s) | `(none)` |
| Output type(s) | `IDepo` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `filename` | path to the input JSON file (required) |
| `jsonpath` | dot-separated path to the deposition array within the JSON structure (default "depos") |
| `model` | recombination model type:name or "electrons" for direct electron-count input (default "electrons") |
| `scale` | scale factor applied when using the "electrons" model (default 1.0) |
