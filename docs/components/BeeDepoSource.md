# BeeDepoSource

Reads energy depositions from one or more Bee-format JSON files (with x, y, z, t, q arrays) and emits them as IDepo objects; by default sends an EOS null after each file's set.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `BeeDepoSource` |
| Concrete class | `WireCell::Sio::BeeDepoSource` |
| Node category | source |
| Primary interface | `IDepoSource` |
| Input type(s) | `(none)` |
| Output type(s) | `IDepo` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `filelist` | array of input Bee JSON file paths (required) |
| `policy` | streaming policy; set to "stream" to suppress per-file EOS markers (default "") |
