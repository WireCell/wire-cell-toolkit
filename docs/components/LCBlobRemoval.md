# LCBlobRemoval

Removes blobs from a cluster graph whose charge value falls below a configurable threshold, pruning low-charge blobs and their associated edges from the output.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `LCBlobRemoval` |
| Concrete class | `WireCell::Img::LCBlobRemoval` |
| Node category | function |
| Primary interface | `IClusterFilter` |
| Input type(s) | `ICluster` |
| Output type(s) | `ICluster` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `blob_value_threshold` | minimum blob charge value required to keep a blob |
| `blob_error_threshold` | minimum blob charge uncertainty threshold |
