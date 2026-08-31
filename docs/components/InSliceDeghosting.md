# InSliceDeghosting

Removes ghost blobs within each time slice by scoring two-view blobs against nearby three-view blobs using wire overlap and charge quality metrics, with configurable rounds of deghosting logic.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `InSliceDeghosting` |
| Concrete class | `WireCell::Img::InSliceDeghosting` |
| Node category | function |
| Primary interface | `IClusterFilter` |
| Input type(s) | `ICluster` |
| Output type(s) | `ICluster` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `dryrun` | if true, pass the input cluster through unchanged (default: false) |
| `good_blob_charge_th` | charge threshold above which a blob is considered "good" |
| `clustering_policy` | string name of the geometric clustering policy for re-clustering after ghost removal |
| `config_round` | which round of deghosting algorithm to apply (1, 2, or 3) |
| `deghost_th` | wire-overlap ratio threshold used in round-2 deghosting |
| `deghost_th1` | per-plane score threshold below which a blob is flagged as potentially bad |
