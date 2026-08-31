# ProjectionDeghosting

Removes ghost clusters by comparing 2D wire-plane projections of candidate clusters; clusters whose projections are fully covered by or equivalent to other clusters in the same connected component are tagged and removed.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `ProjectionDeghosting` |
| Concrete class | `WireCell::Img::ProjectionDeghosting` |
| Node category | function |
| Primary interface | `IClusterFilter` |
| Input type(s) | `ICluster` |
| Output type(s) | `ICluster` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `verbose` | enable verbose debug logging (default: false) |
| `nchan` | number of channel bins used when building 2D projections |
| `nslice` | number of slice bins used when building 2D projections |
| `dryrun` | if true, pass input through unchanged (default: false) |
| `uncer_cut` | uncertainty threshold applied when computing projections and coverage |
| `dead_default_charge` | default charge value assigned to dead channels when computing projections |
| `global_deghosting_cut_values` | flat array of cut values for deciding whether a cluster is saved |
| `judge_alt_cut_values` | array of 4 cut values used in the secondary pairwise coverage comparison |
