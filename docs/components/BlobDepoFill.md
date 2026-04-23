# BlobDepoFill

Joins a cluster with a depo set to produce a new cluster whose blob charge values are filled from the Gaussian-weighted integral of ionisation depositions mapped onto each blob's spatial volume.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `BlobDepoFill` |
| Concrete class | `WireCell::Img::BlobDepoFill` |
| Node category | join |
| Primary interface | `IBlobDepoFill` |
| Input type(s) | `(ICluster, IDepoSet)` |
| Output type(s) | `ICluster` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `time_offset` | time offset added to depo times before associating with blobs (default: 0) |
| `nsigma` | number of sigma of each depo Gaussian to consider (default: 3.0) |
| `speed` | nominal drift speed used to convert a depo's longitudinal extent from distance to time units (required) |
| `pindex` | wire-plane index treated as the "primary" plane for dicing depos onto wires (default: 2) |
