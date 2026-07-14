# PointTreeBuilding

Converts one or two ICluster inputs (live and optionally dead) into a hierarchical point-cloud tree serialised as a tensor set, adding sampled blobs, charge/time projection point clouds, and dead wire-index ranges.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `PointTreeBuilding` |
| Concrete class | `WireCell::Clus::PointTreeBuilding` |
| Node category | fanin |
| Primary interface | `IClusterFaninTensorSet` |
| Input type(s) | `vector<ICluster>` |
| Output type(s) | `ITensorSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `multiplicity` | number of input ICluster streams, 1 (live only) or 2 (live + dead) (default 2) |
| `tags` | array of trace tags, one per input stream |
| `anode` | type/name of the IAnodePlane to use |
| `detector_volumes` | type/name of the IDetectorVolumes service |
| `face` | integer face index on the anode to process (default 0) |
| `samplers` | object mapping point-cloud names to IBlobSampler type/name strings; a "3d" entry is required |
| `datapath` | output tensor datapath template; "%d" replaced by cluster-set ident (default "pointtrees/%d") |
