# PointTreeMerging

Merges multiple per-APA point-cloud trees (live and dead sub-trees) arriving as tensor sets into a single combined tensor set.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `PointTreeMerging` |
| Concrete class | `WireCell::Clus::PointTreeMerging` |
| Node category | fanin |
| Primary interface | `ITensorSetFanin` |
| Input type(s) | `vector<ITensorSet>` |
| Output type(s) | `ITensorSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `multiplicity` | number of input ITensorSet streams to merge |
| `inpath` | base datapath for the input live/dead pc-tree tensors; "%d" replaced by tensor-set ident |
| `outpath` | base datapath for the merged output pc-tree tensors; "%d" replaced by tensor-set ident |
