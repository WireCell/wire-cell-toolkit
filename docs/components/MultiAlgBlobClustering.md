# MultiAlgBlobClustering

Runs a configurable pipeline of blob-clustering algorithms on live/dead point-cloud groupings read from an input tensor set and writes the clustered result back as tensors, optionally exporting Bee visualisation JSON.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `MultiAlgBlobClustering` |
| Concrete class | `WireCell::Clus::MultiAlgBlobClustering` |
| Node category | function |
| Primary interface | `ITensorSetFilter` |
| Input type(s) | `ITensorSet` |
| Output type(s) | `ITensorSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `groupings` | list of grouping names to load and process (default ["live","dead"]) |
| `inpath` | base datapath pattern for input pc-tree tensors; "%d" replaced by tensor-set ident |
| `outpath` | base datapath pattern for output pc-tree tensors; "%d" replaced by tensor-set ident |
| `pipeline` | ordered array of IEnsembleVisitor type/name strings defining the clustering algorithm pipeline |
| `anodes` | array of IAnodePlane type/name strings for the detector anodes to use |
| `detector_volumes` | type/name of the IDetectorVolumes service |
| `bee_zip` | path of the output Bee zip file (default "mabc.zip") |
| `save_deadarea` | if true, save dead-channel area patches to the Bee zip (default false) |
| `dead_live_overlap_offset` | offset in slices for dead/live overlap clustering (default 2) |
| `use_config_rse` | if true, use runNo/subRunNo/eventNo from config instead of tensor-set ident (default false) |
| `runNo` | run number used when use_config_rse is true |
| `subRunNo` | sub-run number used when use_config_rse is true |
| `eventNo` | event number used when use_config_rse is true |
| `cluster_id_order` | cluster ID ordering after each pipeline step: "tree", "size", or "" (default "") |
| `perf` | if true, emit time/memory performance log messages (default false) |
