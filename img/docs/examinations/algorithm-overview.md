# Wire-Cell Imaging Algorithm Overview

## Pipeline Architecture

The `img` package implements the Wire-Cell Toolkit's 3D imaging reconstruction pipeline.
The pipeline transforms 2D wire readout data into 3D spatial objects ("blobs") and
clusters of blobs, then solves for charge and removes ghost artifacts.

The high-level data flow is:

```
Frame -> Slice -> Blob (via Tiling) -> Cluster -> Charge Solve -> Deghost -> Output
```

### Stage 1: Slicing (Frame -> Slices)

Input frames (2D: channel x time) are divided into uniform time windows ("slices"),
each spanning `tick_span` ticks (default 4). Two implementations:

- **MaskSlice** (`MaskSlice.cxx`): Full-featured slicer with adaptive thresholding,
  support for active/dummy/masked planes, and noise-aware filtering. Uses wiener-filtered
  traces for activity detection and charge traces for actual charge values.

- **SumSlice** (`SumSlice.cxx`): Simple accumulation slicer without thresholding.

Each slice contains an activity map: `channel -> (charge, uncertainty)`.

### Stage 2: Tiling (Slices -> Blobs)

**GridTiling** (`GridTiling.cxx`) converts per-slice channel activity into 3D blobs
using the RayGrid framework:

1. Per slice, builds Activity objects for each wire plane layer
2. Calls `RayGrid::make_blobs()` which finds 3D regions where active wires from
   all planes intersect
3. Each blob represents a 3D volume consistent with the observed wire activity

### Stage 3: Blob Operations

- **BlobGrouping** (`BlobGrouping.cxx`): Adds "measure" nodes to the cluster graph.
  Per slice, finds electrically connected groups of blob-channel associations using
  `boost::connected_components`. Each connected component becomes a `SimpleMeasure`.

- **BlobClustering** (`BlobClustering.cxx`): Buffers blob sets frame-by-frame and
  calls `geom_clustering()` to form blob-blob edges between adjacent time slices
  based on geometric overlap (via RayGrid).

- **BlobSetSync/Merge/Fanout**: Stream management nodes for routing blob sets.

### Stage 4: Geometric Clustering

**GeomClusteringUtil** (`GeomClusteringUtil.cxx`) implements blob-to-blob association
across time slices using RayGrid overlap with configurable tolerance. Four policies:

| Policy | max_rel_diff | gap_tol |
|--------|-------------|---------|
| `simple` | 1 | {1: 0} |
| `uboone` | 2 | {1: 2, 2: 1} |
| `uboone_local` | 2 | {1: 2, 2: 2} |
| `dead_clus` | special | uses `adjacent_dead()` |

`max_rel_diff` = max time-slice offset to check. `gap_tol` = wire tolerance per offset.

### Stage 5: Charge Solving

**ChargeSolving** (`ChargeSolving.cxx`) orchestrates LASSO regression to determine
blob charges from wire measurements:

1. **Unpack**: Decompose cluster graph into per-slice bipartite blob-measurement
   subgraphs (`CSGraph::unpack`)
2. **Weight**: Assign blob weights based on cross-slice connectivity. Three strategies:
   - `uniform`: all blobs get weight 9.0
   - `simple`: weight = number of unique connected slices
   - `uboone`: weight 9 / 3^(connections) based on prev/next slice connectivity
3. **Solve**: LASSO regression with optional Cholesky whitening (`CSGraph::solve`)
4. **Prune**: Remove blobs below charge threshold (`CSGraph::prune`)
5. **Repack**: Reconstruct cluster graph with solved charges (`CSGraph::repack`)

Multiple weighting strategies can be applied in sequence, each followed by solve+prune.

### Stage 6: Deghosting

Three deghosting approaches remove spurious 3D artifacts:

- **InSliceDeghosting** (`InSliceDeghosting.cxx`): Local, within-slice ghost removal.
  Identifies blob quality (good/bad) based on charge threshold and cross-slice
  connections, then removes 2-view blobs whose wires are well-covered by 3-view blobs.

- **ProjectionDeghosting** (`ProjectionDeghosting.cxx`): Global deghosting via 2D
  projections. Projects 3D clusters onto each wire plane, compares coverage between
  clusters, and removes clusters whose projections are subsets of higher-quality ones.

- **ShadowGhosting** (`ShadowGhosting.cxx`): Currently a pass-through (incomplete).

### Stage 7: Supporting Components

- **Projection2D** (`Projection2D.cxx`): Builds sparse Eigen matrices representing
  2D (channel x time-slice) projections of 3D blob clusters.

- **FrameQualityTagging** (`FrameQualityTagging.cxx`): Detects noisy periods by
  analyzing activity patterns. Outputs quality flags (0=good through 3=too busy).

- **BlobDepoFill** (`BlobDepoFill.cxx`): Maps simulation deposition charges onto
  blobs using Gaussian-weighted integration.

- **ChargeErrorFrameEstimator** (`ChargeErrorFrameEstimator.cxx`): Estimates charge
  measurement uncertainties based on ROI length and plane-specific fudge factors.

- **FrameMasking** (`FrameMasking.cxx`): Zeros out trace samples within masked time ranges.

- **CMMModifier** (`CMMModifier.cxx`): Modifies channel mask maps for shorted wires,
  veto channels, and dead channels.

---

## Key Data Structures

### ISlice
Time window with activity map: `IChannel -> (charge, uncertainty)`.
Created by MaskSlice or SumSlice.

### IBlob
3D spatial region defined by RayGrid shape (wire intersections).
Has value (charge) and uncertainty. Belongs to one slice and one face.

### Cluster Graph (`cluster_graph_t`)
Boost adjacency_list with variant node types:
- `'b'` (blob): IBlob pointer
- `'s'` (slice): ISlice pointer
- `'m'` (measure): IMeasure pointer (per-plane grouped signal)
- `'w'` (wire): IWire pointer
- `'c'` (channel): IChannel pointer

Edges represent associations: b-s (blob in slice), b-m (blob contributes to measure),
b-w (blob covers wire), w-c (wire is channel), b-b (geometric overlap across slices).

### CS Graph (`CS::graph_t`)
Charge-solving bipartite graph with blob and measurement nodes.
Used internally by CSGraph for LASSO regression. Each node carries a `value_t`
(value + uncertainty) and an `ordering` index.

### Projection2D (`sparse_mat_t`)
Eigen sparse matrix (channel x time-slice) representing 2D projection of a cluster
onto a single wire plane. Uses CSC (compressed sparse column) storage.

---

## Configuration Example (pdhd)

The `pdhd/wct-clustering.jsonnet` configuration wires the pipeline as:

```
ClusterFileSource (active + masked per APA)
    -> per_apa pipeline (BlobSampler, PointTreeBuilding, clustering algorithms)
    -> all_apa pipeline (cross-APA merging, deghosting)
    -> Output
```

Key parameters from `pdhd/clus.jsonnet`:
- `drift_speed`: 1.6 mm/us
- `time_offset`: -250 us
- Per-face and per-APA detector volume definitions with fiducial volume margins
- Pipeline stages: live_dead, extend, regular, parallel_prolong, close, extend_loop,
  separate, connect1, deghost, protect_overclustering, neutrino, isolated, examine_bundles
