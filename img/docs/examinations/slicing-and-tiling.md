# Slicing, Tiling, and Blob Clustering

## 1. MaskSlice (`src/MaskSlice.cxx`, `inc/WireCellImg/MaskSlice.h`)

### Algorithm

MaskSlice converts frames into time slices with three plane categories:

**Active planes** (charge-based):
1. Retrieve wiener, charge, error, and summary tagged traces from input frame
2. For each time bin (tick_span ticks wide), accumulate charge per channel
3. Apply adaptive thresholding: a channel is "active" if either:
   - Its wiener signal exceeds `nthreshold * default_threshold` for that plane, OR
   - Its Gaussian-smoothed signal exceeds 1/3 of a neighboring slice's signal
     AND that neighbor exceeds threshold (lines 200-210 in effect)
4. Active channels contribute charge and error to the slice activity map

**Dummy planes** (fill with constant):
- All channels in dummy planes receive `dummy_charge` (default 0) and `dummy_error`
  (default 1e12, marking as unreliable)

**Masked planes** (bad channel masking):
- Uses channel mask map ("bad" tag) to identify dead/masked channels
- Masked channels receive `masked_charge` (0) and `masked_error` (1e12)

The output TracelessFrame wrapper (lines 24-49) provides a minimal IFrame interface
containing only slices, no traces. This is an internal implementation detail.

### Key Configuration Parameters

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `tick_span` | 4 | Ticks per slice |
| `nthreshold` | [3.6, 3.6, 3.6] | Per-plane threshold multiplier |
| `default_threshold` | [2351.3, 3346.6, 2271.9] | Per-plane base RMS x4, UBooNE-specific |
| `dummy_error` | 1e12 | Large error to mark as unreliable |
| `masked_error` | 1e12 | Large error to mark as unreliable |

---

## 2. SumSlice (`src/SumSlice.cxx`)

### Algorithm

Simple accumulation slicer without noise-aware filtering:
1. Groups all non-zero trace samples into time bins of `tick_span` width
2. Accumulates channel charge directly without thresholding
3. Supports single-slice and per-slice output modes
4. Can pad with empty slices if configured

Much simpler than MaskSlice -- no adaptive thresholding, no plane categories.

---

## 3. GridTiling (`src/GridTiling.cxx`, `inc/WireCellImg/GridTiling.h`)

### Algorithm

Converts per-slice channel activity into 3D blobs using the RayGrid framework:

1. For each slice, iterate over wire plane layers (via anode face)
2. Build `Activity` objects per layer: for each wire in the plane, check if
   its channel has activity above threshold in the slice
3. If insufficient layers have activity (< nplanes), return early (no blobs)
4. Call `RayGrid::make_blobs(face, activities)` -- the core blob-finding algorithm
   that finds 3D regions where active wires from all planes intersect
5. Assign unique incrementing blob IDs (`m_blobs_seen++`)
6. Package blobs into `SimpleBlobSet`

The RayGrid blob-finding works by:
- Treating each wire plane as a set of parallel rays
- Finding pairwise strip intersections to define 2D regions
- Intersecting all planes to find 3D volumes ("blobs")
- The `m_nudge` parameter (default 1e-3) provides floating-point robustness

### Key Configuration Parameters

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `threshold` | 0.0 | Minimum channel activity to include |
| `nudge` | 1e-3 | Floating-point correction factor |

---

## 4. BlobGrouping (`src/BlobGrouping.cxx`)

### Algorithm

Adds "measure" nodes to the cluster graph representing electrically connected
blob-channel groups:

1. Per slice, create per-plane subgraphs of blob-channel connections
2. For each of 3 planes (hard-coded), build a subgraph:
   - Vertices: blobs connected to this slice, and channels in this plane
   - Edges: blob-channel connections (blob is connected to channel via wires)
3. Run `boost::connected_components` on each per-plane subgraph
4. Each connected component becomes one `SimpleMeasure`:
   - Signal = sum of channel activity values in the component
   - PlanID from the wire plane
5. Add measure nodes to the cluster graph with edges to their constituent blobs

This creates the bipartite blob-measure structure needed by ChargeSolving.

---

## 5. BlobClustering (`src/BlobClustering.cxx`)

### Algorithm

Frame-level blob clustering:

1. Buffer incoming blob sets until a new frame is detected (frame ident changes)
2. On frame change, sort buffered blob sets by slice time
3. Build cluster graph: add slice, blob, channel, wire nodes with edges
4. Call `geom_clustering()` to form blob-blob edges between adjacent time slices
5. Output one ICluster per frame

The geometric clustering (in GeomClusteringUtil) creates b-b edges based on
RayGrid overlap with configurable tolerance policy.

---

## 6. GeomClusteringUtil (`src/GeomClusteringUtil.cxx`)

### Algorithm

Blob-to-blob association across time slices:

1. For each pair of adjacent blob sets (within `max_rel_diff` time slices apart):
2. For each blob in set1, find overlapping blobs in set2 using RayGrid:
   - `TolerantVisitor` wraps `RayGrid::overlap()` with wire gap tolerance
   - Tolerance from `map_gap_tol[rel_diff]` allows for small wire gaps
3. Add b-b edges for all overlapping pairs

**Special case: `dead_clus` policy**
- Uses `adjacent_dead()` function with hard-coded offset `4*500*us` (line 34)
- This represents 4 ticks at 500us/tick for dead detector regions
- More permissive matching for dead regions

### `grouped_geom_clustering()` variant
- Same algorithm but operates on pre-defined blob groups rather than sequential pairs
- Used by GlobalGeomClustering and LocalGeomClustering for re-clustering operations
