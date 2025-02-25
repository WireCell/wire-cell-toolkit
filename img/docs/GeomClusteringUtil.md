# GeomClusteringUtil in Wire-Cell Toolkit

## Overview

GeomClusteringUtil is a component in the Wire-Cell Toolkit's imaging module (`WireCellImg`) responsible for grouping blobs together based on their spatial and temporal relationships. This document explains its core functionality, the concept of AnodePlane faces, and how the clustering mechanism works.

## Key Files

- `inc/WireCellImg/GeomClusteringUtil.h` - Header file defining the interfaces
- `src/GeomClusteringUtil.cxx` - Implementation file with the core algorithms

## Core Functions

### `geom_clustering()`

This function establishes geometric connections between blobs across adjacent time slices:

```cpp
void wirecell::img::geom_clustering(
    cluster_indexed_graph_t& grind,
    iblobset::vector::iterator beg, 
    iblobset::vector::iterator end, 
    std::string policy
)
```

Parameters:
- `grind` - The indexed graph that will store connections between blobs
- `beg` and `end` - Iterators defining a range of blob sets to process
- `policy` - String specifying which clustering policy to use

### `grouped_geom_clustering()`

This function is similar to `geom_clustering()` but respects pre-established blob groupings:

```cpp
void wirecell::img::grouped_geom_clustering(
    cluster_graph_t& cg, 
    std::string policy, 
    const std::unordered_map<cluster_vertex_t, int> groups
)
```

Parameters:
- `cg` - The cluster graph to update with new blob-blob connections
- `policy` - String specifying which clustering policy to use
- `groups` - Map of blob vertex descriptors to group IDs. Only blobs in the same group can be connected.

## Clustering Policies

The code supports several pre-defined clustering policies:

1. **"simple"** 
   - Maximum relative time difference: 1 time slice
   - Grid tolerance: 0 (no tolerance for spatial overlap)
   - Used for basic clustering with strict requirements

2. **"uboone"**
   - Maximum relative time difference: 2 time slices
   - Grid tolerance: {1→2, 2→1} (more tolerance for adjacent time slices)
   - Optimized for MicroBooNE detector characteristics

3. **"uboone_local"**
   - Maximum relative time difference: 2 time slices
   - Grid tolerance: {1→2, 2→2} (consistent higher tolerance)
   - Requires time slices to be adjacent in the ordered set of slice times

4. **"dead_clus"**
   - Special handling for regions with dead channels
   - Uses `adjacent_dead()` function to determine time adjacency
   - Grid tolerance: {0→1, 1→1} (tolerance even for same-time slices)

## The "Face" Concept in Wire-Cell Toolkit

### What is a Face?

In the Wire-Cell Toolkit, an AnodePlane represents a physical detector plane, and each AnodePlane can have multiple "faces":

1. A face represents one side of a detector plane where charge can be collected
2. Typical TPCs have two faces (front and back)
3. Each face contains multiple wire planes (typically 3 planes: U, V, and W views)

### Face Structure

A face provides:

1. A coordinate system via its `raygrid`
2. Access to the wire planes through `iface->planes()`
3. A unique identifier via `iface->ident()`

### How Faces Impact Clustering

Faces are fundamental to blob clustering because:

1. Each blob is associated with a specific face
2. The face's coordinate system (raygrid) is used to determine spatial relationships
3. Clustering only occurs between blobs on the same face
4. Wire indices and channel numbers are face-specific

In the code, you can see face-specific handling in the `grouped_geom_clustering()` function where blobs from different slices are connected based on their face-specific coordinates.

## Blob Association & Overlap

The core mechanism for determining if blobs should be connected uses a "tolerant visitor" pattern:

```cpp
struct tolerantvisitor {
    raygrid::grid_index_t tolerance{0};
    bool verbose{false};
    raygrid::blobvec_t operator()(const raygrid::blobref_t& blob, 
                                  const raygrid::blobproj_t& proj, 
                                  raygrid::layer_index_t layer) {
        return raygrid::overlap(blob, proj, layer, tolerance, verbose);
    }
};
```

This visitor is used with `raygrid::associate()` to find overlapping blobs and create edges between them.

The overlap calculation takes into account:
1. The wire index bounds in each plane 
2. The time relationship between slices
3. The policy-specific tolerance for gaps

## Practical Example

When processing detector data:

1. Charge depositions are first converted to "blobs" associated with specific detector faces
2. GeomClusteringUtil connects these blobs across time slices based on their spatial overlap
3. The resulting clusters represent potential particle tracks or showers
4. Different policies allow tuning for different detector configurations or reconstruction goals

## Usage in the Wire-Cell Processing Chain

GeomClusteringUtil is typically used in these components:

1. `LocalGeomClustering` - Groups blobs within local regions
2. `GlobalGeomClustering` - Connects clusters across larger regions
3. `InSliceDeghosting` - Uses clustering to identify and remove "ghost" hits

The output is a graph where vertices represent blobs and edges represent their geometric connections, forming the basis for further analysis and particle identification.