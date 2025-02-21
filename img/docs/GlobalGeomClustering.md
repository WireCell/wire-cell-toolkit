# GlobalGeomClustering

## Overview

`GlobalGeomClustering` is a component in Wire-Cell Toolkit's image processing framework that performs geometric clustering of blobs based on their spatial proximity. The class is designed to create connections between charge blobs that likely belong to the same physical particle track or shower.

## Purpose

The primary purpose of this class is to:

1. Take existing clusters of charge blobs
2. Remove existing blob-to-blob connections
3. Create new connections based on geometric proximity criteria
4. Return an updated cluster with improved blob connectivity

## Class Definition

```cpp
class globalgeomclustering : public aux::logger, 
                           public iclusterfilter, 
                           public iconfigurable {
public:
    globalgeomclustering();
    virtual ~globalgeomclustering();
    virtual void configure(const wirecell::configuration& cfg);
    virtual wirecell::configuration default_configuration() const;
    virtual bool operator()(const input_pointer& in, output_pointer& out);
private:
    std::string m_clustering_policy{"uboone"};
};
```

## Key Concepts

### Anode Plane and Faces

In Wire-Cell, detector geometry is represented by:

- **AnodePlane**: An abstraction of the physical wire planes in a Time Projection Chamber (TPC) detector
- **Face**: A subset of the anode that collects charge drifting from a particular direction
  - In a single-sided TPC: one face
  - In a double-sided TPC: two faces (front and back)
  - Each face has multiple wire planes (typically U, V, W in LArTPCs) at different angles

The `GlobalGeomClustering` operates on blobs within the same face. It does not connect blobs across different faces, as these would represent charges originating from opposite directions.

### Blobs and Clusters

- **Blob**: A contiguous region of charge deposition in 2D (wire-time) space
- **Cluster**: A collection of blobs with various connection types:
  - blob-to-slice: connecting blobs to their time slices
  - wire-to-channel: connecting wires to readout channels
  - blob-to-blob: connecting blobs that are likely part of the same physical object

This class specifically modifies the blob-to-blob connections based on geometric proximity.

## Configuration Parameters

| Parameter | Description | Default Value |
|-----------|-------------|---------------|
| `clustering_policy` | Determines the algorithm and parameters for establishing blob-to-blob connections | "uboone" |

### Available Clustering Policies

1. **"uboone"**: Parameters tuned for the MicroBooNE detector
   - Maximum relative time difference between slices: 2 ticks
   - Gap tolerance for adjacent slices: {1:2, 2:1} (time diff : wire tolerance)

2. **"simple"**: A simpler clustering algorithm
   - Maximum relative time difference: 1 tick
   - Gap tolerance: {1:0}

3. **"uboone_local"**: Similar to "uboone" but with modified wire tolerance
   - Maximum relative time difference: 2 ticks
   - Gap tolerance: {1:2, 2:2}

## Implementation Details

1. The `operator()` method takes an input cluster and produces an output cluster:
   - First, it examines blob-blob connectivity in the original cluster
   - It creates a filtered graph that excludes blob-blob edges
   - It calls `grouped_geom_clustering()` to create new edges based on geometry
   - It returns a new cluster with the updated connectivity

2. The geometric clustering considers:
   - Spatial proximity of blobs
   - Time proximity of slices
   - Wire plane topology

3. Only blobs within the same face can be connected, as cross-face connections would require additional transformation logic.

## Usage Example

```cpp
// Configuration JSON
{
    "clustering_policy": "uboone"
}

// In a WCT configuration file
local clusterer = {
    type: "GlobalGeomClustering",
    data: {
        clustering_policy: "uboone"
    }
}
```

## Internal Operation Flow

1. Receives an input cluster
2. Filters the cluster to remove existing blob-blob edges
3. Creates a new graph without these edges
4. Calls `grouped_geom_clustering()` with the specified policy
5. Returns the cluster with new edges based on geometric proximity

## Limitations

- Only operates on blobs within the same face
- Does not connect blobs across different faces
- Uses predefined policies with fixed parameters
- Does not handle time transformations between faces

## References

1. Wire-Cell Toolkit: [https://github.com/WireCell/wire-cell-toolkit](https://github.com/WireCell/wire-cell-toolkit)
2. MicroBooNE Experiment: [https://microboone.fnal.gov/](https://microboone.fnal.gov/)