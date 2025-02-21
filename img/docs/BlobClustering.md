# BlobClustering Class Documentation

## Overview

The BlobClustering class is a key component in the Wire-Cell Toolkit's image processing framework for particle physics detector data. It takes collections of "blobs" (spatial regions of charge) and organizes them into clusters based on spatial and temporal relationships. This clustering is a crucial step in particle track reconstruction.

## Purpose

The main purpose of BlobClustering is to:

1. Assemble blobs from different time slices into coherent particle tracks
2. Establish relationships between detector components (slices, blobs, wires, channels)
3. Apply geometric clustering to connect blobs that are likely from the same particle

## Key Concepts

### Blob

A blob represents a region of charge in 3D space, detected by wire planes. It has:
- A shape (geometric boundaries)
- Associated slice (time information)
- Connection to a specific detector face
- A charge value and uncertainty

### Detector Geometry

The detector hierarchy is:
- **Anode Plane**: The overall detector component
  - **Faces**: Sub-sections of the anode plane
    - **Wire Planes**: Layers of wires at different angles
      - **Wires**: Individual sensing elements
        - **Channels**: Electronics connected to wires

### Clustering

Clustering establishes relationships between:
- Slices and blobs (temporal connection)
- Blobs and wires (spatial connection)
- Wires and channels (detector structure)
- Blobs and blobs (spatial proximity)

## Implementation Details

### Class Definition

```cpp
class blobclustering : public aux::logger, public iclustering, public iconfigurable {
public:
    blobclustering();
    virtual ~blobclustering();
    virtual void configure(const wirecell::configuration& cfg);
    virtual wirecell::configuration default_configuration() const;
    virtual bool operator()(const input_pointer& blobset, output_queue& clusters);

private:
    std::string m_policy{"uboone"};
    iblobset::vector m_cache;
    void flush(output_queue& clusters);
    bool graph_bs(const input_pointer& newbs);
    bool new_frame(const input_pointer& newbs) const;
    int cur_ident() const;
    int m_count{0};
};
```

### Key Methods

#### `operator()`

This is the main processing function that gets called when new blob sets arrive:
- Handles EOS (End Of Stream) signals
- Detects when a new frame begins
- Caches blob sets until a frame is complete
- Calls `flush()` to process accumulated blob sets

#### `flush()`

Processes all cached blob sets:
1. Sorts blob sets by time
2. Creates graph connections between slices, blobs, wires, and channels
3. Performs geometric clustering to connect related blobs
4. Produces a cluster and clears the cache

#### `add_blobs()` (Helper Function)

Sets up the graph structure:
- Connects slices to blobs
- Connects blobs to wires based on their geometry and the face they belong to
- Establishes the basic topology for clustering

#### `geom_clustering()` (Helper Function)

Creates blob-to-blob connections based on geometric proximity:
- Uses the configured policy (e.g., "simple", "uboone")
- Applies different tolerances for connecting blobs based on their time separation
- Considers spatial overlap in wire plane coordinates

## Configuration

The class accepts these configuration parameters:

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `policy` | string | Geometric clustering policy (e.g., "simple", "uboone") | "uboone" |

## Usage Example

In a Wire-Cell Toolkit configuration:

```jsonnet
local bc = {
  type: "BlobClustering",
  data: {
    policy: "uboone",
  }
};

// Pipeline components
local input = {...};  // source of blob sets
local clustering = bc;
local output = {...};  // next processing step

local pipeline = [
  input,
  clustering,
  output
];
```

## The Role of Anode Face

The anode face is a critical concept for BlobClustering because:

1. **Coordinate System**: Each face provides its own coordinate system via the `raygrid` object
2. **Wire Planes**: The face organizes multiple wire planes that detect charge from different angles
3. **Containment**: Blobs are associated with a specific face, keeping clustering local to detector sections
4. **Geometric Relationships**: Blob-to-wire mapping is done through the face's geometry

When processing blobs, the code accesses the face via `iblob->face()` and uses it to:
- Get the wire planes: `auto wire_planes = iface->planes()`
- Map blob strips to wires in these planes
- Establish connections in the graph

## Typical Processing Flow

1. Blob sets are received, each associated with a time slice
2. Sets are cached until a frame boundary is detected
3. When a frame is complete:
   - Slice-blob-wire-channel connections are established
   - Geometric clustering identifies blob-blob connections
   - A cluster is formed and output
4. The process repeats for the next frame

## Debugging

The class provides detailed debug logging:
- Reports blob counts and graph statistics
- Logs when flushing occurs and how many clusters are produced
- Helps trace the flow of data through the clustering process

## Common Issues

1. **Missing Slice References**: Blobs must have valid slice references
2. **Inconsistent Face References**: Blobs should reference consistent detector faces
3. **Configuration Mismatch**: The clustering policy should match detector geometry

## Relationship with Other Components

- **Upstream**: Receives blob sets from blobbing algorithms
- **Downstream**: Provides clusters to later analysis stages
- **Related Classes**: Works with `GeomClusteringUtil` for blob-blob connections