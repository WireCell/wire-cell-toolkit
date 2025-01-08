
# PointTreeBuilding Class Analysis

## Overview
PointTreeBuilding is a class for converting clusters of blobs into point cloud trees and output tensors, primarily used in wire-cell data processing. It inherits from three base classes:
- `Aux::Logger` - For logging functionality
- `IClusterFaninTensorSet` - For tensor set operations
- `IConfigurable` - For configuration management

## Core Components

### 1. Key Member Variables

```cpp
// Sampling and configuration
size_t m_multiplicity {2};              // Number of input vectors (default 2)
std::vector<std::string> m_tags;        // Tags for input frames
size_t m_count{0};                      // Processing counter

// Physics parameters
double m_tick {0.5*units::us};          // Time tick size
double m_drift_speed {1.101*units::mm/units::us}; // Drift velocity
double m_time_offset {-1600 * units::us}; // Time offset
double m_dead_threshold {1e10};         // Threshold for dead channels

// Angle parameters for projections
double m_angle_u {1.0472};  // 60 degrees
double m_angle_v {-1.0472}; // -60 degrees
double m_angle_w {0};       // 0 degrees

// Core components
IAnodePlane::pointer m_anode;           // Anode plane interface
std::map<std::string, IBlobSampler::pointer> m_samplers; // Blob samplers
std::string m_datapath = "pointtrees/%d"; // Output data path format
```

### 2. Main Processing Flow

The class processes data through its operator() method with the following steps:

1. Input validation and EOS (End Of Stream) checking
2. Processing live clusters through `sample_live()`
3. Processing dead clusters through `sample_dead()` (if multiplicity = 2)
4. Adding channel-time point clouds (CTPC) through `add_ctpc()`
5. Adding dead wire information through `add_dead_winds()`
6. Converting results to tensors for output

## Key Algorithms

### 1. Live Cluster Sampling (sample_live)

```cpp
Points::node_ptr sample_live(const WireCell::ICluster::pointer icluster) const {
    1. Get geometric clusters from the input cluster graph
    2. Create root node for point tree
    3. For each cluster:
        - Sample 3D points using the "3d" sampler
        - Create 2D projections (2dp0, 2dp1, 2dp2)
        - Calculate blob center
        - Add scalar metadata (charge, margins, etc.)
        - Add to point tree
    4. Return root node
}
```

### 2. Dead Cluster Sampling (sample_dead)

```cpp
Points::node_ptr sample_dead(const WireCell::ICluster::pointer icluster) const {
    1. Get geometric clusters from input graph
    2. Create root node
    3. For each cluster:
        - Create scalar dataset with basic metadata
        - Extract corner points from blob shape
        - Add to point tree
    4. Return root node
}
```

### 3. CTPC Addition (add_ctpc)

This algorithm adds channel-time point clouds to the tree:

1. Extract slice information from cluster graph
2. For each slice:
   - Process activity in channels
   - Calculate positions (x,y) for each wire
   - Group data by face and plane
3. Create datasets for each face/plane combination containing:
   - Position coordinates (x,y)
   - Charge information
   - Channel identifiers
   - Wire indices
   - Slice indices

### 4. Dead Wire Processing (add_dead_winds)

This algorithm processes dead wire information:

1. Scan cluster graph for dead channels (charge uncertainty > threshold)
2. For each dead channel:
   - Calculate drift coordinates (xbeg, xend)
   - Group by face and plane
3. Create datasets containing:
   - Start/end positions
   - Wire indices
   - Dead channel information

## Configuration System

The class uses a JSON-based configuration system with these key parameters:

```json
{
    "multiplicity": 2,
    "tags": ["tag1", "tag2"],
    "datapath": "pointtrees/%d",
    "anode": "AnodePlaneType",
    "samplers": {
        "3d": "BlobSamplerType",
        "dead": "DeadBlobSamplerType"
    }
}
```

## Error Handling

The class implements several error checks:
1. Input validation for multiplicity
2. EOS detection
3. Sampler availability verification
4. Configuration parameter validation

## Key Data Structures

### 1. Point Cloud Tree
- Hierarchical structure representing spatial relationships
- Each node contains:
  - 3D point data
  - 2D projections
  - Scalar metadata
  - Channel/time information

### 2. Tensor Output
- Organized by live/dead classification
- Contains position, charge, and metadata information
- Structured for efficient processing downstream

## Performance Considerations

1. Memory Management:
   - Uses smart pointers for tree nodes
   - Efficient data structure sharing
   - Careful handling of large point clouds

2. Computational Efficiency:
   - Organized processing by face/plane
   - Efficient point cloud sampling
   - Structured data access patterns

## Usage Example

```cpp
// Creating and configuring the processor
auto ptb = make_shared<PointTreeBuilding>();
Configuration cfg;
cfg["multiplicity"] = 2;
cfg["anode"] = "AnodePlane";
cfg["samplers"]["3d"] = "BlobSampler";
ptb->configure(cfg);

// Processing clusters
ICluster::pointer live_cluster = /* input live cluster */;
ICluster::pointer dead_cluster = /* input dead cluster */;
input_vector invec = {live_cluster, dead_cluster};
output_pointer tensorset;
bool success = (*ptb)(invec, tensorset);
```


This class serves as a crucial component in the wire-cell data processing pipeline, transforming blob clusters into structured point cloud trees. The key points to understand are:

1. It handles both "live" and "dead" clusters separately, with different processing strategies for each
2. It creates multiple projections (3D and 2D) of the point clouds
3. It maintains detailed metadata about charges, channels, and wire positions
4. It organizes all data into a hierarchical tree structure for efficient processing

The class is particularly sophisticated in its handling of:
- Multiple coordinate systems (3D space, wire planes, time dimensions)
- Different types of data (live/dead clusters, charges, geometrical information)
- Complex metadata relationships
- Efficient data organization for downstream processing

