# LCBlobRemoval Class Documentation

## Overview

The `LCBlobRemoval` class is a component of the Wire-Cell Toolkit's image processing module (`WireCellImg`). It functions as a filter for removing low-charge blobs from a cluster, effectively reducing noise and simplifying downstream processing.

## Class Definition

```cpp
class lcblobremoval : public aux::logger, public iclusterfilter, public iconfigurable {
public:
    lcblobremoval();
    virtual ~lcblobremoval();
    virtual void configure(const wirecell::configuration& cfg);
    virtual wirecell::configuration default_configuration() const;
    virtual bool operator()(const input_pointer& in, output_pointer& out);
private:
    // used to hold measurement and blob values
    // (central+uncertainty).
    using value_t = islice::value_t;
    
    // config: blob_{value,error}_threshold. blob
    // central value less dropped.
    // uncertainty currently considered.
    value_t m_blob_thresh{0,1};
};
```

## Class Inheritance

- **aux::logger**: Provides logging capabilities
- **iclusterfilter**: Defines the interface for filters that process clusters
- **iconfigurable**: Allows the class to be configured via JSON configuration

## Key Components

### Configuration Parameters

The class accepts the following configuration parameters:

- **blob_value_threshold**: The minimum charge value a blob must have to be kept (default: 0)
- **blob_error_threshold**: Related to uncertainty of measurements (default: 1, but not directly used in the filtering logic)

Example configuration:

```json
{
  "blob_value_threshold": 300,
  "blob_error_threshold": 1
}
```

### Default Configuration

```cpp
wirecell::configuration img::lcblobremoval::default_configuration() const {
    wirecell::configuration cfg;
    cfg["blob_value_threshold"] = m_blob_thresh.value();
    cfg["blob_error_threshold"] = m_blob_thresh.uncertainty();
    return cfg;
}
```

### Main Operation

The main operation is performed in the `operator()` method:

```cpp
bool img::lcblobremoval::operator()(const input_pointer& in, output_pointer& out) {
    out = nullptr;
    if(!in) {
        log->debug("eos");
        return true;
    }
    
    const auto in_graph = in->graph();
    dump_cg(in_graph, log);
    
    auto out_graph = prune(in_graph, m_blob_thresh.value());
    dump_cg(out_graph, log);
    
    out = std::make_shared<aux::simplecluster>(out_graph, in->ident());
    return true;
}
```

The key processing is in the `prune` function:

```cpp
cluster_graph_t prune(const cluster_graph_t& cg, float threshold) {
    cluster_graph_t cg_out;
    size_t nblobs = 0;
    std::unordered_map<cluster_vertex_t, cluster_vertex_t> old2new;
    
    // Iterate through all vertices
    for(const auto& vtx : mir(boost::vertices(cg))) {
        const auto& node = cg[vtx];
        if(node.code() == 'b') {
            const auto iblob = get<blob_t>(node.ptr);
            auto bval = iblob->value();
            if(bval < threshold) continue;  // Skip low-charge blobs
        }
        ++nblobs;
        old2new[vtx] = boost::add_vertex(node, cg_out);
    }
    
    // Copy edges between remaining vertices
    for(auto edge : mir(boost::edges(cg))) {
        auto old_tail = boost::source(edge, cg);
        auto old_head = boost::target(edge, cg);
        auto old_tit = old2new.find(old_tail);
        if(old_tit == old2new.end()) {
            continue;
        }
        auto old_hit = old2new.find(old_head);
        if(old_hit == old2new.end()) {
            continue;
        }
        boost::add_edge(old_tit->second, old_hit->second, cg_out);
    }
    
    return cg_out;
}
```

## Cluster Graph and Faces in Wire-Cell

The `LCBlobRemoval` class operates within the Wire-Cell Toolkit's graph-based representation of detector data:

### Cluster Graph Structure

- **Vertices**: Represent different entities including:
  - Blobs ('b'): 2D projections of charge deposits
  - Slices ('s'): Time windows
  - Wires ('w'): Individual detector wires
  - Channels ('c'): Electronics channels
  - Measures ('m'): Measurement nodes

- **Edges**: Represent relationships between entities:
  - Blob-Slice: A blob exists within a time slice
  - Blob-Wire: A blob covers certain wires
  - Wire-Channel: Wires connect to channels
  - Blob-Blob: Spatial or temporal relationships between blobs

### The "Face" Concept

Although `LCBlobRemoval` doesn't directly manipulate face information, understanding faces is important:

1. **Anode Plane Assembly (APA)**:
   - Physical detector component with sensing surfaces

2. **Face**:
   - A specific sensing surface of an APA
   - Contains multiple wire planes (typically U, V, W)
   - Each blob is associated with a specific face

3. **How Faces Are Used**:
   - Blobs are associated with faces: `blob->face()`
   - This association is preserved during filtering
   - The face determines the coordinate system of the blob
   - Different faces might have different geometries or calibrations

## Usage in Processing Pipeline

The `LCBlobRemoval` class is typically used in a processing pipeline where:

1. Detector data is converted to clusters of blobs
2. `LCBlobRemoval` filters out low-charge blobs
3. Subsequent components process the filtered clusters

Example configuration in a pipeline:

```json
{
  "configs": [
    {
      "type": "LCBlobRemoval",
      "data": {
        "blob_value_threshold": 300
      }
    }
  ],
  "connections": [
    {
      "input": "SomeBlobGenerator:output",
      "output": "LCBlobRemoval:input"
    },
    {
      "input": "LCBlobRemoval:output",
      "output": "NextProcessor:input"
    }
  ]
}
```

## Summary

The `LCBlobRemoval` class serves as a simple but essential filter in the Wire-Cell Toolkit's image processing module. It removes blobs with charge values below a configurable threshold, thereby reducing noise and simplifying downstream processing.

While the class itself doesn't directly interact with the "face" concept, it preserves the face associations of blobs that pass the filter. Understanding the role of faces in the Wire-Cell framework is important for comprehending how blob filtering fits into the overall data processing pipeline.