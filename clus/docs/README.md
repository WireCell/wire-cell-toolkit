# Pattern Recognition Documentation

This directory contains documentation for the `clus/` Pattern Recognition (PR)
codebase. Start with the overview, then follow the links to each deep-dive document.

## Documents

| File | Contents |
|---|---|
| [overview.md](overview.md) | Architecture, visitor pattern, data model, pipeline summary |
| [pipeline_stages.md](pipeline_stages.md) | Each of the 6 IEnsembleVisitor stages in detail |
| [steiner_graph.md](steiner_graph.md) | Steiner tree construction from blob point clouds |
| [pattern_recognition.md](pattern_recognition.md) | `find_proto_vertex` and the full PR loop |
| [track_shower_separation.md](track_shower_separation.md) | Classifying particles, shower clustering, kinematics |
| [track_fitting.md](track_fitting.md) | Calorimetry, dE/dx, wire-coordinate fitting |
| [data_structures.md](data_structures.md) | Facade layer, PR::Graph, PR::Vertex, PR::Segment |

## Quick Reference: Pipeline Stages

```
Input: PointTree blobs (wire charge depositions)
    │
    ▼ ClusteringTaggerFlagTransfer   transfer upstream tagger metadata to cluster flags
    ▼ ClusteringRecoveringBundle     undo over-merging from beam-flash association
    ▼ ClusteringSwitchScope          apply T0/space-charge coordinate correction
    ▼ CreateSteinerGraph             build minimal spanning tree skeleton per cluster
    ▼ MakeFiducialUtils              attach fiducial/dead-wire geometry queries
    ▼ TaggerCheckNeutrino            find vertex, classify particles, compute kinematics
    │
    ▼
Output: PR::Graph (vertices + track/shower segments with kinematics)
```

## Key Source Directories

```
clus/inc/WireCellClus/     C++ headers (interfaces, data types)
clus/src/                  C++ implementations
cfg/pgrapher/common/       jsonnet factory (clus.jsonnet defines the cm object)
cfg/pgrapher/experiment/   detector-specific configurations
```
