# Pattern Recognition Codebase Overview

## Introduction

This document provides a high-level overview of the Wire-Cell Pattern Recognition (PR)
codebase in `clus/`. The goal of this system is to take raw charge measurements from a
Liquid Argon Time Projection Chamber (LArTPC) — already organized into 3D "blobs" by
the upstream imaging stage — and reconstruct the physical event: identify the neutrino
interaction vertex, classify outgoing particles as tracks or showers, and compute their
kinematics.

## Entry Point: The jsonnet Pipeline

The processing pipeline is configured in jsonnet. The relevant tail-end section in
`qlport/uboone-mabc.jsonnet` (lines ~1167–1176) reads:

```jsonnet
cm.tagger_flag_transfer("tagger"),
cm.clustering_recovering_bundle("recover_bundle", graph_name="relaxed_pid"),
cm.switch_scope(),
cm.steiner(retiler=improve_cluster_2, perf=perf),
cm.fiducialutils(),
cm.tagger_check_neutrino(...),
```

`cm` is a factory object (defined in `cfg/pgrapher/common/clus.jsonnet`) that maps
method calls to WireCell component type names, detector parameters, and configuration
objects. Each item in the list becomes one entry in the `pipeline` array consumed by
`MultiAlgBlobClustering`.

## Architectural Pattern: IEnsembleVisitor

The entire PR system is organized around a single, simple interface:

```
clus/inc/WireCellClus/IEnsembleVisitor.h
```

```cpp
class IEnsembleVisitor : public IComponent<IEnsembleVisitor> {
public:
    virtual void visit(Facade::Ensemble& ensemble) const = 0;
};
```

`MultiAlgBlobClustering` (the top-level WireCell node) holds a vector of
`IEnsembleVisitor` instances. On each input event it:

1. Deserializes the input tensor set into a `Facade::Ensemble`.
2. Calls `visitor->visit(ensemble)` for each visitor **in order**.
3. Serializes the result back to a tensor set.

This means every algorithm in the pipeline — from geometry correction to full neutrino
reconstruction — is just a class that implements `visit()`.

```
Input TensorSet (PointTrees)
        │
        ▼
MultiAlgBlobClustering
        │
        ├─► [1] ClusteringTaggerFlagTransfer   (transfer tagger flags)
        ├─► [2] ClusteringRecoveringBundle      (recover over-merged clusters)
        ├─► [3] ClusteringSwitchScope           (apply T0 correction)
        ├─► [4] CreateSteinerGraph              (build geometric skeleton)
        ├─► [5] MakeFiducialUtils               (attach fiducial geometry)
        └─► [6] TaggerCheckNeutrino             (full PR & particle ID)
                │
                ▼
        Output TensorSet (annotated clusters + PR graphs)
```

## Core Data Model

Three layers of data nesting represent the event:

```
Ensemble
 └── Grouping  ("live", "dead", ...)
      └── Cluster  (one reconstructed charge deposition)
           ├── Blob[]          (child nodes: individual wire-cell blobs)
           ├── PointCloud[]    (named arrays: "3d", "tagger_info", "isolated", ...)
           ├── Flags           (beam_flash, associated_cluster, ...)
           ├── Scope           (active coordinate frame after T0 correction)
           └── PR::Graph       (named: "relaxed_pid", "steiner", ...)
                ├── PR::Vertex[]   (interaction points)
                └── PR::Segment[]  (tracks / showers between vertices)
```

See [data_structures.md](data_structures.md) for a detailed description of each layer.

## Pipeline Stages Summary

| # | Visitor Type | Purpose |
|---|---|---|
| 1 | `ClusteringTaggerFlagTransfer` | Copy tagger metadata onto cluster flag bits |
| 2 | `ClusteringRecoveringBundle` | Split over-merged beam-flash clusters |
| 3 | `ClusteringSwitchScope` | Apply T0 space-charge correction to coordinates |
| 4 | `CreateSteinerGraph` | Build Steiner tree skeleton for each cluster |
| 5 | `MakeFiducialUtils` | Attach fiducial volume geometry queries |
| 6 | `TaggerCheckNeutrino` | Full pattern recognition: vertex + particle ID |

See [pipeline_stages.md](pipeline_stages.md) for details on each stage.

## Key Source Directories

```
clus/
├── inc/WireCellClus/   C++ headers
│   ├── IEnsembleVisitor.h          the visitor interface
│   ├── MultiAlgBlobClustering.h    top-level WireCell node
│   ├── Facade_*.h                  data model facades
│   ├── PRVertex.h / PRSegment.h    PR graph node types
│   ├── PRGraph.h                   graph construction helpers
│   ├── NeutrinoPatternBase.h       master PR algorithm class
│   ├── TaggerCheckNeutrino.h       neutrino tagger visitor
│   ├── TrackFitting.h              calorimetry / track fitting
│   └── FiducialUtils.h             fiducial geometry queries
└── src/
    ├── ClusteringTaggerFlagTransfer.cxx
    ├── ClusteringRecoveringBundle.cxx
    ├── clustering_switch_scope.cxx
    ├── CreateSteinerGraph.cxx
    ├── SteinerGrapher.cxx          Steiner tree implementation
    ├── make_fiducialutils.cxx
    ├── TaggerCheckNeutrino.cxx
    ├── NeutrinoPatternBase.cxx     main PR loop
    ├── NeutrinoVertexFinder.cxx
    ├── NeutrinoTrackShowerSep.cxx
    ├── NeutrinoKinematics.cxx
    ├── TrackFitting.cxx            ~8,350 lines
    └── Facade_Cluster.cxx / ...
```

## Deep Dive Documents

- [pipeline_stages.md](pipeline_stages.md) — Detailed walkthrough of each visitor
- [steiner_graph.md](steiner_graph.md) — How the Steiner tree skeleton is built
- [pattern_recognition.md](pattern_recognition.md) — `find_proto_vertex` and the PR loop
- [track_shower_separation.md](track_shower_separation.md) — Classifying particles
- [track_fitting.md](track_fitting.md) — Calorimetry and dE/dx
- [data_structures.md](data_structures.md) — Facade layer and PR graph types

## Detector Support

The same C++ code supports multiple detectors via jsonnet configuration:

| Detector | Config location |
|---|---|
| MicroBooNE (uBooNE) | `cfg/pgrapher/experiment/uboone/` |
| SBND | `cfg/pgrapher/experiment/sbnd/` |
| DUNE | `cfg/pgrapher/experiment/dune/` |

Detector-specific parameters (wire pitch, drift velocity, fiducial volume, dead wire
masks) are passed through the jsonnet `cm` factory and land in each component via
`IConfigurable::configure()`.
