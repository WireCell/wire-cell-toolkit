# Pipeline Stages

This document covers each of the six `IEnsembleVisitor` stages that form the
tail end of the pattern recognition pipeline (configured in
`qlport/uboone-mabc.jsonnet`, lines 1167–1176).

---

## Stage 1 — ClusteringTaggerFlagTransfer

**jsonnet call:** `cm.tagger_flag_transfer("tagger")`  
**C++ type:** `ClusteringTaggerFlagTransfer`  
**Source:** `clus/src/ClusteringTaggerFlagTransfer.cxx`

### What it does

Upstream algorithms (cosmic tagger, flash tagger, etc.) encode their decisions
as numeric arrays stored in a point cloud named `"tagger_info"` on each cluster.
This visitor reads those arrays and copies them into the cluster's flag system,
where downstream visitors can query them with simple boolean checks.

### Why it runs first

All subsequent stages inspect cluster flags to decide how to process a cluster
(e.g., skip cosmic-tagged clusters, handle beam-flash clusters specially). The
flags must be populated before anything else runs.

### Key interfaces

- Implements `IEnsembleVisitor`
- Reads: local point cloud `"tagger_info"` on each cluster
- Writes: cluster flag bits

---

## Stage 2 — ClusteringRecoveringBundle

**jsonnet call:** `cm.clustering_recovering_bundle("recover_bundle", graph_name="relaxed_pid")`  
**C++ type:** `ClusteringRecoveringBundle`  
**Source:** `clus/src/ClusteringRecoveringBundle.cxx`

### What it does

The flash-matching stage can over-merge clusters: multiple disconnected charge
depositions get grouped together because they all matched the same beam flash.
This visitor undoes that merging by:

1. Looking at the `"isolated"` and `"perblob"` point cloud arrays on each
   beam-flash–tagged cluster.
2. Using the named graph (`relaxed_pid`) to find the connected-component
   structure within the cluster.
3. Splitting the cluster into separate sub-clusters, one per connected component.
4. Marking recovered sub-clusters with the `associated_cluster` flag.

### Configuration

```jsonnet
clustering_recovering_bundle(name="", graph_name="relaxed") :: {
    type: "ClusteringRecoveringBundle",
    data: dv_cfg + pcts_cfg + scope_cfg + {
        grouping: "live",
        array_name: "isolated",
        pcarray_name: "perblob",
        graph_name: graph_name,
    }
}
```

The `graph_name` parameter selects which pre-existing graph (built earlier in
the pipeline by a `connect_graph_relaxed` visitor) defines the connectivity.
For the neutrino reconstruction pass, `"relaxed_pid"` is used.

### Key interfaces

- Implements `IEnsembleVisitor`, `IConfigurable`
- Uses: `IDetectorVolumes`, `IPCTransformSet`
- Reads: cluster flag `beam_flash`, local PCs `"isolated"` and `"perblob"`
- Writes: new clusters in the `"live"` grouping with `associated_cluster` flag

---

## Stage 3 — ClusteringSwitchScope

**jsonnet call:** `cm.switch_scope()`  
**C++ type:** `ClusteringSwitchScope`  
**Source:** `clus/src/clustering_switch_scope.cxx`

### What it does

LArTPC reconstruction works in two coordinate frames:

- **Raw frame**: wire number + drift time (before space-charge or T0 correction)
- **Corrected frame**: physical (x, y, z) after applying T0 and space-charge
  corrections

This visitor applies the `T0Correction` transform to every cluster in the
`"live"` grouping, adding corrected point-cloud arrays. It then sets the
cluster's *default scope* to the corrected frame, so all downstream geometry
queries automatically use corrected coordinates.

### Why it matters

The Steiner graph and all pattern recognition code operate on 3D Euclidean
distances. Those distances are only meaningful in the corrected coordinate
frame. This stage is the transition point.

### Configuration

```jsonnet
switch_scope(name="", correction_name="T0Correction") :: {
    type: "ClusteringSwitchScope",
    data: {
        correction_name: correction_name,
    } + pcts_cfg + scope_cfg
}
```

### Key interfaces

- Implements `IEnsembleVisitor`, `IConfigurable`
- Uses: `IPCTransformSet`
- Reads: cluster raw scope
- Writes: corrected point cloud arrays; updates cluster default scope

---

## Stage 4 — CreateSteinerGraph

**jsonnet call:** `cm.steiner(retiler=improve_cluster_2, perf=perf)`  
**C++ type:** `CreateSteinerGraph`  
**Source:** `clus/src/CreateSteinerGraph.cxx`, `clus/src/SteinerGrapher.cxx`

### What it does

Builds a **Steiner tree** — the minimal connected subgraph that spans all the
blobs in a cluster. This graph becomes the geometric skeleton that pattern
recognition walks to find tracks and vertices.

The process runs in two phases:

1. **Retiling** (`ImproveCluster_2`): The raw blob set is resampled into a
   denser, more uniform point cloud. This removes holes and improves geometric
   quality before graph construction.

2. **Graph construction** (`SteinerGrapher`): Points in the resampled cloud
   become graph nodes. Edges are added between nearby points, weighted by 3D
   Euclidean distance. A minimum spanning tree (approximating the Steiner
   tree) is then computed. The result is stored as the named graph `"steiner"`
   on each cluster.

### Why a Steiner tree?

A Steiner tree is the topology that best represents branching particle
trajectories: it connects all deposited charge while minimizing total path
length. The tree structure (degree-1 nodes = endpoints, degree-3+ nodes =
vertices) directly maps to the physical event topology.

### Configuration

```jsonnet
steiner(name="", retiler={}, grouping="live", graph="steiner", perf=true) :: {
    type: "CreateSteinerGraph",
    data: {
        grouping: grouping,
        graph: graph,
        retiler: wc.tn(retiler),
        perf: perf,
    } + dv_cfg + pcts_cfg
}
```

### Key interfaces

- Implements `IEnsembleVisitor`, `IConfigurable`
- Uses: `IDetectorVolumes`, `IPCTransformSet`, `IPCTreeMutate` (retiler)
- Reads: cluster point clouds in corrected scope
- Writes: named graph `"steiner"` on each cluster

See [steiner_graph.md](steiner_graph.md) for a detailed walkthrough of the
graph construction algorithm.

---

## Stage 5 — MakeFiducialUtils

**jsonnet call:** `cm.fiducialutils()`  
**C++ type:** `MakeFiducialUtils`  
**Source:** `clus/src/make_fiducialutils.cxx`  
**Header:** `clus/inc/WireCellClus/FiducialUtils.h`

### What it does

Constructs a `FiducialUtils` object from the detector's static geometry (the
fiducial volume definition and dead wire regions) and attaches it to the
`"live"` grouping. The `FiducialUtils` object provides two key queries used
by `TaggerCheckNeutrino`:

- `is_inside_fiducial(point)` — is a 3D point within the active fiducial volume?
- `is_in_dead_region(point)` — does the point project onto a known dead wire?

### Why run it as a separate stage?

The `FiducialUtils` object requires both static data (the detector geometry,
which never changes) and dynamic data (the current live/dead groupings, which
change per event). Building it here, after clustering is complete, allows the
dynamic part to be finalized before the expensive neutrino reconstruction step.

### Configuration

```jsonnet
fiducialutils(name="", live_grouping="live", dead_grouping="dead", target_grouping="live") :: {
    type: "MakeFiducialUtils",
    data: {
        live: live_grouping,
        dead: dead_grouping,
        target: target_grouping,
    } + dv_cfg + fiducial_cfg + pcts_cfg
}
```

### Key interfaces

- Implements `IEnsembleVisitor`, `IConfigurable`
- Uses: `IDetectorVolumes`, fiducial definition config, `IPCTransformSet`
- Reads: `"live"` and `"dead"` groupings
- Writes: `FiducialUtils` pointer attached to the `"live"` grouping

---

## Stage 6 — TaggerCheckNeutrino

**jsonnet call:** `cm.tagger_check_neutrino(...)`  
**C++ type:** `TaggerCheckNeutrino`  
**Source:** `clus/src/TaggerCheckNeutrino.cxx`  
**Header:** `clus/inc/WireCellClus/TaggerCheckNeutrino.h`  
**Algorithm class:** `clus/inc/WireCellClus/NeutrinoPatternBase.h` (~29 KB)

### What it does

This is the main pattern recognition engine. For each cluster in the `"live"`
grouping it:

1. Calls `find_proto_vertex()` to walk the Steiner tree and locate the
   primary interaction vertex.
2. Calls `examine_structure()` to refine the graph topology (merge nearby
   vertices, split poorly fitted segments).
3. Calls `separate_track_shower()` to label each segment as track or shower
   using dE/dx and geometric criteria.
4. Calls `shower_clustering_with_nv()` to group shower segments into EM shower
   objects.
5. Calls `calculate_shower_kinematics()` and track kinematics to assign
   energies and momenta.
6. Optionally calls the deep-learning vertex predictor (`determine_overall_main_vertex_DL`)
   if `dl_weights` is provided, replacing the heuristic vertex with a CNN
   prediction.
7. Tags the cluster as neutrino / cosmic / other.

### Configuration parameters

| Parameter | Purpose |
|---|---|
| `trackfitting_config_file` | YAML file with `TrackFitting` detector parameters |
| `particle_dataset` | Path to particle data tables for particle identification |
| `recombination_model` | WireCell component providing charge→energy conversion |
| `perf` | If `true`, print per-step timing to stdout |
| `dl_weights` | Path to SCN deep-learning vertex `.pth` weights (optional) |
| `dQdx_scale`, `dQdx_offset` | Normalization for DL vertex network input |

### Key interfaces

- Implements `IEnsembleVisitor`, `IConfigurable`
- Mixes in: `NeedDV`, `NeedPCTS`, `NeedRecombModel`, `NeedParticleData`
- Reads: `"live"` grouping, `"steiner"` graph, `FiducialUtils`
- Writes: PR::Graph with labelled vertices and segments on each cluster

See [pattern_recognition.md](pattern_recognition.md) for a deep dive into the
`find_proto_vertex()` loop and all sub-algorithms called from this stage.
