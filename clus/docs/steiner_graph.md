# Deep Dive: Steiner Graph Construction

**Source files:**
- `clus/src/CreateSteinerGraph.cxx` — the IEnsembleVisitor wrapper
- `clus/src/SteinerGrapher.cxx` — the core algorithm (~1,088 lines)
- `clus/src/SteinerGrapher.h` — the `Steiner::Grapher` class
- `clus/src/SteinerFunctions.cxx` / `SteinerFunctions.h` — free helper functions
- `clus/inc/WireCellClus/Graphs.h` — underlying weighted graph type

---

## Purpose

After imaging, each cluster contains a set of 3D "blobs" — small charge
depositions on overlapping wire cells. The blobs are not yet connected in
any meaningful geometric structure. The Steiner graph stage builds a
**minimum spanning tree** over the cluster's point cloud, producing a
branching skeleton that directly represents the underlying particle trajectories.

This skeleton is the input to all subsequent pattern recognition logic.

---

## The Two-Phase Algorithm

### Phase 1 — Retiling (`ImproveCluster_2`)

The raw blob point cloud often has irregular density: dense near the wire
planes, sparse in the drift direction. Direct graph construction on the raw
cloud leads to poor skeletons with many short-range connections.

**Retiling** resamples the cluster into a new point cloud with more uniform
spacing. The `ImproveCluster_2` component (configured via `improve_cluster_2`
in the jsonnet) uses the `IBlobSampler` interface to:

1. Sample each blob to a target density (typically ~3 mm spacing).
2. Optionally apply geometric refinements to fill gaps near dead wire regions.
3. Store the result as a new named point cloud on the cluster.

The retiler is an `IPCTreeMutate` component and runs **in-place** on the
cluster before the graph is built.

### Phase 2 — Graph Construction (`Steiner::Grapher`)

The `Steiner::Grapher` class encapsulates the workspace for one cluster:

```cpp
// clus/src/SteinerGrapher.h
class Grapher {
    struct Config {
        IDetectorVolumes::pointer dv;
        IPCTransformSet::pointer  pcts;
        IPCTreeMutate::pointer    retile;
        bool perf{false};
    };
    // ...
    graph_type& get_graph(const std::string& flavor = "basic");
};
```

The `graph_type` is a Boost adjacency list with floating-point edge weights
(`Graphs::Weighted::graph_type`).

**Steps inside `Grapher`:**

1. **Point extraction**: Read the resampled point cloud. Each 3D point becomes
   a graph vertex. The `blob_vertex_map` records which blob each vertex came
   from.

2. **Edge construction**: For each pair of nearby points (within a configurable
   radius, typically ~6–12 cm), an edge is added with weight equal to the
   3D Euclidean distance between the two points. A k-d tree accelerates the
   neighbor search.

3. **Minimum spanning tree (MST)**: Boost's `kruskal_minimum_spanning_tree` is
   applied to the weighted graph. The MST connects all points with the minimum
   total edge weight — this approximates the Steiner tree for the biological
   topology (true Steiner tree is NP-hard; MST is a good practical proxy).

4. **Post-processing**: Very short edges (< ~1 mm) may be contracted. The
   resulting graph is stored as a named graph `"steiner"` on the cluster via
   `Facade::Cluster::Mixins::Graphs`.

---

## Graph Node and Edge Semantics

In the stored Steiner graph (a `Graphs::Weighted::graph_type`), vertices
and edges carry the following meaning:

| Element | Meaning in physics |
|---|---|
| Degree-1 vertex | Track/shower endpoint (particle stops or exits detector) |
| Degree-2 vertex | Intermediate point along a straight trajectory |
| Degree-3+ vertex | Branching point — candidate for an interaction vertex |
| Edge | Line segment of a track or shower spine |
| Edge weight | 3D Euclidean path length (cm) |

**This degree-based classification is what `find_proto_vertex()` exploits** to
identify candidate interaction vertices: it looks for high-degree nodes in the
Steiner graph as seeds for the primary vertex search.

---

## Relationship to the PR Graph

The Steiner graph is a **low-level geometric** graph (points from the resampled
cloud, stored as a `Graphs::Weighted::graph_type`). It is distinct from the
**PR graph** (`PR::Graph`) which is built later during pattern recognition.

| Property | Steiner graph | PR graph |
|---|---|---|
| Node type | Raw 3D point | `PR::Vertex` (physics object) |
| Edge type | Weighted distance edge | `PR::Segment` (track/shower) |
| Built by | `CreateSteinerGraph` | `TaggerCheckNeutrino` / `find_proto_vertex` |
| Purpose | Geometric skeleton for path finding | Physics result: vertex + particles |

The PR algorithm reads paths **through** the Steiner graph (via
`do_rough_path()`) to initialize the trajectories of `PR::Segment` objects.

---

## Key Free Functions (`SteinerFunctions.h`)

```cpp
// Find the path between two 3D points through the Steiner graph
std::vector<geo_point_t> do_rough_path(const Cluster& cluster,
                                        geo_point_t& first_point,
                                        geo_point_t& last_point);

// Extend a path to find the next branching point in a given direction
std::pair<geo_point_t, size_t> proto_extend_point(const Cluster& cluster,
                                                   geo_point_t& p,
                                                   geo_vector_t& dir,
                                                   geo_vector_t& dir_other,
                                                   bool flag_continue);

// Break a path at a kink point into two sub-paths
bool proto_break_tracks(const Cluster& cluster,
                        const geo_point_t& first_wcp,
                        geo_point_t& curr_wcp,
                        const geo_point_t& last_wcp,
                        std::list<geo_point_t>& wcps_list1,
                        std::list<geo_point_t>& wcps_list2,
                        bool flag_pass_check);
```

These functions are used extensively inside `NeutrinoPatternBase` to extract
geometric paths that seed `PR::Segment` construction.

---

## Performance

The most expensive part is the k-d tree neighbor search for edge construction,
which scales as O(N log N) in the number of points. For a typical MicroBooNE
neutrino event cluster (a few thousand points after retiling), the Steiner
graph stage takes O(10–100 ms). Setting `perf=true` in the configuration
prints per-cluster timing to stdout.
