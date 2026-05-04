# Deep Dive: Pattern Recognition — find_proto_vertex and the PR Loop

**Source files:**
- `clus/inc/WireCellClus/NeutrinoPatternBase.h` — `PatternAlgorithms` class (~29 KB)
- `clus/src/NeutrinoPatternBase.cxx` — core implementations
- `clus/src/NeutrinoVertexFinder.cxx` — vertex candidate finding
- `clus/src/NeutrinoStructureExaminer.cxx` — graph topology refinement
- `clus/src/NeutrinoOtherSegments.cxx` — secondary particle discovery
- `clus/src/TaggerCheckNeutrino.cxx` — the `visit()` entry point

---

## Overview

`TaggerCheckNeutrino::visit()` drives the full reconstruction loop. For each
cluster, it delegates to the `PatternAlgorithms` class (a mixin base,
`NeutrinoPatternBase.h`) which owns every sub-algorithm. The class is
stateless per-event — all working state is passed through function parameters
(`Graph&`, `Cluster&`, `TrackFitting&`, `IDetectorVolumes::pointer`).

The master function is:

```cpp
// NeutrinoPatternBase.h, line 116
bool find_proto_vertex(Graph& graph, Facade::Cluster& cluster,
                       TrackFitting& track_fitter,
                       IDetectorVolumes::pointer dv,
                       bool flag_break_track = true,
                       int nrounds_find_other_tracks = 2,
                       bool flag_back_search = true);
```

---

## The PR Loop: Step by Step

### Step 1 — Initialize the first segment (`init_first_segment`)

```cpp
SegmentPtr init_first_segment(Graph& graph, Facade::Cluster& cluster,
                               Facade::Cluster* main_cluster,
                               TrackFitting& track_fitter,
                               IDetectorVolumes::pointer dv,
                               bool flag_back_search = true);
```

1. Find the cluster's extremal points (farthest-apart pair of points in the
   Steiner graph) by doing a two-pass graph traversal (BFS from one endpoint,
   then BFS from the farthest point found).
2. Call `do_rough_path()` to extract the Steiner path between those endpoints.
3. Fit the path with `TrackFitting` to produce a smooth `PR::Segment`.
4. Add the segment's two endpoint vertices to the `PR::Graph`.

This gives the algorithm a "spine" segment to start from.

### Step 2 — Find candidate vertices (`find_cluster_vertices`)

```cpp
std::vector<VertexPtr> find_cluster_vertices(Graph& graph,
                                              const Facade::Cluster& cluster);
```

Scans the Steiner graph for **high-degree nodes** (degree ≥ 3). Each such
node is a topological branching point and a candidate for a physical
interaction vertex. Returns a sorted list of candidates (sorted by estimated
number of outgoing tracks/showers).

### Step 3 — Break the initial segment at candidate vertices

For each high-degree Steiner node that lies near the initial segment,
`proto_break_tracks()` splits the segment at the kink point into two
sub-segments. This refines the coarse initial spine into a more accurate
representation of the actual trajectory.

Internally this calls:

```cpp
bool proto_break_tracks(const Cluster& cluster,
                         const geo_point_t& first_wcp,
                         geo_point_t& curr_wcp,
                         const geo_point_t& last_wcp,
                         std::list<geo_point_t>& wcps_list1,
                         std::list<geo_point_t>& wcps_list2,
                         bool flag_pass_check);
```

### Step 4 — Examine and refine graph structure (`examine_structure`)

```cpp
void examine_structure(Graph& graph, Facade::Cluster& cluster,
                       TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
```

A composite call that runs four sub-passes:

| Sub-function | What it does |
|---|---|
| `examine_structure_1` | Merge near-collinear segments (remove spurious vertices) |
| `examine_structure_2` | Split curved segments at significant kink angles |
| `examine_structure_3` | Handle isochronous (wire-parallel) segments specially |
| `examine_structure_4` | Resolve topology at specific candidate vertices |

After this step the PR graph has a stable topology.

### Step 5 — Find other segments (`find_other_segments`)

```cpp
void find_other_segments(Graph& graph, Facade::Cluster& cluster,
                          TrackFitting& track_fitter,
                          IDetectorVolumes::pointer dv,
                          bool flag_break_track = true,
                          double search_range = 1.5*units::cm,
                          double scaling_2d = 0.8);
```

The initial segment and its splits cover the longest trajectory, but a
neutrino event has secondary particles. This function:

1. Inspects every vertex in the current graph.
2. Searches the **unused** region of the Steiner graph near each vertex
   (within `search_range`).
3. If a significant unused cluster of points is found, builds a new
   `PR::Segment` from the vertex into that region using `do_rough_path()`.
4. Calls `TrackFitting` to fit the new segment.

This is iterated `nrounds_find_other_tracks` times (default 2) to
progressively recover all secondary tracks.

### Step 6 — Examine segments (`examine_segment`)

```cpp
void examine_segment(Graph& graph, Facade::Cluster& cluster,
                     TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
```

Refines each segment individually:

- `crawl_segment()`: Walks a segment point-by-point looking for significant
  kinks. Splits the segment at any kink exceeding an angular threshold.
- `examine_partial_identical_segments()`: Removes duplicate segments that
  share more than ~80% of their points with another segment.

### Step 7 — Examine vertices (`examine_vertices`)

```cpp
void examine_vertices(Graph& graph, Facade::Cluster& cluster,
                      TrackFitting& track_fitter, IDetectorVolumes::pointer dv,
                      VertexPtr main_vertex = nullptr);
```

Sub-passes:

| Sub-function | What it does |
|---|---|
| `examine_vertices_1` | Merge vertices closer than 0.1 cm |
| `examine_vertices_2` | Remove degree-2 vertices (replace with single segment) |
| `examine_vertices_4` | Handle special topologies (T-junctions, etc.) |
| `examine_vertices_3` | Final cleanup with main-cluster context |

### Step 8 — Break segments at kinks (`break_segments`)

```cpp
bool break_segments(Graph& graph, TrackFitting& track_fitter,
                    IDetectorVolumes::pointer dv,
                    std::vector<SegmentPtr>& remaining_segments,
                    float dis_cut = 0);
```

A final pass that uses the track fitting residuals to identify segments that
are actually two tracks (or a track + shower) bent by a secondary scattering.
Any segment with a large fit residual at a kink point is split.

### Step 9 — Merge nearby vertices (`merge_nearby_vertices`)

```cpp
bool merge_nearby_vertices(Graph& graph, Facade::Cluster& cluster,
                            TrackFitting& track_fitter,
                            IDetectorVolumes::pointer dv);
```

After all splitting is done, vertices that ended up closer than 0.1 cm to each
other are merged into a single vertex. The incident segments are re-fitted to
the merged position.

---

## Determining the Primary Vertex

Once the PR graph is topologically stable, the algorithm must identify which
vertex is the **neutrino interaction vertex**. This is done by
`determine_overall_main_vertex()`:

```cpp
VertexPtr determine_overall_main_vertex(
    Graph& graph,
    ClusterVertexMap map_cluster_main_vertices,
    Facade::Cluster* main_cluster,
    std::vector<Facade::Cluster*>& other_clusters,
    IndexedVertexSet& vertices_in_long_muon,
    IndexedSegmentSet& segments_in_long_muon,
    TrackFitting& track_fitter,
    IDetectorVolumes::pointer dv,
    const ParticleDataSet::pointer& particle_data,
    const IRecombinationModel::pointer& recomb_model,
    bool flag_dev_chain = true);
```

### Scoring logic

For each candidate vertex `v`, a score is computed based on:

1. **Number of outgoing segments** — neutrino vertices typically have 2+
   outgoing prongs.
2. **Fiducial position** — vertices inside the fiducial volume score higher.
3. **dE/dx near the vertex** — a Bragg peak (high dE/dx) at one end of a
   segment implies that end is the stopping end, not the start.
4. **Shower presence** — a vertex with an attached shower is more likely a
   neutrino vertex than a cosmic endpoint.
5. **Conflict map** (`calc_conflict_maps`) — a metric of how well the local
   topology is explained by this vertex being the primary.

The vertex with the highest score is selected and marked
`VertexFlags::kNeutrinoVertex`.

### Deep-learning vertex refinement (optional)

If `dl_weights` is configured, an alternative path is taken:

```cpp
bool determine_overall_main_vertex_DL(
    Graph& graph, ClusterVertexMap& map_cluster_main_vertices,
    Facade::Cluster*& main_cluster, ...,
    const std::string& dl_weights,
    double dl_vtx_cut, double dQdx_scale, double dQdx_offset);
```

This function:
1. Assembles a sparse 3D charge image around each candidate vertex.
2. Normalizes charge values using `dQdx_scale` and `dQdx_offset`.
3. Runs inference on a SCN (Sparse Convolutional Network) loaded from
   `dl_weights`.
4. Returns the vertex the network ranks highest, provided it is within
   `dl_vtx_cut` (default 2 cm) of a heuristic candidate.
5. Returns `false` if the network fails or does not improve on the heuristic
   (triggering fallback to `determine_overall_main_vertex`).

---

## Final Structure Pass (`examine_structure_final`)

After the primary vertex is known, a second round of topology cleanup is
performed with the main vertex as an anchor:

```cpp
bool examine_structure_final(Graph& graph, VertexPtr main_vertex,
                              Facade::Cluster& cluster,
                              TrackFitting& track_fitter,
                              IDetectorVolumes::pointer dv);
```

This handles cases where the topology changes meaning once the primary vertex
is fixed (e.g., a segment that was a "track entering from outside" now looks
like a "track exiting from the vertex").

---

## Output

At the end of `find_proto_vertex()`, the `PR::Graph` on the cluster contains:

- One or more `PR::Vertex` objects, exactly one of which is marked
  `kNeutrinoVertex`.
- Multiple `PR::Segment` objects connecting the vertices, each with:
  - A vector of `WCPoint` waypoints (the fitted trajectory).
  - A `dirsign` (+1/-1/0) indicating whether the segment points away from
    or towards the primary vertex.
  - Segment flags: `kShowerTrajectory`, `kShowerTopology`, `kAvoidMuonCheck`.

This graph is the input to the track/shower separation and kinematics stages.
See [track_shower_separation.md](track_shower_separation.md).
