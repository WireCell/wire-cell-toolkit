# Deep Dive: Neutrino Vertex Determination

**Source files:**
- `clus/src/NeutrinoVertexFinder.cxx` — all vertex functions (~3,534 lines)
- `clus/inc/WireCellClus/NeutrinoPatternBase.h` — `PatternAlgorithms` class declarations
- `clus/src/TaggerCheckNeutrino.cxx` — top-level orchestration

---

## Overview

After `find_proto_vertex()` builds the PR graph for each cluster and
`determine_direction()` assigns segment directions, the algorithm must
identify **which vertex is the neutrino interaction vertex**. This involves
two sequential phases:

1. **Per-cluster vertex determination** (`determine_main_vertex`): for each
   cluster independently, find the best candidate vertex.
2. **Global vertex selection** (`determine_overall_main_vertex` or its DL
   variant): across all clusters in the event, pick the single primary vertex
   and select the "main cluster" to which it belongs.

---

## Per-Cluster Vertex Determination

### `determine_main_vertex`

```cpp
void PatternAlgorithms::determine_main_vertex(
    Graph& graph, Facade::Cluster& cluster,
    VertexPtr& main_vertex,
    IndexedVertexSet& vertices_in_long_muon,
    IndexedSegmentSet& segments_in_long_muon,
    TrackFitting& track_fitter, IDetectorVolumes::pointer dv,
    const ParticleDataSet::pointer& particle_data,
    const IRecombinationModel::pointer& recomb_model);
```

**Step 1 — improve_vertex (first pass)**

```cpp
void PatternAlgorithms::improve_vertex(
    Graph& graph, Facade::Cluster& cluster, VertexPtr& main_vertex,
    ..., bool flag_search_vertex_activity, bool flag_final_vertex);
```

Called before candidates are evaluated. It:
- Runs `fit_vertex()` to geometrically refine the vertex position (least-squares
  fit to the directions of all incident segments).
- If `flag_search_vertex_activity=true`, calls `search_for_vertex_activities()`
  to look for unassociated Steiner terminal points near the vertex and, if
  found, builds new short segments from them.
- After any new segments are added, re-fits the vertex.
- Runs `fix_maps_shower_in_track_out()` to correct direction inconsistencies
  introduced by the new segments.

**Step 2 — scan for "only-showers" topology**

If the cluster has no track segments (only shower-flagged segments), the main
vertex is chosen differently: call `compare_main_vertices_all_showers()` which
selects the vertex with the best shower-start geometry (photon conversion gap,
shower opening angle, etc.).

**Step 3 — `examine_main_vertices_local`**

```cpp
void PatternAlgorithms::examine_main_vertices_local(
    Graph& graph, std::vector<VertexPtr>& candidates,
    const ParticleDataSet::pointer& particle_data,
    const IRecombinationModel::pointer& recomb_model);
```

Filters the candidate list: removes vertices that are clearly wrong (e.g.,
degree-1 vertices at the far end of a long muon with no other activity).

**Step 4 — `compare_main_vertices`** (for clusters with tracks)

```cpp
VertexPtr PatternAlgorithms::compare_main_vertices(
    Graph& graph, Facade::Cluster& cluster,
    std::vector<VertexPtr>& vertex_candidates);
```

Scores each candidate using a `map_vertex_num` score accumulator:

| Score contribution | Condition |
|---|---|
| +1 | Vertex is connected to the longest muon segment |
| +1 | Vertex is connected to a proton segment (outgoing) |
| +1 | Vertex is inside the fiducial volume |
| −1 | Vertex is outside fiducial volume |
| +0.5 | Vertex has ≥2 outgoing segments (neutrino-like topology) |
| −0.5 | Proton segment enters rather than exits the vertex |

The candidate with the highest score is returned. On a tie, the vertex closer
to the "Steiner graph center of mass" is preferred.

### `compare_main_vertices_all_showers`

Used when the cluster contains only shower segments. Evaluates each candidate
vertex for:
- Whether any attached shower has a photon conversion gap (start_connection_type=2).
- The minimum angle between shower axes (smaller = more likely a pair production).
- Total reconstructed shower energy.

Returns the vertex that best matches a photon pair (π0) or single-photon
topology.

---

## Global Vertex Selection

### `determine_overall_main_vertex`

```cpp
VertexPtr PatternAlgorithms::determine_overall_main_vertex(
    Graph& graph,
    ClusterVertexMap map_cluster_main_vertices,
    Facade::Cluster*& main_cluster,
    std::vector<Facade::Cluster*>& other_clusters,
    IndexedVertexSet& vertices_in_long_muon,
    IndexedSegmentSet& segments_in_long_muon,
    TrackFitting& track_fitter, IDetectorVolumes::pointer dv,
    const ParticleDataSet::pointer& particle_data,
    const IRecombinationModel::pointer& recomb_model,
    bool flag_dev_chain = true);
```

`map_cluster_main_vertices` maps each cluster to its best local vertex (from
`determine_main_vertex`). This function selects among these:

1. **`examine_main_vertices`**: loops over all cluster→vertex pairs. For each,
   checks whether another cluster's vertex is spatially closer to the current
   cluster's activity and swaps if so via `swap_main_cluster()`.

2. **`check_switch_main_cluster` / `check_switch_main_cluster_2`**: evaluates
   whether `main_cluster` should change to another cluster based on:
   - Which cluster contains the longest muon.
   - Which cluster's vertex is inside the fiducial volume.
   - Which cluster has more outgoing prongs.

3. **`compare_main_vertices_global`**: final scoring pass across all candidate
   vertices from all clusters. Uses `calc_conflict_maps()` to measure how many
   segment directions are inconsistent with each candidate vertex being the
   primary.

4. The selected vertex is marked `kNeutrinoVertex`:
   ```cpp
   final_vertex->set_flags(VertexFlags::kNeutrinoVertex);
   ```

---

## Deep-Learning Vertex Refinement

### `determine_overall_main_vertex_DL`

```cpp
bool PatternAlgorithms::determine_overall_main_vertex_DL(
    Graph& graph,
    ClusterVertexMap& map_cluster_main_vertices,
    Facade::Cluster*& main_cluster,
    std::vector<Facade::Cluster*>& other_clusters,
    IndexedVertexSet& vertices_in_long_muon,
    IndexedSegmentSet& segments_in_long_muon,
    TrackFitting& track_fitter, IDetectorVolumes::pointer dv,
    const ParticleDataSet::pointer& particle_data,
    const IRecombinationModel::pointer& recomb_model,
    const std::string& dl_weights,
    double dl_vtx_cut,
    double dQdx_scale, double dQdx_offset);
// Returns: true → DL selected a vertex (do not call traditional algorithm)
//          false → DL unavailable or did not improve (fall back to traditional)
```

**Step 1 — Build input point cloud**

Collects `(x, y, z, q)` for all fitted vertex positions and all interior
segment fit points in the graph:

```cpp
for each vertex in graph:
    q = vertex->fit().dQ * dQdx_scale + dQdx_offset
    append (x, y, z, q) to vec_xyzq

for each segment (interior points only):
    q = fit.dQ * dQdx_scale + dQdx_offset
    append (x, y, z, q) to vec_xyzq
```

The charge normalization (`dQdx_scale`, `dQdx_offset`) is configured in the
jsonnet (`dQdx_scale=0.1`, `dQdx_offset=-1000` for MicroBooNE) so that the
network sees charge values in its training range.

**Step 2 — SCN inference**

```cpp
// Guarded by #ifdef HAVE_PYTHON_INC
auto dnn_vtx = WCPPyUtil::SCN_Vertex("SCN_Vertex", "SCN_Vertex",
                                      dl_weights, vec_xyzq, "float32", false);
// Returns std::vector<float> of size 3: {x_pred, y_pred, z_pred} in cm
```

`SCN_Vertex` is a Python/PyTorch Sparse Convolutional Network (SCN) wrapper
from the `WCPPyUtil` package. The network was trained on simulated MicroBooNE
events to regress the 3D neutrino interaction vertex position directly from
the charge point cloud.

The entire block is compiled out if `HAVE_PYTHON_INC` is not set (i.e., if
the toolkit was built without Python support), in which case the function
always returns `false`.

**Step 3 — Nearest candidate matching**

Finds the vertex in `cand_vertices` (all graph vertices) that is closest to
the DL prediction `(x_reg, y_reg, z_reg)`:

```cpp
for (auto vtx : cand_vertices) {
    auto pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
    double dis = (pt - Point(x_reg, y_reg, z_reg)).magnitude();
    if (dis < min_dis) { min_dis = dis; min_vertex = vtx; }
}
```

**Step 4 — Acceptance cut**

```cpp
bool flag_pass = (min_dis < dl_vtx_cut);  // default dl_vtx_cut = 2 cm
```

If the DL prediction is more than `dl_vtx_cut` from any candidate vertex, the
network is considered unreliable and the function returns `false` (fall back
to traditional).

**Step 5 — Switch to DL vertex**

If `flag_pass`:
1. Update `map_cluster_main_vertices` to point to `min_vertex`.
2. If `min_vertex` belongs to a different cluster than `main_cluster`, call
   `swap_main_cluster()`.
3. Re-run `improve_vertex()` and `examine_direction()` with the new vertex.
4. Re-run long-muon cleanup (proton tagging, `find_cont_muon_segment`).
5. Return `true` to suppress the traditional algorithm.

**Fallback**: If `determine_overall_main_vertex_DL` returns `false`, the
caller (`TaggerCheckNeutrino::visit()`) calls `determine_overall_main_vertex`
with the traditional scoring logic.

---

## `vertices_in_long_muon` and `segments_in_long_muon`

These `IndexedVertexSet` / `IndexedSegmentSet` containers track the chain of
segments that form a "long muon" in the event. They are populated during
`examine_direction()` when a sequence of segments satisfies:
- Median `dQ/dx < 1.3 × MIP` (minimum ionizing particle charge density).
- Continuous chain of collinear segments via `find_cont_muon_segment()`.
- Total chain length > 40 cm.

They are passed to vertex determination because:
- Vertices **inside** a long muon chain should be de-prioritized as the
  primary neutrino vertex (they are more likely scatter kinks).
- The shower clustering code uses them to distinguish muon sub-segments
  from genuine EM shower sub-segments.
- The DL vertex function uses them during the cleanup pass after the DL
  vertex is accepted.

---

## Sequence in `TaggerCheckNeutrino::visit()`

```
for each cluster:
    find_proto_vertex()          ← builds PR graph
    separate_track_shower()      ← set shower flags
    determine_direction()        ← set dirsign + PDG
    examine_direction()          ← BFS propagation from local vertex
    determine_main_vertex()      ← choose best per-cluster vertex

if dl_weights not empty:
    success = determine_overall_main_vertex_DL(...)
if !success:
    determine_overall_main_vertex(...)  ← traditional global selection

examine_structure_final()        ← final topology cleanup
shower_clustering_with_nv()      ← group showers with known primary vertex
calculate_shower_kinematics()    ← energy reconstruction
```
