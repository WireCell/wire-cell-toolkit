# Deep Dive: Shower Clustering

**Source files:**
- `clus/src/NeutrinoShowerClustering.cxx` — all shower clustering functions (~3,310 lines)
- `clus/inc/WireCellClus/PRShower.h` — `PR::Shower`, `ShowerData`, comparators
- `clus/src/PRShower.cxx` — `Shower` implementation
- `clus/inc/WireCellClus/PRShowerFunctions.h` — free functions on showers

---

## Purpose

Shower clustering takes a PR graph whose segments already have PDG codes and
direction signs (set by `separate_track_shower` and `determine_direction`) and
groups **shower-flagged segments** into coherent `PR::Shower` objects.  A
`PR::Shower` is the physics-level object for an EM shower: it knows its start
vertex, its start segment, the collection of all sub-segments and sub-vertices
it contains, and its reconstructed kinematics.

---

## Data Structures

### `PR::Shower` (`PRShower.h`)

```cpp
class Shower : public TrajectoryView      // graph view of the shower sub-graph
             , public Flagged<ShowerFlags>// kShower, kKinematics
             , public HasDPCs<Shower>     // associated DynamicPointClouds
```

The `TrajectoryView` base holds a filtered view of the main `PR::Graph`
restricted to the segments and vertices belonging to this shower.

```cpp
struct ShowerData {
    int    particle_type;          // PDG: 11=electron, 22=photon, 13=muon, …
    double kenergy_range;          // kinetic energy from range–energy tables
    double kenergy_dQdx;           // kinetic energy from dE/dx integration
    double kenergy_charge;         // kinetic energy from total charge sum
    double kenergy_best;           // best estimate (chosen from the three above)

    WireCell::Point  start_point;  // shower start position
    WireCell::Point  end_point;    // shower end / furthest point
    WireCell::Vector init_dir;     // initial direction at start_point

    int start_connection_type;     // 1=direct, 2=gap, 3=association, 4=excluded
};
```

Key accessors:
```cpp
auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
SegmentPtr start_seg = shower->start_segment();
void shower->set_start_vertex(VertexPtr v, int connection_type);
void shower->set_start_segment(SegmentPtr s);
void shower->complete_structure_with_start_segment(std::set<SegmentPtr>& used);
void shower->update_particle_type(ParticleDataSet::pointer, IRecombinationModel::pointer);
void shower->calculate_kinematics(TrackFitting&, IDetectorVolumes::pointer, ...);
```

### Determinism infrastructure

Because `PR::Vertex` and `PR::Segment` are owned by `shared_ptr`, a raw
`std::set<VertexPtr>` would order by pointer address — non-deterministic
across runs. The shower clustering code uses index-stable ordered sets:

```cpp
// PRShower.h
struct VertexIndexCmp  { bool operator()(const VertexPtr& a, const VertexPtr& b) const; };
struct SegmentIndexCmp { bool operator()(const SegmentPtr& a, const SegmentPtr& b) const; };
using IndexedVertexSet  = std::set<VertexPtr,  VertexIndexCmp>;
using IndexedSegmentSet = std::set<SegmentPtr, SegmentIndexCmp>;
using IndexedShowerSet  = std::set<ShowerPtr,  ShowerIndexCmp>;
```

The comparators order by `get_graph_index()` — a monotonically assigned
integer that is stable within a single event.

### The four maps

```cpp
// Defined in NeutrinoPatternBase.h
using ShowerVertexMap    = std::map<VertexPtr,  ShowerPtr>;   // non-start vtx → shower
using ShowerSegmentMap   = std::map<SegmentPtr, ShowerPtr>;   // segment → shower
using VertexShowerSetMap = std::map<VertexPtr, IndexedShowerSet>; // start vtx → showers
using ClusterPtrSet      = std::set<Facade::Cluster*>;        // clusters with ≥1 shower seg
```

These maps are the "ownership registry": before adding a segment or vertex to
a shower, the algorithm checks whether it is already in `map_segment_in_shower`
or `map_vertex_in_shower`.

---

## Main Clustering Function

### `shower_clustering_with_nv_in_main_cluster`

```cpp
void PatternAlgorithms::shower_clustering_with_nv_in_main_cluster(
    Graph& graph, VertexPtr main_vertex,
    IndexedShowerSet& showers, /* the four maps */,
    IndexedVertexSet& vertices_in_long_muon,
    IndexedSegmentSet& segments_in_long_muon,
    const ParticleDataSet::pointer& particle_data,
    const IRecombinationModel::pointer& recomb_model);
```

**Algorithm:**

1. **BFS over the track skeleton** starting from `main_vertex`. Uses
   `sorted_out_edges()` for deterministic traversal order. The BFS queue
   holds `(SegmentPtr, VertexPtr)` pairs — the segment to cross and the
   vertex we arrive at.

2. **Shower detection test** at each BFS step:
   ```cpp
   bool is_shower_seg = seg->flags_any(kShowerTrajectory)
                     || seg->flags_any(kShowerTopology)
                     || (seg->has_particle_info() && abs(pdg) == 11);
   bool in_long_muon  = segments_in_long_muon.count(seg) > 0;
   ```

3. **On shower detection**: create a `PR::Shower`, record its start vertex
   (the BFS "parent vertex") and start segment. BFS does **not** descend
   into the shower — the shower's internal structure is built separately by
   `complete_structure_with_start_segment()`.

4. **On track detection**: continue BFS from `daughter_vtx`, expanding the
   track skeleton.

5. After BFS: for each new shower call:
   - `complete_structure_with_start_segment(used_segments)`: internal BFS
     inside the shower sub-tree.
   - `update_particle_type()`: majority-vote PDG over all sub-segments.
   - Force PDG=11 on the start segment if its PDG is 0 (catches single-segment
     showers that were not updated by `update_particle_type`).

**Long muon re-classification**: After all showers are grouped, any shower
whose `particle_type == 13` (started as a long muon) is inspected. If it has
more non-muon sub-segments (by count and length) than muon sub-segments AND
the longest muon segment is < 60 cm, the shower is reclassified as EM (PDG=11)
and the muon sub-segments are removed from `segments_in_long_muon`.

---

## Connecting Showers to the Main Vertex

### `shower_clustering_connecting_to_main_vertex`

Handles showers that were found in a satellite cluster but whose geometry
places them within a few cm of the main vertex. The function:
1. Collects all showers starting at `main_vertex`.
2. For each shower segment endpoint, checks nearby clusters (not yet in
   `used_shower_clusters`).
3. If a nearby cluster is within `search_range` and its closest point to the
   shower connects naturally, links it as a sub-cluster of the shower with
   `start_connection_type = 2` (gap connection).

### `shower_clustering_with_nv_from_vertices`

For satellite clusters not covered by the above: for each cluster in
`other_clusters`, runs `find_proto_vertex()` independently, then calls
`shower_clustering_with_nv_in_main_cluster()` on each satellite cluster
using its local main vertex.  Results are merged into the global shower maps
via `update_shower_maps()`.

### `shower_clustering_in_other_clusters`

Handles clusters that are not connected to any shower yet. A geometric
proximity search links these to the nearest shower endpoint, assigning
`start_connection_type = 3` (association).

---

## Shower Examination and Merging

After initial clustering, several passes examine and potentially merge showers:

### `examine_shower_1` / `examine_showers`

Check topology consistency:
- Showers with `start_connection_type = 4` are excluded from final output.
- Very short showers (< 2 cm) that are adjacent to a long muon track are
  reclassified as a muon segment.
- Two showers starting at the same vertex with similar directions and compatible
  invariant masses may be merged.

### `examine_merge_showers`

Specifically for EM showers: if two showers have an opening angle and energy
consistent with a π0 decay, they are kept separate but linked for kinematic
reconstruction. If they are too similar in direction (conversion pair), they
may be merged.

---

## Kinematics Calculation

After shower structure is finalized:

```cpp
void PatternAlgorithms::calculate_shower_kinematics(
    IndexedShowerSet& showers,
    IndexedVertexSet& vertices_in_long_muon,
    IndexedSegmentSet& segments_in_long_muon,
    Graph& graph, TrackFitting& track_fitter,
    IDetectorVolumes::pointer dv,
    const ParticleDataSet::pointer& particle_data,
    const IRecombinationModel::pointer& recomb_model);
```

For each shower, three energy estimates are computed and the best is chosen:

| Estimate | Method |
|---|---|
| `kenergy_dQdx` | Integrate dQ/dx × (1/W_ion/recombination) along all shower segments |
| `kenergy_charge` | `cal_kine_charge()` using pre-collected 2D charge maps |
| `kenergy_range` | Range–energy table lookup (only meaningful for short, well-contained showers) |
| `kenergy_best` | Weighted combination; `kenergy_charge` preferred for showers |

The shower's `init_dir` is set from `segment_cal_dir_3vector()` at the start
segment, pointing away from `start_vertex`.

---

## `update_shower_maps`

```cpp
void PatternAlgorithms::update_shower_maps(
    IndexedShowerSet& showers,
    ShowerVertexMap& map_vertex_in_shower,
    ShowerSegmentMap& map_segment_in_shower,
    VertexShowerSetMap& map_vertex_to_shower,
    ClusterPtrSet& used_shower_clusters);
```

Rebuilds all four maps from scratch by iterating over all showers and their
`TrajectoryView` contents. Called after any shower modification (creation,
merging, reclassification) to keep the maps consistent.

The distinction between `map_vertex_to_shower` (start vertex → shower) and
`map_vertex_in_shower` (internal vertex → shower) is important: a vertex in
`map_vertex_to_shower` can be the start of multiple showers (e.g., two photons
from π0 decay), whereas a vertex in `map_vertex_in_shower` belongs to exactly
one shower.
