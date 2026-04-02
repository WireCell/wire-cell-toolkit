# Particle Identification: Segment-Level Classification

This document covers how individual `PR::Segment` objects are classified as
tracks or showers, assigned particle types (PDG codes), and given direction
signs. This is purely **segment-level** PID — event-level flavor tagging
(νμ, νe, cosmic, etc.) is a separate, later step.

**Source files:**
- `clus/src/NeutrinoTrackShowerSep.cxx` — `separate_track_shower()`, `determine_direction()`, map fixers
- `clus/src/NeutrinoVertexFinder.cxx` — `examine_direction()` BFS propagation
- `clus/src/NeutrinoShowerClustering.cxx` — `shower_clustering_with_nv_in_main_cluster()`, `update_shower_maps()`
- `clus/inc/WireCellClus/PRSegmentFunctions.h` — all primitive PID functions
- `clus/src/PRSegmentFunctions.cxx` — implementations of the primitives

---

## The "is shower" predicate

Throughout the entire codebase, whether a segment is classified as a shower
is tested with this three-way OR:

```cpp
bool is_shower = seg->flags_any(SegmentFlags::kShowerTrajectory)
              || seg->flags_any(SegmentFlags::kShowerTopology)
              || (seg->has_particle_info() && abs(seg->particle_info()->pdg()) == 11);
```

This mirrors the WCP prototype's `get_flag_shower()`. All BFS traversals,
vertex analysis, and shower clustering use exactly this test.

---

## Stage 1 — Topology and Trajectory Classification (`separate_track_shower`)

```cpp
// NeutrinoTrackShowerSep.cxx
void PatternAlgorithms::separate_track_shower(Graph& graph, Facade::Cluster& cluster);
```

Iterates every segment in the graph. For each:

1. **Topology test first** (`segment_is_shower_topology`):
   - Collects `"associate_points"` cloud (built by `clustering_points()`)
   - Computes the spatial spread (RMS of transverse distances to the segment
     axis) in a rolling window along the segment.
   - If the spread exceeds a threshold (`rms_cut = 0.4 cm`) indicating a cone
     shape → sets `kShowerTopology` flag.
   - Also checks the ratio of above-MIP charge (`dQ/dx > 43000/cm`) vs. total
     length — a high ratio indicates dense EM shower activity.

2. **Trajectory test only if topology is not set** (`segment_is_shower_trajectory`):
   ```cpp
   bool segment_is_shower_trajectory(SegmentPtr seg,
       double step_size = 10*units::cm,
       double mip_dQ_dx = 50000/units::cm);
   ```
   - Steps along the segment in `step_size` chunks.
   - In each chunk computes `dQ/dx` from the fit data.
   - If the profile has a clear rising front followed by a fall (shower maximum
     shape) → sets `kShowerTrajectory` flag.
   - Uses `eval_ks_ratio()` and `do_track_comp()` to compare the observed
     `dQ/dx` profile against tabulated shower and track templates.

**Result**: each segment emerges with zero, one, or two shower flags set.
A segment with neither flag is tentatively a track.

---

## Primitive PID Functions (`PRSegmentFunctions.h`)

These free functions are the building blocks used by both `separate_track_shower`
and `determine_direction`.

### `segment_do_track_pid`

```cpp
// PRSegmentFunctions.h line 93
std::tuple<bool, int, int, double>
segment_do_track_pid(SegmentPtr segment,
                     std::vector<double>& L,
                     std::vector<double>& dQ_dx,
                     const ParticleDataSet::pointer& particle_data,
                     double compare_range = 35*units::cm,
                     double offset_length = 0,
                     bool flag_force = false,
                     double MIP_dQdx = 50000/units::cm);
// Returns: (success, flag_dir, pdg_code, particle_score)
```

- `do_track_comp()`: computes a comparison metric between the observed
  `dQ/dx(L)` profile and the expected profiles from `ParticleDataSet` for
  muon, proton, and pion hypotheses.
- `eval_ks_ratio()`: evaluates the Kolmogorov–Smirnov test ratio between
  two comparison scores to determine which hypothesis wins.
- `flag_dir`: +1 if the Bragg peak is at the end (stopping track going forward),
  −1 if at the beginning (track entering from outside).
- `pdg_code`: 13 (muon), 2212 (proton), 211 (pion), or 0 (undetermined).
- `particle_score`: confidence score; lower = more confident.

### `segment_determine_dir_track`

```cpp
void segment_determine_dir_track(SegmentPtr segment,
    int start_n, int end_n,
    const ParticleDataSet::pointer& particle_data,
    const IRecombinationModel::pointer& recomb_model,
    double MIP_dQdx = 43000/units::cm,
    bool flag_print = false);
```

Called for pure track segments. Steps:
1. Extracts `{L, dQ/dx}` vectors from `segment->fits()`.
2. Calls `segment_do_track_pid()` in both the forward and reverse directions.
3. Selects the direction (and PDG) that gives the better KS score.
4. Sets `seg->dirsign()`, `seg->dir_weak()`, `seg->particle_info()`, and
   `seg->particle_score()`.

### `segment_determine_shower_direction` (topology showers)

```cpp
bool segment_determine_shower_direction(SegmentPtr segment,
    const ParticleDataSet::pointer& particle_data,
    const IRecombinationModel::pointer& recomb_model,
    const std::string& cloud_name = "associate_points",
    double MIP_dQdx = 43000/units::cm,
    double rms_cut = 0.4*units::cm);
```

For `kShowerTopology` segments: uses the transverse spread profile to find
which end is the start of the shower (narrow end). Sets `dirsign()`. Then
assigns PDG=11 (electron) unconditionally.

### `segment_determine_shower_direction_trajectory` (trajectory showers)

```cpp
void segment_determine_shower_direction_trajectory(SegmentPtr segment,
    int start_n, int end_n,
    const ParticleDataSet::pointer& particle_data,
    const IRecombinationModel::pointer& recomb_model,
    double MIP_dQdx = 43000/units::cm,
    bool flag_print = false);
```

For `kShowerTrajectory` segments: same KS-test machinery as track PID but
applied in shower mode. The result sets `dirsign()` and PDG (usually 11).

### `segment_cal_4mom`

```cpp
WireCell::D4Vector<double> segment_cal_4mom(SegmentPtr segment,
    int pdg_code,
    const ParticleDataSet::pointer& particle_data,
    const IRecombinationModel::pointer& recomb_model,
    double MIP_dQdx = 50000/units::cm);
```

Computes the 4-momentum `(E, px, py, pz)` for a given PDG hypothesis:
- `cal_kine_dQdx()`: integrates `dE/dx` via the recombination model → kinetic energy.
- `cal_kine_range()`: range–energy lookup from `ParticleDataSet` tables.
- Best estimate: uses `kine_dQdx` for showers and EM-like segments, `kine_range`
  for stopping tracks (proton/pion).
- Direction from `segment_cal_dir_3vector()` and `dirsign`.

---

## Stage 2 — Per-Segment Direction and PDG (`determine_direction`)

```cpp
// NeutrinoTrackShowerSep.cxx
void PatternAlgorithms::determine_direction(Graph& graph, Facade::Cluster& cluster,
    const ParticleDataSet::pointer& particle_data,
    const IRecombinationModel::pointer& recomb_model);
```

Iterates every segment and dispatches based on its shower flags:

```
kShowerTrajectory  →  segment_determine_shower_direction_trajectory()
kShowerTopology    →  segment_determine_shower_direction()
                      + assign PDG=11, score=100 unconditionally
neither (track)    →  segment_determine_dir_track()
```

After the dispatch, every segment has:
- `dirsign()` set to +1, −1, or 0
- `particle_info()` set with PDG code and 4-momentum
- `particle_score()` set

The debug log line (always emitted at DEBUG level) shows the resulting
classification per segment:
```
determine_direction: Track nfits=42 nwcpts=42 len=18.3cm dirsign=1 dir_weak=0 start_n=3 end_n=1 pdg=13
```

### Special reclassification: `examine_good_tracks`

```cpp
// NeutrinoTrackShowerSep.cxx
void PatternAlgorithms::examine_good_tracks(Graph& graph, Facade::Cluster& cluster,
    const ParticleDataSet::pointer& particle_data);
```

A post-pass over all **non-shower track segments** with a confirmed direction.
Reclassifies a segment as PDG=11 (electron) if:
- It has ≥4 daughter showers downstream, OR (≥2 daughter showers with total
  shower length > 50 cm).
- AND its endpoint has a near-180° angle back-scatter OR is drift-parallel.
- AND the segment itself is short (< 15 cm).

This catches short "stub" segments that are actually the start of an EM shower
but were not flagged by the topology/trajectory tests.

---

## Stage 3 — Direction Consistency Fixers

After per-segment direction is set, several passes enforce physical consistency
at vertices (a particle cannot simultaneously enter and exit from both sides).

### `fix_maps_multiple_tracks_in`

If a vertex has multiple **incoming track** segments (`dirsign` pointing in),
that is unphysical — at most one particle can enter a vertex. Reset `dirsign=0`
and mark `dir_weak=true` for all incoming tracks at such vertices.

### `fix_maps_shower_in_track_out`

If a vertex has incoming showers **and** a strong outgoing non-shower track,
the shower directions are flipped (`dirsign *= -1`, `dir_weak=true`). The
rationale: a hard outgoing track most likely defines the true "out" direction.

### `improve_maps_one_in`

Iteratively: for each vertex that has exactly one incoming segment (regardless
of shower/track type), set all other connected segments' directions outward
relative to that vertex. Repeat until stable.

### `improve_maps_shower_in_track_out` / `improve_maps_no_dir_tracks` / `improve_maps_multiple_tracks_in`

Additional convergence passes for topologies not handled by the above:
- Showers with undetermined direction adjacent to known-direction tracks.
- Segments with `dirsign=0` that can be resolved by their neighbors.

### `examine_maps`

Final sanity check: verifies no vertex still has multiple physical incoming
tracks after the fix passes. Returns `true` if any corrections were made
(triggers another round of fix passes in the caller).

---

## Stage 4 — BFS Direction Propagation (`examine_direction`)

```cpp
// NeutrinoVertexFinder.cxx
bool PatternAlgorithms::examine_direction(Graph& graph,
    VertexPtr vertex, VertexPtr main_vertex,
    IndexedVertexSet& vertices_in_long_muon,
    IndexedSegmentSet& segments_in_long_muon,
    const ParticleDataSet::pointer& particle_data,
    const IRecombinationModel::pointer& recomb_model,
    bool flag_final);
```

Called once per cluster from the vertex determination stage, starting from
the chosen primary vertex. It **re-propagates** directions outward from the
vertex using BFS, which can correct directions that were set locally in
Stage 2 without knowing the vertex.

### BFS mechanics

- Queue: `std::vector<pair<VertexPtr, SegmentPtr>>` — current vertex + the
  segment to cross.
- `used_vertices`, `used_segments` prevent revisiting.
- **Shower segments are boundaries**: BFS does not descend into shower sub-trees.
  This keeps shower sub-topologies independent of the track topology.

### Direction rule

For each queued segment, the BFS "previous vertex" is known. The segment's
start/end are compared to `prev_vtx.wcpt.point` by distance:

```cpp
bool flag_start = (dist(wcpts.front(), prev_vtx) < dist(wcpts.back(), prev_vtx));
seg->dirsign(flag_start ? 1 : -1);
```

### PDG re-assignment rules (in priority order)

| Condition | Result |
|---|---|
| `flag_shower_in && dirsign==0 && !is_shower` | PDG=11 |
| `flag_shower_in && length < 2 cm && !is_shower` | PDG=11 |
| `flag_shower_in && (abs(pdg)==13 OR pdg==0)` | PDG=11 |
| `num_daughter_showers>=4` OR `(≥2 and total>50cm)` + angle check | PDG=11 |
| `pdg==11 && num_daughter_showers<=2 && length>10cm && !flag_only_showers && score<=100 && straight` | Reclassify back to track (PDG=0) |
| Long sustained `dQ/dx > threshold` + `length > 40 cm` | PDG=13 → add to `segments_in_long_muon` |

`flag_shower_in` is true when a **previously-processed** sibling segment at
the same vertex is already classified as a shower pointing inward — this
prevents order-dependent results by only checking already-visited segments.

`flag_only_showers` is a cluster-level flag (set during topology analysis at
the top of `examine_direction`) indicating the cluster has no clear track
segments, which relaxes the electron-to-track reclassification guard.

`flag_final`: when `true`, more aggressively resets `dirsign` on segments
that had a strong direction from Stage 2 but are now inconsistent with the
chosen vertex. Used in the cleanup pass after the primary vertex is finalized.

---

## Stage 5 — Shower Grouping (`shower_clustering_with_nv_in_main_cluster`)

```cpp
// NeutrinoShowerClustering.cxx
void PatternAlgorithms::shower_clustering_with_nv_in_main_cluster(
    Graph& graph, VertexPtr main_vertex,
    IndexedShowerSet& showers,
    ShowerVertexMap& map_vertex_in_shower,
    ShowerSegmentMap& map_segment_in_shower,
    VertexShowerSetMap& map_vertex_to_shower,
    ClusterPtrSet& used_shower_clusters,
    IndexedVertexSet& vertices_in_long_muon,
    IndexedSegmentSet& segments_in_long_muon,
    const ParticleDataSet::pointer& particle_data,
    const IRecombinationModel::pointer& recomb_model);
```

Once all segments have PDG codes and directions, this function groups them
into `PR::Shower` objects. Starting from `main_vertex`:

1. **BFS over the track skeleton** (does not cross shower segments).
2. At each step, if the next segment is a shower (or is in `segments_in_long_muon`),
   a new `PR::Shower` is created with:
   - `set_start_vertex(parent_vtx, connection_type=1)`: the track vertex where
     the shower branches off.
   - `set_start_segment(curr_sg)`: the first shower segment.
3. `complete_structure_with_start_segment(used_segments)`: BFS inside the shower
   sub-tree to collect all shower segments into the `PR::Shower` object.
4. `update_particle_type(particle_data, recomb_model)`: majority vote over all
   segments in the shower to decide the final PDG (electron vs. photon).

### Long-muon re-classification

After initial grouping, any shower object whose start segment has PDG=13 is
inspected. If it contains more non-muon segments (by count and length) than
muon segments, and the longest muon sub-segment is < 60 cm, the entire shower
is reclassified as an EM shower (PDG=11). The muon segments are removed from
`segments_in_long_muon`.

### Map maintenance (`update_shower_maps`)

After any shower modification, the three maps are rebuilt:
- `map_vertex_in_shower`: non-start vertices → their shower (used to detect
  vertices already claimed by a shower)
- `map_segment_in_shower`: segment → shower
- `map_vertex_to_shower`: start vertex → set of showers at that vertex

These maps are the "ownership registry" that prevents segments from being
claimed by multiple showers.

---

## Summary: Classification Flow

```
separate_track_shower()
    ├─ segment_is_shower_topology()    → set kShowerTopology flag
    └─ segment_is_shower_trajectory()  → set kShowerTrajectory flag

determine_direction()
    ├─ kShowerTrajectory → segment_determine_shower_direction_trajectory()
    ├─ kShowerTopology   → segment_determine_shower_direction() + PDG=11
    └─ track             → segment_determine_dir_track()
        └─ segment_do_track_pid() → PDG ∈ {13, 2212, 211, 0}

examine_good_tracks()  ← reclassify short track stubs with many daughter showers

fix_maps_*()           ← enforce direction consistency at vertices
improve_maps_*()       ← propagate directions from known segments

examine_direction()    ← BFS re-propagation from primary vertex
    └─ re-assigns dirsign and PDG based on context

shower_clustering_with_nv_in_main_cluster()
    └─ groups is_shower segments into PR::Shower objects
```
