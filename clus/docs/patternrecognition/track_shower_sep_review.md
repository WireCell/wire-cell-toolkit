# Track/Shower Separation — Port Review

**Date:** 2026-04-11  
**Reviewer:** Claude (automated review + manual verification)  
**Prototype refs:** `NeutrinoID_track_shower.h`, `NeutrinoID_shower_clustering.h`  
**Toolkit refs:** `clus/src/NeutrinoTrackShowerSep.cxx`, `clus/src/NeutrinoPatternBase.cxx`  
**Scope:** Functional equivalence, bugs, determinism, efficiency, multi-APA correctness.

---

## 1. Function Map

| Prototype name | Toolkit name | Toolkit location | Status |
|---|---|---|---|
| `NeutrinoID::separate_track_shower(PR3DCluster*)` | `PatternAlgorithms::separate_track_shower` | `NeutrinoTrackShowerSep.cxx:37` | ✅ Equivalent |
| `NeutrinoID::determine_direction(PR3DCluster*)` | `PatternAlgorithms::determine_direction` | `NeutrinoTrackShowerSep.cxx:68` | ✅ Equivalent (improved) |
| `NeutrinoID::calculate_num_daughter_showers` | `PatternAlgorithms::calculate_num_daughter_showers` | `NeutrinoTrackShowerSep.cxx:164` | ✅ Equivalent |
| `NeutrinoID::examine_good_tracks(int)` | `PatternAlgorithms::examine_good_tracks` | `NeutrinoTrackShowerSep.cxx:327` | ✅ Equivalent (improved) |
| `NeutrinoID::fix_maps_multiple_tracks_in(int)` | `PatternAlgorithms::fix_maps_multiple_tracks_in` | `NeutrinoTrackShowerSep.cxx:440` | ✅ Equivalent |
| `NeutrinoID::fix_maps_shower_in_track_out(int)` | `PatternAlgorithms::fix_maps_shower_in_track_out` | `NeutrinoTrackShowerSep.cxx:505` | ✅ Equivalent |
| `NeutrinoID::improve_maps_one_in(PR3DCluster*, bool)` | `PatternAlgorithms::improve_maps_one_in` | `NeutrinoTrackShowerSep.cxx:573` | ✅ Equivalent (determinism fix) |
| `NeutrinoID::improve_maps_shower_in_track_out(int, bool)` | `PatternAlgorithms::improve_maps_shower_in_track_out` | `NeutrinoTrackShowerSep.cxx:680` | ✅ Equivalent (determinism fix) |
| `NeutrinoID::improve_maps_no_dir_tracks(int)` | `PatternAlgorithms::improve_maps_no_dir_tracks` | `NeutrinoTrackShowerSep.cxx:824` | ✅ Equivalent |
| `NeutrinoID::improve_maps_multiple_tracks_in(int)` | `PatternAlgorithms::improve_maps_multiple_tracks_in` | `NeutrinoTrackShowerSep.cxx:1298` | ✅ Equivalent (bug fixed — see §2.1) |
| `NeutrinoID::judge_no_dir_tracks_close_to_showers(int)` | `PatternAlgorithms::judge_no_dir_tracks_close_to_showers` | `NeutrinoTrackShowerSep.cxx:1385` | ✅ Equivalent (multi-APA enhanced) |
| `NeutrinoID::examine_maps(int)` | `PatternAlgorithms::examine_maps` | `NeutrinoTrackShowerSep.cxx:1463` | ✅ Equivalent |
| `NeutrinoID::examine_all_showers(PR3DCluster*)` | `PatternAlgorithms::examine_all_showers` | `NeutrinoTrackShowerSep.cxx:1542` | ✅ Equivalent (improved) |
| `NeutrinoID::shower_determing_in_main_cluster(PR3DCluster*)` | `PatternAlgorithms::shower_determining_in_main_cluster` | `NeutrinoTrackShowerSep.cxx:1869` | ✅ Equivalent |
| `NeutrinoID::print_segs_info(int, ProtoVertex*)` | `PatternAlgorithms::print_segs_info` | `NeutrinoPatternBase.cxx:1718` | ✅ Equivalent |

---

## 2. Bugs

### 2.1 `improve_maps_multiple_tracks_in`: incomplete `is_shower` check (fixed)

**Status:** Fixed in this review.  
**File:** `NeutrinoTrackShowerSep.cxx`, inside `improve_maps_multiple_tracks_in`.

**Before (buggy):**
```cpp
bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                 sg->flags_any(SegmentFlags::kShowerTopology);
```

**After (fixed):**
```cpp
// matches prototype get_flag_shower() = kShowerTrajectory || kShowerTopology || particle_type==11
bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                 sg->flags_any(SegmentFlags::kShowerTopology) ||
                 (sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 11);
```

**Prototype code (reference):** `NeutrinoID_track_shower.h:861`
```cpp
if (sg->get_flag_shower()) n_in_shower++;
```
where `get_flag_shower()` = `flag_shower_trajectory || flag_shower_topology || (particle_type==11)`.

**Impact:** Without the fix, segments reclassified as electrons by dQ/dx (e.g., by earlier passes of `improve_maps_shower_in_track_out`) would not be recognised as showers inside `improve_maps_multiple_tracks_in`. They would appear as incoming non-shower tracks and be pushed onto `in_tracks`, potentially triggering incorrect further electron promotion cascades at multi-track vertices.

**Why it is isolated to this function:** All 13 other locations that compute `is_shower` in this file use the full three-part check. The two-condition version was a copy-paste omission.

---

## 3. Determinism Improvements (toolkit over prototype)

The prototype iterates `map_vertex_segments` (a `std::map<ProtoVertex*, ...>`) and `std::set<ProtoSegment*>`, both ordered by raw pointer address. The toolkit replaces these with stable, graph-index-ordered traversals.

| Function | Prototype iteration | Toolkit iteration |
|---|---|---|
| `improve_maps_one_in` | `map_vertex_segments` (pointer order) | `ordered_nodes(graph)` |
| `improve_maps_shower_in_track_out` | `map_vertex_segments` (pointer order) | `ordered_nodes(graph)` |
| `improve_maps_multiple_tracks_in` | `map_vertex_segments` (pointer order) | `ordered_nodes(graph)` |
| `judge_no_dir_tracks_close_to_showers` | `std::set<ProtoSegment*>` (pointer order) | `ordered_edges(graph)` + `std::vector` |
| `improve_maps_no_dir_tracks` | `map_segment_vertices` (pointer order) | `boost::edges(graph)` (stable index) |
| `separate_track_shower` | `map_segment_vertices` (pointer order) | `boost::edges(graph)` (stable index) |
| `determine_direction` | `map_segment_vertices` (pointer order) | `boost::edges(graph)` (stable index) |

All iterative update loops (`flag_update = true` loops) in `improve_maps_*` visit vertices in the same deterministic order each outer iteration, removing the run-to-run variation that was present in the prototype.

---

## 4. Multi-APA / Multi-face Correctness

### 4.1 `judge_no_dir_tracks_close_to_showers`: explicit APA/face resolution

The prototype's inner loop calls:
```cpp
(*it1)->get_closest_2d_dis(pts.at(i))
```
which computes 2D wire distances using a single implicit TPC geometry.

The toolkit resolves the APA and face for each test point first:
```cpp
auto test_wpid = dv->contained_by(test_p);
int apa  = test_wpid.apa();
int face = test_wpid.face();
// ...
auto [dist_u, dist_v, dist_w] = segment_get_closest_2d_distances(*it1, test_p, apa, face, "fit");
```

This correctly handles detectors with multiple APAs/faces (e.g., DUNE PDHD). The function also short-circuits early if `apa == -1 || face == -1` (point outside any known volume), which the prototype could not do.

The toolkit also uses `fits()` (track-fitted points) rather than the raw `get_point_vec()` that the prototype uses, which is the appropriate cloud for closeness checks after fitting.

### 4.2 `determine_direction`: distance-based vertex matching

The prototype identifies which vertex is `start_v` vs `end_v` by comparing wire-cell-point indices:
```cpp
if ((*it)->get_wcpt().index == sg->get_wcpt_vec().front().index) start_v = *it;
```

The toolkit uses Euclidean distance instead:
```cpp
double dis_sv_front = ray_length(Ray{start_v->wcpt().point, front_pt});
double dis_sv_back  = ray_length(Ray{start_v->wcpt().point, back_pt});
if (dis_sv_front > dis_sv_back) std::swap(start_v, end_v);
```

Distance-based matching is robust across APAs where global wire indices may not be contiguous or unique across faces.

### 4.3 `examine_good_tracks`: distance-based endpoint detection

Same pattern as `determine_direction` above — the toolkit uses `ray_length` comparisons to determine `start_vertex` / `end_vertex` rather than point-index matching.

### 4.4 `examine_all_showers`: cluster identity via flag not pointer

The prototype compares `temp_cluster == main_cluster` (pointer identity). The toolkit uses:
```cpp
bool is_main_cluster = cluster.get_flag(Facade::Flags::main_cluster);
```

This is the correct approach in a multi-cluster system where "main cluster" status is encoded as a flag on the cluster object rather than by maintaining a separate pointer.

---

## 5. Per-function Detailed Notes

### 5.1 `separate_track_shower`

Functionally identical. Iterates all edges in the graph belonging to the cluster; calls `segment_is_shower_topology` first, then `segment_is_shower_trajectory` only if the segment is not a topology shower. Timing via `m_perf` flag replaces the prototype's unconditional stdout.

Note: The prototype has a second overload `separate_track_shower()` (no args) that loops all clusters and populates `point_flag_showers`. This functionality is handled separately in the toolkit (not part of this function group).

### 5.2 `determine_direction`

Functionally identical logic. Key differences:
- Distance-based `flag_start` determination (see §4.2).
- Topology shower branch: toolkit calls `segment_determine_shower_direction` then immediately sets electron particle info with `particle_score(100)`. This matches what the prototype's `determine_dir_shower_topology` does internally.
- Trajectory shower branch: calls `segment_determine_shower_direction_trajectory` passing `43000/units::cm` as the reference dQ/dx. The prototype passes the same value via `TPCParams`.

### 5.3 `calculate_num_daughter_showers`

BFS traversal is identical. The toolkit's `is_shower` check (`kShowerTrajectory || kShowerTopology || abs(pdg)==11`) correctly maps to the prototype's `get_flag_shower()`.

One structural difference: the prototype inserts `curr_vertex` into `used_vertices` after pushing all its neighbours onto `temp_segments`, while the toolkit inserts it before. Both are correct because `used_vertices` is checked at the start of the outer loop (not per-neighbour push), so duplicates in the queue are harmlessly skipped.

### 5.4 `examine_good_tracks`

Functionally identical. Condition checked:
```
(num_daughter_showers >= 4 || length_daughter_showers > 50cm && num_daughter_showers >= 2)
  && (max_angle > 155 || drift_angle < 15 && min_para_angle < 15 && sum < 25)
  && length < 15cm
```
maps precisely to prototype line 256. On trigger, sets particle to electron, clears direction (`dirsign=0`), marks weak.

### 5.5 `fix_maps_multiple_tracks_in`

Functionally identical. At each vertex with >1 segment: if `n_in > 1 && n_in != n_in_shower`, resets all incoming non-shower tracks to `dirsign=0`, `dir_weak=true`. No particle-type change (intentionally — this is a direction reset, not a type change).

### 5.6 `fix_maps_shower_in_track_out`

Functionally identical. At each vertex: if any incoming segment is a shower AND any outgoing non-shower track has strong direction (`!dir_weak`), flip all incoming showers' `dirsign` by -1 and mark weak. No particle-type change.

### 5.7 `improve_maps_one_in`

Functionally identical (with determinism fix). Propagates direction: if ≥1 incoming segment exists at a vertex (strong, or any depending on `flag_strong_check`), all undirected/weak segments at that vertex are set to point out. The toolkit recalculates 4-momentum for changed segments via `segment_cal_4mom`, matching the prototype's `cal_4mom()` call.

### 5.8 `improve_maps_shower_in_track_out`

Functionally identical (with determinism fix). When an incoming shower exists:
- Outgoing weak tracks → set particle to electron (PDG 11), clear direction.
- No-direction segments → set to electron if not already shower; also recalculates 4-mom for showers with valid energy (matching prototype `if (sg1->get_particle_4mom(3)>0) sg1->cal_4mom()`).

### 5.9 `improve_maps_no_dir_tracks`

All eight labelled cases (A–H) port correctly:

| Case | Condition | Action |
|---|---|---|
| A | >2 showers total or all-shower vertex on both sides, length <5cm | Promote to electron |
| B | Segment between shower-dominated vertices (>2 showers total), length <25cm or no direction | Set direction, promote to electron |
| C/D | Proton with all-shower vertex on one side (≥2 showers), dQ/dx rms threshold met | Set direction, promote to electron |
| E | Muon with all-shower vertex, short/kinked trajectory | Promote to electron |
| F | Muon with many daughter showers and large opening angle | Promote to electron |
| G | Muon with single-segment vertex neighbour that is much longer | Promote to electron |
| H | Unknown type, high median dQ/dx, near showers with small angle | Promote to electron |

One precision note: in cases C and D, the toolkit uses `segment_rms_dQ_dx(sg)` which maps to `get_rms_dQ_dx()`. In case H, `segment_median_dQ_dx(sg)` maps to `get_medium_dQ_dx()`. Both are correct.

### 5.10 `improve_maps_multiple_tracks_in`

After the bug fix (§2.1), functionally identical to prototype. At each vertex with multiple incoming non-shower segments, all incoming tracks are promoted to electron (PDG 11).

### 5.11 `judge_no_dir_tracks_close_to_showers`

Logic identical. For each no-direction (dirsign==0) non-shower track, checks every point in `fits()` against all shower segments. If any point has minimum 2D wire distance to showers >0.6cm in any plane, the track is not promoted. Otherwise it becomes an electron. Enhanced for multi-APA as described in §4.1.

### 5.12 `examine_maps`

Functionally identical. Checks two physics constraints at every vertex with >1 segment:
1. Multiple non-shower particles entering the same vertex → print warning, return false.
2. Showers entering with strong tracks exiting → print warning, return false.

The function is used as a diagnostic/sanity check; its return value is not used to halt processing.

### 5.13 `examine_all_showers`

Functionally identical (with `is_main_cluster` flag fix, §4.4). Decision tree:

1. Single good track: check daughter showers and opening angles from each endpoint; check beam direction against showers at the connected multi-segment vertex. If conditions met, demote the good track.
2. No good tracks, exactly 2 weak tracks ≤35cm: check beam direction of the longer track against its more-connected vertex. If anti-beam, promote that track to shower counts.
3. If `n_good_tracks == 0` after the above: apply length-ratio thresholds to decide whether to convert all remaining tracks to electrons.

The main-cluster branch (`is_main_cluster`) controls a stricter condition in the `length_tracks < 35cm` block, exactly matching the prototype's `temp_cluster == main_cluster` pointer comparison.

### 5.14 `shower_determining_in_main_cluster`

Call sequence is identical to prototype `shower_determing_in_main_cluster` (note: typo fixed in toolkit name):

```
examine_good_tracks
fix_maps_multiple_tracks_in
fix_maps_shower_in_track_out         ← 1st pass
improve_maps_one_in
improve_maps_shower_in_track_out     ← 1st pass (flag_strong_check=true)
improve_maps_no_dir_tracks
improve_maps_shower_in_track_out     ← 2nd pass (flag_strong_check=false)
improve_maps_multiple_tracks_in
fix_maps_shower_in_track_out         ← 2nd pass
judge_no_dir_tracks_close_to_showers
examine_maps
examine_all_showers
```

### 5.15 `print_segs_info`

Functionally identical. Prints per-segment: ID, length (cm), type (S_topo/S_traj/Track), dirsign, PDG, mass (MeV), KE (MeV), dir_weak, in/out status relative to optional vertex. Uses `std::cout` directly (appropriate for debug output), with a `print_segs_info:` prefix added in the toolkit for log parsing.

---

## 6. Efficiency Notes

### 6.1 Repeated `is_shower` pattern

The three-part `is_shower` check:
```cpp
sg->flags_any(SegmentFlags::kShowerTrajectory) ||
sg->flags_any(SegmentFlags::kShowerTopology) ||
(sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 11)
```
appears ~30 times across this file. It is not performance-critical (each check is O(1)), but a helper method `segment_is_shower(SegmentPtr)` could reduce maintenance risk and make future is_shower definition changes easier. Not blocking.

### 6.2 `judge_no_dir_tracks_close_to_showers`: O(S × T × P) scan

For each no-direction track segment, for each of its fit points, all shower segments are scanned. Complexity is O(shower_count × track_count × points_per_track). In practice the counts are small (10s of segments, 100s of points) so this is not a bottleneck. A spatial index (k-d tree on shower fit points) would improve asymptotic cost but is not warranted here.

### 6.3 `improve_maps_no_dir_tracks`: `edge_range1` / `edge_range2` re-use

The two `boost::out_edges` ranges (`edge_range1`, `edge_range2`) are computed once per segment and reused across cases C/D, E/F, and H. This avoids redundant Boost graph lookups and is already efficient.

---

## 7. Summary Table

| Function | Equivalent? | Notes |
|---|---|---|
| `separate_track_shower` | ✅ Yes | — |
| `determine_direction` | ✅ Yes | Distance-based vertex matching (multi-APA safe) |
| `calculate_num_daughter_showers` | ✅ Yes | — |
| `examine_good_tracks` | ✅ Yes | Distance-based endpoint detection |
| `fix_maps_multiple_tracks_in` | ✅ Yes | — |
| `fix_maps_shower_in_track_out` | ✅ Yes | — |
| `improve_maps_one_in` | ✅ Yes | `ordered_nodes` for determinism |
| `improve_maps_shower_in_track_out` | ✅ Yes | `ordered_nodes` for determinism |
| `improve_maps_no_dir_tracks` | ✅ Yes | Cases A–H verified |
| `improve_maps_multiple_tracks_in` | ✅ Yes (fixed) | Bug fixed: missing electron check in `is_shower` |
| `judge_no_dir_tracks_close_to_showers` | ✅ Yes | Multi-APA enhanced (apa/face-aware 2D distances) |
| `examine_maps` | ✅ Yes | — |
| `examine_all_showers` | ✅ Yes | `is_main_cluster` flag replaces pointer comparison |
| `shower_determining_in_main_cluster` | ✅ Yes | Call sequence identical; typo fixed in name |
| `print_segs_info` | ✅ Yes | — |
