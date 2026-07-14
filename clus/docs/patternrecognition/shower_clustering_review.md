# Shower Clustering Review

## Scope

Line-by-line comparison of the prototype shower-clustering pipeline with the
toolkit port.  The prototype lives in two files:

| Prototype file | Functions |
|---|---|
| `prototype_base/pid/src/NeutrinoID_shower_clustering.h` | `shower_clustering_with_nv`, `shower_clustering_with_nv_in_main_cluster`, `shower_clustering_connecting_to_main_vertex`, `shower_clustering_with_nv_from_main_cluster`, `shower_clustering_with_nv_from_vertices`, `examine_merge_showers`, `shower_clustering_in_other_clusters`, `id_pi0_with_vertex`, `id_pi0_without_vertex`, `update_shower_maps` |
| `prototype_base/pid/src/NeutrinoID_em_shower.h` | `examine_showers`, `examine_shower_1` |
| `prototype_base/pid/src/WCShower.cxx` | `complete_structure_with_start_segment` |
| `prototype_base/pid/src/PR3DCluster.cxx` | `fill_2d_charge_dead_chs` |

The toolkit equivalents are in:

| Toolkit file | Functions |
|---|---|
| `clus/src/NeutrinoShowerClustering.cxx` | All `NeutrinoPatternAlgos::` shower clustering functions (3310 lines) |
| `clus/src/PRShower.cxx` | `Shower::complete_structure_with_start_segment` (line 269) |
| `clus/src/TrackFitting.cxx` | `collect_2D_charge` (toolkit equivalent of `fill_2d_charge_dead_chs`) |

---

## Implementation Status

| ID | Issue | Status | Location |
|---|---|---|---|
| B.1 | `shower_clustering_connecting_to_main_vertex`: `calculate_num_daughter_showers` called with `false` instead of `true` | **Fix** | `NeutrinoShowerClustering.cxx:280` |
| B.2 | `id_pi0_without_vertex`: `continue` should be `break` when both showers are short | **Fix** | `NeutrinoShowerClustering.cxx:2969` |
| B.3 | `id_pi0_without_vertex`: missing `pdg==11` in `get_flag_shower` equivalent | **Fix** | `NeutrinoShowerClustering.cxx:2850-2851` |
| B.4 | `shower_clustering_in_other_clusters` sub-pass 1: missing second `update_particle_type()` and `calculate_kinematics()` after shower merge | **Fix** | `NeutrinoShowerClustering.cxx:1461-1463` |
| B.5 | `id_pi0_with_vertex`: direction calculation uses candidate pi0 vertex instead of shower's own start vertex | **Fix** | `NeutrinoShowerClustering.cxx:2592-2597` |
| B.6 | `shower_clustering_in_other_clusters` sub-pass 1: extra electron-forcing block not in prototype | Documented only | `NeutrinoShowerClustering.cxx:1421-1435` |
| L.1 | `shower_clustering_with_nv_in_main_cluster`: hardcoded electron mass `0.511*MeV` | **Fix** | `NeutrinoShowerClustering.cxx:204` |

---

## Per-Function Walk-Through

### 1. `shower_clustering_with_nv` (entry point)

**Prototype**: `NeutrinoID_shower_clustering.h:268-377`
**Toolkit**: `NeutrinoShowerClustering.cxx:3086-3245`

The orchestration function calls 11 sub-functions in the same order as the
prototype.  The toolkit adds:

- `collect_charge_maps(track_fitter)` before `calculate_shower_kinematics` to
  pre-populate 2D charge maps once per event (efficiency improvement).
- `spdlog`-based timing diagnostics guarded by `m_perf`.
- A debug block printing `dirsign`/`dir_weak`/flags for main-cluster segments.

**Verdict**: Correct.  No bugs.

---

### 2. `update_shower_maps`

**Prototype**: `NeutrinoID_shower_clustering.h:1388-1403`
**Toolkit**: `NeutrinoShowerClustering.cxx:30-74`

Clears four maps, iterates showers, populates `map_vertex_to_shower` /
`map_vertex_in_shower` / `map_segment_in_shower` / `used_shower_clusters`.
The toolkit adds null guards on `start_vtx`, `seg`, and `seg->cluster()`.
Uses `IndexedShowerSet` (ident-ordered) instead of pointer-keyed `std::set`.

**Verdict**: Correct.  No bugs.  Improved null-safety.

---

### 3. `shower_clustering_with_nv_in_main_cluster`

**Prototype**: `NeutrinoID_shower_clustering.h:1654-1771`
**Toolkit**: `NeutrinoShowerClustering.cxx:76-218`

BFS from main vertex looking for EM-shower-like sub-trees to convert from
muon to electron.

- **Prototype bug fixed**: line 1688 `fabs(curr_sg->get_particle_type()==13)` is
  a comparison-inside-fabs bug; the toolkit correctly uses
  `std::abs(pdg()) == 13`.
- **Toolkit addition**: `update_particle_type()` and PDG=0 guard (lines 143-158)
  have no prototype equivalent — intentional enhancement for segments that
  never received a PID.
- **L.1**: Hardcoded `0.511 * units::MeV` at line 203 should ideally use
  `particle_data->get_particle_mass(11)` for consistency, but the value is
  correct.
- Determinism: uses `sorted_out_edges`, `IndexedSegmentSet`, max-muon
  tiebreaker by `sg1->id()`.

**Verdict**: Correct.  Prototype bug fixed.  One minor style issue (L.1).

---

### 4. `shower_clustering_connecting_to_main_vertex`

**Prototype**: `NeutrinoID_shower_clustering.h:114-266`
**Toolkit**: `NeutrinoShowerClustering.cxx:221-458`

Creates showers from segments at the main vertex that look like electrons.

- **B.1 (Bug)**: Line 280 calls
  `calculate_num_daughter_showers(graph, main_vertex, sg, false)`.
  The prototype (line 140) passes `true` — meaning "count all segments
  including showers".  The toolkit's default for `flag_count_shower` is `true`
  (header line 134), but the call site explicitly passes `false`.
  **Fix**: remove the `false` argument (or change to `true`).
- Optimization: the `map_segment_in_shower` check is moved before the
  expensive BFS traversal (lines 277-278) — correct since neither has side
  effects.
- Operator-precedence: dQ/dx skip condition (prototype line 146) is correctly
  parenthesized in the toolkit (lines 293-296).
- End-vertex determination uses distance-based matching instead of index
  matching — acceptable.

**Verdict**: One bug (B.1).

---

### 5. `shower_clustering_with_nv_from_main_cluster`

**Prototype**: `NeutrinoID_shower_clustering.h:1775-1939`
**Toolkit**: `NeutrinoShowerClustering.cxx:461-730`

Creates showers from clusters in `map_cluster_main_vertices` (non-main).
Four-step process: initial BFS, segment addition, direction filtering,
iterative segment addition.

- Operator precedence in the `flag_shower` condition (prototype line 1806) is
  correctly made explicit in the toolkit (line 525-526).
- Angle computation with `fabs(dot)` folding into [0,90°] is mathematically
  equivalent to prototype's `TVector3::Angle()` with `abs(angle-90)<5` check.
- `map_shower_dir` iteration sorted by segment ID (lines 649-653) for
  determinism.

**Verdict**: Correct.  No bugs.

---

### 6. `shower_clustering_with_nv_from_vertices`

**Prototype**: `NeutrinoID_shower_clustering.h:995-1386`
**Toolkit**: `NeutrinoShowerClustering.cxx:733-1180`

Largest function (~400 lines).  Creates showers from non-main-vertex segments
and clusters segments from other clusters.

- **Prototype bug fixed**: `get_start_point()` returns `(0,0,0)` in the
  prototype (line 1352) because `set_start_point` is never called before use.
  The toolkit (lines 1154-1155) correctly uses the start vertex fit point.
- All operator-precedence-sensitive conditions correctly parenthesized.
- `np > 0` guard added before division (line 792) — defensive improvement.
- Multi-APA: uses `dv->contained_by()` and `compute_wireplane_params()` for
  point cloud building.

**Verdict**: Correct.  Prototype bug fixed.

---

### 7. `examine_merge_showers`

**Prototype**: `NeutrinoID_shower_clustering.h:380-426`
**Toolkit**: `NeutrinoShowerClustering.cxx:1182-1260`

Merges type-1 (conn_type=1) showers with type-2 (conn_type=2) when direction
angle < 10°.

- Early-exit optimization when either type1 or type2 list is empty (line 1216).
- `claimed` set is `unordered_set` but only used for membership testing — no
  iteration-order dependence.
- Particle-type filter correctly applied during pre-classification.

**Verdict**: Correct.  No bugs.

---

### 8. `shower_clustering_in_other_clusters`

**Prototype**: `NeutrinoID_shower_clustering.h:1442-1651`
**Toolkit**: `NeutrinoShowerClustering.cxx:1263-1593`

Two sub-passes: (1) clusters with vertices in `map_cluster_main_vertices`,
(2) remaining other clusters.

- **B.4 (Bug)**: Sub-pass 1 (toolkit lines 1438-1463): after creating the
  shower and merging with existing showers, the prototype calls
  `update_particle_type()` a second time (line 1555) and
  `calculate_kinematics()` (line 1556) before pushing the shower.
  The toolkit only calls `update_particle_type()` once (line 1438, before
  merge) and omits both the post-merge `update_particle_type()` and
  `calculate_kinematics()`.
  **Fix**: add `shower->update_particle_type(...)` and
  `shower->calculate_kinematics(...)` after the merge loop (after line 1461).
- **B.6 (Documented)**: Sub-pass 1 (toolkit lines 1421-1435) has an
  electron-forcing block that does not exist in the prototype's sub-pass 1.
  It mirrors the pattern in sub-pass 2 and `shower_clustering_with_nv_from_vertices`.
  This is an intentional extension — documented but not removed.
- Sub-pass 2 (lines 1470+) correctly matches the prototype.
- Operator precedence on angle/distance conditions correctly parenthesized.

**Verdict**: One bug (B.4).  One intentional extension (B.6).

---

### 9. `examine_showers`

**Prototype**: `NeutrinoID_em_shower.h:1-334`
**Toolkit**: `NeutrinoShowerClustering.cxx:2021-2375`

Three cases for converting segment-to-shower and shower-to-segment based on
connectivity and dQ/dx.

- Minor logic difference: `abs(pdg) != 13` (line 2113) also excludes
  antimuons (-13) unlike prototype's `particle_type != 13`.  Low risk.
- `add_shower` call uses dereferenced shower (`*shower1`); this matches the
  toolkit's `add_shower(Shower& other)` signature (reference, not pointer).
  Correct.
- Determinism: uses `SegmentIndexCmp` for `map_merge_seg_shower`,
  `shower_cmp` for sorted shower vectors, `ordered_edges` for graph
  traversal.

**Verdict**: Correct.  No bugs.

---

### 10. `examine_shower_1`

**Prototype**: `NeutrinoID_em_shower.h:337-652`
**Toolkit**: `NeutrinoShowerClustering.cxx:1596-2019`

Reclassifies segments near main vertex: converts certain track-like segments
connected to showers into shower segments.  Two-part algorithm.

- **Prototype bug fixed**: line 501 `flag_good_track == true;` is a no-op
  comparison (should be assignment).  Toolkit line 1831 correctly uses `=`.
- Skip-filter condition (prototype line 368) with complex operator precedence
  is correctly decomposed in the toolkit (lines 1653-1662).  After analysis,
  the net effect is equivalent because the `!dir_weak` clause is the first
  filter.
- `is_shower` check at line 1839 checks `kShowerTrajectory || kShowerTopology`
  but not `pdg==11`.  In context, this matches the prototype's
  `get_flag_shower()` because the relevant segments have already been filtered
  to have `dir_weak()==true`, and segments with `pdg==11` would already have
  `kShowerTrajectory` set.
- `map_shower_showers` keyed by `ShowerPtr` — pointer-ordered, but the
  max-energy selection uses a segment-ID tiebreaker (lines 1988-1991).
- `complete_structure_with_start_segment` takes only `used_segments`; the
  toolkit's `Shower` navigates the graph internally, so no external map
  parameter is needed.

**Verdict**: Correct.  Prototype bug fixed.

---

### 11. `id_pi0_with_vertex`

**Prototype**: `NeutrinoID_shower_clustering.h:735-993`
**Toolkit**: `NeutrinoShowerClustering.cxx:2378-2716`

Identifies pi0 candidates from shower pairs sharing a vertex.

- **B.5 (Documented)**: Lines 2588-2592: when `dis < 3cm`, the prototype
  computes the shower direction from the shower's own start vertex
  (`shower_1->get_start_vertex().first->get_fit_pt()`, prototype line 870).
  The toolkit uses the candidate pi0 decay vertex (`vtx_pt`).  These differ
  when the candidate vertex is not the shower's start vertex, but both are
  within 3cm of the shower's start point, so the practical impact is minimal.
  Documented for awareness.
- `flag_start` determination uses distance-based matching instead of index
  matching — acceptable.
- Determinism: `IndexedVertexSet`, `IndexedShowerSet`, `shower_pair_cmp` used
  throughout.

**Verdict**: One documented difference (B.5).

---

### 12. `id_pi0_without_vertex`

**Prototype**: `NeutrinoID_shower_clustering.h:428-732`
**Toolkit**: `NeutrinoShowerClustering.cxx:2718-3084`

Identifies pi0 candidates from disconnected shower pairs using closest-approach
geometry.

- **B.2 (Bug)**: Line 2969: `continue` should be `break`.  When both showers
  in a pair are short (length ≤ 15cm), the prototype (line 614) uses `break`
  to exit the inner loop entirely, skipping all remaining pairs for that
  `shower_1`.  The toolkit uses `continue`, which only skips the current pair
  and keeps checking other `shower_2` candidates.  This changes behavior when
  a short shower is paired with multiple others.
  **Fix**: change `continue` to `break`.
- **B.3 (Bug)**: Lines 2850-2851: the prototype's `get_flag_shower()` returns
  `flag_shower_trajectory || flag_shower_topology || (|particle_type| == 11)`.
  The toolkit only checks `kShowerTrajectory` and `kShowerTopology`, missing
  the `pdg==11` case.  Showers whose start segment has electron PDG but
  neither topology nor trajectory flag set will be incorrectly skipped.
  **Fix**: add `|| std::abs(pdg) == 11` to the condition.
- `pair_less` comparator and sorted `shower_less` maps ensure determinism.
- `std::clamp` and zero-magnitude guards on theta/phi are robustness
  improvements.

**Verdict**: Two bugs (B.2, B.3).

---

### 13. `complete_structure_with_start_segment`

**Prototype**: `WCShower.cxx:703-756`
**Toolkit**: `PRShower.cxx:269-362`

Worklist BFS that adds all segments reachable from the start segment
(excluding the start vertex sub-tree) to the shower.

- The start segment is added to the view by `set_start_segment()` (line 88),
  which calls `TrajectoryView::add_segment(seg)`.
  `complete_structure_with_start_segment` only needs to mark it as used.
- The toolkit navigates the graph internally via `find_vertices()` and
  `sorted_out_edges()`, eliminating the need for external
  `map_vertex_segments` / `map_segment_vertices` parameters.
- Uses `IndexedSegmentSet` for deterministic iteration.

**Verdict**: Correct.  No bugs.

---

### 14. `fill_2d_charge_dead_chs` → `TrackFitting::collect_2D_charge`

**Prototype**: `PR3DCluster.cxx:465-521`
**Toolkit**: `TrackFitting.cxx:778-953`

The prototype iterates blobs in a single cluster, checks each plane for "bad"
status, distributes charge evenly across wires.

Toolkit improvements:
- **Multi-slice support**: handles blobs spanning multiple time slices by
  dividing charge by `num_wires * num_slices`.
- **Flag-guarded accumulation**: only accumulates onto existing entries when
  `flag == 0` (dead channel), preventing dead-channel synthetic charge from
  overwriting live measurements.
- **Multi-APA/face**: correctly extracts `apa`/`face` from each blob's
  `wpid()` and uses `fetch_channel_from_anode(apa, face, plane, wire_index)`.
- **Parameterized uncertainties**: uses `m_params.rel_charge_uncer` and
  `m_params.add_charge_uncer` instead of hardcoded `0.1` and `600`.

**Verdict**: Correct.  Properly multi-APA generalized.

---

## Bugs Summary

### B.1 — `calculate_num_daughter_showers` called with wrong flag

**File**: `NeutrinoShowerClustering.cxx:280`
**Severity**: Medium — changes which segments are filtered out for shower creation.

The prototype passes `true` to `calculate_num_daughter_tracks(main_vertex, sg, true)`,
meaning "count all downstream segments including shower-flagged ones".
The toolkit passes `false`, excluding shower segments from the count.
This can produce a different `pair_result.first` value, affecting the
`medium_dQ_dx_1 > 1.45 && pair_result.first <= 3` proton-skip condition.

**Fix**: Remove the `false` argument (default is `true`).

### B.2 — `continue` instead of `break` in `id_pi0_without_vertex`

**File**: `NeutrinoShowerClustering.cxx:2969`
**Severity**: Medium — changes pi0 candidate selection for short showers.

**Fix**: Change `continue` to `break`.

### B.3 — Missing `pdg==11` in shower flag check

**File**: `NeutrinoShowerClustering.cxx:2850-2851`
**Severity**: Low-Medium — only affects segments that have electron PDG but
neither `kShowerTrajectory` nor `kShowerTopology` flags.

**Fix**: Add `|| std::abs(pdg) == 11` to the condition.

### B.4 — Missing post-merge `update_particle_type` and `calculate_kinematics`

**File**: `NeutrinoShowerClustering.cxx:1461-1463`
**Severity**: Medium — after merging with existing showers, the particle type
is not re-evaluated and kinematics are not recomputed.

**Fix**: Add `shower->update_particle_type(particle_data, recomb_model)` and
`shower->calculate_kinematics(particle_data, recomb_model)` after the merge
loop and before `showers.insert(shower)`.

---

## Documented Differences (not bugs)

### B.5 — Direction calculation in `id_pi0_with_vertex`

**File**: `NeutrinoShowerClustering.cxx:2592-2597`
**Severity**: Low-Medium — affects theta/phi/angle BDT features for pi0.

When `dis < 3cm`, the prototype computes the shower direction from the
shower's own start vertex.  The toolkit was using the candidate pi0 decay
vertex instead.  **Fixed**: now uses each shower's own start vertex point
(`sv1_pt`, `sv2_pt`) to match the prototype.

### B.6 — Extra electron-forcing in `shower_clustering_in_other_clusters` sub-pass 1

**File**: `NeutrinoShowerClustering.cxx:1421-1435`

This block forces PDG=0 or short+weak muon segments to electron before the
majority-vote `update_particle_type()`.  The prototype's sub-pass 1 does not
have this; it only appears in sub-pass 2.  The toolkit applies it
consistently in both sub-passes — an intentional improvement for robustness.

### L.1 — Hardcoded electron mass

**File**: `NeutrinoShowerClustering.cxx:203`

`0.511 * units::MeV` is correct but should ideally use
`particle_data->get_particle_mass(11)` for consistency with the dependency
injection pattern used elsewhere.

---

## Determinism Improvements

| Function | Prototype issue | Toolkit fix |
|---|---|---|
| `update_shower_maps` | Pointer-keyed `std::set<WCShower*>` | `IndexedShowerSet` (ident-ordered) |
| `shower_clustering_connecting_to_main_vertex` | `map_vertex_segments[main_vertex]` pointer order | `ordered_edges(graph)` + `seg_order` |
| `shower_clustering_with_nv_from_main_cluster` | `map_vertex_segments` pointer order | `ordered_edges(graph)` + segment ID tiebreaker |
| `shower_clustering_with_nv_from_vertices` | `map_vertex_segments` / `map_segment_vertices` pointer order | `ordered_edges(graph)` / `ordered_nodes(graph)` / `seg_order` |
| `shower_clustering_in_other_clusters` | `map_vertex_segments` / `map_segment_vertices` pointer order | `ordered_nodes(graph)` / `seg_order` |
| `examine_showers` | Pointer-keyed `map<WCShower*, ...>` | `SegmentIndexCmp` / `shower_cmp` sorted vectors |
| `examine_shower_1` | Pointer-keyed `map<WCShower*, set<WCShower*>>` | Segment-ID tiebreaker for max-energy selection |
| `id_pi0_with_vertex` | Pointer-keyed shower pairs | `IndexedVertexSet` / `IndexedShowerSet` / `shower_pair_cmp` |
| `id_pi0_without_vertex` | Pointer-keyed `map<WCShower*, ...>` | `shower_less` / `pair_less` comparators |

---

## Multi-APA / Face Correctness

Most shower clustering functions operate on graph topology, 3D geometry, and
segment properties — they do not directly touch wire/tick coordinates.
The multi-APA critical functions are:

1. **`fill_2d_charge_dead_chs` → `collect_2D_charge`**: Correctly handles
   multiple APAs/faces via `blob->wpid()` extraction, `fetch_channel_from_anode`,
   and `CoordReadout` keying that includes the `apa` field.

2. **`shower_clustering_with_nv_from_vertices`**: Uses `dv->contained_by()`
   and `compute_wireplane_params()` when building point clouds for multi-APA
   geometry.

3. **Cluster distance queries** (`get_closest_dis`): Facade method handles
   multi-TPC internally.

No multi-APA/face issues were found in any of the reviewed functions.

---

## Prototype Bugs Fixed by Toolkit

| Prototype line | Bug | Toolkit fix |
|---|---|---|
| `NeutrinoID_shower_clustering.h:1688` | `fabs(curr_sg->get_particle_type()==13)` — comparison inside fabs | `std::abs(pdg()) == 13` |
| `NeutrinoID_shower_clustering.h:1352` | `get_start_point()` returns `(0,0,0)` (never set) | Uses start vertex fit point |
| `NeutrinoID_em_shower.h:501` | `flag_good_track == true;` — comparison instead of assignment | `flag_good_track = true;` |

---

## Summary of Changes

- **6 bugs fixed** (B.1–B.5, L.1) in `NeutrinoShowerClustering.cxx`
- **1 documented difference** (B.6) preserved intentionally
- **3 prototype bugs** already fixed in the toolkit port
- **9 determinism improvements** over prototype
- **Multi-APA/face** handling verified correct across all functions
