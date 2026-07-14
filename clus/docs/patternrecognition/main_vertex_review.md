# Main Vertex Determination Review

**Reviewer:** Claude  
**Date:** 2026-04-12  
**Branch:** apply-pointcloud  
**Prototype source:** `NeutrinoID_track_shower.h`, `NeutrinoID.cxx`, `NeutrinoID_improve_vertex.h`  
**Toolkit source:** `clus/src/NeutrinoVertexFinder.cxx`, `clus/src/NeutrinoPatternBase.cxx`

---

## 1. Function Map

| Prototype (NeutrinoID) | Toolkit (PatternAlgorithms) | File | Line | Status |
|---|---|---|---|---|
| `examine_main_vertex_candidate` | `examine_main_vertex_candidate` | `NeutrinoVertexFinder.cxx` | 234 | ✓ Equivalent |
| `determine_main_vertex` | `determine_main_vertex` | `NeutrinoVertexFinder.cxx` | 2362 | ✓ Better |
| `compare_main_vertices_all_showers` | `compare_main_vertices_all_showers` | `NeutrinoVertexFinder.cxx` | 317 | ✓ Equivalent |
| `compare_main_vertices` | `compare_main_vertices` | `NeutrinoVertexFinder.cxx` | 719 | ✓ Better |
| `calc_conflict_maps` | `calc_conflict_maps` | `NeutrinoVertexFinder.cxx` | 503 | ✓ Equivalent |
| `examine_direction` | `examine_direction` | `NeutrinoVertexFinder.cxx` | 1018 | ✓ Better |
| `find_cont_muon_segment` | `find_cont_muon_segment` | `NeutrinoVertexFinder.cxx` | 920 | ✓ Equivalent |
| `eliminate_short_vertex_activities` | `eliminate_short_vertex_activities` | `NeutrinoVertexFinder.cxx` | 1648 | ✓ Better |
| `change_daughter_type` | `change_daughter_type` | `NeutrinoVertexFinder.cxx` | 2583 | ✓ Equivalent |
| `examine_main_vertices(vertices)` | `examine_main_vertices_local` | `NeutrinoVertexFinder.cxx` | 2679 | ✓ Equivalent + bug fixed |
| `examine_main_vertices()` (global) | `examine_main_vertices` | `NeutrinoPatternBase.cxx` | 1922 | ✓ Better |

> **Note on naming:** The prototype's `examine_main_vertices` has two overloads:
> - The zero-argument global version → toolkit `examine_main_vertices` in `NeutrinoPatternBase.cxx`
> - The `(ProtoVertexSelection& vertices)` local candidate-filtering version → toolkit `examine_main_vertices_local` in `NeutrinoVertexFinder.cxx`
>
> `examine_main_vertex_candidate` (singular, returns `tuple<bool,int,int>`) is an internal helper; it is not in the user's list as a standalone but is called by several of the reviewed functions.

---

## 2. Bug Found and Fixed

### 2.1 `examine_main_vertices_local` — non-deterministic `tmp_vertices` set

**File:** `clus/src/NeutrinoVertexFinder.cxx` line 2684 (before fix)

**Problem:** The local accumulator `tmp_vertices` was declared as `std::set<VertexPtr>`. The default comparator for `std::shared_ptr` is pointer-address comparison, which varies across runs and produces non-deterministic iteration order when the result is copied back into the output `vertices` vector.

```cpp
// BEFORE (non-deterministic):
std::set<VertexPtr> tmp_vertices;
```

```cpp
// AFTER (fixed):
IndexedVertexSet tmp_vertices;  // order by stable graph index, not pointer address
```

`IndexedVertexSet` is `std::set<VertexPtr, VertexIndexCmp>` where `VertexIndexCmp` orders by `vtx->get_graph_index()`, a stable per-run integer. This ensures that when `vertices` is rebuilt from `tmp_vertices`, the ordering is the same across runs for the same input.

**Impact:** The selected `main_vertex` (or intermediate vertex added via `find_cont_muon_segment` path) could vary across runs when multiple vertices survive the back-to-back filtering, leading to non-reproducible event classification.

> The prototype's `examine_main_vertices(vertices)` also used `std::set<ProtoVertex*>` (pointer-ordered) and had the same non-determinism. The toolkit's fix is an improvement beyond the prototype.

---

## 3. Determinism Improvements (Toolkit over Prototype)

| Function | Prototype Non-Determinism | Toolkit Fix |
|---|---|---|
| `determine_main_vertex` | Iterates `map_vertex_segments` (pointer-ordered `std::map`) | Uses `ordered_nodes(graph)` (stable graph-index order) |
| `compare_main_vertices` | Iterates `map_segment_vertices` / `map_vertex_segments` | Uses `boost::edges(graph)` + `boost::out_edges(vd,graph)` with consistent index-based scoring |
| `examine_direction` | `flag_shower_in` detection iterates all segments at `prev_vtx` regardless of BFS wave | Only segments in `used_segments` (already-processed earlier BFS levels) are eligible as incoming showers; sibling-order dependence eliminated |
| `examine_main_vertices_local` | `std::set<ProtoVertex*>` (pointer-ordered) | `IndexedVertexSet` (graph-index-ordered) — **this session's bug fix** |
| `examine_main_vertices` (global) | None (different structure) | Sorts candidates by `cluster_id` before processing; uses `sorted_candidates` vector |

### 3.1 Detail: `examine_direction` `flag_shower_in` fix

The prototype's BFS loop checked all segments at `prev_vtx` to determine whether an incoming shower exists:

```cpp
// Prototype: scans all prev_vtx segments including unprocessed siblings
for (auto it1 = map_vertex_segments[prev_vtx].begin(); ... ) {
    if (flag_start && sg->get_flag_dir()==-1 || ...) {
        if (sg->get_flag_shower()) { flag_shower_in = true; break; }
    }
}
```

If two sibling segments were processed in the same BFS wave, neither was yet in `used_segments` when the other's `flag_shower_in` was evaluated. Which sibling was processed first was determined by pointer ordering — non-deterministic.

The toolkit correctly guards with `used_segments.find(sg) == used_segments.end()) continue` (line 1134), ensuring only segments from earlier BFS levels can contribute to `flag_shower_in`. This also matches the semantic intent: a shower is "incoming" only if it was already assigned a direction pointing toward `prev_vtx`.

---

## 4. Multi-APA Correctness

### 4.1 `eliminate_short_vertex_activities` — Case 5 (2D proximity check)

The prototype's Case 5 (points close to existing segments in all three wire planes) called `get_closest_2d_dis(point)` without specifying APA/face, implicitly assuming a single TPC:

```cpp
// Prototype (single-APA implicit):
auto tuple_results = (*it1)->get_closest_2d_dis(pts.at(i));
```

The toolkit resolves the APA/face per point using the detector geometry and passes them explicitly:

```cpp
// Toolkit (multi-APA correct):
auto wpid = dv->contained_by(pt);           // resolve apa/face for this 3D point
if (wpid.face() == -1 || wpid.apa() == -1) continue;
auto [dist_u, dist_v, dist_w] =
    segment_get_closest_2d_distances(existing_sg, pt, wpid.apa(), wpid.face(), "fit");
```

This ensures the 2D wire-plane distances are measured in the correct TPC, preventing spurious matches across APA boundaries.

### 4.2 `examine_main_vertices` (global) — KNN distance queries

The prototype used `ToyPointCloud::get_closest_dis(point)` (linear scan) for measuring how far a vertex is from the main cluster's Steiner point cloud. The toolkit uses:

```cpp
auto knn = main_cluster->kd_steiner_knn(1, vtx_pt, "steiner_pc");
double closest_dis = knn.empty() ? 1e9 : std::sqrt(knn[0].second);
```

This uses an indexed KD-tree, which is both multi-APA-safe (the Steiner PC already stores the correct 3D positions) and algorithmically faster (O(log N) vs O(N)).

### 4.3 `compare_main_vertices` — fiducial volume check

The toolkit queries the fiducial volume via `cluster.grouping()->get_fiducialutils()`, which is aware of the full multi-APA detector geometry. The prototype's `fid->inside_fiducial_volume(pt, offset_x)` was single-APA.

### 4.4 `find_cont_muon_segment` — vertex position

The toolkit uses the fitted vertex position (`vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point`) as the reference point for computing direction vectors, rather than just the Steiner-graph WCPoint. This is correct in multi-APA because the fitted position may differ from the raw WCPoint and is needed for accurate angle computations.

---

## 5. Per-Function Notes

### 5.1 `examine_main_vertex_candidate`

Returns `tuple<bool flag_in, int ntracks, int nshowers>` for a single candidate vertex.

**Prototype:** Index-based start/end matching: `sg->get_wcpt_vec().front().index == vertex->get_wcpt().index`.  
**Toolkit:** Distance-based: `ray_length(Ray{wcps.front().point, vertex->wcpt().point}) < ray_length(Ray{wcps.back().point, ...})`.

Both correctly identify whether a segment has its head at this vertex and is pointing inward (strong non-weak direction). The Michel electron special case (2 segments, 1 track + 1 shower, small daughter shower count) is reproduced identically.

### 5.2 `determine_main_vertex`

High-level orchestrator. Calls (in order):
1. Scans all vertices for `flag_save_only_showers`
2. If not only-showers and main cluster: calls `improve_vertex` + `fix_maps_shower_in_track_out`
3. Builds `map_vertex_track_shower` and `main_vertex_candidates`
4. For only-showers: calls `compare_main_vertices_all_showers`
5. Otherwise: calls `examine_main_vertices_local` → `compare_main_vertices` (if >1 candidate)
6. Calls `examine_structure_final` (non-shower case only)
7. Calls `examine_direction`

Functionally equivalent to the prototype. Determinism improved via `ordered_nodes`.

**Note:** The prototype stored `map_cluster_main_candidate_vertices[temp_cluster] = main_vertex_candidates` as a global map. The toolkit does not maintain this global map — the candidates are local to the function. This is correct since the map was only used for diagnostic printing in the prototype.

### 5.3 `compare_main_vertices_all_showers`

Used when the cluster contains only showers (no strong tracks). Procedure:

1. Collect all intermediate segment points and vertex positions.
2. Run PCA to find the main axis.
3. Project all vertex candidates onto the axis; find min/max endpoints.
4. Run Dijkstra (toolkit: `do_rough_path`) between the two endpoints using the Steiner point cloud.
5. Fit a temporary segment; call `segment_determine_shower_direction`.
6. Pick the vertex at the end the shower points away from.
7. Override with forward-z vertex if segment > 80 cm and z-separation > 40 cm (large shower case).

**Toolkit improvement:** Manages the temporary graph/segment via RAII (local `shared_ptr`; automatically freed on scope exit), whereas the prototype manually `delete`d three raw pointers. No functional difference but eliminates leak risk.

PCA result is order-independent (computed over the full point set), so the `boost::edges`/`boost::vertices` unordered iteration here does not cause reproducibility issues.

### 5.4 `compare_main_vertices`

Scores each candidate vertex on four criteria and picks the maximum-scoring candidate:

| Criterion | Score delta | Notes |
|---|---|---|
| Proton topology | ±(n_in−n_out)/4 | Rewards outgoing protons |
| Z position | −(z−z_min)/(200 cm) | Prefers upstream (small z) vertices |
| Segment multiplicity and type | +1/4 per track, +1/8 per shower, +1/4 for proton, +1/8 for long muon | More connected particles → better |
| Fiducial volume | +0.5 | Inside FV bonus |
| Topology conflicts | −n_conflicts/4 | Calls `calc_conflict_maps` |

**Toolkit improvement:** Uses fitted vertex position (`vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point`) for z-scoring and fiducial volume check. The prototype used `get_fit_pt()`. After `improve_vertex`, fitted positions are available and more accurate than raw WCPoints.

### 5.5 `calc_conflict_maps`

BFS from the assumed-neutrino vertex, building `map_seg_dir[seg] = {start_vtx, end_vtx}` (tree direction). Then counts conflicts:

1. **Segment direction conflicts:** If a segment's recorded direction disagrees with the tree direction → +1 (strong) or +0.5 (weak).
2. **Vertex angle conflicts:** Large bending at a vertex (in/out angle < 110°) penalized. Backward-going particles penalized.
3. **Multiple incoming:** n_in > 1 → +n_in−1 (or halved if all showers).
4. **Shower-in track-out:** min(n_in_shower, n_out_tracks) penalty.

The `is_shower` check in the toolkit uses all three conditions (kShowerTrajectory ‖ kShowerTopology ‖ pdg==11), matching `get_flag_shower()` in the prototype. The angle computation uses `std::acos + std::clamp` instead of `TVector3::Angle`, which is equivalent.

### 5.6 `examine_direction`

Propagates directions and particle types via BFS from `main_vertex`. Also identifies long-muon chains via `find_cont_muon_segment`, and handles Michel electron 4-momentum calculation.

**Key BFS logic (both prototype and toolkit):**
- If segment has no direction, or weak direction, or is a shower, or `flag_final`: set direction from topology, assign particle type.
- If segment has strong direction: only handle the proton-followed-by-shower-in case (reclassify shower as proton/pion).

**Particle type assignment rules (in order of priority):**
1. If incoming shower exists (`flag_shower_in`): set current segment to electron.
2. If many daughter showers: potentially convert to electron (check back-to-back angle).
3. If already electron but few daughters and no shower topology: potentially convert back to muon.
4. If electron + shower trajectory + 1 daughter: check if daughters dominate → convert to muon.
5. Default for undetermined segments from main vertex: dQ/dx > 1.4 → proton, else muon.

**Toolkit extras:**
- Calls `set_default_shower_particle_info(graph, cluster, ...)` at the end (line 1643) as a final pass ensuring every shower segment has `particle_info` set to electron. This mirrors the prototype's `get_particle_type()` which always returned 11 for showers.
- The Michel electron 4-mom check replaces `get_flag_shower_dQdx()` with `kShowerTrajectory || abs(pdg)==11`. `flag_shower_dQdx` was set when a segment was identified as a shower by dQ/dx — in toolkit, `kShowerTrajectory` covers trajectory-based identification; electron PDG covers prior classification. This is a close but not identical mapping; the practical impact is negligible since both conditions indicate a shower endpoint.

**Long-muon chain detection:**
Both prototype and toolkit call `find_cont_muon_segment` iteratively to find collinear segments forming a long muon, then reclassify them all as muons and add to `segments_in_long_muon`. The condition is `total_length > 45 cm && max_length > 35 cm && num_segments > 1`.

### 5.7 `find_cont_muon_segment`

Searches the segments at `vtx` for a continuation of `sg` that is nearly collinear (angle < 10° at 15 cm, or < 10° at 50 cm for longer segments) and has low dQ/dx (< 1.3 × MIP) unless `flag_ignore_dQ_dx`.

Selects the candidate with maximum `length × cos(angle)` (projected length).

**Functional equivalence:** Angle threshold (10°/15°), dQ/dx threshold (1.3), and projected-length selection are all identical. Toolkit uses `vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point` for the reference point.

**Distinction from `find_cont_muon_segment_nue`:** The `_nue` variant uses 12.5° threshold and checks at 30 cm instead of 50 cm for the long-segment angle. It is used in νe and single-photon taggers, not in the functions reviewed here.

### 5.8 `eliminate_short_vertex_activities`

Iteratively removes stub segments too short to be real, in five cases (in priority order):

| Case | Condition | Removal |
|---|---|---|
| 1 | v1 has 1 seg, v2 has ≥3 segs; direct_length < 0.36 cm (or < 0.5 cm if v2 > 3) | Remove seg + v1 |
| 2 | v2 has 1 seg, v1 has ≥3 segs; direct_length < 0.36 cm (or < 0.5 cm if v1 > 3) | Remove seg + v2 |
| 3 | Segment connected to main_vertex; direct_length < 0.1 cm | Remove seg + non-main vertex |
| 4 | Isolated endpoint within 0.36 cm of another segment (or < 0.45 cm near main_vertex) | Remove seg + endpoint |
| 5 | All points of seg within 0.45 cm of existing segments in all three wire planes | Remove seg (+ endpoint if deg-1) |

**Correctness improvement (Case 3):** The prototype's Case 3 always removed `v2`, relying on the `std::set<ProtoVertex*>` ordering to ensure main_vertex was `v1`. The toolkit explicitly computes `VertexPtr to_remove = (v1 == main_vertex) ? v2 : v1`, making the logic correct regardless of set ordering.

**Multi-APA improvement (Case 5):** See §4.1. The toolkit uses `dv->contained_by(pt)` to resolve APA/face before calling `segment_get_closest_2d_distances(..., apa, face)`.

### 5.9 `change_daughter_type`

Recursively reclassifies downstream segments as a given particle type (typically muon, PDG=13). Applied when a back-to-back track pair is found in `examine_main_vertices_local`.

**Logic:** For each segment at `other_vtx` (the far end of `segment`):
1. Skip if same type, shower trajectory, or strong non-weak direction.
2. If large shower topology (length > 40 cm, dirsign==0): check 170° collinearity → reclassify + unset topology flag → recurse.
3. Fall through to 165° check for length > 10 cm → reclassify → recurse.
4. All other segments: skip (`continue`).

**Note:** The 10 cm check at step 3 is only reached for segments that passed the large-shower-topology check in step 2. Non-topology segments and small topology segments hit the `continue` in the else branch. This matches the prototype's structure exactly.

**4-momentum:** The toolkit sets `particle_info()->set_pdg()` and `set_mass()` but does not immediately call `segment_cal_4mom`. The prototype called `cal_4mom()` if 4-momentum was already set. In practice, `examine_direction` is called after `examine_main_vertices_local`, and it recalculates 4-momenta for all modified segments; so the deferred recalculation is correct.

### 5.10 `examine_main_vertices_local` (prototype: `examine_main_vertices(vertices)`)

Filters the candidate vertex list by removing vertices that appear to be throughgoing tracks (back-to-back pairs). Logic:

1. For vertices with 1 segment: always keep.
2. For vertices with >1 segment: find pairs of segments (length > 10 cm) where the opening angle > 165° and at least one is a muon (PDG 13) with length > 30 cm, OR > 170° and both are protons (PDG 2212) with length > 20 cm.
3. If such back-to-back pairs exist and no other significant segment remains (no strong track > 6 cm, no shower > 35 cm total daughter length): skip the vertex.
   - Reclassify the back-to-back segments as muons via `change_daughter_type`.
   - Follow the muon chain with `find_cont_muon_segment` and add the endpoint to keep.
4. Otherwise: keep the vertex.

**Bug fixed:** `tmp_vertices` changed from `std::set<VertexPtr>` to `IndexedVertexSet` (§2.1).

### 5.11 `examine_main_vertices` (global, NeutrinoPatternBase.cxx)

Operates on all clusters across the event (each with its own candidate main vertex). Two passes:

**Pass 1 — Remove short clusters without strong tracks:**
- `length < cluster_length_cut` (min of 0.6 × main_cluster_length and 6 cm):
  - Keep if any connected segment is not a shower, has direction, and has strong direction or high dQ/dx (> 2 × MIP).
  - Additional cut: length < 5 cm but vertex > 100 cm from main cluster Steiner cloud → also remove.
- `length >= cluster_length_cut` but vertex > 200 cm from main cluster Steiner cloud → remove.

**Pass 2 — If main cluster is all-showers:**
- Compare other cluster vertices against the main cluster direction.
- If angle < 10° and cluster is small (< 15 cm) and close (< 40 cm): remove.
- If angle > 160° and cluster is a large shower (> 0.5 × main cluster length) and close (< 10 cm, angle2 < 25°): swap main cluster.

**Toolkit improvements:**
- Uses `kd_steiner_knn` for distance queries (O(log N) vs O(N)).
- Sorts candidates by `cluster_id` for deterministic evaluation order in pass 2.
- The `swap_main_cluster` function correctly updates `main_cluster` pointer and `other_clusters` list.

---

## 6. Efficiency Notes

| Function | Issue | Severity |
|---|---|---|
| `eliminate_short_vertex_activities` Case 5 | Inner loop checks existence of `existing_sg` in graph via O(N_edges) scan every iteration | Low — existing_segments is typically small |
| `compare_main_vertices_all_showers` | Calls `clustering_points_segments` to associate all cluster points to the temp segment | Acceptable — only done for all-shower clusters with >2 Steiner points |
| `calc_conflict_maps` | `std::map<SegmentPtr, pair<...>>` keyed by shared_ptr (pointer hash) | Acceptable — map is small and built per call |

---

## 7. Summary

| Category | Count | Functions |
|---|---|---|
| Bug fixed | 1 | `examine_main_vertices_local` (`tmp_vertices` non-determinism) |
| Determinism improvements | 5 | `determine_main_vertex`, `examine_direction`, `examine_main_vertices_local`, `examine_main_vertices` (global), `compare_main_vertices` |
| Multi-APA improvements | 4 | `eliminate_short_vertex_activities` (Case 5), `examine_main_vertices` (global), `compare_main_vertices`, `find_cont_muon_segment` |
| Correctness improvements | 3 | `eliminate_short_vertex_activities` (Case 3), `examine_direction` (final shower pass), `compare_main_vertices` (fit position) |
| Functionally equivalent | 4 | `examine_main_vertex_candidate`, `compare_main_vertices_all_showers`, `calc_conflict_maps`, `find_cont_muon_segment`, `change_daughter_type` |
