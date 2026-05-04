# Review: Point Cloud Clustering & Proto-Vertex Initialization

**Reviewed by:** Claude Code  
**Date:** 2026-04-11  
**Branch:** `apply-pointcloud`  
**Prototype refs:**
- `prototype_base/pid/src/NeutrinoID_proto_vertex.h` (find_proto_vertex, init_first_segment, init_point_segment, find_vertices)
- `prototype_base/pid/src/PR3DCluster_point_clustering.h` (clustering_points_master, clustering_points)
- `prototype_base/pid/src/NeutrinoID.cxx` lines 1005-1073 (NeutrinoID::clustering_points)
- `prototype_base/pid/src/PR3DCluster_path.h` lines 288-316 (get_local_extension)  
**Function map:** `clus/docs/porting/neutrino_id_function_map.md` lines 33–68

---

## Implementation status

| Item | Description | Status | File |
|---|---|---|---|
| B1 | `proto_extend_point`: duplicated `has_pc("steiner_pc")` null check | ✅ Applied | `NeutrinoPatternBase.cxx` |
| B2 | `clustering_points_segments`: exact float equality between KD-tree global min and linear-scan main_sg min may differ by ULP | ✅ Applied | `PRSegmentFunctions.cxx` |
| E1 | `clustering_points_segments`: global 2D KD-trees replace O(S×P) inner loop | Already implemented | `PRSegmentFunctions.cxx` |
| E2 | `clustering_points_segments`: eager angle-cache population avoids lazy-init bug on faces not represented by first segment in cache | Already implemented | `PRSegmentFunctions.cxx` |
| D1 | `find_vertices`: distance-based vertex ordering replaces non-deterministic pointer-address iteration | Already implemented | `PRGraph.cxx` |
| D2 | `clustering_points_segments`: `ClusterIdCmp` / `SegmentIndexCmp` comparators for all maps | Already implemented | `PRSegmentFunctions.cxx` |
| I1 | `init_first_segment`: local PCA refinement of boundary endpoints via power iteration | Already implemented | `NeutrinoPatternBase.cxx` |
| I2 | `find_proto_vertex`: `has_segment` survivorship check after all merging/cleanup | Already implemented | `NeutrinoPatternBase.cxx` |
| I3 | `find_proto_vertex`: `main_cluster_initial_pair_vertices` is a local variable, not a member | Already implemented | `NeutrinoPatternBase.cxx` |
| M1 | `clustering_points_segments`: per-(plane,apa,face) 2D KD-trees and per-(apa,face) angle caches | Already implemented | `PRSegmentFunctions.cxx` |

---

## 1. Scope

This review audits the initial pattern recognition pipeline for a cluster: establishing the first trajectory, finding additional track segments, assigning cluster points to segments, and the fallback for small clusters.

### 1.1 Function mapping table

| # | Prototype | Toolkit | Toolkit file:line |
|---|---|---|---|
| 1 | `NeutrinoID::find_proto_vertex(cluster, flag_break, nrounds, flag_back)` | `PatternAlgorithms::find_proto_vertex(graph, cluster, tf, dv, ...)` | `NeutrinoPatternBase.cxx:1453` |
| 2 | `NeutrinoID::init_first_segment(cluster, flag_back)` | `PatternAlgorithms::init_first_segment(graph, cluster, main_cluster, tf, dv, flag_back)` | `NeutrinoPatternBase.cxx:200` |
| 3 | `NeutrinoID::find_vertices(sg)` | `PR::find_vertices(graph, seg)` | `PRGraph.cxx:72` |
| 4 | `PR3DCluster::get_local_extension(wcp, flag)` | `PatternAlgorithms::get_local_extension(cluster, wcp)` | `NeutrinoStructureExaminer.cxx:2343` |
| 5 | `NeutrinoID::clustering_points(cluster)` | `PatternAlgorithms::clustering_points(graph, cluster, dv, ...)` | `NeutrinoTrackShowerSep.cxx:11` |
| 6 | `PR3DCluster::clustering_points_master` + `clustering_points` | `clustering_points_segments(segments, dv, ...)` | `PRSegmentFunctions.cxx:1748` |
| 7 | `NeutrinoID::init_point_segment(cluster)` | `PatternAlgorithms::init_point_segment(graph, cluster, tf, dv)` | `NeutrinoPatternBase.cxx:1614` |

### 1.2 Review goals

All five goals addressed:
1. Functional equivalence with the prototype.
2. Bugs in the toolkit.
3. Algorithmic efficiency.
4. Determinism (pointer-keyed / unordered containers).
5. Multi-APA / multi-face correctness.

### 1.3 Caller context

`TaggerCheckNeutrino.cxx` lines 156-234 drive these functions:

- **Main cluster** (always a "long" cluster): `find_proto_vertex(…, flag_break=true, nrounds=2, flag_back=true)` → `clustering_points` → `separate_track_shower` → …
- **Long other clusters** (>6 cm): same pipeline with `flag_back=false`.
- **Short other clusters** (≤6 cm): `find_proto_vertex(…, flag_break=false, nrounds=1, flag_back=false)` → on failure (`return false`), fall back to `init_point_segment` → `clustering_points` → …

---

## 2. Per-function walk-through

### 2.1 `find_proto_vertex`

**Prototype** (`NeutrinoID_proto_vertex.h:69-194`)  
**Toolkit** (`NeutrinoPatternBase.cxx:1453-1611`)

**Functional equivalence:** The top-level orchestration is identical:
1. Return false if no Steiner PC or fewer than 2 points.
2. Call `init_first_segment`; return false on failure.
3. If main cluster, save initial pair of vertices.
4. If `flag_break_track`: call `break_segments`, then `examine_structure`. Else call `do_multi_tracking`.
5. Run `find_other_segments` for `nrounds_find_other_tracks` rounds.
6. If main cluster, call `examine_structure_3`; if it merges, refit with `do_multi_tracking`.
7. Call `examine_vertices`, `examine_partial_identical_segments`.
8. If main cluster and initial pair saved, call `examine_vertices_3`.
9. Final `do_multi_tracking`.
10. Return true.

**Difference — prototype:** Returns true unconditionally after reaching step 10. Does not check whether any segments survive the merging/cleanup passes.

**Toolkit improvement I2:** Lines 1595-1606 add a survivorship check: scan `ordered_edges(graph)` for at least one segment belonging to this cluster. Returns false if none remain, allowing the caller to fall back to `init_point_segment`. This handles the rare edge case where `examine_structure` merges away the sole segment.

**Difference — prototype side-effect:** `main_cluster_initial_pair_vertices` is a member variable of `NeutrinoID`, set at line 90 as a side effect. Future calls or concurrent use could see stale state.

**Toolkit improvement I3:** `main_cluster_initial_pair_vertices` is a local variable (line 1473) passed explicitly to `examine_vertices_3`. No side effects, no shared state.

**Minor difference:** Prototype `wcpts` size check is inside `if (sg1->get_wcpt_vec().size()>1)` block; everything else (break_segments, examine_structure, etc.) is inside that block. Toolkit has an early `return false` at line 1482-1484 with equivalent semantics.

### 2.2 `init_first_segment`

**Prototype** (`NeutrinoID_proto_vertex.h:416-568`)  
**Toolkit** (`NeutrinoPatternBase.cxx:200-487`)

**Functional equivalence of core logic:**
1. Get two boundary points from the Steiner PC via `get_two_boundary_wcps(2)` (prototype) / `get_two_boundary_steiner_graph_idx` (toolkit).
2. Order endpoints: main cluster by z (ascending or descending per `flag_back_search`); non-main cluster so the point closest to `main_cluster`'s Steiner PC comes first.
3. Run Dijkstra shortest path (Steiner graph); return null if path ≤ 1 point.
4. Create vertices and segment; run single-track fit; assign fit data to vertices/segment.
5. Return null if fine tracking path ≤ 1 point.

**Toolkit improvement I1: local PCA refinement of boundary endpoints (lines 245-341).** The global boundary search picks the extreme Steiner point along a fixed projection axis. When the cluster curves near its end, the true tip may be around a corner. The toolkit collects all Steiner terminal nodes within 10 cm of each boundary candidate, computes their local principal component direction via power iteration on the covariance matrix, then selects the point with maximum projection along that direction. Seeds the power iteration with the global outward direction for sign consistency and fast convergence. Falls back to all Steiner neighbors (not just terminals) if there are fewer than 3 terminals. This is a purely additive improvement; if the cluster is straight near its end, the refined endpoint coincides with the original.

**Difference — prototype uses `main_cluster` member variable directly.** Toolkit takes `main_cluster` as an explicit parameter, adding:
- Sanity check at lines 344-349: if `is_main_flag` (Facade flag) disagrees with pointer comparison, logs a warning and uses pointer comparison as authoritative.
- Fallback at lines 383-390: if `main_cluster == nullptr` for a non-main cluster, applies ascending-z ordering and logs a warning (graceful degradation instead of UB).

**Non-main cluster endpoint ordering:** Prototype calls `main_cluster->get_point_cloud_steiner()->get_closest_dis(tp)`. Toolkit calls `main_cluster->kd_steiner_knn(1, pt, "steiner_pc")` and takes the square root of the returned squared distance. Functionally equivalent: both find the distance from each boundary candidate to the nearest Steiner point of the main cluster and swap if the second point is closer.

### 2.3 `find_vertices`

**Prototype** (`NeutrinoID_proto_vertex.h:3224-3240`)  
**Toolkit** (`PRGraph.cxx:72-101`)

**Core purpose:** Given a segment, return its two incident vertices as an ordered pair `(front, back)`.

**Prototype implementation:** Looks up the segment in `map_segment_vertices` (a `std::map<ProtoSegment*, std::set<ProtoVertex*>>`). Iterates the set and assigns first two elements to `v1`, `v2`. Returns `(null, null)` if the segment is not in the map.

**Non-determinism in prototype:** `std::set<ProtoVertex*>` is ordered by pointer address, which changes between runs. The pair `(v1, v2)` has run-to-run non-deterministic ordering. Every downstream caller that uses "front vertex" for direction/orientation gets different results across runs.

**Toolkit determinism fix D1:** Retrieves `source` and `target` vertices from the Boost graph edge descriptor. Computes the distance from each vertex's `wcpt().point` to the front of `seg->wcpts()` (or `seg->fits()` as fallback). Returns the vertex closer to the front as the first element. This is fully deterministic: same segment, same point cloud, same result every run.

**Fallback:** If both `wcpts()` and `fits()` are empty, returns the pair in graph-native order — deterministic given a fixed graph structure.

**Unrelated but nearby helper `find_other_vertex`** (`PRGraph.cxx:103`): uses the same descriptor-based lookup, no ordering concerns.

### 2.4 `get_local_extension`

**Prototype** (`PR3DCluster_path.h:288-316`)  
**Toolkit** (`NeutrinoStructureExaminer.cxx:2343-2400`)

**Prototype signature:** `WCPoint get_local_extension(WCPoint wcp, int flag)` — `flag=1` uses the regular point cloud, `flag=2` uses the Steiner PC.

**Toolkit signature:** `geo_point_t get_local_extension(Cluster& cluster, const geo_point_t& wcp)` — always uses `"steiner_pc"`.

**Simplification is safe:** The only prototype caller is `examine_vertices_3` (in `NeutrinoID_examine_structure.h`), which always passes `flag=2`. The `flag=1` path was never used in the pattern recognition pipeline. The toolkit correctly omits the parameter.

**Return type:** Prototype returns `WCPoint` (carries 3D + 2D index info). Toolkit returns `geo_point_t`. Callers only use the 3D position, so the simplification is safe.

**Functional equivalence:**
1. Compute local direction via Hough transform at 10 cm radius; negate.
2. Compute angle with drift direction x; if within 7.5° of perpendicular (90°), return original point unchanged.
3. Gather all Steiner points within 10 cm radius.
4. Return the point with maximum projection along the local direction from the query point.

**Toolkit improvement:** Explicit `if (dir1_mag == 0) return wcp` at line 2356 before the angle computation. The prototype relied on `TVector3::Angle()` to handle a zero vector implicitly. The toolkit's explicit check is more defensive.

**Multi-APA note (M2):** Hardcoded `drift_dir(1,0,0)`. The perpendicularity check `|angle - 90| < 7.5` uses the absolute angular distance from 90°, which is sign-agnostic. A track nearly perpendicular to drift is correctly identified regardless of whether drift is +x or -x.

### 2.5 `clustering_points` (dispatcher)

**Prototype** (`NeutrinoID.cxx:1005-1073` + `PR3DCluster_point_clustering.h:1-34`)  
**Toolkit** (`NeutrinoTrackShowerSep.cxx:11-35`)

The toolkit wrapper:
1. Collects all segments belonging to this cluster from the graph using `ordered_edges(graph)` (deterministic iteration).
2. If non-empty, calls `clustering_points_segments(segments, dv, cloud_name, search_range, scaling_2d)`.

The prototype's `NeutrinoID::clustering_points` was more complex: it called `clustering_points_master`, then separately assigned points from `point_sub_cluster_ids` to segments, assigned Steiner points from `point_steiner_sub_cluster_ids`, then built per-segment kdtrees. The toolkit unifies all of this inside `clustering_points_segments`.

**Note on Steiner points:** The prototype assigned both regular and Steiner points to segments via `add_associate_point` and `add_associate_point_steiner`. The toolkit's `create_segment_point_cloud` (called from `clustering_points_segments`) operates on the regular point cloud only (the Voronoi graph is `"basic_pid"`). Steiner points are part of `"steiner_pc"` which is used for structural analysis (Steiner graph, boundary finding) but is separate from the `associate_points` cloud built here. This is a deliberate architectural choice: the output `cloud_name` (default `"associate_points"`) collects regular cluster points near each fitted segment trajectory, not Steiner intermediate nodes.

### 2.6 `clustering_points_segments` (core algorithm)

**Prototype** (`PR3DCluster_point_clustering.h:36-281`)  
**Toolkit** (`PRSegmentFunctions.cxx:1748-2050`)

This is the most complex function. Discussion is organized by sub-topic.

#### Architecture

The prototype used a two-pass design:
1. `PR3DCluster::clustering_points()`: Voronoi assignment → store segment IDs in `point_sub_cluster_ids[]`.
2. `NeutrinoID::clustering_points()`: read IDs → assign points to segments → build kdtrees.

The toolkit unifies both into `clustering_points_segments`, eliminating the intermediate ID array and reducing data movement.

#### Terminal selection (Voronoi seeds)

Both implementations build `map_pindex_segment`: for each segment, find the 5 nearest cluster points to each interior fit point, and assign the first uncontested point as a terminal for that segment. Segments with only 2 fit points use the midpoint instead of interior points.

The toolkit uses `seg->fits()` (the fitted trajectory) for terminal selection, matching the prototype's `get_point_vec()`.

#### Voronoi diagram

Both call Dijkstra multi-source shortest paths from the terminal set using the `paal` library's nearest-recorder visitor, yielding a `nearest_terminal` array. Prototype uses `std::vector<int> terminals`; toolkit converts to `std::vector<vertex_type>` (the Boost graph vertex descriptor type). Functionally identical.

#### Ghost removal — flag logic

Both initialize `flag_change = true` then apply seven conditions that set it to `false` (point belongs to main_sg). A dead-channel check can flip `flag_change` back to `true`. Final action: prototype marks ghost points `point_sub_cluster_ids[i] = -1`; toolkit collects non-ghost points into `map_segment_points[main_sg]`. Logically equivalent.

The seven conditions are identical between prototype and toolkit. The toolkit adds explicit parentheses that were missing in the prototype (lines 196-201 of `PR3DCluster_point_clustering.h` use `||` with `&&` without parens; C++ operator precedence makes `&&` tighter, producing `A || (B && C)` — the toolkit's parenthesised form matches this semantics explicitly, making the intent clear).

#### Ghost removal — 2D distance queries

**Prototype (O(P × S × F)):** For each cluster point `i`, iterates `temp_segments` (all segments across all clusters), calling `get_closest_2d_dis(p)` for each. That method itself scans all fit points of the segment. Total: O(points × segments × fit_points_per_segment).

**Toolkit (O(P × log N)) — efficiency improvement E1:** Builds per-(plane,apa,face) 2D KD-trees from ALL segment fit-point projections upfront (lines 1795-1810). For each cluster point, one KD-tree query per plane gives the global minimum squared 2D distance (`res[0].second`, no sqrt taken). For the main_sg distance, the linear scan lambda `get_2d_dist2_fast` also returns squared distances. Both stay in squared-distance space throughout so the equality comparisons are bit-exact. The 2D threshold checks use `sq_2d_thr = (scaling_2d × search_range)²`, hoisted once before the per-point loop (see §B.2 for the original bug and the fix).

**Cross-cluster scope:** The prototype's `temp_segments` was all segments from ALL clusters (line 162-163), not just the current cluster. The toolkit builds global KD-trees from all input segments (including cross-cluster ones) before processing any cluster (lines 1777-1810). The cross-cluster comparison behavior is preserved.

#### Angle cache — correctness fix E2

The toolkit pre-populates an `ang_cache` keyed by `(apa, face)` before the main loop (lines 1816-1830). The old lazy approach could insert `(0,0,0)` if the first segment encountered for a face lacked DPC data, corrupting all subsequent 2D projections for that face. The eager approach iterates ALL segments and ALL faces in `seg_pts2d`, using the first segment that has non-zero angles for each `(apa, face)` key.

#### Dead channel check

Prototype: `ct_point_cloud.get_closest_dead_chs(p, pind)` — single-APA, single-face.  
Toolkit: `grouping->get_closest_dead_chs(gp, ch_range, apa, face, pind)` — per-(apa,face) lookup. Correctly generalizes to multi-APA geometry (M1).

### 2.7 `init_point_segment`

**Prototype** (`NeutrinoID_proto_vertex.h:395-414`)  
**Toolkit** (`NeutrinoPatternBase.cxx:1614-1653`)

**Purpose:** Fallback for small or degenerate clusters where `find_proto_vertex` fails. Uses the regular point cloud (not Steiner) to find a rough path.

**Functional equivalence:**
1. Get two boundary WCPoints using the regular point cloud (`get_two_boundary_wcps(1)` prototype / `get_two_boundary_wcps(false)` toolkit).
2. Run Dijkstra shortest path; return if path ≤ 1 point.
3. Create two vertices and one segment; add to graph.
4. Run `do_multi_tracking` to fit.

**Difference — prototype graph:** Uses `temp_cluster->graph` (the regular graph built and cached during clustering). Toolkit uses the `"relaxed_pid"` named graph.

**Toolkit improvement:** Line 1620 explicitly calls `cluster.graph_algorithms("relaxed_pid", dv, track_fitter.get_pc_transforms())` to ensure the graph is built and cached before `do_rough_path_reg_pc`. The prototype assumed the graph was already available; missing this could cause a crash if `init_point_segment` were called on a cluster whose graph had not been built.

---

## 3. Bugs

### B1 — Duplicated null check in `proto_extend_point` *(LOW)* ✅ Applied

**File:** `NeutrinoPatternBase.cxx`  
**Original lines 493-495:**
```cpp
if (!cluster.has_pc("steiner_pc")) return {p, 0};

if (!cluster.has_pc("steiner_pc")) return {p, 0};
```
**Problem:** Identical guard appears twice. The second is unreachable dead code; harmless but confusing.  
**Fix:** Remove the second occurrence. Applied in this review.

---

### B2 — Exact float equality between KD-tree global min and linear-scan main_sg min *(MEDIUM)* ✅ Applied

**File:** `PRSegmentFunctions.cxx`  
**Problem:** The ghost-removal logic checks whether main_sg achieves the global 2D minimum per plane via `min_2d_dis == closest_2d_dis`. `min_2d_dis` came from `sqrt(knn_result.second)` (KD-tree) and `closest_2d_dis` from `sqrt(dx²+dy²)` (linear scan). Although both compute the Euclidean distance to the same underlying point when main_sg wins, the intermediate `sqrt` operations run on independently accumulated squared distances and may differ by one ULP, causing the equality to fail spuriously and the point to be discarded.

**Direction of error:** `global_min <= linear_scan_min` always holds (global tree includes main_sg). A spurious mismatch leaves `flag_change = true`, discarding the point. This is a **conservative** error: a correctly-assigned point may be dropped; no mis-assignment occurs.

**Fix:** Removed all `sqrt` calls from both the linear-scan lambda (`get_2d_dist2_fast`, formerly `get_2d_dist_fast`) and the KD-tree query (`res[0].second` used directly). All equality comparisons now operate on squared distances. The 2D threshold checks in conditions 5–7 and in the dead-channel block now compare against `sq_2d_thr = (scaling_2d * search_range)²` (computed once outside the per-point loop). The 3D distance `closest_dis_point.first` retains its `sqrt` since it is used only in `< search_range` threshold checks (not in equality comparisons with another `sqrt` path).

**Renamed:** `get_2d_dist_fast` → `get_2d_dist2_fast`; `closest_2d_dis` → `closest_2d_dis2`; `min_2d_dis` → `min_2d_dis2` to make the squared-distance semantics explicit.

---

## 4. Efficiency

### E1 — Global 2D KD-trees in `clustering_points_segments` *(Already implemented)*

**File:** `PRSegmentFunctions.cxx:1795-1810`  
**Complexity improvement:** O(P × S × F) → O((S × F × pts_per_seg × log N) build + P × 3 × log N query)  
where P = cluster points, S = segments, F = planes (3), N = total fit projections.  
For typical event sizes (P ~ 1000, S ~ 10), this is roughly 30× faster ghost removal.

### E2 — Eager angle-cache population *(Already implemented)*

**File:** `PRSegmentFunctions.cxx:1816-1830`  
Prevents an O(1) lazy-init bug where the angle cache could be poisoned by `(0,0,0)` for a face not represented by the first segment in cache-iteration order.

### E3 — Final has_segment check in `find_proto_vertex`

**File:** `NeutrinoPatternBase.cxx:1597-1601`  
`ordered_edges(graph)` scan is O(E). For typical graph sizes (E ~ 10-50 edges), this is negligible. No optimization needed.

---

## 5. Determinism

| Function | Container / Iteration | Status |
|---|---|---|
| `find_proto_vertex` | `ordered_edges(graph)` in has_segment check | Deterministic |
| `init_first_segment` | Endpoint ordering: z-based (main) or KNN distance (non-main) | Deterministic |
| `find_vertices` | Distance-based ordering (vertex closer to `wcpts().front()` is first) | Deterministic (improved over prototype's pointer-order) |
| `get_local_extension` | `kd_steiner_radius` returns all points in radius; max projection is independent of iteration order | Deterministic |
| `clustering_points` | `ordered_edges(graph)` for segment collection | Deterministic |
| `clustering_points_segments` | `ClusterIdCmp` on `map_cluster_segs`; `SegmentIndexCmp` on `seg_dpc_cache`, `seg_pts3d`, `seg_pts2d` | Deterministic |
| `init_point_segment` | No container iteration beyond path-finding (graph-based, deterministic) | Deterministic |

---

## 6. Multi-APA / Multi-face correctness

### M1 — `clustering_points_segments` per-(plane,apa,face) 2D KD-trees and angle caches

**File:** `PRSegmentFunctions.cxx:1771-1836`  

The toolkit uses `ApFaceKey = std::tuple<int,int,int>` (plane, apa, face) as the key for both the 2D KD-trees and the angle cache. Per-point APA/face is obtained from `clus->wire_plane_id(i)` (line 1940). Dead channel checks use `grouping->get_closest_dead_chs(gp, ch_range, apa, face, pind)` (lines 2006-2016). All three of these correctly generalize the prototype's single-TPC geometry.

The angle cache eagerly fills for all `(apa, face)` combinations that appear in the input segments' DPC data, using the first segment that has non-zero angles for each key (lines 1817-1830). A previously observed bug where a lazily-initialized cache could be poisoned by `(0,0,0)` is eliminated.

### M2 — `get_local_extension` hardcoded drift direction

**File:** `NeutrinoStructureExaminer.cxx:2352`  
The drift direction is hardcoded as `(1,0,0)`. The perpendicularity guard `|angle - 90| < 7.5°` uses the absolute angular distance from 90°, which is invariant under sign flip of the drift vector. A nearly-perpendicular track is correctly identified whether drift is +x or -x. No correctness issue for multi-face detectors.

### M3 — `init_first_segment` z-ordering for main cluster

**File:** `NeutrinoPatternBase.cxx:352-364`  
The main cluster is seeded starting from the high-z (upstream) or low-z (downstream) boundary, controlled by `flag_back_search`. Since z is the beam direction in all supported detector configurations (MicroBooNE, SBND, DUNE-ND, PDHD), this ordering is geometry-independent.

### M4 — `init_point_segment` graph construction

**File:** `NeutrinoPatternBase.cxx:1620`  
The `cluster.graph_algorithms("relaxed_pid", dv, pcts)` call passes `dv` (the detector volumes object), which provides per-APA/face geometry. The resulting graph is APA-aware.

---

## 7. Summary

| Function | Key differences from prototype | Classification |
|---|---|---|
| `find_proto_vertex` | Added `has_segment` survivorship check; local (not member) `main_cluster_initial_pair_vertices` | Improvement |
| `init_first_segment` | PCA refinement of boundary endpoints; explicit `main_cluster` parameter; nullptr fallback | Improvement |
| `find_vertices` | Distance-based vertex ordering replaces non-deterministic pointer-address ordering | Determinism fix (important) |
| `get_local_extension` | Hardcodes steiner PC (flag=2 only, safe); explicit zero-magnitude guard | Simplification + improvement |
| `clustering_points` | Dispatcher pattern; `ordered_edges` for deterministic segment collection | Correct |
| `clustering_points_segments` | One-pass architecture; global 2D KD-trees; eager angle cache; `ClusterIdCmp`/`SegmentIndexCmp`; multi-APA dead channels | Major efficiency + correctness improvement |
| `init_point_segment` | Pre-caches `"relaxed_pid"` graph; uses `do_rough_path_reg_pc` | Improvement |

All identified issues have been resolved.
