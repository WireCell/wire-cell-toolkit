# Review: `apply_pointcloud` — Pattern Recognition Functions

**Reviewed by:** Claude Code  
**Date:** 2026-04-11  
**Branch:** `apply-pointcloud`  
**Prototype refs:**  
  - `prototype_base/pid/src/NeutrinoID_proto_vertex.h` (functions 1–9, 13)  
  - `prototype_base/pid/src/PR3DCluster_nu_proto_vertex.h` (functions 11–12)  
**Function map:** `clus/docs/porting/neutrino_id_function_map.md`

---

## Implementation status

| Item | Description | Status | File |
|---|---|---|---|
| B.1 | `proto_extend_point`: duplicate `has_pc("steiner_pc")` guard | ✅ Already clean | `NeutrinoPatternBase.cxx` |
| B.2 | Prototype `check_end_point` dir_p.x triple-assignment bug | ✅ Fixed in toolkit | `NeutrinoOtherSegments.cxx:852` |
| B.3 | `modify_vertex_isochronous` / `modify_segment_isochronous` / `crawl_segment`: outside-TPC step points silently treated as good | Design choice | `NeutrinoOtherSegments.cxx`, `NeutrinoStructureExaminer.cxx` |
| D.1 | `check_end_point` segment fallback uses `boost::edges(graph)` instead of `ordered_edges` | ✅ Applied | `NeutrinoOtherSegments.cxx:906` |
| D.2 | `crawl_segment`: `map_segment_point` is pointer-ordered `std::map<SegmentPtr,...>` | ✅ Applied | `NeutrinoStructureExaminer.cxx:710` |
| E.1 | `find_other_segments` Step 1 tagging: O(N × S × P) loop, no spatial index | Deferred | `NeutrinoOtherSegments.cxx:56–118` |
| E.2 | `check_end_point`: iterates all vertices and segments per call, no spatial index | Deferred | `NeutrinoOtherSegments.cxx:835` |
| E.3 | `examine_partial_identical_segments`: O(degree² × P) per high-degree vertex | Deferred | `NeutrinoStructureExaminer.cxx:2080` |
| M.1 | `find_other_segments`: multi-APA/face tagging via `wpid_array` | ✅ Verified | `NeutrinoOtherSegments.cxx:61–118` |
| M.2 | `crawl_segment`: multi-APA/face via `dv->contained_by()` per step | ✅ Verified | `NeutrinoStructureExaminer.cxx:792` |
| M.3 | `modify_vertex_isochronous` / `modify_segment_isochronous`: multi-APA/face | ✅ Verified | `NeutrinoOtherSegments.cxx:1105,1253` |

---

## 1. Scope

This review audits 13 functions:

| # | Prototype | Toolkit | File:line |
|---|---|---|---|
| 1 | `break_segments` | `PatternAlgorithms::break_segments` | `NeutrinoPatternBase.cxx:841` |
| 2 | `find_other_segments` | `PatternAlgorithms::find_other_segments` | `NeutrinoOtherSegments.cxx:31` |
| 3 | `check_end_point` | `PatternAlgorithms::check_end_point` | `NeutrinoOtherSegments.cxx:835` |
| 4 | `find_vertex_other_segment` | `PatternAlgorithms::find_vertex_other_segment` | `NeutrinoOtherSegments.cxx:949` |
| 5 | `modify_vertex_isochronous` | `PatternAlgorithms::modify_vertex_isochronous` | `NeutrinoOtherSegments.cxx:1053` |
| 6 | `modify_segment_isochronous` | `PatternAlgorithms::modify_segment_isochronous` | `NeutrinoOtherSegments.cxx:1179` |
| 7 | `examine_segment` | `PatternAlgorithms::examine_segment` | `NeutrinoStructureExaminer.cxx:954` |
| 8 | `crawl_segment` | `PatternAlgorithms::crawl_segment` | `NeutrinoStructureExaminer.cxx:701` |
| 9 | `examine_partial_identical_segments` | `PatternAlgorithms::examine_partial_identical_segments` | `NeutrinoStructureExaminer.cxx:2080` |
| 10 | `find_other_vertex` | `PR::find_other_vertex` | `PRGraph.cxx:103` |
| 11 | `proto_extend_point` | `PatternAlgorithms::proto_extend_point` | `NeutrinoPatternBase.cxx:490` |
| 12 | `proto_break_tracks` | `PatternAlgorithms::proto_break_tracks` | `NeutrinoPatternBase.cxx:594` |
| 13 | `add_proto_connection` | Inlined as `PR::add_segment` | `PRGraph.cxx:28` |

Six review goals:
1. Functional equivalence with the prototype (or demonstrably better).
2. Bugs in the toolkit implementation.
3. Algorithmic efficiency.
4. Determinism (pointer-keyed / unordered containers, iteration order).
5. Multi-APA / multi-face correctness (new in toolkit vs. single-TPC prototype).
6. Helper functions called by these functions.

---

## 2. Per-function walk-through

### 2.1 `break_segments` (NeutrinoPatternBase.cxx:841)

**Purpose:** Find kink points along a segment and split it into two sub-segments at the break point. Iterates at most twice (`count < 2`) over the `remaining_segments` list.

**Functional equivalence:** Faithfully ported. The inner loop structure — `segment_search_kink` → `proto_extend_point` → `proto_break_tracks` — matches the prototype exactly. Start/end vertex identification is now done geometrically (distance to front/back wcpt) rather than by wcpt-index equality, which is more robust to index changes.

**Toolkit improvements over prototype:**
- `merge_nearby_vertices` (line 1051) is called as a post-processing step after the main loop. The prototype omits this. It collapses duplicate vertices that arise from the oscillating-break pattern (two consecutive breaks leaving two vertices at the same position), preventing a degenerate graph topology.
- Debug timing information is emitted via `SPDLOG_LOGGER_DEBUG` instead of unconditional `std::cout`.

**Multi-APA:** N/A for this function; delegates to `do_multi_tracking` and `create_segment_for_cluster` which are multi-APA aware.

**Determinism:** The `remaining_segments` vector is caller-provided; internal loop processing is deterministic once that order is fixed. `merge_nearby_vertices` uses `ordered_nodes(graph)`. ✓

**Efficiency:** No structural issues. The outer `count < 2` cap prevents runaway looping.

---

### 2.2 `find_other_segments` (NeutrinoOtherSegments.cxx:31)

**Purpose:** Discover additional track segments in a cluster that are not yet covered by the existing pattern recognition graph. Uses a Voronoi diagram over steiner terminals, MST, connected-components, and quality cuts to find and add new segments.

**Functional equivalence:** The 9-step algorithm (tag → terminals → Voronoi → terminal graph → MST → components → filter → process → break) is faithfully ported. Deviations are all improvements:

1. **Path extension (lines 511–537):** When `find_vertex_other_segment` returns a vertex whose position differs from the rough-path endpoint by more than 0.1 cm, the toolkit re-routes the segment via `do_rough_path` so that `new_seg->wcpts()` spans the actual vertex positions. Without this, `do_multi_tracking` would use a mis-anchored starting trajectory. Not in prototype.

2. **`v1 == v2` corner case (lines 483–505):** When both endpoint searches return the same vertex (e.g. a very short segment where both ends find the same existing vertex), the toolkit replaces the less-appropriate end with a new vertex. The prototype could silently create a self-loop edge (`add_proto_connection(v1, sg, cluster)` twice with the same vertex). Improvement.

3. **VLA removal (line 231):** Prototype used `int ncounts[num]` (a C99 VLA, non-standard in C++). Toolkit uses `std::vector<int> ncounts(num_components, 0)`. Correct.

**Multi-APA (M.1):** Each steiner point carries a `WirePlaneId wpid` (from `wpid_array`). The tagging loop (line 61) reads `wpid.apa()` and `wpid.face()` to call `transform->backward()` with the correct detector coordinates before querying `get_closest_dead_chs`. The prototype used a single `ct_point_cloud` with hardcoded plane indices 0, 1, 2. This is the primary structural improvement for multi-detector support.

**Determinism:** `find_cluster_segments(graph, cluster)` returns segments in insertion order (deterministic). Component numbering follows `boost::connected_components` on an adjacency_list built from the MST, which is deterministic for a fixed MST. `remaining_segments` is a `std::set<int>` keyed by component index (integer), so ordering is deterministic. ✓

**Efficiency (E.1):** Step 1 tagging (lines 56–118) has complexity O(N × S × P) where N = steiner points (~hundreds to thousands), S = existing segments, P = fit points per segment. This is the dominant cost for large clusters. A k-d tree built on all existing-segment fit points would reduce the per-steiner-point search from O(S × P) to O(log(S × P)). Deferred.

---

### 2.3 `check_end_point` (NeutrinoOtherSegments.cxx:835)

**Purpose:** Given a tracking path, find the existing vertex or segment in the graph whose endpoint is nearest to the front (or back) of the path, using a direction-aware distance test.

**Functional equivalence:** Faithfully ported.

**Bug B.2 — verified fixed:** The prototype (NeutrinoID_proto_vertex.h:1804–1806) contains a copy-paste bug in the backward-direction case:
```cpp
// Prototype (WRONG):
dir_p.x = test_p.x - tracking_path.front().x;
dir_p.x = test_p.y - tracking_path.front().y;  // should be dir_p.y
dir_p.x = test_p.z - tracking_path.front().z;  // should be dir_p.z
```
This sets `dir_p.y` and `dir_p.z` to zero, making the backward direction always point purely in x. The bug also appears in the segment search inner loop (prototype lines 1887–1895). The toolkit fixes this correctly via the `tracking_direction()` helper (line 852), which uses `.x()`, `.y()`, `.z()` setters.

**Determinism (D.1):** The vertex search loop (line 859) uses `ordered_nodes(graph)` — deterministic. However, the segment fallback (line 906) uses:
```cpp
auto [ebegin, eend] = boost::edges(graph);
for (auto eit = ebegin; eit != eend; ++eit) { ... }
```
`boost::edges()` on an `adjacency_list` returns edges in an implementation-defined order that may vary between runs (edge descriptors are invalidated and re-assigned during graph mutations). When multiple segments tie for `min_dis`, the one selected is non-deterministic. **Fix:** replace `boost::edges(graph)` with `ordered_edges(graph)` (which iterates by insertion index) at line 906.

**Multi-APA:** N/A — purely geometric search using Euclidean distances and line distances. No wire-plane checks.

**Efficiency (E.2):** Iterates all vertices (O(V)) in the outer loop and all segments × fit points (O(E × P)) in the inner fallback, called up to three times per `find_vertex_other_segment` invocation. No spatial index. For typical graph sizes (< 100 vertices, < 100 segments) this is not a bottleneck; for pathologically large clusters it could be.

---

### 2.4 `find_vertex_other_segment` (NeutrinoOtherSegments.cxx:949)

**Purpose:** Find or create the vertex where one endpoint of a new segment connects to the existing graph. Calls `check_end_point` up to three times with progressively relaxed thresholds. If a segment is found instead of a vertex, breaks the segment at the closest point and creates a new junction vertex.

**Functional equivalence:** Faithfully ported. The three-call pattern with thresholds (default, 1.2/2.5 cm, 1.5/3.0 cm) is identical to the prototype. When breaking a segment (`proto_break_tracks` path), the toolkit uses `create_segment_for_cluster` + `add_segment` in place of prototype's `new ProtoSegment` + `add_proto_connection`, which is the standard toolkit replacement pattern.

**Multi-APA:** N/A — delegates to `check_end_point` and `create_segment_for_cluster`.

**Determinism:** Inherits D.1 from `check_end_point`. Otherwise deterministic.

---

### 2.5 `modify_vertex_isochronous` (NeutrinoOtherSegments.cxx:1053)

**Purpose:** Handle the case where a new isolated segment's endpoint `v1` lies at approximately the same drift time (x-coordinate) as an existing vertex `vtx`. Projects v1 along the track direction to the x-plane of `vtx`, finds the nearest steiner point, validates connectivity, shifts `vtx` to the new position, and connects the new segment.

**Functional equivalence:** Faithfully ported. Projection formula (lines 1065–1070) matches prototype. Connectivity check with `step_size = 0.6 cm` matches prototype.

**Multi-APA (M.3):** The connectivity scan uses `dv->contained_by(step_p)` per step to get the correct APA/face, then `transform->backward()` before calling `is_good_point`. This is multi-APA correct. The prototype called `ct_point_cloud->is_good_point()` directly without APA/face discrimination.

**Bug B.3 — outside-TPC steps:** Lines 1106–1111:
```cpp
auto test_wpid = dv->contained_by(step_p);
if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
    auto temp_p_raw = transform->backward(step_p, cluster_t0, test_wpid.face(), test_wpid.apa());
    if (!grouping->is_good_point(...)) {
        n_bad++;
    }
}
// Outside-TPC: silently skipped, n_bad unchanged
```
When `dv->contained_by(step_p)` returns `face == -1` (the step point is outside any TPC), the block is skipped entirely and `n_bad` is not incremented. This means outside-TPC steps are treated as "good," making the connectivity check more permissive than the prototype (which would have counted them as bad via `ct_point_cloud->is_good_point`). At TPC boundaries, a straight-line path between two valid points may cross outside the active volume; skipping these steps could allow incorrect isochronous connections. **Fix:** add an `else { n_bad++; }` branch.

---

### 2.6 `modify_segment_isochronous` (NeutrinoOtherSegments.cxx:1179)

**Purpose:** Similar to `modify_vertex_isochronous` but snaps to a point on an existing segment `sg1` rather than to an existing vertex. Finds a fit point on `sg1` that is nearly isochronous with `v1`, projects along the track direction, validates, shifts `v1`, creates two new sub-segments from `sg1`'s endpoints to `v1`, and connects the new segment.

**Functional equivalence:** Faithfully ported. Prototype loops over `sg1->get_point_vec()` (the fine tracking path), toolkit loops over `sg1->fits()` — equivalent. Angle/distance cuts match.

**Multi-APA (M.3):** Same `dv->contained_by()` pattern as `modify_vertex_isochronous`.

**Bug B.3:** Same outside-TPC issue at lines 1253–1258. Same fix applies.

**Note on path update:** Prototype (line 1655) calls `temp_cluster->dijkstra_shortest_paths(vtx_new_wcp, 2)` then `temp_cluster->cal_shortest_path(v2->get_wcpt(), 2)` to compute the new path. Toolkit (line 1275) calls `do_rough_path(cluster, vtx_new_pt, v2->wcpt().point)`. Both use Dijkstra on the steiner graph; the toolkit version is more explicit and passes points by value. Functionally equivalent.

---

### 2.7 `examine_segment` (NeutrinoStructureExaminer.cxx:954)

**Purpose:** Four-step cleanup: (1) find short segments (< 4 cm) with high-degree endpoints where angular analysis suggests a mis-placed vertex, and call `crawl_segment`; (2) merge co-located vertices (< 0.01 cm); (3) remove duplicate segments; (4) remove isolated vertices.

**Functional equivalence:** Faithfully ported.

**Step 1:** Uses `ordered_edges(graph)` (line 958) for deterministic iteration. Angle thresholds (`max_angle > 150`, `min_angle > 105`) match prototype.

**Step 2 (vertex merge):** Prototype merges vertices whose `wcpt().index` is identical. Toolkit merges vertices within 0.01 cm geometric distance (line 1051). The toolkit criterion is more robust — it does not depend on point-cloud index stability across runs or after re-tiling. Functionally equivalent for normal operation.

**Step 3 (duplicate segment detection):** Uses `std::unordered_set<std::string>` keyed by canonical `"min_ptr:max_ptr"` (line 1115). Within a single run this correctly identifies same-endpoint pairs. Cross-run pointer values are not reproducible, but that is immaterial here because the set is local to this call. Correct and efficient (O(1) per segment).

**Step 4:** Uses `ordered_nodes(graph)` for deterministic vertex iteration.

**Determinism:** All iteration uses `ordered_edges` / `ordered_nodes`. ✓

**Multi-APA:** N/A — pure graph topology manipulation.

---

### 2.8 `crawl_segment` (NeutrinoStructureExaminer.cxx:701)

**Purpose:** "Crawl" a segment inward from a vertex endpoint, testing each fit point for good connectivity to other segments at that vertex. If a farther-inward point passes all connectivity checks, shorten the segment to that point and update all connected segments to reach the new vertex position.

**Functional equivalence:** Faithfully ported. The key algorithm — test points from vertex end inward, check each against reference points at ~3 cm on other segments, break on first bad step — matches the prototype.

**Multi-APA (M.2):** The good-point check (line 792) uses `dv->contained_by(test_p)` to get the correct APA/face:
```cpp
auto test_wpid = dv->contained_by(test_p);
if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
    auto temp_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
    if (!cluster.grouping()->is_good_point(...)) {
        n_bad++;
        if (n_bad > 0) break;
    }
}
```
The prototype called `ct_point_cloud->is_good_point(test_p, ...)` directly. Multi-APA improvement confirmed.

**Bug B.3:** Same outside-TPC issue (line 793 guard). Outside-TPC step points increment neither `n_bad` nor any counter — they are silently skipped. Same fix as §2.5.

**Determinism (D.2):** The `map_segment_point` container (line 710) is:
```cpp
std::map<SegmentPtr, Facade::geo_point_t> map_segment_point;
```
`SegmentPtr` is `std::shared_ptr<Segment>`. The default `operator<` for `shared_ptr` compares stored raw pointers, which are non-deterministic across process invocations (ASLR). The iteration at line 779 walks this map in pointer order. Different orderings determine which reference segment's connectivity check breaks first, potentially giving different `max_bin` values for the same input geometry.

**Fix:** Replace with `std::vector<std::pair<SegmentPtr, geo_point_t>>` populated in the `boost::out_edges` loop, then sorted by `graph[edge].index` (insertion order):
```cpp
std::vector<std::pair<SegmentPtr, Facade::geo_point_t>> seg_points;
// ... populate ...
std::sort(seg_points.begin(), seg_points.end(),
    [&](const auto& a, const auto& b) {
        return segment_graph_index(graph, a.first) < segment_graph_index(graph, b.first);
    });
```

**Path update:** The prototype builds the updated path for other segments by stepping through the steiner cloud at 1 cm intervals (lines 2296–2306), snapping each step to the nearest steiner point. The toolkit (line 907) calls `do_rough_path(cluster, vtx_new_point, min_wcpt_point)`, which is equivalent (also uses Dijkstra on the steiner graph). The toolkit version is cleaner and more correct (optimal path vs. greedy stepping).

---

### 2.9 `examine_partial_identical_segments` (NeutrinoStructureExaminer.cxx:2080)

**Purpose:** For vertices with degree > 2, detect pairs of outgoing segments that overlap by > 5 cm (fit points within 0.3 cm of each other). If such overlap is found, either merge to an existing nearby vertex (< 0.3 cm) or create a new split vertex at the overlap endpoint.

**Functional equivalence:** Faithfully ported with improved determinism.

**Prototype issues addressed:**
- Iterates `map_vertex_segments` (pointer-ordered) → toolkit uses `ordered_nodes(graph)` (line 2089). ✓
- Prototype iterates `it->second` (pointer-ordered `ProtoSegmentSet`) → toolkit sorts out-edges by `graph[a].index < graph[b].index` (lines 2108–2111). ✓

**Path update for reconnected segments:** Toolkit uses `do_rough_path` + `create_segment_point_cloud` + `remove_segment` + `add_segment` to reconnect `max_sg1` and `max_sg2`. Prototype used `del_proto_connection` + `cal_shortest_path` + `add_proto_connection`. Both recompute paths via Dijkstra on the steiner graph. Equivalent.

**Bridge segment creation:** When creating a new vertex (`vtx2`), the toolkit creates a segment between `vtx2` and `vtx` via `do_rough_path` (line 2277). Prototype does the same via `dijkstra_shortest_paths` + `cal_shortest_path`. Equivalent.

**Multi-APA:** Relies on `do_rough_path` and `create_segment_for_cluster` which are multi-APA aware. No direct wire-plane checks in this function.

**Efficiency (E.3):** Outer loop: O(V) vertices. For each vertex with degree > 2: O(degree² × P) where P = fit points per segment. `segment_get_closest_point` (line 2145) is called per fit point of `sg1` against `sg2`. For high-degree vertices (e.g. the neutrino vertex) with long segments, this is expensive. A bounding-box pre-filter on segment pairs could reduce the constant factor.

---

### 2.10 `find_other_vertex` (PRGraph.cxx:103)

**Purpose:** Given a segment (graph edge) and one of its endpoint vertices, return the other endpoint.

**Functional equivalence:** Correct. Prototype retrieved vertices from `map_segment_vertices[sg]` (a `ProtoVertexSet`). Toolkit uses `boost::source`/`target` on the edge descriptor — more direct.

**Additional safety check:** Line 114 calls `boost::edge(vd1, vd2, graph)` to verify the edge still exists in the graph before deriving the other vertex. Not present in prototype. Minor defensive improvement.

**Determinism / multi-APA:** N/A.

---

### 2.11 `proto_extend_point` (NeutrinoPatternBase.cxx:490)

**Purpose:** Extend a point forward along a given direction through the steiner point cloud, advancing 1 cm per step with a momentum-weighted direction update. Returns the final steiner point and its index.

**Functional equivalence:** Faithfully ported. Forward step sizes (1, 2, 3 cm per iteration), angle thresholds (< 25° for steiner, < 17.5° for regular cloud), direction update formula (`dir = dir2.norm() + dir * 5.0; dir = dir / dir.magnitude()`) all match the prototype.

**Bug B.1:** Lines 493–495 contain a duplicate check:
```cpp
if (!cluster.has_pc("steiner_pc")) return {p, 0};
if (!cluster.has_pc("steiner_pc")) return {p, 0};  // duplicate
```
The second guard is dead code. **Fix:** remove line 495.

**Prototype difference — `dir` units:** In the prototype (PR3DCluster_nu_proto_vertex.h:77), the momentum trick is:
```cpp
dir = dir2 + dir * 5 * units::cm;  // dir2 has units cm, dir has units cm, result has units cm²/cm = cm
dir = dir.Unit();
```
The prototype mixes `TVector3` (magnitude-aware) with spatial units in a dimension-inconsistent way, but `Unit()` normalizes it. In the toolkit (line 545):
```cpp
dir = dir2.norm() + dir * 5.0;  // both dimensionless unit vectors, result is dimensionless
dir = dir / dir.magnitude();
```
The toolkit uses dimensionless unit vectors throughout. The `5.0` scalar weight matches the prototype's `5 * units::cm / units::cm = 5.0`. Functionally equivalent.

**Determinism / multi-APA:** N/A — pure steiner kNN lookup, deterministic for k=1.

---

### 2.12 `proto_break_tracks` (NeutrinoPatternBase.cxx:594)

**Purpose:** Given start, break, and end points, compute two Dijkstra paths (start→break, break→end) and remove overlapping points at the junction. Returns false if either sub-path has ≤ 1 point or both distances < 1 cm.

**Functional equivalence:** Faithfully ported.

**Overlap removal:** Prototype (PR3DCluster_nu_proto_vertex.h:23–36) checks `(*it).index == (*it1).index` (wcpt index equality). Toolkit (lines 636–657) uses geometric distance < 0.01 cm. Both correctly identify the junction overlap; the toolkit version does not depend on point-cloud index stability.

**`do_rough_path` as replacement for `dijkstra_shortest_paths` + `cal_shortest_path`:** The prototype called these as two separate member functions. The toolkit wraps both into `do_rough_path(cluster, from, to)`. Functionally equivalent.

**Determinism / multi-APA:** `do_rough_path` uses Dijkstra on the steiner graph (deterministic). N/A for multi-APA.

---

### 2.13 `add_proto_connection` → `PR::add_segment` (PRGraph.cxx:28)

**Purpose:** Register a segment in the graph between two vertices, maintaining all bookkeeping (vertex cluster assignment, insertion-ordered index).

**Functional equivalence:** The prototype's `add_proto_connection` maintained 6 maps (`map_vertex_cluster`, `map_cluster_vertices`, `map_segment_cluster`, `map_cluster_segments`, `map_vertex_segments`, `map_segment_vertices`). The toolkit Boost graph edge replaces all of these. `add_segment` additionally calls `add_vertex` for each endpoint (no-op if the descriptor is already valid), sets `graph[ed].segment = seg` and `graph[ed].index = num_edge_indices++`. The monotonically increasing edge index is the mechanism for deterministic `ordered_edges(graph)` iteration throughout the codebase.

**Determinism:** Edge indices assigned in call order → deterministic graph iteration. ✓

---

## 3. Bugs

### B.1 — Duplicate steiner-pc guard in `proto_extend_point`

**Status:** ✅ Already clean  
**Location:** `NeutrinoPatternBase.cxx:493`

The duplicate guard observed during review was already removed before the fixes in this session. Only one `has_pc("steiner_pc")` check is present.

---

### B.2 — Prototype `check_end_point` dir_p.x triple-assignment bug (verified fixed)

**Status:** ✅ Fixed in toolkit  
**Location:** `NeutrinoOtherSegments.cxx:852` (via `tracking_direction()` helper)  
**Prototype location:** `NeutrinoID_proto_vertex.h:1804–1806` and `1892–1895`

**Problem (prototype only):** In the backward-direction case, the prototype assigns to `dir_p.x` three times:
```cpp
dir_p.x = test_p.x - tracking_path.front().x;
dir_p.x = test_p.y - tracking_path.front().y;  // BUG: should be dir_p.y
dir_p.x = test_p.z - tracking_path.front().z;  // BUG: should be dir_p.z
```
This leaves `dir_p.y = 0` and `dir_p.z = 0`, making the backward direction always point purely in x. The bug appears in both the vertex search loop and the segment search inner loop.

**Toolkit fix:** The `tracking_direction()` free function (lines 798–829 of NeutrinoOtherSegments.cxx) correctly sets all three components using `.x()`, `.y()`, `.z()` setters.

---

### B.3 — Outside-TPC step points silently treated as good

**Status:** Design choice  
**Location:**
  - `NeutrinoOtherSegments.cxx:1106–1111` (`modify_vertex_isochronous`)
  - `NeutrinoOtherSegments.cxx:1253–1258` (`modify_segment_isochronous`)
  - `NeutrinoStructureExaminer.cxx:793–799` (`crawl_segment`)

**Observation:** All three functions check `dv->contained_by(step_p)` and only evaluate `is_good_point` if `face != -1 && apa != -1`. Points outside any TPC are silently skipped (not counted as bad). The prototype's `ct_point_cloud->is_good_point()` would return `false` for out-of-TPC points in a single-TPC setup. In a multi-TPC detector, a step crossing between two active volumes naturally lands outside both TPCs, so counting it as bad would incorrectly veto valid cross-TPC connectivity. Treating outside-TPC steps as neutral (neither good nor bad) is therefore the correct behaviour for multi-TPC detectors. No fix needed.

---

## 4. Determinism

### D.1 — `check_end_point` segment search uses non-ordered edge iteration

**Status:** ✅ Applied  
**Location:** `NeutrinoOtherSegments.cxx:906`

**Problem:** The original code used `boost::edges(graph)` whose iteration order changes after `remove_segment` / `add_segment` operations, giving non-deterministic `min_sg` selection when two segments tie for minimum distance.

**Fix applied:**
```cpp
// Before:
auto [ebegin, eend] = boost::edges(graph);
for (auto eit = ebegin; eit != eend; ++eit) {
    SegmentPtr candidate_seg = graph[*eit].segment;

// After:
for (const auto& ed : ordered_edges(graph)) {
    SegmentPtr candidate_seg = graph[ed].segment;
```

---

### D.2 — `crawl_segment` uses pointer-ordered map for segment → reference point

**Status:** ✅ Applied  
**Location:** `NeutrinoStructureExaminer.cxx:710`

**Problem:**
```cpp
std::map<SegmentPtr, Facade::geo_point_t> map_segment_point;
```
`SegmentPtr = std::shared_ptr<Segment>`. The default `std::map` key comparator uses `std::less<shared_ptr<T>>`, which compares raw stored pointers. Pointer values depend on heap allocation order (affected by ASLR and prior allocations). Iterating `map_segment_point` at line 779 tested reference segments in an order that varied across runs, potentially yielding different `max_bin` results for the same geometry.

**Fix applied:** Replaced the pointer-ordered map with an edge-index-sorted vector. During the `boost::out_edges` loop, each `(edge_index, segment, ref_point)` triple is collected into a temporary vector, sorted by `graph[ed].index` (monotonically assigned at edge creation time), and then the ordered pairs are stored in `seg_ref_points`:
```cpp
std::vector<std::pair<SegmentPtr, Facade::geo_point_t>> seg_ref_points;
{
    std::vector<std::tuple<int, SegmentPtr, Facade::geo_point_t>> tmp;
    auto edge_range = boost::out_edges(vd, graph);
    for (auto eit = edge_range.first; eit != edge_range.second; ++eit) {
        SegmentPtr sg = graph[*eit].segment;
        if (!sg || sg == seg) continue;
        const auto& fits = sg->fits();
        if (fits.empty()) continue;
        Facade::geo_point_t min_point = fits.front().point;
        double min_dis = 1e9;
        for (size_t i = 0; i < fits.size(); i++) {
            double dis = std::fabs(ray_length(Ray{fits[i].point, vertex->fit().point}) - 3.0 * units::cm);
            if (dis < min_dis) { min_dis = dis; min_point = fits[i].point; }
        }
        tmp.emplace_back(graph[*eit].index, sg, min_point);
    }
    std::sort(tmp.begin(), tmp.end(),
              [](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });
    for (auto& [idx, sg, pt] : tmp) {
        seg_ref_points.emplace_back(sg, pt);
    }
}
```
All three downstream loops that previously iterated `map_segment_point` now iterate `seg_ref_points` with structured bindings (`for (const auto& [ref_sg, ref_pt] : seg_ref_points)`), giving deterministic traversal order tied to graph insertion sequence.

---

## 5. Efficiency

### E.1 — `find_other_segments` Step 1: O(N × S × P) tagging loop

**Location:** `NeutrinoOtherSegments.cxx:56–118`

**Current cost:** For each of N steiner points, iterates over S existing segments, and for each segment calls `segment_get_closest_point` (iterates P fit points) and `segment_get_closest_2d_distances`. Total: O(N × S × P). For a cluster with 1000 steiner points, 5 segments, and 100 fit points each, this is 500,000 distance evaluations per tagging pass.

**Improvement:** Build a k-d tree over all existing-segment fit points once, then for each steiner point query the nearest fit point in O(log(S × P)). Reduces total to O(N × log(S × P) + S × P × build). This would be a significant speedup for clusters with many segments. Deferred.

---

### E.2 — `check_end_point`: full graph scan per call

**Location:** `NeutrinoOtherSegments.cxx:835`

**Current cost:** Iterates all V vertices (O(V)) in the outer `ncount=5` direction loop, then if no vertex found, iterates all E edges × P fit points (O(E × P)). Called up to 3 times per `find_vertex_other_segment` invocation, which is called twice per candidate segment in `find_other_segments`.

**Improvement:** A spatial hash or k-d tree on vertex positions would reduce the vertex search to O(log V). A bounding-box index on segments would reduce the segment search. Deferred given typical graph sizes (< 100 vertices).

---

### E.3 — `examine_partial_identical_segments`: O(degree² × P) per vertex

**Location:** `NeutrinoStructureExaminer.cxx:2113–2161`

**Current cost:** For each vertex with degree D: O(D² × P) where P = fit points per segment. For the neutrino vertex with D ≥ 5 and segments with 50+ fit points, this is ~3000+ `segment_get_closest_point` calls per outer iteration.

**Improvement:** Pre-compute pairwise segment bounding-box overlaps as a pre-filter; only compare detailed fit points for overlapping bounding-box pairs. Deferred.

---

## 6. Multi-APA / Multi-Face Handling

| Function | Multi-APA aware | Mechanism | Notes |
|---|---|---|---|
| `find_other_segments` | ✅ Yes | `wpid_array[i]` per steiner point; `transform->backward(p, t0, face, apa)` | M.1: each steiner point carries its APA/face |
| `crawl_segment` | ✅ Yes | `dv->contained_by(test_p)` per step point | M.2: see B.3 for outside-TPC issue |
| `modify_vertex_isochronous` | ✅ Yes | `dv->contained_by(step_p)` per step | M.3: see B.3 |
| `modify_segment_isochronous` | ✅ Yes | `dv->contained_by(step_p)` per step | M.3: see B.3 |
| `break_segments` | Delegates | Via `do_multi_tracking`, `create_segment_for_cluster` | Those functions are multi-APA aware |
| `examine_segment` | N/A | Pure graph topology | No wire-plane checks |
| `check_end_point` | N/A | Geometric distance only | No wire-plane checks |
| `find_vertex_other_segment` | N/A | Delegates to `check_end_point` and `create_segment_for_cluster` | |
| `examine_partial_identical_segments` | Delegates | Via `do_rough_path`, `create_segment_for_cluster` | |
| `proto_extend_point` | N/A | Steiner kNN only | |
| `proto_break_tracks` | N/A | `do_rough_path` on steiner graph | |
| `find_other_vertex` | N/A | Pure graph lookup | |
| `add_proto_connection` (inlined) | N/A | Pure graph edge addition | |

The four functions that do wire-plane / dead-channel checks are all multi-APA correct in their check mechanisms. The common issue (B.3) is that outside-TPC points are silently skipped rather than counted as bad.

---

## 7. Helper Functions

Key helpers called by the 13 functions reviewed above:

| Helper | Purpose | Notes |
|---|---|---|
| `tracking_direction()` | Compute direction from tracking path at index `i` | Free function in anonymous namespace; fixes prototype's dir_p.x bug (B.2) |
| `segment_search_kink()` | Find a kink point in a segment's fit points | Not reviewed here; called by `break_segments` |
| `segment_get_closest_point(seg, p, "fit")` | Nearest fit point on segment to query point | O(P); used extensively |
| `segment_get_closest_2d_distances(seg, p, apa, face, "fit")` | 2D projected distances on wire planes | Multi-APA aware via apa/face parameters |
| `segment_cal_dir_3vector(seg, pt, range)` | Direction vector at a point along a segment, up to `range` cm | Used by `modify_*_isochronous` and `examine_segment` |
| `do_rough_path(cluster, from, to)` | Dijkstra shortest path on steiner graph | Deterministic; used widely for path computation |
| `create_segment_for_cluster(cluster, dv, pts, dirsign)` | Create a segment with wcpts, "main" point cloud | Multi-APA aware (uses `dv`) |
| `create_segment_point_cloud(seg, pts, dv, "main")` | Rebuild "main" point cloud on existing segment | Used when updating wcpts in-place |
| `merge_vertex_into_another(graph, from, into, dv)` | Re-parent all edges from `from` to `into`, remove `from` | Used in `examine_segment` Step 2 |
| `merge_nearby_vertices(graph, cluster, track_fitter, dv)` | Collapse co-located vertices (< 0.1 cm) | Toolkit-only addition; called by `break_segments` |
| `ordered_nodes(graph)` | Iterate graph vertices in insertion order | Provides deterministic ordering; used throughout |
| `ordered_edges(graph)` | Iterate graph edges in insertion order | Provides deterministic ordering; D.1 fix requires using this |
| `find_vertices(graph, seg)` | Return both endpoint vertices of a segment | Used pervasively |
| `find_segment(graph, v1, v2)` | Find segment connecting two vertices, if any | Used in `examine_partial_identical_segments` |
| `remove_segment(graph, seg)` | Remove edge from graph | Decrements edge count; does not invalidate ordered iteration |
| `remove_vertex(graph, vtx)` | Remove vertex from graph | |

---

## 8. Summary

### Functional correctness

The 13 functions are faithfully ported from the prototype. The toolkit is in all cases equivalent to or better than the prototype:

1. **B.2 prototype bug fixed:** `check_end_point` backward direction now correctly computes y and z components.
2. **VLA removed:** `find_other_segments` uses `std::vector<int>` instead of non-standard `int ncounts[num]`.
3. **Corner cases added:** `v1 == v2` handler in `find_other_segments` prevents self-loop edges.
4. **Path extension:** `find_other_segments` re-anchors segments to actual vertex positions before fitting.
5. **Post-processing:** `break_segments` calls `merge_nearby_vertices` to fix oscillating-break residuals.
6. **Graph robustness:** `find_other_vertex` validates edge existence before returning.

### Applied fixes

| Item | Function(s) | Change | Status |
|---|---|---|---|
| D.1 | `check_end_point` | `boost::edges(graph)` → `ordered_edges(graph)` in segment fallback search | ✅ Applied (`NeutrinoOtherSegments.cxx:906`) |
| D.2 | `crawl_segment` | `std::map<SegmentPtr,geo_point_t>` → edge-index-sorted `vector<pair<...>>` | ✅ Applied (`NeutrinoStructureExaminer.cxx:710`) |

### Deferred / design choices

| Item | Function(s) | Notes |
|---|---|---|
| B.1 | `proto_extend_point` | Duplicate `has_pc` guard was already absent in toolkit — no action needed |
| B.3 | `modify_vertex_isochronous`, `modify_segment_isochronous`, `crawl_segment` | Outside-TPC neutral treatment is correct for multi-TPC detectors — design choice |
| E.1 | `find_other_segments` | O(N × S × P) tagging loop — deferred, spatial index improvement possible |
| E.2 | `check_end_point` | Full graph scan per call — deferred, graph sizes are small in practice |
| E.3 | `examine_partial_identical_segments` | O(degree² × P) per vertex — deferred, bounding-box pre-filter possible |

### Deferred (efficiency)

E.1, E.2, E.3 are noted but not blocking. Spatial indexing for `find_other_segments` (E.1) would give the largest speedup for large clusters.

### Multi-APA: confirmed correct

All four functions performing wire-plane checks (`find_other_segments`, `crawl_segment`, `modify_vertex_isochronous`, `modify_segment_isochronous`) correctly use per-point APA/face information from `wpid_array` and `dv->contained_by()`. The outside-TPC issue (B.3) is an edge case at detector boundaries, not a fundamental APA-awareness failure.
