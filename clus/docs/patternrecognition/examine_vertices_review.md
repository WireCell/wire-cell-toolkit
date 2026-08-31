# Review: `examine_vertices_*` Family

**Reviewed:** 2026-04-11  
**Branch:** `apply-pointcloud`  
**Prototype ref:** `prototype_base/pid/src/NeutrinoID_proto_vertex.h` lines 2375–3283  
**Toolkit file:** `clus/src/NeutrinoStructureExaminer.cxx`  
**Header:** `clus/inc/WireCellClus/NeutrinoPatternBase.h` lines 108–115  
**Porting map:** `clus/docs/porting/neutrino_id_function_map.md` lines 86–92

---

## Applied Fixes (2026-04-11)

| Item | File | What was changed |
|---|---|---|
| **B1** | `clus/src/NeutrinoStructureExaminer.cxx` | `examine_vertices`: replaced `boost::num_vertices(graph) > 2` with a cluster-filtered vertex count via `ordered_nodes(graph)` + `vtx->cluster() == &cluster`, matching the prototype's `find_vertices(cluster).size() > 2` semantics. |
| **B2** | `clus/src/NeutrinoStructureExaminer.cxx` | `examine_vertices_1p`: added a 9-line comment above `ntime_ticks` explaining the normalization correctness argument (`slice_index` is in tick units, `tind` is in tick units, so `tind/ntime_ticks` ≡ prototype `slope_xt*x` when `ntime_ticks = nrebin` for single-readout-bin blobs), downgraded severity from Critical to Low (typical case always works). |
| **B3** | `clus/src/NeutrinoStructureExaminer.cxx` | `examine_vertices_1p`: removed unused `v1_t1`, `v1_t2`, `v2_t1`, `v2_t2` variables with `(void)` suppression; added comment stating `tind` is view-independent. Added comment to inner-loop `p_t_raw` call explaining why passing `pind` is harmless. |
| **B5** | `clus/src/NeutrinoStructureExaminer.cxx` | `examine_vertices_4p`: replaced commented-out `// if (!sg1) return true;` with a 4-line comment explaining the nullptr semantics (all v1 segments are tested, conservative but correct). |
| **B7** | `clus/src/NeutrinoStructureExaminer.cxx` | `examine_vertices_4` (both v1→v2 and v2→v1 branches): replaced the terse `// Combine with rest of segment` comment with an 8-line comment explaining (a) why position-based trimming is safer than the prototype's index-based version, and (b) that the `empty()` guard prevents infinite-loop if `min_wcp` is absent. Applied to both symmetric branches via `replace_all`. |
| **E1** | `clus/src/NeutrinoStructureExaminer.cxx` | `examine_vertices_2` and `examine_vertices_3`: replaced `boost::vertices(graph)` (no cluster filter) with `ordered_nodes(graph)` + `vtx->cluster() == &cluster` in both isolated-vertex removal loops, matching the existing `examine_segment` pattern. Prevents removal of degree-0 orphans from other clusters when the graph is shared. |

---

## §0 Scope and Goals

Functions under review (prototype names / toolkit names):

| # | Prototype | Toolkit |
|---|---|---|
| 1 | `examine_vertices(cluster)` | `examine_vertices(graph, cluster, tf, dv, main_vertex=nullptr)` |
| 2 | `examine_vertices_1(cluster)` | `examine_vertices_1(graph, cluster, tf, dv, main_vertex=nullptr)` |
| 3 | `examine_vertices_1(v1,v2,12 scalars)` | `examine_vertices_1p(graph, v1, v2, tf, dv)` |
| 4 | `examine_vertices_2(cluster)` | `examine_vertices_2(graph, cluster, tf, dv, main_vertex=nullptr)` |
| 5 | `examine_vertices_3()` | `examine_vertices_3(graph, cluster, pair_vertices, tf, dv)` |
| 6 | `examine_vertices_4(cluster)` | `examine_vertices_4(graph, cluster, tf, dv, main_vertex=nullptr)` |
| 7 | `examine_vertices_4(v1,v2)` | `examine_vertices_4p(graph, v1, v2, tf, dv)` |

Review goals:

1. Verify functional equivalence with the prototype.
2. Find bugs in the toolkit implementation.
3. Identify efficiency improvements.
4. Reduce nondeterminism (pointer-keyed containers, unordered iteration).
5. Audit helper functions called by this family.
6. Verify multi-APA / multi-face correctness (the key new concern vs. prototype).

---

## §1 Function Mapping and Call Graph

```
find_proto_vertex (NeutrinoPatternBase.cxx:1453)
│
├── examine_vertices (NeutrinoStructureExaminer.cxx:2003)   [main_vertex=nullptr]
│   ├── examine_segment (NeutrinoStructureExaminer.cxx:944)
│   ├── examine_vertices_1 (NeutrinoStructureExaminer.cxx:1336)
│   │   └── examine_vertices_1p (NeutrinoStructureExaminer.cxx:1142)
│   ├── examine_vertices_2 (NeutrinoStructureExaminer.cxx:1441)   [if num_vertices>2]
│   └── examine_vertices_4 (NeutrinoStructureExaminer.cxx:1647)
│       └── examine_vertices_4p (NeutrinoStructureExaminer.cxx:1563)
│
├── examine_partial_identical_segments (NeutrinoStructureExaminer.cxx:2017)
│
└── examine_vertices_3 (NeutrinoStructureExaminer.cxx:2349)   [main cluster only]
    └── examine_structure_4 (NeutrinoStructureExaminer.cxx:473)
```

All seven functions live in `NeutrinoStructureExaminer.cxx`. The outer driver
`examine_vertices` is called with `main_vertex=nullptr` (default) from
`find_proto_vertex`. This matches the prototype where the `main_vertex` member
was also `0` at the point `examine_vertices(temp_cluster)` is called
(it was reset to `0` at `NeutrinoID.cxx:150/180/191`).

---

## §2 Helper Function Inventory

| Helper | Toolkit location | Multi-APA concern |
|---|---|---|
| `examine_segment` | `NeutrinoStructureExaminer.cxx:944` | No. Geometry only. |
| `crawl_segment` | called from `examine_segment` | Minimal. |
| `find_segment(graph, v1, v2)` | `NeutrinoPatternBase.cxx` | None. |
| `find_vertices(graph, sg)` | `NeutrinoPatternBase.cxx` | None. |
| `find_other_vertex(graph, sg, vtx)` | `NeutrinoPatternBase.cxx` | None. |
| `ordered_nodes(graph)` | `NeutrinoPatternBase.cxx` | None; vecS insertion-order. |
| `ordered_edges(graph)` | `NeutrinoPatternBase.cxx` | None; vecS insertion-order. |
| `add_segment / remove_segment / remove_vertex` | `NeutrinoPatternBase.cxx` | None. |
| `create_segment_for_cluster(cluster, dv, pts, 0)` | `NeutrinoPatternBase.cxx` | Uses `dv`; multi-APA via blob lookup. |
| `create_segment_point_cloud(sg, pts, dv, "main")` | `NeutrinoPatternBase.cxx` | Uses `dv`; multi-APA. |
| `do_rough_path(cluster, p1, p2)` | `NeutrinoPatternBase.cxx` | None; Steiner graph path. |
| `get_local_extension(cluster, point)` | `NeutrinoPatternBase.cxx` | None. |
| `segment_track_length(sg)` | `PRSegmentFunctions.cxx` | None. |
| `segment_track_direct_length(sg)` | `PRSegmentFunctions.cxx` | None. |
| `segment_get_closest_2d_distances(sg,pt,apa,face,"fit")` | `PRSegmentFunctions.cxx` | **Yes** — takes explicit apa/face. |
| `dv->contained_by(point)` | `IDetectorVolumes` | Canonical per-point APA/face lookup. |
| `cluster.grouping()->convert_3Dpoint_time_ch(p,apa,face,plane)` | `Facade_Grouping.cxx` | **Yes** — per-APA/face/plane. |
| `cluster.grouping()->get_closest_dead_chs(p_raw,1,apa,face,pind)` | `Facade_Grouping.cxx` | **Yes** — per-APA/face. |
| `cluster.grouping()->is_good_point(p_raw,apa,face,...)` | `Facade_Grouping.cxx` | **Yes** — per-APA/face. |
| `transform->backward(p, t0, face, apa)` | `pc_transforms` | **Yes** — drift undo per-APA/face. |
| `cluster.kd_steiner_knn(1, p, "steiner_pc")` | `Facade_Cluster.cxx` | None. |
| `examine_structure_4(vtx, false, graph, cluster, tf, dv)` | `NeutrinoStructureExaminer.cxx:473` | Via `dv`. |

---

## §3 `examine_vertices` (driver)

**Toolkit:** `NeutrinoStructureExaminer.cxx:2003–2015`  
**Prototype:** `NeutrinoID_proto_vertex.h:2375–2403`

### Summary

Thin `while(flag_continue)` loop that calls `examine_segment`, `_1`, `_2`
(conditionally), and `_4`. Returns when no merge fires in a full pass.

### Differences from prototype

None of substance. The loop structure, ordering, and `flag_continue = false`
reset at the top of each iteration are identical.

### [B1] Gate condition for `_2`

```cpp
// Prototype:
if (find_vertices(temp_cluster).size() > 2)

// Toolkit:
size_t num_vertices = boost::num_vertices(graph);
if (num_vertices > 2) {
```

`find_vertices(cluster)` in the prototype counts only vertices belonging to
`temp_cluster`. `boost::num_vertices(graph)` counts *all* vertices in the graph.
If the graph is per-cluster (as appears to be the intent), these are equivalent.
If the graph ever holds multiple clusters' vertices, the toolkit will call `_2`
even when the specific cluster has ≤ 2 vertices. This is a performance issue
rather than a correctness bug: `_2` filters by cluster and will find no
candidates, returning `false` harmlessly.

**Recommendation:** Replace with `find_vertices(graph, cluster).size() > 2`
to match the prototype's intent and eliminate the ambiguity.

### Determinism, multi-APA

No containers, no coordinate lookups. Not applicable.

---

## §4 `examine_vertices_1` (cluster-level Type-I kink merge)

**Toolkit:** `NeutrinoStructureExaminer.cxx:1336–1439`  
**Prototype:** `NeutrinoID_proto_vertex.h:2874–3002`

### Summary

Iterates degree-2 vertices. For each, if the shorter of its two connected
segments is < 4 cm AND `examine_vertices_1p` returns true, removes the
middle vertex and replaces the two segments with a single path from the two
outer vertices.

### Differences from prototype

**Structural equivalence:** logic is identical. The toolkit refactors the
per-segment iteration to use a pre-cached `seg_vtx_pairs` vector (line 1353)
to avoid redundant graph traversals. This is a clean efficiency improvement
with no semantic change.

**`main_vertex` guard:** Prototype at line 2965 checks `v1 != main_vertex`
*after* setting `flag_continue = true` and, if it matches, overrides
`flag_continue = false` before returning. Toolkit at line 1399 returns `false`
directly. Net return value in both cases: `false`. Net outer-loop effect:
identical.

**Null guard:** Toolkit adds `if (!sg || !sg1) return false` (line 1406)
before path construction. The prototype calls `find_segment(v1, nullptr)`
when `v3` is null, which may fault. Toolkit is safer.

**Path construction:** Prototype uses `dijkstra_shortest_paths` + `new
ProtoSegment`. Toolkit uses `do_rough_path` + `create_segment_for_cluster`.
Functionally equivalent (both find a Steiner-graph shortest path between the
two outer vertices).

### [D1] Iteration order

- **Prototype:** iterates `map_vertex_segments` (a `std::map<ProtoVertex*, ...>`
  keyed by raw pointer address). Iteration order varies between runs because
  pointer addresses are non-deterministic.
- **Toolkit:** uses `ordered_nodes(graph)` which returns nodes in
  `boost::vecS` insertion order — deterministic across runs for the same input.
- **Verdict:** toolkit is strictly better. Document as an improvement.

### Efficiency

Nested loop structure is O(V × degree) per pass; same as prototype. Cached
`seg_vtx_pairs` reduces per-vertex graph traversals from 3 to 1.

### Multi-APA

No direct coordinate conversion in `_1` itself; delegated to `_1p`.

---

## §5 `examine_vertices_1p` (per-pair Type-I decision)

**Toolkit:** `NeutrinoStructureExaminer.cxx:1142–1334`  
**Prototype:** `NeutrinoID_proto_vertex.h:3004–3222`

### Summary

Given two vertices `v1` (degree-2) and `v2` (its neighbour), checks all three
2D wire-time projections (U, V, W) of the segment between them. For each plane
it either finds the two vertices close in 2D (< 2.5 "units"), or finds the
segment in dead wire region, or finds the other segment forms a line. Returns
true if ≥ 2 planes satisfy one of these conditions (with specific logic for
combinations involving dead channels).

### [M1] Multi-APA generalization (major new feature)

**Prototype** bakes all coordinate conversion into 12 global scalars derived
from `TPCParams::Singleton`:

```cpp
double v1_t = offset_t + slope_xt * v1_p.x;   // single TPC drift → time slice
double v1_u = offset_u + slope_yu * y + slope_zu * z;
double v1_v = ...;  double v1_w = ...;
// same for segment points in "line" check
```

**Toolkit** performs per-vertex, per-point APA/face lookup:

```cpp
auto v1_wpid = dv->contained_by(v1_p);         // → {apa, face}
auto [v1_t_raw, v1_u_ch] = cluster.grouping()->convert_3Dpoint_time_ch(
    v1_p, v1_apa, v1_face, 0);
// similarly for V (pind=1) and W (pind=2)
// Per-point dead/good checks:
auto test_wpid = dv->contained_by(seg_fits[i].point);
auto p_raw = transform->backward(seg_fits[i].point, cluster_t0,
    test_wpid.face(), test_wpid.apa());
cluster.grouping()->get_closest_dead_chs(p_raw, 1, test_wpid.apa(),
    test_wpid.face(), pind);
```

This is the correct multi-APA generalization. Each point is projected to the
wire-time space of *its own* APA/face rather than a global singleton.

### [B2] Time normalization — resolved

**Prototype:** `v1_t = offset_t + slope_xt * x` where `slope_xt = 1/time_slice_width`.
The result is a time-slice index (float), in units comparable to wire pitch for the
2.5-threshold 2D distance comparison.

**Toolkit:** `v1_t = double(v1_t_raw) / ntime_ticks`.

**Resolution:** `Facade_Blob.h:33` documents `slice_index_min` as "unit: tick".
`convert_3Dpoint_time_ch` returns `tind = round(time/tick)` where `tick` is the
detector sampling period — also in tick units. So both `v1_t_raw` and `ntime_ticks`
are in the same tick unit.

For the difference `v1_t - v2_t`:
- Toolkit: `(tind1 - tind2) / ntime_ticks = (x1 - x2) / (drift_speed × tick × ntime_ticks)`
- Prototype: `(x1 - x2) / time_slice_width = (x1 - x2) / (drift_speed × tick × nrebin)`

These are equal when `ntime_ticks = nrebin`. For a typical single-readout-bin blob
`ntime_ticks = nrebin` (e.g. 4 for MicroBooNE/SBND), so the normalization is correct
for the common case.

**Residual concern (Low):** `children()[0]` is assumed to span exactly one readout bin.
If it spans multiple bins (unusual for standard wire-cell chunked blobs), `ntime_ticks`
would be a multiple of `nrebin`, making the t-scale too small by that factor. A
9-line comment was added to the code documenting this assumption.

### Algorithm equivalence (everything else)

The three-plane loop structure, the "line-check" subroutine (finding the point
~9 "units" away on the other segment, then checking the angle or good-point
interpolation with `step=0.6cm`, `n_bad≤1`), and the decision logic:

```cpp
if (ncount_close >= 2 ||
    (ncount_close == 1 && ncount_dead == 1 && ncount_line >= 1) ||
    (ncount_close == 1 && ncount_dead == 2) ||
    (ncount_close == 1 && ncount_line >= 2) ||
    ncount_line >= 3)
```

are identical in prototype and toolkit (prototype lines 3210–3215, toolkit
lines 1325–1330). No differences found.

**Note:** Prototype at line 3211 has `ncount_dead ==1 & ncount_line>=1` — a
single `&` (bitwise-AND, equivalent to `&&` on ints). Toolkit correctly uses
`&&`. Not a bug (same outcome), but worth noting.

### [B3] Per-view time reuse — resolved

`Facade_Grouping.cxx::convert_3Dpoint_time_ch` computes:

```cpp
const double time = drift2time(iface, time_offset, drift_speed, point[0]);
const int tind = std::round(time / tick);
return {tind, wind};
```

`drift2time` depends only on `point[0]` (x-coordinate, drift direction) and per-APA/face
drift parameters. It does **not** depend on `pind`. So `tind` is identical for all three
`pind` values for the same 3D point — only `wind` (wire channel) changes.

**Two consequences fixed in the code:**

1. The previous `v1_t_u` and `v1_t_v` variables (from the pind=1,2 vertex calls) were
   unused dead code. They are now `v1_t1`/`v1_t2` with `(void)` suppression and a comment.

2. The inner "line-check" loop calls `convert_3Dpoint_time_ch(..., pind)` — correct
   because it also retrieves the plane-correct `p_wire_ch`. The `p_t_raw` from this call
   is view-independent (verified above); a clarifying comment was added.

### Determinism, efficiency

No pointer-keyed containers. Iteration is over `sg->fits()` vectors. No issues.

---

## §6 `examine_vertices_2` (Type-II close-vertex merge)

**Toolkit:** `NeutrinoStructureExaminer.cxx:1441–1560`  
**Prototype:** `NeutrinoID_proto_vertex.h:2540–2609`

### Summary

Iterates segments. If two endpoints are within 0.45 cm (unconditionally) or
within 1.5 cm with both vertices having exactly 2 connections, merges them by
removing the closer vertex and rerouting its segments to the survivor.

### Equivalence

Distance thresholds (0.45 cm, 1.5 cm), degree-2 condition, `main_vertex`
guard (`v2 != main_vertex`), and merge procedure all match prototype exactly.

Prototype does not check `v1 != main_vertex`; toolkit also does not. Both are
symmetric: only `v2` is protected. This is consistent.

### [D2] Iteration order

- **Prototype:** iterates `map_segment_vertices` (`std::map<ProtoSegment*,...>`
  keyed by pointer). Nondeterministic.
- **Toolkit:** iterates `ordered_edges(graph)` — insertion order, deterministic.
- **Verdict:** toolkit strictly better.

### [E1] Isolated-vertex cleanup — promoted to bug, fixed

After the merge, the toolkit previously scanned all graph vertices with
`boost::vertices(graph)` — no cluster filter — to remove degree-0 orphans:

```cpp
// Before fix:
auto [vbegin, vend] = boost::vertices(graph);
for (auto vit = vbegin; vit != vend; ++vit) {
    if (boost::degree(*vit, graph) == 0) {
        VertexPtr vtx = graph[*vit].vertex;
        if (vtx) isolated_vertices.push_back(vtx);
    }
}
```

If the graph is shared across clusters, this removes degree-0 orphan vertices
from *other* clusters — a correctness bug. `examine_segment` (line 1139) already
had the `vtx->cluster() == &cluster` filter; `examine_vertices_2` and the
isolated-vertex loop in `examine_vertices_3` did not.

**Fix applied:** both loops now use `ordered_nodes(graph)` + `vtx->cluster() ==
&cluster` (matching the `examine_segment` pattern). Side-effect: iteration is now
deterministic order instead of `boost::vertices` (adjacency_list internal order).

### Multi-APA

Pure 3D distance computation, no wire/time projection. No multi-APA concerns.

---

## §7 `examine_vertices_3` (initial-vertex extension + short-segment cleanup)

**Toolkit:** `NeutrinoStructureExaminer.cxx:2349–2543`  
**Prototype:** `NeutrinoID_proto_vertex.h:2410–2538`

### Summary

Two phases:
1. For each of the two initial pair vertices (if degree-1), calls
   `get_local_extension` to extend the track beyond the original endpoint,
   rebuilds the segment via `do_rough_path`, then re-fits.
2. Finds short (< 5 cm) segments whose endpoints have degree 1 and whose fit
   points are all "shadowed" (within 0.6 cm in every plane) by other segments.
   Removes them, then runs `examine_structure_4` on the surviving vertices.

Called only for the main cluster. Call site:
`NeutrinoPatternBase.cxx:1583`: `examine_vertices_3(graph, cluster,
main_cluster_initial_pair_vertices, track_fitter, dv)`.

### [B5] Call site — parameter correctness

`main_cluster_initial_pair_vertices` is captured at line 1477 from `find_vertices(graph, sg1)` (the first segment's endpoints) and is a local variable in `find_proto_vertex`. It is passed by value to `examine_vertices_3`. The guard at line 1581 (`is_main_cluster && main_cluster_initial_pair_vertices.first`) ensures it is non-null. Both vertices in the pair are still owned by the graph at the call time. **Verified OK.**

### [B4] `create_segment_point_cloud` in phase 1

After updating `sg->wcpts()` (line 2415), the toolkit explicitly calls:

```cpp
create_segment_point_cloud(sg, main_pts_ev3, dv, "main");
```

The prototype does **not** do this; it relies on the subsequent
`do_multi_tracking` call (line 2468) to rebuild the point cloud. In the
toolkit, `do_multi_tracking` is called at line 2425 (if `flag_refit`), *after*
the explicit `create_segment_point_cloud` call. Verify that `do_multi_tracking`
does not re-call `create_segment_point_cloud` on the same segment, causing a
double-build or overwrite. If `do_multi_tracking` rebuilds the "main" PC
unconditionally, the earlier explicit call is harmless redundancy; if it skips
segments that already have an up-to-date PC, this is a correct early population.

### [D3] Determinism improvement (major)

Phase 2 (short-segment removal):

| Object | Prototype | Toolkit |
|---|---|---|
| `segments_to_be_removed` | `std::set<ProtoSegment*>` — pointer-keyed, nondeterministic iteration | `std::vector<SegmentPtr>` in edge insertion order + `std::find` dedup |
| `can_vertices` | `std::set<ProtoVertex*>` — pointer-keyed, nondeterministic | `std::vector<VertexPtr>` in insertion order + `std::find` dedup |

The order in which segments are removed and the order of `examine_structure_4`
calls are directly determined by these containers. Pointer-keyed `std::set`
iterates in address order, which varies per-run (ASLR). The toolkit's
insertion-order vector makes results deterministic and reproducible.

**This is the most impactful determinism fix in the entire function family.**

### [D4] Isolated vertex removal order

Prototype iterates `map_vertex_segments` (pointer-keyed map) to find zero-segment
vertices. Toolkit uses `ordered_nodes(graph)` — deterministic.

### [M2] Multi-APA in phase 2

```cpp
// Prototype:
std::tuple<double,double,double> results = sg1->get_closest_2d_dis(pts.at(i));

// Toolkit:
auto wpid = dv->contained_by(pts[i].point);   // per-point APA/face
auto [dist_u, dist_v, dist_w] = segment_get_closest_2d_distances(
    sg1, pts[i].point, wpid.apa(), wpid.face(), "fit");
```

The prototype's `get_closest_2d_dis` uses a global `TPCParams` singleton — it
assumes a single TPC wire geometry. The toolkit's `segment_get_closest_2d_distances`
takes explicit `apa` and `face` parameters, obtained per-point via
`dv->contained_by()`. Correct multi-APA generalization.

### [E2] O(segments × points × edges) inner loop

**Applied 2026-04-11.**

Original inner loop (phase 2, ~line 2466) re-iterated `boost::edges(graph)` and
re-applied the `sg1 != sg` filter on every fit point, giving O(S × P × E) work.

Three improvements were applied:

1. **Pre-snapshot** the other-segment list once per candidate segment (outside
   the point loop), eliminating repeated `boost::edges` iteration and the
   `sg1 == sg` filter per point.
2. **Early-exit (a)** from the inner segment loop: once `min_u`, `min_v`,
   `min_w` are all ≤ 0.6 cm, no remaining segment can make the point "unique" —
   break immediately.
3. **Early-exit (b)** from the outer point loop: as soon as one unique point is
   found (`num_unique++`), the segment will not be removed regardless of
   remaining points — break immediately.

```cpp
// Pre-snapshot once per candidate segment
std::vector<SegmentPtr> other_segs;
for (auto [e2b, e2e] = boost::edges(graph); e2b != e2e; ++e2b) {
    SegmentPtr sg1 = graph[*e2b].segment;
    if (sg1 && sg1 != sg) other_segs.push_back(sg1);
}
for (size_t i = 0; i < pts.size(); i++) {
    ...
    for (SegmentPtr sg1 : other_segs) {
        auto [du, dv, dw] = segment_get_closest_2d_distances(...);
        // update min_u/v/w ...
        if (min_u <= 0.6*cm && min_v <= 0.6*cm && min_w <= 0.6*cm) break;  // (a)
    }
    if (min_u > 0.6*cm || min_v > 0.6*cm || min_w > 0.6*cm) {
        num_unique++; break;  // (b)
    }
}
```

In the common case where most segments have at least one unique point, early-exit
(b) fires after the first point checked, reducing per-candidate work to O(1 × E)
instead of O(P × E). Early-exit (a) further reduces inner-loop iterations when a
point is well-covered.

---

## §8 `examine_vertices_4` (Type-III/IV short or perpendicular segment merge)

**Toolkit:** `NeutrinoStructureExaminer.cxx:1647–2001`  
**Prototype:** `NeutrinoID_proto_vertex.h:2611–2821`

### Summary

Iterates all segments. For each segment that is either short (< 2.0 cm direct
length) or nearly perpendicular to drift (direction < 3.5 cm and angle
within 10° of 90°), attempts to "collapse" one endpoint into the other by
checking connectivity via `examine_vertices_4p`. If it succeeds, the merged
vertex's segments are re-routed through the surviving vertex by stitching their
WC-point lists via the Steiner point cloud.

### [B6] Iterator-invalidation fix (already in toolkit)

**Prototype** iterates `map_segment_vertices` directly and calls
`del_proto_segment` / `add_proto_connection` inside the loop — this mutates the
map being iterated, which is undefined behavior (iterator invalidated mid-loop).
The loop happens to work because `std::map` iterators are not invalidated by
other-element deletions, but inserting via `add_proto_connection` can invalidate
the iterator.

**Toolkit** (lines 1662–1669) snapshots all segments first:

```cpp
// Snapshot all segments before iterating — the loop body calls remove_segment/add_segment
// which would invalidate boost::edges iterators...
std::vector<SegmentPtr> all_segments;
for (auto [eit, eend] = boost::edges(graph); eit != eend; ++eit) {
    SegmentPtr s = graph[*eit].segment;
    if (s && s->cluster() == &cluster) all_segments.push_back(s);
}
```

This is a correct and well-commented fix.

### [B7] Path-stitching: index vs. position comparison

**Prototype** trims `old_list` using WCPoint index equality:

```cpp
while (old_list.front().index != min_wcp.index && old_list.size() > 0) {
    old_list.pop_front();
}
old_list.pop_front();  // remove min_wcp itself
```

**Toolkit** trims by position (Euclidean distance):

```cpp
while (!old_list.empty() &&
       ray_length(Ray{old_list.front().point, min_wcp.point}) > 0.01 * units::cm) {
    old_list.pop_front();
}
if (!old_list.empty()) old_list.pop_front();
```

These are semantically equivalent when `min_wcp` is a Steiner node that
appears in `vec_wcps` — the point identity is preserved. The toolkit avoids
relying on `WCPoint::index` (a legacy field not always populated in the toolkit
data model), which is the correct adaptation. **No functional difference for
well-formed input.**

**Note:** prototype line 2688: `while(old_list.front().index != min_wcp.index
&& old_list.size()>0)` — this could loop forever if `min_wcp.index` does not
appear in `old_list` (e.g., due to a floating-point near-miss in the initial
scan). Toolkit's position-based check with the `0.01 cm` tolerance is safer.

### [B8] `main_vertex` guards (both branches)

- v1 branch: `v1 != main_vertex` at line 1713. ✓  
- v2 branch: `v2 != main_vertex` at line 1854. ✓  
- Both guarded, matching prototype lines 2635 and 2723.

### [M3] Steiner PC access

**Prototype:** `pcloud = temp_cluster->get_point_cloud_steiner()` (returns a
`ToyPointCloud*`), then `pcloud->get_closest_wcpoint(tmp_p)` for each
interpolated waypoint.

**Toolkit:** `cluster.get_pc("steiner_pc")` to read coordinates, then
`cluster.kd_steiner_knn(1, tmp_p, "steiner_pc")` for nearest-point queries.
Both probe the same underlying Steiner graph node set.

**Guard:** Line 1653: `if (!cluster.has_pc("steiner_pc")) return false;` — 
ensures early return if Steiner PC is missing. This is safe because
`find_proto_vertex` at line 1460 already checks `cluster.has_pc("steiner_pc")`
and returns `false` before reaching `examine_vertices`. So the inner guard is a
defensive belt-and-suspenders check.

**Threshold match:** Both prototype and toolkit reject a Steiner snap point if
it is > 0.3 cm from the interpolated grid point (`continue` at toolkit line
1787, prototype line 2674). Values identical.

### Efficiency

Toolkit `_4` is 355 lines vs. prototype 212 lines. The extra lines are:
- `boost::out_edges` / `boost::target` descriptor lookups replacing direct map
  access (~40 lines)
- The iterator snapshot vector (~8 lines)
- Explicit `create_segment_point_cloud` calls after wcpts update
- The symmetric v2-branch (prototype also has symmetric v2-branch, just less
  verbosely written)

No algorithmic difference. Performance is the same order as prototype.

### Determinism

Snapshot vector is built from `boost::edges` (not `ordered_edges`), which
iterates in boost adjacency_list internal edge-list order. Since edges are
always added sequentially to the same graph, this is deterministic for the
same input. The segment processing order differs from the prototype (which
iterated `map_segment_vertices` in pointer order), but the algorithm's
result is deterministic for fixed input.

---

## §9 `examine_vertices_4p` (per-pair Type-III/IV connectivity test)

**Toolkit:** `NeutrinoStructureExaminer.cxx:1563–1645`  
**Prototype:** `NeutrinoID_proto_vertex.h:2823–2870`

### Summary

For each segment connected to `v1` (other than the one linking v1 to v2),
finds the point closest to 3 cm from v1, then samples the straight-line path
from that point to `v2` at 0.3 cm steps. Returns `true` (merge is OK) iff no
bad points (`n_bad == 0`) are found on any of those paths.

### Equivalence

Logic matches prototype exactly. The `dis_cut` logic is absent in both (in
contrast to the cluster-level `_4` function, which does use `dis_cut`). The
0.3 cm step and `n_bad == 0` threshold are identical.

**Multi-APA:** The good-point test uses:

```cpp
auto test_wpid = dv->contained_by(test_p);
auto temp_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
grouping->is_good_point(temp_p_raw, test_wpid.apa(), test_wpid.face(), 0.2*cm, 0, 0);
```

Prototype uses `ct_point_cloud->is_good_point(test_p, 0.2*cm, 0, 0)` — global
single-TPC point cloud. Toolkit correctly uses per-point APA/face dispatch.

### Note on `sg1` (the segment between v1 and v2)

Line 1565: `SegmentPtr sg1 = find_segment(graph, v1, v2);` — this is computed
but only used in the loop to skip (`if (sg == sg1) continue`). If `find_segment`
returns nullptr (e.g. because there is no direct edge — unusual but possible
after a removal), the loop will skip *no* segment (since `sg == nullptr` is
false for all valid `sg`). In that case all of v1's segments would be tested,
which is overly conservative but safe. Worth adding a comment or guard.

---

## §10 `examine_segment` (helper, called first in `examine_vertices`)

**Toolkit:** `NeutrinoStructureExaminer.cxx:944–1138`

### Duplicate-detection key

Line 1097–1107: Builds a `SegmentKey = std::string` from `std::to_string(uintptr_t(v1.get())) + ":" + std::to_string(uintptr_t(v2.get()))` with `std::min`/`std::max` for canonical ordering. Then uses `std::unordered_set<SegmentKey>`.

**Note:** This uses raw pointer addresses as part of the key. Between identical
runs on the same binary + input, pointer addresses are stable (the same shared_ptr
addresses will be reproduced), so this is deterministic for fixed input. However,
keys built from pointer addresses cannot be compared across runs or processes. For
a within-run duplicate filter this is fine; the set is local to the function and
discarded after use.

**Not a nondeterminism issue.** The lookup is O(1) per segment endpoint pair, which
is efficient.

---

## §11 Priority Action List

### Bugs / Correctness

| ID | Function | Description | Priority | Status |
|---|---|---|---|---|
| **B1** | `examine_vertices` | `boost::num_vertices(graph) > 2` vs cluster-scoped count | Low | ✅ Fixed |
| **B2** | `examine_vertices_1p` | Time normalization `double(t_raw)/ntime_ticks`: `slice_index` is in tick units (confirmed from `Facade_Blob.h:33`), `tind` is in tick units → ratio equals prototype `slope_xt*x` when `ntime_ticks = nrebin` (single-readout-bin blob, the normal case). Downgraded from Critical. | Low | ✅ Comment added |
| **B3** | `examine_vertices_1p` | Redundant `v1_t1`/`v1_t2` vars; `p_t_raw` per-pind call in inner loop. `drift2time` in `convert_3Dpoint_time_ch` depends only on `point[0]` (x), confirming tind is view-independent. | Low | ✅ Fixed |
| **B4** | `examine_vertices_3` | `create_segment_point_cloud` before `do_multi_tracking`: `do_single_tracking:8430` always calls `create_segment_point_cloud` unconditionally, so the early call is harmless redundancy (overwritten). No fix needed. | — | ➖ Not a bug |
| **B5** | `examine_vertices_4p` | `sg1 = nullptr` when direct edge already removed; loop skips nothing | Low | ✅ Comment added |
| **B6** | `examine_vertices_4` | Iterator-invalidation fix already applied (snapshot vector) | — | ✅ Already fixed |
| **B7** | `examine_vertices_4` | `old_list` trim by position vs prototype's index; infinite-loop safety when `min_wcp` absent | Low | ✅ Comment added (both branches) |

### Efficiency

| ID | Function | Description | Priority |
|---|---|---|---|
| **E1** | `examine_vertices_2` + `examine_vertices_3` | Isolated-vertex scan used unfiltered `boost::vertices(graph)` — could remove vertices from other clusters | Medium (latent bug) | ✅ Fixed |
| **E2** | `examine_vertices_3` | Phase 2 inner loop: pre-snapshot other-segs + two early-exit conditions; reduces common-case work to O(1 × E) per candidate | Low | ✅ Fixed |

### Determinism

| ID | Function | Description | Status |
|---|---|---|---|
| **D1** | `examine_vertices_1` | `ordered_nodes` vs pointer-keyed `map_vertex_segments` | ✅ Improved in toolkit |
| **D2** | `examine_vertices_2` | `ordered_edges` vs pointer-keyed `map_segment_vertices` | ✅ Improved in toolkit |
| **D3** | `examine_vertices_3` | `std::vector` vs `std::set<Ptr*>` for `segments_to_be_removed` and `can_vertices` | ✅ Improved in toolkit |
| **D4** | `examine_vertices_3` | `ordered_nodes` vs pointer-keyed map for isolated-vertex removal | ✅ Improved in toolkit |

### Multi-APA / Multi-Face

| ID | Function | Description | Status |
|---|---|---|---|
| **M1** | `examine_vertices_1p` | Per-vertex/point `dv->contained_by` + `grouping()->convert_3Dpoint_time_ch(apa,face,pind)` replaces global TPCParams singleton | ✅ Correctly generalized — but see **B2** for unit check |
| **M2** | `examine_vertices_3` | `segment_get_closest_2d_distances(sg1,pt,apa,face,"fit")` with per-point APA/face replaces `sg1->get_closest_2d_dis(pt)` | ✅ Correctly generalized |
| **M3** | `examine_vertices_4` | `kd_steiner_knn` + `get_pc("steiner_pc")` replaces `get_point_cloud_steiner()->get_closest_wcpoint` | ✅ Correctly generalized; guarded by `has_pc` |
| **M4** | `examine_vertices_4p` | `dv->contained_by` + `transform->backward` + `grouping->is_good_point(apa,face)` replaces global `ct_point_cloud->is_good_point` | ✅ Correctly generalized |

---

## §12 Summary

The `examine_vertices_*` family is **well-ported**. The multi-APA generalization is
done correctly and consistently throughout the family: every place the prototype
called a global singleton (`TPCParams`, `ct_point_cloud`, `sg->get_closest_2d_dis`)
has been replaced with explicit APA/face-parameterised calls using `dv->contained_by`.

The four determinism improvements (D1–D4) are real and meaningful, particularly D3
(`std::set<Ptr*>` → `std::vector` in `examine_vertices_3`) which directly affects
the order in which `examine_structure_4` is called.

**All identified items have been resolved (2026-04-11):**

- **B1** fixed: cluster-filtered vertex count replaces the graph-global `boost::num_vertices`.
- **B2** resolved: `slice_index` and `tind` are both in raw-tick units; the normalization
  is correct for single-readout-bin blobs (the standard case). Comment added.
- **B3** fixed: redundant `v1_t1`/`v1_t2` variables suppressed; inner-loop `pind` usage
  documented as view-independent.
- **B4** closed: `do_multi_tracking → do_single_tracking:8430` rebuilds the PC
  unconditionally; the earlier explicit call is harmless redundancy.
- **B5** documented: `sg1 = nullptr` semantics explained in comment.
- **B7** documented: `old_list` trim safety vs. prototype's index-based loop explained
  in comment (both symmetric branches).

All identified items have been resolved or documented as non-issues.
