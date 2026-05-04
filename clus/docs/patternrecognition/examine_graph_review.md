# Review: `Examine_graph` / `Connect_graph_overclustering_protection`

**Reviewed by:** Claude Code  
**Date:** 2026-04-11  
**Branch:** `apply-pointcloud`  
**Prototype ref:** `/nfs/data/1/xqian/prototype-dev/wire-cell/data/src/PR3DCluster.cxx` lines 1855–2353

---

## Implementation status

| Item | Description | Status | File |
|---|---|---|---|
| B.1 | `wpid_U/V/W_dir` rekeyed to `{apa,face}` pair | ✅ Applied | `connect_graph_relaxed.cxx` |
| B.2 | Out-of-volume path steps counted as bad (main loop + dir1 + dir2) | ✅ Applied | `connect_graph_relaxed.cxx` |
| B.3 | Out-of-volume steps counted as bad in `check_connectivity` | ✅ Applied | `connect_graph_relaxed.cxx` |
| B.4 | Cross-reference comment added to `Separate_overclustering` | ✅ Applied | `clustering_protect_overclustering.cxx` |
| B.5 | `ordered_components.reserve(num)` (not `component.size()`) | ✅ Applied | `clustering_protect_overclustering.cxx` |
| B.6 | Dead if/else for `dis` collapsed (3× in `connect_graph_relaxed`, 3× in `clustering_protect_overclustering`) | ✅ Applied | both files |
| B.2-equiv | Out-of-volume steps counted as bad in `check_path` lambda | ✅ Applied | `clustering_protect_overclustering.cxx` |
| C.2 | Linear wire-interval scan → `lower_bound/upper_bound` | ✅ Applied | `clustering_protect_overclustering.cxx` |
| C.3 | Redundant `num²` initialization loop eliminated in both files | ✅ Applied | both files |
| C.1 | `std::map<Blob*,…,BlobLess>` → `unordered_map` in `connect_graph_closely` | ✅ Applied | `connect_graph_closely.cxx` |
| C.4 | Early-exit Hough when closest distance > 80 cm | ✅ Applied | `connect_graph_relaxed.cxx` |

---

## 1. Scope

This review audits the WCP prototype methods

- `PR3DCluster::Examine_graph(ct_point_cloud)` (line 2311)
- `PR3DCluster::Connect_graph_overclustering_protection(ct_point_cloud)` (line 1855)
- `PR3DCluster::Establish_close_connected_graph()` (called inside both of the above)

against their WCT toolkit counterparts.  Goals are: (1) verify functional equivalence, (2) find bugs, (3) propose efficiency improvements, (4) reduce nondeterminism, (5) audit called helpers, (6) verify multi-APA / multi-face correctness.

---

## 2. Function mapping

| WCP | WCT | File:line |
|---|---|---|
| `Establish_close_connected_graph()` | `Graphs::connect_graph_closely()` | `clus/src/connect_graph_closely.cxx:12` |
| `Connect_graph_overclustering_protection(ct_pc)` | `Graphs::connect_graph_relaxed()` (primary) | `clus/src/connect_graph_relaxed.cxx:14` |
| — | `Separate_overclustering()` (duplicate, see §B.5) | `clus/src/clustering_protect_overclustering.cxx:47` |
| `Examine_graph(ct_pc)` return value | `Cluster::connected_blobs(dv, pcts, "relaxed")` | `clus/src/Facade_Cluster.cxx:2968` |
| Full pipeline | `Graphs::make_graph_relaxed(cluster, dv, pcts)` | `clus/src/make_graphs.cxx:78` |
| `VHoughTrans(p, r, sub_cloud)` | `Cluster::vhough_transform(p, r, theta_phi, sub_cloud, global_idx)` | called from `connect_graph_relaxed.cxx:124` |
| `ct_point_cloud.test_good_point(p)` | `grouping->test_good_point(p_raw, apa, face)` | `connect_graph_relaxed.cxx:182` |
| `ct_point_cloud.is_good_point_wc(p)` | `grouping->is_good_point(p_raw, apa, face)` | `connect_graph_relaxed.cxx:334` |
| `pt_cloud->get_closest_points(other)` | `Simple3DPointCloud::get_closest_points(*other)` | `Facade_Util.cxx` |
| `VHoughTrans` + `MCUGraph` MST | `process_mst_deterministically` | `Facade_Util.cxx:58` |
| *(none — single APA)* | `Facade::get_wireplaneid(...)` (2 overloads) | `Facade_Util.cxx:724,732` |

PID variants of `connect_graph_closely` and `connect_graph_relaxed` exist in the same source files and are omitted from prototype comparison since the prototype lacked them.

### Call graph

```
Cluster::connected_blobs(dv, pcts, "relaxed")
  └─ Cluster::graph_algorithms("relaxed", dv, pcts)        [cached]
       └─ Graphs::make_graph_relaxed(cluster, dv, pcts)
            ├─ Graphs::make_graph_closely(cluster)
            │    └─ Graphs::connect_graph_closely(cluster, graph)
            └─ Graphs::connect_graph_relaxed(cluster, dv, pcts, graph)
                 ├─ Cluster::vhough_transform(...)
                 ├─ Simple3DPointCloud::get_closest_points(...)
                 ├─ Simple3DPointCloud::get_closest_point_along_vec(...)
                 ├─ Facade::get_wireplaneid(...)
                 ├─ grouping->test_good_point / is_good_point
                 └─ Facade::process_mst_deterministically(...)

ClusteringProtectOverclustering visitor  (separate parallel path)
  └─ Separate_overclustering(cluster, dv, pcts, scope)
       ├─ [rebuilds its own close graph internally]
       ├─ [runs its own simple check_path]
       └─ grouping->separate(cluster, b2groupid, ...)
```

---

## 3. §A — Functional equivalence

### A.1 `Examine_graph` → `connected_blobs`

WCP `Examine_graph` does:
1. Delete old graph.
2. Call `Create_point_cloud()` if needed.
3. `graph = new MCUGraph(N)`.
4. `Establish_close_connected_graph()` — fills intra-/inter-blob edges.
5. `Connect_graph_overclustering_protection(ct_pc)` — bridges components across dead channels.
6. `connected_components` → partition blob set by component id.

WCT `connected_blobs` (line 2968) does:
1. Call `graph_algorithms("relaxed", dv, pcts)` which calls `make_graph_relaxed` — steps 3–5.  The result is **cached by flavor** on the cluster, so repeat calls within one event are free.
2. `ga.connected_components()` on the cached graph.
3. Remap point-index → blob-index via `kd3d().major_index(i)`.

**Verdict: functionally correct and superior** — the caching avoids re-running the expensive overclustering-protection scan if `connected_blobs` is called more than once on the same cluster in one event.

### A.2 Main loop of `connect_graph_relaxed` vs `Connect_graph_overclustering_protection`

Lines 110–477 of `connect_graph_relaxed.cxx` are a near-1:1 port of WCP lines 1874–2296:

| Feature | WCP | WCT | Parity |
|---|---|---|---|
| `get_closest_points` every (j,k) pair | line 1912 | line 113 | ✓ |
| Hough probes only for large clouds (`num<100 && >400 \|\| >500`) | lines 1913–1915 | lines 116–118 | ✓ |
| 6-score path test (`scores[0..5]`) | lines 1970–1987 | lines 185–201 | ✓ |
| `flag_strong_check` angle classification | lines 2016–2049 | lines 242–301 | ✓ |
| `3 cm` distance override for MST | line 2251 | line 524 | ✓ |
| `5 cm` → `*1.1` weight bump | lines 2276, 2287 | lines 549, 563 | ✓ |
| dir-MST stores `index_index_dis`, not dir1/dir2 | line 2235,2237 | `process_mst_deterministically` call at line 518 | ✓ (intentional, see A.3) |

### A.3 Known prototype quirks carried over intentionally

**Quirk 1 — WCP dangling-else bug fixed in WCT:**

WCP line 2215–2216:
```cpp
if (std::get<0>(index_index_dis_dir1[j][k])>=0 || std::get<0>(index_index_dis_dir2[j][k])>=0)
  auto edge = add_edge(index1, index2, ..., temp_graph);
```
Because `auto` declarations cannot be the body of an unbraced `if`, this was compiled as-written and the intent was "only add when condition holds". In WCT `connect_graph_relaxed.cxx:506–512`, braces are added and the guard is correct. **This is a genuine prototype bug fix, not a porting error.**

**Quirk 2 — dir-MST uses `index_index_dis` not `dir1`/`dir2` as source:**

In both WCP lines 2234–2238 and the WCT `process_mst_deterministically` call at line 518, the MST of the *direction connectivity graph* records the standard closest-pair tuple as its MST entry, not the directional pair.  This is fine because the dir-MST result is only queried for its index being `>=0` (i.e., as a boolean "these two components are directionally connected via the MST").  The actual directional edge endpoints come from `index_index_dis_dir1/dir2` at lines 545/559.  No action needed, but add a clarifying comment.

### A.4 `Establish_close_connected_graph` → `connect_graph_closely`

The intra-blob graph (lines 68–161 in `connect_graph_closely.cxx`) and the inter-blob cross-time-slice graph (lines 163–466) are correct ports.  The WCT version uses `lower_bound/upper_bound` for the wire-interval lookup (lines 118–126) instead of the linear scan in the prototype — this is an efficiency improvement (see §C.3).  The multi-APA partitioning via `af_time_slices` (lines 163–232) is a correct new addition (see §E.1).

---

## 4. §B — Bugs

### B.1 `connect_graph_relaxed.cxx:32–60` — `wpid_U_dir` keyed by full `WirePlaneId` including layer (**multi-APA bug**) ✅ Fixed

**Location:** `connect_graph_relaxed.cxx` lines 32–60, lookups at lines 212, 218, 224, 352, 357, 367, 438, 443, 452.

**Problem:** The setup loop
```cpp
for (const auto& wpid : wpids) {          // wpid has a specific layer (U/V/W)
    wpid_U_dir[wpid] = geo_point_t(0, cos(angle_u), sin(angle_u));
    wpid_V_dir[wpid] = geo_point_t(0, cos(angle_v), sin(angle_v));
    wpid_W_dir[wpid] = geo_point_t(0, cos(angle_w), sin(angle_w));
}
```
stores three identical U/V/W direction vectors at three *different* WirePlaneId keys (one per layer) for each apa/face.

Later `test_wpid = get_wireplaneid(p1, wpid_p1, p2, wpid_p2, dv)` returns the wpid whose *inner_bounds* intersection with the ray is longer.  That wpid's layer comes from `dv->inner_bounds(wpid1/wpid2)` / `dv->contained_by(point)`.  If the returned layer happens not to match any key in the map for that apa/face (e.g. the map was populated with a U/V/W wpid but `contained_by` returns a W-layer wpid for a different apa), `wpid_U_dir.at(test_wpid)` throws `std::out_of_range`.

**Fix:** Key the three maps on `std::pair<int,int>{apa, face}` instead of the full `WirePlaneId`.

```cpp
// Replace
std::map<WirePlaneId, geo_point_t> wpid_U_dir;
// With
std::map<std::pair<int,int>, geo_point_t> afpair_U_dir;

// In setup loop:
auto af = std::make_pair(wpid.apa(), wpid.face());
if (afpair_U_dir.find(af) == afpair_U_dir.end()) {
    afpair_U_dir[af] = geo_point_t(0, cos(angle_u), sin(angle_u));
    afpair_V_dir[af] = geo_point_t(0, cos(angle_v), sin(angle_v));
    afpair_W_dir[af] = geo_point_t(0, cos(angle_w), sin(angle_w));
}

// Lookups:
auto af = std::make_pair(test_wpid.apa(), test_wpid.face());
double angle1 = tempV1.angle(afpair_U_dir.at(af));
```

This also reduces 3× memory waste when `wpids` has all three layers for each apa/face.

### B.2 `connect_graph_relaxed.cxx:166–204` — out-of-volume path steps not counted (**multi-APA bug**) ✅ Fixed

**Location:** step loop at lines 166–204; same issue at lines 320–343 and 406–429.

**Problem:** When `test_wpid.apa() == -1` (point is between two APA volumes), the entire `if` body is skipped and nothing is added to `num_bad*` counters.  The denominator `num_steps` still includes these skipped steps.  In the prototype this situation never arose because there was only one volume.  In a multi-APA detector, a substantial gap between APAs can make a connection appear much "better" than it really is.

**Fix:** Count out-of-volume steps as bad:
```cpp
if (test_wpid.apa() != -1) {
    // ... existing scoring code ...
} else {
    // step is between APA volumes — treat as all-bad
    num_bad[0]++;  num_bad[1]++;  num_bad[2]++;  num_bad[3]++;
    num_bad1[0]++; num_bad1[1]++; num_bad1[2]++; num_bad1[3]++;
}
```

Apply the same fix to the dir1/dir2 loops (`num_bad++` and `num_bad1++` in the equivalent `else` branch).

### B.3 Same issue in `check_connectivity` (`connect_graph_relaxed.cxx:659–720`) ✅ Fixed

**Location:** line 669 `if (test_wpid.apa() == -1) continue;`

Same gap: out-of-volume steps are silently skipped, making connections between components that pass through dead inter-APA space appear valid.  Apply same fix: increment `num_bad[3]++` when `apa == -1`.

### B.4 `clustering_protect_overclustering.cxx` — duplicate and divergent implementation ✅ Comment added

**Location:** `clustering_protect_overclustering.cxx:47–563`, called by `ClusteringProtectOverclustering::visit()`.

**Problem:** `Separate_overclustering()` re-implements the entire close-graph + overclustering-protection pipeline from scratch.  It diverges from `connect_graph_relaxed` in several ways:

- Uses a simpler `check_path` (lines 352–381) that applies only the hard `num_bad > 7 || num_bad > 0.75*steps` threshold, discarding the full `num_bad[4]` / `num_bad1[4]` / angle classification logic.
- Does **not** apply the `!boost::edge(...)` guard on `add_edge` for the main cluster graph (`*graph`) — see lines 156, 264, 501, 519, 535.  Since `Weighted::Graph` uses `vecS` edge storage, this creates duplicate edges (harmless for `connected_components` but wasteful).
- Does not apply the `3 cm` MST override.
- Does not apply the `5 cm → *1.1` weight bump for directional edges.
- Bug fixes applied to `connect_graph_relaxed` (e.g. B.1, B.2) **will not be inherited**.

**Recommendation:** Keep `clustering_protect_overclustering.cxx` as a standalone visitor but ensure that any bug fix applied to `connect_graph_relaxed.cxx` (especially B.1–B.3 above) is mirrored into `Separate_overclustering`.  The two code paths serve different callers in the pipeline and it is intentional to keep them separate.  Add a comment at the top of `Separate_overclustering` cross-referencing `connect_graph_relaxed.cxx` as the authoritative implementation to keep in sync with.

### B.5 `clustering_protect_overclustering.cxx:279–297` — `ordered_components` over-allocated ✅ Fixed

Minor: `ordered_components.reserve(component.size())` at line 281 reserves `num_vertices` entries but only `num` (connected-component count, typically ≪ num_vertices) are filled.  Should be `reserve(num)`.

### B.6 Redundant `if/else` with identical bodies in `connect_graph_relaxed.cxx:532–538` and `clustering_protect_overclustering.cxx:492–498` ✅ Fixed

Both have:
```cpp
float dis;
if (std::get<2>(index_index_dis_mst[j][k]) > 5 * units::cm) {
    dis = std::get<2>(index_index_dis_mst[j][k]);
} else {
    dis = std::get<2>(index_index_dis_mst[j][k]);
}
```
Both branches are identical.  This is a copy-paste residue from the prototype.  Simplify to:
```cpp
float dis = std::get<2>(index_index_dis_mst[j][k]);
```
(The intent was probably `*1.1` in the `>5 cm` branch, matching the directional edge treatment at line 549/563; whether that was intentionally omitted for the closest-pair MST edge is worth confirming.)

---

## 5. §C — Efficiency improvements

### C.1 `connect_graph_closely.cxx:17–63` — expensive comparator on outer map ✅ Fixed

**Current code:**
```cpp
using mcell_wire_wcps_map_t = std::map<const Blob*, std::map<int, std::set<int>>, BlobLess>;
mcell_wire_wcps_map_t map_mcell_uindex_wcps, ...;
```

`BlobLess` resolves wpid + geometry on every comparison, making each `.at(mcell)` lookup `O(log N_blobs)` with a non-trivial constant.  The maps are only ever accessed by exact key (never iterated in order for output).

**Fix:** Use `std::unordered_map<const Blob*, …>` with `std::hash<const void*>` (or the simpler `std::hash<const Blob*>`).  Iteration order of the outer map is never the source of determinism — all loops are driven by `cluster.children()` and `time_blob_map()`.

```cpp
struct BlobPtrHash {
    std::size_t operator()(const Blob* b) const noexcept {
        return std::hash<const void*>{}(static_cast<const void*>(b));
    }
};
using mcell_wire_wcps_map_t = std::unordered_map<const Blob*, std::map<int, std::set<int>>, BlobPtrHash>;
```

Apply the same change in `connect_graph_closely_pid`.

**Reference pattern:** `clustering_protect_overclustering.cxx:71–92` already does the right thing — it maps blobs through an integer `blob_to_idx` table built once, then uses `std::vector<std::map<int,std::set<int>>>` indexed by integer.  This is the preferred pattern for new code.

### C.2 `clustering_protect_overclustering.cxx:121–131` — linear wire-interval scan ✅ Fixed

**Current code:**
```cpp
for (auto it2 = map_max_index_wcps->begin(); it2 != map_max_index_wcps->end(); it2++) {
    if (std::abs(it2->first - index_max_wire) <= max_wire_interval) { ... }
}
```
This is `O(W)` where W is the number of distinct wire indices in the blob.

**Fix:** Use the `lower_bound/upper_bound` idiom already used in `connect_graph_closely.cxx:118–126`:
```cpp
auto lo = map_max_index_wcps->lower_bound(index_max_wire - max_wire_interval);
auto hi = map_max_index_wcps->upper_bound(index_max_wire + max_wire_interval);
for (auto it2 = lo; it2 != hi; ++it2) max_wcps_set.push_back(&it2->second);
```

### C.3 `connect_graph_relaxed.cxx:87–102` — redundant initialization pass ✅ Fixed (both files)

The five `index_index_dis*` arrays are value-initialized to `{0,0,0.0}` by the vector constructor but immediately overwritten to `{-1,-1,1e9}` in the loop at lines 94–102.  Use direct construction:

```cpp
const auto sentinel = std::make_tuple(-1, -1, 1e9);
const std::vector<std::tuple<int,int,double>> row(num, sentinel);
std::vector<std::vector<std::tuple<int,int,double>>>
    index_index_dis(num, row),
    index_index_dis_mst(num, row),
    index_index_dis_dir1(num, row),
    index_index_dis_dir2(num, row),
    index_index_dis_dir_mst(num, row);
```

This eliminates a 5 × num² write pass.

### C.4 `connect_graph_relaxed.cxx:124–127` — dominant-cost Hough calls not short-circuited ✅ Fixed

For clusters with many disconnected components (noisy events), the double-loop runs O(num²) Hough calls.  Each `vhough_transform` at 30 cm radius is a full k-d tree range query + histogram.

**Cheap early-exit:** before entering the Hough-probe block (line 116), check if the closest-pair distance already fails the hard threshold that will follow (`> 80 cm` or `> 3 cm` with poor path quality) and skip the Hough entirely.  At a minimum, skip when `std::get<2>(index_index_dis[j][k]) > 80 * units::cm` since `get_closest_point_along_vec` won't return anything useful beyond that range anyway.

```cpp
if (/* large-cloud condition */ && std::get<2>(index_index_dis[j][k]) < 80 * units::cm) {
    // ... Hough probes ...
}
```

### C.5 `connect_graph_relaxed.cxx:105–107` — scope-transform objects allocated outside loops (already done)

The `needs_transform`, `ctpc_transform`, and `cluster_t0` are correctly hoisted outside the component-pair loop (lines 105–107).  ✓  No action.

---

## 6. §D — Determinism audit

### D.1 `connect_graph_closely` — map ordering

`std::map<const Blob*, …, BlobLess>` at line 17: `BlobLess` uses `blob_less()` in `Facade_Blob.cxx:248` which sorts by `wpid()` (integer-based), then by geometry.  **Deterministic.** ✓

The outer loops that produce the *output graph edges* are always driven by `cluster.children()` (children insertion order in the n-ary tree) and `time_blob_map()` (sorted `std::map<int, ...>`).  **Deterministic.** ✓

### D.2 `connect_graph_relaxed` — component labeling

`connected_components` is called on the closely-graph (line 67).  Boost BGL assigns component ids in BFS/DFS order starting from vertex 0.  Because the closely graph's edge insertion order is deterministic (D.1), the component labels are deterministic for the same input. ✓

### D.3 `process_mst_deterministically` — MST tie-breaking

`Facade_Util.cxx:58–117` sorts components by `(size desc, min_vertex asc)` and then runs `prim_minimum_spanning_tree` with root `comp_vertices[0]` (the smallest vertex index).  Prim's tie-breaking depends on the adjacency list structure; since `Weighted::Graph` uses `vecS` (adjacency stored as `std::vector`), out-edge iteration order is insertion order — which is deterministic per D.2.  **Deterministic.** ✓

### D.4 Main cluster split (`clustering_examine_bundles.cxx:108–123`)

`overlap_counts` is `std::map<int,int>` (sorted by int key = new group id).  The scan to find `max_overlap` at line 118 picks the **lowest** group id on ties.  **Deterministic.** ✓

### D.5 No raw-pointer-keyed `std::map` with default comparator

All `Blob*` maps use `BlobLess` or integer indexing.  **No pointer-order nondeterminism.** ✓  
*(After the fix in C.1, the `unordered_map` outer map still does not affect output order.)*

---

## 7. §E — Multi-APA / multi-face audit

### E.1 `connect_graph_closely` — correct APA/face partitioning ✓

Lines 164–232 build `af_time_slices[apa][face]` and all `connected_mcells` pairs are within the same `(apa, face)`.  Wire-index comparison only happens between blobs in the same APA, so wire-index space collisions between APAs cannot occur.

`nticks_per_slice.at(apa).at(face)` is used for time-adjacency thresholds (lines 205, 208, 211, 683, 686) — correct per-face generalization of the single-value prototype constant. ✓

### E.2 `connect_graph_relaxed` — cross-APA path checks

Components in the relaxed graph can span multiple APAs (e.g. a track crossing an APA boundary forms one large component after the closely graph is built).  Path checks between components therefore potentially cross APA boundaries.

`get_wireplaneid(test_p, wpid_p1, wpid_p2, dv)` at each path step identifies the correct APA/face for intermediate points. ✓

`grouping->test_good_point(test_p_raw, test_wpid.apa(), test_wpid.face())` consults the per-face dead-channel map. ✓

**But see bug B.1** (wpid_U_dir key mismatch) and **bug B.2** (silent skip when between volumes).

### E.3 Wire-angle computation for cross-APA pairs

When `wpid_p1 != wpid_p2`, `get_wireplaneid(p1, wpid_p1, p2, wpid_p2, dv)` (the 5-arg overload at `Facade_Util.cxx:732`) picks the APA/face whose `inner_bounds` contains the longer segment of the p1→p2 ray.  This is used to select the U/V/W wire angles for the angle classification.

This is a reasonable heuristic for the majority of the connection.  For a pair straddling the exact midpoint of two equally-sized APA volumes the choice is arbitrary but stable (tied to geometry, not pointers). ✓

After fixing B.1, the lookup `afpair_U_dir.at({test_wpid.apa(), test_wpid.face()})` will always succeed as long as the apa/face exists in the grouping's `wpids()`. ✓

### E.4 `grouping->separate()` with multi-APA blobs

`grouping->separate(cluster, b2groupid, ...)` at `clustering_protect_overclustering.cxx:557` splits a cluster by component id.  Each blob carries its own wpid, so blobs from different APAs in the same component remain correctly associated after the split.  No multi-APA issue. ✓

### E.5 `connected_blobs` remapping

`kd3d().major_index(i)` returns the blob-level index in `cluster.children()` order.  Since `cluster.children()` is insertion-order and blobs from all APAs appear there, the mapping is correct for multi-APA clusters. ✓

---

## 8. §F — Open questions / recommendations

1. **`Separate_overclustering()` is kept** (see B.4).  The `ClusteringProtectOverclustering` visitor and `ClusteringExamineBundles` are separate pipeline stages and their implementations remain distinct.  The required action is to keep them synchronized when bugs are fixed in `connect_graph_relaxed`.

2. **Out-of-volume steps are counted as bad** (B.2–B.3 fix confirmed).  Counting them as bad (all six scores zero) is conservative, matches the prototype spirit, and correctly prevents connections through large inter-APA dead regions.

3. **Is the `>5 cm → *1.1` weight bump intentionally absent for the closest-pair MST edge?** (B.6)  The prototype has it for the directional edges but not the closest-pair edge.  WCT preserves this asymmetry.  Confirm it is intentional.

4. **C.4 short-circuit distance threshold** — `80 cm` is suggested but should be confirmed against the maximum useful range of `get_closest_point_along_vec` as called in this context.

---

## 9. Verification checklist

- [ ] Build toolkit after B.1 fix: `./wcb` passes.
- [ ] Run `clus/test/test-porting/pdhd/clus.jsonnet` — diff cluster output (blob counts per CC, CC count per cluster) before/after B.1 fix.
- [ ] Run `clus/test/test-porting/fgval/stage2.jsonnet` (multi-APA) — verify no `std::out_of_range` from `wpid_U_dir.at()`.
- [ ] Run same event twice with same settings and diff output — confirm determinism (D.1–D.5).
- [ ] After B.4 refactor: rerun both pdhd and fgval and confirm numerical equivalence of cluster split counts.
- [ ] After C.1 (unordered_map) — rerun determinism diff to confirm no ordering change.
- [ ] Profile `connect_graph_relaxed` wall time before/after C.3–C.4 on the 1000-cluster pdhd sample.
