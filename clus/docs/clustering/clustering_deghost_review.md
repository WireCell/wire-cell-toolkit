# Review: `ClusteringDeghost` / `clustering_deghost()`

**Reviewed:** 2026-04-10  
**Reviewer:** Claude (assisted)  
**Toolkit file:** `clus/src/clustering_deghost.cxx` (795 lines)  
**Prototype file:** `prototype-dev/wire-cell/2dtoy/src/ToyClustering_deghost.h` (1104 lines, two overloads)

## Review axes

1. Functional parity with prototype  
2. Bug detection  
3. Efficiency / algorithmic improvements  
4. Determinism (pointer-keyed associative containers)  
5. Called-helper review  
6. Multi-APA / multi-face correctness

---

## 1. Files and functions covered

| Artifact | Path |
|---|---|
| Facade class | `clus/src/clustering_deghost.cxx:78-104` |
| Entry point | `clus/src/clustering_deghost.cxx:119-795` |
| Prototype overload A (no CTPC) | `ToyClustering_deghost.h:6-538` |
| Prototype overload B (with CTPC) | `ToyClustering_deghost.h:541-1104` |
| `DynamicPointCloud` helpers | `clus/src/DynamicPointCloud.cxx` |
| `merge_clusters` | `clus/src/ClusteringFuncs.cxx:48-120` |
| `extract_geometry_params` | `clus/src/ClusteringFuncs.cxx:12-46` |
| `NaryTreeFacade::destroy_child` | `util/inc/WireCellUtil/NaryTreeFacade.h:239` |
| `DetUtils::face_parameters` | `clus/src/DetUtils.cxx:15-44` |

Called prototype callers:
- `ToyClustering.cxx:374` — calls **Overload B** (`use_ctpc=true` path)  
- `ToyClustering.cxx:566` — calls **Overload A** (`use_ctpc=false` path)

---

## 2. Prototype → toolkit stage map

| Stage | Prototype A lines | Prototype B lines | Toolkit lines | Notes |
|---|---|---|---|---|
| Sort clusters by length (desc) | 11-21 | 546-556 | 180-182 | Toolkit: `std::sort` with lambda; prototype: pair-vec + `std::sort`. Both unstable on ties → minor non-determinism on equal-length clusters. |
| TPC/geometry params | 28-35 | 563-570 | 127-170 | Prototype: `TPCParams` singleton. Toolkit: iterates `dv_wpids()` (one `kAllLayers` wpid per (apa,face)); builds `wpid_params` map. Correct. |
| Create `global_point_cloud` + `global_skeleton_cloud` | 38-39 | 573-574 | 190-192 | Prototype: `DynamicToyPointCloud(angle_u,v,w)`. Toolkit: `DynamicPointCloud(wpid_params)`. Correct. |
| Declare `to_be_removed_clusters`, merge bookkeeping | 41-42 | 576-577 | 195-204 | Prototype: `std::set<std::pair<Cluster*,Cluster*>> to_be_merged_pairs`. Toolkit: `cluster_connectivity_graph_t g` + `merge_clusters()`. Better. |
| Seed cluster (i==0) | 45-52 (A), 581-590 (B) | — | 216-232 | See §3-B. |
| Per-point U/V/W dead/close scan | 71-180 | 609-742 | 249-426 | Ghost predicate, thresholds, early break: all match. |
| Ghost branch — max cluster selection | 207-224 | 771-788 | 463-480 | **B1 determinism bug** — see §4. |
| Ghost branch — three-way merge/remove decisions | 231-301 | 795-865 | 484-581 | Logic matches. |
| Non-ghost / secondary-heuristic branch | 306-413 | 870-977 | 582-738 | Logic matches. Max cluster selection at 592-609: **same B1 bug**. |
| Add surviving clusters to global clouds | 418-425 | 982-989 | 742-751 | Match. |
| Above-`length_cut` path (unconditional add) | 426-433 | 990-997 | 753-762 | Match. |
| Merge via graph / union-find | 440-476 | 1004-1041 | 769 (`merge_clusters`) | Toolkit replaces prototype's manual union-find with boost CC. Better. |
| Destroy removed clusters | 504-509 | 1069-1075 | 773-777 | Match. |

---

## 3. Functional parity

### 3-A. Ghost predicate (L435-455 vs prototype L186-200): **MATCH**

The 14-clause `&&`/`||` predicate at toolkit L435-455 matches the prototype
token-for-token including the unsuppressed operator-precedence warning (suppressed
by `#pragma GCC diagnostic ignored "-Wparentheses"` at L108-109). No logic
difference.

### 3-B. Seed cluster skeleton handling: **PARTIAL DIFFERENCE**

The toolkit always adds a short (≤30 cm) seed cluster to `global_skeleton_cloud`
using the full point cloud as fallback (L227-230). This matches **Overload B**
(prototype L586-589). When `use_ctpc=false` (Overload A path), the prototype
does NOT add short seed clusters to `global_skeleton_cloud`, but the toolkit
does. Effect: with `use_ctpc=false`, the skeleton cloud is slightly more dense
for short leading clusters, which can only make ghost detection slightly stricter
(more reference points to match against). Since SBND (the only active user)
runs with `use_ctpc=true`, this difference is currently inactive.

### 3-C. Long-cluster skeleton construction: **MATCH**

Prototype Overload B: `cluster->Construct_skeleton(ct_point_cloud)` then
`AddPoints(cluster, flag=1)` (skeleton points only).  
Toolkit: `get_path_wcps(*cluster, dv, pcts, use_ctpc)` followed by
`make_points_cluster_skeleton(..., path_wcps, flag_wrap=true)`.  
`flag_wrap=true` adds multi-face wrapped 2D projections, which is the
correct multi-TPC generalization of the single-TPC skeleton. Semantically
equivalent.

### 3-D. Merge bookkeeping: **IMPROVED**

The prototype builds `std::set<std::pair<Cluster*,Cluster*>>` and then
runs a manual union-find loop (lines 440-476). The toolkit uses a boost
`cluster_connectivity_graph_t` (adjacency list) + `boost::connected_components`
inside `merge_clusters()` (`ClusteringFuncs.cxx:48`). This is deterministic
(all int-keyed), correct, and avoids the O(M²) prototype union-find cost.

### 3-E. `assert(live == nullptr)` after `destroy_child`: **CORRECT**

`NaryTreeFacade::destroy_child(child_type*& kid)` sets `kid = nullptr`
(`NaryTreeFacade.h:239-242`). The assert is valid.

---

## 4. Bugs found

### BUG-1 — Non-deterministic max-cluster selection (high severity)

**Location:** `clustering_deghost.cxx:245, 463-480, 592-609`

```cpp
std::map<const Cluster *, int> map_cluster_num[3];   // L245
...
for (auto it = map_cluster_num[0].begin(); it != map_cluster_num[0].end(); it++) {
    if (it->second > max_value_u) { ... }             // L463-468
}
```

`std::map<const Cluster*, int>` is sorted by pointer address. The arg-max
loop uses strict `>`, so when two clusters have equal hit counts, the one
that appears *first* in address order wins — a different cluster across
runs, memory layouts, or compiler versions.

The same pattern is repeated in the "else" branch at L592-609.

**Fix:** Replace `std::map` with `std::unordered_map` and break ties
deterministically using the insertion index from `map_cluster_index`.
Factor the three duplicated arg-max scans into a single lambda to reduce
code. See §6 for the applied fix.

**Prototype:** Has the same non-determinism (`std::map<PR3DCluster*,int>`
iterated at prototype L207-224). The toolkit fix improves on the prototype.

### BUG-2 — `std::sort` tie-break on equal-length clusters (low severity)

**Location:** `clustering_deghost.cxx:180-182`

```cpp
std::sort(live_clusters.begin(), live_clusters.end(), [](const Cluster *a, const Cluster *b) {
    return a->get_length() > b->get_length();
});
```

`std::sort` is not stable. Two clusters of identical length will appear in
implementation-defined order, causing downstream ghost decisions to differ
across runs. The prototype has the same issue (also uses `std::sort`).

**Fix:** Change to `std::stable_sort`, or add a secondary sort key
(cluster ident, or position in `live_grouping.children()`). See §6.

### BUG-3 — Member type mismatch: `use_ctpc_` declared as `double` (bug, low severity)

**Location:** `clustering_deghost.cxx:102`

```cpp
double use_ctpc_{true};   // should be bool
double length_cut_{0};
```

`use_ctpc_` is documented and used as a boolean flag but is declared
`double`. This silently promotes `true` to `1.0`. The code still works
because any non-zero double is truthy in `if (use_ctpc_)`, but the type
declaration is wrong and misleading.

**Fix:** Change `double use_ctpc_{true};` to `bool use_ctpc_{true};`.

---

## 5. Efficiency

### E-1 — Duplicated arg-max scan (low severity)

Lines 463-480 and 592-609 are verbatim copies of the three-map arg-max
loop. Any future change to the selection criterion must be applied in
both places. Factor into a local helper (done in the BUG-1 fix, §6).

### E-2 — Per-point kd-tree queries: O(N·P·log M) — expected

For each cluster processed, every point issues up to 6 kd-tree NN queries
(3 planes × 2 clouds). This is the same scaling as the prototype and cannot
be improved without rearchitecting the global cloud accumulation. No action
needed.

### E-3 — `dv_wpids()` vs `wpids()` usage (informational)

The toolkit uses `live_grouping.dv_wpids()` at L127 (returns all (apa,face)
pairs known to the detector volume, keyed on `kAllLayers` WirePlaneIds) to
populate `wpid_params` and `af_dead_{u,v,w}_index`. Sibling stages such as
`ClusteringConnect1` use `grouping.wpids()` (only (apa,face) pairs present
in the live clusters). Using `dv_wpids()` is intentional here: it ensures
geometry params and dead-wire maps are ready for any (apa,face) a point
might belong to, even if no cluster currently spans that face. No action
needed; document the intent with a comment.

### E-4 — TODO: DetUtils helpers (low-effort cleanup)

In-file TODOs at L136-138 and L185-187 note that the `wpid_params`
construction block and the `DynamicPointCloud` creation could be replaced
by `DetUtils::face_parameters(dv)` + `DetUtils::make_dynamicpointcloud(dv)`
(`clus/src/DetUtils.cxx:15-49`). This would shorten the function by ~40
lines. Left as a follow-up; not addressed in this review.

---

## 6. Determinism

Summary of all pointer-keyed containers:

| Container | Line | Iterated? | Risk | Action |
|---|---|---|---|---|
| `map_cluster_num[3]` — `std::map<const Cluster*,int>` | 245 | Yes (L463-480, 592-609) | **HIGH — non-deterministic arg-max** | Fixed (see below) |
| `map_cluster_index` — `std::map<const Cluster*,int>` | 199 | No — only lookup | None | Optional: change to `unordered_map` |
| `live_clusters` — `std::vector<Cluster*>` | 177 | Iterated in length order | Low (equal-length tie-break) | Fixed with `stable_sort` |

**Applied fix for `map_cluster_num`:** Replaced with
`std::unordered_map<const Cluster*, int>` and introduced a deterministic
arg-max helper lambda that breaks ties by the cluster's index in
`map_cluster_index` (which is built from `live_grouping.children()` in
deterministic insertion order). Also changed `map_cluster_index` to
`std::unordered_map` to make the intent clear.

---

## 7. Multi-APA / multi-face correctness

### 7-A. How geometry is obtained per (apa,face)

`live_grouping.dv_wpids()` returns one `kAllLayers` WirePlaneId per
(apa,face). The loop at L140-170 runs once per face and correctly builds:

- `wpid_params[wpid_kAllLayers] = (drift_dir, angle_u, angle_v, angle_w)`  
  Used by `DynamicPointCloud` and `make_points_cluster{,_skeleton}`.
- `af_dead_{u,v,w}_index[apa][face]` — per-face dead wire maps.

`BlobSampler.cxx:337` confirms that per-point `wpid` values stored in the
point cloud are also `kAllLayers` WirePlaneIds, so `cluster->wire_plane_id(j)`
and the `wpid_params` keys are consistent.

### 7-B. Per-point apa/face routing

Every kd-tree query passes the point's `test_wpid.face()` and
`test_wpid.apa()` through to `get_closest_2d_point_info(point, plane,
face, apa)` (L273, 286, 332, 353, 386, 400). Dead-wire lookups use
`af_dead_u_index.at(test_wpid.apa()).at(test_wpid.face())` (L261, 306, 375).
Multi-face within one APA is handled correctly.

### 7-C. Multi-APA guard

The function enforces `apas.size() == 1` at L172-174 with `raise<ValueError>`.
The module comment at L118 says "This can handle entire APA (including all
faces)" but does not mention the multi-APA rejection. This is intentional:
the algorithm is applied globally over one APA's live clusters at a time.
Since `SBND` (and UBoone) have one APA, this is not limiting now.

The `visit()` method (L96-99) does not loop over APAs — it passes the whole
`live` grouping. If a future detector provides multi-APA groupings to this
function, the guard will throw. The needed change would be either:
(a) iterate per-APA in `visit()`, or (b) remove the guard after verifying
cross-APA cluster matching is safe.

**Updated comment in code** to clarify the single-APA restriction.

### 7-D. `merge_clusters` is APA-agnostic

The boost-CC merge step (L769) operates over the connectivity graph built
from all live clusters regardless of face/APA. Connected clusters across
faces can be merged, which is correct: a ghost spread across two faces
should be merged into the real cluster it ghosts.

### 7-E. `DynamicPointCloud::get_closest_2d_point_info` — correct

`DynamicPointCloud.cxx:283-321` creates a `kAllLayers` wpid from the query
`(face, apa)`, looks up the matching entry in `m_wpid_params`, selects the
correct wire angle for the requested plane, and builds a per-(apa,face)
per-layer kd-tree on first access. Fully multi-face aware.

---

## 8. Called-helper review

### `merge_clusters` (ClusteringFuncs.cxx:48-120)

- Uses `std::map<int, std::set<int>>` — both int-keyed. **Deterministic.**
- Children iterated via `auto orig_clusters = grouping.children()` — order
  matches insertion order. OK.
- Each CC creates a fresh `Cluster` with merged content, taking children from
  the original clusters and destroying originals via `destroy_child`. Correct.

### `extract_geometry_params` (ClusteringFuncs.cxx:12-46)

Single-face: `break`s after the first valid wpid. Used by `ClusteringConnect1`
(not by deghost, which hand-rolls its own loop). **Not a multi-APA concern.**

### `DynamicPointCloud::get_closest_2d_point_info` (DynamicPointCloud.cxx:283)

Correctly per-(apa,face). Returns `(distance, cluster*, point_index)`.
When the kd-tree for that (apa,face,plane) is empty (no points yet), returns
`(-1.0, nullptr, size_t(-1))`. **Callers must guard distance ≤ threshold,
which they do** (L277, 288 check distance against `dis_cut/3.` and
`dis_cut*2.0` respectively, so null pointer is never dereferenced directly).

### `make_points_cluster` / `make_points_cluster_skeleton` (DynamicPointCloud.cxx)

`make_points_cluster` (L449): looks up per-point `kAllLayers` wpid in
`wpid_params`; raises if not found. `flag_wrap=true` activates multi-face
wrapped 2D projections.

`make_points_cluster_skeleton` (L692): same pattern; also caches angles per
wpid ident to avoid repeated tuple unpacking. Correct.

### `Cluster::get_closest_points(other)` (Facade_Cluster.cxx)

Used 12 times in merge decisions (L491, 505, 522, 536, 553, 567, 635, 650,
665, 686, 706, 726). Returns the minimum 3D distance between two cluster
point-sets via kd-tree. Distance threshold of `20 * units::cm` is hardcoded
throughout, matching the prototype.

### `NaryTreeFacade::destroy_child` (NaryTreeFacade.h:239)

Takes `child_type*& kid` (reference), calls `remove_child`, then `kid = nullptr`.
The `assert(live == nullptr)` pattern in deghost (L776) and in `merge_clusters`
(ClusteringFuncs.cxx:107) is correct.

---

## 9. Action items → applied fixes

| ID | Severity | Status | Description |
|---|---|---|---|
| BUG-1 | High | **Fixed** | `map_cluster_num[3]` → `unordered_map`; deterministic arg-max lambda; duplicate scan factored out. |
| BUG-2 | Low | **Fixed** | `std::sort` → `std::stable_sort` for length-sorted cluster order. |
| BUG-3 | Low | **Fixed** | `double use_ctpc_{true}` → `bool use_ctpc_{true}`. |
| 3-B (parity) | Info | Documented | Short seed cluster added to skeleton cloud even with `use_ctpc=false`. Matches Overload B. Acceptable. |
| E-3 (dv_wpids) | Info | Documented | `dv_wpids()` usage is intentional; added clarifying comment. |
| E-4 (DetUtils TODO) | Info | Deferred | Replace hand-rolled wpid loop with `DetUtils::face_parameters` + `make_dynamicpointcloud`. |
| 7-C (multi-APA comment) | Info | **Fixed** | Module comment at L118 updated to state single-APA restriction explicitly. |
