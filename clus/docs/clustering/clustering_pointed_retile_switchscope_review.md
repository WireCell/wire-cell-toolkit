# Review: ClusteringPointed, ClusteringRetile, ClusteringSwitchScope

**Date:** 2026-04-10  
**Files examined:**
- `clus/src/clustering_pointed.cxx`
- `clus/src/clustering_retile.cxx`
- `clus/src/clustering_switch_scope.cxx`
- `clus/src/Facade_Cluster.cxx` (especially `add_corrected_points`, `set_default_scope`, `default_scope_from`, `separate`-related methods)
- `clus/src/Facade_Grouping.cxx` (especially `separate`, `from`)
- `clus/inc/WireCellClus/Facade_Grouping.h`
- `clus/inc/WireCellClus/Facade_Cluster.h`
- `clus/inc/WireCellClus/IPCTransform.h`
- `clus/inc/WireCellClus/ClusteringFuncsMixins.h`
- `util/inc/WireCellUtil/NaryTree.h` (`separate`, `children`)

These three components have no prototype counterparts — they are toolkit-native. They form a preprocessing pipeline that operates on the "live" grouping before downstream clustering:

```
ClusteringPointed → ClusteringSwitchScope → ClusteringRetile
```

---

## 1. ClusteringPointed

### Purpose
Removes empty blobs (no points) from clusters and then removes empty clusters (no children) from a configured set of groupings. Effectively a cleanup/pruning pass before further processing.

### Logic summary
```
for grouping_name in m_groupings:
    collect doomed_blobs per cluster
    destroy doomed_blobs
    if cluster has no children → collect for doom
    destroy doomed clusters
```

### Bugs / Correctness Issues

**None found.** The two-pass pattern (collect then destroy) is safe: `grouping->children()` and `cluster->children()` both return vectors by value (copies, per `NaryTree::Node::children()`), so modifying the tree during the outer loop does not invalidate iterators.

**Minor: silent skip on missing grouping.** If a configured grouping name is not in the ensemble, the code silently continues (line 47: `if (got.empty()) { continue; }`). This is reasonable but worth noting — a misconfigured grouping name produces no warning.

### Efficiency

- **OK.** Both `grouping->children()` and `cluster->children()` are called once per level and their return vectors are iterated linearly. No caching concerns.
- `blob->npoints()` on a leaf blob returns the size of the local `"3d"` point cloud array — this is O(1) and correct.

### TPC / APA / Face handling

Not applicable. `ClusteringPointed` is geometry-agnostic — it operates purely on tree structure (blobs with points vs without). No APA, face, or wire-plane information is needed or used.

---

## 2. ClusteringSwitchScope

### Purpose
Applies a named coordinate correction (currently only `"T0Correction"`) to each cluster: the corrected x-coordinate (`x_t0cor`) is appended to each blob's `"3d"` point cloud, and then each cluster is split into two sub-clusters based on a per-blob filter (did any corrected point fall inside the detector active volume?). The new "correction scope" (`{"3d", {"x_t0cor", "y", "z"}}`) is then set as the default scope of all resulting clusters.

### Logic summary (T0Correction path)
```
for cluster in live_grouping.children() [snapshot]:
    filter_results = cluster->add_corrected_points(pcts, "T0Correction")
        # adds x_t0cor to each blob's 3d PC; returns 0/1 per blob
    correction_scope = cluster->get_scope("T0Correction")   # {"3d",{"x_t0cor","y","z"}}
    cluster->set_default_scope(correction_scope)
    cluster->set_scope_transform(correction_scope, "T0Correction")
    separated = live_grouping.separate(cluster, filter_results, remove=true)
        # → {0: cluster_fail, 1: cluster_pass}  (or only one of these)
    for id, new_cluster in separated:
        new_cluster->set_default_scope(correction_scope)   # redundant
        new_cluster->set_scope_filter(correction_scope, id==1)
        new_cluster->set_scope_transform(correction_scope, "T0Correction")  # redundant
```

### Bugs / Correctness Issues

**B1 — Dead parameter: `default_scope` is never used.**  
The static helper `clustering_switch_scope()` accepts a `const Tree::Scope& default_scope` parameter (line 55 of `clustering_switch_scope.cxx`) but the parameter is never referenced anywhere in the function body. The NeedScope mixin is configured, and `m_scope` is passed to the helper, but it is silently ignored. If the intent was to use it as the pre-correction scope for reference or validation, that logic is missing.

**B2 — Silent no-op for unknown correction names.**  
The outer `if (correction_name == "T0Correction")` block in the static helper has no else branch. If any correction name other than `"T0Correction"` is configured, the entire per-cluster loop body is skipped silently. The function returns with the live grouping unchanged. Note that `add_corrected_points()` itself *does* raise a `RuntimeError` for unknown names — but it is never reached because the check occurs first.  
**Recommendation:** add a `raise<RuntimeError>` or at minimum a `spdlog::warn` for unrecognised `correction_name` values.

**B3 — Dead code: unused `filter_result_set` (lines 102–106).**  
```cpp
std::set<int> filter_result_set(filter_results.begin(), filter_results.end());
```
This set is constructed but never read after its construction. It was used by commented-out `info()` calls. The object construction is harmless but wastes CPU for every cluster. Remove it.

**B4 — Redundant calls inside the post-separation loop.**  
After `Grouping::separate()` calls `c->from(*cluster)` (which calls `default_scope_from()` and thus `set_default_scope(correction_scope)` and `set_scope_transform(correction_scope, ...)`), the explicit loop at lines 119 and 124 calls those exact same setters again. `set_default_scope` calls `clear_cache()` every time it is called — so the cache is cleared twice per separated cluster. Only the `set_scope_filter()` call at lines 121/123 is non-redundant.  
**Recommendation:** simplify the post-separation loop to only set the scope filter:
```cpp
for (auto& [id, new_cluster] : separated_clusters) {
    new_cluster->set_scope_filter(correction_scope, id == 1);
}
```

**B5 — `add_corrected_points`: `children()` called O(N) times in the loop.**  
In `Facade_Cluster.cxx` lines 206–222:
```cpp
for (size_t iblob = 0; iblob < this->children().size(); ++iblob) {
    Blob* blob = this->children().at(iblob);
```
`NaryTree::Node::children()` constructs and returns a new `std::vector` of raw pointers on every call (see `NaryTree.h` lines 391–400). The loop above calls `children()` **twice per iteration**: once for `size()` in the condition, once for `at(iblob)` in the body. This allocates and fills a new vector 2×N_blobs times per cluster.  
**Recommendation:** hoist the children vector out of the loop:
```cpp
const auto blobs = this->children();
blob_passed.resize(blobs.size(), 0);
for (size_t iblob = 0; iblob < blobs.size(); ++iblob) {
    Blob* blob = blobs[iblob];
    ...
}
```

**B6 — `add_corrected_points`: scope name is hardcoded as `"T0Correction"` inside the function.**  
Line 225 of `Facade_Cluster.cxx`:
```cpp
m_scopes["T0Correction"] = {"3d", {"x_t0cor", "y", "z"}};
```
Even though the function receives `correction_name` as a parameter, the scope is registered under the literal string `"T0Correction"`. If a caller passes a different correction name that reaches the T0Correction branch (not currently possible but a maintenance hazard), the scope would be registered under the wrong key.  
**Recommendation:** use `m_scopes[correction_name] = ...`.

### Efficiency

- Remove `filter_result_set` (B3): removes one `std::set<int>` construction per cluster.
- Hoist `children()` in `add_corrected_points` (B5): eliminates 2×N_blobs vector allocations per cluster.
- The redundant `set_default_scope` / `set_scope_transform` calls (B4) each clear the entire cluster cache. Removing them avoids O(redundant × N_clusters) cache clears.

### TPC / APA / Face handling

**Correct.** `add_corrected_points` applies the T0 correction on a **per-blob basis**, using each blob's `wpid().face()` and `wpid().apa()`. This correctly handles multi-APA and multi-face clusters: each blob gets the correction appropriate to its own wire-plane ID.

**Coarseness of filter per blob.** The filter result per blob is 1 if *any* point in the blob's corrected point cloud falls inside the detector active volume, 0 otherwise. This is intentionally coarse. A blob that is mostly outside the volume but has one edge point inside will be placed in the "pass" sub-cluster. This is a design choice and should be understood by callers.

**`y_t0cor` and `z_t0cor` discarded.** The `pct->forward()` call produces `{"x_t0cor", "y_t0cor", "z_t0cor"}` but only `x_t0cor` is added to the blob's `"3d"` PC (line 212). This is intentional for a drift-direction (x) T0 correction — y and z coordinates are unchanged. However, the `filter()` call still receives `{"x_t0cor", "y_t0cor", "z_t0cor"}` from the forward-transform output, which is correct because `y_t0cor == y` and `z_t0cor == z` for a T0 correction.

**The scope assigned** (`{"3d", {"x_t0cor", "y", "z"}}`) uses the corrected x and the original y and z. All downstream operations that use the default scope will use `x_t0cor` for the x coordinate. This is the intended behavior.

---

## 3. ClusteringRetile

### Purpose
For each cluster in "live" that passes the configured scope filter (i.e., `get_scope_filter(m_scope) == true`), delegates retiling to an `IPCTreeMutate` component. Retiled clusters are inserted into a newly created "shadow" grouping. Non-retiled clusters (failing the filter or the mutate) are not placed in shadow.

### Logic summary
```
orig_grouping = ensemble["live"]
shadow_grouping = ensemble.make_grouping("shadow")
shadow_grouping.from(orig_grouping)    # copies anodes + dv, not children

for cluster in orig_grouping.children():
    if not cluster->get_scope_filter(m_scope): continue
    if cluster->get_default_scope() != m_scope:
        cluster->set_default_scope(m_scope)
    shad_node = m_retiler->mutate(*cluster->node())
    if shad_node:
        shad_root->insert(shad_node)
```

### Bugs / Correctness Issues

**B7 — Scope hash mismatch with ClusteringSwitchScope output (critical).**  
`NeedScope` initialises `m_scope` to `{"3d", {"x", "y", "z"}}` by default. However, after `ClusteringSwitchScope` runs, `set_scope_filter()` is set on each cluster for the *correction scope* `{"3d", {"x_t0cor", "y", "z"}}`. These two scopes have **different hashes**. Therefore, unless `ClusteringRetile` is explicitly configured with:
```json
{ "pc_name": "3d", "coords": ["x_t0cor", "y", "z"] }
```
`get_scope_filter(m_scope)` will return `false` for every cluster (missing key → default false), and **no cluster will be retiled**. The shadow grouping will be empty.

This is not a code bug per se (the user is expected to configure the scope), but:
- The `default_configuration()` method does not include `pc_name` or `coords` defaults, giving no hint to users.
- The `NeedScope` default differs from the expected post-SwitchScope value.
- There is no warning or diagnostic when all clusters are filtered out.

**Recommendation:**
1. Add `pc_name` and `coords` to `default_configuration()`.
2. Add a spdlog::warn if `shadow_grouping.children().empty()` after the loop.

**B8 — "live" and "shadow" grouping names are hardcoded.**  
Line 84 has `// fixme: make grouping names configurable`. Both `"live"` and `"shadow"` are literals. This is a known issue already flagged by the developer. No other review action needed beyond noting it here.

**B9 — `ensemble.with_name("live").at(0)` crashes on missing grouping.**  
If no grouping named `"live"` exists in the ensemble, `with_name()` returns an empty vector and `.at(0)` throws `std::out_of_range`. The same issue exists in `ClusteringSwitchScope` (line 41).  
**Recommendation:** check the return value:
```cpp
auto live_vec = ensemble.with_name("live");
if (live_vec.empty()) {
    spdlog::warn("ClusteringRetile: no 'live' grouping found, skipping");
    return;
}
auto& orig_grouping = *live_vec.at(0);
```

**B10 — Non-retiled (filtered-out) clusters are absent from shadow.**  
Clusters where `get_scope_filter(m_scope) == false` (outside detector volume) and clusters where `m_retiler->mutate()` returns null are silently dropped from the shadow grouping. If downstream algorithms treat shadow as a complete replacement of live, these clusters are permanently lost. If this is intentional (shadow only contains "good" retiled clusters), it should be documented explicitly.

**B11 — `set_default_scope` on a live cluster mid-loop.**  
Lines 98–100 update `orig_cluster->set_default_scope(m_scope)` in-place on the live cluster. `set_default_scope` calls `clear_cache()`. If the loop also queries any cached values on the original cluster after this point (e.g., through the retiler examining the original cluster's cache), the cache has been cleared. Currently the retiler receives the node directly via `mutate(*orig_cluster->node())` and does not call back into the cluster facade, so this is probably harmless. But it is a side effect on the live grouping that is worth documenting.

### Efficiency

**E1 — `m_scope.hash()` recomputed every iteration.**  
`get_scope_filter(m_scope)` and the hash comparison `orig_cluster->get_default_scope().hash() != m_scope.hash()` both compute `m_scope.hash()` on every cluster. If `Scope::hash()` is not trivially O(1) (it likely involves hashing strings), precomputing it once before the loop would help:
```cpp
const auto scope_hash = m_scope.hash();
for (auto* orig_cluster : orig_grouping.children()) {
    if (!orig_cluster->get_scope_filter(m_scope)) continue;
    if (orig_cluster->get_default_scope().hash() != scope_hash) {
        orig_cluster->set_default_scope(m_scope);
    }
    ...
}
```

### TPC / APA / Face handling

`ClusteringRetile` itself does not perform any APA/face-specific logic — this is fully delegated to `m_retiler->mutate()`. The scope filter is at the cluster level (a single boolean per scope hash). Multi-APA clusters are handled by the retiler implementation.

---

## 4. Cross-cutting issues

### Scope chain dependency
These three components form a chain with strict scope dependencies:

| Stage                  | Scope used / set                           |
|------------------------|--------------------------------------------|
| ClusteringPointed      | N/A (geometry-free pruning)               |
| ClusteringSwitchScope  | reads original scope; sets `{"3d",{"x_t0cor","y","z"}}` as new default |
| ClusteringRetile       | must be configured with `{"3d",{"x_t0cor","y","z"}}` |

The scope value `{"3d", {"x_t0cor", "y", "z"}}` is a magic constant that must agree between `add_corrected_points` (line 225 of `Facade_Cluster.cxx`), the pipeline configuration of `ClusteringSwitchScope`, and the pipeline configuration of `ClusteringRetile`. It is currently not centralised.

**Recommendation:** define this scope in a shared header (e.g., `Facade_Cluster.h`) as a `static constexpr` or `inline` so all three points of use share a single definition.

### Stale-pointer safety after `separate()`
In `clustering_switch_scope`, after `live_grouping.separate(cluster, filter_results, true)`, the pointer `cluster` is invalid (the node has been destroyed by `destroy_child`). The code does not dereference `cluster` after this point (correction_scope was captured as a value copy before separate). This is correct but fragile — future modifications that add a `cluster->something()` call after line 110 would be a use-after-free. Consider adding a comment:
```cpp
auto separated_clusters = live_grouping.separate(cluster, filter_results, true);
// 'cluster' is destroyed at this point; do not dereference it after here.
```

### No logging / diagnostics
All three components produce no log output. `ClusteringSwitchScope` has many commented-out `info()` and `std::cout` calls — these should either be restored at `debug` level or permanently removed. Dead commented-out code accumulates and makes the codebase harder to read.

---

## Summary table

| # | Component / File | Severity | Issue |
|---|---|---|---|
| B1 | clustering_switch_scope.cxx:55 | Low | `default_scope` parameter never used |
| B2 | clustering_switch_scope.cxx:77 | Medium | Unknown correction names silently do nothing |
| B3 | clustering_switch_scope.cxx:102 | Low | `filter_result_set` constructed but never read |
| B4 | clustering_switch_scope.cxx:119,124 | Low | `set_default_scope`/`set_scope_transform` redundant after `from()`; causes extra `clear_cache()` |
| B5 | Facade_Cluster.cxx:206-207 | Medium | `children()` called 2× per loop iteration, allocating vectors each time |
| B6 | Facade_Cluster.cxx:225 | Low | Scope key hardcoded as `"T0Correction"` rather than using `correction_name` |
| B7 | clustering_retile.cxx:93 | **Critical** | Default scope `{"x","y","z"}` ≠ SwitchScope output `{"x_t0cor","y","z"}`; no cluster passes filter without explicit config |
| B8 | clustering_retile.cxx:84 | Low | Hardcoded grouping names (known fixme) |
| B9 | clustering_retile.cxx:86; clustering_switch_scope.cxx:41 | Medium | `.at(0)` crash if "live" grouping absent |
| B10 | clustering_retile.cxx:93-95 | Design | Filtered-out clusters silently absent from shadow — needs documentation |
| B11 | clustering_retile.cxx:99 | Low | `set_default_scope` clears live cluster cache in-place |
