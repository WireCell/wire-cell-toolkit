# ClusteringExtend Analysis

Comparison of `clus/src/clustering_extend.cxx` (toolkit) against
`prototype_base/2dtoy/src/ToyClustering_extend.h` (prototype).

---

## Functional Equivalence

The core algorithm is identical between toolkit and prototype:

1. Compute `length_1_cut` based on `flag` and `num_try`
2. Outer loop over all clusters: skip if `length < length_1_cut`
3. Four modes controlled by `flag`:
   - **flag=1** (prolong): find earliest/latest points, check if direction is near-parallel to a wire plane, call `Clustering_4th_prol` for each candidate
   - **flag=2** (parallel): find highest/lowest points, check if direction is near-perpendicular to drift, call `Clustering_4th_para`
   - **flag=3** (regular): find highest/lowest and earliest/latest extremes, call `Clustering_4th_reg` from each endpoint
   - **flag=4** (dead): if cluster has `live_dead` flag set, call `Clustering_4th_dead` against all candidates
4. Build boost graph of merge edges; call `merge_clusters` for transitive closure

**Intentional divergences from prototype (not bugs):**

- **D1 (multi-detector)**: Wire angles (`angle_u`, `angle_v`, `angle_w`) are read from
  `dv->wire_direction()` per APA/face, stored in `wpid_U_dir / wpid_V_dir / wpid_W_dir` maps.
  Prototype hardcodes MicroBooNE angles: `U_dir(0,cos60°,sin60°)`, `V_dir(0,cos60°,-sin60°)`, `W_dir(0,1,0)`.
  This is a correct generalization for PDHD.

- **D2 (detector-specific region removed)**: The early-return in `Clustering_4th_dead` that checks
  hardcoded z/y coordinate ranges (a MicroBooNE dead-channel region) is commented out in the toolkit.
  Correct — those coordinates are specific to one detector configuration.

- **D3 (equivalent refactor)**: Prototype collects merge pairs in
  `std::set<std::pair<PR3DCluster*, PR3DCluster*>>` and applies a manual union-find.
  Toolkit builds a `boost::adjacency_list` and calls `boost::connected_components` via `merge_clusters`.
  Results are equivalent; toolkit approach is more robust.

---

## Bugs Found

None. The missing `used_clusters.find()` check in the flag=2 lowest-point inner loop exists
identically in both prototype and toolkit — not a toolkit-specific bug.

---

## Efficiency and Randomness Improvements (in this session)

### R1 — Sort `live_clusters` for deterministic outer loop order

**Before:** `live_clusters` was a const-reference to `live_grouping.children()`, whose order
depends on the order clusters were inserted into the grouping after the previous merge step.

```cpp
const auto& live_clusters = live_grouping.children();
```

**Fix:** Take a sorted copy before the main loops:
```cpp
auto live_clusters = live_grouping.children();  // sorted copy for deterministic order
sort_clusters(live_clusters);
```

**Impact:** The outer-loop order determines which cluster absorbs which short cluster into
`used_clusters`. With an unsorted order, results could vary between runs depending on
allocation order of prior clustering steps.

### E1/R2 — `map_cluster_index` changed to `unordered_map`

**Before:**
```cpp
std::map<const Cluster*, int> map_cluster_index;
```
O(log N) lookup; key ordering is by pointer address (ASLR-dependent, though in this
context only used for lookup, not iteration).

**Fix:**
```cpp
std::unordered_map<const Cluster*, int> map_cluster_index;
```
O(1) average lookup; no pointer-address dependence.

### E2/R3 — `used_clusters` changed to `unordered_set`

**Before:** `cluster_set_t used_clusters` where `cluster_set_t = std::set<const Cluster*>`.
O(log N) `find` and `insert`; ordered by pointer address.

**Fix:**
```cpp
std::unordered_set<const Cluster*> used_clusters;
```
O(1) average `find` and `insert`; no pointer-address ordering.

---

## Notes

- The `wpid_*_dir` maps (`std::map<WirePlaneId, ...>`) are keyed by `WirePlaneId` (a struct
  with APA/face/layer fields), not by pointer. Their ordering is deterministic and
  hardware-topology-based. No change needed.
- The `ilive2desc` map was already `std::unordered_map<int, int>` before this session.
- `ClusteringExtendLoop` wraps `clustering_extend` in a loop with specific flag/length_cut
  combinations. Its `num_try` check (reduces to 1 try for >1100 clusters) is also present
  in the prototype's caller (`ToyClustering.cxx` lines ~293–330).
