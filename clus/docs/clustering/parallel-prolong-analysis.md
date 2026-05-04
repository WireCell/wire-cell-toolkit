# ClusteringParallelProlong Analysis

Comparison of `clus/src/clustering_parallel_prolong.cxx` (toolkit) against
`prototype_base/2dtoy/src/ToyClustering_para_prol.h` (prototype).

---

## Functional Equivalence

The core algorithm (`Clustering_2nd_round`) is identical between toolkit and prototype:

1. Bail out early if both clusters are < 10 cm
2. Find closest points `p1`, `p2` via `Find_Closest_Points`
3. If within distance threshold, compute local average positions and check:
   - **Parallel case**: connecting vector near-perpendicular to drift + average-position vector near-perpendicular to drift + Hough directions near-perpendicular to drift → check parallel-U and parallel-V sub-cases
   - **Prolonged case**: connecting vector (YZ projection) near-parallel to one of U/V/W wire directions → check Hough directions for alignment
4. Return false if none of the merge criteria fire

The outer `clustering_parallel_prolong` function is equivalent: all-pairs loop `j=i+1`
over live clusters, edges added to boost graph, `merge_clusters` for transitive closure.

**Intentional divergences (not bugs):**

- **D1 (multi-detector)**: Wire directions (`U_dir`, `V_dir`, `W_dir`) and angles
  are read from `dv->wire_direction()` per APA/face via `wpid_U_dir` / `wpid_V_dir` /
  `wpid_W_dir` maps. Prototype hardcodes MicroBooNE angles:
  `U_dir(0,cos60°,sin60°)`, `V_dir(0,cos60°,-sin60°)`, `W_dir(0,1,0)`.

- **D2 (removed dead prototype variables)**: Prototype declares `flag_para_U` and
  `flag_para_V` (lines 159-160 of prototype) but assigns them only once and never reads
  them for any decision. Toolkit omits them entirely — correct.

- **D3 (merge mechanism)**: Prototype uses a manual union-find with
  `std::set<std::pair<PR3DCluster*, PR3DCluster*>>`. Toolkit builds a
  `boost::adjacency_list` and calls `merge_clusters`. Results are equivalent.

---

## Bugs Found

None.

---

## Dead Code Note

`wpid_params` (line 271) and `apas` (line 275) are filled by `compute_wireplane_params`
but never read afterward in this function. Unlike the inline loops removed from
`clustering_extend.cxx` and `clustering_regular.cxx`, these are now output parameters
of the shared helper `compute_wireplane_params`. The caller signature cannot be changed
without affecting all other callers (some of which may read these values). Not changed.

---

## Efficiency and Randomness Improvements (in this session)

### E1/R1 — `map_cluster_index` changed to `unordered_map`

**Before:**
```cpp
std::map<const Cluster*, int> map_cluster_index;
```
O(log N) lookup; key ordering is by pointer address (ASLR-dependent).

**Fix:**
```cpp
std::unordered_map<const Cluster*, int> map_cluster_index;
```
O(1) average lookup; no pointer-address dependence.

### R2 — Sort `live_clusters` for deterministic outer loop order

**Before:** `const auto& live_clusters = live_grouping.children()` — reference to
grouping children in insertion order, which depends on prior merge-step allocation.

**Fix:** Sorted copy before the main loops:
```cpp
auto live_clusters = live_grouping.children();  // sorted copy for deterministic order
sort_clusters(live_clusters);
```

**Note:** As in `clustering_regular`, the all-pairs inner loop (`j=i+1`) visits every
pair regardless of order. The sort ensures consistent vertex numbering in the boost
graph and consistent `map_cluster_index` assignments across runs.

---

## Notes

- The `wpid_U_dir` / `wpid_V_dir` / `wpid_W_dir` maps remain as
  `std::map<WirePlaneId, ...>` — keyed by `WirePlaneId` (APA/face/layer struct), not
  pointer. Their ordering is deterministic and hardware-topology-based. No change needed.
- `ilive2desc` was already `std::unordered_map<int, int>`.
- The setup block uses `compute_wireplane_params()` (extracted shared helper in
  `Clustering_Util.cxx`) instead of the inline loop present in `clustering_extend.cxx`
  and `clustering_regular.cxx` before their cleanup. This is a pre-existing refactor.
