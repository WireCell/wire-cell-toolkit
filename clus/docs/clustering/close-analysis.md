# ClusteringClose Analysis

Comparison of `clus/src/clustering_close.cxx` (toolkit) against
`prototype_base/2dtoy/src/ToyClustering_close.h` (prototype).

---

## Functional Equivalence

The core algorithm (`Clustering_3rd_round`) is identical between toolkit and prototype:

1. Find closest points `p1`, `p2` via `Find_Closest_Points`
2. If `dis < 0.5 cm` → merge unconditionally
3. If `dis < 1.0 cm` and both clusters short (< 12 cm) → merge
4. If `dis < 2.0 cm` and at least one cluster long (≥ 12 cm): compute Hough
   directions, dipole counts, nearby-point counts; merge if endpoint is
   at the extremity of a long straight track
5. If `dis < length_cut` and at least one cluster long: check average-position
   vectors and Hough directions for alignment → merge

The outer `clustering_close` function is equivalent: skip clusters < 1.5 cm,
skip already-absorbed clusters, mark short (< 5 cm) clusters as used after first
merge to prevent re-absorption.

**One toolkit improvement over prototype** (not a bug):
The prototype leaves `num_p1`, `num_p2`, `num_tp1`, `num_tp2` uninitialized when
`dis ≥ 2 cm` (technically UB). The toolkit initializes them to `0` via brace-init.
With 0-values the `(num_p1 > 25 || ...)` guard prevents `dir1`/`dir2` from being
accessed before they are computed, making the behavior well-defined. Correct fix.

**Intentional divergences:**
- **D1 (merge mechanism)**: Prototype uses manual union-find with
  `std::set<std::pair<...>>`. Toolkit builds a `boost::adjacency_list` and calls
  `merge_clusters`. Results are equivalent.
- **D2 (flag_print)**: Prototype sets `flag_print = true` (enables timing output).
  Toolkit had `flag_print = false` — dead code, removed (see DC1 below).

---

## Bugs Found

None.

---

## Dead Code Removed (in this session)

### DC1 — `flag_print`, `ExecMon em`, `if (flag_print)` guards (`clustering_close.cxx`)

**Before:** Three dead lines at the top of `Clustering_3rd_round`:
```cpp
bool flag_print = false;   // always false — all dependent branches are dead
ExecMon em("starting");    // only referenced inside dead if (flag_print) blocks
```
Plus two guarded print calls that never fired:
```cpp
if (flag_print) std::cout << em("Hough Transform") << std::endl;
if (flag_print) std::cout << em("Get Number Points") << std::endl;
```
And the now-unused include:
```cpp
#include "WireCellUtil/ExecMon.h"
```

**Fix:** Removed all of the above. `ExecMon.h` include replaced by
`<unordered_map>` and `<unordered_set>` (needed for efficiency changes below).

---

## Efficiency and Randomness Improvements (in this session)

### E1/R1 — `used_clusters` changed to `unordered_set`

**Before:** `cluster_set_t used_clusters` where `cluster_set_t = std::set<const Cluster*>`.
O(log N) `find` and `insert`; ordered by pointer address (ASLR-dependent).

**Fix:**
```cpp
std::unordered_set<const Cluster*> used_clusters;
```
O(1) average `find` and `insert`; no pointer-address ordering.

### E2/R2 — `map_cluster_index` changed to `unordered_map`

**Before:**
```cpp
std::map<const Cluster*, int> map_cluster_index;
```
O(log N) lookup; key ordering by pointer address.

**Fix:**
```cpp
std::unordered_map<const Cluster*, int> map_cluster_index;
```
O(1) average lookup; no pointer-address dependence.

### R3 — Sort `live_clusters` for deterministic outer loop order

**Before:** `const auto& live_clusters = live_grouping.children()` — reference to
grouping children in insertion order (depends on prior merge-step allocation).

**Fix:** Sorted copy before the main loops:
```cpp
auto live_clusters = live_grouping.children();  // sorted copy for deterministic order
sort_clusters(live_clusters);
```

**Impact:** The outer-loop order controls which short cluster (< 5 cm) gets absorbed
first into `used_clusters`. Sorting ensures consistent behavior across runs.

---

## Notes

- `ilive2desc` was already `std::unordered_map<int, int>`.
