# ClusteringLiveDead Analysis

Comparison of `clus/src/clustering_live_dead.cxx` (toolkit) against
`prototype_base/2dtoy/src/ToyClustering_dead_live.h` (prototype).

---

## Functional Equivalence

The core algorithm is identical between toolkit and prototype:

1. Build `dead → [live clusters]` mapping via `is_connected` (overlap check)
2. For each dead cluster connecting ≥2 live clusters, test each pair
3. Convergence loop: alternating nearest-neighbour search to find closest points between two clusters
4. Compute Hough-transform directions, distance, angle metrics
5. Four merge-criteria branches based on cluster length (both-short / one-short / both-long)
6. Transitively merge all pairs decided for merging

**One toolkit enhancement not in prototype** (D1): Direction flip protection at the
Hough-transform step. After computing `dir1 = vhough_transform(mcell1_center, 30cm)`, the
toolkit also computes `dir5 = vhough_transform(p1, 30cm)` and flips `dir1` if the two
disagree by more than 120°. This prevents an ambiguous Hough result from producing a
backwards direction vector. The prototype uses `dir1` directly without this check.
This is a legitimate improvement; it is intentional divergence from the prototype.

**Prototype bug fixed in toolkit** (B-proto): In the `length_2 > 12cm && length_1 ≤ 12cm`
merge-criteria branch, the prototype has a copy-paste error:
```cpp
// prototype (wrong):
if ((dis <= 3*cm) && ((angle_diff2 <= 45 || angle_diff2<=45) ...
//                      ^^^^^^^^^^^^^^^^^^ same variable twice — should be angle_diff1
```
The toolkit correctly uses `(angle_diff1 <= 45 || angle_diff2 <= 45)`.

---

## Bugs Fixed (in this session)

### B1 — `cluster_less` PCA center self-comparison (`Facade_Cluster.cxx`)

**Lines 2558, 2560, 2562** (before fix):
```cpp
if (bc[0] < bc[0]) return false;   // always false — dead code
if (bc[1] < bc[1]) return false;   // always false — dead code
if (bc[2] < bc[2]) return false;   // always false — dead code
```
Should be `bc[i] < ac[i]`. The PCA center tie-breaking step of `cluster_less` was
completely non-functional, causing both `sort_clusters` (used everywhere) and
`ClusterLess` comparators to fall through to the pointer-address fallback
(`return a < b`) more often than necessary.

**Fix applied:**
```cpp
if (bc[0] < ac[0]) return false;
if (bc[1] < ac[1]) return false;
if (bc[2] < ac[2]) return false;
```

**Impact:** Affects all `sort_clusters` calls across all clustering methods. Reduces
frequency of pointer-based ordering (which varies with ASLR between runs).

---

## Efficiency and Randomness Improvements (in this session)

### R2 / B2 — Inconsistent sort of live clusters (`clustering_live_dead.cxx`)

**Before:** `live_clusters` was sorted with a simple lambda (sort key: `get_length()`
only; unstable), while `dead_clusters` used `sort_clusters()` which has full
tie-breaking (length → nblobs → npoints → wire ranges → PCA center → pointer).

```cpp
// before — unstable, single-key sort:
std::sort(live_clusters.begin(), live_clusters.end(), [](const Cluster* a, const Cluster* b) {
    return a->get_length() > b->get_length();
});
// sort_clusters(live_clusters);   // was commented out
```

**Fix:** Use `sort_clusters` for live clusters, consistent with dead clusters.
```cpp
sort_clusters(live_clusters);
```

### E1 / R3 — `tested_pairs` set replaced with index-based unordered_set

**Before:** `std::set<std::pair<const Cluster*, const Cluster*>>` — O(log N) lookup,
comparison by pointer address.

Two insertions per tested pair `(c1,c2)` were needed:
```cpp
tested_pairs.insert({c1, c2});
tested_pairs.insert({c2, c1});   // prevent reverse pair from being re-tested
```

**Fix:** Build an integer index for each cluster (`map_cluster_index`), then store a
single symmetric key in an `unordered_set<size_t>`:
```cpp
std::unordered_map<const Cluster*, int> map_cluster_index;  // O(1) lookup
// ...
std::unordered_set<size_t> tested_pairs;
auto pair_key = [&](const Cluster* a, const Cluster* b) -> size_t {
    int ia = map_cluster_index.at(a);
    int ib = map_cluster_index.at(b);
    if (ia > ib) std::swap(ia, ib);         // normalize: key(a,b) == key(b,a)
    return (size_t)ia * nlive + ib;
};
// One insert, one lookup:
if (tested_pairs.insert(pair_key(c1, c2)).second) { ... }
```

Benefits:
- O(1) average lookup instead of O(log N)
- Halves the number of set entries (symmetric key, single insertion per pair)
- Removes pointer-address dependence from the data structure entirely
- `map_cluster_index` changed from `std::map` (O(log N)) to `std::unordered_map` (O(1))

### E2 — Collapsed double map lookup for `dead_live_cluster_mapping`

**Before:** Two map traversals per connection found:
```cpp
if (dead_live_cluster_mapping.find(dead) == dead_live_cluster_mapping.end()) { ... }
dead_live_cluster_mapping[dead].push_back(live);
```

**Fix:** Use `emplace` return value:
```cpp
auto [it, inserted] = dead_live_cluster_mapping.emplace(dead, std::vector<Cluster*>{});
if (inserted) { dead_cluster_order.push_back(dead); }
it->second.push_back(live);
```

---

## Remaining notes

- **D2 (`get_scope_filter` skip)**: `live->get_scope_filter(m_scope)` returns false only
  if the cluster's "3d" xyz scope has been explicitly disabled. By default all clusters
  have this scope enabled (`m_scope_3d_raw = {"3d", {"x","y","z"}}`). In the per-face
  pipeline, `live_dead` always runs before any scope-switching, so the check is
  effectively always true. Not changed.

- **Pointer fallback in `cluster_less`**: After fixing B1, the PCA center comparison
  is now active. A true pointer-address last resort remains at line 2572 as intended
  ("randomness is the better choice" per the existing comment), but will now be
  reached much less frequently.
