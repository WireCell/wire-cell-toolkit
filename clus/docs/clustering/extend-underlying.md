# ClusteringExtend / ClusteringExtendLoop — Underlying Function Audit

Functions called by `clustering_extend` and `ClusteringExtendLoop::visit` in
`clus/src/clustering_extend.cxx`.

---

## Underlying Functions

| Function | Source | Called where |
|---|---|---|
| `Find_Closest_Points` | `clustering_extend.cxx` (defined at bottom) | multiple helpers |
| `get_strategic_points` | `clustering_extend.cxx` | inside `Find_Closest_Points` |
| `get_hull` | `Facade_Cluster.cxx` | inside `get_strategic_points` |
| `get_closest_point_blob` | `Facade_Cluster.cxx` | inside `Find_Closest_Points` + helpers |
| `sort_clusters` | `Facade_Cluster.cxx` | line 473 |
| `merge_clusters` | `ClusteringFuncs.cxx` | line 504 |

---

## Randomness Issues Found and Fixed (this session)

### R1 — `get_strategic_points`: `std::set<std::pair<geo_point_t, const Blob*>>`

**Before:**
```cpp
std::set<std::pair<geo_point_t, const Blob*>> unique_points;
// ...
unique_points.emplace(p, blob);
// ...
std::vector<...> points;
points.insert(points.end(), unique_points.begin(), unique_points.end());
```
`std::pair` comparison falls through to `const Blob*` pointer comparison when two entries
share the same `geo_point_t` (e.g., two hull points mapped to the same cluster point via
different blobs). Pointer ordering is ASLR-dependent. The resulting vector order fed into
`Find_Closest_Points` can differ run-to-run; with the early-exit at `dis < 0.5 cm` this
can alter `p1_save`/`p2_save`.

**Fix:** Replace with a `std::vector`, sort by `geo_point_t` only (using `operator<`,
which is deterministic content-based), then deduplicate via `std::unique`:
```cpp
std::vector<std::pair<geo_point_t, const Blob*>> points;
// ... emplace_back ...
std::sort(points.begin(), points.end(),
    [](const auto& a, const auto& b) { return a.first < b.first; });
points.erase(std::unique(points.begin(), points.end(),
    [](const auto& a, const auto& b) {
        return !(a.first < b.first) && !(b.first < a.first);
    }), points.end());
```
No pointer value is ever used for ordering or deduplication.

---

## Dead Code Removed (this session)

### DC1 — Commented-out original `Find_Closest_Points` implementation (~94 lines)

The original two-pass convergence algorithm (starting from `get_first_blob()` /
`get_last_blob()`) was replaced by the `get_strategic_points`-based multi-start approach
but left in place as a commented-out block after `return dis_save;`. Removed.

---

## Efficiency

No additional improvements found.

- `get_hull` uses `quickhull::QuickHull` with `std::set<int>` for index deduplication
  (integer keys, deterministic) — no issues.
- `get_closest_point_blob` wraps a KD-tree `kd_knn(1, point)` — O(log N), already
  optimal.
- `Find_Closest_Points` iterates over all strategic points with a convergence inner loop
  and early-exit at 0.5 cm. No further algorithmic improvement identified.

---

## `cluster_less` / `sort_clusters`

`cluster_less` has a pointer-comparison fallback (`return a < b`) as the absolute last
tiebreaker after exhausting all content-based comparisons (length, nblobs, npoints, wire
ranges, PCA center). This is intentional and documented in the source: pointer comparison
is preferred over returning `false` (which would incorrectly treat distinct clusters as
equal in ordered containers). No action needed.

---

## `merge_clusters`

- `desc2id`: `unordered_map<int,int>` — fine.
- `id2desc`: `map<int, set<int>>` — integer keys, ascending component-ID iteration is
  deterministic.
- Vertex numbering is driven by `sort_clusters`-ordered `live_clusters` → already fixed
  in prior session. `merge_clusters` is deterministic.
