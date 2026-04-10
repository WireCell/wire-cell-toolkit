# Review: `clustering_connect1()` / `ClusteringConnect1`

**Toolkit file:** `clus/src/clustering_connect.cxx`  
**Prototype reference:** `prototype_base/2dtoy/src/ToyClustering_connect.h`  
**Reviewer:** Xin Qian (via Claude Code)  
**Date:** 2026-04-10

---

## A. Logic fidelity — overall faithful

Section-by-section comparison confirms the toolkit and prototype produce the same control flow,
distance cuts, merge conditions, and second-pass PCA alignment check.

API-level differences (intentional, not bugs):

| Prototype | Toolkit |
|-----------|---------|
| `TVector3` | `geo_point_t` |
| `cluster_length_map[cluster]` | `cluster->get_length()` (cached) |
| Hard-coded 60° U/V/W angles | `extract_geometry_params(grouping, dv)` — improvement, still reproduces MicroBooNE values |
| `global_cloud.pts[k].index_u` repurposed as `dis_cut` (integer) | `DPCPoint::dist_cut[plane]` (double, `std::floor`-rounded to same integer values) |
| `std::set<std::pair<...>>` + hand-rolled union-find (×2) | `boost::adjacency_list` + `merge_clusters()` helper (×2) — functionally equivalent |
| `cluster->Calc_PCA()` called inside inner loop | `cluster->get_pca()` reads cached result |
| `cluster->Create_point_cloud()` per cluster | Point cloud built eagerly in Facade, no explicit call needed |

### Prototype quirk faithfully preserved

Toolkit lines 364–367 (after this review) reproduce prototype line 225–226:
```cpp
if (cluster->get_length() < 100 * units::cm ||
    fabs(dir2.angle(...)) < 5 * ... &&
    fabs(dir1.angle(...)) < 5 * ... &&
    cluster->get_length() < 200 * units::cm)
```
C++ `&&` binds tighter than `||`, so this parses as `A || (B && C && D)`, which is the intended
semantics: enter the merge-test block for short clusters, OR for drift-parallel long clusters up to
200 cm. **Do not parenthesize to change behaviour.**  A comment was added to mark this.

---

## B. Bugs found and fixed

### B1. Nondeterministic tie-break in max-cluster selection (FIXED)

`map_cluster_num[0/1/2]` and `temp_clusters` were `std::map<const Cluster*, int>` /
`std::set<const Cluster*>`, which iterate in pointer-address order.  When two clusters had the
same overlap count the tie was resolved by allocator-dependent heap order, making results
non-reproducible across runs.

**Fix:** Changed to `std::map<const Cluster*, int, ClusterLess>` and
`std::set<const Cluster*, ClusterLess>` throughout.  `ClusterLess` orders by
`cluster->get_graph_index()`, which is stable across runs.  Applied to `map_cluster_dir1/dir2`
and `map_cluster_index` as well for consistency.

### B2. Double map lookup (FIXED)

Pattern `if (map.find(k) != end()) { ... map[k] ... }` causes two hash/tree lookups per access.
Replaced everywhere with a single `find()` and `it->second`.  Affects the dead-wire flag test
(3 planes × all points) and the cross-plane count lookups inside the max-finding loops.

### B3. Repeated `point3d(j)` / `point3d_raw(j)` calls (FIXED)

Each was called twice (or three times) per loop iteration.  Cached into
```cpp
const auto p3  = cluster->point3d(j);
const geo_point_t test_point(p3.x(), p3.y(), p3.z());
const double raw_x = cluster->point3d_raw(j).x();
```

### B4. Dead multi-face scaffolding (REMOVED)

`std::set<int> apas`, `af_dead_u/v/w_index` (3 nested maps), and their population inside the
`for (const auto& wpid : wpids)` loop were computed but **never read**.  They were removed to
avoid misleading readers into thinking multi-face is already supported.  The single-face
`dead_u/v/w_index` references (lines 74-76) and the `wpid_params` map remain, as both are used.

### B5. Scope-filtered clusters added to connectivity graphs (FIXED)

The loops that built `g` and `g2` added every child of `live_grouping`, including those that
fail `get_scope_filter(scope)`.  Added `if (!live->get_scope_filter(scope)) continue;` before
`boost::add_vertex` in both loops.  The orphaned isolated vertices were harmless but inflated the
graph and confused `merge_clusters`.

---

## C. Efficiency improvements

### C1. Three-plane code triplication eliminated

The inner per-point body was ~200 lines of three near-identical U/V/W blocks differing only in
plane index.  Replaced with a single `process_plane` lambda capturing `[&]`:

```cpp
auto process_plane = [&](int plane,
                         const std::map<int,std::pair<double,double>>& dead_index,
                         int wire_idx, double raw_x,
                         const geo_point_t& test_point,
                         int& num_unique_ref) { ... };

for (int j = 0; j != num_total_points; j++) {
    const auto p3 = cluster->point3d(j);
    const geo_point_t test_point(p3.x(), p3.y(), p3.z());
    const double raw_x = cluster->point3d_raw(j).x();
    process_plane(0, dead_u_index, winds[0][j], raw_x, test_point, num_unique[0]);
    process_plane(1, dead_v_index, winds[1][j], raw_x, test_point, num_unique[1]);
    process_plane(2, dead_w_index, winds[2][j], raw_x, test_point, num_unique[2]);
}
```

Logic preserved exactly: non-dead points use `skel_pts[gidx].dist_cut[plane]`;
dead points use `loose_dis_cut / 3. * 2.` — matching prototype behaviour.

### C2. Max-finding loops simplified

The three `for (it = map_cluster_num[k].begin() ...)` blocks for max-count detection were
refactored to range-for with structured bindings and a reusable `lookup_count` lambda, removing
the 6 double-lookup pairs.

---

## D. Determinism summary

After the B1 and B5 fixes the function is fully deterministic:

- `ClusterLess` ordering → stable map/set iteration
- `boost::connected_components` with integer-indexed adjacency_list → stable merge grouping
- `DynamicPointCloud::get_2d_points_info` KD-tree search → index-deterministic
- Scope-filter applied before graph-vertex assignment → consistent vertex numbering

---

## E. Multi-APA/face (multi-TPC) — known limitation

`clustering_connect1` is explicitly single-face (comment at line 54; `ValueError` raised at
line 61–66 if `wpids().size() > 1`).  Four load-bearing reasons this constraint must remain
until a deeper refactor:

1. **`extract_geometry_params`** (`ClusteringFuncs.cxx:12-46`) walks `wpids` and `break`s after
   the first entry — returns one angle triplet for the whole grouping.
2. **`make_points_linear_extrapolation`** (`DynamicPointCloud.cxx:855`) hard-codes
   `*(cluster->grouping()->wpids().begin())` and stamps that single wpid on every synthetic point.
3. **Dead-wire lookup**: `winds[0][j]` is a raw U wire index with no face/apa tag.  In a
   multi-face grouping a point from face B would be looked up in face A's `dead_u_index`,
   silently giving wrong dead-region hits.
4. **`vhough_transform`** pools all faces into one 3D KD tree — usually fine for direction
   estimation but imprecise if a cluster straddles faces.

The **ingest side** is already multi-face-correct:
- `make_points_cluster` reads `cluster->wire_plane_id(ipt)` per point and routes to the right
  per-face 2D KD tree.
- `DynamicPointCloud::get_2d_points_info(..., face, apa)` correctly filters the `(plane,face,apa)`
  KD tree.

So `DynamicPointCloud` is ready for multi-face, but `clustering_connect1` itself is not.
The `ValueError` guard is load-bearing — keep it until helpers 1–4 are updated.
The previously dead `af_dead_u/v/w_index` scaffolding was removed (B4) to avoid implying
multi-face support that does not exist.

---

## Callees reviewed

| Callee | Location | Finding |
|--------|----------|---------|
| `extract_geometry_params` | `ClusteringFuncs.cxx:12-46` | Picks first wpid, break — single-face; safe given caller's guard |
| `merge_clusters` | `ClusteringFuncs.cxx:48-120` | Boost `connected_components` + index grouping; correctly equivalent to prototype union-find |
| `DynamicPointCloud::add_points` | `DynamicPointCloud.cxx:81-187` | Multi-face-correct; per-wpid 2D KD trees |
| `DynamicPointCloud::get_2d_points_info` | `DynamicPointCloud.cxx:218-267` | Filters by `(plane,face,apa)` correctly |
| `make_points_cluster` | `DynamicPointCloud.cxx:449-516` | Multi-face-correct; reads `wire_plane_id(ipt)` per point |
| `make_points_linear_extrapolation` | `DynamicPointCloud.cxx:842-925` | Single-face hard-coded at line 855; safe given caller's single-face guard |
| `DynamicPointCloud::vhough_transform` | `DynamicPointCloud.cxx:441-445` | 3D KD tree pools all faces; acceptable for direction, not face-ambiguous in practice |
| `Cluster::wire_indices()` | `Facade_Cluster.cxx:1084-1091` | Raw wire index, no face tag; single-face-only valid |
| `ClusterLess` / `cluster_set_t` | `ClusteringFuncs.h:80-86` | Used here for determinism; `cluster_set_t` is an unordered alias — use explicit `ClusterLess` comparator |
