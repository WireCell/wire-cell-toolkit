# ClusteringRegular Analysis

Comparison of `clus/src/clustering_regular.cxx` (toolkit) against
`prototype_base/2dtoy/src/ToyClustering_reg.h` (prototype).

---

## Functional Equivalence

The core algorithm (`Clustering_1st_round`) is identical between toolkit and prototype.
The outer `clustering_regular` function is also equivalent; the prototype wraps its inner
loop in `for (int kk=0;kk!=1;kk++)` (a single-iteration no-op) which the toolkit omits.

**Intentional divergences (not bugs):**

- **D1 (multi-detector)**: Wire directions (`U_dir`, `V_dir`, `W_dir`) are read from
  `dv->wire_direction()` per APA/face via `wpid_U_dir` / `wpid_V_dir` / `wpid_W_dir` maps.
  Prototype hardcodes MicroBooNE angles: `U_dir(0,cos60°,sin60°)`, `V_dir(0,cos60°,-sin60°)`,
  `W_dir(0,1,0)`.

- **D2 (PDHD cross-APA addition)**: Lines 283–286 add an extra merge condition for very
  long tracks (>100 cm each, <5 cm apart) spanning different APAs or faces. Not in prototype;
  intentional PDHD geometry fix.

---

## Bugs Found

None. The prototype declared but never used `angle5` / `angle5_1` variables — the toolkit
already removed them.

---

## Dead Code Removed (in this session)

### DC1 — `wpid_params`, `apas`, `drift_dir`, `face_dirx` (`clustering_regular.cxx`)

**Before:** Four declarations that fed only into `wpid_params`, which was never read:
```cpp
std::map<WirePlaneId , std::tuple<geo_point_t, double, double, double>> wpid_params;
std::set<int> apas;
...
for (const auto& wpid : wpids) {
    apas.insert(apa);
    int face_dirx = dv->face_dirx(wpid_u);
    geo_point_t drift_dir(face_dirx, 0, 0);
    ...
    wpid_params[wpid] = std::make_tuple(drift_dir, angle_u, angle_v, angle_w);
```

**Fix:** Removed `wpid_params` declaration and its only assignment, `apas` and its
`insert`, `face_dirx` and `drift_dir` (inside the for loop) that fed nothing else.

### DC2 — Same dead code in `clustering_extend.cxx`

Identical `wpid_params` / `apas` / `drift_dir` / `face_dirx` dead-code block was also
present in `clustering_extend.cxx` and removed in the same session.

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

**Before:** `const auto& live_clusters = live_grouping.children()` — reference to grouping
children in insertion order, which depends on prior merge-step allocation.

**Fix:** Sorted copy before the main loops:
```cpp
auto live_clusters = live_grouping.children();  // sorted copy for deterministic order
sort_clusters(live_clusters);
```

**Note:** Unlike `clustering_extend` where `used_clusters` made the result directly
order-dependent, here the all-pairs inner loop (`j=i+1`) visits every pair regardless of
order. The sort ensures consistent vertex numbering in the boost graph and consistent
`map_cluster_index` assignments across runs.

---

## Notes

- The `wpid_U_dir` / `wpid_V_dir` / `wpid_W_dir` maps remain as `std::map<WirePlaneId, ...>`
  — keyed by `WirePlaneId` (APA/face/layer struct), not pointer. Their ordering is
  deterministic and hardware-topology-based. No change needed.
- The `ilive2desc` map was already `std::unordered_map<int, int>` before this session.
