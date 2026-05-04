# ClusteringSeparate — Analysis

Comparison of `clustering_separate.cxx` against the prototype
`ToyClustering_separate.h`. Covers correctness, bugs, randomness, and efficiency.

Prototype functions examined:
- `Clustering_separate` (lines 1827–2193, and 2358+)
- `JudgeSeparateDec_1` (lines 11–34)
- `JudgeSeparateDec_2` (lines 36–554)
- `Separate_1` (lines 556–1182, and ctpc variant 1183–1826)
- `Separate_2` (lines 2197–2355)

---

## Summary Table

| Issue | Location | Type | Severity |
|---|---|---|---|
| `del_clusters.push_back(cluster)` stores dangling pointer | `clustering_separate` line ~346 | Bug | Low (dead code, never dereferenced) |
| `new_clusters` / `del_clusters` never used after loop | `clustering_separate` lines ~102–370 | Dead code | Cosmetic |
| `scope_transform` retrieved but not applied after `Separate_2` | `clustering_separate` lines ~324, ~361 | Potential Bug | Medium |
| Sort lacks stable tiebreaker for equal-length clusters | `clustering_separate` line ~70 | Randomness | Low |
| `sqrt` in 30+ distance comparisons in `JudgeSeparateDec_2` | `JudgeSeparateDec_2` | Efficiency | Medium |

---

## 1. Functional Equivalence — Overview

The toolkit reproduces the prototype's three-level separation cascade faithfully:
```
Level 0:  Separate_1(cluster)           → cluster1 + cluster2
Level 1:  Separate_1(cluster2)          → cluster3 + cluster4  [if length>100cm]
Level 2:  Separate_1(cluster4)          → cluster5 + cluster6  [if length>100cm]
Final:    Separate_1(final_sep_cluster) → strip one more time  [if length>60cm]
Last:     Separate_2(final_sep_cluster)                        [remaining]
```
And the "stripping" path (`else if length < 6m`):
```
Separate_1(cluster) → cluster1 + cluster2
Separate_2(cluster2)
```

The toolkit uses `grouping->separate(cluster, b2id, true)` in-place rather than the prototype's
explicit `new_clusters` / `del_clusters` / `delete` management. This is architecturally
correct: the grouping owns the lifetime of clusters.

---

## 2. JudgeSeparateDec_1

### temp_angle1 computation

**Prototype:**
```cpp
double temp_angle1 = asin(cluster->get_num_time_slices()*time_slice_length/length)/3.1415926*180.;
```

**Toolkit:**
```cpp
auto points = cluster->get_earliest_latest_points();
double temp_angle1 = asin(fabs(points.first.x() - points.second.x()) / length) / 3.1415926 * 180.;
```

These are semantically equivalent: `num_time_slices * tick_width ≈ Δx`. The
toolkit uses the actual x-extent between the earliest and latest 3D points, which
is a marginally more precise estimate. **Not a bug.**

### PDHD guard (new addition)
```cpp
if (angle1 < 5 && ratio2 < 0.05 && cluster->get_length() > 300*units::cm) return false;
```
Intentional new guard for PDHD geometry (longer tracks). Not in the prototype. **Correct.**

---

## 3. JudgeSeparateDec_2

### flag_outx generalization (correct)
Prototype hardcodes `hx > 257cm || lx < -1cm` (MicroBooNE cathode/wire-planes).
Toolkit uses `det_FV_xmax + det_FV_xmax_margin` / `det_FV_xmin - det_FV_xmin_margin`
from detector metadata. **Correct generalization.**

### lx_points used inside hx_points loop (not a bug)
Both the prototype (lines 336–348) and the toolkit (lines 757–774) intentionally use
`lx_points.at(j)` for `independent_surfaces` classification inside the `hx_points` loop.
Since `hx_points` and `lx_points` are always exactly size 1 (never pushed, only at(0)
updated in-place), `j` is always 0 and `lx_points.at(0)` is always valid. The logic
reflects: when the highest-x point exits the detector, classify the surface using the
lowest-x point's coordinates (prototype design choice). **Faithfully reproduced. Not a bug.**

---

## 4. Separate_2 — Multi-APA Extension (correct)

The prototype has a single `time_cells_set_map` (one APA/face). The toolkit extends
this via `af_time_slices[apa][face]` to support multiple APAs/faces. The blob
connectivity logic (same-slice overlap + adjacent-slice lookahead ±1 or ±2) is
identical. The `mcell_index_map` uses `std::map<const Blob*, int>` but is accessed
only by pointer key lookup — its internal ordering does not affect results.
**Correct multi-APA extension.**

---

## 5. Bug — `del_clusters.push_back(cluster)` with Dangling Pointer

### Location
`clustering_separate`, lines ~198–202 (commented out) and ~346 (active).

### Problem
`grouping->separate(Cluster*& cluster, ...)` takes `cluster` by reference and sets it to
`nullptr` after destroying the cluster. This happens inside `Separate_1` for the local
parameter copy. The outer `cluster` pointer in `clustering_separate` is passed **by
value** to `Separate_1`, so it is NOT nulled out — it remains a dangling pointer to
the now-destroyed cluster object.

```cpp
// After Separate_1(cluster, ...) returns, outer `cluster` is dangling:
del_clusters.push_back(cluster);  // BUG: dangling pointer
```

Since `del_clusters` is never processed after the loop, no crash occurs. But the
code is incorrect and misleading.

### Fix
Remove `del_clusters.push_back(cluster)` at line 346 (and the entire `del_clusters`
and `new_clusters` mechanism — see next issue).

---

## 6. Dead Code — `new_clusters` / `del_clusters`

### Location
`clustering_separate`, lines ~102–370 (declarations and all push_backs).

### Problem
Both vectors are populated throughout the main loop but **never read after the loop
ends**. The grouping is modified in-place via `Separate_1` (which calls
`grouping->separate()` internally) and `live_grouping.separate()` (called directly
for `Separate_2` results). The prototype's explicit `new_clusters`/`del_clusters`
management was needed because the prototype owns raw pointers, but the toolkit
grouping owns cluster lifetimes automatically.

### Fix
Remove both vector declarations and all `push_back` calls.

---

## 7. Potential Bug — `scope_transform` Retrieved but Not Applied After `Separate_2`

### Locations
`clustering_separate` lines ~324 and ~361:
```cpp
auto scope_transform = final_sep_cluster->get_scope_transform(scope);
auto final_sep_clusters = live_grouping.separate(final_sep_cluster, b2id, true);
assert(final_sep_cluster == nullptr);
// scope_transform never used!
```

### Problem
After `Separate_2` + `live_grouping.separate()`, the resulting clusters have their
scope set via the grouping's internal inheritance. However, the `scope_transform`
retrieved from `final_sep_cluster` is never applied to any of the resulting clusters.
In `Separate_1`, the equivalent is done explicitly:
```cpp
clusters_step0[0]->set_scope_transform(scope, scope_transform);
```

Whether the `separate()` call automatically propagates `scope_transform` needs
verification. If it does not, the clusters from these two `Separate_2` calls may lack
a scope_transform. **Needs investigation.**

---

## 8. Randomness — Sort Lacks Stable Tiebreaker

### Location
`clustering_separate`, lines ~70–72:
```cpp
std::sort(live_clusters.begin(), live_clusters.end(), [](const Cluster *c1, const Cluster *c2) {
    return c1->get_length() > c2->get_length();
});
```

### Problem
When two clusters have equal `get_length()`, `std::sort` may order them differently
across runs (ASLR-dependent heap addresses, non-stable sort).

### Fix
Add pointer address as a deterministic tiebreaker:
```cpp
std::sort(live_clusters.begin(), live_clusters.end(), [](const Cluster *c1, const Cluster *c2) {
    if (c1->get_length() != c2->get_length()) return c1->get_length() > c2->get_length();
    return c1 < c2;
});
```

---

## 9. Efficiency — sqrt in Distance Comparisons

### Location
`JudgeSeparateDec_2`, approximately 30 occurrences of:
```cpp
double dis = sqrt(pow(a.x() - b.x(), 2) + pow(a.y() - b.y(), 2) + pow(a.z() - b.z(), 2));
if (dis < threshold) ...
```

### Problem
`sqrt` is called inside nested loops (over `boundary_points` × `hy/ly/hz/lz/hx/lx_points`
+ `independent_points`). `boundary_points` can contain hundreds of hull vertices.
`JudgeSeparateDec_2` is called for every cluster > 100 cm long, and called up to 5
times per cluster in the worst case (levels 0–2 + final + Separate_2 check).

### Fix
Replace all `sqrt(pow(a,2)+pow(b,2)+pow(c,2)) < threshold` with
`pow(a,2)+pow(b,2)+pow(c,2) < threshold*threshold`. Thresholds used: 25 cm and 15 cm.
Pre-compute the squared threshold constants.

Similarly in `clustering_separate` line ~500:
```cpp
double dis = sqrt(pow(hy_points.at(k).x() - boundary_points.at(j).x(), 2) + ...);
```
All can be squared.

---

## Execution Plan

### Step 1: Fix dead code and dangling pointer
- Remove `new_clusters` and `del_clusters` declarations
- Remove all `new_clusters.push_back(...)` calls
- Remove all `del_clusters.push_back(cluster)` calls

### Step 2: Fix sort tiebreaker
- Add `c1 < c2` as secondary sort key

### Step 3: Optimize sqrt in JudgeSeparateDec_2
- Replace `sqrt(pow()+pow()+pow()) < threshold` with squared comparison in all
  occurrences (approximately 30 sites, two threshold values: 25 cm and 15 cm)

### Step 4: Build and verify
