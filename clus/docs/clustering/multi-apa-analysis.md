# Multi-APA/Face Support Analysis

Examination of how the six clustering functions extend the prototype's single-TPC
MicroBooNE geometry to support multiple APAs and faces.

---

## Summary Table

| Function | Multi-APA mechanism | Status |
|---|---|---|
| `ClusteringLiveDead` | Single-APA only; hard exception if grouping has >1 wpid | **BUG** |
| `ClusteringExtend` | Per-wpid `wpid_U/V/W_dir` maps; per-cluster wpid lookup | Correct |
| `ClusteringRegular` | Per-wpid `wpid_U/V/W_dir` maps; `get_wireplaneid` for cross-APA pairs | Correct |
| `ClusteringParallelProlong` | Per-wpid maps via `compute_wireplane_params`; `wpid_ps` dominant-APA | Correct |
| `ClusteringClose` | No wire-direction geometry (angle checks only by distance/length) | N/A |
| `ClusteringExtendLoop` | Inherits from `ClusteringExtend` | Correct |

---

## The Multi-APA Extension Pattern

The prototype hardcoded MicroBooNE wire geometry:
```cpp
// prototype (ToyClustering_extend.h)
TVector3 U_dir(0, cos(60°), sin(60°));
TVector3 V_dir(0, cos(60°), -sin(60°));
TVector3 W_dir(0, 1, 0);
```

The toolkit generalizes this via per-APA/face maps.  For each wpid in the
grouping, `compute_wireplane_params` (or inline equivalent) builds:
- `wpid_U_dir[wpid]` = U-wire direction + angle for that APA/face
- `wpid_V_dir[wpid]` = V-wire direction + angle
- `wpid_W_dir[wpid]` = W-wire direction + angle

For a cluster pair (c1, c2), the lookup key is the wpid of the relevant
point on that cluster, obtained via `cluster.wpid(point)`.

Cross-APA pairs are resolved by `get_wireplaneid(p1, wpid_p1, p2, wpid_p2, dv)`,
which returns the wpid whose APA bounding box has the longer ray intersection.

---

## Bug: `ClusteringLiveDead` — Hard Exception for Multi-APA

### Location
`clus/src/clustering_live_dead.cxx`, lines 71–81.

### Code
```cpp
// Check that groupings has less than one wpid
if (live_grouping.wpids().size() > 1 || dead_grouping.wpids().size() > 1) {
    for (const auto& wpid : live_grouping.wpids()) {
        std::cout << "Live grouping wpid: " << wpid.name() << std::endl;
    }
    for (const auto& wpid : dead_grouping.wpids()) {
        std::cout << "Dead grouping wpid: " << wpid.name() << std::endl;
    }
    raise<ValueError>("Live %d > 1, Dead %d > 1",
                      live_grouping.wpids().size(), dead_grouping.wpids().size());
}
auto [drift_dir, angle_u, angle_v, angle_w] = extract_geometry_params(live_grouping, m_dv);
```

### Problem

1. **Hard crash for any multi-APA grouping.** `grouping.wpids()` returns all
   distinct blob WirePlaneIds across all clusters. Any real multi-APA detector
   produces groupings with `wpids().size() > 1`, causing a `ValueError` at
   runtime before any geometry is processed.

2. **Single-APA geometry for all clusters.** `extract_geometry_params` reads
   only the first wpid and returns a single `(drift_dir, angle_u, angle_v,
   angle_w)` tuple. The `angle_u/v/w` values are then used in
   `is_angle_consistent` calls for both cluster 1 and cluster 2 regardless of
   which APA they belong to.

### Contrast with Other Functions

`clustering_extend`, `clustering_regular`, and `clustering_parallel_prolong` all:
- Do NOT raise an exception for multi-APA groupings.
- Build per-wpid maps from all grouping wpids.
- Look up geometry per-cluster using `cluster.wpid(point)`.
- Apply cluster-1's angles to cluster-1's direction check and cluster-2's
  angles to cluster-2's direction check.

### Fix Required

Replace lines 71–81 with per-wpid map building via `compute_wireplane_params`:

```cpp
// Build per-APA/face wire geometry maps (supports multiple APAs/faces)
const auto& all_wpids = live_grouping.wpids();
std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> wpid_angle_params;
std::map<WirePlaneId, std::pair<geo_point_t, double>> wpid_U_dir_ld, wpid_V_dir_ld, wpid_W_dir_ld;
std::set<int> apas_ld;
compute_wireplane_params(all_wpids, m_dv,
    wpid_angle_params, wpid_U_dir_ld, wpid_V_dir_ld, wpid_W_dir_ld, apas_ld);
```

Then, in the inner merge-decision loop (lines ~221–265), replace the single
`angle_u, angle_v, angle_w` with per-cluster lookups:

```cpp
// Look up per-cluster APA geometry
auto wpid_1 = cluster_1->wpid(mcell1_center);
auto wpid_2 = cluster_2->wpid(mcell2_center);
const auto& [drift_dir_1, angle_u_1, angle_v_1, angle_w_1] = wpid_angle_params.at(wpid_1);
const auto& [drift_dir_2, angle_u_2, angle_v_2, angle_w_2] = wpid_angle_params.at(wpid_2);

// drift angles (sign-insensitive for fabs(angle-π/2) checks, but use correct APA)
angle1 = dir1.angle(drift_dir_1);
angle2 = dir2.angle(drift_dir_1);
angle3 = dir3.angle(drift_dir_2);
```

And replace each `is_angle_consistent(dir1, ..., angle_u, angle_v, angle_w, ...)` call
with `angle_u_1/v_1/w_1` (cluster 1's APA) and each `is_angle_consistent(dir3, ..., angle_u, angle_v, angle_w, ...)` with `angle_u_2/v_2/w_2` (cluster 2's APA).

---

## Non-Bugs (Confirmed Correct)

### `drift_dir_abs(1, 0, 0)` hardcoding in `Clustering_4th_reg`, `Clustering_4th_dead`

Both functions use `geo_point_t drift_dir_abs(1, 0, 0)` exclusively in checks of the
form `fabs(angle - π/2) < threshold`. Because `fabs(dir.angle(-X) - π/2) = fabs(dir.angle(+X) - π/2)`,
these checks are sign-insensitive with respect to drift direction. The hardcoding
does not cause incorrect results for faces with -X drift. **Not a bug.**

### `get_wireplaneid` cross-APA choice in `clustering_regular` / `clustering_extend`

When `wpid_p1 != wpid_p2`, `get_wireplaneid` returns the APA whose bounding box
has the longer intersection with the p1→p2 ray. Both clusters' wpids are always
present in the geometry maps (built from all grouping wpids), so the fallback
cannot crash and returns the geometrically dominant APA. Reasonable approximation
for cross-APA pairs. **Not a bug.**

### Single `wpid_ps` in `clustering_parallel_prolong` for both direction checks

For the parallel/prolong case, both clusters must be nearly perpendicular to drift
and nearly co-planar. In practice, cross-APA parallel candidates are geometrically
unlikely (they would be at the same Y-Z location but opposite sides of a cathode).
Using one dominant `wpid_ps` for both direction angle checks is a safe simplification.
**Not a bug.**

### `clustering_close` has no per-APA geometry

`Clustering_3rd_round` makes merge decisions purely by distance thresholds, cluster
lengths, Hough directions, and local point counts — none of which depend on wire
plane orientation. Multi-APA extension is therefore automatic and requires no
per-wpid maps. **Not applicable / already correct.**

---

## `is_connected` Cross-APA Behavior in `ClusteringLiveDead`

`live->is_connected(*dead, dead_live_overlap_offset_)` checks whether a live cluster
blob overlaps in wire+time with any blob in the dead cluster. This uses integer
wire-index ranges and does not depend on APA geometry — it is already correct for
multi-APA groupings. The only APA-specific logic is in the subsequent angle/direction
checks, which is exactly the part affected by the bug above.
