# Multi-APA/Face Support Analysis

Examination of how the six clustering functions extend the prototype's single-TPC
MicroBooNE geometry to support multiple APAs and faces.

---

## Summary Table

| Function | Multi-APA mechanism | Status |
|---|---|---|
| `ClusteringLiveDead` | Per-wpid `wpid_params` map; per-cluster wpid lookup (same pattern as Extend) | Fixed |
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

## `ClusteringLiveDead` Multi-APA Support — **Fixed**

### What Was Fixed

An earlier version of `clustering_live_dead.cxx` raised a `ValueError` for any
grouping with `wpids().size() > 1` and used a single `extract_geometry_params`
call that returned one `(drift_dir, angle_u, angle_v, angle_w)` tuple shared
across all clusters regardless of their APA/face.

This has been fixed. Current `clustering_live_dead.cxx` lines 71–78:

```cpp
// Build per-APA/face wire geometry maps (supports multiple APAs/faces)
const auto& all_wpids = live_grouping.wpids();
std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> wpid_params;
std::map<WirePlaneId, std::pair<geo_point_t, double>> wpid_U_dir;
std::map<WirePlaneId, std::pair<geo_point_t, double>> wpid_V_dir;
std::map<WirePlaneId, std::pair<geo_point_t, double>> wpid_W_dir;
std::set<int> apas;
compute_wireplane_params(all_wpids, m_dv, wpid_params, wpid_U_dir, wpid_V_dir, wpid_W_dir, apas);
```

The inner merge-decision loop (lines 217–221) now uses per-cluster lookups:

```cpp
auto wpid_1 = cluster_1->wpid(mcell1_center);
auto wpid_2 = cluster_2->wpid(mcell2_center);
const auto& [drift_dir_1, angle_u_1, angle_v_1, angle_w_1] = wpid_params.at(wpid_1);
const auto& [drift_dir_2, angle_u_2, angle_v_2, angle_w_2] = wpid_params.at(wpid_2);
```

Each `is_angle_consistent` call uses `angle_u_1/v_1/w_1` for cluster-1 direction
checks and `angle_u_2/v_2/w_2` for cluster-2 direction checks — identical to the
pattern in `clustering_extend`, `clustering_regular`, and `clustering_parallel_prolong`.

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
