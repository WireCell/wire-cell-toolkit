# Code Review: `shower_to_wall`, `broken_muon_id`, `gap_identification`

**Date:** 2026-04-08  
**Reviewer:** Claude Code  
**Scope:** Logic fidelity vs prototype, bug hunting, efficiency, determinism.

---

## Files Examined

| Role | File |
|---|---|
| Toolkit implementations | `clus/src/NeutrinoTaggerNuE.cxx` |
| Prototype implementations | `prototype_pid/src/NeutrinoID_nue_tagger.h` |
| Translation table / API map | `NeutrinoTaggerNuE.cxx:30-45` |
| FiducialUtils tolerance note | `clus/src/FiducialUtils.cxx:75` |

---

## Functions Reviewed

| Function | Toolkit Location | Prototype Location |
|---|---|---|
| `broken_muon_id` | `NeutrinoTaggerNuE.cxx:1368` | `NeutrinoID_nue_tagger.h:1010` |
| `shower_to_wall` | `NeutrinoTaggerNuE.cxx:2294` | `NeutrinoID_nue_tagger.h:1219` |
| `gap_identification` | `NeutrinoTaggerNuE.cxx:2689` | `NeutrinoID_nue_tagger.h:1438` |

---

## Bugs Found and Fixed

| Bug | Location | Fix |
|---|---|---|
| `stw_4_v_flag` uses accumulated `flag_bad4_save` instead of per-element flag | `shower_to_wall`, `NeutrinoTaggerNuE.cxx:2495` (pre-fix) | Introduced per-element `bool flag_bad4`; push `!flag_bad4`, accumulate separately into `flag_bad4_save` |

---

## `broken_muon_id`

### Purpose

Checks if the shower is a broken (gap-crossing) muon track by greedily walking forward from
the shower stem, chaining nearly-collinear segments (first shower-internal, then across cluster
gaps), and testing whether the accumulated track properties are more consistent with a muon
than an EM shower.

### Findings

#### ✅ `map_vtx_segs` / `map_seg_vtxs` → `fill_sets` + boost graph

The prototype uses shower-internal maps for Step A (connected continuation) and Step B (gap
crossing):

```cpp
Map_Proto_Vertex_Segments& map_vtx_segs = shower->get_map_vtx_segs();
Map_Proto_Segment_Vertices& map_seg_vtxs = shower->get_map_seg_vtxs();
```

The toolkit replaces both with:
```cpp
shower->fill_sets(shower_vtxs, shower_segs, false);
```
and traverses `boost::out_edges` filtered to `shower_segs` for Step A, and iterates
`shower_segs` directly for Step B. Semantically equivalent. ✅

#### ✅ Step A (connected continuation) — identical collinearity check

```cpp
// Both: back-to-back within 15° AND length > 6cm
if (180 - dir1.angle(dir2)/π*180 < 15 && length > 6*cm)
```
✅

#### ✅ Step B (gap crossing) — collinearity conditions correctly parenthesized

Prototype (one long expression with implicit `&&`>`||` grouping):
```cpp
if ( (std::min(angle1, angle2) < 10 && angle1+angle2 < 25 || angle3 < 15 && dis < 5*cm) && dis < 25*cm ||
     std::min(angle1, angle2) < 15 && angle3 < 30 && dis > 30*cm && length > 25*cm && dis < 60*cm ||
     (std::min(angle1, angle2) < 5 && angle1+angle2 < 15 || angle3 < 10 && dis < 5*cm) && dis < 30*cm )
```

Toolkit decomposes into named booleans:
```cpp
bool close_collinear  = ((min(a1,a2) < 10 && a1+a2 < 25) || (a3 < 15 && dis < 5*cm)) && dis < 25*cm;
bool far_collinear    = min(a1,a2) < 15 && a3 < 30 && dis > 30*cm && length > 25*cm && dis < 60*cm;
bool strict_collinear = (min(a1,a2) < 5 && a1+a2 < 15) || (a3 < 10 && dis < 5*cm);
bool passes = close_collinear || (strict_collinear && dis < 30*cm) || far_collinear;
```
Prototype clauses 1, 3, 2 map to `close_collinear`, `strict_collinear && dis < 30*cm`,
`far_collinear` respectively. Parenthesization matches prototype's implicit grouping. ✅

#### ✅ Vertex selection after gap crossing — farther endpoint chosen

```cpp
// Both: pick the endpoint of min_seg that is farther from curr_muon_vertex
if (dis4 > dis3) min_vtx = pair_vertices.second;
else             min_vtx = pair_vertices.first;
```
✅

#### ✅ `tmp_ids` (cluster IDs) → `tmp_clusters` (Cluster*) — equivalent

`tmp_ids.insert(seg->get_cluster_id())` → `tmp_clusters.insert(seg->cluster())`.  
Both count unique clusters and test cluster membership. ✅

#### ✅ `add_length` — correctly dropped

Prototype accumulates `add_length += min_dis` (gap distance) but uses it only in a `flag_print`
debug `std::cout` at the end. Toolkit drops it, documented with a comment. ✅

#### ✅ Connected length (`connected_length`) — equivalent traversal

Both sum lengths of shower-internal segments that belong to clusters in the muon track:
- Prototype: `if (tmp_ids.find(sg1->get_cluster_id()) != tmp_ids.end())`
- Toolkit: `if (tmp_clusters.count(sg1->cluster()))`
✅

#### ✅ `7022_42_2123` parallel-segment addition — equivalent

Both add shower-internal segments (not yet in muon_segments, same main cluster as sg) whose
minimum alignment angle to `dir_sg` is < 10°. ✅

#### ✅ Cut condition — identical

```cpp
// Both:
if (muon_segments.size() > 1 &&
    (Ep > Eshower*0.55 || acc_length > 0.65*total || connected > 0.95*total) &&
    n_clusters > 1 && acc_direct_length > 0.94*acc_length && Eshower < 350*MeV) {
    if (main_segs <= 3 && main_segs - num_muon_main < 2 &&
        (n_segs < n_muon_segs+6 || acc_length > connected*0.9 || acc_length > 0.8*total))
        flag_bad = true;
}
```
✅

#### ✅ All 11 fills match

| Field | Prototype | Toolkit |
|---|---|---|
| `brm_n_mu_segs` | `muon_segments.size()` | same ✅ |
| `brm_Ep` | `Ep/MeV` | same ✅ |
| `brm_energy` | `Eshower/MeV` | same ✅ |
| `brm_acc_length` | `acc_length/cm` | same ✅ |
| `brm_shower_total_length` | `shower->get_total_length()/cm` | same ✅ |
| `brm_connected_length` | `connected_length/cm` | same ✅ |
| `brm_n_size` | `tmp_ids.size()` | `tmp_clusters.size()` ✅ |
| `brm_acc_direct_length` | `acc_direct_length/cm` | same ✅ |
| `brm_n_shower_main_segs` | `shower->get_num_segments()` | same ✅ |
| `brm_n_mu_main` | `num_muon_main` | same ✅ |
| `brm_flag` | `!flag_bad` | same ✅ |

---

## `shower_to_wall`

### Purpose

Rejects events where the shower's backward direction points toward a detector wall (stw_1),
another non-electron shower (stw_2), another cluster's vertex (stw_3), or a forward shower
whose end is within 3 cm of the wall (stw_4).

### Findings

#### ✅ `stm_tol_vec` wall-walk tolerance — consistent with prior review

The prototype passes `stm_tol_vec = {-1.5*cm, ...}` to `fid->inside_fiducial_volume` for
the backward wall-walk and the stw_4 forward walk. The toolkit passes the same vector to
`fiducial_utils->inside_fiducial_volume`. As documented in `FiducialUtils.cxx:75`, the
toolkit's implementation ignores the tolerance vector entirely, delegating straight to
`m_sd.fiducial->contained(p)`. This is the same behavior established in the `angular_cut`
review. ✅

#### ✅ `n_other_shower` — dead variable in both

Both prototype and toolkit compute `n_other_shower` in the stw_1 vertex-shower loop (counting
companion showers with `E > 60*MeV`) but never reference it in any cut condition or fill.
Correctly retained as a no-op in the toolkit. ✅

#### ✅ `dir1.magnitude() == 0` guard moved before `norm()` — improvement

Prototype: `dir1 = dir1.Unit(); if (dir1.Mag()==0) continue;` — check after normalization,
so it never triggers (`.Unit()` always produces magnitude 1).  
Toolkit: `if (dir1.magnitude() == 0) continue; dir1 = dir1.norm();` — correctly guards
before normalization. ✅

#### ✅ stw_3 vertex loop — `map_vertex_segments` → `graph_nodes`

Prototype iterates `map_vertex_segments` (global vertex→segments map), checks
`vtx1->get_cluster_id() == vertex->get_cluster_id()` and `it->second.size()==1`.  
Toolkit iterates `graph_nodes(ctx.graph)`, checks `vtx1->cluster() == vertex->cluster()` and
`boost::out_degree(vd, ctx.graph) == 1`. Equivalent filters. ✅

#### ✅ stw_3 vertex point — `get_wcpt()` vs `wcpt().point`

Both use the raw wire-cell grid point (not the fit point) for the direction vector:
- Prototype: `TVector3 dir1(vtx1->get_wcpt().x - vertex_point.x, ...)`
- Toolkit: `Vector dir1 = vtx1->wcpt().point - vertex_point;`
✅

#### 🐛 FIXED: `stw_4_v_flag` used accumulated `flag_bad4_save` instead of per-element flag

**Root cause:** The prototype resets `flag_bad4 = false` at the top of each loop iteration and
uses `!flag_bad4` for the per-element push:

```cpp
// Prototype:
for (...) {
    flag_bad4 = false;   // reset per element
    ...
    if (dis1 < 3*cm && shower_energy < 500*MeV && flag_single_shower) flag_bad4 = true;
    if (flag_fill) tagger_info.stw_4_v_flag.push_back(!flag_bad4);
    if (flag_bad4) flag_bad4_save = true;
}
```

The toolkit instead wrote directly to the accumulated `flag_bad4_save` and pushed
`!flag_bad4_save` — meaning every entry after the first triggering shower would also carry
`flag=0` (bad) even if that entry did not satisfy the condition:

```cpp
// Pre-fix toolkit (WRONG):
if (dis1 < 3*cm && ...) flag_bad4_save = true;   // accumulates
ti.stw_4_v_flag.push_back(!flag_bad4_save);       // wrong: not per-element
```

**Fix applied:** Introduced a per-element `bool flag_bad4 = false` inside the loop:

```cpp
// Post-fix toolkit (CORRECT):
bool flag_bad4 = false;
if (dis1 < 3*cm && shower_energy < 500*MeV && flag_single_shower) flag_bad4 = true;
if (flag_bad4) flag_bad4_save = true;
ti.stw_4_v_flag.push_back(!flag_bad4);   // per-element, as in prototype
```

The overall rejection decision (`flag_bad4_save` feeding into the final `flag_bad`) was
correct before and after the fix. Only the per-element `stw_4_v_flag` vector was wrong.

#### ✅ All stw_1 fills (7 fields), stw_2 per-element fills (6 fields), stw_3 per-element fills (5 fields), stw_4 per-element fills (4 fields), and `stw_flag` match

---

## `gap_identification`

### Purpose

Checks for signal gaps in the shower stem near the interaction vertex by generating 3
sub-points for each consecutive pair of fit points within 2.4 cm of the vertex and querying
the charge point cloud for wire signal in each plane (U, V, W). A sub-point with fewer than
3 contributing planes (signal + dead channel + prolongation exemption) increments `n_bad`.

### Findings

#### ✅ `check_direction` inlining — correct algorithm

Prototype: `std::vector<bool> flag = main_cluster->check_direction(dir)` where `dir` is the
vector from `vertex_point` to the fit point closest to 3 cm along the stem.

Toolkit inlines the same logic:
1. Find fit point closest to 3 cm from vertex → `dir_cd = closest_p − vertex_point`
2. Get wire angles for U, V, W planes via `grouping->wire_angles(ctx.apa, ctx.face)`
3. For each wire plane: project `dir_cd` into the YZ plane, compute the angle to the wire
   direction, then form the effective 2D vector vs drift axis — flag prolongation if angle <
   12.5° (≈ `cut1`)
4. `flag_parallel`: check if `|dir_cd.angle(drift_abs) − π/2| < 10°` (≈ `cut2`)

The `check_prolong` lambda reproduces the same computation. ✅

#### ✅ `closest_p` initialization — both branches identical

In both the `flag_start=true` and `flag_start=false` branches the prototype initializes
`closest_p = pts.back()` then searches all `pts`. The toolkit initializes `closest_p =
vertex_point` and searches all `sg_fits`. The search overwrites `closest_p` with the
nearest-to-3cm element either way; the initialization only matters if `sg_fits` is empty
(in which case the toolkit leaves `closest_p = vertex_point` and `dir_cd = zero`, making all
prolong flags false — a safe default). ✅

#### ✅ Sub-point scan — equivalent loop structure

The toolkit's `query_sub_points(i_start, i_end, i_step)` lambda reproduces both the
forward (`flag_start=true`: i=0…n-2, inner j=0,1,2) and backward
(`flag_start=false`: i=n-1…1, inner j=0,1,2) traversals.

Sub-point formula: `test_p = sg_fits[i].point + (sg_fits[i_next].point - sg_fits[i].point) * (j/3.0)` —
matches prototype's `pts[i] + j/3.*(pts[next]−pts[i])`. ✅

Break condition: both break when the next fit point exceeds 2.4 cm from `vertex_point`. ✅

#### ✅ `get_closest_dead_chs` API change — documented

Prototype: `ct_point_cloud->get_closest_dead_chs(test_p, plane, 0)` — third arg is a mode flag.  
Toolkit: `grouping->get_closest_dead_chs(test_p, 1, q_apa, q_face, plane)` — second arg is
radius (1 wire). Documented in the translation notes at the top of the function. ✅

#### ✅ `n_bad` reset logic — identical

```cpp
// Both:
if (num_connect + num_bad_ch + num_spec == 3) {
    if (n_bad == n_points && n_bad <= 5) n_bad = 0;
} else {
    ++n_bad;
}
++n_points;
```
✅

#### ✅ Cut conditions — identical

| Condition | Prototype | Toolkit |
|---|---|---|
| `E > 900*MeV`, `!single` `!parallel`, `E > 1200*MeV`: `n_bad > 2/3 * n_pts` | ✅ | ✅ |
| `E > 900*MeV`, `!single` `!parallel`, else: `n_bad > 1/3 * n_pts` | ✅ | ✅ |
| `E > 900*MeV`, `parallel` `!single`: `n_bad > 0.5 * n_pts` (`1/2.` in prototype) | ✅ | ✅ |
| `150 < E ≤ 900*MeV`, `!single`, `parallel`: `n_bad > 4` | ✅ | ✅ |
| `150 < E ≤ 900*MeV`, `!single`, `!parallel`: `n_bad > 1` | ✅ | ✅ |
| `150 < E ≤ 900*MeV`, `single`: `n_bad > 2` | ✅ | ✅ |
| `E ≤ 150*MeV`, `!single`, `parallel`: `n_bad > 3` | ✅ | ✅ |
| `E ≤ 150*MeV`, `!single`, `!parallel`: `n_bad > 1` | ✅ | ✅ |
| `E ≤ 150*MeV`, `single`: `n_bad > 2` | ✅ | ✅ |
| 7021: `n_bad >= 6 && E < 1000*MeV` | ✅ | ✅ |
| `E ≤ 900*MeV && n_bad > 1` | ✅ | ✅ |

#### ✅ All 11 fills match

| Field | Prototype | Toolkit |
|---|---|---|
| `gap_flag` | `!flag_gap` | same ✅ |
| `gap_flag_prolong_u` | `flag_prolong_u` | same ✅ |
| `gap_flag_prolong_v` | `flag_prolong_v` | same ✅ |
| `gap_flag_prolong_w` | `flag_prolong_w` | same ✅ |
| `gap_flag_parallel` | `flag_parallel` | same ✅ |
| `gap_n_points` | `n_points` | same ✅ |
| `gap_n_bad` | `n_bad` | same ✅ |
| `gap_energy` | `E_shower/MeV` | same ✅ |
| `gap_num_valid_tracks` | `num_valid_tracks` | same ✅ |
| `gap_flag_single_shower` | `flag_single_shower` | same ✅ |
| `gap_filled` | `1` | same ✅ |

---

## Summary Table

| Check | Result | Action |
|---|---|---|
| `broken_muon_id` — `fill_sets` replaces `map_seg_vtxs`/`map_vtx_segs` | ✅ Equivalent | — |
| `broken_muon_id` — Step A back-to-back collinearity | ✅ Identical | — |
| `broken_muon_id` — Step B gap-crossing conditions parenthesized | ✅ Identical grouping | — |
| `broken_muon_id` — farther-endpoint vertex selection | ✅ Identical | — |
| `broken_muon_id` — `tmp_ids` → `tmp_clusters` | ✅ Equivalent | — |
| `broken_muon_id` — `add_length` dropped | ✅ Debug-print-only | — |
| `broken_muon_id` — `connected_length` and `num_muon_main` | ✅ Identical | — |
| `broken_muon_id` — cut condition and 11 fills | ✅ Identical | — |
| `shower_to_wall` — `stm_tol_vec` tolerance behavior | ✅ Consistent with FiducialUtils.cxx:75 | — |
| `shower_to_wall` — `n_other_shower` dead variable | ✅ Kept as no-op (matches prototype) | — |
| `shower_to_wall` — `dir1.magnitude() == 0` before `norm()` | ✅ Improvement over prototype | — |
| `shower_to_wall` — stw_3 `graph_nodes` vs `map_vertex_segments` | ✅ Equivalent | — |
| `shower_to_wall` — stw_3 `vtx1->wcpt().point` vs `get_wcpt().x` | ✅ Same point | — |
| **`shower_to_wall` — stw_4 `stw_4_v_flag` used accumulated flag** | **🐛 Bug** | **Fixed: per-element `flag_bad4`** |
| `shower_to_wall` — stw_1/stw_2/stw_3/stw_4 fills | ✅ Match (after fix) | — |
| `gap_identification` — `check_direction` inlining | ✅ Correct algorithm | — |
| `gap_identification` — `closest_p` initialization | ✅ Same point found | — |
| `gap_identification` — sub-point scan loop (forward and backward) | ✅ Equivalent | — |
| `gap_identification` — `get_closest_dead_chs` API change | ✅ Documented | — |
| `gap_identification` — `n_bad` reset logic | ✅ Identical | — |
| `gap_identification` — all 11 cut conditions | ✅ Identical | — |
| `gap_identification` — all 11 fills | ✅ Match | — |

---

## Changes Made

**`clus/src/NeutrinoTaggerNuE.cxx` — `shower_to_wall` stw_4 loop**

Introduced a per-element `bool flag_bad4 = false` inside the stw_4 loop. The push to
`stw_4_v_flag` now uses `!flag_bad4` (per-element, matching prototype semantics), and
accumulation into `flag_bad4_save` happens explicitly after the per-element check.
The overall rejection signal (`flag_bad4_save`) was already correct; only the per-element
vector fill was wrong.
