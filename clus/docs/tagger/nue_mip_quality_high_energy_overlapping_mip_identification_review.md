# Code Review: `mip_quality`, `high_energy_overlapping`, `mip_identification`

**Date:** 2026-04-08  
**Reviewer:** Claude Code  
**Scope:** Logic fidelity vs prototype, bug hunting, efficiency, determinism.

---

## Files Examined

| Role | File |
|---|---|
| Toolkit implementations | `clus/src/NeutrinoTaggerNuE.cxx` |
| Prototype implementations | `prototype_pid/src/NeutrinoID_nue_tagger.h` |

---

## Functions Reviewed

| Function | Toolkit Location | Prototype Location |
|---|---|---|
| `mip_quality` | `NeutrinoTaggerNuE.cxx:1591` | `NeutrinoID_nue_tagger.h:1669` |
| `high_energy_overlapping` | `NeutrinoTaggerNuE.cxx:1781` | `NeutrinoID_nue_tagger.h:2284` |
| `mip_identification` | `NeutrinoTaggerNuE.cxx:3806` | `NeutrinoID_nue_tagger.h:1823` |

---

## Bugs Found and Fixed

### `high_energy_overlapping`: `flag_all_showers` incorrectly cleared on zero-direction segment

**Location:** `NeutrinoTaggerNuE.cxx:1829` (before fix)

**Prototype behaviour** (`NeutrinoID_nue_tagger.h:2329`):

```cpp
if ((*it)->get_particle_type()==11 || ...) {
    TVector3 dir2 = (*it)->cal_dir_3vector(vtx_point, 5*units::cm);
    if (dir2.Mag() == 0) continue;   // just skip — flag_all_showers unchanged
    ...
} else {
    flag_all_showers = false;   // only non-pdg11/muon segments clear this
}
```

`flag_all_showers` is intended to be `true` iff every vertex segment (other than `sg`) is a
pdg-11 electron shower or a short weak-direction muon. When a pdg-11/muon segment happens to have
zero computed direction the prototype simply skips it — the segment's particle type still qualifies
it, so `flag_all_showers` stays `true`.

**Toolkit bug** (before fix):

```cpp
if (is_pdg11 || is_weak_muon) {
    ...
    if (dir2.magnitude() == 0) { flag_all_showers = false; continue; }  // WRONG
```

Setting `flag_all_showers = false` here was incorrect: it treated a pdg-11 segment with a
degenerate direction vector as if it were a hadronic track, potentially suppressing the two
`flag_overlap1` conditions that check `flag_all_showers`.

**Fix applied:**

```cpp
if (dir2.magnitude() == 0) { continue; }
```

The `flag_all_showers = false` was removed. The segment is still skipped for the
angle/length update (as in the prototype), but no longer changes the shower-composition flag.

---

## `mip_quality`

### Purpose

Checks two rejection modes for the shower candidate:
- `flag_overlap`: first ≤3 stem fit points lie in 2D wire proximity of another shower-internal
  segment — likely a crossing track misreconstructed as a shower branch.
- `flag_split`: exactly two pdg-11 showers at the vertex with no tracks — the event looks like
  a single shower split across a reconstruction gap.

### Findings

#### ✅ `n_protons` correctly dropped

Prototype computes `n_protons` inside the vertex-loop track count:

```cpp
if (sg1->get_particle_type()==2212) n_protons ++;
```

This variable is never used in any condition or fill. Toolkit omits it. ✅

#### ✅ Loop variable `shower` shadow resolved

Inside the prototype's `map_vertex_to_shower` loop, the loop variable is named `shower`,
shadowing the outer `shower` function parameter:

```cpp
WCPPID::WCShower *shower = *it1;  // shadows outer parameter
WCPPID::ProtoSegment *sg1 = shower->get_start_segment();  // loop's shower
```

The toolkit renames the loop variable to `shower1`. ✅

#### ✅ `map_vertex_segments[main_vertex].find(sg1)` → boost::out_edges on `ctx.main_vertex`

Both check whether the start segment `sg1` of the candidate shower connects to the main
interaction vertex. Equivalent. ✅

#### ✅ Overlap check: `map_seg_vtxs` → `shower_segs` from `fill_sets`

Prototype iterates `shower->get_map_seg_vtxs()` to build the set of other shower segments.
Toolkit calls `shower->fill_sets(shower_vtxs_tmp, shower_segs, false)`.
Both produce the same set of shower-internal segments. ✅

#### ✅ `nconnected` via `boost::out_degree` — equivalent

Prototype: `map_vertex_segments[other_vertex].size()`  
Toolkit: `boost::out_degree(other_vtx->get_descriptor(), ctx.graph)`  
Same count. ✅

#### ✅ Vertex-end / other-end exception logic — equivalent

Prototype exception conditions (implicit test-loop index `i`):
- `i == 0 && u==0 && v==0 && w==0`: first test point at vertex end — shared vertex, skip.
- `i+1 == n && u==0 && v==0 && w==0 && nconnected==2`: only triggers when n (segment size)
  ≤ 3, i.e., the far end of a short segment is at a shared two-way vertex — skip.

Toolkit: `k == 0` (first test point, vertex end) and `is_other_end = (orig_idx == n_fits-1)`
when vertex_at_front or `(orig_idx == 0)` when vertex_at_back. Both conditions fire at
the same points as the prototype. ✅

#### ✅ `flag_bad` precedence — explicitly parenthesized

Prototype:
```cpp
flag_bad = (Eshower < 800*units::MeV) && flag_overlap || (Eshower < 500*units::MeV) && flag_split;
```
With C++ precedence (`&&` tighter than `||`) this is two &&-clauses ||'d together.

Toolkit:
```cpp
bool flag_bad = ((Eshower < 800*units::MeV) && flag_overlap) ||
                ((Eshower < 500*units::MeV) && flag_split);
```
Explicit parentheses, same semantics. ✅

#### ✅ `shortest_length > 10*units::m` → 1000 cm fill

Prototype: `if (shortest_length > 10*units::m) tagger_info.mip_quality_shortest_length = 1000;`  
Toolkit: ternary with `1000.0f`. Same. ✅

#### ✅ All 12 fills match

| Field | Prototype | Toolkit |
|---|---|---|
| `mip_quality_flag` | `!flag_bad` | same ✅ |
| `mip_quality_energy` | `Eshower/MeV` | same ✅ |
| `mip_quality_overlap` | `flag_overlap` | same ✅ |
| `mip_quality_n_showers` | `n_showers` | same ✅ |
| `mip_quality_n_tracks` | `n_tracks` | same ✅ |
| `mip_quality_flag_inside_pi0` | `flag_inside_pi0` | same ✅ |
| `mip_quality_n_pi0_showers` | `tmp_pi0_showers.size()` | same ✅ |
| `mip_quality_shortest_length` | `1000` or `length/cm` | same ✅ |
| `mip_quality_acc_length` | `shortest_acc_length/cm` | same ✅ |
| `mip_quality_shortest_angle` | `shortest_angle` | same ✅ |
| `mip_quality_flag_proton` | `flag_proton` | same ✅ |
| `mip_quality_filled` | `1` | same ✅ |

---

## `high_energy_overlapping`

### Purpose

Two complementary overlap checks for high-energy showers at a type-1 vertex:
- `flag_overlap1`: no valid tracks at vertex, with a small minimum angle to another
  electron/weak-muon segment — likely an overlapping cluster.
- `flag_overlap2`: consecutive stem fit points are within 0.6 cm of the most collinear
  vertex segment, and that segment has elevated dQ/dx in the overlap region.

### Findings

#### ❌ → ✅ Bug fixed: `flag_all_showers` on zero-direction segment (see above)

#### ✅ `flag_start` / `vtx_point` — equivalent to prototype

Prototype: compares `vtx->get_wcpt().index == sg->get_wcpt_vec().front().index`.  
Toolkit: compares `ray_length(Ray{vtx_fit_pt(vtx), sg_fits.front().point})` ≤ back distance.  
Both determine which end of `sg` is near the vertex and set `vtx_point` accordingly. ✅

#### ✅ `conn_type == 1` guard — both flag_overlap1 and flag_overlap2 inside

Prototype uses two separate `if (pair_result.second == 1)` blocks. Toolkit combines them into a
single `if (conn_type == 1)` block. Equivalent control flow. ✅

#### ✅ `n_valid_tracks` / `num_showers` loop — identical

All conditions, thresholds, and `n_valid_tracks` increment rules match prototype exactly. ✅

#### ✅ `medium_dQ_dx` end detection — proximity replaces wcpt index comparison

Prototype determines which end of `min_sg` is near the vertex by comparing wcpt indices of `sg`
and `min_sg`:

```cpp
// flag_start == true:
if (sg->get_wcpt_vec().front().index == min_sg->get_wcpt_vec().front().index)
    medium_dQ_dx = min_sg->get_medium_dQ_dx(0, ncount) / ...;
else
    medium_dQ_dx = min_sg->get_medium_dQ_dx(n-1-ncount, n-1) / ...;
```

Toolkit uses direct proximity:

```cpp
bool min_front_near = (ray_length(Ray{vtx_point, min_fits.front().point}) <=
                       ray_length(Ray{vtx_point, min_fits.back().point}));
if (min_front_near) medium_dQ_dx = segment_median_dQ_dx(min_sg, 0, ncount) / ...;
else                medium_dQ_dx = segment_median_dQ_dx(min_sg, n_min-1-ncount, n_min-1) / ...;
```

Functionally equivalent: if `min_sg` connects to the same vertex as `sg`, its near-end fit
point is close to `vtx_point`. ✅

#### ✅ Dead `min_dQ_dx` loop variable dropped

Prototype declares `double min_dQ_dx = 0;` inside the forward ncount loop but never reads it.
Toolkit omits it. ✅

#### ✅ All 11 fills match

| Field | Prototype | Toolkit |
|---|---|---|
| `hol_1_n_valid_tracks` | `n_valid_tracks` | same ✅ |
| `hol_1_min_angle` | `min_angle` | same ✅ |
| `hol_1_energy` | `Eshower/MeV` | same ✅ |
| `hol_1_flag_all_shower` | `flag_all_showers` | same ✅ |
| `hol_1_min_length` | `min_length/cm` | same ✅ |
| `hol_1_flag` | `!flag_overlap1` | same ✅ |
| `hol_2_min_angle` | `min_angle` | same ✅ |
| `hol_2_medium_dQ_dx` | `medium_dQ_dx` | same ✅ |
| `hol_2_ncount` | `ncount` | same ✅ |
| `hol_2_energy` | `Eshower/MeV` | same ✅ |
| `hol_2_flag` | `!flag_overlap2` | same ✅ |
| `hol_flag` | `!flag_overlap` | same ✅ |

---

## `mip_identification`

### Purpose

Classifies the shower stem as EM-like (`1`), ambiguous (`0`), or MIP/track-like (`-1`) using
the per-fit-point dQ/dx profile, shower energy, single-shower topology, and various geometric
and quality checks. Returns the integer classification.

### Findings

#### ✅ Dead shower-loop variables correctly dropped

Prototype computes in the indirect-shower loop:

```cpp
double E_direct_max_energy = 0, E_direct_total_energy = 0;
int n_direct_showers = 0;
int n_indirect_showers = 0;
double E_indirect_total_energy = 0;
```

None of these five variables appear in any fill. Only `E_indirect_max_energy` is used (two
conditions + one fill). Toolkit computes only `E_indirect_max_energy`. ✅

#### ✅ Loop variable `sg` shadow resolved

Prototype indirect-shower loop: `WCPPID::ProtoSegment *sg = shower1->get_start_segment()`
shadows the outer function parameter `sg`. C++ resolves correctly to the local, but the
toolkit renames to `sg1` to eliminate the shadowing. ✅

#### ✅ `vec_dQ_dx` size guard — defensive improvement

Prototype directly accesses `.front()` and `.at(0)/.at(1)` before any size check. Toolkit
adds:

```cpp
while (vec_dQ_dx.size() < 2) vec_dQ_dx.push_back(0.0);
```

This prevents UB when `get_stem_dQ_dx` returns an empty or size-1 vector. The additional `0.0`
values are below every threshold so logic is unchanged. ✅

#### ✅ `vec_dQ_dx` fill-padding — equivalent for indices [0..19]

Prototype appends exactly 20 × `3.0` regardless of current size
(so a size-20 input becomes size-40). Toolkit pads to minimum size-20.
Only indices [0..19] are written to `ti.mip_vec_dQ_dx_*`, so both produce
identical fill values for all `get_stem_dQ_dx` outputs (which are capped at 20). ✅

#### ✅ `median_dQ_dx` computed before padding

Both prototype and toolkit compute the median on the unpadded vector.
Padding with `3.0` after the median preserves the quality-check result. ✅

#### ✅ `n_first_mip` / `n_first_non_mip_1` / `n_first_non_mip_2` scan loops — identical

All three loops (L3861-3887 toolkit; L1880-1908 prototype) produce the same integer results.
The toolkit uses `vec_dQ_dx[i]` and `vec_threshold[i]` direct indexing; prototype uses
`.at(i)`. Equivalent. ✅

#### ✅ `early_mip` extracted to named bool — same semantics

Prototype: `if (n_first_non_mip_2 - n_first_mip >= 2 && (n_first_mip <= 2 || ...)) mip_id = 1; else mip_id = -1;`

Toolkit: `int run = n_first_non_mip_2 - n_first_mip; bool early_mip = ...; if (run >= 2 && early_mip) mip_id = 1; else mip_id = -1;`

Same condition, clearer naming. ✅

#### ✅ `flag_strong_check` override — identical

All thresholds and sub-conditions in the strong-check block match, including the
`run == 3 && max(dQ_dx[0], dQ_dx[1]) > 3.3` exception. ✅

#### ✅ `vtx_n_segs` via `boost::out_degree` — equivalent

Prototype: `map_vertex_segments[vertex].size()`  
Toolkit: `boost::out_degree(vertex->get_descriptor(), ctx.graph)`  
Same count. Used in single/multi-vertex branch conditions and `n_lowest <= 2` cut. ✅

#### ✅ `n_other_vertex` and other-end check — equivalent

Prototype: `map_vertex_segments[other_vertex].size() > 2`  
Toolkit: `other_deg > 2` via `boost::out_degree`  

Prototype: `sg->get_point_vec().size() <= n_first_mip + 1`  
Toolkit: `n_fits <= (size_t)(n_first_mip + 1)` where `n_fits = sg->fits().size()`  
Equivalent. ✅

#### ✅ `min_dis` computation — equivalent

Prototype iterates `map_seg_vtxs` (off-cluster segments) × `map_vtx_segs` (shower vertices
in main cluster). Toolkit uses `fill_sets` and filters `shower_segs` / `shower_vtxs` for
off-cluster vs main-cluster. Same pairwise distance computation, same `min_dis`. ✅

#### ✅ `shower->get_total_length(vertex->get_cluster_id())` → `shower->get_total_length(vertex->cluster())`

`vertex->cluster()` is the cluster pointer corresponding to `vertex->get_cluster_id()`. ✅

#### ✅ `iso_angle` formula — identical

Prototype: `fabs(pi/2 - dir_shower.Angle(dir_drift))/pi*180`  
Toolkit: `std::fabs(M_PI/2.0 - dir_shower.angle(drift_dir)) / M_PI * 180.0`  
Same. ✅

#### ✅ All 32 fills match (12 scalar + 20 `mip_vec_dQ_dx_*`)

All scalar fills (`mip_flag`, `mip_energy`, `mip_n_*`, `mip_length_*`, `mip_*_dQ_dx`, etc.)
and all 20 per-stem-point dQ/dx fills (`mip_vec_dQ_dx_0` … `mip_vec_dQ_dx_19`) match the
prototype in name, value expression, and order. ✅

---

## Summary Table

| Check | Result | Action |
|---|---|---|
| `mip_quality` — `n_protons` dead variable dropped | ✅ Correctly dropped | — |
| `mip_quality` — loop `shower` shadow resolved | ✅ Renamed to `shower1` | — |
| `mip_quality` — `main_vertex` connectivity via boost | ✅ Equivalent | — |
| `mip_quality` — overlap check via `fill_sets` shower_segs | ✅ Equivalent to `map_seg_vtxs` | — |
| `mip_quality` — vertex-end exception (k==0, is_other_end) | ✅ Equivalent | — |
| `mip_quality` — `flag_bad` precedence | ✅ Explicitly parenthesized | — |
| `mip_quality` — all 12 fills | ✅ Match | — |
| `high_energy_overlapping` — `flag_all_showers` on zero-direction pdg11 | ❌ Bug: set false instead of continue | **Fixed** |
| `high_energy_overlapping` — `flag_start` / `vtx_point` | ✅ Proximity equivalent to wcpt index | — |
| `high_energy_overlapping` — `conn_type == 1` combined block | ✅ Equivalent | — |
| `high_energy_overlapping` — `medium_dQ_dx` end detection | ✅ Proximity replaces wcpt index | — |
| `high_energy_overlapping` — dead `min_dQ_dx` dropped | ✅ Correctly dropped | — |
| `high_energy_overlapping` — all 12 fills | ✅ Match | — |
| `mip_identification` — dead `E_direct_*`, `n_direct/indirect` dropped | ✅ Correctly dropped | — |
| `mip_identification` — loop `sg` shadow resolved | ✅ Renamed to `sg1` | — |
| `mip_identification` — `vec_dQ_dx` size guard (< 2) | ✅ Defensive improvement | — |
| `mip_identification` — fill-padding (to 20) | ✅ Equivalent for indices [0..19] | — |
| `mip_identification` — median before padding | ✅ Correct order | — |
| `mip_identification` — `early_mip` bool extraction | ✅ Same semantics | — |
| `mip_identification` — `flag_strong_check` override | ✅ Identical | — |
| `mip_identification` — `vtx_n_segs` via `boost::out_degree` | ✅ Equivalent | — |
| `mip_identification` — `n_other_vertex` connectivity check | ✅ Equivalent | — |
| `mip_identification` — `min_dis` via `fill_sets` | ✅ Equivalent to `map_seg_vtxs` | — |
| `mip_identification` — all 32 fills | ✅ Match | — |

---

## Changes Made

### `high_energy_overlapping` — `flag_all_showers` fix

**File:** `clus/src/NeutrinoTaggerNuE.cxx:1829`

```cpp
// Before (WRONG):
if (dir2.magnitude() == 0) { flag_all_showers = false; continue; }

// After (CORRECT):
if (dir2.magnitude() == 0) { continue; }
```

A pdg-11 / weak-muon segment with a degenerate direction vector should not clear
`flag_all_showers`. The prototype's `continue` skips the angle update without modifying the
flag. The two `flag_overlap1` conditions that gate on `flag_all_showers` could have been
incorrectly suppressed in the toolkit.
