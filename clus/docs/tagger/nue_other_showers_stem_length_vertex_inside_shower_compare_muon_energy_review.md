# Code Review: `other_showers`, `stem_length`, `vertex_inside_shower`, `compare_muon_energy`

**Date:** 2026-04-08  
**Reviewer:** Claude Code  
**Scope:** Logic fidelity vs prototype, bug hunting, efficiency, determinism.

---

## Files Examined

| Role | File |
|---|---|
| Toolkit implementations | `clus/src/NeutrinoTaggerNuE.cxx` |
| Prototype implementations | `prototype_pid/src/NeutrinoID_nue_functions.h` |
| Translation table / API map | `NeutrinoTaggerNuE.cxx:30-45` |

---

## Functions Reviewed

| Function | Toolkit Location | Prototype Location |
|---|---|---|
| `stem_length` | `NeutrinoTaggerNuE.cxx:189` | `NeutrinoID_nue_functions.h:379` |
| `compare_muon_energy` | `NeutrinoTaggerNuE.cxx:362` | `NeutrinoID_nue_functions.h:625` |
| `other_showers` | `NeutrinoTaggerNuE.cxx:963` | `NeutrinoID_nue_functions.h:209` |
| `vertex_inside_shower` | `NeutrinoTaggerNuE.cxx:1132` | `NeutrinoID_nue_functions.h:407` |

---

## No bugs found. No changes made.

---

## `stem_length`

### Purpose

Rejects events where the shower's stem segment is too long for an EM shower at the given
energy, modulo an exception for segments with many daughter tracks (muon with kinks).

### Findings

#### ✅ `calculate_num_daughter_tracks` signature — correctly adapted

Prototype: `calculate_num_daughter_tracks(vertex, sg, true)` — no graph parameter (member).  
Toolkit: `ctx.self.calculate_num_daughter_tracks(ctx.graph, vertex, sg, true, 0)` — standard toolkit API with explicit graph and mode=0.  
Same semantic: count daughter tracks from `vertex` beyond start segment `sg`, including shower segments. ✅

#### ✅ `kAvoidMuonCheck` flag — correct translation

`sg->get_flag_avoid_muon_check()` → `sg->flags_any(SegmentFlags::kAvoidMuonCheck)`. ✅

#### ✅ Cut conditions — identical

```cpp
// Both:
if (energy < 500*MeV && sg_length > 50*cm && !avoid_muon_check) {
    flag_bad = true;
    if (pair_result.first > 6 && sg_length < 55*cm) flag_bad = false;
}
```

#### ✅ All 6 fills match

| Field | Prototype | Toolkit |
|---|---|---|
| `stem_len_energy` | `max_energy/MeV` | `energy/MeV` ✅ |
| `stem_len_length` | `sg->get_length()/cm` | `sg_length/cm` ✅ |
| `stem_len_flag_avoid_muon_check` | `get_flag_avoid_muon_check()` | `flags_any(kAvoidMuonCheck)` ✅ |
| `stem_len_num_daughters` | `pair_result.first` | same ✅ |
| `stem_len_daughter_length` | `pair_result.second/cm` | same ✅ |
| `stem_len_flag` | `!flag_bad` | same ✅ |

---

## `compare_muon_energy`

### Purpose

Rejects events where a muon-range energy estimate is comparable to or larger than the shower
energy, indicating the shower candidate is more likely a muon track.

### Findings

#### ✅ `flag_numuCC` / `neutrino_type` dropped

Prototype computes:
```cpp
bool flag_numuCC = (neutrino_type >> 2) & 1U;
```
and uses it only inside a `flag_print` debug `std::cout` — never in any cut condition or fill.
The toolkit translation table (`NeutrinoTaggerNuE.cxx:45`) documents this: *"neutrino_type flag →
dropped (only appeared in debug prints)"*. **Correctly omitted.** ✅

#### ✅ `dir_shower` selection — identical to `single_shower`

Prototype and toolkit both apply the same three-way conditional:
1. `sg->get_length() > 12*cm` → segment direction at 15 cm
2. else → shower direction at 15 cm
3. Override: if `|angle_to_drift − 90| < 10` or `energy > 800*MeV` → shower direction at 25 cm

Then `.Unit()` / `.norm()`. ✅

#### ✅ Muon range function — equivalent API

Prototype: `TPCParams::get_muon_r2ke()` (ROOT `TGraph::Eval`).  
Toolkit: `ctx.particle_data->get_range_function("muon")->scalar_function(length/cm)`.  
Both map cm→MeV using the same muon range-to-kinetic-energy table. ✅

#### ✅ E_muon update loop — iterates `main_vertex`, not shower vertex

Both prototype and toolkit loop over segments at `main_vertex` (the global neutrino interaction
vertex, not the shower start vertex):

```cpp
// Prototype: map_vertex_segments[main_vertex]
// Toolkit:   boost::out_edges(ctx.main_vertex)
```

The `has_particle_info()` null guard before `pdg()` access is a defensive addition in the
toolkit; the muon/proton pdg check (`pdg == 13 || pdg == 2212`) and the dQ_dx_cut formula
`0.8866 + 0.9533 * pow(18*cm/length, 0.4234)` are identical. ✅

#### ✅ Cut conditions — identical (explicit parentheses added)

Prototype (implicit `&&` > `||` grouping):
```cpp
if (E_muon > max_energy && max_energy < 550*MeV ||
    muon_length > tmp_shower_total_length ||
    muon_length > 80*cm ||
    muon_length > 0.6*tmp_shower_total_length && max_energy < 500*MeV)
```

Toolkit (explicit parentheses):
```cpp
if ((E_muon > energy && energy < 550*MeV) ||
    muon_length > tmp_shower_total_length ||
    muon_length > 80*cm ||
    (muon_length > 0.6 * tmp_shower_total_length && energy < 500*MeV))
```
Same grouping, same logic. ✅

#### ✅ All 6 fills match

| Field | Prototype | Toolkit |
|---|---|---|
| `cme_mu_energy` | `E_muon/MeV` | same ✅ |
| `cme_energy` | `max_energy/MeV` | `energy/MeV` ✅ |
| `cme_mu_length` | `muon_length/cm` | same ✅ |
| `cme_length` | `tmp_shower_total_length/cm` | same ✅ |
| `cme_angle_beam` | `dir_beam.Angle(dir_shower)/π*180` | same formula ✅ |
| `cme_flag` | `!flag_bad` | same ✅ |

---

## `other_showers`

### Purpose

Evaluates whether other electron-like showers (direct conn_type=1, or indirect conn_type=2
within 72 cm) are energetically or geometrically inconsistent with the main shower being a
primary nue electron.

### Findings

#### ✅ Inner loop shadow `sg` → `sg1` — correctly resolved

The prototype's second loop declares a local variable `sg` that shadows the outer `sg`:

```cpp
for (...){
    WCPPID::WCShower *shower1 = *it1;
    WCPPID::ProtoSegment *sg = shower1->get_start_segment();  // shadows outer sg
    if (sg->get_particle_type()!=11) continue;
    ...
}
```

The outer `sg` is the main shower's start segment (set earlier as
`shower->get_start_segment()`). The inner `sg` is `shower1`'s start segment — a different
object. In all uses inside the loop the local `sg` is intended, so the shadowing happens to
be harmless in the prototype. The toolkit renames it to `sg1` to eliminate shadowing. ✅

#### ✅ `has_particle_info()` guard — defensive improvement

Prototype accesses `sg1->get_particle_type()` unconditionally (both loops). The toolkit adds:
```cpp
if (!sg1 || !sg1->has_particle_info() || sg1->particle_info()->pdg() != 11) continue;
```
Same filter outcome; the guard prevents a crash on segments without particle info. ✅

#### ✅ `flag_direct_max_pi0` — dead variable correctly dropped

The prototype computes `flag_direct_max_pi0` inside the `conn_type==1` max-energy update:
```cpp
if (E_shower1 > E_direct_max_energy && shower1->get_start_vertex().first == main_vertex) {
    E_direct_max_energy = E_shower1;
    if (map_shower_pio_id.find(shower1) != ...) flag_direct_max_pi0 = true;
    else flag_direct_max_pi0 = false;
}
```
`flag_direct_max_pi0` is never referenced in any cut condition or fill. The toolkit correctly
omits it. ✅

#### ✅ `n_indirect_showers` — dead counter correctly dropped (comment preserved)

Prototype increments `n_indirect_showers` for each indirect shower with `E_shower1 > 80*MeV`,
but the variable is never used in cuts or fills. The toolkit replaces the increment with a
comment:
```cpp
if (E_shower1 > 80*units::MeV) { /* n_indirect_showers — not filled, keep for logic */ }
```
✅

#### ✅ `flag_indirect_max_pi0` tracking — equivalent

Prototype:
```cpp
if (map_shower_pio_id.find(shower1) != map_shower_pio_id.end())
    flag_indirect_max_pi0 = true;
else flag_indirect_max_pi0 = false;
```
Toolkit: `flag_indirect_max_pi0 = ctx.map_shower_pio_id.count(shower1) > 0;` — equivalent. ✅

#### ✅ All cut conditions — identical

| Condition | Prototype | Toolkit |
|---|---|---|
| `flag_single_shower && max_energy > Eshower` | ✅ | ✅ |
| `flag_single_shower && Eshower < 150*MeV && total_other_energy > 0.27*Eshower` | ✅ | ✅ |
| `max_energy_1 > Eshower * 0.75` | ✅ | ✅ |
| `E_indirect_max_energy > Eshower + 350*MeV \|\| E_direct_max_energy > Eshower` | ✅ | ✅ |
| `Eshower < 1000*MeV && n_direct_showers > 0 && E_direct_max_energy > 0.33*Eshower` | ✅ | ✅ |
| `Eshower >= 1000*MeV && n_direct > 0 && E_direct_max > 0.33*Eshower && E_direct_total > 900*MeV` | ✅ | ✅ |
| `flag_indirect_max_pi0` branch (6748 + 7004 pattern) | ✅ | ✅ |

#### ✅ All 12 fills match

| Field | Prototype | Toolkit |
|---|---|---|
| `mgt_flag_single_shower` | `flag_single_shower` | same ✅ |
| `mgt_max_energy` | `max_energy/MeV` | same ✅ |
| `mgt_energy` | `Eshower/MeV` | same ✅ |
| `mgt_total_other_energy` | `total_other_energy/MeV` | same ✅ |
| `mgt_max_energy_1` | `max_energy_1/MeV` | same ✅ |
| `mgt_e_indirect_max_energy` | `E_indirect_max_energy/MeV` | same ✅ |
| `mgt_e_direct_max_energy` | `E_direct_max_energy/MeV` | same ✅ |
| `mgt_n_direct_showers` | `n_direct_showers` | same ✅ |
| `mgt_e_direct_total_energy` | `E_direct_total_energy/MeV` | same ✅ |
| `mgt_e_indirect_total_energy` | `E_indirect_total_energy/MeV` | same ✅ |
| `mgt_flag_indirect_max_pio` | `flag_indirect_max_pi0` | same ✅ |
| `mgt_flag` | `!flag_bad` | same ✅ |

---

## `vertex_inside_shower`

### Purpose

Detects two failure modes — `flag_bad1` (vis_1): a segment at the vertex is nearly anti-parallel
to the shower and similarly short (kink/broken-track topology); `flag_bad2` (vis_2): a segment
is nearly collinear with the shower direction combined with weak direction or high dQ/dx
(crossing track / broken muon).

### Findings

#### ⚠️ Prototype bug `max_sg = sg` — faithfully reproduced

The prototype's Block 1 vertex-segment loop:
```cpp
if (angle > max_angle && sg1->get_length() > 1.0*cm) {
    max_angle = angle;
    max_sg = sg;    // BUG: assigns shower start segment, not sg1
}
```

`max_sg` is always set to `sg` (the shower's start segment), regardless of which `sg1` achieved
the maximum angle. Consequently:

```cpp
double tmp_length1 = max_sg->get_length();  // = sg->get_length(), always
double tmp_length2 = shower->get_start_segment()->get_length();  // = sg->get_length(), also
// → tmp_length1 == tmp_length2 always
```

And `vis_1_particle_type` is `max_sg->get_particle_type()` = `sg->get_particle_type()`.

The toolkit reproduces this exactly:
```cpp
max_sg = sg;    // prototype assigns sg (start seg), not sg1
```
documented in the function header with the note:
> *"Prototype line 456 assigns `max_sg = sg` instead of `max_sg = sg1` — this is a prototype bug
> that is faithfully reproduced here to preserve BDT input values."*

Changing it to `max_sg = sg1` would alter `tmp_length1`, `vis_1_tmp_length1`, and
`vis_1_particle_type`, breaking BDT reproducibility. ✅

#### ✅ `max_angle1` / `max_medium_dQ_dx` in Block 2 — dead max-tracker variables omitted

The prototype tracks two extra variables alongside the max-angle segment:
```cpp
double max_angle1 = 0;     // fabs(π/2 - dir1.Angle(drift_dir))/π*180
double max_medium_dQ_dx = 0;
```
Neither `max_angle1` nor `max_medium_dQ_dx` appears in any cut condition or fill in either
Block 2 or anywhere else. The toolkit omits both from the max-update block. ✅

#### ✅ Block 1 shower-direction — equivalent formulation

Prototype recomputes `vertex_point` from wcpt-index comparison for each reference:
```cpp
if (vertex->get_wcpt().index == shower->get_start_segment()->get_wcpt_vec().front().index)
    dir1 = shower->cal_dir_3vector(front_point, 30*cm);
else
    dir1 = shower->cal_dir_3vector(back_point, 30*cm);
```

Toolkit pre-computes `vertex_point` once via `seg_endpoint_near(sg, vtx_fit_pt(vertex))` and
passes it directly:
```cpp
Vector dir1 = shower_cal_dir_3vector(*shower, vertex_point, 30*units::cm);
```
Same geometric result. ✅

#### ✅ Block 2 `dir_shower` selection — identical to `single_shower`

Same three-way conditional (segment length > 12 cm, drift-parallel override, `norm()`). ✅

#### ✅ Block 2 cut conditions — all 5 correctly parenthesized

All five `flag_bad2` conditions match the prototype's implicit `&&`>`||` grouping. Key
re-parenthesizations:

Condition 1 (6090_85_4300):
```cpp
// Prototype (implicit): (min_angle < 25 && min_weak_track==1 || min_angle < 20)
// Toolkit (explicit):   ((min_angle < 25 && min_weak_track == 1) || min_angle < 20)
```

Condition 3 (7003_1740_87003):
```cpp
// Prototype (implicit):
// (min_angle<15 && dQ_dx<2.1 || min_angle<17.5 && length<5cm && dQ_dx<2.5) &&
// (min_weak==1 && max_angle>120 || length<6cm && max_angle>135 && min_angle<12.5 && max_weak==1)
//
// Toolkit: explicit parentheses around each && sub-clause within each ||
```
All groupings match. ✅

#### ✅ `iso_angle1` pre-computation — correct

`fabs(π/2 − drift_dir.Angle(dir_shower))/π*180` appears three times in the prototype's five
conditions. The toolkit pre-computes it as `iso_angle1` (and as `angle_beam` for the beam
angle) and reuses them, making the code clearer without changing values. ✅

#### ✅ `vis_flag` return — equivalent

Prototype: `flag_bad = flag_bad1 || flag_bad2; vis_flag = !flag_bad;`  
Toolkit: `ti.vis_flag = !(flag_bad1 || flag_bad2); return flag_bad1 || flag_bad2;`  
Same result. ✅

#### ✅ Block 1 vis_1 fills (9 fields)

| Field | Prototype | Toolkit |
|---|---|---|
| `vis_1_filled` | `true` | `1.0f` ✅ |
| `vis_1_n_vtx_segs` | `map_vertex_segments[vertex].size()` | `boost::out_degree` ✅ |
| `vis_1_energy` | `Eshower/MeV` | same ✅ |
| `vis_1_num_good_tracks` | `num_good_tracks` | same ✅ |
| `vis_1_max_angle` | `max_angle` | same ✅ |
| `vis_1_max_shower_angle` | `max_shower_angle` | same ✅ |
| `vis_1_tmp_length1` | `max_sg->get_length()/cm` (= `sg->get_length()/cm`) | `segment_track_length(max_sg)/cm` ✅ |
| `vis_1_tmp_length2` | `sg->get_length()/cm` | `segment_track_length(sg)/cm` ✅ |
| `vis_1_particle_type` | `max_sg->get_particle_type()` (= `sg`'s type) | `sg->particle_info()->pdg()` ✅ |
| `vis_1_flag` | `!flag_bad1` | same ✅ |

#### ✅ Block 2 vis_2 fills (13 fields)

| Field | Prototype | Toolkit |
|---|---|---|
| `vis_2_filled` | `true` | `1.0f` ✅ |
| `vis_2_n_vtx_segs` | `map_vertex_segments[vertex].size()` | `boost::out_degree` ✅ |
| `vis_2_min_angle` | `min_angle` | same ✅ |
| `vis_2_min_weak_track` | `min_weak_track` | same ✅ |
| `vis_2_angle_beam` | `beam_dir.Angle(dir_shower)/π*180` | `angle_beam` ✅ |
| `vis_2_min_angle1` | `min_angle1` | same ✅ |
| `vis_2_iso_angle1` | `fabs(π/2−drift.Angle(shower))/π*180` | `iso_angle1` ✅ |
| `vis_2_min_medium_dQ_dx` | `min_medium_dQ_dx` | same ✅ |
| `vis_2_min_length` | `min_length/cm` | same ✅ |
| `vis_2_sg_length` | `sg->get_length()/cm` | `segment_track_length(sg)/cm` ✅ |
| `vis_2_max_angle` | `max_angle` | same ✅ |
| `vis_2_max_weak_track` | `max_weak_track` | same ✅ |
| `vis_2_flag` | `!flag_bad2` | same ✅ |

---

## Summary Table

| Check | Result | Action |
|---|---|---|
| `stem_length` — `calculate_num_daughter_tracks` API adaptation | ✅ Equivalent | — |
| `stem_length` — `kAvoidMuonCheck` flag translation | ✅ Identical | — |
| `stem_length` — cut and exception | ✅ Identical | — |
| `stem_length` — all 6 fills | ✅ Match | — |
| `compare_muon_energy` — `flag_numuCC`/`neutrino_type` dropped | ✅ Debug-print-only; correctly omitted | — |
| `compare_muon_energy` — `dir_shower` selection | ✅ Identical to `single_shower` | — |
| `compare_muon_energy` — muon range function | ✅ Equivalent API | — |
| `compare_muon_energy` — E_muon loop at `main_vertex` | ✅ Correct vertex; `has_particle_info()` guard added | — |
| `compare_muon_energy` — cut conditions parenthesized | ✅ Identical | — |
| `compare_muon_energy` — all 6 fills | ✅ Match | — |
| `other_showers` — inner loop shadow `sg` → `sg1` | ✅ Semantics preserved | — |
| `other_showers` — `has_particle_info()` guards | ✅ Defensive improvement | — |
| `other_showers` — `flag_direct_max_pi0` dead var dropped | ✅ Correctly omitted | — |
| `other_showers` — `n_indirect_showers` dead counter dropped | ✅ Correctly omitted (comment kept) | — |
| `other_showers` — all 7 cut conditions | ✅ Identical | — |
| `other_showers` — all 12 fills | ✅ Match | — |
| `vertex_inside_shower` — `max_sg = sg` prototype bug reproduced | ✅ Intentional; BDT reproducibility preserved | — |
| `vertex_inside_shower` — `max_angle1`, `max_medium_dQ_dx` dead vars omitted | ✅ Correctly omitted | — |
| `vertex_inside_shower` — Block 2 five conditions re-parenthesized | ✅ Identical logic | — |
| `vertex_inside_shower` — `iso_angle1`, `angle_beam` pre-computed | ✅ Equivalent, cleaner | — |
| `vertex_inside_shower` — Block 1 vis_1 fills (10 fields) | ✅ Match | — |
| `vertex_inside_shower` — Block 2 vis_2 fills (13 fields) | ✅ Match | — |
| `vertex_inside_shower` — `vis_flag` | ✅ Match | — |

---

## Changes Made

None. No bugs or logic divergences found.
