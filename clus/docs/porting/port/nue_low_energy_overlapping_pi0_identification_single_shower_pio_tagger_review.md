# Code Review: `low_energy_overlapping`, `pi0_identification`, `single_shower_pio_tagger`

**Date:** 2026-04-09  
**Reviewer:** Claude Code  
**Scope:** Logic fidelity vs prototype, bug hunting, efficiency, determinism.

---

## Files Examined

| Role | File |
|---|---|
| Toolkit implementations | `clus/src/NeutrinoTaggerNuE.cxx` |
| Prototype implementations | `prototype_pid/src/NeutrinoID_nue_tagger.h` |
| Translation table / API map | `NeutrinoTaggerNuE.cxx:30-45` |
| Units definitions | `util/inc/WireCellUtil/Units.h` |

---

## Functions Reviewed

| Function | Toolkit Location | Prototype Location |
|---|---|---|
| `low_energy_overlapping` | `NeutrinoTaggerNuE.cxx:1972` | `NeutrinoID_nue_tagger.h:2485` |
| `pi0_identification` | `NeutrinoTaggerNuE.cxx:548` | `NeutrinoID_nue_tagger.h:2643` |
| `single_shower_pio_tagger` | `NeutrinoTaggerNuE.cxx:2529` | `NeutrinoID_nue_tagger.h:2756` |

---

## No bugs found. No changes made.

One prototype typo documented (`800*units::cm` → kept as `800*units::MeV` in toolkit).

---

## `low_energy_overlapping`

### Purpose

Three complementary checks for low-energy shower overlap:
1. `flag_overlap_1` — shower-internal vertex with 2 shower segments forming a narrow opening angle.
2. `flag_overlap_2` — short collinear muon or weak segment at the shower start vertex.
3. `flag_overlap_3` — backward/isolated low-energy shower with outward-pointing hits.

### Findings

#### ✅ `nseg` counting — equivalent

Prototype iterates `map_seg_vtxs` (segment→vertices map, shower-internal), counts segments in the
main cluster:

```cpp
int nseg = 0;
for (auto it = map_seg_vtxs.begin(); it!= map_seg_vtxs.end(); it++){
    if (it->first->get_cluster_id() == sg->get_cluster_id()) nseg ++;
}
```

Toolkit iterates `shower_segs` from `fill_sets`:

```cpp
for (SegmentPtr sg1 : shower_segs) {
    if (sg1->cluster() == sg->cluster()) ++nseg;
}
```

`shower_segs` is the same set of shower-internal segments as `map_seg_vtxs` keys. ✅

#### ✅ `n_valid_tracks` and `min_angle` — equivalent

Prototype iterates `map_vertex_segments[vtx]` (all segments at start vertex), skips `sg`,
checks `!is_dir_weak() || particle_type==2212 || length > 20cm`. Toolkit uses `boost::out_edges`
with the same skip and the same three conditions:

```cpp
bool is_proton = sg1->has_particle_info() && sg1->particle_info()->pdg() == 2212;
if ((!sg1->dir_weak() || is_proton || segment_track_length(sg1) > 20*units::cm) && !seg_is_shower(sg1))
    ++n_valid_tracks;
```

`segment_track_length` → `sg->get_length()` per translation table L33. ✅

#### ✅ `n_out / n_sum` — equivalent

Prototype iterates `map_vtx_segs` (shower vertices) and `map_seg_vtxs` (shower segments),
using interior points only (`for i=1; i+1<pts.size()`). Toolkit iterates `shower_vtxs` and
`shower_segs`, with `for (size_t i = 1; i + 1 < fits1.size(); ++i)`. Same interior-points-only
logic. Both use the 15 cm direction vector from the stem. ✅

#### ✅ `flag_overlap_1` loop — open-angle argument commutative

Prototype uses `*it->second.begin()` and `*it->second.rbegin()` on an ordered `std::set` of
shower segments at that vertex. Toolkit collects shower segments into a `std::vector vtx_ss` and
uses `vtx_ss.front()` / `vtx_ss.back()`. The ordering differs, but the angle computation
`dv1.angle(dv2)` is symmetric, so the result is identical when `vtx_ss.size() == 2`. ✅

Both only fill for main-cluster vertices (`cluster == sg->cluster()`). ✅

#### ✅ `n_vtx_segs_global` — `boost::out_degree` replaces `map_vertex_segments[vtx].size()`

Used in `flag_overlap_1` (cut: `== 1`), `flag_overlap_2` (cut: `> 1` and `== 2`), and
`flag_overlap_3` (cuts: `> 1` and `== 1`). `boost::out_degree(vtx, graph)` gives the same
count as `map_vertex_segments[vtx].size()`. ✅

#### ✅ `lol_2_v_type` fill — `pdg()` with null guard

Prototype: `(*it)->get_particle_type()` — returns 0 when no particle type is set.
Toolkit: `sg1->has_particle_info() ? sg1->particle_info()->pdg() : 0` — same default. ✅

#### ✅ All fills match

| Field | Prototype | Toolkit |
|---|---|---|
| `lol_1_v_energy` | `Eshower/MeV` | same ✅ |
| `lol_1_v_vtx_n_segs` | `map_vertex_segments[vtx].size()` | `n_vtx_segs_global` ✅ |
| `lol_1_v_nseg` | `nseg` | same ✅ |
| `lol_1_v_angle` | `dir1.Angle(dir2)/π*180` | `dv1.angle(dv2)/π*180` ✅ |
| `lol_1_v_flag` | `!flag_overlap_1` | `!flag_ov1` ✅ |
| `lol_2_v_flag` | `!flag_overlap_2` | `!flag_ov2` ✅ |
| `lol_2_v_length` | `(*it)->get_length()/cm` | `len1/cm` ✅ |
| `lol_2_v_angle` | `dir1.Angle(dir2)/π*180` | `ang2` ✅ |
| `lol_2_v_type` | `get_particle_type()` | `pdg()` or 0 ✅ |
| `lol_2_v_vtx_n_segs` | `map_vertex_segments[vtx].size()` | `n_vtx_segs_global` ✅ |
| `lol_2_v_energy` | `Eshower/MeV` | same ✅ |
| `lol_2_v_shower_main_length` | `shower->get_total_length(cluster)/cm` | `main_len/cm` ✅ |
| `lol_2_v_flag_dir_weak` | `is_dir_weak()` | `dir_weak()` ✅ |
| `lol_3_flag` | `!flag_overlap_3` | same ✅ |
| `lol_3_angle_beam` | `angle_beam` | same ✅ |
| `lol_3_min_angle` | `min_angle` | `min_angle_vtx` ✅ |
| `lol_3_n_valid_tracks` | `n_valid_tracks` | same ✅ |
| `lol_3_vtx_n_segs` | `map_vertex_segments[vtx].size()` | `n_vtx_segs_global` ✅ |
| `lol_3_energy` | `Eshower/MeV` | same ✅ |
| `lol_3_shower_main_length` | `shower->get_total_length(cluster)/cm` | `main_length/cm` ✅ |
| `lol_3_n_sum` | `n_sum` | same ✅ |
| `lol_3_n_out` | `n_out` | same ✅ |
| `lol_flag` | `!flag_overlap` | same ✅ |

---

## `pi0_identification`

### Purpose

Determines whether a shower is likely the photon from a pi0 decay, via two branches:

- `flag_pi0_1`: shower is in the reconstructed pi0 map (`map_shower_pio_id`); checks mass and
  energy balance of the pair.
- `flag_pi0_2`: shower is NOT in the pi0 map; looks for another cluster in the anti-parallel
  direction at ≤36 cm (back-to-back geometry).

### Findings

#### ✅ `used_vertices` collection — equivalent

Prototype iterates `map_vtx_segs` of each known pi0 shower to collect vertices:

```cpp
for (auto it1 = map_shower_pio_id.begin(); ...) {
    Map_Proto_Vertex_Segments& map_vtx_segs = shower1->get_map_vtx_segs();
    for (auto it2 = map_vtx_segs.begin(); ...) used_vertices.insert(it2->first);
}
```

Toolkit uses `fill_sets(vtxs, segs, false)`:

```cpp
for (auto& [shower1, pio_id] : ctx.map_shower_pio_id) {
    IndexedVertexSet vtxs; IndexedSegmentSet segs;
    shower1->fill_sets(vtxs, segs, false);
    used_vertices.insert(vtxs.begin(), vtxs.end());
}
```

`fill_sets` populates `vtxs` with exactly the same shower-internal vertex set as
`map_vtx_segs` keys. ✅

#### ✅ `flag_pi0_1` branch — mass/energy conditions identical

Both prototype and toolkit:
1. Check `|mass - 135 MeV| < 35 MeV && type==1` OR `|mass - 135 MeV| < 60 MeV && type==2`.
2. If true: set `flag_pi0_1 = true` if `min(E1,E2) > 15 MeV && |E1-E2|/(E1+E2) < 0.87`.
3. Also set `flag_pi0_1 = true` if `min(E1,E2) > max(10 MeV, threshold) && max(E1,E2) < 400 MeV`.
4. Veto (7049_875_43775): clear `flag_pi0_1` if the pair is asymmetric and far-separated. ✅

#### ✅ `dis1` / `dis2` — equivalent

Prototype uses explicit `sqrt(pow(...))`. Toolkit uses `ray_length(Ray{start_point, vertex_point})`.
Both compute the Euclidean distance from shower start point to vertex. ✅

#### ✅ `flag_pi0_2` branch — global vertex iteration equivalent

Prototype iterates `map_vertex_segments` (all event vertices), skipping same-cluster vertices and
`used_vertices`. Toolkit iterates `graph_nodes(ctx.graph)`, applying the same two filters. ✅

#### ✅ `acc_length` per other cluster — equivalent

Prototype inner loop iterates all `map_segment_vertices` entries, summing `sg1->get_length()` for
segments in `vtx1`'s cluster. Toolkit precomputes `cluster_acc_length` using `boost::edges` and
`segment_track_length` (`= sg->get_length()` per translation table L33). ✅

#### ✅ Back-angle computation — identical

Both compute `back_angle = 180 - dir1.angle(dir2) / π * 180`. ✅

#### ✅ `pio_1_flag` placement — inside `flag_pi0_1` branch only

Prototype sets `pio_1_flag = !flag_pi0_1` inside the `it != map_shower_pio_id.end()` block
(L2749). Toolkit does the same at L618. The field retains its default value when the shower is
not in the pi0 map. ✅

#### ✅ `pio_2_v_*` fills — guarded by `dis2 > 0` in both branches

Prototype fills both the "is pi0" and "is not pi0" branches only when `dir2.Mag() > 0`.
Toolkit: `if (dis2 <= 0) continue;` guards both branches identically. ✅

#### ✅ All fills match

| Field | Prototype | Toolkit |
|---|---|---|
| `pio_flag_pio` | `(it != end)` | same ✅ |
| `pio_1_mass` | `mass_pair.first/MeV` | same ✅ |
| `pio_1_pio_type` | `mass_pair.second` | same ✅ |
| `pio_1_energy_1` | `Eshower_1/MeV` | same ✅ |
| `pio_1_energy_2` | `Eshower_2/MeV` | same ✅ |
| `pio_1_dis_1` | `dis1/cm` | same ✅ |
| `pio_1_dis_2` | `dis2/cm` | same ✅ |
| `pio_1_flag` | `!flag_pi0_1` (in pi0-map branch) | same ✅ |
| `pio_2_v_flag` | `false`/`true` | `0.0f`/`1.0f` ✅ |
| `pio_2_v_dis2` | `dir2.Mag()/cm` | `dis2/cm` ✅ |
| `pio_2_v_angle2` | `180 - angle/π*180` | `back_angle` ✅ |
| `pio_2_v_acc_length` | `acc_length/cm` | same ✅ |

---

## `single_shower_pio_tagger`

### Purpose

Identifies events where the main shower is likely a pi0 photon rather than a nue electron, via
two checks:

- `flag_bad1`: another pdg-11 shower at the same vertex (type==2) points away from the vertex
  in an angle < 30° with the shower direction (suggesting an off-vertex gamma).
- `flag_bad2`: the farthest-forward shower vertex has a MIP-like dQ/dx, while the shower start
  is clean — suggesting a broken track + hadronic tail.

### Findings

#### ✅ `flag_bad1` shower loop — identical

Prototype iterates `map_vertex_to_shower[vtx]`. Toolkit uses `ctx.map_vertex_to_shower.find(vtx)`.
Both:
1. Skip non-pdg-11 start segments.
2. Skip `shower1 == shower`.
3. Skip type > 2.
4. Fill only for type == 2. ✅

#### ✅ `dir1` and `dir2` computation — equivalent

Prototype:
```cpp
TVector3 dir1(sp1.x - vtx->get_fit_pt().x, ...);
TVector3 dir2 = shower1->cal_dir_3vector(sp1, 15*cm);
```
Toolkit:
```cpp
Vector dir1 = sp1 - vtx_fit_pt(vtx);
Vector dir2 = shower_cal_dir_3vector(*shower1, sp1, 15*units::cm);
```
✅

#### ✅ `max_vtx` search — farthest forward vertex in main cluster

Prototype iterates `map_vtx_segs`, projects onto `dir` (shower direction at 15 cm):

```cpp
double dis = dir1.Dot(dir);  // dir1 = vtx - vertex_point, dir = shower_dir
if (dis > max_dis) { max_dis = dis; max_vtx = tmp_vtx; }
```

Toolkit:
```cpp
double dis = dir1.dot(shower_dir);
if (dis > max_dis) { max_dis = dis; max_vtx = tmp_vtx; }
```

Same forward-projection criterion. ✅

#### ✅ `max_sg` selection — most anti-aligned shower segment at `max_vtx`

Prototype iterates `map_vtx_segs[max_vtx]` (shower-internal segments at `max_vtx`). Toolkit
iterates `boost::out_edges(max_vtx, graph)` filtered to `shower_segs`. Both maximize
`dir1.Angle(dir)` / `dir1.angle(shower_dir)`. ✅

#### ✅ `medium_dQ_dx` at far end — proximity replaces wcpt-index comparison

Prototype detects near-end of `max_sg` by comparing `wcpt` indices. Toolkit uses geometric
proximity of fit endpoints to `vtx_fit_pt(max_vtx)`. Functionally equivalent (same translation
strategy as used in `high_energy_overlapping`). ✅

#### ✅ `start_dQ_dx` — per-fit-point dQ/dx max over first 3 points

Prototype iterates `vec_dQ`/`vec_dx` arrays. Toolkit computes `it->dQ / (it->dx + 1e-9)` per
fit entry in the same direction. Both cap at 3 points (`ncount > 2` break). ✅

#### ⚠️ `800*units::cm` — prototype typo; toolkit has correct `800*units::MeV`

Prototype line 2882:
```cpp
if (Eshower < 800*units::cm && flag_single_shower)
```

Toolkit line 2654:
```cpp
if (Eshower < 800*units::MeV && flag_single_shower) {
```

In WireCell's unit system: `units::mm = 1` (base), `units::cm = 10`, `units::MeV = 1`.
Therefore:
- Prototype threshold: `800 × 10 = 8000` in raw units → fires for Eshower < 8 GeV (essentially always).
- Toolkit threshold: `800 × 1 = 800` in raw units → fires for Eshower < 800 MeV.

The prototype's `units::cm` is almost certainly a typo for `units::MeV`. The surrounding context
(function is a low-energy pi0 tagger, cut at 250/500 MeV in the preceding block) and the
comment "7020_1108_55428" confirm the intent is an energy cut. The toolkit's `800*units::MeV`
is the physically correct interpretation.

**Decision: Keep the toolkit's `800*units::MeV`.** Adopting the prototype's bug would cause the
dQ/dx 2.0 cut to fire on all events regardless of energy, which is unlikely to be intended and
would differ from the 800 MeV intent visible in the surrounding conditions.

#### ✅ All fills match

| Field | Prototype | Toolkit |
|---|---|---|
| `sig_1_v_angle` | `dir1.Angle(dir2)/π*180` | `ang` ✅ |
| `sig_1_v_flag_single_shower` | `flag_single_shower` | same ✅ |
| `sig_1_v_energy` | `Eshower/MeV` | same ✅ |
| `sig_1_v_energy_1` | `Eshower1/MeV` | `E1/MeV` ✅ |
| `sig_1_v_flag` | `!flag_bad1` | same ✅ |
| `sig_2_v_energy` | `Eshower/MeV` | same ✅ |
| `sig_2_v_shower_angle` | `shower_angle` | same ✅ |
| `sig_2_v_flag_single_shower` | `flag_single_shower` | same ✅ |
| `sig_2_v_medium_dQ_dx` | `medium_dQ_dx` | same ✅ |
| `sig_2_v_start_dQ_dx` | `start_dQ_dx` | same ✅ |
| `sig_2_v_flag` | `!flag_bad2` | same ✅ |
| `sig_flag` | `!flag_bad` | same ✅ |

---

## Summary Table

| Check | Result | Action |
|---|---|---|
| `low_energy_overlapping` — `nseg` via `shower_segs` filter | ✅ Equivalent to `map_seg_vtxs` count | — |
| `low_energy_overlapping` — `n_valid_tracks` / `min_angle` via `boost::out_edges` | ✅ Identical logic | — |
| `low_energy_overlapping` — `n_out/n_sum` interior-points-only | ✅ Identical logic | — |
| `low_energy_overlapping` — `flag_overlap_1` open angle (symmetric, 2-segment case) | ✅ Order-independent | — |
| `low_energy_overlapping` — `n_vtx_segs_global` via `boost::out_degree` | ✅ Equivalent | — |
| `low_energy_overlapping` — `lol_2_v_type` null guard | ✅ Default 0 matches prototype | — |
| `low_energy_overlapping` — all 22 fills | ✅ Match | — |
| `pi0_identification` — `used_vertices` via `fill_sets` | ✅ Equivalent to `map_vtx_segs` traversal | — |
| `pi0_identification` — `flag_pi0_1` mass/energy/veto conditions | ✅ Identical | — |
| `pi0_identification` — `dis1`/`dis2` via `ray_length` | ✅ Equivalent | — |
| `pi0_identification` — global vertex iteration for `flag_pi0_2` | ✅ Equivalent | — |
| `pi0_identification` — `acc_length` via precomputed `cluster_acc_length` | ✅ Equivalent | — |
| `pi0_identification` — `pio_1_flag` placement inside `flag_pi0_1` branch | ✅ Matches prototype | — |
| `pi0_identification` — `pio_2_v_*` fills guarded by `dis2 > 0` | ✅ Identical | — |
| `pi0_identification` — all 12 fills | ✅ Match | — |
| `single_shower_pio_tagger` — `flag_bad1` pdg-11 type-2 shower loop | ✅ Identical | — |
| `single_shower_pio_tagger` — `max_vtx` forward-projection search | ✅ Identical | — |
| `single_shower_pio_tagger` — `max_sg` max-angle selection via `shower_segs` | ✅ Identical | — |
| `single_shower_pio_tagger` — `medium_dQ_dx` near-end via proximity | ✅ Equivalent (same design as `high_energy_overlapping`) | — |
| `single_shower_pio_tagger` — `start_dQ_dx` max of first 3 per-fit values | ✅ Identical | — |
| `single_shower_pio_tagger` — `800*units::cm` prototype typo | ⚠️ Prototype bug; toolkit correct at `800*units::MeV` | Kept as `800*units::MeV` |
| `single_shower_pio_tagger` — all 12 fills | ✅ Match | — |

---

## Changes Made

None. No toolkit bugs found. One prototype typo (`800*units::cm` for `800*units::MeV`) noted and
kept corrected in the toolkit.
