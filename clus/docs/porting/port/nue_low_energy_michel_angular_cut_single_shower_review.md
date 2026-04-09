# Code Review: `low_energy_michel`, `angular_cut`, `single_shower`

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
| Fiducial volume check impl | `clus/src/FiducialUtils.cxx:72-77` |

---

## Functions Reviewed

| Function | Toolkit Location | Prototype Location |
|---|---|---|
| `low_energy_michel` | `NeutrinoTaggerNuE.cxx:124` | `NeutrinoID_nue_tagger.h:266` |
| `angular_cut` | `NeutrinoTaggerNuE.cxx:232` | `NeutrinoID_nue_tagger.h:436` |
| `single_shower` | `NeutrinoTaggerNuE.cxx:693` | `NeutrinoID_nue_tagger.h:320` |

---

## No bugs found. No changes made.

---

## `low_energy_michel`

### Purpose

Rejects a shower candidate that looks like a low-energy Michel electron: too short
relative to its topology, or with charge dominated by shower multiplicity rather than
MIP dQ/dx.

### Findings

#### ✅ Unused `E_range` and `n_segs` in prototype correctly omitted

The prototype computes two dead variables:

```cpp
double E_range = shower->get_kine_range();
if (E_range == 0) { E_range = shower->get_start_segment()->cal_kine_range(...); }
// ... E_range never referenced again in any condition or fill

int n_segs = shower->get_num_main_segments();
// ... n_segs not used in any condition (fill uses shower->get_num_main_segments() directly)
```

The toolkit omits both. Only the variables actually used in conditions
(`E_dQdx`, `E_charge`, `total_length`, `main_length`, `n_3seg`) and the seven fill
fields are retained. **Correct.**

#### ✅ `n_3seg` counting — equivalent logic

Prototype iterates `shower->get_map_vtx_segs()` (shower-internal vertex→segments map),
filters to the main cluster, and counts entries with ≥3 shower segments:

```cpp
Map_Proto_Vertex_Segments& map_vtx_segs = shower->get_map_vtx_segs();
for (auto it1 = map_vtx_segs.begin(); it1 != map_vtx_segs.end(); it1++) {
    if (vtx1->get_cluster_id() != shower->get_start_segment()->get_cluster_id()) continue;
    if (it1->second.size() >= 3) n_3seg++;
}
```

Toolkit uses `shower->fill_sets(shower_vtxs, shower_segs, false)` then traverses the
boost graph, counting out-edges whose segment is in `shower_segs`:

```cpp
shower->fill_sets(shower_vtxs, shower_segs, false);
for (VertexPtr vtx1 : shower_vtxs) {
    if (!vtx1->cluster() || vtx1->cluster() != start_cl) continue;
    if (!vtx1->descriptor_valid()) continue;
    int deg = 0;
    for (auto [eit, eend] = boost::out_edges(vtx1->get_descriptor(), ctx.graph); ...)
        if (shower_segs.count(ctx.graph[*eit].segment)) ++deg;
    if (deg >= 3) ++n_3seg;
}
```

`fill_sets` populates `shower_segs` with exactly the shower's segments (same set as
`map_vtx_segs` values). The filtered boost traversal counts the same shower-segment
degree per vertex. **Equivalent.** The `descriptor_valid()` guard is a defensive
addition for robustness only.

#### ✅ Short-shower criterion — identical

```cpp
// Both prototype and toolkit:
if ((total_length < 25*cm && main_length > 0.75*total_length && n_3seg == 0) ||
    (total_length < 18*cm && main_length > 0.75*total_length && n_3seg >  0))
    flag_bad = true;
```

#### ✅ Low-charge MIP criterion — identical

```cpp
if (E_charge < 100*MeV && E_dQdx < 0.7*E_charge &&
    shower->get_num_segments() == shower->get_num_main_segments())
    flag_bad = true;
```

#### ✅ All fills match

| Field | Prototype | Toolkit |
|---|---|---|
| `lem_shower_total_length` | `total_length/cm` | same ✅ |
| `lem_shower_main_length` | `get_total_length(cluster)/cm` | `main_length/cm` ✅ |
| `lem_n_3seg` | `n_3seg` | same ✅ |
| `lem_e_charge` | `E_charge/MeV` | same ✅ |
| `lem_e_dQdx` | `E_dQdx/MeV` | same ✅ |
| `lem_shower_num_segs` | `get_num_segments()` | same ✅ |
| `lem_shower_num_main_segs` | `get_num_main_segments()` | same ✅ |
| `lem_flag` | `!flag_bad` | same ✅ |

---

## `angular_cut`

### Purpose

Rejects events where track material is predominantly backward relative to the shower
direction (suggesting a hadronic interaction product, not an electron), or where
shower vertices lie outside the fiducial volume.

### Findings

#### ✅ Segment-loop variable shadowing — correctly resolved

Prototype has a loop-local variable also named `angle` that shadows the function
parameter `angle`:

```cpp
bool WCPPID::NeutrinoID::angular_cut(..., double angle, ...) {
    for (...) {
        double angle = dir1.Angle(dir_beam)/3.1415926*180.;  // shadows outer `angle`
        if (angle > 90) ...  // uses local `angle` — CORRECT in prototype
    }
    // ... uses outer `angle` for cut conditions
}
```

The prototype happens to work because the inner `angle` is in a nested block and
C++ lookup resolves to the local first. The toolkit renames the loop variable to
`seg_angle` to eliminate the shadowing, which preserves identical semantics. ✅

#### ✅ Cut conditions — identical after re-parenthesization

The prototype has one long `if (A || B || C || D || E)` expression across five lines.
With C++ operator precedence (`&&` tighter than `||`), each clause groups naturally.
The toolkit splits these into five named boolean variables (`cut_1`…`cut_5`) with
explicit parentheses matching the prototype's implicit grouping exactly. ✅

Specifically `cut_4`:

```cpp
// Prototype (implicit grouping):
energy < 650*MeV &&
    (acc_forward_length==0 || ... < 0.03*acc_backward_length) &&
    acc_backward_length > 0 &&
    (acc_backward_length - shower_main_length > acc_forward_length && angle > 90 || angle <= 90)

// Toolkit (explicit):
bool cut_4 = (energy < 650*MeV &&
              (acc_forward_length == 0 || ...) &&
              acc_backward_length > 0 &&
              ((acc_backward_length - shower_main_length > acc_forward_length && angle > 90) ||
               angle <= 90));
```

`shower->get_total_length(shower->get_start_segment()->get_cluster_id())` in the
prototype → `shower_main_length` in the toolkit (same expression, pre-computed). ✅

#### ✅ `dir_shower` — full-shower direction at 30 cm — correct

Prototype: `shower->cal_dir_3vector(vertex_point, 30*cm)`  
Toolkit: `shower_cal_dir_3vector(*shower, vertex_point, 30*cm)` ✅

Note: unlike `single_shower` which has a conditional (start segment length > 12 cm),
`angular_cut` always uses the shower-level direction for `dir_shower`. Toolkit matches. ✅

#### ✅ Vertex-segment loop — no start-segment exclusion (matches prototype)

`angular_cut` iterates all segments at the vertex, including the shower's start segment.
This differs from `single_shower` (which skips `sg`). Both toolkit implementations
correctly match their respective prototypes.

#### ✅ Fiducial check — `stm_tol_vec` has no effect

Prototype:
```cpp
std::vector<double> stm_tol_vec = {-1.5*cm, -1.5*cm, -1.5*cm, -1.5*cm, -1.5*cm};
fid->inside_fiducial_volume(vtx1->get_fit_pt(), offset_x, &stm_tol_vec)
```

Toolkit:
```cpp
fiducial_utils->inside_fiducial_volume(vtx_fit_pt(vtx1))  // no tolerance
```

The toolkit's `FiducialUtils::inside_fiducial_volume` implementation ignores the
`tolerance_vec` parameter entirely (`FiducialUtils.cxx:75`: *"currently tolerance
vector is not used"*) and delegates straight to `m_sd.fiducial->contained(p)`.
The toolkit's no-tolerance call is therefore numerically identical to the prototype's
`stm_tol_vec` call. ✅

#### ✅ Fiducial check vertex iteration — equivalent

Prototype iterates `shower->get_map_vtx_segs()` keys, skipping start vertex and
off-cluster vertices. Toolkit iterates `fv_vtxs` from `fill_sets` with the same
filters. Both produce the same vertex set. ✅

#### ✅ All fills match

| Field | Prototype | Toolkit |
|---|---|---|
| `anc_energy` | `energy/MeV` | same ✅ |
| `anc_angle` | `angle` | same ✅ |
| `anc_max_angle` | `max_angle` | same ✅ |
| `anc_max_length` | `max_length/cm` | same ✅ |
| `anc_acc_forward_length` | `acc_forward_length/cm` | same ✅ |
| `anc_acc_backward_length` | `acc_backward_length/cm` | same ✅ |
| `anc_acc_forward_length1` | `acc_forward_length1/cm` | same ✅ |
| `anc_shower_main_length` | `main_length/cm` | `shower_main_length/cm` ✅ |
| `anc_shower_total_length` | `total_length/cm` | `shower_total_length/cm` ✅ |
| `anc_flag_main_outside` | `flag_main_outside` | same ✅ |
| `anc_flag` | `!flag_bad` | same ✅ |

---

## `single_shower`

### Purpose

Evaluates whether a single-shower topology passes geometric and dQ/dx quality cuts
for the nue BDT. The `flag_single_shower` parameter selects between two distinct cut
sets — `true` for a truly isolated shower (no vertex-connected tracks), `false` for a
shower with vertex-connected tracks present.

### Findings

#### ✅ `Eshower` selection — identical

```cpp
// Both: kine_best if non-zero, else kine_charge
Eshower = (shower->get_kine_best() != 0) ? shower->get_kine_best()
                                          : shower->get_kine_charge();
```

#### ✅ `max_dQ_dx` stem scan — identical

Both iterate `get_stem_dQ_dx(vertex, sg, 20)`, break after index 2, take the maximum.
The values are already normalized by the shower API (no `/(43e3/cm)` factor needed),
consistent with the dimensionless comparison thresholds (3.0, 2.4). ✅

#### ✅ `dir_shower` conditional — identical

Prototype and toolkit both select the direction source based on start-segment length:

```cpp
if (segment_track_length(sg) > 12*cm)   // sg->get_length() in prototype
    dir_shower = segment direction at 15cm;
else
    dir_shower = shower direction at 15cm;
// Override if nearly orthogonal to drift or high energy:
if (|angle_to_drift - 90| < 10 || Eshower > 800*MeV)
    dir_shower = shower direction at 25cm;
dir_shower = dir_shower.norm();   // .Unit() in prototype
```

#### ✅ `num_valid_tracks` loop — identical

Prototype iterates `map_vertex_segments[vertex]`, toolkit iterates boost out-edges.
Both skip the shower start segment (`sg`). For each other segment:

1. Non-shower path (`!get_flag_shower()` / `!seg_is_shower()`): count if not dir-weak,
   or dir-weak with sufficient length, or sufficient `medium_dQ_dx`. ✅
2. Shower-flagged fallback: count if `medium_dQ_dx > dQ_dx_cut` with the same
   `0.8866 + 0.9533 * pow(18*cm/length, 0.4234)` formula. ✅

`segment_median_dQ_dx(sg1)/(43e3/cm)` matches prototype's
`sg1->get_medium_dQ_dx()/(43e3/cm)` — same normalization. ✅

#### ✅ `flag_single_shower == true` branch — identical

Three sub-conditions (6572_18_948, near-vertical/horizontal dQ/dx, drift-parallel
dQ/dx) match prototype exactly in thresholds and logical structure. ✅

#### ✅ `flag_single_shower == false` branch — identical

```cpp
// Both:
if (num_valid_tracks == 0 && angle_beam > 60 && n_vtx_segs <= 3)  flag_bad = true;
if (Eshower < 200*MeV && angle_vertical < 10 && angle_drift < 5 &&
    max_length < 5*cm && num_valid_tracks <= 1)                     flag_bad = true;
```

`n_vtx_segs = boost::out_degree(vertex, graph)` replaces prototype's
`map_vertex_segments[vertex].size()`. Both include the shower start segment in the
count (no exclusion here). ✅

#### ✅ All fills match

| Field | Prototype | Toolkit |
|---|---|---|
| `spt_flag_single_shower` | `flag_single_shower` | same ✅ |
| `spt_energy` | `Eshower/MeV` | same ✅ |
| `spt_shower_main_length` | `get_total_length(cluster)/cm` | `shower_main_length/cm` ✅ |
| `spt_shower_total_length` | `get_total_length()/cm` | `shower_total_length/cm` ✅ |
| `spt_angle_beam` | `angle_beam` | same ✅ |
| `spt_angle_vertical` | `angle_vertical` | same ✅ |
| `spt_max_dQ_dx` | `max_dQ_dx` | same ✅ |
| `spt_angle_beam_1` | `dir_shower1.Angle(dir_beam)*180/pi` | same formula ✅ |
| `spt_angle_drift` | `|pi/2-angle_to_drift|*180/pi` | same ✅ |
| `spt_angle_drift_1` | same with dir_shower1 | same ✅ |
| `spt_num_valid_tracks` | `num_valid_tracks` | same ✅ |
| `spt_n_vtx_segs` | `map_vertex_segments[vertex].size()` | `boost::out_degree` ✅ |
| `spt_max_length` | `max_length/cm` | same ✅ |
| `spt_flag` | `!flag_bad` | same ✅ |

---

## Summary Table

| Check | Result | Action |
|---|---|---|
| `low_energy_michel` — unused `E_range`, `n_segs` omitted | ✅ Correctly dropped | — |
| `low_energy_michel` — `n_3seg` via boost graph filtered by shower_segs | ✅ Equivalent to prototype map | — |
| `low_energy_michel` — short-shower / MIP criteria | ✅ Identical | — |
| `low_energy_michel` — all 8 fills | ✅ Match | — |
| `angular_cut` — loop-local `angle` shadow → renamed `seg_angle` | ✅ Semantics preserved | — |
| `angular_cut` — 5 cut conditions, re-parenthesized | ✅ Identical logic | — |
| `angular_cut` — `dir_shower` = shower direction at 30 cm (no segment condition) | ✅ Matches prototype | — |
| `angular_cut` — no start-segment exclusion in vertex loop | ✅ Matches prototype | — |
| `angular_cut` — `stm_tol_vec` ignored by `inside_fiducial_volume` impl | ✅ No-tolerance call equivalent | — |
| `angular_cut` — all 11 fills | ✅ Match | — |
| `single_shower` — `Eshower` selection | ✅ Identical | — |
| `single_shower` — `max_dQ_dx` stem scan (first 3, max) | ✅ Identical | — |
| `single_shower` — `dir_shower` conditional on segment length | ✅ Identical | — |
| `single_shower` — `num_valid_tracks` loop, skips `sg`, same normalization | ✅ Identical | — |
| `single_shower` — both flag_single_shower branches | ✅ Identical | — |
| `single_shower` — `n_vtx_segs` via `boost::out_degree` | ✅ Equivalent to map size | — |
| `single_shower` — all 14 fills | ✅ Match | — |

---

## Changes Made

None. No bugs or logic divergences found.
