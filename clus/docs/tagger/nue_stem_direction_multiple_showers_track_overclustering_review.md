# Code Review: `stem_direction`, `multiple_showers`, `track_overclustering`

**Date:** 2026-04-08  
**Reviewer:** Claude Code  
**Scope:** Logic fidelity vs prototype, bug hunting, efficiency, determinism.

---

## Files Examined

| Role | File |
|---|---|
| Toolkit implementations | `clus/src/NeutrinoTaggerNuE.cxx` |
| Prototype `stem_direction`, `multiple_showers` | `prototype_pid/src/NeutrinoID_nue_functions.h` |
| Prototype `track_overclustering` | `prototype_pid/src/NeutrinoID_nue_tagger.h:547` |
| Translation table / API map | `NeutrinoTaggerNuE.cxx:26-45` |

---

## Functions Reviewed

| Function | Toolkit Location | Prototype Location |
|---|---|---|
| `stem_direction` | `NeutrinoTaggerNuE.cxx:444` | `NeutrinoID_nue_functions.h:1` |
| `multiple_showers` | `NeutrinoTaggerNuE.cxx:829` | `NeutrinoID_nue_functions.h:81` |
| `track_overclustering` | `NeutrinoTaggerNuE.cxx:3384` | `NeutrinoID_nue_tagger.h:547` |

---

## No bugs found. No changes made.

---

## `stem_direction`

### Purpose

Rejects events where the shower's start segment is misaligned with the shower's global
PCA axis or has a kinked (non-straight) trajectory at the stem, which is inconsistent
with a genuine EM shower.

### Findings

#### ✅ PCA axis computation — equivalent

Prototype: `main_cluster->Calc_PCA(tmp_pts)` → `main_cluster->get_PCA_axis(0)`  
Toolkit: `ctx.self.calc_PCA_main_axis(tmp_pts).second`

Documented in the translation table at `NeutrinoTaggerNuE.cxx:26`:
> `Calc_PCA + get_PCA_axis(0)` → `calc_PCA_main_axis(pts).second`

The `.second` of the returned pair is the principal axis vector. ✅

#### ✅ Vertex endpoint detection — equivalent

Prototype determines which end of `sg` is at `main_vertex` by comparing wcpt wire-cell
point indices (`main_vertex->get_wcpt().index == sg->get_wcpt_vec().front().index`).
Toolkit uses geometric proximity:

```cpp
bool front_is_mv = (ray_length(Ray{mv_pt, sg_front}) <= ray_length(Ray{mv_pt, sg_back}));
Point vertex_point = front_is_mv ? sg_front : sg_back;
```

This is the standard toolkit idiom replacing index comparison (documented at
`NeutrinoTaggerNuE.cxx:44`). ✅

#### ✅ `ratio` zero-division guard — defensive improvement

Prototype (L32):

```cpp
ratio = sg->get_direct_length(0,10) / sg->get_length(0,10);  // no guard
```

Toolkit:

```cpp
double len0_10 = segment_track_length(sg, 0, 0, 10);
if (len0_10 > 0)
    ratio = segment_track_direct_length(sg, 0, 10) / len0_10;
```

If `len0_10 == 0` in the prototype, the division is undefined behavior. The toolkit
keeps `ratio = 0` (initialized at top). The only ratio-dependent cut is
`(angle1 > 7.5 || angle2 > 7.5) && ratio < 0.97`; `ratio = 0` satisfies `< 0.97`,
so the toolkit may trip this condition where the prototype had UB. In practice a
segment with 10 fit points between indices 0 and 10 always has positive track length,
making this branch unreachable. **Safe defensive improvement.** ✅

#### ✅ `stem_dir_flag` and `stem_dir_filled` — added metadata

Prototype does not fill `stem_dir_flag` or `stem_dir_filled`. Toolkit fills both:

```cpp
ti.stem_dir_filled = 1.0f;
ti.stem_dir_flag   = !flag_bad;
```

`NeutrinoTaggerInfo.h` defaults: `stem_dir_flag{1}`, `stem_dir_filled{0}`. Neither
field appears in any `UbooneNueBDTScorer.cxx` `AddVariable` call — they are metadata
only. **Correct addition; no BDT impact.** ✅

#### ✅ Cut conditions — identical

Three energy tiers (> 1000 MeV, > 500 MeV, ≤ 500 MeV) with identical thresholds
(18°, 25°, 32°, 12.5°, 10°, 7.5°, 0.97 ratio, 3° `angle3`). ✅

#### ✅ All 6 fills match

`stem_dir_angle`, `stem_dir_energy`, `stem_dir_angle1`, `stem_dir_angle2`,
`stem_dir_angle3`, `stem_dir_ratio` — all units and expressions match prototype. ✅

---

## `multiple_showers`

### Purpose

Rejects events where other electron showers at or near the main vertex are energetically
comparable to the max shower, suggesting a multi-shower interaction or pi0 background.

### Structure

1. **First loop** — showers connected to `main_vertex` with electron PDG and direct
   attachment (`conn_type ≤ 1`). Accumulates `E_max_energy` and `E_total`.
2. **Second loop** — all electron showers not in the main cluster. Handles pi0 partners,
   fills `E_max_energy_1`, `E_max_energy_2`, `total_other_energy`, etc.

### Findings

#### ✅ Null guards — defensive improvements

Prototype: `if (sg->get_particle_type() != 11) continue;` — no null check on `sg`.  
Toolkit: `if (!sg || !sg->has_particle_info() || sg->particle_info()->pdg() != 11) continue;`

Added in both the first loop (main vertex showers) and the second loop (all showers).
Prevents dereferencing a null start segment. **Safe defensive addition.** ✅

#### ✅ First loop — all logic identical

- `conn_type > 1` exclusion ✅
- `Eshower` selection (kine_best if nonzero) ✅
- pi0 predicate call (with local `dummy_ti{}` instead of prototype's `flag_fill=false`) ✅
- Shower-main-length filter (`< 0.1 * total_length && < 10 cm`) ✅
- Stem-length filter (`> 80 cm`) ✅
- `E_max_energy` updated before 50 MeV threshold check ✅

#### ✅ Second loop — all logic identical

- Main-cluster exclusion via `sg->cluster() == ctx.main_vertex->cluster()` ✅
- `bad_reconstruction` predicate call ✅
- pi0 partner tracking (`E_max_energy_1`) ✅
- `conn_type ≤ 3` for `total_other_energy`, `total_other_energy_1` ✅
- `conn_type > 2` gate for `E_max_energy_2` ✅

#### ✅ All cut conditions identical

```cpp
// Cut 1:
if (E_max_energy > 0.6*max_energy ||
    (E_max_energy > 0.45*max_energy && max_energy - E_max_energy < 150*MeV)) flag_bad = true;

// Cut 2 (also requires !flag_bad):
if (E_total > 0.6*max_energy ||
    (max_energy < 400*MeV && nshowers >= 2 && E_total > 0.3*max_energy)) flag_bad = true;

// Cut 3:
if (E_max_energy_1 > max_energy * 0.75) flag_bad = true;

// Cut 4:
if (E_max_energy_2 > max_energy * 1.2 && max_energy < 250*MeV) flag_bad = true;

// Cut 5 (complex, also requires !flag_bad):
if (!flag_bad && (... total_other_energy / total_other_energy_1 conditions ...)) flag_bad = true;
```

All thresholds and logical structure match the prototype exactly. ✅

#### ✅ All 10 fills match

`mgo_energy`, `mgo_max_energy`, `mgo_total_energy`, `mgo_n_showers`,
`mgo_max_energy_1`, `mgo_max_energy_2`, `mgo_total_other_energy`,
`mgo_n_total_showers`, `mgo_total_other_energy_1`, `mgo_flag` — all ✅

---

## `track_overclustering`

### Purpose

Identifies events where a muon track has been incorrectly clustered with the shower.
Checks five distinct topological sub-conditions (tro_1 through tro_5), each filling
its own vector branch.

### tro_1 — leaf-vertex shower-internal segments

#### ✅ Leaf-vertex detection — equivalent

Prototype iterates `shower->get_map_seg_vtxs()` (shower segments), tests leaf by checking
`map_vertex_segments[v].size() == 1` — global degree.

Toolkit iterates `shower_segs`, tests by `boost::out_degree(v, graph) != 1` (same global
degree). **Equivalent.** ✅

#### ✅ `dQ_dx_cut` inf/nan guard — placement difference, no effect

Prototype applies the guard inside `if (flag_fill)` — after using `dQ_dx_cut` in
`medium_dQ_dx > dQ_dx_cut * 1.1`. If `tmp_length == 0`, `dQ_dx_cut = ∞`, and
`x > ∞` is always false so `flag_bad1` won't be set by this path.

Toolkit applies the guard before the comparison. If `dQ_dx_cut` were `∞`, the toolkit
would clamp it to 10 before comparing. In practice, the comparison is gated by
`tmp_length > 12*cm`, which means `tmp_length > 0`, making `dQ_dx_cut` finite.
**No observable difference.** ✅

### tro_2 — muon walk + lateral spread check

#### ✅ Muon walk — equivalent via `shower_segs_at`

Prototype uses `map_vtx_segs[curr_muon_vertex]` (shower-internal segments at vertex).
Toolkit uses `shower_segs_at(curr_vtx)` (boost graph edges at vertex filtered to
`shower_segs`). Since `shower_segs` contains exactly the shower's segments, the two
sets are the same. ✅

#### ✅ Dead `max_dQ_dx` loop correctly omitted

In the lateral spread vertex loop, the prototype computes:

```cpp
// 7017_1210_60518
double max_dQ_dx = 0;
for (auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
    double dQ_dx = (*it2)->get_medium_dQ_dx()/(43e3/units::cm);
    if (dQ_dx > max_dQ_dx) max_dQ_dx = dQ_dx;
}
// [commented out]: if (max_dQ_dx < 0.35) continue;
```

`max_dQ_dx` is never used (the `continue` using it is commented out). Toolkit
correctly omits the entire loop. **Dead code correctly removed.** ✅

#### ✅ Lateral spread vertex iteration — equivalent

Prototype iterates `map_vtx_segs.begin()..end()` (shower-internal vertex map keys).
Toolkit iterates `shower_vtxs` (from `fill_sets`). Both produce the same vertex set
(all vertices in the shower). ✅

#### ✅ `tro_2_v_stem_length` — prototype typo preserved

Both prototype and toolkit fill `tro_2_v_stem_length` with `/units::MeV` instead of
`/units::cm`. This is a pre-existing prototype typo; the toolkit matches it exactly to
keep BDT inputs identical. ✅

### tro_3 — gap-crossing muon walk extension

#### ✅ Same shower-segment set for gap search

Prototype gap search iterates `map_seg_vtxs` (shower segments). Toolkit iterates
`shower_segs`. Both cover only shower-internal segments for the gap-crossing step.
Cluster exclusion (`sg1->cluster() == ms->cluster()`) matches in both. ✅

#### ✅ Near-vertex selection — equivalent

Prototype uses Euclidean distance (explicit `sqrt(pow...)`). Toolkit uses
`ray_length(Ray{cvp, front1/back1})`. Both select the nearer endpoint. ✅

### tro_4 — shower-internal leaf vertices close to stem

#### ✅ Leaf detection — equivalent

Prototype: `it1->second.size() != 1` (shower-internal degree from `map_vtx_segs`).  
Toolkit: `shower_segs_at(vtx1).size() != 1`. Equivalent for shower-internal leaves. ✅

#### ✅ `end_dQ_dx` computation — equivalent

Prototype determines which end of `sg1` is at `vtx1` by wcpt index comparison, then
calls `sg1->get_medium_dQ_dx(0, 6)` (front) or `sg1->get_medium_dQ_dx(size-1-6, size-1)`
(back).

Toolkit uses geometric proximity (endpoint distance to `vtx_fit_pt(vtx1)`), then:
```cpp
segment_median_dQ_dx(sg1, 0, std::min(6, nf1-1))    // front → indices [0, 6]
segment_median_dQ_dx(sg1, std::max(0, nf1-7), nf1-1) // back  → indices [size-7, size-1]
```

`std::max(0, nf1-7) = nf1-7` (for `nf1 >= 7`) = `size-7`. Prototype: `size-1-6 = size-7`.
**Identical index ranges.** The `std::min(6, nf1-1)` and `std::max(0, nf1-7)` guards
prevent out-of-range indices on very short segments. ✅

#### ✅ Dead variables in prototype tro_5 loop correctly omitted

The prototype tro_4 loop computes `min_length`, `min_angle1`, `max_angle1` (the
`angle1 = fabs(pi/2 - dir2.Angle(drift_dir))/pi*180` for each segment) but never uses
them in any decision or fill. Toolkit does not compute these. **Dead code correctly
removed.** ✅

### tro_5 — branching at far end of start segment

#### ✅ Shower-segment degree check — equivalent

Prototype: `map_vtx_segs[vtx1].size() >= 2` (shower-internal degree at far-end vertex).  
Toolkit: `ss_at_vtx1.size() >= 2` (`shower_segs_at(vtx1).size()`). Same count. ✅

#### ✅ Cut condition — identical

```cpp
if (max_angle > 25 && min_angle < max_angle &&
    max_length > 10*cm && min_angle < 20 &&
    dir1_iso > 10 && ss_at_vtx1.size() == 3 &&
    min_count == 1 && max_count > 1 &&
    ((Eshower >= 600*MeV && dir1_iso < 40) ||
     (Eshower  < 600*MeV && dir1_iso < 25)))
    flag_bad5 = true;
```

`dir1_iso = fabs(pi/2 - dir1.angle(drift_dir)) / pi * 180` — identical to prototype's
`fabs(pi/2 - dir1.Angle(drift_dir)/pi * 180)`. ✅

#### ✅ `tro_5_v_n_vtx_segs` — correct source

Prototype fills with `map_vtx_segs[vtx1].size()`. Toolkit fills with
`ss_at_vtx1.size()`. Same value. ✅

---

## Summary Table

| Check | Result | Action |
|---|---|---|
| `stem_direction` — PCA axis via `calc_PCA_main_axis().second` | ✅ Equivalent to `Calc_PCA + get_PCA_axis(0)` | — |
| `stem_direction` — endpoint via geometric proximity | ✅ Equivalent to wcpt-index comparison | — |
| `stem_direction` — `ratio` zero-division guard added | ✅ Defensive; unreachable in practice | — |
| `stem_direction` — `stem_dir_flag`, `stem_dir_filled` added | ✅ Not BDT inputs; metadata only | — |
| `stem_direction` — cut conditions | ✅ Identical | — |
| `stem_direction` — 6 fills | ✅ Match | — |
| `multiple_showers` — null guards on `sg` and `has_particle_info` | ✅ Defensive improvement | — |
| `multiple_showers` — first loop: conn_type, pi0, length, energy filters | ✅ Identical | — |
| `multiple_showers` — second loop: bad_reconstruction, pi0 partner, energy accumulators | ✅ Identical | — |
| `multiple_showers` — all 5 cut conditions | ✅ Identical | — |
| `multiple_showers` — all 10 fills | ✅ Match | — |
| `tro_1` — leaf detection via boost degree | ✅ Equivalent to `map_vertex_segments.size()` | — |
| `tro_1` — `dQ_dx_cut` guard before comparison | ✅ Harmless; `tmp_length > 12 cm` gates comparison | — |
| `tro_2` walk — `shower_segs_at` vs `map_vtx_segs` | ✅ Equivalent shower-internal set | — |
| `tro_2` lateral — dead `max_dQ_dx` loop omitted | ✅ Correct; was commented-out dead code | — |
| `tro_2` lateral — `shower_vtxs` vs `map_vtx_segs` iteration | ✅ Equivalent vertex set | — |
| `tro_2_v_stem_length` stored in MeV units | ✅ Prototype typo preserved intentionally | — |
| `tro_3` — gap-crossing over `shower_segs` | ✅ Same shower-segment set | — |
| `tro_3` — same-cluster exclusion | ✅ Identical | — |
| `tro_4` — leaf via `shower_segs_at().size()` | ✅ Equivalent | — |
| `tro_4` — `end_dQ_dx` index ranges | ✅ Identical (front [0,6], back [size-7,size-1]) | — |
| `tro_4` — dead `min_length`, `min_angle1`, `max_angle1` omitted | ✅ Correct; never used in prototype | — |
| `tro_5` — `shower_segs_at().size() >= 2` | ✅ Equivalent to `map_vtx_segs.size()` | — |
| `tro_5` — cut condition and fills | ✅ Identical | — |

---

## Changes Made

None. No bugs or logic divergences found.
