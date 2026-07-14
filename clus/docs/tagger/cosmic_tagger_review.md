# Code Review: `cosmic_tagger()` and daughter functions

**Date:** 2026-04-08  
**Reviewer:** Claude Code  
**Scope:** Logic fidelity vs prototype, bug hunting, efficiency, determinism.

---

## Files Examined

| Role | File |
|---|---|
| Toolkit implementation | `clus/src/NeutrinoTaggerCosmic.cxx` |
| Toolkit daughter functions | `clus/src/NeutrinoTrackShowerSep.cxx:164-325` |
| Prototype cosmic tagger | `prototype_pid/src/NeutrinoID_cosmic_tagger.h` |
| Prototype daughter functions | `prototype_pid/src/NeutrinoID_track_shower.h:688-763, 2372-2441` |
| Prototype bad_reconstruction | `prototype_pid/src/NeutrinoID_nue_tagger.h:3450` |

---

## Functions Reviewed

| Function | Toolkit Location |
|---|---|
| `bad_reconstruction()` | `NeutrinoTaggerCosmic.cxx:101-453` |
| `cosmic_tagger()` | `NeutrinoTaggerCosmic.cxx:475-1329` |
| `calculate_num_daughter_showers()` | `NeutrinoTrackShowerSep.cxx:164-215` |
| `calculate_num_daughter_tracks()` | `NeutrinoTrackShowerSep.cxx:221-269` |
| `find_cont_muon_segment_nue()` | `NeutrinoTrackShowerSep.cxx:275-325` |

---

## BUG FIXED — `shower_energy()` redundant helper

**Location:** `NeutrinoTaggerCosmic.cxx:85-88` (pre-fix), with 6 call sites.

**Problem:**

```cpp
static inline double shower_energy(ShowerPtr shower) {
    double e = shower->get_kine_best();
    return (e != 0) ? e : shower->get_kine_charge();
}
```

`PRShower::get_kine_best()` (PRShower.h:172-173) already performs this fallback internally:

```cpp
double get_kine_best(){
    if (data.kenergy_best != 0) return data.kenergy_best;
    else return data.kenergy_charge;
}
```

So `shower_energy(shower)` was always identically equal to `shower->get_kine_best()`. The helper was dead code, identical to the same redundancy fixed in `NeutrinoKinematics.cxx`.

**Fix:** Removed the `shower_energy()` helper and replaced all 6 call sites with `shower->get_kine_best()`:
- `bad_reconstruction` L112 (pre-fix)
- `cosmic_tagger` L611, L783, L851, L961, L1223 (pre-fix)

---

## `bad_reconstruction()`

### Structure (both prototype and toolkit)

Three independent sub-checks (`flag_bad_shower_1/2/3`), any of which marks the shower as bad.

| Sub-check | Condition |
|---|---|
| 1 | Stem length > 80 cm, or non-shower stem > 10 cm with low energy/few segments |
| 2 | Longest muon-like track extension (via `find_cont_muon_segment_nue`) exceeds energy-dependent threshold |
| 3 | Long straight track near start-segment far end, no shower topology, single main segment |

### Findings

#### ✅ Logic fidelity — sub-check 1

Toolkit condition:
```cpp
if (start_type == 1 && vtx_degree(vtx, graph) == 1 &&
    Eshower < 120 * units::MeV && (int)shower_segs.size() <= 3) {
    if (!topo && !traj && sg_length > 10 * units::cm)
        flag_bad_shower_1 = true;
}
if (sg_length > 80 * units::cm)
    flag_bad_shower_1 = true;
```

Matches prototype (`NeutrinoID_nue_tagger.h:3465-3486`) exactly. ✅

#### ✅ Logic fidelity — sub-check 2 (find_cont_muon_segment_nue calls)

The loop over `shower_segs` calls `find_cont_muon_segment_nue` at both endpoints of each segment (sv1, sv2), skipping `main_vertex`. Extension lengths are accumulated; then a length cut depending on `Eshower` and `n_connected` is applied. All thresholds match the prototype. ✅

#### ✅ Logic fidelity — sub-check 3

Near-end geometry test: checks for anti-parallel track segments within 5 cm of `other_vtx`, with bridge-segment search for the "neither end close" case. All thresholds and angular cuts match the prototype. ✅

#### ✅ Determinism — `IndexedSegmentSet shower_segs`

The prototype uses raw `map_seg_vtxs` / `map_vtx_segs`. The toolkit uses `shower->fill_sets(shower_vtxs, shower_segs, false)` which populates `IndexedSegmentSet` (sorted by stable integer index). Iteration over `shower_segs` is deterministic. **Improvement over prototype.**

---

## `cosmic_tagger()`

Ten independent cosmic-rejection flags (flag_1 through flag_10).

### Findings

#### ✅ Flag 1 — vertex outside fiducial volume

Toolkit:
```cpp
double tmp_dis = ray_length(Ray{mv_pt, mv_wcp_pt});
Point test_p = (tmp_dis > 5 * units::cm) ? mv_wcp_pt : mv_pt;
if (fiducial_utils && !fiducial_utils->inside_fiducial_volume(test_p))
    flag_cosmic_1 = true;
```

Prototype applies an extra `stm_tol_vec` shrinkage of 1.5 cm on all FV boundaries. The toolkit's `FiducialUtils::inside_fiducial_volume` does not yet support this per-call shrinkage. A `TODO` comment in the toolkit documents this limitation. The behavior difference is conservative (toolkit accepts slightly more events near the FV boundary). **Documented as known limitation; not a logic error in the porting.**

#### ✅ Flags 2–5 — single muon / long muon tests

The dominant muon/long-muon candidate selection loop, valid-tracks counting, connected-shower counting, and the dQ/dx front/end orientation logic all match the prototype. Angle conversions (`Angle()` in ROOT → `dir.angle() / M_PI * 180.0`) are correct. ✅

#### ✅ Flags 6–8 — stopped muon + Michel / secondary tests

The three-pass loop structure (muon+michel candidate search, muon_2nd search, valid_tracks counting) matches the prototype. The `valid_tracks` refinement block (angle-based decrements) is correctly gated on `(muon || long_muon) && michel_ele && muon_2nd`. ✅

#### NOTE — `Emi < 25 * units::cm` (prototype quirk, faithfully reproduced)

**Location:** toolkit L880 / prototype L398.

Both toolkit and prototype compare `Emi` (a shower energy in WireCell internal units) to `25 * units::cm`:

```cpp
double Emi = michel_energy;  // energy from shower->get_kine_best()
if (Emi < 25 * units::cm && ...)   // prototype has identical threshold
```

In WireCell, `units::MeV = 1` and `units::cm ≈ 10` (mm-based internal units), so this threshold equals `~250 MeV` rather than `25 MeV`. This appears to be a threshold accidentally written using length units instead of energy units in the prototype. **The toolkit faithfully reproduces this behavior.** Changing it would be a physics decision outside the scope of this porting review.

#### ✅ Flag 9 — global cluster-direction PCA analysis

The cluster loop, PCA/centroid direction computation, angle thresholds, and the `flagp_cosmic` → `flag_cosmic_9` decision all match the prototype.

#### NOTE — `pca.center.y()` vs `pts.front().y` (deliberate improvement)

**Location:** toolkit L1082-1086 / prototype L612-615.

For small clusters (length ≤ 3 cm), the prototype checks:
```cpp
if (it->first->get_point_cloud()->get_cloud().pts.front().y > 50*units::cm)
```
using the **first** point in the point cloud (an implementation-order artifact).

The toolkit uses:
```cpp
const auto& pca = cl->get_pca();
if (pca.center.y() > 50 * units::cm)
```
using the **PCA centroid**, which is a more stable and representative measure of where the cluster sits in y.

**This is a deliberate improvement.** Both implementations check whether the small cluster is in the upper part of the detector; the centroid is more robust than an arbitrary first point. The logic intent is preserved.

#### ✅ Flag 10 — front-face vertex / beam-direction check

Per-vertex loop, inside-FV check, z < 15 cm, segment angle to beam direction, and the `flag_cosmic_10_save` accumulation all match the prototype. ✅

---

## `calculate_num_daughter_showers()`

### Findings

#### ✅ Logic fidelity

BFS starting from `(vertex, segment)`, counting segments based on shower-flag and `flag_count_shower` parameter. Matches prototype `NeutrinoID_track_shower.h:688-722` exactly. ✅

#### ⚠ Determinism — raw pointer sets

Both `used_vertices` and `used_segments` are `std::set<VertexPtr>` / `std::set<SegmentPtr>` (pointer-address ordered). However, since the function only accumulates a **count and sum** (commutative operations), the BFS traversal order does not affect the result. For a tree-like graph, each segment/vertex is visited exactly once regardless of order.

**Action:** No code change required for correctness. Upgrading to `IndexedVertexSet` / `IndexedSegmentSet` would be consistent with the toolkit's determinism goal but is not urgent.

---

## `calculate_num_daughter_tracks()`

### Findings

#### ✅ Logic fidelity

Same BFS structure as `calculate_num_daughter_showers`, but counts non-shower segments above `length_cut`. Matches prototype `NeutrinoID_track_shower.h:724-763`. ✅

#### ⚠ Determinism

Same raw pointer set issue as `calculate_num_daughter_showers`; same reasoning applies. No impact on correctness.

---

## `find_cont_muon_segment_nue()`

### Findings

#### ✅ Logic fidelity

Iterates adjacent segments at `vtx`, computes opening angle (180° − `dir1.angle(dir2)`) using both 15 cm and 30 cm direction vectors for long segments, checks `angle_ok` and `dqdx_ok`, selects the best candidate by maximum projected length `length * cos(angle)`. Matches prototype `NeutrinoID_track_shower.h:2372-2441`. ✅

#### ✅ Efficiency — `dir3` precomputed

The toolkit computes:
```cpp
WireCell::Vector dir3 = (sg_length > 30 * units::cm)
                            ? segment_cal_dir_3vector(sg, vtx_pt, 30 * units::cm)
                            : dir1;
```
once outside the loop, rather than recomputing it per candidate as the prototype does. Minor efficiency improvement. ✅

---

## Summary Table

| Check | Result | Action |
|---|---|---|
| `shower_energy()` redundant helper | ⚠ Dead code | **Fixed** (6 sites) |
| `bad_reconstruction` sub-checks 1/2/3 | ✅ Equivalent | — |
| `IndexedSegmentSet` in `bad_reconstruction` | ✅ Improved | — |
| Flags 1–10 logic fidelity | ✅ Equivalent | — |
| `segs_at_vtx` lambda (deterministic BGL iteration) | ✅ Improved vs prototype | — |
| `pca.center.y()` vs `pts.front().y` in flag 9 | ✅ Deliberate improvement | — |
| `Emi < 25 * units::cm` threshold | ℹ Faithful prototype quirk | — |
| FV stm tolerance (flag 1) | ℹ Known limitation (TODO) | — |
| `calculate_num_daughter_showers/tracks` raw pointer sets | ⚠ Non-deterministic traversal order | No correctness impact; optional cleanup |
| `find_cont_muon_segment_nue` fidelity | ✅ Equivalent | — |
| `find_cont_muon_segment_nue` dir3 precompute | ✅ Minor efficiency improvement | — |

---

## Changes Made

**File:** `clus/src/NeutrinoTaggerCosmic.cxx`

1. Removed the `shower_energy()` static helper function (L85-88 pre-fix) — it was dead code duplicating `PRShower::get_kine_best()`.
2. Replaced all 6 call sites with `shower->get_kine_best()` or `michel_ele->get_kine_best()`.
