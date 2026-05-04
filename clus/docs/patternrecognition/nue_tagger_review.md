# NuE Tagger / SinglePhoton Tagger / SSM Tagger Review

**Date**: 2026-04-12
**Scope**: `nue_tagger`, `singlephoton_tagger`, `ssm_tagger`, and all sub-functions
**Files reviewed**:
- Toolkit: `clus/src/NeutrinoTaggerNuE.cxx` (4402 lines), `clus/src/NeutrinoTaggerSinglePhoton.cxx` (2509 lines), `clus/src/NeutrinoTaggerSSM.cxx` (1517 lines)
- Prototype: `pid/src/NeutrinoID_nue_tagger.h`, `pid/src/NeutrinoID_nue_functions.h`, `pid/src/NeutrinoID_singlephoton_tagger.h`, `pid/src/NeutrinoID_ssm_tagger.h`, `pid/src/NeutrinoID_nue_bdts.h`

**Review criteria**: (1) functional equivalence with prototype, (2) bugs, (3) efficiency, (4) determinism, (5) multi-APA support

---

## Summary of Issues

| ID | Severity | File | Function | Line | Description | Status |
|----|----------|------|----------|------|-------------|--------|
| B1 | Medium | NeutrinoTaggerNuE.cxx | `angular_cut` | 322 | Missing `-1.5cm` fiducial tolerance vector | FIXED |
| ~~B2~~ | ~~Medium~~ | ~~NeutrinoTaggerNuE.cxx~~ | ~~`bad_reconstruction_1`~~ | ~~2173~~ | ~~Missing `kine_charge` fallback~~ | NOT A BUG — `get_kine_best()` already falls back internally |
| ~~B3~~ | ~~Medium~~ | ~~NeutrinoTaggerNuE.cxx~~ | ~~`bad_reconstruction_3`~~ | ~~2868~~ | ~~Missing `kine_charge` fallback~~ | NOT A BUG — same reason |
| ~~B4~~ | ~~Medium~~ | ~~NeutrinoTaggerNuE.cxx~~ | ~~`bad_reconstruction_2`~~ | ~~3084~~ | ~~Missing `kine_charge` fallback~~ | NOT A BUG — same reason |
| B5 | **Critical** | NeutrinoTaggerSSM.cxx | `ssm_tagger` | 623 | Missing `/(43e3/units::cm)` dQ/dx normalization | FIXED |
| B6 | Medium | NeutrinoTaggerSSM.cxx | `ssm_tagger` | 1128-1130 | `nu_all` includes `mom_pi0` — prototype does not | FIXED |
| B7 | Medium | NeutrinoTaggerSSM.cxx | `ssm_tagger` | 1022-1028 | Missing `medium_dq_dx_bp = medium_dq_dx` reset in degenerate break_point | FIXED |
| B8 | Medium | NeutrinoTaggerSinglePhoton.cxx | `bad_reconstruction_2_sp` | 729-739 | br3_3 `angle>105` check and TaggerInfo fills incorrectly guarded by `dir1.magnitude()>10cm` | FIXED |
| B9 | Low | NeutrinoTaggerSinglePhoton.cxx | `high_energy_overlapping_sp` | 1626 | Incorrect `flag_all_showers=false` when `dir2.magnitude()==0` | FIXED |
| B10 | Medium | NeutrinoTaggerSinglePhoton.cxx | `singlephoton_tagger` | 2434-2463 | No SCE correction on vertex/shower positions for distance calc | OPEN (may be intentional) |
| M1 | Info | NeutrinoTaggerSSM.cxx | `ssm_tagger` | 580-583 | Hardcoded MicroBooNE beam/target/absorber directions | NOTED |
| M2 | Info | — | `cal_bdts_xgboost` | — | Not ported to toolkit | NOTED |

> **Note on B2/B3/B4:** The toolkit's `PRShower::get_kine_best()` already contains the fallback: `if (kenergy_best != 0) return kenergy_best; else return kenergy_charge;`. The explicit fallback in the prototype is redundant when ported to the toolkit. No code change needed.

---

## NeutrinoTaggerNuE.cxx — Detailed Findings

### `low_energy_michel` (line 125 vs prototype line 266)
- **Equivalence**: OK. All thresholds (25cm, 18cm, 100MeV, 0.7) match. `n_3seg` counting via `boost::out_edges` is equivalent to prototype's `map_vtx_segs[vtx].size()`. Unused prototype variables (`E_range`, `n_segs`) correctly omitted.
- **Bugs**: None
- **Efficiency**: OK
- **Determinism**: OK — `IndexedVertexSet` iteration is index-ordered
- **Multi-APA**: OK

### `stem_length` (line 190 vs prototype line 379)
- **Equivalence**: OK. All thresholds (500MeV, 50cm, 55cm) and exception logic identical.
- **Bugs**: None
- **Efficiency**: OK
- **Determinism**: OK
- **Multi-APA**: OK

### `angular_cut` (line 233 vs prototype line 436)
- **Equivalence**: DIFFERENCE — fiducial tolerance omitted
- **Bug B1** (line 322): The prototype passes `stm_tol_vec = {-1.5*units::cm, ...}` (5 elements) to `inside_fiducial_volume` at line 504/511, which shrinks the fiducial volume by 1.5 cm on each face. The toolkit calls `inside_fiducial_volume(vtx_fit_pt(vtx1))` with no tolerance. This makes the toolkit's `flag_main_outside` check less strict — vertices within 1.5cm of boundaries that prototype flags as outside will be considered inside.
- **Fix**: Add `{-1.5*units::cm, -1.5*units::cm, -1.5*units::cm, -1.5*units::cm, -1.5*units::cm}` as second arg to `inside_fiducial_volume` at line 322.
- **Efficiency**: OK
- **Determinism**: OK
- **Multi-APA**: OK

### `compare_muon_energy` (line 363 vs prototype line 625)
- **Equivalence**: OK. All thresholds and logic identical. Muon range function call equivalent.
- **Bugs**: None
- **Efficiency**: OK
- **Determinism**: OK
- **Multi-APA**: OK

### `stem_direction` (line 445 vs prototype line 1)
- **Equivalence**: OK. All angle/ratio computations, thresholds, and branching logic match.
- **Bugs**: None
- **Efficiency**: OK
- **Determinism**: OK
- **Multi-APA**: OK

### `pi0_identification` (line 549 vs prototype line 2643)
- **Equivalence**: OK. Mass conditions and veto logic match.
- **Bugs**: None
- **Efficiency**: IMPROVEMENT — precomputed `cluster_acc_length` avoids O(V*S) per-vertex recomputation.
- **Determinism**: OK — `IndexedVertexSet` used
- **Multi-APA**: OK

### `single_shower` (line 694 vs prototype line 320)
- **Equivalence**: OK. All conditionals and thresholds match.
- **Bugs**: None
- **Efficiency**: OK
- **Determinism**: OK
- **Multi-APA**: OK — direction vectors `dir_beam(0,0,1)`, `dir_drift(1,0,0)`, `dir_vertical(0,1,0)` are coordinate conventions

### `multiple_showers` (line 830 vs prototype line 81)
- **Equivalence**: OK. All thresholds and logic match.
- **Bugs**: None
- **Efficiency**: OK
- **Determinism**: OK — iterates `ctx.showers` (ordered container)
- **Multi-APA**: OK

### `other_showers` (line 964 vs prototype line 209)
- **Equivalence**: OK. Benign omission of unused variables `flag_direct_max_pi0` and `n_indirect_showers` (only used in debug prints in prototype).
- **Bugs**: None
- **Efficiency**: OK
- **Determinism**: OK
- **Multi-APA**: OK

### `vertex_inside_shower` (line 1133 vs prototype line 407)
- **Equivalence**: OK. All conditionals, thresholds, and block 1/2 structure match.
- **Bugs**: None
- **Efficiency**: OK
- **Determinism**: OK
- **Multi-APA**: OK — `dir_drift(1,0,0)` and `dir_beam(0,0,1)` are physics constants

### `broken_muon_id` (line 1369 vs prototype line 1010)
- **Equivalence**: OK. Angular cuts, distance thresholds, and muon walk algorithm match. `add_length` correctly dropped (only used in debug prints).
- **Bugs**: None
- **Efficiency**: OK
- **Determinism**: `std::set<SegmentPtr>` for `muon_segments` uses pointer ordering, but only for membership checks and accumulation — no order-dependent results. Acceptable.
- **Multi-APA**: OK

### `mip_quality` (line 1592 vs prototype line 1669)
- **Equivalence**: OK. Overlap check and split check match.
- **Bugs**: None
- **Efficiency**: OK
- **Determinism**: OK
- **Multi-APA**: OK — 2D distance calls properly take `ctx.apa` and `ctx.face`

### `high_energy_overlapping` (line 1782 vs prototype line 2284)
- **Equivalence**: OK. All flag_overlap1/2 cuts match.
- **Bugs**: None. The `medium_dQ_dx` endpoint logic uses geometric proximity instead of index comparison — semantically equivalent, negligible risk.
- **Efficiency**: OK
- **Determinism**: OK
- **Multi-APA**: OK

### `low_energy_overlapping` (line 1973 vs prototype line 2485)
- **Equivalence**: OK. All three sub-flags match. Prototype's never-set `flag_overlap_4/5` correctly omitted.
- **Bugs**: None
- **Efficiency**: OK
- **Determinism**: OK
- **Multi-APA**: OK

### `bad_reconstruction_1` (line 2167 vs prototype line 3854)
- **Equivalence**: OK — `get_kine_best()` already contains the `kine_charge` fallback internally in the toolkit's `PRShower` class, so the explicit ternary in the prototype is redundant.
- **Bugs**: None
- **Efficiency**: OK
- **Determinism**: OK
- **Multi-APA**: OK

### `shower_to_wall` (line 2295 vs prototype line 1219)
- **Equivalence**: OK. All four sub-flags (flag_bad1-4) and their cuts match.
- **Bugs**: None
- **Efficiency**: OK
- **Determinism**: OK — `graph_nodes` returns deterministic vertex descriptors
- **Multi-APA**: OK — uses `fiducial_utils->inside_fiducial_volume()` with tolerance vec

### `single_shower_pio_tagger` (line 2530 vs prototype line 2756)
- **Equivalence**: OK + IMPROVEMENT. Prototype line 2882 has `Eshower < 800*units::cm` (wrong units); toolkit line 2655 correctly uses `units::MeV`. This is a prototype bug fix.
- **Bugs**: None (improvement)
- **Efficiency**: OK
- **Determinism**: OK — `IndexedSegmentSet`/`IndexedVertexSet` used
- **Multi-APA**: OK

### `gap_identification` (line 2692 vs prototype line 1438)
- **Equivalence**: OK + IMPROVEMENT. Prototype omits `flag_prolong_w` in the front path (line 1515-1516) but includes it in the back path. Toolkit consistently includes it in both paths (line 2780). This fixes a prototype bug.
- **Bugs**: None (improvement)
- **Efficiency**: OK
- **Determinism**: OK
- **Multi-APA**: OK — uses `ctx.dv->contained_by()` for APA/face lookup

### `bad_reconstruction_3` (line 2862 vs prototype line 2907)
- **Equivalence**: OK — `get_kine_best()` already contains the `kine_charge` fallback internally.
- **Bugs**: None
- **Efficiency**: OK
- **Determinism**: OK
- **Multi-APA**: OK

### `bad_reconstruction_2` (line 3076 vs prototype line 3145)
- **Equivalence**: OK — `get_kine_best()` already contains the `kine_charge` fallback internally.
- **Bugs**: None
- **Efficiency**: OK
- **Determinism**: OK
- **Multi-APA**: OK

### `track_overclustering` (line 3387 vs prototype line 547)
- **Equivalence**: OK. All five sub-checks (tro_1 through tro_5) match.
- **Note**: Line 3562 faithfully reproduces prototype typo `stem_length_1 / units::MeV` (should be `/units::cm`). Since this is a BDT input trained on prototype outputs, changing it would break the model.
- **Bugs**: None
- **Efficiency**: OK
- **Determinism**: `std::set<SegmentPtr>` for `muon_segs` uses pointer ordering but only for membership/accumulation — acceptable.
- **Multi-APA**: OK

### `mip_identification` (line 3807 vs prototype line 1823)
- **Equivalence**: OK. All dQ/dx cuts, threshold logic, early_mip classification, and strong_check overrides match.
- **Bugs**: None
- **Efficiency**: OK
- **Determinism**: OK
- **Multi-APA**: OK

### `nue_tagger` (line 4221 vs prototype line 2)
- **Equivalence**: OK. Calling sequence matches prototype exactly. All sub-taggers called in same order with same flag logic.
- **Bugs**: None
- **Efficiency**: OK
- **Determinism**: OK
- **Multi-APA**: OK — `NuEContext` carries `apa` and `face`

---

## NeutrinoTaggerSinglePhoton.cxx — Detailed Findings

### `bad_reconstruction_sp` (line 136 vs prototype line 3766)
- **Equivalence**: OK. All three sub-checks (br1_1, br1_2, br1_3) match.
- **Bugs**: None
- **Efficiency**: OK — lambdas consolidate repeated if-else chains
- **Determinism**: OK — `IndexedSegmentSet` used
- **Multi-APA**: OK

### `bad_reconstruction_1_sp` (line 482 vs prototype line 4170)
- **Equivalence**: OK. All angle computations and PCA calls match.
- **Bugs**: None
- **Efficiency**: OK
- **Determinism**: OK
- **Multi-APA**: OK

### `bad_reconstruction_2_sp` (line 612 vs prototype line 3461)
- **Equivalence**: DIFFERENCE — br3_3 scoping error
- **Bug B8** (lines 729-739): In the prototype (line 3588-3597), the `angle > 90` accumulation check (`acc_length += length`) is guarded by `dir1.Mag() > 10*units::cm`, but `angle > 105 && len1 > 15cm` check at line 3590 is NOT guarded by `dir1.Mag()`. The toolkit places BOTH checks inside the `if (dir1.magnitude() > 10*units::cm)` guard at line 729. Additionally, the TaggerInfo vector fills (lines 735-739) are inside this guard, while the prototype fills for ALL segments regardless of `dir1.Mag()`.
- **Impact**: (1) Segments with short endpoint-to-endpoint vectors but long track lengths (>15cm, possible for curvy segments) will not be flagged. (2) BDT input vectors will have fewer entries than prototype, potentially breaking inference.
- **Fix**: Move `angle > 105` check and vector pushes outside the `dir1.magnitude()` guard. Keep the `acc_length` accumulation and `angle > 150` check inside it.

### `bad_reconstruction_3_sp` (line 905 vs prototype line 3223)
- **Equivalence**: OK. All thresholds and sub-checks match.
- **Bugs**: None
- **Efficiency**: OK
- **Determinism**: OK
- **Multi-APA**: OK

### `mip_identification_sp` (line 1114 vs prototype line 2102)
- **Equivalence**: OK. All dQ/dx logic, thresholds, and fill values match.
- **Bugs**: None
- **Efficiency**: OK
- **Determinism**: OK
- **Multi-APA**: OK

### `high_energy_overlapping_sp` (line 1579 vs prototype line 2596)
- **Equivalence**: DIFFERENCE — `flag_all_showers` handling
- **Bug B9** (line 1626): When `dir2.magnitude() == 0` for an electron/weak-muon segment, the prototype's `continue` does NOT modify `flag_all_showers`. The toolkit sets `flag_all_showers = false` before continuing. This makes the overlap cut slightly more conservative.
- **Fix**: Remove `flag_all_showers = false;` from the zero-magnitude guard.
- **Efficiency**: OK
- **Determinism**: OK
- **Multi-APA**: OK

### `low_energy_overlapping_sp` (line 1766 vs prototype line 2797)
- **Equivalence**: OK. All three sub-checks match.
- **Bugs**: None
- **Efficiency**: OK
- **Determinism**: OK
- **Multi-APA**: OK

### `pi0_identification_sp` (line 1957 vs prototype line 2955)
- **Equivalence**: OK. Both pio_1/pio_2 paths match.
- **Bugs**: None
- **Efficiency**: IMPROVEMENT — precomputed cluster lengths
- **Determinism**: OK
- **Multi-APA**: OK

### `low_energy_michel_sp` (line 2094 vs prototype line 545)
- **Equivalence**: OK. Both length and charge criteria match. Unused `E_range` correctly omitted.
- **Bugs**: None
- **Efficiency**: OK
- **Determinism**: OK
- **Multi-APA**: OK

### `singlephoton_tagger` (line 2167 vs prototype line 2)
- **Equivalence**: DIFFERENCE — no SCE correction
- **Bug B10** (lines 2434-2463): Prototype applies `func_pos_SCE_correction` to `nu_vtx`, `trk_vtx`, and `shw_vtx_pt` before computing distances. Toolkit uses raw positions (`ctx.geom_helper` is set to `nullptr`). This affects `shw_vtx_dis` and `max_shw_dis` values and the final `max_shw_dis < 2` cut.
- **Note**: This may be intentional if SCE is handled differently in the toolkit pipeline, but should be verified.
- **Efficiency**: OK
- **Determinism**: `std::set<ShowerPtr>` for `good_showers`/`ok_showers` uses pointer ordering but only for membership checks — acceptable.
- **Multi-APA**: OK — derives `apa` and `face` from vertex position

---

## NeutrinoTaggerSSM.cxx — Detailed Findings

### `ssm_tagger` (line 562 vs prototype `NeutrinoID_ssm_tagger.h`)

- **Equivalence**: MULTIPLE DIFFERENCES

#### Bug B5 (CRITICAL) — Missing dQ/dx normalization (line 623-626, 918-925, 965-989)

The prototype normalizes all dQ/dx values by `/(43e3/units::cm)` before computing differences and comparing thresholds (prototype line 438):
```cpp
double last_dq_dx = vec_dQ.at(0)/vec_dx.at(0)/(43e3/units::cm);
```

The toolkit uses raw dQ/(dx/cm) without the normalization factor:
```cpp
double last = fits[0].dQ / (fits[0].dx/units::cm);  // line 623
```

**Impact**: The threshold `0.7` for vertex activity detection was tuned for normalized (dimensionless ~1.0) values. Raw dQ/dx values are ~43000x larger, so:
- Phase A: vertex activity will fire on nearly every segment (false positives)
- Phase B: break_point detection will be incorrect
- All `ssm_dq_dx_fwd_*`, `ssm_dq_dx_bck_*`, `ssm_d_dq_dx_*`, `ssm_max_dq_dx_*`, `ssm_max_d_dq_dx_*` BDT inputs will be on wrong scale

**Fix**: Divide dQ/dx by `(43e3/units::cm)` everywhere it is computed. E.g. line 623:
```cpp
double last = fits[0].dQ / (fits[0].dx/units::cm) / (43e3/units::cm);
```

#### Bug B6 — `nu_all` includes `mom_pi0` (line 1128-1130)

Prototype line 1722: `nu_dir = mom + mom_prim_track1 + ... + mom_offvtx_track1 + mom_offvtx_shw1;` — does NOT include `mom_pi0`. Toolkit line 1128-1130 includes `mom_pi0` in the sum.

**Impact**: `ssm_nu_angle_z/target/absorber/vertical` values will differ when pi0 is present.

**Fix**: Remove `+ mom_pi0` from line 1130.

#### Bug B7 — Missing `medium_dq_dx_bp` reset (line 1022-1028)

When break_point is degenerate (covers all points), the prototype resets `medium_dq_dx_bp = medium_dq_dx` (line 750-756). The toolkit does not.

**Fix**: Add `medium_dq_dx_bp = medium_dq_dx;` inside the degenerate break_point block.

#### Intentional Fixes (toolkit improvements over prototype)

1. **`dQ_dx_cut` formula** (line 1030): Prototype has spurious `units::cm` in denominator (`18/length/units::cm`). Toolkit correctly uses unitless `18.0/length`.

2. **`max_d_dq_dx` swap** (line 1000): Prototype has copy-paste bug swapping `fwd_3<->bck_5` instead of `fwd_3<->bck_3`. Toolkit correctly swaps `fwd_3<->bck_3` and `fwd_5<->bck_5`. **Warning**: If BDT models were trained on prototype's buggy outputs, the toolkit's correct values will cause a model-data mismatch.

#### Determinism
- IMPROVEMENT: Uses `SegmentIndexCmp` for ordered map (line 593) instead of raw pointer comparison. This is correct and deterministic.

#### Multi-APA
- Hardcoded MicroBooNE beam/target/absorber directions (lines 580-583). Same as prototype — not a regression, but should be configurable for other detectors.

---

## `cal_bdts_xgboost` — Not Ported

The prototype's `NeutrinoID_nue_bdts.h` contains `cal_bdts_xgboost()` (line 3) which uses TMVA readers for nue BDT scoring. No corresponding function exists in the toolkit's `clus/src/` directory. This function is not called from the current tagger code and remains unported.

---

## Improvements in Toolkit Over Prototype

1. **pi0_identification / pi0_identification_sp**: Precomputed `cluster_acc_length` avoids O(V*S) per-vertex recomputation
2. **single_shower_pio_tagger**: Fixed `units::cm` -> `units::MeV` (prototype bug at line 2882)
3. **gap_identification**: Added missing `flag_prolong_w` in front path
4. **ssm_tagger**: Fixed `dQ_dx_cut` formula units and `max_d_dq_dx` swap bug
5. **Determinism**: `IndexedSegmentSet`/`IndexedVertexSet`/`SegmentIndexCmp` used throughout instead of pointer-keyed containers
6. **Multi-APA**: `NuEContext` carries `apa`/`face`, `gap_identification` uses `dv->contained_by()` for APA/face lookup
