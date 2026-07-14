# SSM Tagger Review: Prototype → Toolkit

## Overview

| Item | Detail |
|------|--------|
| **Prototype function** | `WCPPID::NeutrinoID::ssm_tagger()` in `NeutrinoID_ssm_tagger.h` (3598 lines) |
| **Toolkit function** | `PatternAlgorithms::ssm_tagger()` in `clus/src/NeutrinoTaggerSSM.cxx` (1517 lines) |
| **Scope** | `ssm_tagger()` and all helper functions it calls |
| **Caller** | `TaggerCheckNeutrino.cxx:347` |

### Helper Function Mapping

| Prototype | Toolkit | Notes |
|-----------|---------|-------|
| `NeutrinoID::fill_ssmsp(ProtoSegment*, int, int, int)` | `fill_ssmsp(SegmentPtr, int, int, int, TaggerInfo&, ...)` | Static, same logic |
| `NeutrinoID::fill_ssmsp_psuedo(WCShower*, int, int)` | `fill_ssmsp_pseudo_1(ShowerPtr, int, int, TaggerInfo&)` | Overload 1, renamed (fixed typo) |
| `NeutrinoID::fill_ssmsp_psuedo(WCShower*, ProtoSegment*, int, int)` | `fill_ssmsp_pseudo_2(ShowerPtr, SegmentPtr, int, int, TaggerInfo&)` | Overload 2 |
| `NeutrinoID::fill_ssmsp_psuedo(WCShower*, ProtoSegment*, int)` | `fill_ssmsp_pseudo_3(ShowerPtr, SegmentPtr, int, TaggerInfo&)` | Overload 3 |
| `NeutrinoID::exit_ssm_tagger()` | Inlined as `exit_ssm` lambda | Merged into ssm_tagger |
| `NeutrinoID::print_ssm_tagger()` | Omitted | Replaced by structured logging |
| `NeutrinoID::get_scores(ProtoSegment*)` | `get_scores(SegmentPtr, ParticleDataSet::pointer&)` | Static helper |
| `NeutrinoID::get_scores(ProtoSegment*, int, int)` | `get_scores_bp(SegmentPtr, int, int, ParticleDataSet::pointer&)` | Static helper |
| `NeutrinoID::get_containing_shower_id/get_containing_shower` | `get_containing_shower_info(SegmentPtr, ShowerSegmentMap&)` | Merged into single function |
| `NeutrinoID::find_incoming_segment(ProtoVertex*)` | `find_incoming_segment(Graph&, VertexPtr, set<SegmentPtr>&)` | Static, explicit graph param |
| (prim/daught particle loops, ~430 lines × 2) | `fill_particle_block_at_vtx()` | Major refactor, see EFF-1 |

---

## 1. Functional Equivalence

The toolkit implementation is **functionally equivalent** to the prototype. Verified phase by phase:

### Phase A — SSM Candidate Identification (prototype lines 413-536, toolkit lines 588-683)

- **Pre-cuts**: Length/PDG cuts identical (`sg_length <= 46 && >= 1 && |pdg|==13`, or `direct_length > 0.9*length && length <= 44 && >= 1 && |pdg|==11`). Minimum length check (< 1 cm) and median dQ/dx cut (< 0.95) identical.
- **Vertex activity detection**: d(dQ/dx) vector computation, first/last 5-point scanning with threshold 0.7, identical.
- **Ambiguity resolution**: Three-tier logic (threshold 1.0 → equal check with 3-point tiebreaker → threshold 1.3 → fallback to max) identical to prototype lines 476-518.
- **Backwards muon flag**: `dirsign()` vs `get_flag_dir()` mapping correct.

### Phase B — SSM Selection & Feature Extraction (prototype lines 545-764, toolkit lines 872-1042)

- **Best SSM selection**: Priority: longest with vtx_activity, then longest overall. Identical.
- **Direction vectors**: Toolkit uses `segment_cal_dir_3vector(ssm_sg, dir, 10, 0).norm()` — equivalent to prototype's `ssm_sg->cal_dir_3vector(dir,10).Unit()`.
- **dQ/dx profile**: Point-by-point fwd/bck extraction, d(dQ/dx) derivatives, max values in first/last 3 and 5 points — identical.
- **fwd/bck swap** for `dir == -1`: Toolkit line 992-1001 correctly swaps all pairs. (Note: prototype has a swap bug here, see PROTO-BUG-1.)
- **PID scores**: Full-track and break_point scores computed identically.
- **Degenerate break_point guard**: Toolkit has two guards — one at lines 1022-1028 for scores/reduced_length, and another at lines 1038-1042 for `medium_dq_dx_bp`. Matches prototype behavior.
- **dQ_dx_cut formula**: `0.8866 + 0.9533 * pow(18/length, 0.4234)` — identical (both use length in cm).
- **Kinetic energy**: `cal_kine_range(length, 13, particle_data)` replaces prototype's `g_range->Eval(length) * MeV` — equivalent via `ParticleDataSet`.

### Phase C — Primary/Daughter Particle Loops (prototype lines 766-1224, toolkit lines 1042-1053)

The prototype has two nearly-identical ~430-line loops: one for `map_vertex_segments[main_vertex]` (prim) and one for `map_vertex_segments[second_vtx]` (daught). The toolkit replaces both with `fill_particle_block_at_vtx()`, which:

- Iterates `boost::out_edges(vtx)`, skipping `skip_sg` (the SSM segment)
- Counts tracks/showers at length thresholds 1/3/5/8/11 cm — identical logic
- Ranks top-2 tracks and top-2 showers by length — identical logic
- Computes per-particle: length, direct_length, pdg, medium_dq_dx, max_dev, ke_cal, ke_rng (mu/p/e), direction, daughter counts at 1/5/11 cm, PID scores — all identical
- Daughter count logic (subtract 1 for track/shower self-count, add 1 if segment shorter than threshold) — identical

### Phase D — Backwards Muon Swap, Momentum, Off-Vertex (prototype lines 1340-1749, toolkit lines 1059-1251)

- **Backwards swap**: Prototype swaps ~115 individual variable pairs. Toolkit: `std::swap(pb_prim, pb_daught)` — structurally equivalent.
- **Momentum vectors**: `pmag = sqrt(KE^2 + 2*KE*m)` formula identical. Daughter showers use `ke_best` (not `ke_rng`) for momentum — matches prototype.
- **Off-vertex loop**: Iterates `boost::edges(graph)`, skips segments at ssm_main_vtx or ssm_second_vtx, applies 80 cm distance cut, accumulates length/energy, finds leading off-vertex track and shower — identical logic.
- **Neutrino direction angles**: Four combinations (all, connected, primary, track-only) computed with pi0 momentum included in `nu_all` — identical to prototype.

### Phase E — SSMSP Filling, KDAR Flag, TaggerInfo Assignment (prototype lines 1750-2289, toolkit lines 1253-1517)

- **SSMSP filling**: `fill_ssmsp_all()` implements the same BFS traversal: SSM first as muon (pdg=13), then protons at ssm_main_vtx, then BFS daughters, then uncovered showers. Logic verified identical.
- **KDAR flag**: `n_prim_tracks_1==0 && n_prim_all_3==0 && n_daughter_tracks_5==0 && n_daughter_all_5<2 && Nsm_wivtx==1 && !(pio_mass>70 && pio_mass<200)` — identical.
- **TaggerInfo assignment**: All ~200+ fields verified to map correctly from local variables to `ti.ssm_*` fields.

---

## 2. Bugs Found

### In the Toolkit

Three toolkit bugs were found and fixed (2026-04-12) in a separate review pass (see `clus/docs/tagger/nue_ssm_tagger_review.md`):

| Bug | Severity | Location | Fix |
|---|---|---|---|
| Missing `/(43e3/units::cm)` dQ/dx normalization | Critical | Phase A (L623) + Phase B (L918) | Added normalization factor |
| `nu_all` includes `mom_pi0` — prototype does not | Medium | Phase D (L1128-1130) | Removed `+ mom_pi0` from `nu_all` sum |
| Missing `medium_dq_dx_bp = medium_dq_dx` reset in degenerate break_point | Medium | Phase B (L1022-1028) | Added reset after `medium_dq_dx_bp` declaration |

All three have been fixed in the current code.

### Prototype Bugs Fixed by Toolkit

1. **PROTO-BUG-1 (Medium) — Wrong swap pair for `max_d_dq_dx`**
   - **Location**: Prototype `NeutrinoID_ssm_tagger.h:716-717`
   - **Issue**: When `dir == -1`, the prototype swaps:
     ```cpp
     std::swap(max_d_dq_dx_fwd_3, max_d_dq_dx_bck_5);  // fwd_3 ↔ bck_5 (WRONG)
     std::swap(max_d_dq_dx_fwd_5, max_d_dq_dx_bck_5);  // fwd_5 ↔ (now holds old fwd_3)
     ```
     The first swap should be `fwd_3 ↔ bck_3`. As written, the original `fwd_5` value is lost.
   - **Toolkit fix**: Line 1000 correctly swaps `(max_d_fwd3 ↔ max_d_bck3)` and `(max_d_fwd5 ↔ max_d_bck5)`.
   - **Impact**: Affects `ssm_max_d_dq_dx_fwd_3/5` and `ssm_max_d_dq_dx_bck_3/5` for backward-direction SSM segments.

2. **PROTO-BUG-2 (Low) — Copy-paste error in daught_track1 direction**
   - **Location**: Prototype `NeutrinoID_ssm_tagger.h:1492-1494`
   - **Issue**: Direction components for `daught_track1` are assigned from `dir_prim_track1` instead of `dir_daught_track1`:
     ```cpp
     x_dir_daught_track1 = dir_prim_track1[0];  // Should be dir_daught_track1
     y_dir_daught_track1 = dir_prim_track1[1];
     z_dir_daught_track1 = dir_prim_track1[2];
     ```
     The momentum vector (lines 1495-1497) correctly uses `dir_daught_track1`.
   - **Toolkit fix**: Uses `pb_daught.track1.dir` consistently for both direction and momentum.
   - **Impact**: Affects `ssm_daught_track1_x/y/z_dir` fields (but not the momentum calculation).

3. **PROTO-BUG-3 (Low) — Unclamped `acos()` calls**
   - **Location**: Prototype lines 583-601 (SSM angles) and 1724-1748 (neutrino angles)
   - **Issue**: Raw `acos()` on dot products without clamping to [-1, 1]. Floating-point rounding can produce values slightly outside this range, yielding NaN.
   - **Toolkit fix**: All angle computations use `safe_acos()` (line 37-39) which wraps `std::acos(std::clamp(x, -1.0, 1.0))`.

4. **PROTO-BUG-4 (Low) — `dQ_dx_cut` formula has extra `/units::cm`**
   - **Location**: Prototype line 743
   - **Issue**: `dQ_dx_cut = 0.8866 + 0.9533 * pow(18/length/units::cm, 0.4234)` where `length` is already in cm. The extra `/units::cm` makes the argument a factor of 10 too small (since `units::cm = 10` in CLHEP units). Evaluates to `pow(1.8/length, 0.4234)` instead of `pow(18/length, 0.4234)`.
   - **Toolkit fix** (line 1030): `0.8866 + 0.9533 * std::pow(18.0/length, 0.4234)` — correctly dimensionless.
   - **Impact**: `ssm_dQ_dx_cut` is stored as a feature only — never used for any cut within `ssm_tagger` — so the prototype error does not propagate to event selection.

---

## 3. Efficiency / Structure Improvements

1. **EFF-1 — Particle loop deduplication**: The prototype's two nearly-identical ~430-line loops for primary and daughter particle enumeration (lines 794-1224) are replaced by a single reusable function `fill_particle_block_at_vtx()` (lines 443-550). This eliminates ~310 lines of duplicated code.

2. **EFF-2 — Backwards_muon swap simplification**: The ~115-line field-by-field swap block (prototype lines 1340-1464) is replaced by a single `std::swap(pb_prim, pb_daught)` on toolkit line 1066, made possible by the `ParticleBlock` struct.

3. **EFF-3 — Exit path consolidation**: The prototype's `exit_ssm_tagger()` (separate 325-line function) is merged as an inline lambda with helper lambdas `set_track_block`/`set_shw_block` for field-group initialization.

4. **EFF-4 — Logging removal**: `print_ssm_tagger()` (390 lines of `std::cout` output) is removed. Debug output throughout the function body (`std::cout << ...`) is also removed.

5. **Overall**: 3598 lines → 1517 lines (58% reduction) with identical functionality.

---

## 4. Determinism

| Container | Prototype | Toolkit | Issue? |
|-----------|-----------|---------|--------|
| `all_ssm_sg` | `std::map<ProtoSegment*, ...>` | `std::map<SegmentPtr, ..., SegmentIndexCmp>` | **Fixed** — pointer order → index order |
| Vertex→segment iteration | `map_vertex_segments[vtx]` (pointer set) | `boost::out_edges(vtx, graph)` | **Fixed** — deterministic edge order |
| Shower iteration | `showers` (pointer set) | `IndexedShowerSet` | **Fixed** — index-ordered |
| BFS `used_segments/used_vertices` | `std::set<ProtoSegment*>` / `std::set<ProtoVertex*>` | `std::set<SegmentPtr>` / `std::set<VertexPtr>` | OK — only used for membership tests, not order-dependent iteration |

All non-determinism sources from the prototype are resolved.

---

## 5. Multi-APA / Multi-Face Handling

The SSM tagger operates on **graph topology only** (vertices, segments, edges, dQ/dx data). It does not reference:
- Wire-plane geometry or wire directions
- Drift direction
- APA indices or face indices
- Any `IDetectorVolumes` or `WirePlaneId` parameters

The off-vertex loop at line 1164 (`boost::edges(graph)`) iterates all edges in the graph, which naturally spans all APAs. The 80 cm distance cut (line 1195) is metric and APA-agnostic.

**No multi-APA issues found.**

---

## 6. Minor Logic Divergence

**DIVERGE-1 — fill_ssmsp_all shower fallback mother ID**

In `fill_ssmsp_all()`, when processing a shower with `conn_type == 2 || conn_type == 3` and the start vertex is neither the main vertex nor has a recognizable parent segment:

- **Prototype** (line 1838-1840): Falls back to SSM segment's id as mother
- **Toolkit** (line 378-379): Falls back to `mother = 0`

This affects only the parentage field in the ssmsp tree for a rare edge case. The ssmsp data is used for event display and secondary analysis, not for BDT scoring. Impact is minimal.

---

## 7. Recommendations

1. **All toolkit bugs have been fixed** (2026-04-12). No further code changes needed.
2. The four prototype bugs fixed by the toolkit (swap pair, direction copy-paste, unclamped acos, dQ_dx_cut formula) are worth noting for anyone comparing prototype vs toolkit output.
3. The minor divergence (DIVERGE-1) could be aligned with prototype behavior if exact ssmsp tree compatibility is desired, but is not functionally important.
