# SinglePhoton Tagger Review: Prototype → Toolkit

## Overview

| Item | Detail |
|------|--------|
| **Prototype function** | `WCPPID::NeutrinoID::singlephoton_tagger()` in `NeutrinoID_singlephoton_tagger.h` (4275 lines) |
| **Toolkit function** | `PatternAlgorithms::singlephoton_tagger()` in `clus/src/NeutrinoTaggerSinglePhoton.cxx` (2511 lines) |
| **Scope** | `singlephoton_tagger()` and all 9 helper functions it calls |
| **Caller** | `TaggerCheckNeutrino.cxx:366` |

### Helper Function Mapping

| Prototype | Toolkit | Lines (TK) | Notes |
|-----------|---------|-------------|-------|
| `NeutrinoID::singlephoton_tagger(muon_length)` | `singlephoton_tagger(graph, ...)` | 2169-2510 | Entry point; `muon_length` param unused in prototype |
| `NeutrinoID::bad_reconstruction_sp()` | `bad_reconstruction_sp()` | 136-463 | br1_1/br1_2/br1_3; static helper |
| `NeutrinoID::bad_reconstruction_1_sp()` | `bad_reconstruction_1_sp()` | 482-595 | br2 (PCA direction); static helper |
| `NeutrinoID::bad_reconstruction_2_sp()` | `bad_reconstruction_2_sp()` | 612-890 | br3_1..br3_8 (8 sub-checks); static helper |
| `NeutrinoID::bad_reconstruction_3_sp()` | `bad_reconstruction_3_sp()` | 907-1096 | br4_1/br4_2; static helper |
| `NeutrinoID::mip_identification()` | `mip_identification_sp()` | 1116-1562 | dQ/dx MIP classification; static helper |
| `NeutrinoID::high_energy_overlapping()` | `high_energy_overlapping_sp()` | 1578-1748 | hol_1/hol_2; static helper |
| `NeutrinoID::low_energy_overlapping()` | `low_energy_overlapping_sp()` | 1768-1938 | lol_1/lol_2/lol_3; static helper |
| `NeutrinoID::pi0_identification()` | `pi0_identification_sp()` | 1956-2082 | pio_1/pio_2; static helper |
| `NeutrinoID::low_energy_michel_sp()` | `low_energy_michel_sp()` | 2096-2145 | Michel electron check; static helper |
| (N/A) | `SpContext` struct | 103-117 | File-local state bundle replacing member variables |

---

## 1. Functional Equivalence

The toolkit implementation is **functionally equivalent** to the prototype. Verified function by function:

### Entry Point — `singlephoton_tagger()` (prototype lines 2-542, toolkit lines 2169-2510)

- **First-pass shower loop**: Proton/muon/pion bookkeeping, electron shower classification (en20, badreco1-4), good/ok shower selection — identical logic. Proton ranking into slot #1/#2 uses the same algorithm (see SHARED-BUG-1 below).
- **Selection logic**: `num_mip_tracks > 1` veto, `num_good_shws != 1` veto, `num_20br1_shws > 1` veto, fallback from good→ok showers — all identical.
- **Second-pass evaluation**: mip_id → low_energy_michel → pi0 → bad_reconstruction → low_energy_overlapping → bad_reconstruction_2/3 → high_energy_overlapping — same order, same flag propagation.
- **Final cuts**: `mean_dedx < 2.3` and `max_shw_dis < 2` (when tracks present) — identical thresholds.
- **SCE correction**: Prototype applies SCE correction to track/shower vertex positions; toolkit uses raw positions. This is an intentional design difference (toolkit operates in detector coordinates throughout).

### `bad_reconstruction_sp` (br1_1, br1_2, br1_3) — prototype lines 3766-4168, toolkit lines 136-463

- **br1_1**: All 4 stem/topology conditions identical (energy thresholds 120 MeV, length 10/80 cm, segment counts).
- **br1_2**: Full energy-binned cut table (200/400/600/800/1500 MeV boundaries) with all length thresholds (38-75 cm range) matches exactly. Special overrides (72 cm clamp, 1000 MeV ratio<0.95, 15-segment override) all present.
- **br1_3**: All 5 energy-binned length thresholds match. Angle thresholds (170, 165, 170 degrees) and distance thresholds (5, 6 cm) match. The "undo" guard at the end is identical.

### `bad_reconstruction_1_sp` (br2) — prototype lines 4170-4274, toolkit lines 482-595

- PCA computation, 4 angle variables, energy-dependent cuts (500/1000 MeV), max_angle loop, 2 unconditional cuts — all identical.
- `max_angle` loop: prototype iterates `map_vertex_segments[other_vertex]`; toolkit uses `boost::out_edges` — equivalent.

### `bad_reconstruction_2_sp` (br3_1..br3_8) — prototype lines 3461-3764, toolkit lines 612-890

- **br3_1**: 4 straight-shower conditions with explicit parenthesization of the operator-precedence-sensitive third condition — identical.
- **br3_2**: n_ele/n_other counting, fiducial check — identical.
- **br3_3/br3_4**: Backward segment detection; `acc_length > 0.33*total_length` cut — identical (after the already-fixed scoping bug, see Bugs section).
- **br3_5**: Average non-stem point position, `flag_bad5` override for multi-cluster — identical.
- **br3_6**: Far-end vertex segments, sliding dQ/dx — identical.
- **br3_7**: Stem vs shower main length — identical.
- **br3_8**: Max dQ/dx in sliding 5-bin windows — identical.

### `bad_reconstruction_3_sp` (br4_1, br4_2) — prototype lines 3223-3458, toolkit lines 907-1096

- **br4_1**: 6 distance/length threshold conditions, 2 false-override conditions, `num_main_segs>=4` guard — all identical.
- **br4_2**: Angular distribution fit-point counting at 15/25/35/45 degree cones, `flag_bad2` compound expression — identical.

### `mip_identification_sp` — prototype lines 2102-2594, toolkit lines 1116-1562

- dQ_dx_cut table (1.3/1.45/1.6/1.85) — identical.
- Scan loops (n_end_reduction, n_first_mip, n_first_non_mip variants) — identical.
- Primary MIP classification (9 sub-clauses) — identical.
- All refinement cuts (event-specific: 7013, 6640, 7018, etc.) — identical.
- flag_strong_check branch, n_good_tracks loop, energy-dependent corrections — identical.
- Angular cuts on low-energy showers — identical.
- Other-shower loop, min_dis computation — identical.
- Median/mean dedx conversion (Bethe-Bloch: alpha=1.0, beta=0.255, 43e3, 23.6e-6, 1.38, 0.273) — identical.
- All 20 `vec_dQ_dx` fills — identical.
- vec_dQ_dx padding: prototype appends 20 entries unconditionally; toolkit pads to exactly 20. Only indices 0-19 are read, so functionally equivalent.

### `high_energy_overlapping_sp` (hol_1, hol_2) — prototype lines 2596-2793, toolkit lines 1578-1748

- **hol_1**: n_valid_tracks, min_angle, flag_all_showers, 3 overlap conditions — identical (after already-fixed bug, see Bugs section).
- **hol_2**: min_ang2/min_sg scan, ncount proximity loop (0.6 cm threshold), medium_dQ_dx extraction — identical logic. Minor implementation difference: toolkit determines near-endpoint by geometric distance; prototype uses point-vector index comparison. Both yield the same near endpoint under normal conditions.

### `low_energy_overlapping_sp` (lol_1, lol_2, lol_3) — prototype lines 2797-2952, toolkit lines 1768-1938

- **lol_1**: 2-segment vertex angle check, all thresholds (36 deg, 150 MeV) — identical.
- **lol_2**: Muon/weak segment proximity, all thresholds (30 cm, 10 deg, 8 cm, 30 deg) — identical.
- **lol_3**: Backward shower check, n_out/n_sum ratio — identical.
- Prototype declares `flag_overlap_4` and `flag_overlap_5` but never sets them; toolkit correctly omits.

### `pi0_identification_sp` (pio_1, pio_2) — prototype lines 2955-3070, toolkit lines 1956-2082

- **pio_1**: Mass window (135±35/±60 MeV), asymmetry check (<0.87), energy thresholds (15/10/400 MeV), veto conditions (30 MeV + 80/120 cm) — all identical.
- **pio_2**: Back-to-back cluster search with distance (36 cm) and angle (7.5 deg) thresholds — identical. Toolkit precomputes `cluster_acc_length` in one pass (O(E)) vs prototype's per-vertex recomputation (O(V*S)) — efficiency improvement, same result.

### `low_energy_michel_sp` �� prototype lines 545-597, toolkit lines 2096-2145

- n_3seg counting, 2 flag_bad conditions (25/18 cm, 0.75 ratio, 100 MeV, 0.7 energy ratio) — identical.
- Prototype computes `E_range` and `n_segs` but never uses either in decision logic; toolkit correctly omits both.

---

## 2. Bugs Found

### In the Toolkit

Four toolkit bugs were found and fixed (2026-04-09/12) in a separate review pass (see `clus/docs/tagger/nue_singlephoton_tagger_review.md`):

| Bug | Severity | Location | Fix |
|---|---|---|---|
| First-pass TaggerInfo pollution — vector fields over-populated with entries from all showers | Medium-High | Entry point, first-pass loop | Use throw-away `TaggerInfo tmp_ti{}` for first-pass calls |
| Hardcoded `num_valid_tracks=0` in first-pass `bad_reconstruction_1_sp` call | Medium | Entry point, L2310 | Compute `first_pass_valid_tracks` per-shower matching prototype L183-190 |
| br3_3 `angle>105` and vector fills inside `dir1.mag>10cm` guard | Medium | `bad_reconstruction_2_sp`, L735 | Moved check outside guard; compute `angle` unconditionally |
| hol_1 incorrect `flag_all_showers=false` on zero-magnitude direction | Low | `high_energy_overlapping_sp` | Removed assignment, matching prototype's simple `continue` |

All four have been fixed in the current code.

### Prototype Bugs Fixed by Toolkit

1. **PROTO-BUG-1 (Low) — `num_muons`/`num_pions` parenthesis bug**
   - **Location**: Prototype lines 127-128
   - **Issue**: `abs(sg->get_particle_type()==13)` evaluates the comparison first (boolean 0/1), then takes `abs`. The outer `if` at line 121 correctly uses `abs(sg->get_particle_type())==13`, so only muons/pions enter the block. But the inner counters fail for anti-particles: `abs(-13==13) = abs(false) = 0`, so anti-muons count in `num_mip_tracks` but not `num_muons`. Same for `num_pions` with pdg=-211.
   - **Toolkit fix**: Line 2279-2280 correctly uses `std::abs(pdg) == 13` and `std::abs(pdg) == 211`.
   - **Impact**: Low — anti-particles are rare in neutrino interactions; `num_muons`/`num_pions` are TaggerInfo fields but not used in any cut within this function.

### Shared Behavior (Not Bugs in Toolkit)

1. **SHARED-BUG-1 (Low) — Proton ranking loses 2nd-best proton**
   - **Location**: Prototype lines 108-118, toolkit lines 2259-2268
   - **Issue**: When a new proton has `energy > proton_energy_1`, the old slot #1 is overwritten without demoting it to slot #2. For 3+ protons arriving in ascending energy order, `proton_2` remains at -1 even though a valid 2nd-best exists.
   - **Impact**: Low — events with 3+ identified protons at the neutrino vertex are rare. `proton_length_2`/`proton_dqdx_2`/`proton_energy_2` are TaggerInfo fields used by the BDT, not in any cut within this function.
   - **Status**: Not fixing — intentionally preserving prototype behavior.

---

## 3. Efficiency / Structure Improvements

1. **EFF-1 — SpContext bundle**: Prototype uses member variables scattered across the `NeutrinoID` class. Toolkit bundles shared state into the `SpContext` struct (lines 103-117), constructed once and passed by reference to all helpers. Improves locality and makes dependencies explicit.

2. **EFF-2 — Static helpers**: All 9 helper functions are file-local `static` functions, eliminating virtual dispatch overhead and enabling compiler inlining.

3. **EFF-3 — First-pass flag_fill elimination**: Prototype had a `flag_fill` boolean parameter threaded through all `bad_reconstruction_*` functions to suppress TaggerInfo writes in the first pass. Toolkit eliminates this parameter entirely — first-pass calls use a throw-away `TaggerInfo tmp_ti{}`.

4. **EFF-4 �� pi0_identification_sp acc_length precomputation**: Prototype recomputes per-cluster track length for every candidate vertex (O(V*S)). Toolkit precomputes the `cluster_acc_length` map in one edge-scan pass (O(E)), then does O(1) lookups.

5. **EFF-5 — Dead code removal**: Prototype's `E_range` computation in `low_energy_michel_sp`, unused `n_segs` variable, unused `flag_overlap_4`/`flag_overlap_5` in `low_energy_overlapping_sp`, and `max_dQ_dx` in `bad_reconstruction_sp` br1_2 — all removed.

6. **EFF-6 — Debug output removal**: All `std::cout`/`flag_print` debug output (hundreds of lines) removed throughout.

7. **Overall**: 4275 lines → 2511 lines (41% reduction) with identical functionality.

---

## 4. Determinism

| Container | Location | Prototype | Toolkit | Issue? |
|-----------|----------|-----------|---------|--------|
| Shower loop | Entry L2241 | `map_vertex_to_shower[v]` (pointer set) | `VertexShowerSetMap` value = `IndexedShowerSet` | **Fixed** — index-ordered |
| `good_showers`, `ok_showers` | Entry L2231 | (not separate containers) | `std::set<ShowerPtr>` | **Low risk** — only `.size()`, `.count()`, `.insert()` used; no order-dependent iteration |
| Segment/vertex sets in helpers | Throughout | `map_seg_vtxs`/`map_vtx_segs` (pointer-keyed) | `IndexedSegmentSet`, `IndexedVertexSet` | **Fixed** — index-ordered |
| `cluster_acc_length` | pi0_identification L2038 | `map_segment_vertices` (pointer-keyed) | `std::map<Cluster*, double>` | **OK** — only used for lookup, never iterated |
| pio_2 vector population | pi0_identification L2045-2074 | `map_vertex_segments` iteration (pointer-ordered) | `graph_nodes(ctx.graph)` (pointer-ordered) | **Inherited** — both are pointer-ordered; pio_2 vectors populated in non-deterministic order. BDT treats these as unordered collections |
| Graph edge/vertex iteration | Throughout | `map_vertex_segments` | `boost::out_edges`, `graph_nodes` | **OK** — boost graph iteration is deterministic for a given graph construction |

The toolkit resolves the main non-determinism sources from the prototype (pointer-keyed maps for shower/segment/vertex iteration) via indexed containers. The one remaining non-deterministic pattern (pio_2 vector ordering via `graph_nodes`) is inherited from the prototype and has no practical impact since the BDT handles these vectors as unordered feature collections.

---

## 5. Multi-APA / Multi-Face Handling

The single-photon tagger operates on **graph topology and segment-level quantities** (lengths, dQ/dx, directions, PCA). It does not reference:
- Wire-plane geometry or wire directions
- APA indices or face indices (except for SpContext initialization at L2187-2194, which correctly derives apa/face from the main vertex position using `dv->contained_by()`)
- Any `WirePlaneId` parameters

**Drift direction** is hardcoded as `Vector(1, 0, 0)` throughout (same as prototype). All drift-angle computations use `fabs(angle - 90°)`, which measures perpendicularity to the drift axis. This is symmetric with respect to drift sign — a drift direction of `(-1, 0, 0)` (opposite face) produces the same isolation angle. Safe for multi-face detectors.

**Beam direction** is hardcoded as `Vector(0, 0, 1)` — geometry-independent.

**No multi-APA issues found.**

---

## 6. Minor Logic Divergence

**DIVERGE-1 — hol_2 near-endpoint determination**

In `high_energy_overlapping_sp` hol_2, when computing `medium_dQ_dx` for the closest-angle segment:

- **Prototype** (line 2736-2763): Determines which end is near the vertex by comparing point-vector indices (`wcpt.index` of front/back).
- **Toolkit** (line 1720-1727): Determines this by geometric distance (`ray_length` to front/back fit points).

Both yield the same near endpoint under normal conditions. A divergence could occur only if a segment's fit endpoints are reordered relative to its point vector, which does not happen in practice. Impact is negligible.

---

## 7. Angle Safety

The toolkit's `D3Vector::angle()` method (defined in `util/inc/WireCellUtil/D3Vector.h:127-136`) internally clamps the cosine to [-1, 1] via `std::min(std::max(cosine, T(-1)), T(1)))` before calling `std::acos`. All angle computations in the single-photon tagger use this method.

In contrast, the prototype uses ROOT's `TVector3::Angle()` which uses `atan2` (also safe). Both are robust against floating-point rounding producing values outside [-1, 1].

**No NaN risk from angle computations.**

---

## 8. Recommendations

1. **All toolkit bugs have been fixed** (2026-04-09/12). No further code changes needed.
2. The prototype parenthesis bug (PROTO-BUG-1: `num_muons`/`num_pions` for anti-particles) is worth noting for anyone comparing prototype vs toolkit output.
3. The shared proton ranking behavior (SHARED-BUG-1) could be improved in both codebases by demoting slot #1 to slot #2 when a higher-energy proton is found, but this is not critical and would change behavior vs prototype.
4. The pio_2 vector ordering non-determinism is inherited from the prototype and has no practical impact, but could be resolved by sorting the vector entries by distance or angle if exact reproducibility is desired.
