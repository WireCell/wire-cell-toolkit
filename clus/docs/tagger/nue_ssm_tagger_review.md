# Code Review: `ssm_tagger`

**Date:** 2026-04-09  
**Reviewer:** Claude Code  
**Scope:** Logic fidelity vs prototype, bug hunting, efficiency, determinism.

---

## Files Examined

| Role | File |
|---|---|
| Toolkit implementation | `clus/src/NeutrinoTaggerSSM.cxx` |
| Prototype implementation | `prototype_pid/src/NeutrinoID_ssm_tagger.h` |

---

## Functions Reviewed

| Function | Toolkit Location | Prototype Location |
|---|---|---|
| `ssm_tagger` | `NeutrinoTaggerSSM.cxx:415` | `NeutrinoID_ssm_tagger.h:1` |

Supporting helpers reviewed: `get_scores` (L42), `get_scores_bp` (L70), `get_containing_shower_info` (L101), `fill_ssmsp` (L117), `fill_ssmsp_pseudo_1/2/3` (L161–264), `find_incoming_segment` (L267), `fill_ssmsp_all` (L281).

---

## No bugs in the toolkit. Three prototype bugs correctly fixed.

---

## Structure

The function is large (~1500 lines in toolkit, ~3600 in prototype). Both share the same five-phase layout:

1. **Phase A** — Scan all segments at `main_vertex`; identify SSM candidates (length/PDG/dQ/dx cut, vertex-activity detection).
2. **Phase B** — Select best SSM (longest with vertex activity wins); compute SSM properties (direction, dQ/dx profile, PID scores, energy).
3. **Phase C** — Loop primary segments (at `main_vertex`, excluding SSM) and daughter segments (at `second_vtx`); rank top-2 tracks and top-2 showers in each category.
4. **Phase D** — If `backwards_muon`, swap all prim↔daughter variables; compute momentum vectors; run off-vertex loop.
5. **Phase E** — Fill `ssmsp` space-point vectors; evaluate `flag_st_kdar`; assign all TaggerInfo fields; return.

---

## Findings

### ✅ Phase A — SSM candidate scan

Pre-cuts (PDG, length, dQ/dx), d(dQ/dx) vector construction, vertex-activity scan of first/last 5 points, and the multi-stage ambiguity resolution (prototype lines 476–518) all match exactly.

Prototype uses `map_vertex_segments[main_vertex]` (ordered set); toolkit uses `boost::out_edges(main_vertex, graph)`. Order is nondeterministic in both (set order depends on pointer address). The `all_ssm_sg` map uses `SegmentIndexCmp` (by graph index) in the toolkit to give deterministic iteration — a correctness improvement over the prototype. ✅

---

### ✅ Phase B — SSM selection and property computation

- Selection policy (longest-with-vtx-activity wins) — identical. ✅
- `dir = ssm_sg->dirsign(); if (backwards_muon) dir = -dir` — matches `sg->get_flag_dir()` / negation in prototype. ✅
- Direction vectors at 5/10/15/20 cm — identical computation. ✅
- Angles filled from 10 cm direction (`init_dir_10`) — matches prototype's `angle_to_z_10` etc. ✅

---

### 🐛 Prototype bug fixed: `max_d_dq_dx_fwd_3` swap (prototype line 716)

**Prototype (when `dir == -1`):**

```cpp
std::swap(max_dq_dx_fwd_3,  max_dq_dx_bck_3);      // line 714 ✓
std::swap(max_dq_dx_fwd_5,  max_dq_dx_bck_5);      // line 715 ✓
std::swap(max_d_dq_dx_fwd_3, max_d_dq_dx_bck_5);  // line 716 ✗ (bck_5 should be bck_3)
std::swap(max_d_dq_dx_fwd_5, max_d_dq_dx_bck_5);  // line 717 (bck_5 is now corrupted)
```

Line 716 is a copy-paste typo: `max_d_dq_dx_bck_5` instead of `max_d_dq_dx_bck_3`. As a result, when `dir == -1`, the final state is:

| Variable | Prototype result | Correct result |
|---|---|---|
| `max_d_dq_dx_fwd_3` | original `bck_5` | original `bck_3` |
| `max_d_dq_dx_fwd_5` | original `fwd_3` | original `bck_5` |
| `max_d_dq_dx_bck_3` | original `bck_3` (unswapped) | original `fwd_3` |
| `max_d_dq_dx_bck_5` | original `fwd_5` | original `fwd_5` |

**Toolkit (correct):**

```cpp
std::swap(max_d_fwd3, max_d_bck3);   // correct
std::swap(max_d_fwd5, max_d_bck5);   // correct
```

These fields (`ssm_max_d_dq_dx_fwd_3/5/bck_3/5`) are BDT input features, so the bug affects BDT scores only when the SSM direction is backward. **The toolkit produces the physically correct values.**

---

### 🐛 Prototype bug fixed: `dQ_dx_cut` formula (prototype line 743)

**Prototype:**

```cpp
dQ_dx_cut = 0.8866 + 0.9533 * pow(18 / length / units::cm, 0.4234);
```

where `length` is already in cm. This evaluates to `pow(1.8 / length_cm, 0.4234)` — a factor-of-10 too small in the argument.

**Toolkit:**

```cpp
double dQ_dx_cut = 0.8866 + 0.9533 * std::pow(18.0 / length, 0.4234);
```

This evaluates to `pow(18 / length_cm, 0.4234)`, which is consistent with the `single_shower` formula `pow(18*cm / length_raw, 0.4234)`. **The toolkit is correct.**

`ssm_dQ_dx_cut` is stored as a feature only — never used for any cut within `ssm_tagger` — so the prototype error does not propagate to any event-selection decision within this function.

---

### 🐛 Prototype bug fixed: `x/y/z_dir_daught_track1` copy-paste error (prototype lines 1492–1494)

**Prototype:**

```cpp
if (length_daught_track1 > 0) {
    ...
    x_dir_daught_track1 = dir_prim_track1[0];   // ✗ wrong vector
    y_dir_daught_track1 = dir_prim_track1[1];   // ✗ wrong vector
    z_dir_daught_track1 = dir_prim_track1[2];   // ✗ wrong vector
    mom_daught_track1[...] = dir_daught_track1[...] * ...;  // momentum is correct
}
```

`x/y/z_dir_daught_track1` is read from `dir_prim_track1` instead of `dir_daught_track1` — a copy-paste error from the `prim_track1` block immediately above.

**Toolkit (correct):**

```cpp
ti.ssm_daught_track1_x_dir = (float)dir_daught_track1.x();
ti.ssm_daught_track1_y_dir = (float)dir_daught_track1.y();
ti.ssm_daught_track1_z_dir = (float)dir_daught_track1.z();
```

Uses `dir_daught_track1` throughout. **The toolkit is correct.**

---

### ✅ Phase C — Primary/daughter loops

- `n_prim_tracks_*` / `n_prim_all_*` / `n_daughter_tracks_*` / `n_daughter_all_*` counting — matches prototype exactly. ✅
- Top-2-by-length selection for prim_track1/2, prim_shw1/2, daught_track1/2, daught_shw1/2 — identical ranking logic. ✅
- `daughter_counts` helper subtracts self from track count, does not subtract from shower track count (for shower segments), adds back if segment too short — exactly mirrors prototype's `calculate_num_daughter_tracks` pattern. ✅
- `kine_energy_best_*` for showers: uses `shower->get_kine_best()` / `get_kine_charge()` fallback — matches prototype. ✅

---

### ✅ Phase D — Backwards-muon swap, momentum, off-vertex loop

#### Backwards-muon swap

All prim↔daughter swaps (counts, lengths, energies, directions, scores, segment pointers) match the prototype's 120-line block exactly. The toolkit also includes the segment pointer swaps (`prim_track1_sg ↔ daught_track1_sg` etc.) that the prototype does at line 1460–1463. ✅

#### Momentum vectors

`p = sqrt(KE² + 2·KE·m)`, `m` from `particle_data->get_particle_mass()` — matches prototype's `sg->get_particle_mass()`. ✅

#### Off-vertex loop

- Prototype iterates `map_segment_vertices` (all segments), skipping those at `ssm_main_vtx` or `ssm_second_vtx`. Toolkit iterates `boost::edges(graph)` with the same filter. ✅
- Distance filter: `sep > 80 cm` — matches prototype's `sep_dist > 80`. ✅
- Track vs. shower classification: `pdg == 22 || == 11` → shower; else → track — matches prototype. ✅
- Single best off-vertex track and single best off-vertex shower by length. ✅

---

### ✅ `flag_st_kdar` condition

```cpp
// Both prototype and toolkit:
n_prim_tracks_1 == 0 && n_prim_all_3 == 0 && n_daughter_tracks_5 == 0 &&
n_daughter_all_5 < 2 && Nsm_wivtx == 1 &&
!(pio_mass > 70 && pio_mass < 200)
```

✅ Matches exactly.

---

### ✅ Exit path (exit_ssm / exit_ssm_tagger)

All `ssm_*` fields set to −999, `Nsm`/`Nsm_wivtx` always written before exit. Matches prototype's `exit_ssm_tagger`. ✅

---

### ✅ SSMSP space-point filling

`fill_ssmsp_all` — BFS from start vertex, SSM segment first, then all connected segments, then showers by connection type (1/2/3). All three `fill_ssmsp_pseudo` variants (`vtx→shower_start`, `shower_start→parent_sg`, `mother_sg→daughter_shower_start`) match the prototype's three `fill_ssmsp_psuedo` overloads. ✅

---

## Summary Table

| Check | Result | Action |
|---|---|---|
| Phase A — candidate scan: PDG/length/dQ/dx pre-cuts | ✅ Match | — |
| Phase A — d(dQ/dx) vertex-activity detection and ambiguity resolution | ✅ Match | — |
| Phase A — `all_ssm_sg` iteration determinism (SegmentIndexCmp) | ✅ Toolkit improvement | — |
| Phase B — SSM selection (longest with vtx activity) | ✅ Match | — |
| Phase B — direction vectors at 5/10/15/20 cm | ✅ Match | — |
| Phase B — angles stored from 10 cm direction | ✅ Match | — |
| Phase B — break-point detection and reduced length | ✅ Match | — |
| Phase B — dQ/dx fwd/bck point-by-point values | ✅ Match | — |
| Phase B — `max_d_dq_dx_fwd/bck_3/5` swap when dir==-1 | 🐛 Prototype typo (bck_5 vs bck_3) | Toolkit correct |
| Phase B — PID scores and dir-based swap | ✅ Match | — |
| Phase B — `score_*_bp` with break_point | ✅ Match | — |
| Phase B — degenerate break_point catch | ✅ Match | — |
| Phase B — `dQ_dx_cut` formula | 🐛 Prototype has extra `/units::cm` | Toolkit correct |
| Phase B — `medium_dq_dx_bp` | ✅ Match | — |
| Phase B — kinetic energy from range (muon, PDG=13) | ✅ Match | — |
| Phase C — prim/daughter count accumulators | ✅ Match | — |
| Phase C — top-2 primary tracks and showers selection | ✅ Match | — |
| Phase C — top-2 daughter tracks and showers selection | ✅ Match | — |
| Phase C — `daughter_counts` subtract/add-back logic | ✅ Match | — |
| Phase C — `kine_energy_best_*` shower lookup | ✅ Match | — |
| Phase D — backwards_muon swap (all variable pairs incl. segment ptrs) | ✅ Match | — |
| Phase D — `x/y/z_dir_daught_track1` direction source | 🐛 Prototype uses `dir_prim_track1` | Toolkit correct |
| Phase D — momentum vector computation | ✅ Match | — |
| Phase D — `ssm_main_vtx`/`ssm_second_vtx` assignment when backwards | ✅ Match | — |
| Phase D — off-vertex loop: distance filter, classification, single best | ✅ Match | — |
| Phase E — `flag_st_kdar` cut conditions | ✅ Match | — |
| Phase E — `fill_ssmsp_all` BFS + shower handling | ✅ Match | — |
| Phase E — all TaggerInfo fills (~200 fields) | ✅ Match | — |
| Exit path — all ssm_* set to -999 | ✅ Match | — |

---

## Changes Made

None. All three divergences from the prototype are prototype bugs that the toolkit already has correct. No toolkit code was modified.
