# Code Review: `init_tagger_info()` and `fill_kine_tree()`

**Date:** 2026-04-08  
**Reviewer:** Claude Code  
**Scope:** Logic fidelity vs prototype, bug hunting, efficiency, and determinism.

---

## Files Examined

| Role | File |
|---|---|
| Toolkit implementation | `clus/src/NeutrinoKinematics.cxx` |
| Toolkit struct definitions | `clus/inc/WireCellClus/NeutrinoTaggerInfo.h` |
| Toolkit shower class | `clus/inc/WireCellClus/PRShower.h` |
| Toolkit call site | `clus/src/TaggerCheckNeutrino.cxx:315-318` |
| Prototype implementation | `prototype_pid/src/NeutrinoID.cxx:2217-2550` (`init_tagger_info`) |
| Prototype implementation | `prototype_pid/src/NeutrinoID_kine.h:1-283` (`fill_kine_tree`) |

---

## `init_tagger_info()`

### Prototype (NeutrinoID.cxx:2217)

Approximately 300 lines of explicit per-field assignments on the `tagger_info` member, covering every field with its default value.

### Toolkit (NeutrinoKinematics.cxx:18-21)

```cpp
void PatternAlgorithms::init_tagger_info(TaggerInfo& ti)
{
    ti = TaggerInfo{};
}
```

The toolkit relies on C++ default-member-initializers in `NeutrinoTaggerInfo.h`. Value-initialising the struct (`TaggerInfo{}`) fires all those initializers in one step.

### Findings

| Check | Result | Detail |
|---|---|---|
| Logic fidelity | ✅ Pass | All non-zero defaults in `NeutrinoTaggerInfo.h` have been cross-checked against the prototype's explicit assignments. Every sentinel (-999), flag (1), and other non-zero value matches. |
| Bugs | ✅ None | — |
| Efficiency | ✅ Optimal | Single aggregate copy; no per-field assignments needed. |
| Randomness | ✅ N/A | No container iteration involved. |

### Notable non-zero defaults verified

| Group | Fields with non-trivial defaults |
|---|---|
| Cosmic / gap / MIP flags | `cosmic_flag{1}`, `gap_flag{1}`, `mip_quality_flag{1}`, `mip_flag{1}` |
| MIP dQ/dx sentinels | `mip_n_first_non_mip{19}`, `mip_n_first_non_mip_1{19}`, `mip_n_first_non_mip_2{19}`, `mip_n_below_threshold{19}` |
| MIP energy sentinels | `mip_vec_dQ_dx_0{1}`, `mip_vec_dQ_dx_1{1}`, `mip_max_dQ_dx_sample{1}`, `mip_lowest_dQ_dx{1}`, etc. |
| MIP counts | `mip_n_vertex{1}`, `mip_n_lowest{1}`, `mip_n_highest{1}`, `mip_n_other_vertex{2}`, `mip_n_stem_size{20}` |
| SSM sentinels | All `ssm_*` fields initialised to `-999` matching prototype |
| Single-photon | `shw_sp_flag{1}`, all dQ/dx sentinels matching MIP group |

**Conclusion:** `init_tagger_info()` is correct. No action required.

---

## `fill_kine_tree()`

### Overview of structure (both prototype and toolkit)

1. Collect SCE-corrected neutrino vertex position.
2. Mark all shower-internal vertices/segments as used.
3. Build a `map_sg_shower` (start-segment → shower).
4. **First pass**: iterate edges from `main_vertex`; dispatch each to shower or track handling.
5. **BFS loop**: traverse remaining track graph; accumulate energies and mass corrections.
6. **Remaining showers**: handle showers not reached in BFS (secondary showers with `vtx_type ≤ 3`).
7. Sum `kine_reco_Enu`; copy π⁰ kinematics.

### Findings

---

#### ✅ Logic fidelity — determinism

The prototype uses raw-pointer `std::set<ProtoVertex*>` and `std::map<ProtoVertex*, ...>`, whose iteration order depends on pointer addresses (non-deterministic across runs). The toolkit uses:

- `IndexedVertexSet`, `IndexedSegmentSet` — sorted by stable integer index via `VertexIndexCmp`, `SegmentIndexCmp`
- `IndexedShowerSet` — sorted by `ShowerIndexCmp`
- `std::map<SegmentPtr, ShowerPtr, SegmentIndexCmp>` for `map_sg_shower`

This eliminates pointer-address-dependent ordering throughout the function. **This is a correctness improvement over the prototype.**

---

#### ✅ Logic fidelity — particle type access

The prototype accesses `shower->get_start_segment()->get_particle_type()` at every push site.  
The toolkit uses `shower->get_particle_type()` (returns `data.particle_type`).

Cross-check: the only call site for `shower->set_particle_type()` in the toolkit is `NeutrinoShowerClustering.cxx:118`:

```cpp
shower->set_particle_type(curr_sg->particle_info()->pdg());
```

`curr_sg` is the start segment at that point, so `shower->get_particle_type()` is always identical to `shower->start_segment()->particle_info()->pdg()`. **Equivalent to prototype.**

---

#### ✅ Logic fidelity — shower duplication guard in BFS

The prototype's BFS loop processes showers found during traversal **without** checking `used_showers` first. If a shower's start segment is reachable from multiple BFS paths, the prototype would double-count it.

The toolkit adds an explicit guard (NeutrinoKinematics.cxx:233):

```cpp
if (!used_showers.count(shower)) {
    push_shower_kine(shower);
    ktree.kine_energy_included.push_back(1);
    used_showers.insert(shower);
}
```

**This is a correctness improvement over the prototype.**

---

#### ✅ Logic fidelity — `kine_energy_included` for remaining showers

Toolkit (NeutrinoKinematics.cxx:286):
```cpp
ktree.kine_energy_included.push_back(vtx_type != 3 ? 1 : vtx_type);
```

Prototype (NeutrinoID_kine.h:233-237):
```cpp
if (pair_vertex.second != 3)
    ktree.kine_energy_included.push_back(1);
else
    ktree.kine_energy_included.push_back(pair_vertex.second);
```

Both produce: push `1` unless `vtx_type == 3`, in which case push `3`. **Equivalent.**

---

#### ✅ Logic fidelity — rest-mass correction for remaining showers

Prototype (NeutrinoID_kine.h:240-245): only proton binding energy with a 5 cm length cut; non-electron non-proton showers receive no correction (consistent with the comment that they are not present in the remaining set).

Toolkit (NeutrinoKinematics.cxx:289-294): same logic — only proton binding energy with `segment_track_length(start_sg) > 5 * units::cm`. **Equivalent.**

---

#### ✅ Logic fidelity — order of guards in remaining showers loop

Toolkit checks `used_showers` **before** calling `get_start_vertex_and_type()`.  
Prototype checks `vtx_type` first.

Both conditions must be false to proceed, so correctness is identical. The toolkit order is actually **more efficient**: it skips the `get_start_vertex_and_type()` call for showers that are already used (the common case after the BFS).

---

#### BUG FIXED — redundant `kine_best` fallback in `push_shower_kine` and remaining showers

**Location:** NeutrinoKinematics.cxx:105 and NeutrinoKinematics.cxx:272 (pre-fix)

**Problem:**

```cpp
double kine_best = shower->get_kine_best();
if (kine_best == 0) kine_best = shower->get_kine_charge();  // dead code
```

`PRShower::get_kine_best()` (PRShower.h:172-173) already performs this fallback internally:

```cpp
double get_kine_best(){
    if (data.kenergy_best != 0) return data.kenergy_best;
    else return data.kenergy_charge;
}
```

So `kine_best` can never be 0 here unless `kine_charge` is also 0, in which case `get_kine_charge()` returns the same 0. The `if` branch was dead code.

**Fix:** Removed both redundant lines (L105 and L272).

The prototype needed this check because the prototype's `WCShower::get_kine_best()` did **not** perform the internal fallback. The toolkit's `PRShower::get_kine_best()` was designed to include it, making the external check obsolete.

---

#### FIXED — wasted dQdx computation in `push_segment_kine`

**Location:** NeutrinoKinematics.cxx:144 (pre-fix)

**Problem:**

```cpp
(void)segment_cal_kine_dQdx(seg, recomb_model); // computed but not yet filled in ktree
```

The prototype also computes `kine_dQdx` at every track site but never uses it in the `kine_energy_info` decision — the ternary only compares `kine_best` against `kine_charge` and `kine_range`, with dQdx as an implicit fallback (`push 0`). The `(void)` cast with comment revealed this was already known to be unused.

If `segment_cal_kine_dQdx` had side effects needed here, those would have been documented. It has no required side effects in this context.

**Fix:** Removed the call entirely.

---

### Summary table for `fill_kine_tree()`

| Check | Result | Action |
|---|---|---|
| Determinism (ordered containers) | ✅ Improved over prototype | — |
| Particle-type access | ✅ Equivalent | — |
| BFS shower duplication guard | ✅ Improved over prototype | — |
| `kine_energy_included` logic | ✅ Equivalent | — |
| Rest-mass correction (remaining showers) | ✅ Equivalent | — |
| Redundant `kine_best` fallback | ⚠️ Dead code | **Fixed** |
| Wasted dQdx computation | ⚠️ Wasted CPU | **Fixed** |

---

## Changes Made

**File:** `clus/src/NeutrinoKinematics.cxx`

1. Removed `if (kine_best == 0) kine_best = shower->get_kine_charge();` from `push_shower_kine` lambda (was L105) — dead code because `get_kine_best()` already handles this.
2. Removed `if (kine_best == 0) kine_best = shower->get_kine_charge();` from remaining-showers loop (was L272) — same reason.
3. Removed `(void)segment_cal_kine_dQdx(seg, recomb_model);` from `push_segment_kine` lambda (was L144) — result was always discarded; no side effects required.
