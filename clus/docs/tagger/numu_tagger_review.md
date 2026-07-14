# Code Review: `numu_tagger()`, `count_daughters()`, `find_cont_muon_segment()`

**Date:** 2026-04-08  
**Reviewer:** Claude Code  
**Scope:** Logic fidelity vs prototype, bug hunting, efficiency, determinism.

---

## Files Examined

| Role | File |
|---|---|
| Toolkit implementation | `clus/src/NeutrinoTaggerNuMu.cxx` |
| Toolkit helper | `clus/src/NeutrinoVertexFinder.cxx:920` |
| Prototype | `prototype_pid/src/NeutrinoID_numu_tagger.h` |
| Prototype helper | `prototype_pid/src/NeutrinoID_track_shower.h:2304` |

---

## Functions Reviewed

| Function | Toolkit Location | Prototype Location |
|---|---|---|
| `count_daughters` (segment overload) | `NeutrinoTaggerNuMu.cxx:81` | `NeutrinoID_numu_tagger.h:264` |
| `count_daughters` (shower overload) | `NeutrinoTaggerNuMu.cxx:119` | `NeutrinoID_numu_tagger.h:294` |
| `numu_tagger` | `NeutrinoTaggerNuMu.cxx:162` | `NeutrinoID_numu_tagger.h:1` |
| `find_cont_muon_segment` | `NeutrinoVertexFinder.cxx:920` | `NeutrinoID_track_shower.h:2304` |

> **Note:** `find_cont_muon_segment_nue` (a related but distinct function) was already reviewed in `cosmic_tagger_review.md`. This review covers `find_cont_muon_segment`, which is called by `numu_tagger` during the flag-3 muon extension check.

---

## BUG FIXED — Unused structured-binding variable in `numu_tagger` flag 3

**Location:** `NeutrinoTaggerNuMu.cxx` L327-332 (pre-fix)

**Problem:**

```cpp
auto [ext_sg1, ext_vtx1] = find_cont_muon_segment(graph, sg, pv1);
if (ext_sg1) tmp_length += segment_track_length(ext_sg1);
...
auto [ext_sg2, ext_vtx2] = find_cont_muon_segment(graph, sg, pv2);
if (ext_sg2) tmp_length += segment_track_length(ext_sg2);
```

`ext_vtx1` and `ext_vtx2` are named bindings that are never used. The prototype returns a `pair` and uses only `.first`; the structured bindings create named-but-unused variables that produce compiler warnings. No correctness impact, but a hygiene issue.

**Fix:** Replaced with explicit `.first` access:

```cpp
SegmentPtr ext_sg1 = find_cont_muon_segment(graph, sg, pv1).first;
if (ext_sg1) tmp_length += segment_track_length(ext_sg1);
...
SegmentPtr ext_sg2 = find_cont_muon_segment(graph, sg, pv2).first;
if (ext_sg2) tmp_length += segment_track_length(ext_sg2);
```

---

## `count_daughters` (segment overload)

### Structure

Counts track daughters and all daughters at the **far end** of a given muon segment:
1. Find both vertices of the segment.
2. Identify which vertex is closer to `main_vertex`.
3. Use `calculate_num_daughter_tracks` from that close vertex through the segment to count what is reachable at the far end.
4. Subtract the muon segment itself from the counts.

### Findings

#### ✅ Logic fidelity

Both toolkit and prototype find the two endpoint vertices, pick the one closer to `main_vertex` (by distance), call `calculate_num_daughter_tracks` from that vertex, and subtract 1 from each count for the muon segment itself. All conditions match. ✅

#### ✅ Unused `muon_dir` removed

Prototype (L265, L272-281): declares and populates `TVector3 muon_dir` but never uses it. Toolkit correctly omits `muon_dir`. **Dead code removal.**

#### ✅ Null guard improvement

Prototype (L268): `if (max_muon != 0)` — if null, returns `{0, 0}` implicitly. Toolkit (L84): early-return `if (!max_muon || !main_vertex) return {0, 0};` — also guards `main_vertex`, which the prototype assumes is always valid. **Defensive improvement.**

---

## `count_daughters` (shower overload)

### Structure

Same semantics as the segment overload, but for a long-muon shower chain:
1. Get the last segment of the chain via `get_last_segment_vertex_long_muon`.
2. Pick the endpoint of that segment closer to `main_vertex`.
3. Count daughters via `calculate_num_daughter_tracks`.
4. Subtract 1 from each count.

### Findings

#### ✅ Logic fidelity

Both toolkit and prototype call `get_last_segment_vertex_long_muon`, find the near/far endpoint, and apply `calculate_num_daughter_tracks`. All conditions match. ✅

#### ✅ Dead `muon_dir` computation removed

Prototype (L300): `muon_dir = max_long_muon->cal_dir_3vector(main_vertex->get_fit_pt(), 30*units::cm);` — computed but never used at any point in the function. Pure dead code. Toolkit correctly omits this. **Dead code removal.**

#### ✅ Null guard on `last_sg`

Toolkit adds `if (!last_sg) return {0, 0};` after `get_last_segment_vertex_long_muon` returns. Prototype would crash if `last_sg` is null. **Defensive improvement.**

---

## `numu_tagger`

### Structure

Three independent checks:

| Flag | Description |
|---|---|
| flag 1 (`numu_cc_1`) | Muon-like segment directly at `main_vertex`: PDG=13, length > 5 cm, dQ/dx below threshold, directness or length cut, not too many daughters |
| flag 2 (`numu_cc_2`) | Long-muon shower in main cluster: PDG=13 start segment, total track length > 18 cm, not too many daughters |
| flag 3 (`numu_cc_3`) | Muon-like segment NOT at `main_vertex`, length > 25 cm (or 30 cm if no other tracks), not too many daughters, not a pion |

Returns `{flag_long_muon, max_muon_length}` where `flag_long_muon = (max_muon_length > 100 cm || max_length_all > 120 cm) && numu_cc`.

### Findings

#### ✅ Logic fidelity — flag 1 cut condition

Prototype (L46):
```cpp
abs(sg->get_particle_type())==13 && length > dis_cut && medium_dQ_dx < dQ_dx_cut * 43e3/units::cm
  && (length > 40*units::cm || length <= 40*units::cm && direct_length > 0.925 * length)
```

`&&` binds tighter than `||`, so this parses as:
```
A && B && C && (D || (E && F))
```
i.e., `length > 40cm` OR `(length <= 40cm AND straight)`.

Toolkit (L208-212):
```cpp
std::abs(pdg) == 13 &&
length > dis_cut &&
medium_dQ_dx < dQ_dx_cut * 43e3/units::cm &&
(length > 40*units::cm || (length <= 40*units::cm && direct_length > 0.925 * length))
```

Explicitly parenthesized to the same meaning. ✅

#### ✅ Logic fidelity — flag 1 TaggerInfo fields

All eight fields (`numu_cc_flag_1`, `numu_cc_1_particle_type`, `numu_cc_1_length`, `numu_cc_1_direct_length`, `numu_cc_1_medium_dQ_dx`, `numu_cc_1_dQ_dx_cut`, `numu_cc_1_n_daughter_tracks`, `numu_cc_1_n_daughter_all`) match the prototype with correct unit conversions. ✅

#### ✅ Logic fidelity — flag 2

Prototype iterates `showers` (raw pointer container). Toolkit iterates `IndexedShowerSet& showers`. Same muon-pdg + same-cluster-id + total-track-length + daughter-count conditions. ✅

Null safety: toolkit guards `start_sg` null (via `if (!start_sg) continue`) and `cluster()` null. Prototype dereferences directly. **Defensive improvement.**

#### ✅ Logic fidelity — flag 3 shower condition simplification

Prototype (L138):
```cpp
if (sg->get_flag_shower() && (!sg->get_flag_shower_topology()) || (!sg->get_flag_shower()) || length > 50*units::cm)
```

Operator precedence (`&&` before `||`): `(A && !B) || !A || C`.  
Boolean simplification: `(A && !B) || !A = !(A && B) = !A || !B`.  
So the condition is `!(shower && topo) || length > 50cm`.

Since `kShowerTopology` implies a shower segment — a segment cannot be topology-flagged without also being shower-flagged — `!(shower && topo)` is equivalent to `!topo`. Hence:

```
!is_shower_topo || length > 50cm
```

Toolkit (L309): `if (!is_shower_topo || length > 50*units::cm)` — same result, with an explanatory comment. ✅

#### ✅ Logic fidelity — flag 3 acc_track_length / max_length_all

Prototype updates `max_length_all` (L120) and `acc_track_length` (L125-127) **before** the `main_vertex` skip at L128. Toolkit updates both at L292-294 before the `seg_at_main_vertex` skip at L297. Order preserved. ✅

#### ✅ Logic fidelity — flag 3 final evaluation

```
(max_length > 25cm && acc_track_length > 0 ||
 max_length > 30cm && acc_track_length == 0) &&
(max_length_all - max_muon_length <= 100cm) &&
!(n3_tracks > 1 || n3_all - n3_tracks > 2)
```

Both prototype and toolkit: same conditions, same pion exclusion (`tmp_particle_type != 211`). ✅

#### ✅ Unused variables removed from flag 3

Prototype (L109-110): declares `max_direct_length` and `max_medium_dQ_dx`, updates them in the loop, but uses them only in the commented-out debug print (L211). Toolkit correctly omits both. **Dead code removal.**

#### ✅ Determinism — BGL traversal

Flag 1: prototype uses `map_vertex_segments[main_vertex]` (raw pointer set). Toolkit uses `boost::out_edges(vd, graph)`.  
Flag 2: prototype iterates raw pointer shower container. Toolkit uses `IndexedShowerSet` (sorted by index).  
Flag 3: prototype iterates `map_segment_vertices` (raw pointer map over all segments). Toolkit uses `boost::edges(graph)`.  
**All three improved to deterministic traversal.**

#### ✅ `neutrino_type` bit-setting omitted

Prototype (L252-255) sets `neutrino_type |= 1UL<<2` (numu) or `1UL<<3` (nc). Toolkit does not manage neutrino_type bits — the return value carries `flag_long_muon` and caller is responsible. This is a deliberate architectural decision documented in the file header. ✅

---

## `find_cont_muon_segment`

### Structure

Given a segment `sg` and one of its endpoint vertices `vtx`, find the adjacent segment (at `vtx`) that best continues `sg` as a muon — most collinear and consistent dQ/dx. Returns the best candidate and its far vertex, or `{nullptr, nullptr}` if none qualifies.

### Findings

#### ✅ Logic fidelity — angle computation

Prototype (L2328): `angle = (3.1415926 - dir1.Angle(dir2)) / 3.1415926 * 180.`  
Toolkit (L964): `angle = (M_PI - std::acos(cos_angle)) / M_PI * 180.0`

Same formula. ✅

#### ✅ Logic fidelity — angle_ok condition

Prototype (L2340):
```cpp
(angle < 10. || angle1 < 10 || sg_length < 6*units::cm && (angle < 15 || angle1 < 15))
```
Parses as: `angle<10 || angle1<10 || (sg_length<6cm && (angle<15 || angle1<15))`.

Toolkit (L982):
```cpp
(angle < 10.0 || angle1 < 10.0 || (sg_length < 6*units::cm && (angle < 15.0 || angle1 < 15.0)))
```
Explicitly parenthesized — same meaning. ✅

#### ✅ Logic fidelity — selection criterion

Both select the candidate with maximum `length * cos(angle)` (projected length along the continuing direction). ✅

#### ✅ Efficiency improvement — `std::clamp` prevents NaN

Prototype (L2328): `dir1.Angle(dir2)` uses ROOT's built-in which clamps internally.  
Toolkit (L963): `std::clamp(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()), -1.0, 1.0)` before `std::acos` — explicit clamp prevents NaN from floating-point rounding that pushes the dot product slightly outside `[-1, 1]`. **Minor robustness improvement.**

#### ✅ Safety — `descriptor_valid()` guard

Toolkit (L938): returns early if `!vtx->descriptor_valid()`. Prototype dereferences directly. **Defensive improvement.**

#### ✅ Determinism — BGL traversal

Prototype (L2318): `map_vertex_segments[vtx]` (raw pointer set). Toolkit (L944-946): `boost::out_edges(vd, graph)`. **Improvement.**

#### NOTE — `dir1` recomputed each iteration (matches prototype)

Both prototype and toolkit compute `dir1 = segment_cal_dir_3vector(sg, vtx_point, 15*units::cm)` inside the loop, even though `sg` is fixed. This is a minor inefficiency. Since the function has few adjacent segments in typical events and this matches the prototype exactly, it is not changed here.

---

## Summary Table

| Check | Result | Action |
|---|---|---|
| Unused `ext_vtx1`/`ext_vtx2` in `numu_tagger` flag 3 | ⚠ Compiler warning / hygiene | **Fixed** |
| `count_daughters` (seg): unused `muon_dir` | ✅ Correctly omitted | — |
| `count_daughters` (shower): dead `muon_dir` computation | ✅ Correctly omitted | — |
| `count_daughters` null guards | ✅ Improvement | — |
| `numu_tagger` flag 1 cut conditions | ✅ Equivalent | — |
| `numu_tagger` flag 1 TaggerInfo fields | ✅ Equivalent | — |
| `numu_tagger` flag 2 logic | ✅ Equivalent | — |
| `numu_tagger` flag 3 shower condition simplification | ✅ Equivalent (correctly simplified) | — |
| `numu_tagger` flag 3 `acc_track_length`/`max_length_all` ordering | ✅ Equivalent | — |
| `numu_tagger` flag 3 final evaluation | ✅ Equivalent | — |
| `numu_tagger` flag 3 unused `max_direct_length`, `max_medium_dQ_dx` | ✅ Correctly omitted | — |
| `numu_tagger` `neutrino_type` bits omitted | ✅ Deliberate architectural decision | — |
| Determinism — BGL traversal throughout | ✅ Improvement over prototype | — |
| `find_cont_muon_segment` angle computation | ✅ Equivalent | — |
| `find_cont_muon_segment` `std::clamp` for acos NaN prevention | ✅ Improvement | — |
| `find_cont_muon_segment` `descriptor_valid()` guard | ✅ Improvement | — |
| `find_cont_muon_segment` `dir1` recomputed per loop iteration | ℹ Minor inefficiency (matches prototype) | — |

---

## Changes Made

**File:** `clus/src/NeutrinoTaggerNuMu.cxx`

1. Replaced unused structured bindings `auto [ext_sg1, ext_vtx1]` and `auto [ext_sg2, ext_vtx2]` in `numu_tagger` flag-3 muon extension block with explicit `.first` access:
   ```cpp
   SegmentPtr ext_sg1 = find_cont_muon_segment(graph, sg, pv1).first;
   SegmentPtr ext_sg2 = find_cont_muon_segment(graph, sg, pv2).first;
   ```
   This eliminates potential compiler warnings for unused variables `ext_vtx1`/`ext_vtx2`.
