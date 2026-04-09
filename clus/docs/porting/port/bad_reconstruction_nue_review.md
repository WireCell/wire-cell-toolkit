# Code Review: `bad_reconstruction_1`, `bad_reconstruction_2`, `bad_reconstruction_3` (NuE tagger)

**Date:** 2026-04-08  
**Reviewer:** Claude Code  
**Scope:** Logic fidelity vs prototype, bug hunting, efficiency, determinism.

---

## Files Examined

| Role | File |
|---|---|
| Toolkit implementation | `clus/src/NeutrinoTaggerNuE.cxx` |
| Prototype | `prototype_pid/src/NeutrinoID_nue_tagger.h` |

---

## Functions Reviewed

| Function | Toolkit Location | Prototype Location | TaggerInfo Fields |
|---|---|---|---|
| `bad_reconstruction_1` | `NeutrinoTaggerNuE.cxx:2166` | `NeutrinoID_nue_tagger.h:3854` | `br2_*` |
| `bad_reconstruction_3` | `NeutrinoTaggerNuE.cxx:2860` | `NeutrinoID_nue_tagger.h:2907` | `br4_1_*`, `br4_2_*` |
| `bad_reconstruction_2` | `NeutrinoTaggerNuE.cxx:3077` | `NeutrinoID_nue_tagger.h:3145` | `br3_1_*` … `br3_8_*` |
| `bad_reconstruction` (br1) | `NeutrinoTaggerCosmic.cxx:101` | already reviewed in `cosmic_tagger_review.md` | — |

> **Note:** `bad_reconstruction` (br1) is not re-implemented in `NeutrinoTaggerNuE.cxx`. The NuE tagger calls `ctx.self.bad_reconstruction(...)` which dispatches to the implementation in `NeutrinoTaggerCosmic.cxx`, reviewed separately.

---

## BUGS FIXED

### 1. Redundant `get_kine_best()` fallback — three sites

**Locations (pre-fix):**
- `bad_reconstruction_1` L2172-2173
- `bad_reconstruction_3` L2866-2867
- `bad_reconstruction_2` L3085-3086

**Problem** (identical in all three):

```cpp
double Eshower = (shower->get_kine_best() != 0)
                 ? shower->get_kine_best() : shower->get_kine_charge();
```

`PRShower::get_kine_best()` already performs this fallback internally:

```cpp
double get_kine_best(){
    if (data.kenergy_best != 0) return data.kenergy_best;
    else return data.kenergy_charge;
}
```

Same class of dead code removed in `NeutrinoKinematics.cxx` and `NeutrinoTaggerCosmic.cxx` in previous sessions.

**Fix:** All three simplified to `double Eshower = shower->get_kine_best();`

---

### 2. Dead accumulation loop in `bad_reconstruction_2`, br3_5

**Location:** `bad_reconstruction_2` L3218-3226 (pre-fix)

**Problem:**

```cpp
Point ave_p(0, 0, 0); int num_p = 0, n_seg = 0;
double side_total_length = 0;
for (SegmentPtr sg1 : shower_segs) {                    // first loop — DEAD
    if (sg1->cluster() != sg->cluster() || sg1 == sg) continue;
    for (const auto& fit : sg1->fits()) {
        ave_p = ave_p + fit.point * (1.0 / 1.0);       // adds fit.point to ave_p
        ave_p.x(ave_p.x() + fit.point.x());             // then adds x AGAIN — corrupt
    }
    // Actually simpler: just sum components
}
// Redo with clean accumulation
ave_p = Point(0, 0, 0); num_p = 0; n_seg = 0; side_total_length = 0;
for (SegmentPtr sg1 : shower_segs) {                    // second loop — correct
    ...
}
```

The first loop double-counted the x component (`ave_p.x` gets `fit.point.x` twice), left `num_p`, `n_seg`, `side_total_length` all at zero, and was immediately undone by the reset at line 3229. **It was a leftover from an incomplete refactor** — dead code and wasted CPU.

The prototype (`NeutrinoID_nue_tagger.h:3305-3322`) has a single clean loop.

**Fix:** Removed the dead first loop; kept only the correct "redo" loop.

---

### 3. Minor efficiency — `ti.br4_1_n_vtx_segs` re-computed via lambda

**Location:** `bad_reconstruction_3` L2954-2957 (pre-fix)

**Problem:**

```cpp
// n_vtx_segs already computed at L2938:
size_t n_vtx_segs = start_vtx && start_vtx->descriptor_valid()
                    ? boost::out_degree(start_vtx->get_descriptor(), ctx.graph) : 0;
if (n_vtx_segs == 1 && ...) flag_bad1 = false;

// Then re-computed unnecessarily via lambda:
ti.br4_1_n_vtx_segs = (start_sg ? [&](){
    VertexPtr sv = shower->get_start_vertex_and_type().first;  // re-call
    return sv && sv->descriptor_valid() ? (int)boost::out_degree(sv->get_descriptor(), ctx.graph) : 0;
}() : 0);
```

The lambda re-calls `shower->get_start_vertex_and_type()` and re-queries `boost::out_degree` for a value already held in `n_vtx_segs`.

**Fix:** Replaced with `ti.br4_1_n_vtx_segs = (int)n_vtx_segs;`

---

## `bad_reconstruction_1`

### Structure

Checks whether the shower stem direction is inconsistent with the shower's PCA axis — a shower that is actually a mis-classified track. Fills `br2_*` fields.

Four angle variables:
- `angle` — PCA axis vs stem direction (at vertex end)
- `angle1` — shower 30 cm direction vs drift (proxy for near-horizontal showers)
- `angle2` — PCA axis vs drift (same proxy but global)
- `angle3` — global shower direction vs stem direction

Three cut blocks: energy-gated misalignment, unconditional large misalignment, and long shower-trajectory stem cut.

### Findings

#### ✅ Logic fidelity — cut conditions

All three blocks match prototype exactly:

| Cut block | Prototype | Toolkit |
|---|---|---|
| High-energy gate (>1 GeV) | no cut | no cut |
| 500–1000 MeV: `(angle1>10 \|\| angle2>10) && angle>30 && angle3>3` | ✅ | ✅ |
| <500 MeV: `(angle>25 && n_main_segs>1 \|\| angle>30) && (angle1>7.5 \|\| angle2>7.5)` | ✅ | ✅ |
| Unconditional: `angle>40 && (angle1>7.5 \|\| angle2>7.5) && max_angle<100` | ✅ | ✅ |
| Long traj stem: `angle>20 && ... && length>21cm && E<600MeV && kShowerTrajectory` | ✅ | ✅ |

#### ✅ Logic fidelity — vertex end determination

Prototype: `vertex->get_wcpt().index == sg->get_wcpt_vec().front().index`  
Toolkit: distance comparison `ray_length({vtx_fit_pt(ctx.main_vertex), sg_fits.front().point}) <= ...`

Distance comparison is more robust and avoids relying on wire-point index coincidence. **Improvement.**

#### ✅ Logic fidelity — PCA computation

Prototype: `main_cluster->Calc_PCA(tmp_pts)` then `get_PCA_axis(0)`.  
Toolkit: `ctx.self.calc_PCA_main_axis(tmp_pts).second` — same PCA, with an added `!tmp_pts.empty()` guard. **Improvement.**

#### ✅ Logic fidelity — `max_angle` loop

Prototype uses `map_vertex_segments[other_vertex]` (raw pointer map, all segments at vertex).  
Toolkit uses `boost::out_edges(ov_vd, ctx.graph)` — same set of segments, deterministic BGL order.

#### ⚠ Redundant `get_kine_best()` fallback | **Fixed** (see above)

---

## `bad_reconstruction_3`

### Structure

Two sub-checks (`flag_bad1`, `flag_bad2`), corresponding to `br4_1` and `br4_2` TaggerInfo blocks.

- **br4_1**: main-cluster fraction is small while off-cluster segments are far from the main-cluster tip (shower crosses clusters).
- **br4_2**: angular distribution of shower hits/vertices around the shower direction is too narrow (track-like fan).

### Findings

#### ✅ Logic fidelity — br4_1 thresholds

All six fraction/distance cuts match the prototype exactly (`0.40/40cm`, `0.25/33cm`, `0.16/23cm`, `0.10/18cm`, `0.05/8cm`, `8cm/0.1`). All three exception cases (`kAvoidMuonCheck`, `n_vtx_segs==1`, `n_main_segs>=4`) match. ✅

#### ✅ Logic fidelity — br4_1 `acc_close_length` refinement

Prototype (L2972):
```cpp
if (acc_close_length > 10*units::cm || num_close >=3 && acc_close_length > 4.5*units::cm || shower->get_start_segment()->get_flag_avoid_muon_check()) min_dis = min_dis1;
```

Toolkit (L2916-2919):
```cpp
if (acc_close_length > 10*units::cm ||
    (num_close >= 3 && acc_close_length > 4.5*units::cm) ||
    (start_sg && start_sg->flags_any(SegmentFlags::kAvoidMuonCheck)))
    min_dis = min_dis1;
```

The prototype has an operator precedence ambiguity (`&&` before `||` without parentheses). The toolkit is correctly parenthesized to the same meaning. ✅

#### ✅ Logic fidelity — br4_2 angular distribution

Prototype iterates `map_seg_vtxs` (interior points) and then `it->second` (endpoint vertices).  
Toolkit iterates `shower_segs` (interior fit points) then `find_vertices(ctx.graph, sg1)` (endpoint vertices).  
Both cover the same point set. ✅

#### ✅ Determinism — `IndexedSegmentSet` / `IndexedVertexSet`

Prototype uses `Map_Proto_Segment_Vertices` and `Map_Proto_Vertex_Segments` (raw pointer keyed maps, non-deterministic). Toolkit uses `shower->fill_sets(shower_vtxs, shower_segs, ...)` with indexed sets. **Improvement.**

#### ⚠ Redundant `get_kine_best()` fallback | **Fixed** (see above)
#### ⚠ `ti.br4_1_n_vtx_segs` lambda re-computing `n_vtx_segs` | **Fixed** (see above)

---

## `bad_reconstruction_2`

### Structure

Eight sub-checks (`br3_1` … `br3_8`):

| Sub-check | Description |
|---|---|
| br3_1 | Low-energy straight single-segment shower (track masquerade) |
| br3_2 | Segment type composition in main cluster looks track-like |
| br3_3 | Per-segment backward direction relative to stem |
| br3_4 | Accumulated backward length fraction |
| br3_5 | Average position of side segments vs stem far-end (Michel topology) |
| br3_6 | Anti-parallel segments at far end of stem |
| br3_7 | Short stem relative to main-cluster length |
| br3_8 | Sliding-window dQ/dx peak across main-cluster segments |

### Findings

#### ✅ Logic fidelity — br3_1 (four conditions)

All four low-energy straight-shower conditions match the prototype thresholds (100 MeV, 200 MeV, 0.95/0.85/0.925 ratios, 25 cm, topology/trajectory flags). ✅

#### ✅ Logic fidelity — br3_2 segment type classification

The `n_ele` / `n_other` classification (`kShowerTopology`, `kShowerTrajectory && med<1.3`, `ratio<0.92`) and the cut conditions (`n_ele==0 && n_other>0`, `n_ele==1 && n_ele<n_other && n_other<=2`, `n_ele==1 && n_other==0 && !other_fid`) match the prototype. ✅

#### ✅ Logic fidelity — br3_3/br3_4 variable naming improvement

Prototype reuses `total_length` as a loop accumulator, shadowing its earlier value:
```cpp
double total_length = shower->get_total_length(); // outer value
...
total_length = 0;                                  // overwritten for br3_3
for (...) total_length += length;                  // now = total main cluster length
```

Toolkit uses a separate `total_main_len2` for the br3_3 accumulation, keeping `total_length` unmodified. **Clarity improvement; functionally equivalent.**

#### ✅ Logic fidelity — br3_5 cut conditions

After removing the dead loop, the remaining logic exactly matches the prototype:
- `dir1.magnitude() > 3 cm || side_total_length > 6 cm`
- `(!avoid_check || n_seg > 1)`
- `dir_stem.angle(dir1) > 60°`
- `length > 10 cm && Eshower < 250 MeV`
- Exception: `num_main_segs + 6 < num_segs && main_length < 0.7*total_length` ✅

The prototype also computes `angle1` (max angle to drift from dir/dir1) but never uses it in any cut — only in debug print. The toolkit omits this unused variable. ✅ **Correct omission.**

#### ✅ Logic fidelity — br3_6 cut

```
angle > 150° && angle1 > 10° && !kShowerTrajectory && dir/len > 0.9
    && len > 7.5 cm && n_other_vtx_segs <= 4 && Eshower < 600 MeV
```
Matches prototype. ✅

#### ✅ Logic fidelity — br3_7 and br3_8

br3_7: `Eshower<200 && min_angle>60 && length < 0.2*shower_main_len` ✅  
br3_8: sliding window of 5 fit points, `max_dQ_dx > 1.85 && Eshower<150 && n_main_segs<=2 && main_frac>0.8` ✅

#### ⚠ Dead accumulation loop in br3_5 | **Fixed** (see above)
#### ⚠ Redundant `get_kine_best()` fallback | **Fixed** (see above)

#### ✅ Determinism — `IndexedSegmentSet` throughout

All segment/vertex iteration uses `IndexedSegmentSet` / `IndexedVertexSet`. **Improvement over prototype's raw pointer maps.**

---

## Summary Table

| Check | Result | Action |
|---|---|---|
| `get_kine_best()` fallback — `bad_reconstruction_1` | ⚠ Dead code | **Fixed** |
| `get_kine_best()` fallback — `bad_reconstruction_3` | ⚠ Dead code | **Fixed** |
| `get_kine_best()` fallback — `bad_reconstruction_2` | ⚠ Dead code | **Fixed** |
| Dead br3_5 accumulation loop | ⚠ Wrong + wasted CPU | **Fixed** |
| `ti.br4_1_n_vtx_segs` lambda re-computation | ⚠ Wasted CPU | **Fixed** |
| `bad_reconstruction_1` logic fidelity | ✅ Equivalent | — |
| `bad_reconstruction_1` vertex-end distance comparison | ✅ Improvement | — |
| `bad_reconstruction_1` PCA guard | ✅ Improvement | — |
| `bad_reconstruction_3` br4_1 thresholds and exceptions | ✅ Equivalent | — |
| `bad_reconstruction_3` operator precedence fix | ✅ Improvement | — |
| `bad_reconstruction_3` br4_2 angular counts | ✅ Equivalent | — |
| `bad_reconstruction_2` br3_1…br3_8 logic fidelity | ✅ Equivalent | — |
| `bad_reconstruction_2` `total_length` shadowing | ✅ Named `total_main_len2` (improvement) | — |
| `bad_reconstruction_2` unused `angle1` in br3_5 | ✅ Correctly omitted | — |
| `IndexedSegmentSet` / `IndexedVertexSet` throughout | ✅ Improvement over prototype | — |
| BGL `out_edges` (deterministic) vs raw `map_vertex_segments` | ✅ Improvement | — |

---

## Changes Made

**File:** `clus/src/NeutrinoTaggerNuE.cxx`

1. Removed redundant `(shower->get_kine_best() != 0) ? ... : shower->get_kine_charge()` from `bad_reconstruction_1` — replaced with `shower->get_kine_best()`.
2. Same fix in `bad_reconstruction_3`.
3. Same fix in `bad_reconstruction_2`.
4. Removed dead first accumulation loop in `bad_reconstruction_2` br3_5 (was L3218–3226 pre-fix) — a leftover from an incomplete refactor that corrupted `ave_p` before the correct "redo" loop reset it.
5. Replaced `ti.br4_1_n_vtx_segs` lambda (which re-called `get_start_vertex_and_type()` and `boost::out_degree`) with the already-computed `(int)n_vtx_segs`.
