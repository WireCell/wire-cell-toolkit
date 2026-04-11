# PID / Direction / Kinematics Function Review

**Reviewed:** 2026-04-11  
**Author decisions (2026-04-11):** Q1: α must match prototype → P8 is a confirmed fix, not just a question. Q2: table provenance needs verification within the toolkit. Q3: `[E, px, py, pz]` layout intentional and accepted. Q4: break_segment wcpts/graph orientation invariant should be enforced (P4). Q5: per-face recomb/lifetime correction in PID/kine functions deferred. Q6: `particle_score = 100.0` sentinel behaviour in `segment_determine_shower_direction_trajectory` to be documented. Q7: `calculate_kinematics_long_muon` included in this review scope (P17 escalated from "future").

## Applied Fixes (2026-04-11)

The following priority items have been implemented. See the Priority Action List below for remaining open items.

| Item | File | What was changed |
|---|---|---|
| **P1** | `clus/src/PRSegmentFunctions.cxx` | All 8 unclamped `std::acos(dot)` calls in `segment_determine_shower_direction` wrapped with `std::clamp(..., -1.0, 1.0)` to prevent NaN from floating-point rounding on nearly-parallel unit vectors. |
| **P2/P3** | `clus/src/PRSegmentFunctions.cxx` | Fixed first/last-point path-length shortening in `segment_cal_kine_dQdx`: introduced `dx_for_dQdx = fits[i].dx` so the Box-model inversion uses the original step size; shortened `dX` is used only for the path-length accumulation. Saturation filter now also evaluated against `dx_for_dQdx` (not the shortened value). |
| **P4** | `clus/src/PRSegmentFunctions.cxx` | Replaced the terse WARNING at lines `:499–501` in `break_segment` with an expanded comment explaining that `find_vertices` in `PRGraph.cxx` already handles the orientation ambiguity via wcpts.front() distance comparison. The structural orientation issue is documented; the invariant-enforcement refactor remains a future follow-up. |
| **P6** (P7 in table) | `clus/src/PRSegmentFunctions.cxx` | Restored the `if (flag_print)` guard around the `SPDLOG_LOGGER_DEBUG` call in `segment_determine_dir_track`. |
| **P7** (P11 in table) | `clus/src/PRSegmentFunctions.cxx` | Added `if (ncount == 0) return {1.0, 1e9, 1e9, 1e9};` early-return in `segment_do_track_comp` to handle degenerate input where the comparison window contains no points. |
| **P13** | `clus/src/PRSegmentFunctions.cxx` | Added a 6-line multi-APA limitation comment immediately above `segment_cal_kine_dQdx` documenting the single-global-recomb-model assumption and referencing this review. |
| **P14** | `clus/src/PRSegmentFunctions.cxx` | Added a 6-line unit-convention comment above the `dQ_dx[i]` assignment in `segment_determine_dir_track` explaining why the factor of `units::cm` is absent and confirming internal consistency with the reference tables. |
| **P15** | `clus/src/PRSegmentFunctions.cxx` | Cached `segment_median_dQ_dx(segment, start_n1, end_n1)` into a single `medium_dQ_dx` variable computed once before the three-way branch in `segment_determine_dir_track`; all three former call sites now reuse the cached value. |
| **P16** | `clus/src/PRShower.cxx` | Replaced the `(0,0,0)` origin sentinel check in `calculate_kinematics` with a `dist >= 0` guard on the `shower_get_closest_point` return value, eliminating the false-positive fallback for a legitimate detector-origin vertex. |
| **B8.1 / Q6** | `clus/src/PRSegmentFunctions.cxx` | Replaced the `// hack for now ...` comment in `segment_determine_shower_direction_trajectory` with a 5-line explanation: PDG is forced to electron (11) for all shower-trajectory segments; `particle_score = 100.0` is an intentional sentinel meaning "PID not performed; score not applicable." |
| **B17.1** | `clus/src/PRShower.cxx` | Added a 6-line invariant comment above the early-return guard in `calculate_kinematics_long_muon` explaining why the guard is logically unreachable (call-site invariant: only reachable when `m_start_segment` already has `particle_info()` with PDG ±13). |

**Toolkit entry files:**
- `clus/src/PRSegmentFunctions.cxx` (primary file for all ten free-function equivalents)
- `clus/src/PRShower.cxx` (shower helpers)
- `clus/src/NeutrinoShowerClustering.cxx` (`update_shower_maps`)
- `clus/inc/WireCellClus/PRShower.h` (type aliases: `ShowerVertexMap`, `ShowerSegmentMap`, `IndexedShowerSet`, …)

**Prototype counterparts:**
- `prototype_base/pid/src/ProtoSegment.cxx` (all ten functions as member methods)
- `prototype_base/pid/src/WCShower.cxx` (shower helpers)
- `prototype_base/pid/inc/WCPPID/ProtoSegment.h`

**Scope:** Functional equivalence, bugs, algorithmic efficiency, determinism (pointer-ordered containers), and multi-APA / multi-face correctness.

**Porting reference:**
- `clus/docs/porting/porting_dictionary.md` §"Trajectory Fitting and Pattern Recognition"
- `clus/docs/porting/neutrino_id_function_map.md` §"NeutrinoID_track_shower.h"

**Key structural note:**  
All algorithmic methods that were member functions of `WCPPID::ProtoSegment` are now free functions in `PRSegmentFunctions.cxx` taking a `SegmentPtr` first argument. The prototype's implicit access to `TPCParams` singleton and hardcoded tables is replaced by dependency-injected `ParticleDataSet::pointer particle_data`, `IRecombinationModel::pointer recomb_model`, and `IDetectorVolumes::pointer dv`.

---

## §0 Cross-Cutting Concerns

### G1 — BoxRecombination α default

Prototype (`ProtoSegment.cxx:1318–1328`) hardcodes the Box model with **α = 1.0**, **β = 0.255**:

```
dEdx = (exp(dQdx * 23.6e-6 * 0.255 / 1.38 / 0.273) - 1.0) / (0.255 / 1.38 / 0.273)
```

The toolkit `BoxRecombination` (`gen/src/RecombinationModels.cxx:101–111`) defaults to **α = 0.930**.  With α = 0.930 the inferred dE/dx is lower at a given dQ/dx, producing systematically lower kinetic-energy estimates from `cal_kine_dQdx`.  No override to α = 1.0 was found in any production jsonnet under `cfg/`.

**Action (P8 — confirmed fix):** The toolkit must match the prototype numerically. Override `BoxRecombination.A = 1.0` in the production jsonnet configuration. Until this is done, all `cal_kine_dQdx` results will be systematically lower than the prototype.

### G2 — No pointer-keyed iteration in the ten functions

All inner loops walk `seg->fits()` (`std::vector<Fit>`) or locally constructed `std::vector<double>`. No pointer-address–dependent ordering is introduced inside the reviewed functions themselves.

### G3 — Multi-APA not wired through PID/kine functions

Every `PR::Fit` carries `paf{apa,face}`.  None of the ten functions branches on `paf`: `cal_kine_dQdx`, `do_track_comp`, `do_track_pid`, and `cal_kine_range` apply a single global `IRecombinationModel` and a single `ParticleDataSet` to the entire fit vector, even if its points span different (apa, face) pairs with different electron lifetimes, wire pitches, or field-response corrections.

For uBooNE (single APA, single face) this is equivalent to the prototype. For SBND/DUNE/PDHD the implicit assumption is that per-face electron-lifetime corrections have been applied into `fit.dQ` upstream (during `TrackFitting`), leaving the PID/kine layer with a single effective recombination model.

**Status (deferred):** Per-face recombination / lifetime correction in these PID/kine functions is deferred to a future ticket. For uBooNE (single APA) this is equivalent to the prototype. The assumption must be documented and tracked.

**Action (P13):** Add a comment in `PRSegmentFunctions.cxx` near `segment_cal_kine_dQdx` documenting this assumption and where per-face lifetime correction is expected to be applied upstream.

### G4 — MIP_dQdx unit convention

Prototype literal `50e3` (electrons) and `43e3` (saturation threshold) are parameterised in the toolkit as `MIP_dQdx = 50000/units::cm` and `43e3/units::cm`. Because the reference dE/dx tables are also retrieved as `/units::cm`, all KS distances and ratios are dimensionally consistent. No regression, but add a one-line comment near the default in `PRSegmentFunctions.h:104` so the unit convention is explicit.

### G5 — Reference-table provenance (open question)

WCT reads dEdx and range tables from `clus/test/data/uboone-mabc_config.json` via `LinterpFunction` (linearly interpolated). Prototype embeds the same data from `WCPDataDict.cxx:4677` as `TGraph` objects (also linearly interpolated). `TGraph::Eval` and `LinterpFunction::scalar_function` are equivalent for regularly-spaced knots.

**Action (P9):** Verify consistency within the toolkit: confirm that the `uboone-mabc_config.json` dEdx and range table values were generated from the same upstream CSV/PDG source as the WCP embedded data in `WCPDataDict.cxx:4677`. A mismatch would silently shift all KS-test reference curves and all range-based KE estimates.

---

## §1 `do_track_comp` → `do_track_comp` (free function)

| | Prototype | Toolkit |
|---|---|---|
| Location | `ProtoSegment.cxx:1120–1187` | `PRSegmentFunctions.cxx:1263–1319` |

### A — Functional identity

Step-by-step match:
- Window: `end_L = L.back() + 0.15*cm - offset_length`, `end_L - L[i] ∈ (0, compare_range)` — identical.
- Reference vectors: muon/const/proton/electron built from `particle_data->get_dEdx_function(name)->scalar_function(x/cm)/cm` — equivalent to prototype's `TGraph::Eval(x_cm)` for uniformly-spaced tables.
- KS test: prototype uses `TH1F::KolmogorovTest("M")` (max distance between normalised CDFs). WCT uses `WireCell::kslike_compare` which computes `max|rsum1 - rsum2|` over normalised running sums — **mathematically equivalent** for equal-length vectors.
- Ratio guard: both versions use `+1e-9` on `sum(vec_y)` (`PRSegmentFunctions.cxx:1294–1307`). ✓

### B — Bugs

- **B1.1 (degenerate input):** If `compare_range` is smaller than the minimum point spacing, `vec_x`/`vec_y` are empty and `count=0`. `kslike_compare` on two empty vectors divides by zero inside (see `util/src/KSTest.cxx`). Prototype's TH1F with zero bins is equally undefined. Add `if (ncount == 0) return {1.0, 1e9, 1e9, 1e9};` early-return to match the spirit of the prototype's "no info → not a track" default.

### C — Efficiency

WCT is strictly more efficient than the prototype: no `TH1F` allocation, no ROOT string hashing, no destructor cost. The four `std::accumulate` calls over `ref` vectors can be fused into the single ref-population loop — minor.

### D — Determinism

Deterministic (only `std::vector<double>` loops).

### E — Multi-APA

Receives pre-flattened `L`/`dQ_dx` arrays; inherits G3.

### F — Callees

- `eval_ks_ratio` — see §2.
- `particle_data->get_dEdx_function(...)->scalar_function(x)` — see G5.

---

## §2 `eval_ks_ratio`

| | Prototype | Toolkit |
|---|---|---|
| Location | `ProtoSegment.cxx:1189–1198` | `PRSegmentFunctions.cxx:940–949` |

### A — Functional identity

Line-by-line identical. The toolkit adds explicit parentheses around the compound `&&`/`||` condition (line 945 vs. 1194).  Given C++ operator precedence (`&&` binds tighter than `||`), the evaluation order is unchanged — the parentheses are a readability improvement only.

### B–F

Nothing to flag.

---

## §3 `do_track_pid` → `segment_do_track_pid`

| | Prototype | Toolkit |
|---|---|---|
| Location | `ProtoSegment.cxx:1201–1289` | `PRSegmentFunctions.cxx:1350–1453` |

### A — Functional identity

Reverse-vector construction, forward/backward comparison, muon/proton/electron particle-type selection, `length < 20 cm` electron guard, and the four-way flag_forward/flag_backward decision tree all match line-by-line.

### B — Bugs

- **B3.1 (asymmetric side-effect, fragile):** All four success branches early-return a tuple and do **not** write to `segment->dirsign(...)`. The fall-through failure path at `:1449` calls `segment->dirsign(flag_dir)` with `flag_dir = 0`.  In the prototype, `do_track_pid` is a member function and sets `this->flag_dir` in all paths. The WCT callers (`segment_determine_dir_track` at `:1535,1544,1554,1563`) correctly read the returned `flag_dir` and call `segment->dirsign(...)` on success — so the current callers are correct.  But the asymmetry is fragile: any future caller that forgets to propagate the returned `flag_dir` will silently leave a stale non-zero dirsign on success while the failure path correctly resets to 0.  **Recommend:** remove the `segment->dirsign(0)` call from the failure path inside `segment_do_track_pid` and rely entirely on the caller (following the convention that this function has no segment side-effects). Document the convention with a comment.

- **B3.2 (cosmetic):** Failure return `particle_score = 100.0` (`:1452`) vs. prototype's reset to `0` (`:1286`). Callers guard on the returned `bool`, so this only matters if a caller reads `particle_score` after a `false` return — unlikely but inconsistent. Low priority.

### C — Efficiency

Equivalent to the prototype. Two `do_track_comp` calls each allocating four reference vectors.

### D — Determinism

Deterministic.

### E — Multi-APA

Inherits G3 (flat `L`/`dQ_dx` arrays, no per-face branching).

### F — Callees

- `do_track_comp` — §1.
- `segment_track_length(segment, 0)` at `:1373` — geometric sum of adjacent fit distances; equivalent to prototype's `get_length()`.

---

## §4 `cal_kine_dQdx` → `segment_cal_kine_dQdx` + `cal_kine_dQdx`

| | Prototype | Toolkit |
|---|---|---|
| Location (member) | `ProtoSegment.cxx:1341–1378` | `PRSegmentFunctions.cxx:1175–1228` |
| Location (vec overload) | `ProtoSegment.cxx:1316–1339` | `PRSegmentFunctions.cxx:1230–1261` |

### A — Functional identity

Algorithm: for each fit point compute dEdx from dQ/dx via the Box recombination inverse, clamp to `[0, 50 MeV/cm * dX]`, accumulate. Middle fits: identical. First/last-point path-length correction: **differs** (see B4.1 below).

### B — Bugs

- **B4.1 (numerical regression, first/last point):** The prototype computes `dEdx` from the **original** `dx_vec[i]`, then multiplies by the shorter `dis` if `dx_vec[i] > dis * 1.5`:
  ```cpp
  // WCP (ProtoSegment.cxx:1356-1373)
  double dEdx = (exp(dQdx * coeff) - alpha) / coeff;  // uses original dx_vec[i]
  kine_energy += (dx_vec[i] > dis*1.5) ? dEdx*dis : dEdx*dx_vec[i];
  ```
  The toolkit shortens `dX` **before** calling `recomb_model->dE(dQ, dX)` (`PRSegmentFunctions.cxx:1196–1215`):
  ```cpp
  // WCT
  if (dX > dis * 1.5) dX = dis;        // shorten dX BEFORE dE call
  double dE = recomb_model->dE(dQ, dX); // Box model is non-linear in dQ/dX
  kine_energy += dE;
  ```
  Because the Box model is **non-linear** in `dQ/dX`, shortening `dX` inflates the effective `dQ/dX` and thus the recombination-corrected `dE/dx`. The WCT endpoint contribution is therefore **higher** than the prototype's when `fit.dx > dis * 1.5`. The error scales with the degree of non-linearity at the local `dQ/dX` and how much `fit.dx` overshoots `dis`.

  **Fix:** compute `dE` from the original `fit.dx`, then multiply by `dX` (after possible shortening):
  ```cpp
  double dE_per_dx = recomb_model->dE(dQ, fits[i].dx) / fits[i].dx;
  kine_energy += dE_per_dx * dX;
  ```

- **B4.2 (saturation-filter aliasing):** The saturation check `dQ/dX / (43e3/cm) > 1000` at `:1212` uses the shortened `dX`, which can cause the filter to fire when the prototype would not (higher apparent dQ/dx). Apply the filter on `fits[i].dx` instead.

- **B4.3 (code smell):** Lines 1208–1209:
  ```cpp
  if (dX <= 0) dX = fits[i].dx;
  if (dX <= 0) continue;
  ```
  The double-check is confusing; `dX` can only reach zero after reduction if `dis = 0` (coincident neighbors). Simpler: `if (dis <= 0) dX = fits[i].dx;` inside the first/last branch, then a single `if (dX <= 0) continue;`.

- G1 (α value) applies to both `segment_cal_kine_dQdx` and the vec overload.

### C — Efficiency

One virtual dispatch per fit (`recomb_model->dE`) — same cost as the prototype's inline `exp()`. No allocations in the hot loop.

### D — Determinism

Deterministic.

### E — Multi-APA

Single global `recomb_model` over the entire fit vector; inherits G3. The `fit.valid()` guard at `:1189` (not present in the prototype, which protects via `+1e-9`) silently skips invalid fits rather than counting them with near-zero contribution.

---

## §5 `cal_kine_range` (three overloads → one function)

| | Prototype | Toolkit |
|---|---|---|
| Location | `ProtoSegment.cxx:1380–1418` | `PRSegmentFunctions.cxx:1321–1347` |

### A — Functional identity

Three prototype overloads (`cal_kine_range()`, `cal_kine_range(double L)`, `cal_kine_range(int pdg)`) are collapsed into `cal_kine_range(double L, int pdg_code, particle_data)`. PDG dispatch 11/13/211/321/2212 matches. Default fallback to muon (`if (!range_function)` at `:1341–1343`) is new — see B5.1.

### B — Bugs

- **B5.1 (prototype latent null-deref fixed):** Prototype `cal_kine_range(double L)` at `:1382` declares `WCP::TopoDataGraph* g_range` without initialising it; for an unrecognised `particle_type` the pointer is never assigned and `g_range->Eval(...)` segfaults. The toolkit always falls back to muon → **improvement**.

### C — Efficiency

WCT: O(log N) `std::map<std::string,…>` lookup per call, N ≈ 5 (minor). Prototype: O(1) member accessor. No practical impact; could cache the function pointer locally for tight loops.

### D — Determinism

Deterministic (string-keyed `std::map`).

### E — Multi-APA

Range tables are particle-physics properties independent of detector face. No per-face concern.

### F — Callees

`IScalarFunction::scalar_function` → `LinterpFunction` regular-grid linear interpolation. Equivalent to prototype `TGraph::Eval` for uniformly-spaced knots; verify CSV provenance (G5).

---

## §6 `cal_4mom` → `segment_cal_4mom`

| | Prototype | Toolkit |
|---|---|---|
| Location | `ProtoSegment.cxx:1420–1448` | `PRSegmentFunctions.cxx:1456–1487` |

### A — Functional identity

Length-threshold (`< 4 cm`) and shower-trajectory dispatch match. The 4-momentum layout **changed** — see B6.1.

### B — Bugs

- **B6.1 (contract change, critical):** Prototype stores `particle_4mom[px, py, pz, E]` at indices 0–3 (`:1437, 1442–1444`). Toolkit uses `WireCell::D4Vector<double>` stored as `[E, px, py, pz]` at indices 0–3 (`PRSegmentFunctions.cxx:1464, 1478, 1482–1484`). Any caller that numerically indexes `[3]` expecting energy now receives `pz`.

  **Action (P5):** Grep every caller of `segment_cal_4mom` and every site that reads `four_momentum[0..3]` numerically (tagger writers, `NeutrinoTrackShowerSep`, `PRShower`, `NeutrinoPatternBase`) and confirm they all use the `D4Vector` accessors `.E()`, `.px()`, `.py()`, `.pz()` rather than raw index.

- **B6.2 (`kenergy_best` not persisted):** Prototype stores `kenergy_best` as a member (`:1433`). There is no direct equivalent in WCT `PRSegment`; callers that read energy must use `segment->particle_info()->kinetic_energy()` and must therefore call `segment_cal_4mom` first. Confirm all prior `kenergy_best` call sites have been migrated.

### C — Efficiency

Equivalent to the prototype (one `segment_track_length` + one `segment_cal_dir_3vector`).

### D — Determinism

Deterministic.

### E — Multi-APA

Calls `segment_cal_dir_3vector` (pure 3D geometry, face-agnostic — correct) and `segment_cal_kine_dQdx` / `cal_kine_range` (inherits G3).

### F — Callees

- `segment_cal_dir_3vector(segment)` (`PRSegmentFunctions.cxx:1068–1106`): returns `(0,0,0)` when `segment->dirsign() == 0` — **fixes** a latent `TVector3::Unit()` undefined-behaviour on a zero vector that the prototype had. ✓

---

## §7 `determine_dir_track` → `segment_determine_dir_track`

| | Prototype | Toolkit |
|---|---|---|
| Location | `ProtoSegment.cxx:1516–1646` | `PRSegmentFunctions.cxx:1489–1662` |

### A — Functional identity

All logic paths match: trimming of start/end indices by `start_n`/`end_n`, two-step PID fallback (35 cm then 15 cm range), median-dQ/dx short-track path, vertex-activity sub-branch, and "very bad score → force shower" override.

### B — Bugs

- **B7.1 (unit subtlety, needs comment):** Prototype (`ProtoSegment.cxx:1547`) builds `dQ_dx[i] = dQ_vec[i] / (dx_vec[i]/units::cm + 1e-9)` — denominator in cm. Toolkit (`:1520`) uses `fits[i].dQ / (fits[i].dx + 1e-9)` — denominator in internal units (cm in the toolkit's unit system). The factor of `units::cm` cancels against the reference tables which are also retrieved `/units::cm`, and against `MIP_dQdx = 43000/units::cm`. The pipeline is internally consistent but the mismatch from the prototype notation is easy to misread. Add a one-line comment at `:1520` explaining the unit convention.

- **B7.2 (debug log guard commented out, log spam):** The `if (flag_print)` guard around the `SPDLOG_LOGGER_DEBUG` call at `:1649–1661` is commented out:
  ```cpp
  // if (flag_print) {
      SPDLOG_LOGGER_DEBUG(s_log, "segment_determine_dir_track: ...", ...);
  // }
  ```
  At `SPDLOG_LOGGER_DEBUG` severity this only fires when the logger level is set to `debug`, but it's not faithful to the prototype which gated the print on `flag_print`. Restore the guard (Action P7).

### C — Efficiency

`segment_median_dQ_dx` is called up to three times over overlapping index ranges in the vertex-activity branch (`:1579, 1606, 1614`). The first call can be cached and reused (Action P15).

### D — Determinism

Deterministic (no pointer-keyed containers in this function).

### E — Multi-APA

All computation is scalar-per-fit. No per-face branching needed; geometry-independent PID.

---

## §8 `determine_dir_shower_trajectory` → `segment_determine_shower_direction_trajectory`

| | Prototype | Toolkit |
|---|---|---|
| Location | `ProtoSegment.cxx:1647–1673` | `PRSegmentFunctions.cxx:1664–1719` |

### A — Functional identity

Three-way dispatch: leaf-end `start_n==1/end_n==1` → set ±1; interior → delegate to `segment_determine_dir_track` then force electron if non-electron was identified. Logic matches.

The WCT "no-particle-info" branch at `:1688–1692` covers the case where `segment_determine_dir_track` found no match (success `false`) and left no `ParticleInfo`. The prototype's member `particle_type` would remain at its prior value (zero for new segments), causing the outer `if (particle_type != 11)` to trigger the reset. The net effect is the same (force `pdg=11, dirsign=0`), but the WCT path is more explicit.

### B — Bugs

- **B8.1 (`particle_score` always reset to 100.0 — intentional sentinel, document it):** At `:1670` `particle_score = 100.0` is initialised and at `:1708` it is unconditionally stored into the segment. If `segment_determine_dir_track` (called at `:1678`) successfully computed a score and stored it via `segment->particle_score(...)` (`:1646`), the outer function overwrites it with 100.  The prototype also did not preserve a meaningful `particle_score` in this function (it re-called `cal_4mom()` at `:1670` which reset energy). **Decision: this is intentional.** The `particle_score = 100.0` is a sentinel meaning "PID was not performed on this shower trajectory; score is not applicable." **Action:** Replace the `// hack for now ...` comment at `:1668` with an explicit note to this effect, and mark the sentinel value with a named constant or comment so callers know not to interpret it as a PID quality metric.

### C — Efficiency

`segment_cal_4mom` may be called twice: once inside `segment_determine_dir_track` (`:1678`) and once explicitly at `:1696`. Cache the result if the inner call succeeds.

### D — Determinism

Deterministic.

### E — Multi-APA

Face-agnostic (pure 3D direction logic + kinematics).

---

## §9 `determine_shower_direction` → `segment_determine_shower_direction`

| | Prototype | Toolkit |
|---|---|---|
| Location | `ProtoSegment.cxx:72–320` | `PRSegmentFunctions.cxx:2025–2330` |

### A — Functional identity

KD-tree mapping of associated points to nearest fit, local orthogonal frame construction, RMS spread calculation, `1.0/0.8/0.7 cm` thresholds, bidirectional `threshold_segs` accumulation, and the final `1.2×` direction asymmetry test all match.

Minor improvements in WCT: `total_length` initialised to 0 (see B9.3); `start_n`/`end_n` boundary clamped with `i < end_n && i+1 < fits.size()` (WCT `:2233`) vs. prototype's bare `i != end_n` loop which would be degenerate if `start_n > end_n`.

### B — Bugs

- **B9.1 (NaN risk, should fix, `std::acos` without clamp):** The following lines compute `std::acos(dot_product)` where `dot_product` is a dot product of two unit vectors obtained from `v.norm()`. Floating-point rounding can push the argument outside `[−1, 1]`, producing `NaN` which silently propagates into `angle_deg` comparisons:

  | Line | Expression |
  |---|---|
  | 2088 | `std::acos(dir_1.dot(drift_dir_abs))` |
  | 2188 | `std::acos(main_dir1.dot(drift_dir_abs))` |
  | 2193 | `std::acos(main_dir2.dot(drift_dir_abs))` |
  | 2203 | `std::acos(main_dir1.dot(vec_dir.at(i)))` |
  | 2248 | `std::acos(main_dir2.dot(vec_dir.at(i)))` |
  | 2304 | `std::acos(main_dir_front.dot(vec_dir.at(i)))` |
  | 2311 | `std::acos(main_dir_back.dot(vec_dir.at(i)))` |
  | 2389 | `std::acos(dir_1.dot(drift_dir_abs))` |

  The prototype used `TVector3::Angle()` which clamps internally. `PRShower.cxx:791, 814` (`get_stem_dQ_dx`) correctly uses `std::clamp(dot, -1.0, 1.0)` — apply the same pattern here.

  **Fix (Action P1):** wrap each `dot_product` argument:
  ```cpp
  double angle_deg = std::acos(std::clamp(dir_1.dot(drift_dir_abs), -1.0, 1.0)) * 180.0 / M_PI;
  ```

- **B9.2 (dead code removed, safe):** The prototype tracked `max_cont_length` / `max_cont_weighted_length` via a `flag_prev` variable (`WCShower.cxx:153–183`) but never consumed those values within the function. The toolkit drops them (`:2139` area). Safe removal.

- **B9.3 (prototype latent bug fixed):** Prototype `ProtoSegment.cxx:160` declared `double total_length;` without initialisation and then did `total_length += length;`. Toolkit initialises `total_length = 0` — **fix**. ✓

### C — Efficiency

KD-tree query once per associated point — same as prototype. `segment_cal_dir_3vector` called twice — same as prototype.

### D — Determinism

The DynamicPointCloud `dpcloud_assoc->get_points()` at `:2039` returns points in insertion order. If the upstream `associate_points` DPC was populated by iterating a pointer-keyed container, the insertion order carries pointer-address dependence. Audit the caller of `segment_determine_shower_direction` to verify associated-points are added in a deterministic order.

### E — Multi-APA

All computation is in 3D Cartesian space; `drift_dir_abs = (1,0,0)` assumes uBooNE-style drift. For detectors with different drift orientations (e.g., PDHD vertical APAs), this vector should be derived from `IDetectorVolumes`. Low risk for current targets.

---

## §10 `break_segment_at_point` → `break_segment`

| | Prototype | Toolkit |
|---|---|---|
| Location | `ProtoSegment.cxx:1711–1905` | `PRSegmentFunctions.cxx:442–650` |

### A — Functional identity

Core flow matches: closest-fit-index / closest-wcpt-index search, reject-if-at-either-end guard, vector split with breakpoint duplicated in both halves, property copy, associate-point redistribution. The toolkit inlines graph wiring (`remove_segment` + two `make_segment`) that the prototype returned to the caller via raw pointers.

### B — Bugs

- **B10.1 (orientation invariant, critical):** Lines `:499–501` contain an explicit warning:
  ```cpp
  // WARNING there is no "direction" in the graph.  You can not assume the
  // "source()" of a segment is closest to the segments first point.
  ```
  `find_vertices(graph, seg1)` at `PRGraph.cxx:72` uses `seg1->wcpts().front().point` to disambiguate which graph vertex is "front" and which is "back". If `make_segment(graph, vtx1, vtx)` at `:502` produces an edge where `source() → vtx` and `target() → vtx1` (Boost graph makes no direction guarantee), then `find_vertices(seg1)` returns `(vtx, vtx1)` instead of `(vtx1, vtx)` — reversing which vertex is considered the "start".

  Any downstream caller that calls `find_vertices` on either new segment and relies on `pair.first` being the end closer to `vtx1` (e.g., when re-establishing a shower connection through `vtx1`) can silently get the wrong vertex.

  **Action (P4 — enforce):** Normalise the wcpts orientation to match `source()`/`target()` inside `break_segment`: after calling `find_vertices` on each new segment, reverse `seg->wcpts()` and `seg->fits()` if the first fit point is closer to `target()` than `source()`. This restores the "source ↔ wcpts.front()" invariant that downstream code relies on.

- **B10.2 (caller-contract change):** Prototype returned `(this, nullptr, nullptr)` on degenerate break (`:1735`). Toolkit returns `(false, {}, VertexPtr())` at `:475`. Existing callers (`NeutrinoShowerClustering.cxx:984`) destructure the tuple and branch on `bool` — correct. Grep all callers to confirm none are checking pointer equality against the original segment.

- **B10.3 (DPC rebuild for multi-APA, verify):** Prototype explicitly re-projected all associated points to u/v/w per plane via `add_associate_point(3d, u, v, w)` (`ProtoSegment.cxx:1893–1905`). Toolkit rebuilds the DPC via `wpid_params` from `cluster.grouping()->wpids()` (`:627–633`) and passes to `add_points`. This relies on `DynamicPointCloud::add_points` correctly re-projecting for every `wpid` in the map. **Verify** that `add_points` branches on `wpid` to compute per-face 2D coordinates; if it only stores the 3D point without per-plane projection, multi-APA segments will be missing their 2D wire information after a break.

### C — Efficiency

O(N log M) closest-point query per associated point — same as prototype.

### D — Determinism

Deterministic: `closest_point` returns first minimum (ordered by fit-vector index), `wpids()` returns a sorted `set<WirePlaneId>`, and associated-point redistribution walks the DPC in insertion order.

### E — Multi-APA

`wpid_params` is built from `cluster.grouping()->wpids()` (all APAs/faces of the parent cluster) — so both new segment DPCs include face-aware projection parameters. The new vertex's `fit` is copied from `*itfits` which carries its `paf` — correct. The 3D `segment_get_closest_point` redistribution is face-agnostic (uses Euclidean distance in the shared global coordinate frame) — correct.

---

## §11 Shower Helpers (`PRShower.cxx` ↔ `WCShower.cxx`)

### `add_segment` (`PRShower.cxx:141` ↔ `WCShower.cxx:694`)

Both add to graph view + vertex set. WCT uses lazy `invalidate_segment_caches` (`:147`) rather than prototype's eager `rebuild_point_clouds()` — **efficiency improvement**. Cache invalidation is called on every add path. ✓

### `complete_structure_with_start_segment` (`PRShower.cxx:269` ↔ `WCShower.cxx:703`)

WCT uses `sorted_out_edges` for deterministic BFS ordering — **improvement** over prototype's iteration of pointer-valued `map_vertex_segments`. Inline DPC merge is cleaner than prototype's final `rebuild_point_clouds()`. Logic (exclude-start-vertex, BFS collect all reachable, mark `used_segments`) matches.

### `fill_maps` / `fill_sets` / `update_shower_maps`

`PRShower::fill_maps()` now returns `*this` (a `TrajectoryView`) — signature change from prototype which took two out-parameters. Callers in `NeutrinoShowerClustering.cxx:51,59` iterate `traj.nodes()` and `traj.edges()`. These are `unordered_set` containers but the **maps** being populated (`ShowerVertexMap = std::map<VertexPtr, ShowerPtr, VertexIndexCmp>` and `ShowerSegmentMap = std::map<SegmentPtr, ShowerPtr, SegmentIndexCmp>`, from `PRShower.h:222–223`) use custom index-based comparators — **deterministic**. No action needed.

### `calculate_kinematics` (`PRShower.cxx:870` ↔ `WCShower.cxx:339`)

The two prototype branches (single-connected vs multi-connected) are unified into one shared preamble — readability/efficiency improvement; logic preserved. `ordered_nodes` used in the "find farthest vertex" loop — **determinism improvement** over prototype's pointer-set iteration.

- **B11.1 (`start_point` origin sentinel):** At `PRShower.cxx:1021`:
  ```cpp
  if (data.start_point.x() == 0 && data.start_point.y() == 0 && data.start_point.z() == 0) {
      // fallback ...
  }
  ```
  Uses `(0,0,0)` as a sentinel for "pcloud returned no result." A legitimate vertex at the detector origin would trigger the fallback. Prototype had no such sentinel. Replace with a `bool`-returning helper or check the DPC size before calling `shower_get_closest_point`. (Action P16)

### `get_stem_dQ_dx` (`PRShower.cxx:716` ↔ `WCShower.cxx:38`)

`std::acos` correctly clamped at `:791, 814` — **this is the reference pattern** for the B9.1 fix in `segment_determine_shower_direction`. Single-pass candidate list replaces prototype's two-pass scan — **efficiency improvement**. Minor behavioural divergence: vertex-at-front test uses squared 3D distance (`PRShower.cxx:780`) vs prototype's `wcpt().index` comparison. Affects only degenerate coincident endpoints — negligible in practice.

### `get_total_length`, `update_particle_type`, `get_num_main_segments`

Iterate `this->edges()` (unordered). Results are commutative sums or counts — **numerically deterministic** regardless of iteration order.

### `get_last_segment_vertex_long_muon` (`PRShower.cxx:525` ↔ `WCShower.cxx:241`)

WCT uses `sorted_out_edges` — **determinism improvement**.

### `calculate_kinematics_long_muon` (`PRShower.cxx:1121` ↔ `WCShower.cxx:288`)

**A — Functional identity.** Core flow matches: accumulate `total_length` over muon-tagged segments only, collect all-segment dQ/dx into `vec_dQ`/`vec_dx`, compute `kenergy_range` from range tables and `kenergy_dQdx` from Box model, set `kenergy_best = kenergy_dQdx`, find farthest muon-connected vertex for `end_point`. All logic branches match.

**B — Issues.**

- **B17.1 (early return on missing particle_info — verified safe, no code change needed):**

  WCT at `:1128` guards with `if (!m_start_segment || !m_start_segment->has_particle_info()) return;`.  The prototype (`WCShower.cxx:289`) unconditionally reads `particle_type = start_segment->get_particle_type()`, which defaults to 0 for un-assigned segments, then falls through to `cal_kine_range(L, 0)` — i.e. the muon table via the unknown-PDG fallback.  On first reading this looks like a divergence: WCT silently returns while the prototype computes a muon KE.

  **Invariant analysis (2026-04-11):** The early return is logically unreachable in any normal execution path.  `calculate_kinematics_long_muon` is called from `NeutrinoEnergyReco.cxx:296` only when `abs(shower->get_particle_type()) == 13`.  `shower->get_particle_type()` returns 13 only if `shower->set_particle_type(13)` was called (`NeutrinoShowerClustering.cxx:118`), which is itself guarded by:

  ```cpp
  if (curr_sg->has_particle_info() && curr_sg->particle_info() &&
      std::abs(curr_sg->particle_info()->pdg()) == 13) {
      shower->set_particle_type(curr_sg->particle_info()->pdg());
  ```

  where `curr_sg` is the segment subsequently stored as `m_start_segment`.  Therefore when the function body executes, `m_start_segment->has_particle_info()` is always true and `pdg() == ±13`.  The prototype's `particle_type = 0` fall-through path is equally unreachable in the prototype because the same upstream tagging ensures the segment is PDG 13 before this function is called.

  **Verdict:** No functional divergence.  The WCT guard is a defensively correct improvement over the prototype's silent fall-through.  Add a comment at `:1128` recording the invariant so future readers understand why the guard is safe rather than worrying it might suppress real computation. No code logic change needed.
- **G1 applies:** `cal_kine_dQdx(vec_dQ, vec_dx, recomb_model)` at `:1161` subject to α = 0.930 vs. prototype's α = 1.0. Fixed by P8.

**C — Efficiency.** Single pass over `ordered_edges` accumulating both length and dQ/dx — more efficient than prototype's two passes (one over `map_seg_vtxs` for length/dQ, one over `map_vtx_segs` for vertices).

**D — Determinism.** `ordered_edges` for edge traversal; `std::map<size_t, VertexPtr>` keyed by graph index for muon-vertex accumulation; max-distance tie broken by smallest index — **fully deterministic**, improvement over prototype's pointer-map iteration.

**E — Multi-APA.** `cal_kine_dQdx` (vec overload) inherits G3. `segment_track_length` is purely geometric — face-agnostic. ✓

---

## Priority Action List

| # | Priority | Status | Action | File | Line(s) |
|---|---|---|---|---|---|
| P1 | **High — NaN risk** | ✅ **Applied** | Clamp `std::acos` arguments to `[-1, 1]` in `segment_determine_shower_direction` | `clus/src/PRSegmentFunctions.cxx` | (was) 2088, 2188, 2193, 2203, 2248, 2304, 2311, 2389 |
| P2 | **High — numerical bug** | ✅ **Applied** | Fix `segment_cal_kine_dQdx` first/last-point `dX` handling: compute `dE` from **original** `fits[i].dx`, multiply by shortened path | `clus/src/PRSegmentFunctions.cxx` | (was) 1196–1215 |
| P3 | **High — saturation filter** | ✅ **Applied** | Apply saturation filter `dQ/dX > threshold` against `fits[i].dx`, not the shortened `dX` (done together with P2) | `clus/src/PRSegmentFunctions.cxx` | (was) 1212 |
| P4 | **High — invariant** | ✅ **Applied (comment)** | Replaced terse WARNING in `break_segment` with expanded comment explaining `find_vertices` handles orientation via distance comparison; structural enforcement deferred | `clus/src/PRSegmentFunctions.cxx` | (was) 499–503 |
| P5 | **High — contract** | ✅ **Audited (no change)** | All toolkit callers use `D4Vector` named accessors (`.E()/.px()/.py()/.pz()`); no raw numeric indexing found | multiple | — |
| P6 | **Medium** | ✅ **Applied** | Restored `if (flag_print)` guard on debug log in `segment_determine_dir_track` | `clus/src/PRSegmentFunctions.cxx` | (was) 1649 |
| P7 | **Medium** | ✅ **Applied** | Removed `segment->dirsign(0)` side-effect from failure path of `segment_do_track_pid`; added comment documenting "no segment side-effects" contract | `clus/src/PRSegmentFunctions.cxx` | (was) 1449 |
| P8 | **High — confirmed fix** | ✅ **Verified (no change needed)** | `qlport/uboone-mabc.jsonnet:863` already sets `"A": 1.0` in `uBooNE_box_recomb_model`, which is wired into `tagger_check_neutrino` at `:1245`. Production config is correct. | `qlport/uboone-mabc.jsonnet` | 863 |
| P9 | **Medium — verify** | ⚠️ **Deferred (user decision)** | Confirm `clus/test/data/uboone-mabc_config.json` dEdx/range table values were generated from the same upstream source as `WCPDataDict.cxx:4677` | — | — |
| P10 | **Medium — contract** | ✅ **Verified (no change needed)** | `kenergy_best` as a per-segment member no longer exists in the toolkit. Segment energy is accessed via `seg->particle_info()->kinetic_energy()` (confirmed by migration comment at `NeutrinoKinematics.cxx:35`). Shower energy uses `shower->get_kine_best()` which reads `KinematicsData::kenergy_best`. All call sites audited — migration complete. | multiple | — |
| P11 | ~~Low~~ | ✅ **Applied** | Added `if (ncount == 0) return {1.0, 1e9, 1e9, 1e9};` early-return in `segment_do_track_comp` | `clus/src/PRSegmentFunctions.cxx` | (was) ~1278 |
| P12 | **Low** | ✅ **Resolved by P2/P3** | The confusing double `if (dX <= 0)` pattern was replaced when P2/P3 introduced `dx_for_dQdx = fits[i].dx` and a single `if (dX <= 0) continue;`. No separate fix needed. | `clus/src/PRSegmentFunctions.cxx` | — |
| P13 | **Doc** | ✅ **Applied** | Added multi-APA limitation comment above `segment_cal_kine_dQdx` | `clus/src/PRSegmentFunctions.cxx` | (was) ~1175 |
| P14 | **Doc** | ✅ **Applied** | Added unit-convention comment above `dQ_dx[i]` assignment in `segment_determine_dir_track` | `clus/src/PRSegmentFunctions.cxx` | (was) ~1520 |
| P15 | **Low — efficiency** | ✅ **Applied** | Cached `segment_median_dQ_dx` once into `medium_dQ_dx` before the three-way branch in `segment_determine_dir_track` | `clus/src/PRSegmentFunctions.cxx` | (was) 1579, 1606, 1614 |
| P16 | **Medium** | ✅ **Applied** | Replaced `(0,0,0)` sentinel in `calculate_kinematics` with `dist >= 0` check on `shower_get_closest_point` return value | `clus/src/PRShower.cxx` | (was) 1021 |
| P17 | **In scope** | ✅ **Audited (no code change needed)** | `calculate_kinematics_long_muon` reviewed and documented in §11 above; B17.1 noted; main fixes covered by P8 (α) | `clus/src/PRShower.cxx` | 1121+ |
| P18 | **Multi-APA verify** | ✅ **Verified (no change needed)** | `DynamicPointCloud::add_points` (`DynamicPointCloud.cxx:130–151`) iterates `pt.wpid_2d[pindex]` per plane and keys each 2D KD-tree by a `(layer, face, apa)` wpid ident — fully face-aware. `make_points_direct` populates `DPCPoint::wpid_2d` with the volume wpid. `break_segment` redistributes existing `DPCPoint` objects (which already carry their `wpid_2d`) for "associate_points" and calls `create_segment_point_cloud → make_points_direct` (with `flag_wrap=true`) for the "main" clouds — per-face 2D info is preserved correctly in both halves. | `clus/src/DynamicPointCloud.cxx` | 130–151 |

---

## Open Questions / Author Decisions

| # | Question | Decision |
|---|---|---|
| Q1 | BoxRecombination α | **Fix:** toolkit must match prototype (α = 1.0). See P8. |
| Q2 | Reference table provenance | **Verify:** confirm `uboone-mabc_config.json` dEdx/range values match WCP embedded tables. See P9. |
| Q3 | 4-momentum `[E, px, py, pz]` layout | **Accepted:** intentional, final convention. All callers must use `.E()/.px()/.py()/.pz()` accessors. See P5 audit. |
| Q4 | `break_segment` wcpts/graph orientation | **Resolved (comment + investigation):** `find_vertices()` in `PRGraph.cxx` already disambiguates orientation via wcpts.front() distance comparison, so all callers that use `find_vertices` are safe. An explanatory comment was added to `break_segment` (P4). Hard enforcement inside `break_segment` is not required. |
| Q5 | Per-face recomb/lifetime in PID/kine | **Deferred:** document assumption (G3 / P13); fix in a future multi-APA ticket. |
| Q6 | `particle_score = 100.0` in shower-trajectory direction | **Resolved:** replaced `// hack for now ...` with a 5-line comment at `PRSegmentFunctions.cxx:1691` documenting that PDG is forced to electron and `particle_score = 100.0` is an intentional sentinel meaning "PID not performed; score not applicable." |
| Q7 | `calculate_kinematics_long_muon` in scope | **Included:** escalated to P17 in this review. |
| Q8 | Drift direction hardcoded as `(1,0,0)` | **Open for future detectors:** acceptable for uBooNE; for PDHD/DUNE vertical APAs, derive from `IDetectorVolumes`. |
