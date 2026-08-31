# BDT Scorers Examination

> **Bug Fix Status**: Bugs 1, 3, 7 fixed. Bug 2 (numu_3_score not fed to xgboost) deferred — needs verification against training XML. Bugs 4-6 deferred — low severity / correct but fragile.
>
> **Efficiency Fix Status**: Efficiency 1 (TMVA::Reader per-event XML parsing) fixed — readers and BookMVA moved to configure() for both nue and numu scorers. All 16 nue readers (15 vector sub-BDTs + 1 xgboost) and all 5 numu readers (4 sub-BDTs + 1 xgboost) are now created once. Expected 10-100x speedup in BDT scoring. Efficiencies 2-4 deferred — subsumed by Efficiency 1 or low impact.

Files examined:
- `root/src/UbooneNueBDTScorer.cxx` (1707 lines) + `root/inc/WireCellRoot/UbooneNueBDTScorer.h`
- `root/src/UbooneNumuBDTScorer.cxx` (401 lines) + `root/inc/WireCellRoot/UbooneNumuBDTScorer.h`

---

## 1. Potential Bugs

### Bug 1: Division-by-zero / infinity in log-odds transformation (both files)

- **File**: `UbooneNueBDTScorer.cxx` line 1703; `UbooneNumuBDTScorer.cxx` line 400
- **Description**: The log-odds transformation `std::log10((1.0 + val1) / (1.0 - val1))` produces `+inf` when `val1 == 1.0`, `-inf` when `val1 == -1.0`, and `NaN` when `|val1| > 1.0`. TMVA BDT outputs are typically in `(-1, 1)`, but edge cases (e.g. a degenerate tree or extreme input values) can push `val1` to exactly `+/-1.0` or beyond. There is no clamping or guard before the transformation.
- **Severity**: Medium
- **Why it is a bug**: An infinite or NaN `nue_score` / `numu_score` can propagate downstream, causing undefined behavior in subsequent cuts or selections. The prototype code (using `TMath::Log10`) would have the same issue, so this may be a known limitation, but a clamp like `val1 = std::clamp(val1, -0.9999, 0.9999)` would be safer.

### Bug 2: numu_3_score computed but not fed to final xgboost model

- **File**: `UbooneNumuBDTScorer.cxx` lines 281, 385-387
- **Description**: `ti.numu_3_score` is computed at line 281 via `cal_numu_3_bdt()`, but it is never added as a variable to the final xgboost `reader` (lines 288-389). Only `cosmict_10_score`, `numu_1_score`, and `numu_2_score` are added (lines 385-387). If the XGBoost XML model was trained with `numu_3_score` as an input, the model will receive a wrong/missing feature; if it was intentionally omitted from the training, then the computation at line 281 is wasted.
- **Severity**: Medium -- needs verification against the training XML to determine whether this is a missing feature or dead code. If the XML expects the variable, this is a high-severity silent data corruption bug.

### Bug 3: cal_numu_3_bdt unconditionally evaluates without a fill gate

- **File**: `UbooneNumuBDTScorer.cxx` lines 245-264
- **Description**: Unlike the other sub-BDTs which have fill condition checks (e.g. `if (!ti.numu_cc_1_particle_type.empty())` at line 174, `if (!ti.cosmict_10_length.empty())` at line 124), `cal_numu_3_bdt` unconditionally calls `reader_numu_3.EvaluateMVA("MyBDT")` at line 261, ignoring the `default_val` parameter entirely. The `default_val` of `-0.2f` passed at line 281 is never returned. If the numu_cc_3 scalar fields happen to contain uninitialized or meaningless values (because the flag-3 check was never triggered for this event), the BDT will score garbage.
- **Severity**: Medium
- **Why it is a bug**: The function signature accepts `default_val` suggesting it should be used when the data is not filled, but the return path at line 263 always overwrites `val` with the BDT output. A guard like `if (ti.numu_cc_flag_3 == 0) return default_val;` would match the pattern of the other sub-BDTs.

### Bug 4: NaN guards placed after BookMVA, not before EvaluateMVA

- **File**: `UbooneNumuBDTScorer.cxx` lines 392-395
- **Description**: The NaN checks for `cosmict_4_angle_beam`, `cosmict_7_angle_beam`, `cosmict_7_theta`, and `cosmict_7_phi` are placed after `reader.BookMVA()` at line 389 but before `reader.EvaluateMVA()` at line 397. This ordering works correctly. However, the variables are bound to `ti` fields via pointer. The concern is that TMVA reads the float* at evaluation time, so the NaN fix at lines 392-395 does take effect. This is correct but fragile -- if someone reorders the code, the guards could be missed.
- **Severity**: Low (correct but fragile)

### Bug 5: Inconsistent vector size assumptions in vector sub-BDTs

- **File**: `UbooneNumuBDTScorer.cxx` lines 126-131; `UbooneNueBDTScorer.cxx` (all vector sub-BDTs)
- **Description**: In `cal_cosmict_10_bdt` (lines 124-138), the loop iterates over `ti.cosmict_10_length.size()` but accesses `ti.cosmict_10_vtx_z.at(i)`, `ti.cosmict_10_flag_shower.at(i)`, etc. If these vectors have different sizes (due to a fill bug upstream), `.at(i)` will throw `std::out_of_range`. The same pattern is used in every vector sub-BDT in both files. While `.at()` is bounds-checked (will throw rather than corrupt memory), no defensive size check is performed.
- **Severity**: Low -- the upstream fill code is expected to maintain consistent vector sizes, and `.at()` provides a safe crash rather than silent corruption. Still, a debug-level assertion would help catch upstream bugs.

### Bug 6: 15 scalar sub-BDT methods are dead code in nue scorer

- **File**: `UbooneNueBDTScorer.cxx` lines 182-908 (functions `cal_mipid_bdt` through `cal_vis_2_bdt`)
- **Description**: The 15 scalar sub-BDT methods (`cal_mipid_bdt`, `cal_gap_bdt`, `cal_hol_lol_bdt`, `cal_cme_anc_bdt`, `cal_mgo_mgt_bdt`, `cal_br1_bdt`, `cal_br3_bdt`, `cal_stemdir_br2_bdt`, `cal_trimuon_bdt`, `cal_br4_tro_bdt`, `cal_mipquality_bdt`, `cal_pio_1_bdt`, `cal_stw_spt_bdt`, `cal_vis_1_bdt`, `cal_vis_2_bdt`) are declared and fully implemented but never called from `cal_bdts_xgboost` or from `visit()`. The final xgboost model uses the raw scalar features directly (lines 1439-1696), not the intermediate BDT scores from these methods. These 15 methods were part of the older `cal_bdts()` variant (not ported), but remain as dead code.
- **Severity**: Low (no correctness impact; code maintenance burden and binary size impact only)

### Bug 7: Header include-guard typo

- **File**: `UbooneNueBDTScorer.h` line 13
- **Description**: The include guard is `WIRECELLROOT_UBOONENUEBTDSCORER_H` -- note `BTD` instead of `BDT`. This is `UBOONENUEBTDSCORER` rather than `UBOONENUEBDTSCORER`. The guard still works (it just needs to be unique), but it is inconsistent with the class name and could cause confusion if someone tries to match the pattern.
- **Severity**: Low (cosmetic, no functional impact)

---

## 2. Efficiency Improvements

### Efficiency 1: TMVA::Reader and BookMVA called on every visit (both files)

- **File**: `UbooneNueBDTScorer.cxx` lines 186-218, 238-249, 269-286, etc. (every sub-BDT and the final scorer); `UbooneNumuBDTScorer.cxx` lines 115-122, 163-172, etc.
- **Current approach**: Every time `visit()` is called (once per event), each sub-BDT function creates a `TMVA::Reader`, calls `AddVariable` for each feature, calls `BookMVA` (which parses the XML weight file from disk), evaluates, then destroys the reader. For the nue scorer, this means up to 16 XML files are parsed per event (15 vector sub-BDTs + 1 final xgboost). For the numu scorer, up to 5 XML files are parsed per event.
- **Suggested improvement**: Create the `TMVA::Reader` objects once during `configure()` and store them as member variables. `BookMVA` is the expensive operation (XML parsing, tree construction). The variables can be bound to member floats that are updated before each `EvaluateMVA` call. This is the standard TMVA usage pattern.
- **Expected impact**: High -- XML parsing and reader construction dominate the per-event cost. A typical BDT XML file is hundreds of KB to several MB. Eliminating redundant parsing could reduce BDT scoring time by 10-100x, depending on file sizes.

### Efficiency 2: Vector sub-BDTs loop with per-element EvaluateMVA

- **File**: `UbooneNueBDTScorer.cxx` lines 473-481, 525-539, 575-586, etc. (all 15 vector sub-BDTs); `UbooneNumuBDTScorer.cxx` lines 126-137, 176-190.
- **Current approach**: For vector-feature sub-BDTs, the code loops over all elements, evaluating the BDT for each element and tracking the min (or max). Each `EvaluateMVA` call traverses the full decision tree ensemble. This is inherent to the per-element scoring and cannot be avoided.
- **Suggested improvement**: Given that Efficiency 1 (persistent readers) is implemented, these loops become cheap. No further optimization is needed beyond Efficiency 1. However, if elements are independent and the number is large, batch evaluation or early termination (e.g. if a score already hits a threshold) could help in rare high-multiplicity events.
- **Expected impact**: Low (the loop body is fast once Efficiency 1 is done; element counts are typically small, O(1-10))

### Efficiency 3: Redundant variable registration in final xgboost reader

- **File**: `UbooneNueBDTScorer.cxx` lines 1438-1697 (the final xgboost reader setup in `cal_bdts_xgboost`)
- **Current approach**: The final xgboost reader registers ~160 variables. Many of these are the same fields that the 15 scalar sub-BDTs also register. This double registration is not a correctness issue (different readers) but contributes to the per-event overhead when combined with Efficiency 1.
- **Suggested improvement**: With persistent readers (Efficiency 1), this registration happens only once. No separate improvement needed.
- **Expected impact**: Low (subsumed by Efficiency 1)

### Efficiency 4: `const` correctness on TaggerInfo

- **File**: All sub-BDT methods in both `.h` and `.cxx` files
- **Current approach**: All sub-BDT methods take `Clus::PR::TaggerInfo& ti` (non-const reference), even though most sub-BDTs only read from `ti`. The non-const reference is needed because `TMVA::Reader::AddVariable` takes `float*` (not `const float*`), so the code cannot use a const reference. However, some methods also write NaN-guard fixes directly to `ti` fields (e.g., `UbooneNumuBDTScorer.cxx` lines 133, 186, 392-395; `UbooneNueBDTScorer.cxx` lines 1429-1432).
- **Suggested improvement**: For the NaN guards that modify `ti` directly, consider using local copies of the values to avoid mutating the shared `TaggerInfo` state. This prevents subtle side effects where a NaN-guarded value in one sub-BDT affects a different sub-BDT that reads the same field.
- **Expected impact**: Low (correctness/clarity improvement rather than performance)

---

## 3. Algorithm and Code Explanation

### 3.1 UbooneNueBDTScorer

#### General Purpose

`UbooneNueBDTScorer` implements a multi-stage Boosted Decision Tree (BDT) scoring pipeline for electron neutrino (nue) charged-current event selection in the MicroBooNE liquid argon TPC experiment. It is an `IEnsembleVisitor` that runs after `TaggerCheckNeutrino` has filled the `TaggerInfo` and `KineInfo` data structures. Its output is a single float, `ti.nue_score`, written back into `TaggerInfo`.

#### Key Data Structures

- **`Clus::PR::TaggerInfo` (ti)**: Contains ~300+ physics variables computed upstream by various tagger algorithms. These include:
  - Scalar flags and features for MIP identification (`mip_*`), gap analysis (`gap_*`), shower topology (`br1_*`, `br3_*`, `br4_*`), pi0 reconstruction (`pio_*`), track/shower discrimination (`stw_*`, `sig_*`, `lol_*`, `tro_*`), and more.
  - Vector features for per-segment or per-candidate evaluations (e.g., `br3_3_v_energy`, `tro_1_v_particle_type`, etc.).
  - Fill-condition flags (`mip_filled`, `gap_filled`, `br_filled`, `pio_filled`, `vis_1_filled`, `vis_2_filled`, `mip_quality_filled`) indicating whether a particular tagger successfully ran.
  - Output score fields (`br3_3_score`, `nue_score`, etc.).

- **`Clus::PR::KineInfo` (ki)**: Contains kinematic reconstruction results, primarily `kine_reco_Enu` (reconstructed neutrino energy).

- **TMVA::Reader**: ROOT/TMVA framework class that loads XML weight files and evaluates BDT/XGBoost models. Variables are bound by pointer; the reader reads the pointed-to floats at evaluation time.

#### Algorithm Walkthrough

The `cal_bdts_xgboost` method (lines 1408-1707) implements a two-tier scoring architecture:

**Tier 1 -- Vector Sub-BDT Scoring (lines 1412-1426)**:
15 sub-BDTs are evaluated, each operating on vector features (per-segment or per-candidate). For each sub-BDT:
1. A `TMVA::Reader` is created and configured with the relevant feature variables.
2. The XML weight file is loaded via `BookMVA`.
3. The code loops over all elements in the vector, evaluating the BDT for each element.
4. The minimum score across all elements is taken (most background-like element dominates).
5. If the vector is empty, a hand-tuned default value is returned (e.g., 0.3 for br3_3, 0.42 for br3_5).

The 15 vector sub-BDTs and their default values:
| Sub-BDT | Default | Purpose |
|---------|---------|---------|
| br3_3 | 0.3 | Bad-reconstruction type 3, sub-check 3 |
| br3_5 | 0.42 | Bad-reconstruction type 3, sub-check 5 |
| br3_6 | 0.75 | Bad-reconstruction type 3, sub-check 6 |
| pio_2 | 0.2 | Pi0 type 2 (non-vertex-attached) |
| stw_2 | 0.7 | Shower-to-wall type 2 |
| stw_3 | 0.5 | Shower-to-wall type 3 |
| stw_4 | 0.7 | Shower-to-wall type 4 |
| sig_1 | 0.59 | Single-shower pi0 type 1 |
| sig_2 | 0.55 | Single-shower pi0 type 2 |
| lol_1 | 0.85 | Low-energy overlapping type 1 |
| lol_2 | 0.7 | Low-energy overlapping type 2 |
| tro_1 | 0.28 | Track overclustering type 1 |
| tro_2 | 0.35 | Track overclustering type 2 |
| tro_4 | 0.33 | Track overclustering type 4 |
| tro_5 | 0.5 | Track overclustering type 5 |

**Tier 1b -- Variable Protection (lines 1429-1432)**:
Four specific variables receive NaN/overflow guards:
- `mip_min_dis` clamped to max 1000
- `mip_quality_shortest_length` clamped to max 1000
- `mip_quality_shortest_angle` NaN -> 0
- `stem_dir_ratio` NaN -> 1.0

**Tier 2 -- Final XGBoost Model (lines 1438-1706)**:
A single XGBoost model (loaded from `XGB_nue_seed2_0923.xml`) takes ~160 input features:
- `match_isFC` (fiducial containment flag)
- `kine_reco_Enu` (reconstructed neutrino energy)
- All raw scalar features from the various taggers (~130 features)
- The 15 vector sub-BDT scores from Tier 1
- `mip_vec_dQ_dx_0` through `mip_vec_dQ_dx_19` (20 dQ/dx profile bins)

The model is evaluated only if `ti.br_filled == 1` (line 1701). The raw output is transformed via log-odds: `nue_score = log10((1 + val) / (1 - val))` (line 1703). If `br_filled != 1`, the score defaults to `-15.0`.

**Note on the 15 scalar sub-BDTs**: The methods `cal_mipid_bdt` through `cal_vis_2_bdt` (lines 182-908) are fully implemented but never called. These were part of the older `cal_bdts()` approach (using a cascade of intermediate BDTs). The xgboost approach feeds raw features directly, making these methods dead code.

#### Event Selection Criteria

- The scorer only runs if `m_nue_xgboost_xml` is configured (line 150).
- It requires a grouping with a `TrackFitting` object (lines 155-166).
- The final BDT is only evaluated if `ti.br_filled == 1` (line 1701), meaning the bad-reconstruction tagger must have successfully run.
- The `pio_2_bdt` has an additional gate: `pio_filled == 1 && pio_flag_pio == 0` (line 795).

### 3.2 UbooneNumuBDTScorer

#### General Purpose

`UbooneNumuBDTScorer` implements a BDT scoring pipeline for muon neutrino (numu) charged-current event selection. It follows the same `IEnsembleVisitor` pattern as the nue scorer. Its output is `ti.numu_score`.

#### Key Data Structures

Same `TaggerInfo` and `KineInfo` as the nue scorer, but focuses on different feature subsets:
- **numu_cc_1_***: Direct-muon track features (per-track vector: particle type, length, dQ/dx, daughter counts)
- **numu_cc_2_***: Long-muon shower features (per-shower vector: lengths, daughter counts)
- **numu_cc_3_***: Indirect-muon features (scalars: particle type, various lengths, daughter counts)
- **cosmict_{2-9}_***: Cosmic-ray tagger features (scalars for tagger flags 2-9)
- **cosmict_10_***: Upstream-dirt cluster features (per-cluster vector)

#### Algorithm Walkthrough

The `cal_numu_bdts_xgboost` method (lines 275-401) implements a two-tier architecture:

**Tier 1 -- Sub-BDT Scoring (lines 279-282)**:
Four sub-BDTs are evaluated:

1. **cal_numu_1_bdt** (lines 151-194): Scores direct-muon candidates. Vector BDT over `numu_cc_1_*` features. Returns the maximum score (most signal-like candidate wins). Default: -0.4. Guards `inf` in `dQ_dx_cut` by replacing with 10.

2. **cal_numu_2_bdt** (lines 204-235): Scores long-muon shower candidates. Vector BDT over `numu_cc_2_*` features. Returns maximum score. Default: -0.1.

3. **cal_numu_3_bdt** (lines 245-264): Scores indirect-muon features. Scalar BDT on `numu_cc_3_*` features. Always evaluates (no fill gate -- see Bug 3). Default: -0.2 (unused).

4. **cal_cosmict_10_bdt** (lines 104-141): Scores upstream-dirt clusters. Vector BDT over `cosmict_10_*` features. Returns minimum score (most cosmic-like element wins). Default: 0.7. Guards NaN in `angle_beam` by replacing with 0.

**Tier 2 -- Final XGBoost Model (lines 284-401)**:
A single XGBoost model (loaded from `numu_scalars_scores_0923.xml`) takes ~70 input features:
- `numu_cc_flag_3` and all `numu_cc_3_*` scalar features
- `cosmict_flag_{2-9}` and all associated scalar features
- Top-level flags: `cosmic_flag`, `cosmic_filled`, `cosmict_flag`, `numu_cc_flag`, `cosmict_flag_1`
- `kine_reco_Enu` and `match_isFC`
- Three of the four sub-BDT scores: `cosmict_10_score`, `numu_1_score`, `numu_2_score`
- Notably missing: `numu_3_score` (see Bug 2)

NaN guards are applied to `cosmict_4_angle_beam`, `cosmict_7_angle_beam`, `cosmict_7_theta`, `cosmict_7_phi` (lines 392-395).

The raw output is transformed via the same log-odds formula: `numu_score = log10((1 + val) / (1 - val))` (line 400).

#### Key Differences from Nue Scorer

1. **Fewer sub-BDTs**: 4 vs. 15 (only the ones needed by xgboost).
2. **No fill gate on final evaluation**: The numu xgboost model is always evaluated (no `br_filled` check), unlike the nue scorer which checks `br_filled == 1`.
3. **Max vs. min aggregation**: The numu vector sub-BDTs use both max (numu_1, numu_2: looking for the best signal candidate) and min (cosmict_10: looking for the worst cosmic contamination). The nue vector sub-BDTs all use min (looking for the most background-like element).
4. **No dead code**: Unlike the nue scorer, the numu scorer does not carry unused sub-BDT methods.

### 3.3 Common Patterns

Both scorers follow identical structural patterns:
- **Factory registration** via `WIRECELL_FACTORY` macro
- **Configuration** via `configure()` / `default_configuration()` with XML weight file paths resolved through `Persist::resolve()`
- **Visitor pattern**: `visit()` extracts `TaggerInfo` and `KineInfo` from the ensemble's grouping, then delegates to the top-level xgboost scorer
- **TMVA variable binding**: Features are bound to `TMVA::Reader` by pointer (either to `ti` member fields or to local float copies)
- **Log-odds output transformation**: Both convert raw BDT output to a log-likelihood ratio via `log10((1+v)/(1-v))`
- **Default values**: Vector sub-BDTs use a sentinel initialization (`1e9` or `-1e9`) and fall back to hand-tuned defaults when vectors are empty
