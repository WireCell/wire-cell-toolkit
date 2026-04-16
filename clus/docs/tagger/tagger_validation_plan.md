# Tagger Validation Plan

## 1. Overview and Goals

The neutrino ID taggers (cosmic, numu, nue, ssm, single photon) and BDT scoring code have been ported from the prototype (`WCPPID::NeutrinoID`) to the toolkit (`PatternAlgorithms`). Validation against prototype output is needed to establish confidence in the ported code.

The fundamental challenge is that pattern recognition upstream of the taggers produces **functionally equivalent but not bit-identical** results between prototype and toolkit. Different neutrino vertices may be selected, different segments identified, or different shower clustering applied. Therefore:

- Tagger input features will differ event-by-event.
- BDT scores will not match exactly.
- The comparison must be **statistical and distributional**, not a direct value match.

**Goals:**

1. **Reasonableness check (toolkit alone):** Verify all toolkit tagger output variables are within physically sensible ranges and have internally consistent fill-gate behavior.
2. **Functional equivalence (toolkit vs prototype):** Verify that distributions of BDT features and final scores are statistically compatible between the two implementations.

---

## 2. Toolkit Side: Getting Tagger Output into ROOT Files

### 2.1 Current State

`TaggerInfo` (~500+ fields) and `KineInfo` (~20 fields) are defined in `clus/inc/WireCellClus/NeutrinoTaggerInfo.h` and stored in-memory on the `TrackFitting` object after the tagger chain runs:

```
TaggerCheckNeutrino::visit()     -> fills TaggerInfo + KineInfo -> TrackFitting
UbooneNumuBDTScorer::visit()     -> reads TaggerInfo, writes numu_score + sub-scores
UbooneNueBDTScorer::visit()      -> reads TaggerInfo + KineInfo, writes nue_score + 30 sub-scores
```

The existing `UbooneMagnifyTrackingVisitor` only writes spatial/charge data (`T_rec_charge`, `T_proj_data`) to ROOT. **No visitor currently writes `TaggerInfo` or `KineInfo` to disk.**

### 2.2 New Decoupled Output Visitor

Create a new, standalone `IEnsembleVisitor` component — **`UbooneTaggerOutputVisitor`** — that writes tagger variables to a ROOT file. This visitor:

- Is completely independent of `UbooneMagnifyTrackingVisitor`; it lives in separate files.
- Accesses `TaggerInfo` and `KineInfo` only through the `TrackFitting` public API (`get_tagger_info()`, `get_kine_info()`), accessible via `grouping->get_track_fitting()`.
- Can be enabled or disabled via Jsonnet configuration without affecting any other output.
- If ROOT storage is later replaced (e.g., HDF5 or tensor output), only this visitor needs to change.

**New files:**
- `root/inc/WireCellRoot/UbooneTaggerOutputVisitor.h`
- `root/src/UbooneTaggerOutputVisitor.cxx`

No build system changes are needed — the `root/` package's `wscript_build` uses `bld.smplpkg()` which auto-discovers all `.cxx` files in `src/`.

### 2.3 Implementation

**Status: DONE.** The visitor has been implemented in:
- `root/inc/WireCellRoot/UbooneTaggerOutputVisitor.h`
- `root/src/UbooneTaggerOutputVisitor.cxx`

The visitor opens the existing tracking ROOT file (created earlier by `UbooneMagnifyTrackingVisitor`) in **UPDATE** mode and adds two new trees. No changes to the existing tracking visitor are needed.

Key implementation notes:
- Uses `TFile::Open(filename, "UPDATE")` — the file must already exist when this runs.
- `visit()` is `const`, so `TaggerInfo`/`KineInfo` are copied as local variables before passing branch addresses.
- All branches are registered explicitly (no loop-based pointer arithmetic) for auditability.
- Build system auto-discovers the new `.cxx` — no `wscript_build` changes needed.

### 2.4 ROOT Tree Structure

Inspection of the prototype output file (`prototype_base/nue_5384_130_6528.root`) confirmed the actual tree names are **`T_tagger`** and **`T_kine`** (not `T_eval`/`T_BDTvars`/`T_KINEvars` as initially assumed).

#### T_tagger (one entry per event)
All fields of `TaggerInfo` (~1200 branches) plus vertex position from `KineInfo`:

| Branch | Source | Notes |
|---|---|---|
| `nu_x/F`, `nu_y/F`, `nu_z/F` | KineInfo (`kine_nu_x/y/z_corr`) | vertex position |
| `cosmic_flag/F` … | TaggerInfo | all cosmic tagger fields |
| `numu_cc_flag/F`, `numu_score/F` … | TaggerInfo | all numu fields + BDT scores |
| `nue_score/F` … | TaggerInfo | all nue fields + 30 sub-BDT scores |
| `photon_flag/F` … | TaggerInfo | single photon fields |
| `ssm_flag/F` … | TaggerInfo | SSM fields (use `-999` sentinel when not filled) |
| `match_isFC/F` | TaggerInfo | fiducial flag (not in prototype file, but included in toolkit output) |

Prototype has 1216 branches. Toolkit output is expected to match closely; any extra branches (like `match_isFC`) are additions.

#### T_kine (one entry per event)
All 21 fields of `KineInfo`:
- Scalar: `kine_nu_x/y/z_corr`, `kine_reco_Enu`, `kine_reco_add_energy`, `kine_pio_mass`, `kine_pio_flag`, `kine_pio_vtx_dis`, `kine_pio_energy_1/2`, `kine_pio_theta_1/2`, `kine_pio_phi_1/2`, `kine_pio_dis_1/2`, `kine_pio_angle`
- Vector (STL): `kine_energy_particle`, `kine_energy_info`, `kine_particle_type`, `kine_energy_included`

### 2.5 Jsonnet Configuration

**Status: DONE.** Added to `cfg/pgrapher/common/clus.jsonnet`:

```jsonnet
// Write T_tagger and T_kine trees into the existing tracking output ROOT file.
// Must run AFTER numu_bdt_scorer and nue_bdt_scorer (BDT scores must be filled).
// Must run AFTER UbooneMagnifyTrackingVisitor (file must already exist to UPDATE).
tagger_output(name="", output_filename="tracking_proj.root") :: {
    type: "UbooneTaggerOutputVisitor",
    name: prefix + name,
    data: {
        grouping: "live",
        output_filename: output_filename,
    }
},
```

Wired into the pipeline in `qlport/uboone-mabc.jsonnet` — `tagger_output_visitor` is conditionally appended after `tracking_visitor`, gated on `tracking_output != ""`:

```jsonnet
local tagger_output_visitor = cm.tagger_output(output_filename=tracking_output);
...
+ (if tracking_output != "" then [tracking_visitor, tagger_output_visitor] else []);
```

The visitor is **on by default** whenever `tracking_output` is set (same gate as the tracking visitor itself). In the `qlport/uboone-mabc.jsonnet` test configuration, the output filename is `track_com_<run>_<event>.root`.

---

## 3. Prototype Side: What to Extract

### 3.1 No Code Changes Needed

For the prototype side, no new code is required. The output ROOT files already exist. The steps are:

1. **Locate prototype output files** for a representative event sample. These should be the same input events used by the toolkit integration tests (e.g., `nuselEval_5384_130_6501.root` from the BATS tests).

2. **Verify tree and branch names.** Open a prototype ROOT file and inspect the tree structure:
   ```bash
   root -l prototype_output.root
   .ls        # list top-level objects
   T_eval->Print()      # list branches
   T_BDTvars->Print()
   T_KINEvars->Print()
   ```
   Document any naming differences from what the toolkit will produce.

   Any example file can be found @ /home/xqian/work/scratch_wcgpu1/toolkit-dev/toolkit/prototype_base/nue_5384_130_6528.root

3. **Note any unit or normalization differences.** Check whether:
   - Energies are in the same units (MeV vs GeV)
   - Lengths are in the same units (cm vs mm)
   - Any fields are sign-flipped or offset

4. **Prepare an event list** with `(runNo, subRunNo, eventNo)` tuples that exist in both prototype and toolkit outputs for matched comparison.

---

## 4. Comparison Methodology

### 4.1 Phase 1: Sanity Checks on Toolkit Output Alone

Run these checks before any cross-comparison to confirm the toolkit output is self-consistent.

#### 4.1.1 No NaN or Inf
Every scalar field in `TaggerInfo` and `KineInfo` must be finite for all events. This catches uninitialized-variable bugs or divide-by-zero.

```python
for branch_name in scalar_branches:
    values = np.array(tree[branch_name].array())
    assert np.all(np.isfinite(values[values != -999])), f"{branch_name} has non-finite values"
```

#### 4.1.2 Physical Range Checks
Verify each field is within its expected physical range:

| Field category | Expected range |
|---|---|
| Energy fields (`kine_reco_Enu`, `kine_pio_energy_*`, etc.) | >= 0 (or == -999 sentinel) |
| Angle fields | [0, π] or [-π, π] as appropriate |
| Length fields | >= 0 (or == -999 sentinel) |
| Flag fields (`cosmic_flag`, `numu_cc_flag`, etc.) | {0, 1} |
| `match_isFC` | {0, 1} |
| BDT scores (`numu_score`, `nue_score`) | Check for Inf if log-transform is applied; expect finite values after clamp |
| Sub-BDT scores (`*_score`) | Typically [-1, 1] for TMVA or [0, 1] for probability output; verify against BDT scorer code |
| SSM sentinel fields | -999 when SSM not triggered |

#### 4.1.3 Fill-Gate Consistency
The tagger structs use `*_filled` flags to indicate whether a sub-tagger ran. When `filled==1`, the associated feature fields should have non-default values:

| If this field == 1 | These should be non-default |
|---|---|
| `cosmic_filled` | `cosmic_n_solid_tracks`, `cosmic_energy_main_showers`, etc. |
| `gap_filled` | `gap_n_bad`, `gap_n_points`, `gap_energy`, etc. |
| `mip_quality_filled` | `mip_quality_energy`, `mip_quality_overlap`, etc. |
| `mip_filled` | `mip_*` fields including `mip_vec_dQ_dx_0`...`_19` |
| `pio_filled` | `pio_*` fields |
| `br_filled` | `br1_*`, `br2_*`, `br3_*`, `br4_*` fields |

Also check vector-valued features for consistent lengths per event:
- `numu_cc_1_length`, `numu_cc_1_fraction`, `numu_cc_1_muon_length`, ... should all be the same length.
- `tro_1_v_particle_type`, `tro_1_v_flag_*`, ... should all be the same length.

#### 4.1.4 Event Fraction Checks
Check that the fraction of events where each major sub-tagger fires is plausible:

| Sub-tagger | Expected fraction (approximate, for neutrino-enhanced sample) |
|---|---|
| `cosmic_filled == 1` | ~100% |
| `numu_cc_flag > 0` | 20-40% (numu CC fraction) |
| `nue_score > 0.7` | 1-5% (nue CC signal region) |
| `ssm_flag_st_kdar >= 0` | Rare; depends on sample |
| `match_isFC == 1` | ~50-70% (FV fraction) |

These are approximate; the key check is that no sub-tagger fires on 0% or 100% of events when it should not.

---

### 4.2 Phase 2: Multi-Event Distributional Comparison

Run both prototype and toolkit on the same sample of at least 100 events (ideally 500-1000 for stable tail statistics). Match events by `(runNo, subRunNo, eventNo)`.

#### 4.2.1 1D Distribution Comparison
For every scalar field in `TaggerInfo` and `KineInfo`:
- Overlay histograms (toolkit in blue, prototype in red)
- Compute Kolmogorov-Smirnov (KS) test p-value between the two distributions
- Compute mean and RMS for each; report the fractional difference in mean

Flag any variable where:
- KS p-value < 0.01 (distributions are statistically inconsistent)
- Mean fractional difference > 20%

**Important:** Exclude sentinel values (-999) before computing statistics, or check them separately.

For the SSM fields (many have -999 default), compare separately:
- Fraction of events where the field is non-sentinel (i.e., SSM fired)
- Distribution among events where SSM fired

#### 4.2.2 BDT Score Comparison (Priority)
The final BDT scores are the most important output. For `numu_score` and `nue_score`:
- Overlay distributions and compute KS p-value
- Compute efficiency at representative working points (e.g., nue_score > 7.0, numu_score > 0.9)
- Compare the efficiency difference: expect < 5% relative difference at signal working points

For the 30+ sub-BDT scores, use the same overlay+KS approach.

#### 4.2.3 Event-by-Event Correlation (2D Scatter)
For matched events, scatter plot toolkit value vs prototype value for each scalar variable:
- Ideal case: points along the y = x diagonal
- Compute Pearson correlation coefficient *r*
- Compute the median absolute deviation from the diagonal as a measure of scatter

Expected correlation levels:
- High-level flags (`match_isFC`, `cosmic_flag`): should agree very often (*r* > 0.95)
- Final BDT scores (`numu_score`, `nue_score`): *r* > 0.9 is the target
- Individual BDT input features: accept *r* > 0.7 for most, since they depend directly on upstream PR variations

#### 4.2.4 Vector-Valued Features
For vector branches (e.g., `numu_cc_1_length`, `tro_1_v_particle_type`):
1. Compare the **vector length distribution** (number of candidates per event). Large differences indicate different numbers of track/shower candidates from PR.
2. Compare the **aggregated reduced value** (the BDT score derived from the vector, e.g., `numu_1_score = max over vector`). This is ultimately what matters.
3. Do **not** attempt direct element-wise comparison of vector entries; the ordering within the vector may differ even when the physical content is equivalent.

---

### 4.3 Phase 3: Event-Level Diagnosis for Outliers

For events identified as outliers in Phase 2 (large deviation from diagonal, or one side fires a tagger when the other does not):

1. **Check the neutrino vertex.** If the toolkit and prototype selected different clusters as the primary vertex, all downstream tagger outputs will naturally differ. This is expected and not a bug.

2. **Check the pattern recognition graph.** Use `PatternDebugIO` (see `clus/inc/WireCellClus/PatternDebugIO.h`) to serialize the input state for that event to JSON, then inspect:
   - Number of vertices, segments, showers in the PR graph
   - Position of the main neutrino vertex
   - Segment and shower counts

3. **Classify each outlier** as one of:
   - **Class A (expected):** Different neutrino vertex selected upstream
   - **Class B (expected):** Different number of track/shower candidates (segment count differs)
   - **Class C (investigate):** Same vertex, same segment structure, but BDT feature differs — this indicates a porting bug

Only Class C events require code-level debugging.

---

### 4.4 Acceptance Criteria

| Category | Criterion | Status |
|---|---|---|
| Sanity | No NaN or Inf in any scalar field | Pass / Fail |
| Sanity | All flag fields in {0, 1} | Pass / Fail |
| Sanity | Fill-gate consistency: 100% of events | Pass / Fail |
| Sanity | Fraction of events per sub-tagger in plausible range | Pass / Fail |
| Distribution | KS p-value > 0.01 for >= 95% of scalar variables | Pass / Fail |
| Distribution | `numu_score` mean difference < 10% relative | Pass / Fail |
| Distribution | `nue_score` mean difference < 10% relative | Pass / Fail |
| Distribution | Nue efficiency at working point differs < 5% relative | Pass / Fail |
| Correlation | `numu_score` Pearson *r* > 0.90 | Pass / Fail |
| Correlation | `nue_score` Pearson *r* > 0.90 | Pass / Fail |
| Correlation | >= 85% of scalar BDT input features have *r* > 0.70 | Pass / Fail |
| Diagnosis | All Class C outliers investigated and resolved | Pass / Fail |

---

## 5. Unit Test Strategy

### 5.1 Doctest Tests (Build-time, `clus/test/`)

#### 5.1.1 TaggerInfo Default Initialization Test
File: `clus/test/doctest_tagger_info.cxx`

```cpp
TEST_CASE("TaggerInfo default initialization") {
    WireCell::Clus::PR::TaggerInfo ti{};
    // Check non-zero sentinels are correct
    CHECK(ti.cosmic_flag == Approx(1.0f));
    CHECK(ti.gap_flag == Approx(1.0f));
    CHECK(ti.mip_quality_flag == Approx(1.0f));
    CHECK(ti.mip_flag == Approx(1.0f));
    // ... all non-zero defaults from init_tagger_info_review.md
    // Check zero-init fields
    CHECK(ti.cosmic_filled == Approx(0.0f));
    CHECK(ti.numu_cc_flag == Approx(0.0f));
    // Check SSM sentinels
    CHECK(ti.ssm_kine_energy == Approx(-999.0f));
    // ... etc.
}
```

This catches any accidental changes to default values in `NeutrinoTaggerInfo.h`.

#### 5.1.2 KineInfo Default Initialization Test
File: same `doctest_tagger_info.cxx`

```cpp
TEST_CASE("KineInfo default initialization") {
    WireCell::Clus::PR::KineInfo ki{};
    CHECK(ki.kine_reco_Enu == Approx(0.0f));
    CHECK(ki.kine_pio_flag == 0);
    CHECK(ki.kine_energy_particle.empty());
    // ...
}
```

### 5.2 Round-Trip Serialization Test (ROOT)

File: `root/test/doctest_tagger_output.cxx` (or a BATS test if ROOT is not available at doctest time)

```
1. Construct TaggerInfo with known non-default values for all fields
2. Write to a TTree via UbooneTaggerOutputVisitor::write_T_BDTvars()
3. Read back from the TTree
4. Check every field matches the original value within float precision
```

This tests the branch registration code is complete and correct — any field added to `TaggerInfo` but not registered as a branch will be caught when the round-trip produces the wrong value.

### 5.3 BATS Integration Tests

#### 5.3.1 End-to-End Tagger Output Test
Extend `clus/test/test-porting.bats` with a new test function:

```bash
@test "tagger output file is produced" {
    do_prep "tagger_output"
    # Download test input file (same as qlport/steiner tests)
    local url="$raw_url/qlport/rootfiles/nuselEval_5384_130_6501.root"
    local dat="$(download_file "$url")"
    local cfg="$(relative_path test-porting/tagger_output/main.jsonnet)"
    
    # Run pipeline with tagger output enabled
    run_idempotently -s "$cfg" -s "$dat" -t "tagger_output.root" -t "$log" -- \
        bash -c "wire-cell -A infiles=$dat $cfg > $log 2>&1"
    
    # Verify output file exists and has expected content
    test -f tagger_output.root
    
    # Verify trees exist and have entries
    root -l -b -q "check_tagger_output.C(\"tagger_output.root\")"
}
```

The helper macro `check_tagger_output.C` checks:
- `T_eval`, `T_BDTvars`, `T_KINEvars` trees all exist
- Each tree has >= 1 entry
- `numu_score` and `nue_score` are finite (not 0 or NaN -- confirms BDT scorers ran)
- `match_isFC` is 0 or 1

#### 5.3.2 Regression Test (Stability)
Once a known-good output is established, add a digest comparison:
- Extract `numu_score`, `nue_score`, `kine_reco_Enu` for the test event
- Save as a historical baseline file (following the existing `saveout -c history` pattern)
- Future runs diff against the baseline; any change triggers investigation

### 5.4 Comparison Script

File: `clus/test/compare_tagger_outputs.py` (or `root/scripts/`)

A Python script using `uproot` (preferred for HEP CI environments) or `ROOT` that:

```
Usage: python compare_tagger_outputs.py toolkit_out.root prototype_out.root [--report report.html]

Steps:
1. Open both files, find matching events by (runNo, subRunNo, eventNo)
2. Phase 1: Sanity checks on toolkit file
3. Phase 2: For each scalar variable:
   a. Overlay histograms (matplotlib or ROOT)
   b. KS test (scipy.stats.ks_2samp)
   c. Correlation coefficient (scipy.stats.pearsonr)
4. Phase 2b: BDT score efficiency comparison at working points
5. Output: summary table (pass/fail per variable) + PDF/HTML with plots
6. Exit code 0 if all acceptance criteria pass, non-zero otherwise
```

Key Python dependencies: `uproot`, `numpy`, `scipy`, `matplotlib` (all standard in LArSoft / HEP Python environments).

---

## 6. Additional Considerations

### 6.1 Events with No Neutrino Vertex

When `TaggerCheckNeutrino` finds no neutrino vertex (e.g., all clusters are purely cosmic), it may not call any tagger. Both `TaggerInfo` and `KineInfo` stay at defaults. Verify that:
- The output ROOT file still gets an entry for these events (all-default values)
- The defaults are the same between prototype and toolkit (check against `init_tagger_info_review.md`)

### 6.2 SSM Sentinel Values (-999)

Many SSM fields use -999 as a sentinel (not -1 or 0). These sentinel values will dominate any histogram if included. In the comparison script:
- Always check and exclude sentinel values separately
- Report the **fraction of events where the field is non-sentinel** as a first-level metric
- Compare the distribution only among events where the SSM tagger actually fired

### 6.3 Log-Transform for Final BDT Scores

The `nue_score` and `numu_score` in the prototype are stored as `log10((1+v)/(1-v))` where `v` is the raw BDT output. If `v` is very close to +/-1, this diverges. Verify that:
- The toolkit BDT scorers (`UbooneNueBDTScorer`, `UbooneNumuBDTScorer`) clamp `v` before the log transform
- The clamping matches the prototype
- The output ROOT branch stores the **post-transform** value (to match prototype convention)

See `root/src/UbooneNueBDTScorer.cxx` for the implementation details.

### 6.4 Vector Branch Length Mismatches

For vector branches (`numu_cc_1_length`, `tro_1_v_particle_type`, etc.), the number of entries per event (the vector length) can differ between toolkit and prototype if PR found different numbers of segments or showers. This is expected. In the comparison:
- Do not fail on different vector lengths
- Compare the **reduced scalar score** (what gets passed to the final BDT)
- Optionally, compare the distribution of vector lengths as a diagnostic

### 6.5 Float Precision

`TaggerInfo` uses `float` (32-bit). ROOT `TTree` branches should be declared as `Float_t` (not `Double_t`) to avoid precision loss. Check especially the SSM fields that include position coordinates (`ssm_vtxX/Y/Z`) — these should remain float.

### 6.6 `singlephoton_tagger` Overlap with `nue_tagger`

The `shw_sp_*` fields are filled by both `nue_tagger()` and `singlephoton_tagger()` — the functions share field name prefixes but fill different sub-sets. Ensure the ROOT tree branches capture all `shw_sp_*` fields correctly, and that the comparison script distinguishes between:
- `shw_sp_*` fields filled by `nue_tagger` (shower spine features)
- `shw_sp_*` fields filled by `singlephoton_tagger` (single-photon hypothesis features)

Reference: `clus/docs/tagger/nue_singlephoton_tagger_review.md`.

### 6.7 Run/Subrun/Event Number Injection

For the BATS integration tests that run on a single downloaded input file, the event number comes from the input ROOT file metadata. `UbooneNumuBDTScorer.cxx` reads `runNo` from the `Trun` tree in the input. The new `UbooneTaggerOutputVisitor` should obtain these the same way — either via the input tree or via Jsonnet configuration — and write them as the event identifier branches in T_eval.

### 6.8 Future: Tensor Serialization Alternative

The toolkit's primary output uses tensor serialization (`.tar.gz` via `MultiAlgBlobClustering`). Adding `TaggerInfo`/`KineInfo` to the tensor output is a viable future path that would not require ROOT and would integrate naturally into the WCT data model. The ROOT output visitor described here is a validation tool; the long-term production path should be discussed separately.

---

## 7. Implementation Sequence

| Step | Task | Notes |
|---|---|---|
| 1 | Verify prototype tree names and branch names | Open prototype ROOT file, run `T_eval->Print()` etc. |
| 2 | Implement `UbooneTaggerOutputVisitor` | Follow `UbooneMagnifyTrackingVisitor` pattern |
| 3 | Add Jsonnet `tagger_output()` to `cfg/pgrapher/common/clus.jsonnet` | Default OFF per project convention |
| 4 | Write validation-specific Jsonnet config | Enables `tagger_output` visitor |
| 5 | Run on single test event | Confirm ROOT file produced with expected trees |
| 6 | Run Phase 1 sanity checks | Use comparison script or manual ROOT inspection |
| 7 | Prepare multi-event sample (100+ events) | Same input files as BATS tests |
| 8 | Run Phase 2 distributional comparison | Use `compare_tagger_outputs.py` |
| 9 | Investigate any Class C outliers | Use PatternDebugIO for PR graph inspection |
| 10 | Add doctest for TaggerInfo/KineInfo defaults | `clus/test/doctest_tagger_info.cxx` |
| 11 | Add round-trip serialization test | `root/test/doctest_tagger_output.cxx` |
| 12 | Add BATS end-to-end and regression tests | Extend `test-porting.bats` |
| 13 | Document results | Add `tagger_validation_results.md` under this directory |

---

## 8. Key File Reference

| File | Role |
|---|---|
| `clus/inc/WireCellClus/NeutrinoTaggerInfo.h` | TaggerInfo and KineInfo struct definitions |
| `clus/src/TaggerCheckNeutrino.cxx` | Primary producer: fills TaggerInfo + KineInfo |
| `clus/src/NeutrinoTaggerCosmic.cxx` | Fills `cosmict_*` features |
| `clus/src/NeutrinoTaggerNuMu.cxx` | Fills `numu_cc_*` features |
| `clus/src/NeutrinoTaggerNuE.cxx` | Fills nue BDT features (largest tagger, ~4400 lines) |
| `clus/src/NeutrinoTaggerSSM.cxx` | Fills `ssm_*` features |
| `clus/src/NeutrinoTaggerSinglePhoton.cxx` | Fills `shw_sp_*` features |
| `clus/src/NeutrinoKinematics.cxx` | `fill_kine_tree()` fills KineInfo |
| `clus/inc/WireCellClus/TrackFitting.h` (lines 188-194) | `get_tagger_info()`, `get_kine_info()` accessors |
| `root/src/UbooneNumuBDTScorer.cxx` | Numu BDT scorer (writes numu_score) |
| `root/src/UbooneNueBDTScorer.cxx` | Nue BDT scorer (writes nue_score, 30 sub-scores) |
| `root/src/UbooneMagnifyTrackingVisitor.cxx` | Reference visitor: ROOT file output pattern |
| `cfg/pgrapher/common/clus.jsonnet` | Pipeline configuration for all visitors |
| `clus/test/test-porting.bats` | Existing BATS integration tests |
| `clus/inc/WireCellClus/PatternDebugIO.h` | JSON serialization for PR state debugging |
| `clus/docs/tagger/init_tagger_info_review.md` | Reference for correct default values |
| `clus/docs/tagger/numu_bdts_xgboost_review.md` | Numu BDT implementation review |
| `clus/docs/tagger/nue_cal_bdts_xgboost_review.md` | Nue BDT implementation review |
