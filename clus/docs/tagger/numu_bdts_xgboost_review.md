# Code Review: `cal_numu_bdts_xgboost()` and `match_isFC`

**Date:** 2026-04-08  
**Reviewer:** Claude Code  
**Scope:** Logic fidelity vs prototype, bug hunting, efficiency, determinism.

---

## Files Examined

| Role | File |
|---|---|
| Toolkit BDT scorer | `root/src/UbooneNumuBDTScorer.cxx` |
| Toolkit `match_isFC` computation | `clus/src/TaggerCheckNeutrino.cxx:377-383` |
| Toolkit FC check helper | `clus/src/Clustering_Util.cxx:57-…` |
| Toolkit FC check in STM | `clus/src/TaggerCheckSTM.cxx:2146` |
| Prototype BDT scorer | `prototype_pid/src/NeutrinoID_numu_bdts.h:1-98` |
| Prototype `match_isFC` origin | `prototype_pid/apps/wire-cell-prod-nue.cxx:232-236` |

---

## Functions Reviewed

| Function | Toolkit Location | Prototype Location |
|---|---|---|
| `cal_numu_bdts_xgboost` | `UbooneNumuBDTScorer.cxx:263` | `NeutrinoID_numu_bdts.h:1` |
| `cal_numu_1_bdt` | `UbooneNumuBDTScorer.cxx:139` | `NeutrinoID_numu_bdts.h:275` |
| `cal_numu_2_bdt` | `UbooneNumuBDTScorer.cxx:192` | `NeutrinoID_numu_bdts.h:321` |
| `cal_numu_3_bdt` | `UbooneNumuBDTScorer.cxx:233` | `NeutrinoID_numu_bdts.h:353` |
| `cal_cosmict_10_bdt` | `UbooneNumuBDTScorer.cxx:92` | `NeutrinoID_numu_bdts.h:238` |
| `match_isFC` fill | `TaggerCheckNeutrino.cxx:377` | constructor arg from external tree |

---

## No bugs found. No changes made.

---

## `match_isFC` — Design and Fidelity

### Prototype design

In the prototype, `match_isFC` is **not computed inside `NeutrinoID`**. It is:

1. Computed by an earlier processing step (cluster-characterisation, before the STM step), which writes it to a ROOT `T_eval` tree as a `bool`.
2. Read back from that tree by the NuE/NuMu apps:
   ```cpp
   bool match_isFC = false;
   T_eval->SetBranchAddress("match_isFC", &match_isFC);
   ```
3. Passed as a constructor argument to `NeutrinoID`:
   ```cpp
   WCPPID::NeutrinoID *neutrino = new WCPPID::NeutrinoID(..., match_isFC, ...);
   ```
4. Stored as `float match_isFC` member in `NeutrinoID` and used directly in the BDT `AddVariable`.

The prototype's boundary-check logic that produces `match_isFC` is equivalent to the toolkit's `cluster_fc_check` function — both perform the same two-round steiner boundary check with FV, signal-processing, and dead-volume tests.

### Toolkit design

The toolkit collapses this two-step (external tree → constructor arg) into a single inline computation:

```cpp
// TaggerCheckNeutrino.cxx:380-383
if (main_cluster) {
    auto fc_result = Facade::cluster_fc_check(*main_cluster, m_dv);
    tagger_info.match_isFC = fc_result.is_fc ? 1.0f : 0.0f;
}
```

`cluster_fc_check` is the **same function** used by `TaggerCheckSTM`:

```cpp
// TaggerCheckSTM.cxx:2146
auto fc_result = Facade::cluster_fc_check(cluster, m_dv);
if (fc_result.is_fc) { return false; }  // FC cluster → not STM
```

So `match_isFC = 1` iff the cluster is fully contained — consistent between both toolkit users and with the prototype's boolean semantics (`true` → 1.0f, `false` → 0.0f). ✅

### Findings

#### ✅ FC check logic is consistent with TaggerCheckSTM

Both `TaggerCheckSTM` and `TaggerCheckNeutrino` call `Facade::cluster_fc_check` with the same cluster and DV pointer. The `is_fc` field has the same meaning in both contexts. **Consistent.**

#### ✅ FC check result stored correctly as float

The prototype declares `float match_isFC` (NeutrinoID.h:1974) and the tree stores a `bool` which ROOT auto-converts. The toolkit stores `1.0f`/`0.0f` in `TaggerInfo::match_isFC` (float, default 0). TMVA `AddVariable("match_isFC", &ti.match_isFC)` receives the correct value. ✅

#### ✅ Null guard on `main_cluster`

Toolkit wraps the FC check in `if (main_cluster)`, leaving `match_isFC` at its default `0` if no main cluster was found. Correct — a missing cluster is not FC. ✅

---

## `cal_numu_bdts_xgboost`

### Structure

1. Calls four sub-BDTs to fill intermediate scores in `ti`.
2. Registers ~74 scalar `TaggerInfo`/`KineInfo` fields with a TMVA reader.
3. Guards four fields against NaN before `EvaluateMVA`.
4. Stores `log10((1+val1)/(1-val1))` in `ti.numu_score`.

### Findings

#### ✅ Sub-BDT default values match

| Sub-BDT | Default | Prototype | Toolkit |
|---|---|---|---|
| `numu_1_score` | −0.4 | ✅ | ✅ |
| `numu_2_score` | −0.1 | ✅ | ✅ |
| `numu_3_score` | −0.2 | ✅ | ✅ |
| `cosmict_10_score` | +0.7 | ✅ | ✅ |

#### ✅ TMVA variable list matches prototype exactly

All 74 `AddVariable` calls match the prototype in name, order, and target member. Key alias preserved:

```cpp
// Training variable name differs from struct member name:
reader.AddVariable("numu_cc_3_track_length", &ti.numu_cc_3_acc_track_length);
```

This matches prototype line 14: `reader.AddVariable("numu_cc_3_track_length", &tagger_info.numu_cc_3_acc_track_length)`. ✅

#### ✅ `numu_3_score` computed but not fed into main model

Prototype computes `numu_3_score` (line 6) but does not add it to the main xgboost `AddVariable` list (only `numu_1_score` and `numu_2_score` appear). Toolkit does the same: `ti.numu_3_score` is set but not registered with the main reader. **Intentional, matches prototype.** ✅

#### ✅ NaN guards — correct placement and scope

Both prototype and toolkit apply four NaN guards **after `BookMVA` but before `EvaluateMVA`**:

```cpp
if (std::isnan(ti.cosmict_4_angle_beam)) ti.cosmict_4_angle_beam = 0;
if (std::isnan(ti.cosmict_7_angle_beam)) ti.cosmict_7_angle_beam = 0;
if (std::isnan(ti.cosmict_7_theta))      ti.cosmict_7_theta = 0;
if (std::isnan(ti.cosmict_7_phi))        ti.cosmict_7_phi = 0;
```

Only these four fields are guarded in the prototype; the toolkit does not over-guard. ✅

#### ✅ Log10 formula — `TMath::Log10` vs `std::log10`

Prototype: `val = TMath::Log10( (1+val1)/(1-val1) );`  
Toolkit: `ti.numu_score = static_cast<float>(std::log10((1.0 + val1) / (1.0 - val1)));`

`TMath::Log10` and `std::log10` are identical for positive finite arguments. The argument `(1+val1)/(1-val1)` is positive when `|val1| < 1`, which is guaranteed by TMVA's bounded output range. ✅

#### ✅ `kine_reco_Enu` local copy necessary

Prototype: `&kine_info.kine_reco_Enu` — direct address of a `float` member (mutable).  
Toolkit: `ki` is a `const KineInfo&`, so `&ki.kine_reco_Enu` would be `const float*`, which TMVA `AddVariable` rejects. The local `float kine_reco_Enu = static_cast<float>(ki.kine_reco_Enu)` provides a mutable address. Since the type is already `float` (KineInfo.h:32), `static_cast<float>` is a no-op but harmless. **Correct design.** ✅

#### ✅ Return value vs void — architectural difference

Prototype returns `float val` (which the caller stores). Toolkit is `void` and writes directly to `ti.numu_score`. The caller (the `visit()` method) already owns the `TaggerInfo&`. **Correct architectural adaptation.** ✅

---

## `cal_numu_1_bdt`

### Findings

#### ✅ `numu_cc_flag_1` — loaded but not in BDT input

Prototype (L301): declares local `float numu_cc_flag_1` and fills it from the vector, but never passes it to `reader_numu_1.AddVariable`. It is dead — the variable is unused after assignment.

Toolkit (L165): `(void)ti.numu_cc_flag_1.at(i);` — explicitly marks it as intentionally loaded-but-unused and still performs the bounds check. ✅

#### ✅ `isinf` guard on `numu_cc_1_dQ_dx_cut`

Both prototype and toolkit clamp infinite `dQ_dx_cut` to 10 before scoring. ✅

#### ✅ Iterates vector, returns maximum

Both initialize `val = -1e9` (when non-empty) and take `max(val, tmp_bdt)` over the per-segment vector. ✅

---

## `cal_numu_2_bdt`

### Findings

#### ✅ Logic fidelity

Iterates `numu_cc_2_length` vector, takes maximum. Four variables. Default −0.1. All match prototype. ✅

---

## `cal_numu_3_bdt`

### Findings

#### ✅ Logic fidelity — scalar fields, single call

Flag-3 has scalar (not per-segment-vector) TaggerInfo fields. Both prototype and toolkit pass `&ti.numu_cc_3_*` addresses directly without iteration — single `EvaluateMVA` call. ✅

---

## `cal_cosmict_10_bdt`

### Findings

#### ✅ Iterates vector, returns minimum

Both initialize `val = 1e9` (when non-empty) and take `min(val, tmp_bdt)` — the most cosmic-like element wins. ✅

#### ✅ NaN guard on `cosmict_10_angle_beam`

Both prototype and toolkit clamp NaN angle to 0 per element inside the loop. ✅

#### ✅ Default value 0.7 (most FC-like) when empty

No upstream-dirt clusters → default 0.7, which makes the event look like a neutrino candidate. Matches prototype. ✅

---

## Summary Table

| Check | Result | Action |
|---|---|---|
| `match_isFC` consistency with TaggerCheckSTM | ✅ Same `cluster_fc_check`, same semantics | — |
| `match_isFC` null guard on `main_cluster` | ✅ Improvement vs prototype | — |
| `match_isFC` float encoding (1.0/0.0) | ✅ Correct | — |
| `cal_numu_bdts_xgboost` variable list (74 vars) | ✅ Matches prototype exactly | — |
| `numu_cc_3_track_length` alias preserved | ✅ Matches prototype | — |
| `numu_3_score` not fed to main model | ✅ Intentional, matches prototype | — |
| NaN guards — 4 fields, correct placement | ✅ Matches prototype | — |
| `TMath::Log10` → `std::log10` | ✅ Equivalent | — |
| `kine_reco_Enu` local float copy | ✅ Necessary (const ref → mutable ptr) | — |
| Return void vs float | ✅ Correct architectural adaptation | — |
| `cal_numu_1_bdt` unused `numu_cc_flag_1` | ✅ Documented via `(void)` | — |
| `cal_numu_1_bdt` `isinf` guard | ✅ Matches | — |
| `cal_numu_2_bdt` logic | ✅ Matches | — |
| `cal_numu_3_bdt` scalar fields, single call | ✅ Matches | — |
| `cal_cosmict_10_bdt` minimum score, NaN guard | ✅ Matches | — |

---

## Changes Made

None. No bugs or logic divergences found.
