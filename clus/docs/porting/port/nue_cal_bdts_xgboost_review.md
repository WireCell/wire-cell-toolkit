# Code Review: `cal_bdts_xgboost` (NuE BDT scorer)

**Date:** 2026-04-09  
**Reviewer:** Claude Code  
**Scope:** Logic fidelity vs prototype, bug hunting, efficiency, determinism.

---

## Files Examined

| Role | File |
|---|---|
| Toolkit implementation | `root/src/UbooneNueBDTScorer.cxx` |
| Prototype implementation | `prototype_pid/src/NeutrinoID_nue_bdts.h:3` |
| TaggerInfo defaults | `clus/inc/WireCellClus/NeutrinoTaggerInfo.h:1397` |

---

## Functions Reviewed

| Function | Toolkit Location | Prototype Location |
|---|---|---|
| `cal_bdts_xgboost` | `UbooneNueBDTScorer.cxx:1396` | `NeutrinoID_nue_bdts.h:3` |

---

## One bug fixed.

**`nue_score` default when `br_filled != 1` — fixed from 0 to −15.**

---

## Structure

Both prototype and toolkit follow the same four-step layout:

1. Compute 15 vector sub-BDT scores.
2. Apply variable protection (clamping and NaN guards).
3. Build TMVA reader with ~160 features; call `BookMVA`.
4. Evaluate and transform via `log10((1+val)/(1-val))`.

---

## Findings

### ✅ Step 1 — 15 vector sub-BDT scores match

All 15 sub-BDT calls match in name, order, and default threshold:

| Score | Threshold | Prototype | Toolkit |
|---|---|---|---|
| `br3_3_score` | 0.30 | ✅ | ✅ |
| `br3_5_score` | 0.42 | ✅ | ✅ |
| `br3_6_score` | 0.75 | ✅ | ✅ |
| `pio_2_score` | 0.20 | ✅ | ✅ |
| `stw_2_score` | 0.70 | ✅ | ✅ |
| `stw_3_score` | 0.50 | ✅ | ✅ |
| `stw_4_score` | 0.70 | ✅ | ✅ |
| `sig_1_score` | 0.59 | ✅ | ✅ |
| `sig_2_score` | 0.55 | ✅ | ✅ |
| `lol_1_score` | 0.85 | ✅ | ✅ |
| `lol_2_score` | 0.70 | ✅ | ✅ |
| `tro_1_score` | 0.28 | ✅ | ✅ |
| `tro_2_score` | 0.35 | ✅ | ✅ |
| `tro_4_score` | 0.33 | ✅ | ✅ |
| `tro_5_score` | 0.50 | ✅ | ✅ |

---

### ✅ Step 2 — Variable protection: 4 guards match

```cpp
// Both prototype and toolkit:
if (mip_min_dis > 1000)                           mip_min_dis = 1000;
if (mip_quality_shortest_length > 1000)           mip_quality_shortest_length = 1000;
if (std::isnan(mip_quality_shortest_angle))        mip_quality_shortest_angle = 0;
if (std::isnan(stem_dir_ratio))                   stem_dir_ratio = 1.0;
```

✅ All four guards match. Only these four fields are protected — no over-guarding.

**Protection placement:** Prototype applies guards after `BookMVA` but before `EvaluateMVA`.
Toolkit applies them before building the reader. Since TMVA reads variable values via registered
pointers only at `EvaluateMVA` time, the ordering relative to `AddVariable` / `BookMVA` has no
effect on results. ✅

---

### ✅ Step 3 — TMVA variable list: ~160 variables match

All variables match in name, order, and target `TaggerInfo` / `KineInfo` member.

Key structural notes:

#### ✅ `kine_reco_Enu` local mutable copy

Prototype: `&kine_info.kine_reco_Enu` (mutable member).  
Toolkit: `float nue_kine_reco_Enu = static_cast<float>(ki.kine_reco_Enu);` — local mutable copy
because `ki` is a `const KineInfo&`. Same approach as `cal_numu_bdts_xgboost`. ✅

#### ✅ `match_isFC` consistent with TaggerCheckNeutrino

Both prototype (`float match_isFC` NeutrinoID member) and toolkit (`ti.match_isFC` TaggerInfo
member) are mutable floats, holding 1.0 (FC) or 0.0 (not FC). ✅

#### ✅ Two "naming issue" aliases preserved

Prototype notes two fields where the TMVA training name differs from the struct member name or
field placement is unexpected:

```cpp
// Prototype comment: "naming issue"
reader.AddVariable("hol_1_flag_all_shower", &tagger_info.hol_1_flag_all_shower);
// ... (field appears after lol_3 block rather than with other hol_1 fields)

reader.AddVariable("brm_acc_direct_length", &tagger_info.brm_acc_direct_length);
// ... (field appears after lem block rather than with other brm fields)

reader.AddVariable("br4_1_n_main_segs", &tagger_info.br4_1_n_main_segs);
// ... (field appears after tro_3 block rather than with other br4_1 fields)
```

All three are reproduced in the toolkit in the same out-of-group positions. ✅

---

### ✅ Step 4 — Log10 transform — `TMath::Log10` vs `std::log10`

Prototype: `val = TMath::Log10((1+val1)/(1-val1));`  
Toolkit: `ti.nue_score = static_cast<float>(std::log10((1.0 + val1) / (1.0 - val1)));`

`TMath::Log10` and `std::log10` are identical for positive finite arguments. The argument is
positive when `|val1| < 1`, guaranteed by TMVA's bounded output. ✅

---

### 🐛 `nue_score` default when `br_filled != 1` — **fixed**

**Prototype:**
```cpp
float default_val = -15;  // background like
double val = 0;
// ...
if (tagger_info.br_filled == 1) {
    double val1 = reader.EvaluateMVA("MyBDT");
    val = TMath::Log10((1+val1)/(1-val1));
} else {
    val = default_val;   // returns -15
}
return val;             // caller stores: tagger_info.nue_score = cal_bdts_xgboost()
```

**Toolkit (before fix):**
```cpp
if (ti.br_filled == 1) {
    float val1 = reader.EvaluateMVA("MyBDT");
    ti.nue_score = static_cast<float>(std::log10((1.0 + val1) / (1.0 - val1)));
}
// else: ti.nue_score stays at TaggerInfo default = 0
```

**Issue:** `TaggerInfo::nue_score` initialises to `0` (NeutrinoTaggerInfo.h:1397). The prototype
returns `−15` (declared as `// background like`) when `br_filled != 1`. A neutral score of 0
instead of −15 would make events with missing `br_filled` appear less background-like than they
should.

**Fix applied** (`UbooneNueBDTScorer.cxx:1692`):
```cpp
} else {
    ti.nue_score = -15.f;  // background-like default; matches prototype default_val = -15
}
```

---

## Summary Table

| Check | Result | Action |
|---|---|---|
| 15 vector sub-BDT scores — names, order, thresholds | ✅ All match | — |
| Variable protection — 4 guards, correct fields | ✅ Match | — |
| Protection placement vs `BookMVA` — order irrelevant | ✅ No effect on output | — |
| TMVA variable list — ~160 vars, name + order + member | ✅ Match | — |
| `kine_reco_Enu` local mutable copy | ✅ Correct (const ref → mutable ptr) | — |
| `match_isFC` float encoding | ✅ Correct | — |
| Three "naming issue" out-of-group fields preserved | ✅ Match | — |
| `TMath::Log10` → `std::log10` | ✅ Equivalent | — |
| `nue_score` default when `br_filled != 1` | 🐛 Toolkit had 0; prototype uses −15 | Fixed |

---

## Changes Made

**`root/src/UbooneNueBDTScorer.cxx:1692`** — added `else { ti.nue_score = -15.f; }` so that the
background-like default matches the prototype when `br_filled` is not set.
