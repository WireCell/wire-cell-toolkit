# Charge–Light χ² and Flags: QLMatching (toolkit) vs ToyMatching (MicroBooNE prototype)

## 1. Purpose & scope

This is a **guideline and map** for the planned update of the χ² / light-matching logic in the
toolkit `QLMatching`, following the MicroBooNE prototype (`ToyMatching` / `FlashTPCBundle`) as a
reference.

Two facts motivate it:

1. In the current toolkit, the geometric/topology flags
   (`flag_close_to_PMT`, `flag_at_x_boundary`, `flag_spec_end`, `flag_window_truncated`) are
   **computed but inert** — none of them feed the per-bundle χ² or the LASSO weights.
2. The χ² construction itself differs from MicroBooNE, both in the error denominator and in the
   `flag_high_consistent` acceptance logic.

The prototype is **MicroBooNE**: 32 PMTs, a single anode plane (`low_x_cut = 0`),
`high_x_cut ≈ 256 cm`. The toolkit target is **SBND**: ~120 channels across **two drift volumes**
(two anodes, a shared central cathode). Many prototype constants are absolute-PE or single-anode
quantities and **cannot be copied** — they must be re-derived for SBND. This document tags each one
(§9).

This document changes no code. It covers the five requested topics: (1) χ² construction,
(2) flag usage, (3) regularization strength, (4) light errors, (5) χ² in the first vs second round —
plus how the prototype reconstructs the per-channel PE error upstream of matching (§7.1) and a
complete reference for the toolkit's jsonnet config parameters (§11).

## 2. File map

**Toolkit (this repo, `match/`)**

| What | File:line |
|---|---|
| Per-bundle χ² + `flag_high_consistent` | `match/src/TimingTPCBundle.cxx:186-201` |
| LASSO round 1 (setup, weights, errors) | `match/src/QLMatching.cxx:670-751` |
| LASSO round 2 (KS shape penalty, best-per-cluster) | `match/src/QLMatching.cxx:753-849` |
| Where endpoint flags are SET | `match/src/QLMatching.cxx:1032-1155` |
| PE error model (floor/frac/knee) | `match/inc/WireCellMatch/Opflash.h:15-19`, `match/src/Opflash.cxx:8-29` |
| Config defaults | `match/inc/WireCellMatch/QLMatching.h:84,120-148` |

**Prototype (MicroBooNE, separate checkout under `/nfs/data/1/xqian/prototype-dev/wire-cell/`)**

| What | File:line |
|---|---|
| Per-bundle χ², flag use, consistency ladder | `data/src/FlashTPCBundle.cxx:446-603` |
| Per-PMT hit PE + PE_err reconstruction | `data/src/COphit.cxx:33-61`, `data/inc/WCPData/COphit.h` |
| Flash PE_err aggregation (per type) | `data/src/Opflash.cxx:7-119`, `data/inc/WCPData/Opflash.h` |
| Flags SET | `2dtoy/src/ToyMatching.cxx:279-290` |
| LASSO round 1 | `2dtoy/src/ToyMatching.cxx:1377-1548` |
| LASSO round 2 | `2dtoy/src/ToyMatching.cxx:1552-1663` |

## 3. Per-bundle χ² construction

This is the χ² stored on each bundle (one flash × one cluster), used to gate
`flag_high_consistent`. It is **separate** from the LASSO normal equations (§8); the two use
different error models in the prototype but the same one in the toolkit (§7).

### Toolkit — `TimingTPCBundle.cxx:189-195`

```cpp
for (int j = 0; j < m_nchan; ++j) {
    if (opdet_mask[j] == 0) continue;
    ++nvalidopdets;
    if (pe[j] < m_qp.pe_ndf_knee && pred_pe[j] < m_qp.pe_ndf_knee) { /* no-op */ }
    else ndf++;
    chi2 += std::pow(pred_pe[j] - pe[j], 2) / (pe[j] + std::pow(pe_err[j], 2));
}
```

```
        Σ_j  (pred_j − pe_j)²
χ² =   ─────────────────────────          (j over masked opdets only)
          pe_j  +  pe_err_j²
```

The denominator is **Poisson variance** (`pe_j` ≈ measured counts) **plus systematic variance**
(`pe_err_j²`). NDF counts opdets where `pe` **or** `pred` ≥ `pe_ndf_knee` (=1.0). No flags, no
outlier handling, no cosmic veto.

### Prototype — `FlashTPCBundle.cxx:459-502`

```cpp
for (int j=0;j!=32;j++){                                  // cosmic veto on the prediction:
  h1->SetBinContent(j+1,pe[j]);
  if ((pred_pe[j] < cos_pe_low[j] ||
       (pred_pe[j] < cos_pe_mid[j]*1.1 && pe[j]==0)) && flash->get_type()==1){
    pred_pe[j] = 0;                                       // zero out channels too dim to be cosmic
  }
  h2->SetBinContent(j+1,pred_pe[j]);
}
...
for (int j=0;j!=32;j++){
  double cur_chi2 = 0;
  if (flag_close_to_PMT){
    if (pe[j]-pred_pe[j]>350 && pe[j]>pred_pe[j]*1.3){    // big measured excess near the PMTs
      cur_chi2 = pow(pred_pe[j]-pe[j],2)/(pow(pe_err[j],2)+pow(pe[j]*0.5,2));  // inflate denom
    }else{
      cur_chi2 = pow(pred_pe[j]-pe[j],2)/pow(pe_err[j],2);
    }
  }else{
    cur_chi2 = pow(pred_pe[j]-pe[j],2)/pow(pe_err[j],2);
  }
  chi2 += cur_chi2;
  if (cur_chi2 > max_chi2){ max_chi2 = cur_chi2; max_bin = j; }
  if (pe[j]==0 && pred_pe[j]==0){}else{ ndf++; }
}
if (pe[max_bin] == 0 && pred_pe[max_bin]>0)              // tolerate one inefficient PMT
  chi2 -= max_chi2-1;
```

```
        Σ_j  (pred_j − pe_j)²
χ² =   ─────────────────────                (default denom)
            pe_err_j²
```

with three behaviors the toolkit lacks:

- **Cosmic veto**: for cosmic-type flashes (`get_type()==1`), predicted PE below the cosmic
  expectation (`cos_pe_low`, `cos_pe_mid`) is zeroed before χ². This is a per-channel mask on the
  *prediction*, externally supplied.
- **`flag_close_to_PMT` error inflation**: a channel with a large *measured* excess over the
  prediction (`pe−pred > 350` **and** `pe > 1.3·pred`) gets its denominator widened by `(pe·0.5)²`,
  softening near-PMT over-response.
- **One-inefficient-PMT subtraction**: the single worst channel, if it is `pe==0 && pred>0`, is
  effectively removed (`chi² −= max_chi² − 1`).

### What this means for the port

The two `(pred − pe)²/denom` denominators **both contain Poisson + systematic** — they just fold
the Poisson term in at *different stages*:

- **Toolkit**: `pe_err` is systematic-only (floor/frac/knee), and the Poisson variance is added
  *here* in the χ², as the explicit `+ pe` term (`pe + pe_err²`).
- **Prototype**: `pe_err` already **has the Poisson variance baked in** at flash-reconstruction
  time (the `2·PE` / `PE + 2·noise` statistical terms in COphit/Opflash — see §7.1), so the χ²
  needs no extra `+ pe`; `pe_err²` alone is Poisson+systematic.

So this is **not** "toolkit has Poisson, prototype doesn't." The real porting hazard is
**double-counting**: if the port adopts the prototype's reconstructed `pe_err` (Poisson inside)
*and* keeps the toolkit's `+ pe`, Poisson is counted twice. The open decision (§10) is therefore
*where* Poisson is folded in, not *whether*. Re-deriving the `350`-PE inflation threshold for SBND
is mandatory regardless (§9).

## 4. `flag_high_consistent` consistency gate

`flag_high_consistent` is a pre-LASSO quality flag: a high-consistency bundle survives, and within
its cluster the *other* (inconsistent) bundles are dropped before the fit
(`QLMatching.cxx` filtering around the consistent-bundle pass).

### Toolkit — `TimingTPCBundle.cxx:197-200` (single condition)

```cpp
flag_high_consistent = false;
if (ks_dis < m_qp.highconsist_ks_max && ndf >= m_qp.highconsist_min_ndf &&
    chi2 < ndf * nvalidopdets)
    flag_high_consistent = true;
```

One branch: `ks < 0.06`, `ndf ≥ 3`, `chi² < ndf · nvalidopdets`. **No boundary relaxation.**

### Prototype — `FlashTPCBundle.cxx:504-522` (multi-branch ladder)

```cpp
if (ks_dis < 0.06 && ndf >=3 && chi2 < ndf * 36){ flag_high_consistent = true; }
else if (ks_dis<0.05 && ndf >=6 && chi2 < ndf * 45){ ... }
else if (ks_dis < 0.12 && ndf >=3 && chi2 < ndf * 25){ ... }
else if (flag_at_x_boundary && ndf >=2 && chi2 < 9 * ndf && ks_dis < 0.12){ ... }   // boundary
else if (flag_at_x_boundary && ndf >=1 && chi2 < 3 * ndf && ks_dis < 0.12){ ... }   // boundary
else if (chi2 < 4 * ndf && ndf >=3 && ks_dis < 0.15){ ... }
else if (chi2 < 1.5 * ndf && ks_dis < 0.2 && ndf >= 3){ ... }
else if (ks_dis < 0.12 && ndf >=5 && chi2 < ndf * 55 && flag_close_to_PMT){ ... }   // near-PMT
else if (ks_dis < 0.14 && ndf >=3 && chi2< ndf * 6){ ... }
```

Several branches, including **boundary-relaxed** ones: `flag_at_x_boundary` permits `χ² < 9·ndf`
(ndf≥2) or `χ² < 3·ndf` (ndf≥1); `flag_close_to_PMT` permits `χ² < 55·ndf` (ndf≥5). Note the
toolkit's `chi2 < ndf·nvalidopdets` cut (nvalidopdets ≈ 100+ for SBND) is structurally *much looser*
than the prototype's `ndf · 36`; this differs because SBND has many more channels than 32.

A further prototype gate beyond `flag_high_consistent` (`FlashTPCBundle.cxx:547-602`) uses the
"fired fraction" (`nfired/ntot` of channels with `pred > 0.33`) to set `flag_potential_bad_match`
and to reject bundles outright — also boundary-aware. The toolkit has no equivalent.

## 5. Flag usage

| Flag | Toolkit: where SET | Prototype: how USED | Toolkit status |
|---|---|---|---|
| `flag_close_to_PMT` | `QLMatching.cxx:1128` (anode-end, within `m_anode_ext2`) | χ² denom inflation `+ (pe·0.5)²` for big excess (`FlashTPCBundle.cxx:480-485`); relaxed consistency `χ²<55·ndf` (`:518`); LASSO weight 1.0→0.2 via paired `flag_at_x_boundary` | **inert** (debug log only) |
| `flag_at_x_boundary` | `QLMatching.cxx:1129` (anode), `:1155` (cathode, CPA 3-D fiducial) | relaxed consistency `χ²<9·ndf`/`3·ndf` (`FlashTPCBundle.cxx:510-513`); LASSO weight 1.0→0.2 (`ToyMatching.cxx:1437-1441`, rd1; `1603-1607`, rd2) | **inert** |
| `flag_spec_end` | `QLMatching.cxx:1080,1104,1124` | informational; tail-alignment bookkeeping | **inert** |
| `flag_window_truncated` | `QLMatching.cxx:1032-1037` (slice within `m_window_edge_ticks` of readout window edge) | *no direct prototype analog* — see note below | **inert** |
| `flag_potential_bad_match` | (toolkit sets via `ks==1` drop, `QLMatching.cxx:550`) | fired-fraction reject, boundary-aware (`FlashTPCBundle.cxx:577-599`) | partial |

Prototype flag-set site (`ToyMatching.cxx:281-290`): `flag_close_to_PMT` + `flag_at_x_boundary` are
set when the cluster's **anode-end** position falls in the low-x boundary window; `flag_at_x_boundary`
alone when the **cathode-end** falls in the high-x window. Three mechanisms downstream:
(a) χ² error inflation, (b) LASSO weight 1.0→0.2, (c) relaxed acceptance thresholds.

**`flag_window_truncated` (user note #2).** This toolkit flag marks a bundle whose leading/trailing
time slice sits at the raw readout-window edge — its light may be partially outside the recorded
window, so the *measured* PE is an underestimate. It is the toolkit's near-analog of
`flag_close_to_PMT`'s intent (tolerate measured-vs-predicted mismatch), but triggered by the **time**
window rather than the **drift** boundary, and it should carry **slightly less tolerance** than
`flag_close_to_PMT`: window truncation removes part of the signal (one-sided), whereas near-PMT
over-response adds signal. A reasonable port would let `flag_window_truncated` relax the high-χ²
cut and/or down-weight in LASSO, but with a tighter multiplier than the `flag_close_to_PMT`
branch — to be tuned (§9, §10).

## 6. Regularization strength & knob cross-map

Both implementations solve a non-negative LASSO `min ‖y − Xs‖² + λ‖s‖₁` in normal-equation form,
with **λ = 0.1** in both rounds.

| Concept | Toolkit name / value | Prototype name / value |
|---|---|---|
| LASSO λ | `m_lasso_lambda = 0.1` | `lambda = 0.1` |
| Charge/track constraint (Σ strengths per cluster → 1) | `delta_charge = 0.01` (row `1/δ`) | `delta_track = 0.01` |
| Flash background absorption (rd1 only) | `delta_light = 0.025` | `delta_flash = 0.025` |
| KS shape penalty scale (rd2 only) | `delta_shape = 0.01` | *(none — toolkit-specific)* |
| Background column weight (rd1) | `bkg_weight = 0.5` | *(implicit in flash DOF)* |
| Per-bundle weight floor / knee | `pe_mismatch_floor = 0.3`, `pe_mismatch_knee = 0.3` | (no floor; weight=mismatch directly) |
| Boundary down-weight | *(not applied — flags inert)* | `0.2` for `flag_at_x_boundary` |
| Strength keep-threshold | `m_strength_cutoff = 0.05` | `0.05` |

The λ value and `delta_*` constants are **dimensionless ratios** (the prototype comment says "the
coefficient is all around 1"), so they are good candidates to **port as-is**; the boundary
down-weight `0.2` is likewise dimensionless. See §9.

## 7. Light error model — two distinct sites

### 7.1 Upstream: where the prototype's `pe_err` is born (flash reconstruction)

(User note #1.) Everything downstream — the χ² denominator (§3) and the LASSO error (§7.2) — reads
the flash's per-channel `flash->get_PE_err(k)`. That array is **not** computed in the matcher; it is
filled during **optical flash reconstruction**, in `COphit` (per-PMT hit) and `Opflash` (the flash
that aggregates hits). The error model **depends on the flash type**, set by which `Opflash`
constructor built it:

| Flash type | Built by | Source | File:line |
|---|---|---|---|
| `type==1` (cosmic / L1) | `Opflash(COphitSelection&)` | per-PMT `COphit` integrals | `data/src/Opflash.cxx:24-52` |
| `type==2` (beam) | `Opflash(TH1F**, …)` | summed waveform histogram | `data/src/Opflash.cxx:54-119` |
| (replay) | `Opflash(int type, …, temp_PE_err, …)` | passthrough of stored arrays | `data/src/Opflash.cxx:7-22` |

**Per-PMT hit (`COphit.cxx:33-61`).** A hit's PE comes from the waveform integral over gain
(`PE = integral/gain · 2`, the ×2 restoring the 0.6 µs window). Its error is built in quadrature:

```cpp
// good_baseline (COphit.cxx:38-42):
PE_err = sqrt( pow(PE*gain_err/gain, 2)      // gain calibration uncertainty
             + 2*PE                          // STATISTICAL (Poisson) — Var(2N)=4N=2·PE after ×2 scaling
             + pow(40/sqrt(3.)/gain*2, 2)    // baseline: 1/√3 ADC over 40 ticks, ×2 window
             + pow(PE*0.03, 2) );            // 3% relative (baseline-guess) systematic
// bad baseline (COphit.cxx:60):  PE_err = 2*PE;   // flat 200% — channel is untrustworthy
```

**Cosmic flash (`Opflash.cxx:24-52`).** Each channel is initialized to a floor
`PE_err = 6.4 PE` (`= 11/√3`, the per-channel noise quantum) and then **overwritten by the
COphit's `PE_err`** for fired channels. So the cosmic path inherits the COphit quadrature above.

**Beam flash (`Opflash.cxx:54-119`).** PE is summed from the waveform histogram (sub-0.2-PE bins
zeroed, a scaled noise pedestal subtracted). The error starts from a `0.2 PE` floor and adds its
own statistical + systematic terms:

```cpp
// Opflash.cxx:103-106:
PE_err[i] = sqrt( pow(PE_err[i], 2)            // 0.2 PE base
                + (PE[i] + pe_noise_scaled*2)  // STATISTICAL (Poisson) term
                + pow(PE[i]*0.02, 2) );        // 2% base systematic
```

**The key takeaway for the port** (and the correction to the old "pure systematic" framing below):
in **both** prototype flash types the reconstructed `pe_err` already includes a **Poisson/statistical
term** (`2·PE` for cosmic, `PE + 2·noise` for beam). The toolkit's `pe_err` (the floor/frac/knee
struct, §7.2) carries **systematic only**, and Poisson is added later in the χ²/LASSO via the
explicit `+ pe`. Both reach Poisson+systematic; only the *stage* differs — hence the
double-counting hazard in §3/§10. The prototype's absolute-PE constants here (`6.4`, `0.2` floors;
`40/√3` baseline; `2·PE`, `0.03`, `0.02` fractions) are MicroBooNE gain/noise quantities and are
**RE-TUNE for SBND** candidates (§9). SBND flashes arrive pre-reconstructed via the canonical PCs,
so the toolkit does not re-run this code — but if SBND ever needs a Poisson-in-`pe_err` model to
match the prototype, this is the formula to reproduce.

### 7.2 In the matcher — two distinct sites

The error denominator is then used in **two different places**, and the two implementations treat
them differently. **Do not conflate them.**

| Site | Toolkit | Prototype |
|---|---|---|
| Per-bundle χ² (`examine_bundle`) | `sqrt(pe + pe_err²)` → denom `pe + pe_err²` (`TimingTPCBundle.cxx:194`) | `pe_err²` (raw), `flag_close_to_PMT` may add `(pe·0.5)²` (`FlashTPCBundle.cxx:482-487`) |
| LASSO normal equations | **same** `sqrt(pe + pe_err²)` (`QLMatching.cxx:698,707,782,791`) | **different**: `sqrt((pe_err·f2)² + (pe·f1)² + (pe·rel_light_yield_err)²)` (`ToyMatching.cxx:1431`, rd1; `1585/1597`, rd2) |

Toolkit base `pe_err` (the floor/frac/knee model, `Opflash.h:15-19`, `Opflash.cxx:24`):

```cpp
struct PEErr { double floor = 0.3; double frac = 0.3; double knee = 1.0; };
PE_err[i] = (PE[i] < knee) ? floor : frac * PE[i];     // 0.3 floor, else 30% of PE
```

Prototype LASSO error (`ToyMatching.cxx:1431`):

```cpp
double pe_err = sqrt(pow(flash->get_PE_err(k)*fudge_factor2,2)   // f2 = 1.0
                   + pow(pe*fudge_factor1,2)                     // f1 = 0.06 (rd1) / 0.05 (rd2)
                   + pow(pe*rel_light_yield_err,2));             // run-dependent yield calib
```

Key differences:

- The toolkit applies **one consistent error model** (`sqrt(pe + pe_err²)`) at *both* sites; the
  prototype uses the **reconstructed `pe_err`** (Poisson already inside, §7.1) in the χ², but a
  **quadrature-augmented** error in the LASSO that adds *further* fractional terms on top.
- The toolkit has **no `rel_light_yield_err`** term (run-dependent light-yield calibration
  uncertainty). For SBND MC this may be unnecessary; for data it is a deliberate decision (§9).
- The Poisson term lives in **different places**: the toolkit adds it explicitly as `+ pe` in the
  denominator; the prototype folds it into `pe_err` upstream (§7.1). The prototype's extra LASSO
  `(pe·f1)²` fractional term (f1 = 0.06/0.05) has no toolkit counterpart — that one *is* a genuine
  difference, not just a relocation.

## 8. First vs second round

### Toolkit (`QLMatching.cxx:670-751`, `753-849`)

| | Round 1 | Round 2 |
|---|---|---|
| Columns | `nbundle + nflash` (one background light column per flash, `PF = 1/delta_light`) | `nbundle` only (no background) |
| Per-bundle weight | PE-total mismatch: `|pred_tot − meas_tot|/meas_tot` if `> 0.3·meas_tot`, else floor `0.3` (`:713-715`) | same base **plus KS shape penalty** `base + delta_shape·nopdet·ks_dis/lambda` (`:798-801`) |
| Background weight | `bkg_weight = 0.5` (`:721`) | — |
| Post-fit | drop bundles with strength ≤ `0.05` (`:743`) | drop ≤ `0.05`; then keep **best flash per cluster** by strength (`:837-846`) |

### Prototype (`ToyMatching.cxx:1377-1548`, `1552-1663`)

| | Round 1 | Round 2 |
|---|---|---|
| Columns | `num_unknowns + nflash` ("flash-alone" DOF, `RF = 1/delta_flash`) | `num_unknowns` only (flash DOF dropped) |
| Error fudge | `f1 = 0.06` | `f1 = 0.05` |
| Weight | `0.2` if `flag_at_x_boundary` else `1.0` (`:1437-1441`) | `0.2` if `flag_at_x_boundary` else `1.0` (`:1603-1607`) |
| Post-fit | drop bundles with `beta < 0.05` (`:1518`) | keep **best beta per cluster** (`:1648-1663`) |

Structural takeaways:

- Both have the **same two-round skeleton**: round 1 has an extra per-flash DOF / background
  column that round 2 removes; both prune at strength `0.05`; both finish best-per-cluster.
- The **toolkit's KS shape penalty in round 2** (`delta_shape·nopdet·ks_dis/lambda`) is
  **toolkit-specific** — the prototype has no shape term in its LASSO weights (it enforces shape
  earlier, via the `ks_dis` cuts in `flag_high_consistent`).
- The **prototype's boundary down-weight (`0.2`)** is **absent in the toolkit** (flags inert) — this
  is the single most direct thing the port should add, in *both* rounds.
- The prototype lowers `f1` from 0.06→0.05 between rounds; the toolkit has no per-round error change.

## 9. SBND vs MicroBooNE: port-as-is vs re-tune

The flags are **already computed but inert** in the toolkit; the prototype shows the χ²/weight
pathways they *should* feed. Tomorrow's job is **wiring existing toolkit flags into those pathways
and re-tuning the SBND-sensitive constants** — not building flags from scratch. The SBND cathode
3-D CPA fiducial that sets `flag_at_x_boundary` already exists and is inert (see memory
`project_ql_cathode_fiducial`).

Geometry difference to keep in mind throughout: prototype has **one anode** (`low_x_cut = 0`,
`high_x_cut ≈ 256 cm`) and **32 PMTs**; SBND has **two anodes + a shared central cathode** and
**~120 channels**. So `flag_close_to_PMT` (anode-side, both anodes) and `flag_at_x_boundary`
(cathode-side) map onto a two-sided drift geometry, and the cathode is in the *middle*, not at one
end.

| Constant | Prototype value | Port verdict | Reason |
|---|---|---|---|
| LASSO `λ` | 0.1 | **as-is** | dimensionless; toolkit already 0.1 |
| `delta_track`/`delta_charge` | 0.01 | **as-is** | dimensionless ratio; matches toolkit |
| `delta_flash`/`delta_light` | 0.025 | **as-is** | dimensionless; matches toolkit |
| boundary down-weight | 0.2 | **as-is (try first)** | dimensionless; re-tune only if SBND boundary purity differs |
| strength cutoff | 0.05 | **as-is** | dimensionless; matches toolkit |
| close-to-PMT excess `350` PE | 350 | **RE-TUNE** | absolute PE — depends on SBND light yield + per-channel response |
| close-to-PMT ratio `1.3×`, inflation `0.5` | 1.3 / 0.5 | **review** | ratios; likely transferable but validate against SBND near-anode response |
| fired threshold `0.33` PE, `pred>1.0` | 0.33 / 1.0 | **RE-TUNE** | absolute PE thresholds — scale with light yield |
| consistency `χ² < ndf·N` multipliers (36/45/25/9/3/55/6) | per branch | **RE-TUNE** | calibrated to 32-PMT χ² scale; SBND has ~120 channels |
| KS cuts (0.05–0.20) | per branch | **review** | KS is shape-normalized → channel-count-robust, but re-validate |
| `f1` 0.06/0.05, `f2` 1.0 | quadrature | **DECISION** | toolkit currently has neither; choose error model (§7) |
| `rel_light_yield_err` | run-dependent | **DECISION** | data-calibration term; may be 0 for SBND MC |
| cosmic veto `cos_pe_low/cos_pe_mid` | external arrays | **DECISION** | toolkit has no cosmic-prediction veto; needs SBND cosmic model if adopted |
| one-inefficient-PMT `chi2 -= max_chi2-1` | — | **DECISION** | reasonable with 32 PMTs; reconsider weight with ~120 SBND channels |

## 10. Open decisions for tomorrow (checklist)

- [ ] **χ² denominator**: keep toolkit `pe + pe_err²` (Poisson added in the χ²), or adopt the
      prototype's reconstructed `pe_err` (Poisson folded in upstream, §7.1) + LASSO quadrature
      `f1`/`yield`? **Do not do both** — that double-counts Poisson. Decide *where* it is folded in,
      not *whether*. (§3, §7.1, §7.2)
- [ ] **Wire flags into LASSO weights**: add the `0.2` boundary down-weight in both rounds for
      `flag_at_x_boundary` (and decide whether `flag_close_to_PMT` shares it). (§5, §8)
- [ ] **Wire flags into the consistency gate**: replace the single toolkit branch with a
      boundary-relaxed ladder; re-derive the `χ² < ndf·N` multipliers for ~120 SBND channels. (§4, §9)
- [ ] **`flag_window_truncated`**: relax χ² / down-weight with *tighter* tolerance than
      `flag_close_to_PMT` (one-sided signal loss). (§5)
- [ ] **`flag_close_to_PMT` error inflation**: re-tune the `350`-PE excess threshold for SBND. (§3, §9)
- [ ] **Per-round error**: adopt the prototype's `f1` 0.06→0.05 step, or keep one model? (§7, §8)
- [ ] **Optional prototype features**: decide on cosmic veto, one-inefficient-PMT subtraction,
      `rel_light_yield_err`. (§3, §7)
- [ ] Keep all changes **jsonnet-togglable, default OFF**, so existing production configs stay
      bit-identical (project convention `feedback_toggleable_behavior_changes`).

## 11. QLMatching config-file parameter reference

(User note #2.) Complete list of the jsonnet keys consumed by `QLMatching::configure`
(`match/src/QLMatching.cxx`), with defaults from `QLMatching.h:55-155`. The **jsonnet key** column
is the name to use in a config; the **§§** column points to where in *this* document the parameter's
role in the χ²/flags/LASSO is analyzed (the cross-maps in §6 and §9 are the *tuning* view of the
same numbers — this table is the *complete* view, including knobs unrelated to χ²). Defaults equal
the historical hard-coded literals (bit-identical if omitted) except the §A cushions, which follow
the MicroBooNE convention and intentionally differ from the old SBND literals.

### 11.1 χ² / flag / LASSO parameters (the subject of this document)

**§C — LASSO weights** (see §6, §8)

| jsonnet key | default | meaning |
|---|---|---|
| `lasso_lambda` | `0.1` | LASSO λ (‖s‖₁ penalty), both rounds |
| `delta_charge` | `0.01` | per-cluster charge/track constraint (row `1/δ`) |
| `delta_light` | `0.025` | flash background-absorption DOF (round 1) |
| `delta_shape` | `0.01` | KS shape-penalty scale (round 2, toolkit-specific) |
| `bkg_weight` | `0.5` | background-column weight (round 1) |
| `pe_mismatch_knee` | `0.3` | PE-mismatch weight knee (fraction of measured PE) |
| `pe_mismatch_floor` | `0.3` | PE-mismatch weight floor |
| `strength_cutoff` | `0.05` | drop bundles with LASSO strength below this, each round |

**§G — flash PE-error model** (forwarded to `Opflash`; see §3, §7)

| jsonnet key | default | meaning |
|---|---|---|
| `pe_err_floor` | `0.3` | systematic `pe_err` floor when `PE < knee` |
| `pe_err_frac` | `0.3` | systematic `pe_err` = `frac·PE` when `PE ≥ knee` |
| `pe_err_knee` | `1.0` | PE below which the floor applies |
| `flash_pe_threshold` | `0.0` | Opflash "fired" channel threshold (PE) |

**§F — bundle-quality thresholds** (forwarded to `TimingTPCBundle`; see §3, §4)

| jsonnet key | default | meaning |
|---|---|---|
| `bundle_ks_merge_max` | `0.2` | KS ceiling for bundle add/merge |
| `bundle_chi2ndf_merge_max` | `20` | χ²/ndf ceiling for bundle add/merge |
| `bundle_addmerge_exponent` | `0.8` | exponent in the add/merge acceptance |
| `highconsist_ks_max` | `0.06` | `flag_high_consistent` KS cut (§4) |
| `highconsist_min_ndf` | `3` | `flag_high_consistent` min ndf (§4) |
| `bundle_pe_ndf_knee` | `1.0` | PE below which a channel does not count toward ndf (§3) |

**§A — active-volume cushions** (set the endpoint/boundary flag windows; see §5, §9). Signed, in the
per-TPC anode→cathode drift coordinate u (u=0 at the anode/PMT plane).

| jsonnet key | default | meaning |
|---|---|---|
| `anode_ext1` | `-2.0 cm` | PE-inclusion window edge below the anode (`low_x_cut_ext1`) |
| `anode_ext2` | `4.0 cm` | anode flag-window outer edge — sets `flag_close_to_PMT` (`low_x_cut_ext2`) |
| `cathode_ext1` | `1.2 cm` | PE-inclusion window edge beyond the cathode (`high_x_cut_ext1`) |
| `cathode_ext2` | `-2.0 cm` | cathode flag-window inner edge — sets `flag_at_x_boundary` (`high_x_cut_ext2`) |
| `y_cushion` | `0.0 cm` | signed inward(+)/outward(−) shift of each \|y\| edge |
| `z_cushion` | `0.0 cm` | signed inward(+)/outward(−) shift of each z edge |

**§E — out-of-beam QA cuts**

| jsonnet key | default | meaning |
|---|---|---|
| `outbeam_ks_max` | `0.2` | KS ceiling for an out-of-beam match |
| `outbeam_chi2ndf_max` | `20` | χ²/ndf ceiling for an out-of-beam match |
| `outbeam_pe_frac` | `0.5` | out-of-beam total-PE mismatch fraction |

**§D — pre-selection / bad-match gate**

| jsonnet key | default | meaning |
|---|---|---|
| `mc_saturation_pe` | `5000` | total-flash-PE trigger for the MC saturated-PMT mask |

**§H — raw readout-window truncation flag** (always computed; sets `flag_window_truncated`, §5)

| jsonnet key | default | meaning |
|---|---|---|
| `readout_window_ticks` | `3427` | window end in raw ticks (SBND `daq.nticks`) |
| `window_edge_ticks` | `4` | edge-proximity threshold (~one live slice) |

### 11.2 Geometry / IO / model parameters (not χ²- or flag-related)

Listed for completeness; these select the detector, the optical model, and the time/PE
pre-selection — they do not enter the χ² or the flags.

| jsonnet key | default | meaning |
|---|---|---|
| `anode` | *(tn, required)* | `IAnodePlane` type-name |
| `detector_volumes` | *(tn, required)* | `IDetectorVolumes` type-name (active-volume bounds) |
| `cathode_fiducial` | *(tn, optional)* | `IFiducial` for the SBND cathode 3-D CPA flag (`project_ql_cathode_fiducial`) |
| `semimodel_file` | `sbnd/photodet/semi-analytical-sbnd.json` | VUVHits/VISHits/geometry/OpDet table |
| `nchan` | `312` | per-flash PE-vector length |
| `VUVEfficiency` / `VISEfficiency` | 312-entry SBND arrays | per-OpDet VUV/VIS efficiency |
| `QtoL` | `0.5` | charge→light scale |
| `drift_speed` | `1.563 mm/µs` | drift speed for the flash-T0 x-offset (pass `params.lar.drift_speed`) |
| `pmts` | `true` | apply the OpDet-type mask |
| `active_opdet_types` | `[1]` | OpDet types kept by the mask (1 = dome PMT) |
| `ch_mask` | `[]` | per-channel disable list |
| `data` | `true` | data vs MC path |
| `beamonly` | `false` | restrict the flash window to `[beam_mintime, beam_maxtime]` |
| `flash_minPE` | `500` | minimum total flash PE |
| `flash_mintime` / `flash_maxtime` | `∓1.5 ms` | accepted flash-time window |
| `beam_mintime` / `beam_maxtime` | `∓5 µs` | beam window (used when `beamonly`) |
| `cluster_t0` | `-1e12` | fixed cluster T0 override (disabled by default) |
| `inpath` / `outpath` | `pointtrees/%d` | tensor-set in/out paths |
| `require_containment` | `false` | discard bundles whose cluster leaves the TPC box at the flash T0 (SBND-on; `project_ql_tpc_containment`) |
