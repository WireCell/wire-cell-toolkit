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
plus how the prototype reconstructs the per-channel PE error upstream of matching (§7.1), how the
LASSO penalty weights (`delta_shape`/`bkg_weight`/`pe_mismatch_*`) enter the fit (§6.1), the bundle
add/merge thresholds (§4.1), the MC saturation mask and out-of-beam QA cuts (§12), the
light-prediction model that builds the predicted PE the χ²/LASSO consume (charge→PE, §13), the
two-block LASSO matrix architecture (§14), and a complete reference for the toolkit's jsonnet config
parameters (§11).

### 1.1 End-to-end matching pipeline (the five stages)

Both codebases share the **same five-stage skeleton**. The mental model "bundles are formed, then
fit, then merged" is correct — with the key facts that (a) the *only* thing happening **between the
two LASSO rounds is a weak-bundle prune** (no merging), and (b) **merging is strictly post-fit**.

| Stage | Toolkit `QLMatching.cxx` | Prototype `ToyMatching.cxx` |
|---|---|---|
| **1. Pre-filter** (geometry + light + consistency cull) | TPC-containment cull (`require_containment`, SBND-on, §15); KS==1 degenerate guard at formation (`:549-552`); **over-prediction light cull** (`reject_overpred`, SBND-on, §15); then for each cluster with a consistent bundle, drop its other non-consistent bundles (`:642-653`) | TPC-containment gate (`flag_good_bundle`); for each flash with a consistent bundle, drop remaining bundles no consistent bundle can absorb (`examine_bundle` merge-test, `:1307-1358`); fired-fraction reject (`:547-602`) |
| **2. R1 fit** | LASSO with per-flash background/light DOF; post-fit prune strength ≤ `0.05` (`:670-751`) | LASSO with per-flash "flash-alone" DOF; post-fit prune `beta < 0.05` (`:1377-1548`) |
| *(between rounds)* | **prune only — no merging** (`:743-750`) | **prune only — no merging** (`:1516-1547`) |
| **3. R2 fit** | background DOF dropped; KS-shape term added to weights; prune + keep best flash per cluster (`:753-849`) | flash DOF dropped; `f1` 0.06→0.05; keep best beta per cluster (`:1552-1663`) |
| **4. Examine / merge bundles** | `organize_bundles` merges split clusters lit by the same flash (`examine_bundle`/`add_bundle` gates, `:1198-1262`) | post-LASSO `examine_bundle_rank`/`add_bundle` merge of non-best bundles into the per-flash best (`:1839-1841`) |
| **5. Post-matching** | **remove only**: out-of-beam QA drops bundles failing `ks>0.2 \|\| chi2/ndf>20 \|\| \|meas−pred\|/meas>0.5`; dropped cluster is unmatched, no recovery (`:1264-1284`) | **aggressively re-match** (`organize_matched_bundles`, `:1744`): orphaned clusters re-tried against *other* flashes over a 3-round cascade (`:1786-1944`) — full algorithm in §12.3 |

Two placement nuances:

- In the **toolkit**, stages 4 and 5 both live inside `organize_bundles` (merge `:1198-1262`, then QA
  `:1264-1284`) — one function, called once after R2 (`:864`).
- In the **prototype**, stages 4 and 5 are **interleaved per flash** in the post-LASSO region: pick the
  best bundle for a flash, merge what it can absorb, push the rest to the orphan re-match rounds.

The **single biggest divergence is stage 5**: the prototype recovers orphaned clusters (re-matches
them to other flashes); the toolkit only *removes* failing bundles and never retries — a known port
gap (see §12.2).

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

### 3.1 The ported χ² relaxation (implemented, `chi2_relax`)

The prototype's two per-bundle χ² relaxations (§3) are now ported behind
`BundleQualityParams::chi2_relax` (`TimingTPCBundle::examine_bundle`; **C++ default OFF =
bit-identical; SBND-on**), keeping the toolkit's `pe + perr²` Poisson base:

- **Close-to-PMT denominator inflation**: when `flag_close_to_PMT` and a channel shows a big measured
  excess (`pe − pred > chi2_pmt_excess` PE **and** `pe > chi2_pmt_ratio·pred`), the denominator is
  widened by `(pe·chi2_pmt_inflate)²`. `chi2_pmt_excess` is absolute PE and was **re-validated for
  SBND**: the prototype's `350` is kept, but only after confirming the SBND close-to-PMT measured
  excess sits at a *median ~2300 PE* (p90 ~7000), so 350 cleanly admits genuine over-response while
  excluding sub-350 noise. `chi2_pmt_ratio = 1.3`, `chi2_pmt_inflate = 0.5` (prototype values).
- **One-inefficient-PMT subtraction**: the single worst-χ² channel, if `pe == 0 && pred > 0`, is
  dropped (`chi² −= max_χ² − 1`).

**Live but benign.** On the 20 hand-scan events the inflation *fires* (the excess condition is met)
but changes **no** `flag_high_consistent` ladder or end-to-end matching decision: the ladder is
KS-led and its B4 branch (§4.1b) already loosens the χ² ceiling for close-to-PMT/boundary/truncated
bundles, so reducing those bundles' χ² further does not move the gate. It is kept as the faithful
prototype port (it does change χ² *values*, so it is not a no-op for any future χ²-keyed consumer).

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
and to reject bundles outright — also boundary-aware. The toolkit ports the *over-prediction*
direction of this as a simpler, data-tuned **light prefilter** (`reject_overpred`, §15): a bundle is
dropped before the fit when its predicted light hugely exceeds the measured light (the prototype's
"way more light than the flash shows" case), boundary-exempt like the prototype.

### 4.1b The ported SBND ladder (implemented, `highconsist_ladder`)

The single toolkit branch is replaced — **behind `highconsist_ladder` (C++ default OFF =
single-branch, bit-identical; SBND-on)** — by a flag-aware ladder that ports the prototype's
*structure* but **re-derives every number from the 10 hand-scan data events**. The prototype's
χ²-per-ndf multipliers (36/45/55) do **not** transfer: after the SBND PE-error recipe a good
match has χ²/ndf ≈ 1–2 (not MicroBooNE's tens), so those multipliers would be absurdly loose here.
The empirics that set the design (`work/ql_recipe` dumps vs `work/ql_labels/data`): **KS is the
purity lever, χ²/ndf is not** — in every flag category true matches sit at ks ~0.09–0.12 while
background sits at ks ~0.44–0.68, but χ²/ndf is ~15 for *both*. So the ladder is KS-led, with
χ²/ndf ceilings only fencing the tail. Branches (`TimingTPCBundle.cxx`, OR; `c2n = chi2/ndf`):

| # | branch | condition | maps to user idea |
|---|---|---|---|
| B1 | clean very-good | `ndf≥3 && ks<0.06 && c2n<6` | "very good KS + low χ²" |
| B2 | general good | `ndf≥3 && ks<0.09 && c2n<4` | (clean, moderate) |
| B3 | two_boundary | `flag_two_boundary && ndf≥3 && ks<0.10 && c2n<8` | "flag_2_boundary: both reasonable" |
| B4 | x-bnd / close-PMT / window-trunc | `(flag_at_x_boundary‖flag_close_to_PMT‖flag_window_truncated) && ndf≥5 && ks<0.08 && c2n<60` | "these flags: KS reasonable, χ² can be higher (missing charge)" |

This is the **TIGHT** operating point, chosen because `flag_high_consistent` gates the pre-LASSO
cull (`cull_inconsistent`, §1.1) and so **must be pure**; lower efficiency is acceptable. B4
deliberately relaxes the χ² ceiling, *not* the KS — the simpler analogue of the prototype's
per-channel `+(pe·0.5)²` denominator inflation (§3), which is **not** ported. The eight ceilings
+ `hc_miss_min_ndf` are jsonnet-exposed (`BundleQualityParams::hc_*`) for retuning; `flag_two_boundary`
is now computed in the matching path whenever the ladder is on (previously calib-dump-only).

**Measured impact** (the 10 hand-scan data events; `work/ql_ladder` vs the single-branch baseline):
flag recall **18%→44%**, purity **82%→88%**, true-matches culled **2→1**. End-to-end, true-match
agreement in the matched output rises **89→93 / 100** (+4). MC (validation only, *not* retuned):
purity 94%, recall 53% (cleaner sim → looser), end-to-end **87→92 / 113** (+5). No regressions.

`flag_spec_end` (prototype demotion gate, `ToyMatching.cxx:790-804`) is **not ported**: it is
never set in SBND (0/99 selected, 0/762 background), so it is inert here.

### 4.1 Bundle add/merge thresholds (`§F` knobs)

`flag_high_consistent` (above) is one of several thresholds carried in the `BundleQualityParams`
struct (`qp`, built at `QLMatching.cxx:435-437`) and forwarded to **every** `TimingTPCBundle`. The
other three govern **bundle merging** — combining a candidate cluster into an existing flash↔cluster
bundle when a flash plausibly matches more than one cluster.

**When/where merging happens.** Merging is **not** part of forming bundles, and **not** part of the
χ² fit — it is a *post-fit consolidation*. The flow is three stages: (1) **formation** — each
flash×cluster pairing becomes a bundle and `examine_bundle()` (no-arg) computes its own KS + χ²/ndf
(the §3 χ², `QLMatching.cxx:548`); (2) **the LASSO fit** (rounds 1–2, `:670-831`) assigns strengths
and prunes, after which a single flash may still retain *several* surviving bundles; (3) **merge** in
`organize_bundles` (called at `:864`, body `:1198-1262`) — for each flash, the highest-strength
bundle becomes the anchor and the flash's other survivors are folded into it
(`examine_bundle(candidate)` then `add_bundle`, `:1240-1241`). **Why:** the upstream clustering can
split one physical object (a long track, or a multi-prong interaction) into several clusters all lit
by the same flash; the LASSO may keep several, and merging recombines those co-illuminated clusters
into one consolidated `flash ↔ (multi-cluster)` match. The cuts below ensure only a cluster whose
added light keeps the combined photon-pattern match good is absorbed.

- **`bundle_ks_merge_max` (0.2), `bundle_chi2ndf_merge_max` (20)** — the *merge-acceptance* gate
  (`TimingTPCBundle::examine_bundle(candidate)`, `TimingTPCBundle.cxx:109-112`). A candidate is
  merged only if the *combined* prediction both improves the fit **and** stays within quality:

  ```cpp
  if ((candidate_ks_dis < ks_dis || candidate_chi2 < chi2) &&        // merging helps, and
      candidate_ks_dis < ks_merge_max && candidate_chi2/ndf < chi2ndf_merge_max)  // stays good
      return true;   // accept the merge
  ```

  So a merge that would push the combined KS above 0.2 or χ²/ndf above 20 is rejected even if it
  nominally lowers one of them. These are the merge-time analog of the consistency cuts in §4.

- **`bundle_addmerge_exponent` (0.8)** — the *main-cluster tie-break* in `add_bundle`
  (`TimingTPCBundle.cxx:118-121`). Once a merge is accepted, the two clusters compete for the
  "main" role by the combined metric `ks_dis · (chi2/ndf)^0.8`; the smaller wins and becomes the
  bundle's `main_cluster` (the other is demoted to `other_clusters`). The exponent sets how much
  χ²/ndf weighs against KS in that ranking (0.8 < 1 → KS dominates slightly).

- **`bundle_pe_ndf_knee` (1.0)** — channels where *both* `pe` and `pred` fall below 1.0 PE do not
  count toward `ndf` (`TimingTPCBundle.cxx:104,192`), so empty-on-both-sides OpDets neither inflate
  χ²/ndf nor the consistency denominators.

**Prototype equivalent.** The prototype merges at the same post-LASSO stage, via
`examine_bundle_rank(other)` + `add_bundle(other)` (`ToyMatching.cxx:1839-1841`, `1888-1891`;
`FlashTPCBundle.cxx:153,350`). Two differences in form: (a) its merge-acceptance gate uses
**relative-improvement** bounds (`temp_ks < ks+0.06 && temp_ks < ks·1.2 && temp_chi2 < chi2+5·ndf &&
temp_chi2 < chi2·1.21`, with a looser OR-branch; `FlashTPCBundle.cxx:106-112`) rather than the
toolkit's improves-OR + **absolute ceilings** `0.2`/`20` — so `bundle_ks_merge_max`/
`bundle_chi2ndf_merge_max` are a toolkit *simplification*, not direct prototype constants; (b) the
`0.8` tie-break exponent **is** a direct port, but the prototype normalizes by predicted light —
`ks·(chi2/ndf)^0.8 / pred_light` (`FlashTPCBundle.cxx:356`) — a `/pred_light` factor the toolkit
drops. (Separately, the prototype also uses `examine_bundle(other)` *pre*-LASSO at `:1329` as a
keep/drop consistency cull — a different use, not a merge.) See §9 for the re-tune view.

## 5. Flag usage

| Flag | Toolkit: where SET | Prototype: how USED | Toolkit status |
|---|---|---|---|
| `flag_close_to_PMT` | `QLMatching.cxx:1128` (anode-end, within `m_anode_ext2`) | χ² denom inflation `+ (pe·0.5)²` for big excess (`FlashTPCBundle.cxx:480-485`); relaxed consistency `χ²<55·ndf` (`:518`); LASSO weight 1.0→0.2 via paired `flag_at_x_boundary` | **active**: LASSO down-weight (`lasso_flag_weight`, §6.2) + χ² denom inflation (`chi2_relax`, §3.1) + ladder B4 (§4.1b) |
| `flag_at_x_boundary` | `QLMatching.cxx:1129` (anode), `:1155` (cathode, CPA 3-D fiducial) | relaxed consistency `χ²<9·ndf`/`3·ndf` (`FlashTPCBundle.cxx:510-513`); LASSO weight 1.0→0.2 (`ToyMatching.cxx:1437-1441`, rd1; `1603-1607`, rd2) | **active**: LASSO down-weight (`lasso_flag_weight`, §6.2) + ladder B4 (§4.1b) |
| `flag_spec_end` | `QLMatching.cxx:1080,1104,1124` | informational; tail-alignment bookkeeping | **inert** (not set in SBND; §4.1b) |
| `flag_window_truncated` | `QLMatching.cxx:1032-1037` (slice within `m_window_edge_ticks` of readout window edge) | *no direct prototype analog* — see note below | **active**: LASSO down-weight (`lasso_flag_weight`, §6.2) + ladder B4 (§4.1b) + cross-TPC scenario 2 (§16) |
| `flag_potential_bad_match` | toolkit sets via `ks==1` drop **and the over-prediction light prefilter** (`reject_overpred`, §15) | fired-fraction reject, boundary-aware (`FlashTPCBundle.cxx:577-599`) | **active** (over-prediction direction; §15) |

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
| Background column weight (rd1) | `bkg_weight = 0.5` | `1.0` (fixed, `ToyMatching.cxx:1449`) |
| Per-bundle L1 weight content | PE-total mismatch (`pe_mismatch_floor = 0.3`, `pe_mismatch_knee = 0.3`) + KS shape (rd2) | **boundary flag only** — `0.2`/`1.0`, *no* PE-mismatch or shape term (`ToyMatching.cxx:1437-1441`, `1603-1607`) |
| Boundary down-weight | *(not applied — flags inert)* | `0.2` for `flag_at_x_boundary` (this **is** the prototype's per-bundle weight) |
| Strength keep-threshold | `m_strength_cutoff = 0.05` | `0.05` |

The λ value and `delta_*` constants are **dimensionless ratios** (the prototype comment says "the
coefficient is all around 1"), so they are good candidates to **port as-is**; the boundary
down-weight `0.2` is likewise dimensionless. See §9.

### 6.1 How `delta_shape`, `bkg_weight`, `pe_mismatch_*` enter the fit

These four are **not** part of the per-bundle χ² of §3 — that χ² (`examine_bundle`,
`TimingTPCBundle.cxx:194`) is built before the fit and only gates `flag_high_consistent`. The four
parameters here build the **LASSO `weights` vector**, i.e. the *L1 penalty* of the regression, a
separate term from the data χ². The solver (`util/src/LassoModel.cxx:11`) minimizes

```
        ½‖y − Xβ‖²        +        N · λ · Σ_k  weights(k) · |β_k|
        └── data term ──┘          └──────── L1 penalty ────────┘
         (the χ²: y,X hold          (weights(k) is a per-column
          pe/pred over pe_err)       penalty multiplier)
```

`β_k` is the strength of column *k* (one per candidate bundle, plus one background column per flash
in round 1). `y`/`X` carry the `pe/pe_err` and `pred/pe_err` data — that is where the χ² of §3/§7
lives. The four parameters only scale `weights(k)`: a **larger** `weights(k)` shrinks `β_k` harder
toward zero (`_soft_thresholding(…, λ·weights(k))`, `LassoModel.cxx:125`), so the bundle is more
likely to fall below the `strength_cutoff = 0.05` keep-threshold and be dropped. They do not change
any `(pred − pe)²/σ²` value; they change *which columns survive* the joint fit.

**`pe_mismatch_knee` / `pe_mismatch_floor` — per-bundle penalty** (`QLMatching.cxx:713-715` rd1,
`:798-800` rd2). For each bundle column the base penalty is the *relative total-PE mismatch*:

```cpp
weights(k) = (|pred_tot − meas_tot| > pe_mismatch_knee · meas_tot)
               ? |pred_tot − meas_tot| / meas_tot   // mismatch beyond the knee: penalize ∝ mismatch
               : pe_mismatch_floor;                 // within the knee: flat floor
```

So a bundle whose predicted total light is within `knee = 0.3` (30%) of the measured flash PE pays
only the floor `0.3`; beyond 30% it pays a penalty equal to its fractional mismatch (e.g. a 2×
overshoot → weight 1.0). `knee` is the dead-band half-width; `floor` is the minimum penalty any
bundle pays so even a perfect-PE bundle is still L1-regularized.

**`bkg_weight` — background-column penalty** (round 1 only, `QLMatching.cxx:721`). Round 1 adds one
"background light" column per flash (`PF = 1/delta_light`) that can absorb flash PE no bundle
explains. `weights(nbundle+k) = bkg_weight = 0.5` is the L1 price of dumping light into that
background instead of into a real bundle: lower `bkg_weight` makes the background cheaper (more
flash PE explained away as background, bundles compete harder); higher makes real bundles
preferred. Round 2 has no background column, so `bkg_weight` is round-1-only.

**`delta_shape` — KS shape penalty** (round 2 only, `QLMatching.cxx:801`). Round 2 augments each
bundle's base penalty with a shape term:

```cpp
weights(k) = base + delta_shape · nopdet · ks_dis / lambda;   // base = the pe_mismatch term above
```

A bundle whose predicted *pattern* across OpDets matches the measured one poorly (large `ks_dis`)
pays extra L1 penalty, scaled by the channel count `nopdet` and divided by `lambda` (so the
`N·λ·weights` product is λ-independent in this term). This is **toolkit-specific** — the prototype
has no shape term in its LASSO and instead enforces shape upstream through the `ks_dis` cuts of
`flag_high_consistent` (§4). `delta_shape = 0.01` is small by design: it nudges the round-2 ranking
toward better-shaped bundles without overriding the PE-mismatch term.

**Contrast with the prototype.** None of these three are prototype features. The prototype's
per-bundle L1 weight is the **boundary flag only** — `0.2` if `flag_at_x_boundary` else `1.0`
(`ToyMatching.cxx:1437-1441`, `1603-1607`), with no PE-mismatch and no shape dependence — and its
background column has a fixed weight of `1.0` (`:1449`). So the two codes put entirely different
content into the per-column penalty: prototype = drift geometry (the boundary down-weight, currently
inert in the toolkit, §5); toolkit = PE-total mismatch + KS shape. From the likelihood view above
these are just different choices of *per-candidate prior plausibility* — geometry-based vs
agreement-based — fed into the same weighted-LASSO machinery.

### 6.2 The ported flag down-weight (implemented, `lasso_flag_weight`)

The prototype's boundary down-weight is now ported behind `lasso_flag_weight` (**C++ default OFF =
factor 1.0, bit-identical; SBND-on**). `QLMatching::lasso_flag_factor(bundle)` returns
`lasso_boundary_weight (0.2)` when the bundle is `flag_at_x_boundary || flag_close_to_PMT ||
flag_window_truncated` (the ladder-B4 group, generalizing the prototype's `flag_at_x_boundary`-only
0.2), else 1.0. It **multiplies** the existing per-column base (`pe_mismatch` floor/knee in round 1;
`base + KS-shape` in round 2) in `fit_round1`/`fit_round2`, so the toolkit's agreement-based content
is kept and only *scaled down* for boundary bundles. A down-weighted bundle pays less L1 penalty,
is shrunk less, and so survives the `strength_cutoff` instead of being killed — the point being that
a boundary/truncated bundle's *measured* light is an underestimate (§5, §13.4), so a raw PE-mismatch
penalty would wrongly reject a real match.

**Measured impact** (10 data + 10 MC hand-scans, true-match agreement): **DATA +2 (92→94)**, MC
net-neutral (one event loses two marginal `ks > 0.25` boundary matches, offset by gains elsewhere).
Purity-first, net-positive, no per-event data regression — kept.

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
| `delta_shape` | `0.01` | KS shape-penalty scale (round 2, toolkit-specific; §6.1) |
| `bkg_weight` | `0.5` | background-column L1 penalty (round 1; §6.1) |
| `pe_mismatch_knee` | `0.3` | PE-mismatch penalty knee (dead-band, fraction of measured PE; §6.1) |
| `pe_mismatch_floor` | `0.3` | PE-mismatch penalty floor (minimum L1 weight; §6.1) |
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
| `bundle_ks_merge_max` | `0.2` | KS ceiling for the merge-acceptance gate (§4.1) |
| `bundle_chi2ndf_merge_max` | `20` | χ²/ndf ceiling for the merge-acceptance gate (§4.1) |
| `bundle_addmerge_exponent` | `0.8` | exponent in the main-cluster tie-break `ks·(chi2/ndf)^p` (§4.1) |
| `highconsist_ks_max` | `0.06` | `flag_high_consistent` KS cut (§4) |
| `highconsist_min_ndf` | `3` | `flag_high_consistent` min ndf (§4) |
| `bundle_pe_ndf_knee` | `1.0` | PE below which a channel does not count toward ndf (§3, §4.1) |

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

**§E — out-of-beam QA cuts** (post-fit cleanup of out-of-beam matches; see §12.2)

| jsonnet key | default | meaning |
|---|---|---|
| `outbeam_ks_max` | `0.2` | KS ceiling for an out-of-beam match |
| `outbeam_chi2ndf_max` | `20` | χ²/ndf ceiling for an out-of-beam match |
| `outbeam_pe_frac` | `0.5` | out-of-beam total-PE mismatch fraction |

**§D — pre-selection / bad-match gate** (see §12.1)

| jsonnet key | default | meaning |
|---|---|---|
| `mc_saturation_pe` | `5000` | total-flash-PE trigger for the MC saturated-PMT mask (MC only) |

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

## 12. Selection gates outside the χ²: MC saturation mask & out-of-beam QA

Two toolkit gates act on the data *around* the fit rather than inside the χ². They are listed in §11
(§D, §E) but their *when/where/why* is here.

### 12.1 `mc_saturation_pe` — per-flash channel mask (pre-fit, MC only)

**Where** (`QLMatching.cxx:455-461`): just before the bundles for a flash are built, the per-TPC
`opdet_mask` is copied to a per-flash `flash_opdet_mask`, and for that copy:

```cpp
if (flash->get_total_PE() > m_mc_saturation_pe && pe_det == 0 && m_data == false)
    flash_opdet_mask[idet] = 0;   // drop this channel for this flash only
```

**What/why** (per the code's own comment at `QLMatching.cxx:455`, "also catches simulated saturated
PMTs in MC"): in **MC** the optical simulation represents a **saturated** PMT by emitting `pe = 0`
for that channel. In a **bright** flash (total PE above `mc_saturation_pe = 5000`) a healthy PMT
should see plenty of light, so a `pe = 0` channel is taken to be one of those simulated saturated
PMTs. Left in, it would contribute a large bogus `(pred − 0)²` to every bundle's χ² for that flash.
Masking it (`flag = 0`) excludes it from the χ² (§3) and the LASSO data term for that flash only; it
also drops out of the KS — but only because the channel is `(measured 0, predicted 0)`, not because
the mask is applied to the KS loop (which it is not — see §12.4). The gate is doubly guarded — **MC only** (`m_data == false`) and **bright only** — so
dim flashes and all data are untouched (data handles saturated/dead channels by other means, e.g.
`ch_mask`). `5000` is an absolute-PE threshold and is a **RE-TUNE for SBND** candidate (it scales
with light yield, like the §9 entries). The behavior is verified from this toolkit code and its
comment; no corresponding bright-flash zero-mask was found in the MicroBooNE prototype matcher
(whose dead-PMT handling is a separate data-run channel veto), so this is treated as toolkit/SBND
logic rather than a port.

### 12.2 Out-of-beam QA cuts (post-fit cleanup)

**Where** (`QLMatching.cxx:1264-1284`): the **final** pass, after the LASSO rounds, the
best-per-cluster selection, and bundle merging have produced `results_bundles`.

**What/why**: a matched bundle whose **flash time is outside the beam window**
(`< beam_mintime || > beam_maxtime`) is held to a **tighter** standard than in-beam matches, because
an out-of-beam (cosmic) match has no beam-timing prior backing it up and is more likely spurious. It
is dropped if **any** of:

```cpp
bundle->get_ks_dis()  > outbeam_ks_max            // 0.2  — shape too far off
bundle->get_chi2()/bundle->get_ndf() > outbeam_chi2ndf_max   // 20  — χ²/ndf too high
|flash_total_PE − pred_light| > outbeam_pe_frac · flash_total_PE   // 0.5 — >50% total-PE mismatch
```

In-beam bundles bypass this entirely (the whole block is inside the out-of-beam time test). So these
three are **acceptance thresholds on the already-computed** `ks_dis`, `chi2/ndf`, and total-PE
mismatch — they prune final cosmic matches, they do not feed back into the χ² or the fit. The KS and
χ²/ndf ceilings mirror the `§F` merge ceilings (0.2 / 20) but are applied at a different stage and
only to out-of-beam flashes.

**What happens to the dropped cluster — it is unmatched, with no recovery.** This is the *terminal*
step (`organize_bundles` is the last processing stage; afterwards `:866-871` only builds the output
map). When a bundle is erased here, its cluster — plus any clusters merged into it (§4.1) — has **no
entry left** in `results_bundles`: no flash, no T0. The toolkit does **not** retry the orphaned
cluster against any other (e.g. in-beam) flash. Note a dropped cluster is actually worse off than a
*never-matched* one: a cluster that simply found no match keeps a `flash=nullptr` strength-0
placeholder (`:857-861`), whereas a QA-dropped cluster is removed outright.

**Prototype equivalent — none directly; the prototype is the opposite philosophy.** The toolkit's
cull is keyed on a **flash-time window**; the prototype keys on **flash type** (1 = cosmic,
2 = beam) and its nearest analogous code runs the *reverse* polarity:

- `examine_beam_bundle()` (`FlashTPCBundle.cxx:383-443`, called `ToyMatching.cxx:1277`) is a
  **promotion**, not a rejection: for a *beam* flash left with no consistent bundle, it promotes the
  best bundle to consistent if it passes `(ks < 0.1 || ks₁ < 0.05) && (χ² < 12·ndf || …)`. It rescues
  beam matches rather than culling cosmic ones.
- For clusters orphaned by the post-fit merge, the prototype runs a **multi-round re-matching**
  (`ToyMatching.cxx:1870-1944`): a cluster that lost its bundle is retried against *other* flashes
  (2nd round), then still-other flashes (3rd round), and is only dropped (`set_flash(0)`, bundle
  deleted) once all flashes are exhausted.

So §12.2 is **toolkit-specific**, and the larger design gap is the missing **orphan re-matching**:
the prototype gives a cluster several chances at different flashes before abandoning it, whereas the
toolkit drops it on the first failed out-of-beam QA. Whether SBND wants the prototype's recovery
rounds is a port decision (§10-adjacent; not yet on that checklist). The full algorithm is §12.3.

**Resolution (2026-06-04 gap review).** The high-value, light-separable slice of the orphan
re-matching IS now ported as `empty_rescue` (§4.4a): each LASSO-emptied flash adopts its best
light-quality candidate (`ks·(χ²/ndf)^0.8`) from the pre-fit snapshot, one-flash-per-cluster.
The residual full 3-round loop is empirically inert for SBND — of ~5 hand-scan misses only 1 is
light-separable (recovered at zero regression); the rest are timing/drift-degenerate (no light bar
separates recover-from-steal) or cross-TPC (handled by the `xtpc` machinery). The one genuinely
**live, un-ported** piece is `examine_beam_bundle`'s beam-flash *promotion*: the toolkit ladder has
no "beam flash with no consistent bundle → promote its best" fallback, so the toolkit is strictly
*more permissive* there (it leaves more candidates in the LASSO; it cannot drop a true match).
Low SBND impact; addable behind a default-OFF / SBND-on knob if a beam-recall need ever shows up.
Verdict: **post-fit / merge / re-matching port is complete for SBND**; only this minor promotion
remains optional.

### 12.3 Prototype stage-5 algorithm — orphan re-matching (`organize_matched_bundles`)

This is the prototype's entire post-fit stage (`ToyMatching.cxx:1744`, called at `:1724`). Its input
`results_bundles` is the post-R2 state: each cluster carries **at most one** matched bundle (its
best flash, from the best-per-cluster step), and unmatched clusters carry a `flash = 0` placeholder.
The routine groups the matched bundles by flash (`:1748-1762`) and runs **three cascading rounds**.

**The merge test used throughout** is `examine_bundle_rank(candidate)`
(`FlashTPCBundle.cxx:153-232`): it builds the *combined* (anchor + candidate) predicted light, then
accepts only if the combined KS and χ² do **not degrade beyond relative bounds** — a 3-way OR
(`:221-229`), e.g. `temp_ks < ks+0.06 && temp_ks < 1.2·ks && temp_chi2 < chi2+5·ndf && temp_chi2 <
1.21·chi2`, or a looser pair, or simply `temp_ks·temp_chi2 < ks·chi2`. In words: *merging this
cluster in must not make the joint fit meaningfully worse.*

**Round 1 — merge within each flash (`:1786-1868`).** For each flash owning ≥1 bundle:
1. Pick an **anchor** `best_bundle`: the one minimizing the score
   `value = ks_dis · (chi2/ndf)^0.8` (`:1807`), with flag down-weights that deliberately *favor*
   boundary/PMT bundles as the anchor — `×0.8` if `flag_at_x_boundary`, a further `×0.8` if also
   `flag_close_to_PMT`, and `×0.75` more if the running best is poor (`ks>0.2 && chi2>60·ndf`)
   (`:1809-1818`).
2. For every other bundle of that flash, test `best_bundle->examine_bundle_rank(bundle)`. **Pass** →
   `add_bundle` folds that cluster's light into the anchor (`:1839-1841`) and the bundle leaves the
   results. **Fail** → the cluster is an **orphan**, pushed to `second_round_bundles` (`:1849`).
3. If the flash is beam (`type == 2`), `best_bundle->examine_merge_clusters()` does an extra
   beam-only cluster merge (`:1863-1866`).

**Round 2 — re-home orphans onto already-matched flashes (`:1870-1916`).** For each orphan cluster:
scan the existing result bundles (reverse order, skipping bundles that are themselves orphans). If an
already-matched flash had **this cluster as a candidate** (the `(flash, cluster)` pair exists in
`fc_bundles_map`) *and* that flash's anchor accepts it via `examine_bundle_rank`, then `add_bundle`
merges the cluster onto that flash (`:1885-1894`); the flash is marked **tried**. If no flash
absorbs it → `third_round_bundles` (`:1899`).

**Round 3 — last chance on untried flashes (`:1918-1944`).** For each still-orphan cluster: collect
every candidate bundle `(cluster, flash)` from `fc_bundles_map` whose flash was **not yet tried**
(`:1922-1927`). If exactly one exists, accept it (`:1928-1930`); if several, accept the one
minimizing `ks_dis · (chi2/ndf)^0.8` (`:1931-1942`). The orphan's original bundle then has its flash
zeroed (`set_flash(0)`, `:1943`) — so if *no* untried candidate exists, the cluster is finally left
unmatched.

**Cleanup (`:1947-1953`).** Any `fc_bundles_map` bundle not in the final results is deleted.

**Net effect:** a cluster gets **up to three chances** at a flash — its own best (R1), any
already-matched flash that wanted it (R2), then any untried flash (R3) — and is abandoned only after
all are exhausted. Note the asymmetry the toolkit would need to replicate: rounds 2–3 reach into
`fc_bundles_map` (the *full* (flash, cluster) candidate set, kept alive from before the LASSO), not
just the surviving matched bundles — the prototype keeps every candidate pairing available for
recovery. The toolkit's `organize_bundles` only merges (§4.1) and then *removes* out-of-beam failures
(§12.2); it has no analogue of any of these three rounds.

#### 12.3a Toolkit status — the organize result is vestigial; the empty-flash rescue (SBND-on)

Two findings from porting this stage (validated on 10 data + 10 MC hand-scans):

1. **`organize_bundles`' output is discarded.** `fit_round2` builds a `results_flash_bundles_map`
   from the organized `results_bundles` and then never assigns it back to `run.flash_bundles_map` —
   the matched output is just the strength-cutoff survivors. So the toolkit's per-flash best-pick
   *removals* and the §12.2 out-of-beam QA are dead w.r.t. output (the same-flash `add_bundle` *merge*
   mutation IS observable, via the shared `shared_ptr`s). See `qlmatching-code.md` §4.4.

2. **The recoverable misses are timing/drift-degenerate, not light-recoverable.** The single piece of
   the prototype's recovery that the hand-scan misses actually need is the *light-quality, flash-centric*
   pick (rounds 2–3 reach `fc_bundles_map` and rank by `ks·(chi2/ndf)^0.8`, **strength-independent**).
   The ported `rescue_empty_flashes` (§4.4a, `empty_rescue`) does exactly that for **emptied** flashes,
   enforcing one-flash-per-cluster by reassignment. But the misses where a cluster lost its flash to a
   neighbour are degenerate: the cluster fits its **wrong** flash as well as the correct one (a correct
   match can have light metric 25.9 while an empty flash out-fits it at 0.68), so **no light bar
   separates recover-vs-steal**. Only **one** of ~5 hand-scan misses is light-separable (MC evt11
   `(10,8)`: 0.13 at the correct flash vs 6.37 at the wrong one). The conservative bar
   `rescue_metric_max = 0.5` (below the lowest data-regression metric 0.68) recovers exactly that one
   at **zero regression** (DATA 95→95, MC 98→99). The rest need a drift/timing discriminator the
   standalone matcher lacks; the cross-TPC ones belong to §16's `xtpc` machinery. Default OFF ⇒
   bit-identical.

### 12.4 Masking and the KS shape metric (a toolkit asymmetry)

Per-bundle quality (§3) is two numbers: the χ² and the KS shape distance. The channel **mask**
(`opdet_mask` — built from opdet type, `ch_mask`, the opposite-TPC side, and the per-flash
`mc_saturation` bit; §12.1) is applied consistently to the χ² and the LASSO fit, and **historically
not to the KS** — an asymmetry now closed behind the `bundle_mask_ks` toggle (see **Fix** below). This
subsection records how masked detectors enter the KS in each codebase, the small piece of KS
arithmetic that makes the `mc_saturation` case harmless but the geometry/`ch_mask` case not, and the
fix. The description below is of the **OFF** (legacy) behaviour.

**Does a `(measured = 0, predicted = 0)` channel change the KS?** No — it is **identical to dropping
it**. The KS here is the maximum absolute difference of the two normalised cumulative distributions
(toolkit `calc_ks_test`, `TimingTPCBundle.cxx:14-30`; prototype `TH1::KolmogorovTest(…,"M")`). A
channel that is zero in *both* distributions (i) adds 0 to both normalisation totals, so every other
channel's normalised value is unchanged, and (ii) leaves both cumulatives flat across it
(`cum_m[k]=cum_m[k−1]`, `cum_p[k]=cum_p[k−1]`), so `|Δcum|` there merely repeats the adjacent value
and can never be a new maximum. The KS is therefore moved **only** by channels where measured and
predicted differ — in particular where `predicted = 0` but `measured > 0`.

**Toolkit — the KS skips the mask (`TimingTPCBundle.cxx:156-184`, merge variant `:68-99`).** The KS
loop fills `measured_dist[j]=pe[j]` and `predicted_dist[j]=pred_pe[j]` over **all `m_nchan`** with
**no `opdet_mask` test**; the mask is consulted only in the χ²/ndf loop (`:190`) and in the LASSO fit
(`opdet_idx_v`, `QLMatching.cxx:664,697,781`). For a masked channel the **prediction is 0** (the
light loop skips masked dets, `:530`) but the **measured PE is not zeroed** — `Opflash::init` stores
the full input vector `PE = flash.pes(nchan)` with no per-TPC/active zeroing (`Opflash.cxx:19-20`), so
`get_PE` returns real light on the opposite-TPC and `ch_mask` channels. Those channels are thus
`(measured > 0, predicted = 0)` → by the rule above they **shift the KS**, even though χ² and the fit
correctly exclude them. The existence of the opposite-TPC mask (`:447-448`) is itself evidence those
channels carry PE — otherwise it would be unnecessary. So this is a genuine **inconsistency**: the KS
does not apply the exclusion that χ² and the LASSO both apply.

- **Exempt:** the `mc_saturation` mask is conditioned on `pe_det == 0` (`:459`), so those channels are
  `(0, 0)` and — by the arithmetic above — invisible to the KS regardless. The asymmetry bites only
  the geometry/`ch_mask` channels that saw light.

**Prototype — no such asymmetry, by construction.** Single TPC, 32 PMTs, no per-channel mask vector;
the KS (`FlashTPCBundle.cxx:459-470`) and χ² (`:477-500`) both run over all 32 identically. "Masking"
is done by **zeroing the prediction** — `norm_factor[17]=0` for the dead PMT (`ToyMatching.cxx:354`)
and the type-1 cosmic veto (`FlashTPCBundle.cxx:461-464`) — and the dead PMT's measured PE is ~0 too,
so its bin is `(0, 0)` and drops out of the KS cleanly.

**Fix (implemented).** The toolkit KS now gates on `opdet_mask`, matching the χ²/LASSO paths. In both
`examine_bundle` overloads (`TimingTPCBundle.cxx:68-114`, `:156-202`) the loop that builds
`measured_dist`/`predicted_dist` zeroes a masked channel's measured PE (the prediction is already 0
there), so it becomes `(0, 0)` and — by the arithmetic above — drops out of the KS, exactly as if it
had been skipped. Following the codebase convention this is a behaviour change behind a knob,
`bundle_mask_ks` (struct `BundleQualityParams::mask_ks`), **C++ default OFF** so every other config
stays bit-identical, and set **`true` for SBND** in `cfg/pgrapher/experiment/sbnd/qlmatching.jsonnet`
(mirrors the `require_containment` precedent). Toggling it OFF recovers the old (asymmetric) KS for
A/B comparison; the magnitude of the change depends on how much real PE the masked channels actually
carry, which the toggle is there to let you measure.

## 13. Light prediction: charge → predicted PE per channel

Everything above consumes `pred_flash` (the per-channel predicted PE) — it is the `pred_pe` in the χ²
(§3), the columns of `P` in the LASSO (§6), and the shape compared in the KS (§12.4). But the two
codebases **build that vector by fundamentally different methods**, and this is the single most
important technical area the rest of the document took as given. The measured side (`pe`, `pe_err`)
is covered in §7; this section is the **predicted** side.

### 13.1 Prototype — voxelized photon library (table lookup)

Predicted light comes from a **pre-computed optical library** (GEANT/optical sim baked in); the
matcher only does lookups (`ToyMatching.cxx:305-360`):

- Each merged cell's charge is split over its sampling points: `charge = mcell->get_q()/pts.size()`
  (`:315`); each point is quantised to a voxel id (`convert_xyz_voxel_id`, ~`:20-40`, a few-cm grid).
- The library returns, per voxel, a list of `(OpChannel, visibility)` pairs
  (`photon_library->at(voxel_id)`, `:327`). Predicted PE accumulates
  `pred += charge · visibility / elifetime_ratio` (`:330`), with the library channel remapped to the
  PMT index via `map_lib_pmt`.
- 32 PMTs; one detector technology.

### 13.2 Toolkit — `SemiAnalyticalModel` (computed at runtime)

There is **no pre-baked library**; visibilities are computed analytically per `(point, OpDet)` from
detector geometry (`QLMatching.cxx:497-540`, `SemiAnalyticalModel.cxx`):

- Per blob point, `q = blob->charge()/blob->npoints()` (`:505`), positions converted **mm→cm** (`:523`).
- **Direct (VUV)** visibility = geometric solid angle × `exp(−d/vuv_absorption_length)` × a
  Gaisser-Hillas border correction (`SemiAnalyticalModel.cxx:265,413-421`). The solid angle uses a
  **dome** model for PMTs (`Omega_Dome_Model`, `:534`) and a **rectangle** model for (X)Arapuca flat
  PDs (`Rectangle_SolidAngle`, `:465`).
- **Reflected (VIS)** visibility = a two-hop model: point→cathode "hotspot", then hotspot→OpDet
  (`:310-358`).
- Predicted PE: `pred += q·QtoL·(dir_vis·VUVEff + ref_vis·VISEff)` (`QLMatching.cxx:535-536`), summed
  over the whole cluster **group** (main + associated sub-clusters, `:497-540`).

### 13.3 Three prediction differences that matter for the port

1. **Electron-lifetime correction: present in the prototype, absent in the toolkit matcher.** The
   prototype divides the predicted contribution by `elifetime_ratio` (`ToyMatching.cxx:323,330`) to
   recover the *original* ionization — charge measured at the wire is post-drift-attenuation, but
   scintillation is emitted pre-drift, so the light scales with the un-attenuated charge. A grep of
   toolkit `match/` finds **no** lifetime/attenuation term; `QLMatching` feeds `blob->charge()`
   straight in. Either the toolkit's input blob charge is already lifetime-corrected upstream
   (clustering/imaging) **or this correction is unported** — must be confirmed for SBND, since a
   missing correction biases predicted light low for far-drift charge (tens of %, drift-dependent).
2. **Light-magnitude normalization differs.** Prototype: one `scaling_light_mag` (≈0.015) × a
   per-**run** `yield_ratio` (TGraph, ~0.62–1.0, data only), with `rel_light_yield_err` folded into
   the PE error (§7). Toolkit: an absolute `QtoL` (default **0.5**) × **per-OpDet**
   `VUVEfficiency`/`VISEfficiency` arrays, and **no** run-dependent yield-drift term. These are the
   §9 "DECISION / RE-TUNE" rows — §13 is where they physically enter.
3. **Two detector technologies (SBND).** The toolkit models dome PMTs (`type 1`) and rectangular
   (X)Arapucas (`type 0`) separately, selecting which to use via `active_opdet_types` (default `{1}` =
   PMTs only, reproducing the historical SBND mask, `QLMatching.cxx:256-263`). The prototype is
   PMT-only. Turning Arapucas on adds rectangle-solid-angle channels to every prediction and hence to
   the χ²/KS — a deliberate, not automatic, extension.

### 13.4 The two TPCs are assumed *completely* optically separated (a hard cathode cut)

The SBND model treats the cathode (`x = 0`) as an **opaque optical wall**: a charge deposit in one
TPC produces **exactly zero** predicted PE on every OpDet of the *other* TPC — not a small number, an
identical zero, by an early `continue` before any solid-angle computation. Both visibility routines
open with the same line (a port of larsim's `SBNDOpticalPath_tool`):

```cpp
// SemiAnalyticalModel.cxx:233 (direct/VUV) and :355 (reflected/VIS)
// only OpDets on the same X-sign as the scintillation point are visible.
if ((scintPoint.x() < 0.) != (od.center.x() < 0.)) continue;
```

`vis` is pre-zeroed, so opposite-side OpDets keep `0`. This holds for **both** light components:
the reflected path places the cathode "hotspot" at `plane_depth` with the **same sign** as the
deposit's `x` (`:317, :352`) and re-applies the same-side cut in the per-PD sum (`:355`), so reflected
light never crosses the cathode either.

- **Geometry confirming the split** (`wire-cell-data/sbnd/photodet/semi-analytical-sbnd.json`):
  `cathode_x = 0`; **312 OpDets split 156 / 156** across it (x ≈ ±213.5 cm, behind each anode);
  192 Arapucas (`type 0`) + 120 dome PMTs (`type 1`). So the x-sign *is* the TPC label.
- **Redundant with the matcher mask.** Even without this cut, QLMatching's per-TPC mask
  (`QLMatching.cxx:445-449`) zeros the opposite-side OpDets in each TPC hypothesis. Predicted far-side
  PE is therefore doubly guaranteed to be 0.
- **This is the physics behind the §12.4 KS masking.** The prediction is strictly one-sided, but a
  reconstructed flash carries measured PE on **both** TPCs' PMTs (the opflash is not split per-TPC).
  The opposite-TPC channels are thus `(measured > 0, predicted = 0)` and *must* be dropped from
  χ²/LASSO/KS (the per-TPC mask + `bundle_mask_ks`), else a one-sided prediction is scored against a
  two-sided measurement.
- **Modeling caveat.** "No cross-cathode light" is a *hard* assumption inherited from
  `SBNDOpticalPath_tool`, not a tunable. If SBND's TPB-coated reflective cathode actually leaks/reflects
  a non-negligible fraction across, the model cannot represent it (it would under-predict the far side);
  the only available remedy is masking, not modeling. The prototype (MicroBooNE, single TPC) has no
  analogue — its photon library is a single-volume lookup with no inter-TPC question.

## 14. LASSO matrix architecture (the two-block normal-equation system)

§6/§8 cover the *columns* and the *weights*; this section records the **row** structure — the soft
constraints that §6 only names ("Σ strengths per cluster → 1"). Both codes assemble the **same**
stacked system and solve it in normal-equation form
(`QLMatching.cxx:727-730`; prototype `ToyMatching.cxx:1481-1485`):

```
   y = Pᵀ·M + PFᵀ·MF        X = Pᵀ·P + PFᵀ·PF        solve  min ½‖…‖² + N·λ·Σ wₖ|βₖ|
```

Two row blocks are stacked before forming the normal equations:

- **Block 1 — measurement rows** (`nopdet × nflash` of them). `M(i·nopdet+j) = pe/pe_err` is the
  measured light; `P` holds, in each bundle column, `pred_pe/pe_err`, plus (round 1 only) one
  **background-light column per flash** holding the raw `pe/pe_err` so unexplained flash PE has
  somewhere to go (`QLMatching.cxx:699-708,717`; prototype `ToyMatching.cxx:1415-1432`). This block
  *is* the χ²/data term.
- **Block 2 — regularization rows.** One **charge-conservation row per cluster**:
  `PF(cluster, bundle) = 1/delta_charge` with target `MF(cluster) = 1/delta_charge`
  (`QLMatching.cxx:722-724`; prototype `delta_track`, `ToyMatching.cxx:1468-1477`). Because every
  bundle anchored on a cluster writes `1/delta_charge` into that cluster's single row against a target
  of `1/delta_charge`, the fit is **pulled toward Σ(strengths on that cluster) = 1** — i.e. *each
  cluster is used about once*, which is what stops one cluster's charge from being split across many
  flashes. Round 1 adds a **per-flash light row** (`PF(flash, bkg_col) = 1/delta_light`; prototype
  `delta_flash`) regularizing the background columns.

The small `delta_*` (0.01–0.025) act as **large** penalty coefficients (`1/delta` ≈ 40–100), making
these near-hard constraints. The per-column L1 `weights` (§6.1) ride on top of this and decide which
columns survive; the `strength_cutoff = 0.05` then converts the soft solution into hard matches (§8).
This architecture is shared; only the per-column `weights` content (§6.1) and the §13 light predictor
differ between the codes.

## 15. The bundle prefilter (what reaches the χ² fit)

Before the LASSO/χ² fit (§8), each (flash × cluster) bundle passes through a prefilter in
`build_bundles` (`QLMatching.cxx`). Two of its stages are **geometry** then **light**, both ported
from the MicroBooNE prototype, both **default OFF in C++ / SBND-on in jsonnet** (convention
`feedback_toggleable_behavior_changes`). Surviving bundles are what the fit sees.

### 15.1 Geometry — `require_containment` (already present; SBND-on)

`require_containment` drops a bundle whose charge cluster is not contained in the TPC drift box once
the flash-T0 x-offset is applied (`QLMatching.cxx`, cull at `if (m_require_containment && !contained)
continue;`). `contained` is the 4-part box test in `compute_endpoint_flags`, a direct port of the
prototype `flag_good_bundle` gate (`ToyMatching.cxx:272-275`) with the same cushions (−2 cm anode,
+1.2 cm cathode — "not too far outside the box"). Validated against the 20 hand-scans: **0 of the
214 ground-truth bundles are dropped** by containment (`sbnd_xin/ql_prefilter_tune.py` step 1).

### 15.2 Light — `reject_overpred` (new; SBND-on)

A bundle is dropped when its **predicted** light is much larger than the **measured** light: a charge
cluster emitting far more light than the flash actually shows cannot be the match. It is
**one-directional** — under-prediction is allowed (a flash may be lit by several clusters; near-PMT
over-response and window truncation make *measured* an underestimate). This ports the over-prediction
direction of the prototype fired-fraction reject (`FlashTPCBundle.cxx:547-602`) as a simpler,
data-tuned cut. Two metrics, both over the **same `opdet_mask` the χ² uses** (PMT, same TPC, active):

```
R_total = Σ pred_masked / Σ meas_masked
R_max   = pred[ch*] / max(meas[ch*], 1),   ch* = argmax(pred_masked)
reject if (R_total > overpred_total_ratio  OR  R_max > overpred_maxch_ratio)
```

Boundary/truncated bundles (`close_to_PMT | window_truncated | at_x_boundary`) are **exempt**
(measured underestimates there), mirroring the prototype's boundary-aware reject and the flash drops
in `ql_pe_error.py`. A rejected bundle gets `flag_potential_bad_match` and is skipped before
`pre_bundles`, exactly like the KS==1 guard, so it never enters the fit.

**Tuning & validation** (`sbnd_xin/ql_prefilter_tune.py`, parity `ql_prefilter_parity.py`). Cuts were
set on the **10 data hand-scans** as (worst non-boundary ground-truth) × 1.5 margin and confirmed on
the **10 MC hand-scans** (the GT baseline is not 1.0 — true matches over-predict ~1.3 — and shifts
data↔MC via `data_qtol=0.86` vs `sim=1.0`, so the data→MC transfer was checked explicitly):

| | worst non-boundary GT `R_total` | worst GT `R_max` | GT removed | non-GT culled |
|---|---|---|---|---|
| data (tune) | 1.92 | 2.89 | 0 | 26% |
| MC (validate) | 1.11 | 2.81 | 0 | 32% |

Chosen SBND values (`cfg/pgrapher/experiment/sbnd/qlmatching.jsonnet`): `overpred_total_ratio = 2.9`,
`overpred_maxch_ratio = 4.3`. Post-build re-dump confirms **full C++/Python parity** (the C++ flags
exactly the bundles the Python metric predicts: 190/190 data, 355/355 MC) and **0 GT rejected**.

Config knobs (defaults inert in C++): `reject_overpred` (false), `overpred_total_ratio` (1e9),
`overpred_maxch_ratio` (1e9).

## 16. Cross-TPC cathode-crossing cull (`flag_xtpc_consistent`)

A cathode-crossing cosmic is reconstructed as two halves — one per TPC — lit by **one coincident
flash group**. If the per-TPC matcher picked the *correct* cluster in both TPCs, the two halves form
one continuous track across the cathode. That cross-TPC geometry is an independent constraint the
per-TPC χ²/KS ladder (§4) cannot see (it is single-TPC). It is used as a **pre-fit cull**: a pair of
candidate bundles confirmed as one cathode-crosser drops their competing rivals before the LASSO, so
fewer bundles enter the fit (the §1.1 stage-1 idea, extended across TPCs). It is a **separate** flag,
not overloaded onto `flag_high_consistent`.

*(History: this was originally a post-matching observation-only confirm-stamp that also wrote a
per-cluster `xtpc_consistent` output scalar. It is now a pre-fit cull that DOES change matching, and
the output scalar was removed — only the `-calib` bundle field remains. The geometry/cuts below are
unchanged except for the new scenario-2 distance ceiling.)*

**Where.** The per-APA pipeline `run_one_apa` is split at the LASSO boundary into
`run_one_apa_prefit` (bundles + the per-TPC `cull_inconsistent`) and `run_one_apa_fit` (the two LASSO
rounds + output). `operator()` runs the prefit for **all** APAs, then
`QLMatching::cull_cross_tpc(runs)` (multi-APA, `xtpc_flag`), then the fit for all APAs. The OFF path
is the historical single-loop order, bit-identical. `cull_cross_tpc` reuses the `-cathode-diag`
geometry primitives (`xtpc_pair_consistent`: full-cluster closest-point pair +
`Cluster::vhough_transform` local dirs + `conn`, with the per-TPC transverse offset `ApaRun.dy/dz`
applied), but pairs **candidate** main-cluster bundles (every coincident-flash bundle, not just
matched mains). A wrong-flash pairing carries the wrong T0 x-offset ⇒ large `d` ⇒ self-vetoing. A
confirmed pair sets `flag_xtpc_consistent` on both bundles; then per run, a cluster owning an
xtpc-consistent bundle drops its other bundles that are neither high- nor xtpc-consistent.

**The two scenarios:** for each TPC0×TPC1 candidate-main pair with coincident flash times
(`|t0−t1| < flash_group_window`), in the T0-corrected frame,

- **Scenario 1 — closest distance** (cathode end present): `flag = d < xtpc_dmax (5 cm)`.
- **Scenario 2 — three-vector collinearity** (a half is `flag_window_truncated`, cathode end missing
  so `d` is large): `conn`, `dir0`, `dir1` mutually collinear (`a01,a0c,a1c < xtpc_angle_max (20°)`)
  **AND** `d < xtpc_dmax2 (300 cm)`. The distance ceiling is load-bearing in the candidate
  population — see below.

Combined: `flag_xtpc_consistent = (d < 5 cm) OR (window_truncated AND d < 300 cm AND a01,a0c,a1c < 20°)`.

**Cuts — tuned on the 10 SBND hand-scan data events only** (MC validation; reconstruct every
coincident **candidate**-main pair, TRUE iff both halves hand-scan-selected, on the real C++ `vhough`
values via a temporary all-pairs log). Scenario 1 is clean: no FALSE pair has `d < 326 cm`. Scenario
2 is the subtle one: pairing *candidate* (not just matched) bundles admits far-apart collinear
truncated pairs the post-fit version never saw — a data false at `d=473 cm` and an MC false at
`d=406/326 cm`, all with angles `< 14°`. **Angle cannot separate them** (the MC false at `3.6°` is
*more* collinear than the data true evt1302 at `2.8°`); **distance can** — every true scenario-2 pair
sits at `d ≤ 264 cm`, every false at `d ≥ 326 cm`, hence `xtpc_dmax2 = 300 cm`. With the ceiling:
**flag purity 100% (data 16/16, MC 30/30 flagged true).** End-to-end (true-match agreement vs the
post-fit baseline): the cull removes ~46 (data) / ~139 (MC) rival bundles pre-fit and takes **MC
recall 91→97 (+6), DATA flat 92**, with no per-event data regression. (Tune the vhough angles/`d` on
the real C++ values, never a python PCA proxy — the PCA proxy mis-estimates the angle.)

Config (default OFF = bit-identical; SBND-on): `xtpc_flag`, `xtpc_dmax`, `xtpc_dmax2`,
`xtpc_angle_max`, `xtpc_hough_radius` (§11 table in `qlmatching-code.md`). Output: the
`flag_xtpc_consistent` bundle field in the `-calib` dump (no root-node scalar).

## 17. QLMatching runtime cost (timing & memory)

**Measured on our own locally-imaged clusters** (file 1 of `input-3files-lan-reco2`, 50
real BNB+cosmic data events, SBND-on, joint two-TPC matching). Active **and** dead blobs
were imaged in-toolkit from the SP frames (`wct-img-all.jsonnet` + in-graph masked
imaging) — **no foreign `*-active.npz`** (see [[feedback_own_imaging_only]]). Per-event
cost is instrumented directly in `QLMatching::operator()` via the existing `ExecMon em`,
which emits one line per event:

```
QLMatching timing: ident <N> took <ms> ms, proc RSS <MB> MB (delta <MB> MB)
```

`took` is the full `operator()` wall (prefit + cross-TPC cull + fit + output build).
`proc RSS` is the **whole matching process** resident set at that point (size→resident is
`MemUsage::current().second`); `delta` is QLMatching's own incremental RSS over its entry
baseline. The chain is single-threaded (`Pgrapher`), so no other node runs inside the
measured window — the only logging overhead inside it is QLMatching's own ~150 debug
lines/event (tens of ms vs seconds), i.e. the numbers are compute, not log I/O.

### 17.1 Timing — QLMatching is seconds/event, not "sub-second"

Per-event `took`, 49 events (stable across three independent runs):

| metric | value |
|---|--:|
| median | **5.2 s** |
| mean | 7.5 s |
| p99 / max | **31 s** (event 60933, the busiest cosmic pile-up) |
| total (49 evt) | **369 s** |

For context, the whole matching process (dead imaging + per-APA clustering + all-APA
clustering + QLMatching) ran ~658 s wall over this file, so **QLMatching alone is ~56 % of
the matching-process wall**. This corrects the earlier "QLMatching is sub-second" note
(which was measured against yuhw's sparser larsoft active.npz; see
[[project_clustering_examinebundles_hotspot]]): on our richer multi-3view active clusters
the cost scales with cluster/point count — hence the long tail (median 5 s, max 31 s).

> **Per-phase drill-down (2026-06-05, supersedes the cause speculated here).** A later
> per-phase profile of event 60933 (`qlmatching-perf-event60933.md`) found the dominant
> cost is **not** PCA/`vhough`/`get_extreme_wcps` or the two-round LASSO (the LASSO solve
> is ~1.4 ms). **91 % of the wall is `cull_cross_tpc`** — its `xtpc_pair_consistent`
> closest-approach step is a brute-force O(npts0·npts1) double loop whose per-point
> `point3d()` re-resolves a scoped k-d view (Scope hash + hashtable lookup), ~100 ns/pair
> over 2.7×10⁸ pairs. It scales linearly with that point-pair count (84–91 % at both
> median and tail).
>
> **Fixed** by hoisting the scope resolution out of the loop (bind the raw coordinate
> arrays once; the T0 shift stays a scalar add — bit-identical). Re-running **all 150
> lan-reco2 events** post-fix: `xtpc_cull` max **27 s → 0.4 s**, and QLMatching `took`
> drops from median 5.2 s / max 31 s to **median 1.1 s / max 5.6 s**. The numbers above are
> the pre-fix baseline. Post-fix, `build_bundles`/`SemiAnalyticalModel` (median 0.9 s) is
> the dominant QLMatching phase, and the per-event reco wall is dominated by imaging and
> `ClusteringExamineBundles`, not QLMatching. See `qlmatching-perf-event60933.md`.

### 17.2 Memory — QLMatching's own footprint is negligible

| metric | value |
|---|--:|
| matching-process RSS at QL end (median) | **691 MB** |
| matching-process RSS at QL end (max) | 707 MB |
| QLMatching incremental ΔRSS (median / max) | **≈ 0 MB** (−6 .. +0) |

The process resident set (~0.7 GB, dominated by the clustering point clouds) is already
fully resident when QLMatching runs; the matcher reuses that pctree in place and allocates
essentially nothing of its own (the small negative deltas are allocator give-back during
matching). The whole-process **peak** RSS (`VmHWM`, polled externally) was **985 MB** for
the matching job and **1131 MB** for the upstream active-imaging job. So matching is not a
memory bottleneck; the per-event *time* is the cost to watch.

> **Caveat (open, 2026-06-05).** The full 150-event local-imaging reprocess is **blocked by
> an intermittent clustering heap-corruption** newly exposed by the richer local active
> imaging (segfaults in `clus_all_apa` on a layout-dependent ~50–70 % of runs; clean under
> gdb/valgrind layout, so it is a memory-corruption heisenbug, not the already-fixed
> `fill_wrap_points` off-by-one). The timings above come from the runs that completed all
> 49–50 events; root-causing is in progress (valgrind memcheck). The numbers are
> representative per-event, but the end-to-end batch is not yet crash-clean.
