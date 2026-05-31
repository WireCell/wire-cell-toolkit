# The Charge–Light (Q–L) Matching Algorithm

This document describes the charge–light (Q–L) matching algorithm implemented in
`toolkit/match/`. It explains how the likelihood / cost function is built, how candidate
charge–flash pairs are filtered, the special treatments applied, and the choice of key
parameters. It also analyses **why the matching currently runs per‑APA (per‑TPC)** and what,
algorithmically, stands in the way of matching the two TPCs of an SBND drift volume *together*
— the motivation for a planned future improvement. A final section catalogues the hard‑coded
constants that should later be made configurable.

All file/line references are to the source as of this writing; they are pointers, not exact
guarantees if the code moves.

---

## 1. Overview and data flow

The matcher is a single Wire‑Cell Toolkit component, `WireCell::Match::QLMatching`
(`match/src/QLMatching.cxx`, `match/inc/WireCellMatch/QLMatching.h`), implementing
`ITensorSetFilter` + `IConfigurable`.

The SBND graph runs **one chain per APA** (`cfg/pgrapher/experiment/sbnd/qlmatching.jsonnet`):

```
TensorFileSource  ->  FlashTensorToOpticalPCs  ->  QLMatching
 (opflash_apa<n>)      (matrix -> flash PCs)        (per-APA Q-L match)
```

- `TensorFileSource` reads the per‑APA optical archive `opflash_apa<n>.tar.gz`.
- `Aux::FlashTensorToOpticalPCs` expands that flash matrix into the canonical
  `flash` / `light` / `flashlight` point clouds on the **live root node** of the cluster tree
  (`aux/src/FlashTensorToOpticalPCs.cxx`) — the same schema MicroBooNE's
  `root/UbooneClusterSource` uses, so all `clus` facade tooling interoperates.
- `QLMatching::operator()` (`QLMatching.cxx:140`) runs **once per APA event**: it reads the
  live/dead pctree at `inpath/{live,dead}`, extracts that APA's flashes and clusters, runs the
  two‑round LASSO match, and writes back, per matched cluster, a `cluster_t0` and a `flash`
  scalar (the matched flash row index) so `Clus::Facade::Cluster::get_flash()` reflects the
  match (`QLMatching.cxx:637-660`).

The component never sees more than one APA's data per call. This single fact drives most of
the cross‑TPC limitations discussed in §8.

---

## 2. Flash formation (optical reconstruction)

**Flashes are not built inside `match/`.** They arrive *pre‑reconstructed* as a 2‑D matrix
tensor of shape `[nflash, 1 + nchan]`:

- column 0 = flash time (ns);
- columns 1..nchan = per‑channel photo‑electron (PE) counts (`nchan = 312` for SBND).

`FlashTensorToOpticalPCs` expands this matrix (`FlashTensorToOpticalPCs.cxx:73-144`) into:

- a `flash` PC — one row per flash: `time`, `tmin`, `tmax`, `value` (= total PE),
  `ident` (= row index), `type`;
- a `light` PC — sparse, one row per fired channel: `ident` (= channel), `time`, `value`
  (= PE), `error`;
- a `flashlight` PC — the join table linking flash rows to light rows.

### Flashes are per‑APA / per‑TPC, not global

Each APA has its **own** optical archive (`opflash_apa0.tar.gz`, `opflash_apa1.tar.gz`, …) and
its own node chain. Flash row indices (`ident`) are therefore unique only *within an APA*. To
keep the matched‑flash association unambiguous after the per‑APA trees are later merged into an
all‑APA tree, `QLMatching` synthesises a **globally‑unique flash id**:

```
gid = anode_ident * kFlashGidStride + per_APA_flash_index   // kFlashGidStride = 1,000,000
```

(`QLMatching.cxx:26-30`, applied at `:650-651` and `:676`.) The per‑APA measured light is also
persisted into a merge‑safe `opflash` PC keyed by `gid` (`QLMatching.cxx:669-688`) so the
downstream all‑APA MABC Bee dump can show every flash, matched or not.

**Consequence for the planned improvement:** because flashes are physically grouped by APA in
the input data, a true two‑TPC joint match would need both APAs' flashes available in one
matching call (see §8).

---

## 3. The cost / likelihood function

### 3.1 Predicted light from charge

For each candidate (flash, cluster) pair, the algorithm predicts the PE that the cluster's
charge would produce on every optical detector (`QLMatching.cxx:268-311`). For each 3‑D point
of each blob in the cluster:

```
q                 = blob->charge() / blob->npoints()           // charge per point
x                 = point.x() + flash_x_offset                 // drift-corrected x
pred[idet] += q * QtoL * dir_vis[idet] * VUVeff[idet]
            + q * QtoL * ref_vis[idet] * VISeff[idet]
```

- `flash_x_offset = sign_offset * flash_time * drift_speed` is the drift correction that places
  the charge at the x implied by the flash time (`QLMatching.cxx:245`);
- `dir_vis` / `ref_vis` are the **direct (VUV)** and **cathode‑reflected (VIS)** visibilities
  from the SBND `SemiAnalyticalModel` (`match/src/SemiAnalyticalModel.cxx`), evaluated with the
  point converted mm→cm (`QLMatching.cxx:297`);
- `VUVeff` / `VISeff` are per‑OpDet quantum/collection efficiency arrays
  (`QLMatching.h:77-82`, 312 entries each);
- `QtoL` is the charge‑to‑light scale factor.

The `SemiAnalyticalModel` is a port of `larsim`'s `PhotonPropagation/SemiAnalyticalModel`:
Gaisser‑Hillas + solid‑angle parameterisation for direct light, a position‑parameterised
visibility for cathode‑reflected light. It is loaded from a JSON table
(`semimodel_file`, default `sbnd/photodet/semi-analytical-sbnd.json`).

### 3.2 χ² between predicted and measured PE

The per‑bundle goodness of fit is computed in `TimingTPCBundle::examine_bundle()`
(`TimingTPCBundle.cxx:151-195`):

```
chi2 = Σ_j  (pred_pe[j] - pe[j])^2 / (pe[j] + pe_err[j]^2)      // sum over masked-in OpDets
```

with the PE error rule (`Opflash.cxx:24`):

```
pe_err = (PE < 1) ? 0.3 : 0.3 * PE                              // 30% rule, floor 0.3
```

`ndf` counts the masked‑in OpDets that have either measured **or** predicted PE ≥ 1
(`TimingTPCBundle.cxx:187-188`). A bundle is flagged `high_consistent` if `ks_dis < 0.06`,
`ndf ≥ 3`, and `chi2 < ndf * nvalidopdets` (`:193`).

### 3.3 Kolmogorov–Smirnov shape distance

In addition to χ² (an absolute‑scale comparison), a **shape** comparison is computed as the
maximum gap between the normalised predicted and measured PE cumulative distributions
(`TimingTPCBundle.cxx:11-30`, `:173-179`). This `ks_dis ∈ [0,1]` is used as a penalty in the
round‑2 LASSO weights (§4) and in the quality gates (§5).

### 3.4 It is not a plain χ² minimisation

The actual matching does **not** minimise χ² pair‑by‑pair. It solves a global, regularised,
non‑negative least‑squares problem (§4) that *simultaneously* enforces two constraints:

- a **light constraint** — predicted PE ≈ measured PE on every OpDet of every flash;
- a **charge constraint** — each cluster's summed bundle strength ≈ 1, i.e. each cluster
  should be explained by (essentially) one flash. This is encoded by the `MF`/`PF` blocks with
  coefficient `1/delta_charge` (`QLMatching.cxx:486-489`, `:569-571`).

χ² and KS feed the fit only through the **per‑bundle weights** and the **pre/post‑fit gates**.

---

## 4. The fit: LASSO (yes), in two rounds

**Lasso (L1‑regularised regression) is used.** The solver is `Ress::solve(...)` with
`params.model = Ress::lasso` (`util/src/Ress.cxx:8-20`, `util/inc/WireCellUtil/Ress.h:29-42`),
backed by `WireCell::LassoModel` — non‑negative coordinate descent minimising

```
min_s  || y - X s ||^2  +  λ || s ||_1            (s ≥ 0)
```

where `s` is the vector of per‑bundle "strengths" (one entry per surviving (flash, cluster)
candidate, plus, in round 1, one per‑flash background term). Solver knobs:
`max_iter = 100000`, `tolerance = 1e-3`, `non_negative = true` (`Ress.h:38-40`).

The system is built in the normal‑equation form `X = PᵀP + PFᵀPF`, `y = PᵀM + PFᵀMF`
(`QLMatching.cxx:491-494`, `:574-577`), where:

- `P` / `M` are the **light** block: `nopdet * nflash` rows, each row a measured/predicted PE
  divided by `pe_err = sqrt(PE + PE_err²)`;
- `PF` / `MF` are the **charge** block: one row per cluster, coefficient `1/delta_charge`,
  enforcing Σ(strengths of that cluster) ≈ 1.

### Round 1 — with per‑flash background columns (`QLMatching.cxx:435-515`)

The design matrix has `nbundle + nflash` columns: one per candidate bundle, **plus one
"background" column per flash** that can absorb unexplained light (coefficient `1/delta_light`
on the flash‑light constraint, `:481`). Per‑bundle weights:

```
w_bundle = (|pred_tot - meas_tot| > 0.3*meas_tot) ? |pred_tot - meas_tot|/meas_tot : 0.3
w_background = 0.5
```

(`QLMatching.cxx:477-485`). After solving, bundles with `strength ≤ strength_cutoff` are
pruned (unless `beamonly`) (`:504-510`).

### Round 2 — bundle columns only, with a KS shape penalty (`QLMatching.cxx:517-626`)

The matrix is rebuilt without background columns over the surviving bundles. The per‑bundle
weight now adds a shape term:

```
w = base + delta_shape * nopdet * ks_dis / λ          // base as in round 1
```

(`QLMatching.cxx:562-565`), penalising bundles whose PE *shape* disagrees with the flash. After
solving and pruning by `strength_cutoff`, the algorithm keeps, **per cluster, only the flash
with the largest strength** (`:600-610`), then organises the surviving bundles (§5) and writes
the matches.

### Post‑match organisation (`organize_bundles`, `QLMatching.cxx:760-827`)

For each flash, the strongest bundle is chosen, and other clusters on the same flash may be
*merged into* it via `examine_bundle(candidate)` / `add_bundle` (`TimingTPCBundle.cxx:67-149`),
which accepts a merge only if it improves KS or χ² and satisfies `ks_dis < 0.2`,
`chi2/ndf < 20`. A final out‑of‑beam QA pass is then applied (§5).

---

## 5. How Q–L pairs are filtered

Filtering happens at several stages, progressively narrowing the candidate set:

**Flash gates** (`QLMatching.cxx:196-202`):
- keep flash only if `flash_mintime ≤ time ≤ flash_maxtime`;
- keep flash only if `total_PE ≥ flash_minPE`.

SBND has **312 OpDets: 120 PMTs** (`type=1`, 60/TPC) **+ 192 X‑Arapucas** (`type=0`,
96/TPC). The matching uses **PMTs only**.

- a per‑channel type mask, derived from the injected `OpDets` table: channel on iff
  its `type` is in `active_opdet_types` (config, default `[1]` = dome PMTs only) — this
  selects the 120 PMTs and drops all 192 Arapucas;
- a static channel mask `ch_mask` disables known‑bad channels (config). For SBND this is
  19 channels, **all PMTs**, leaving **101 of 120 PMTs active**;
- a **per‑TPC split by OpDet position** — an OpDet belongs to TPC0/TPC1 if `center.x`
  is on the low/high side of `cathode_x` (reproduces the old even/odd index split for
  the SBND PMTs);
- a per‑flash MC saturation mask — if `total_PE > 5000` and a channel reads exactly 0 PE in MC
  (`data == false`), that channel is masked (`:249-253`).

**Per‑point spatial gates** (within bundle construction, `QLMatching.cxx:282-293`):
- drift‑corrected `x` must lie within the TPC's drift bounds `[lo_x_bound, hi_x_bound]`;
- `|y| ≤ 2000`, `0 ≤ z ≤ 5000` (mm);
- if more than **25 %** of a cluster's points drift outside the bounds, the whole bundle is
  flagged `potential_bad_match` and dropped (`:293`, `:316-319`).

**Bundle pre‑selection** (`QLMatching.cxx:321-330`):
- drop if total predicted PE `< 10`;
- drop if `ks_dis == 1` (uncomputable — all‑zero prediction);
- drop if `chi2/ndf > 1e4`.

**LASSO pruning:** `strength ≤ strength_cutoff` (0.05) after each round (`:507`, `:591`).

**Out‑of‑beam QA** (`organize_bundles`, `QLMatching.cxx:807-826`): for a matched flash *outside*
the beam window (`time < beam_mintime` or `> beam_maxtime`), reject if
`ks_dis > 0.2`, or `chi2/ndf > 20`, or total‑PE mismatch `> 50 %`.

---

## 6. Special treatments

- **Deterministic ordering.** Flashes, clusters and bundles are iterated in id/index‑sorted
  order (`QLMatching.cxx:363-410`, the `flash_iter_order` / `cluster_iter_order` helpers and the
  inner sort by cluster index) so the LASSO design matrix is assembled identically run‑to‑run,
  independent of allocator/pointer ordering. (This matters because other parts of the stack are
  known to be sensitive to non‑determinism.)
- **`beamonly` mode.** When set, the strength‑cutoff pruning after each round is bypassed and
  the flash time window is collapsed to the beam window (`:59-60`, `:507`, `:591`).
- **Per‑TPC drift sign and bounds.** `sign_offset`, `lo_x_bound`, `hi_x_bound` are keyed off
  `m_anode->ident()` (TPC 0 vs 1) (`:232-235`).
- **Per‑flash background light.** Round 1 only, to absorb unmodelled light before the cleaner
  round‑2 fit.
- **Cluster merging onto a flash.** `add_bundle` lets several clusters share one flash when the
  combined prediction is a better fit (`TimingTPCBundle.cxx:115-149`).
- **Unmatched clusters** are represented by a bundle with a null flash and strength 0
  (`QLMatching.cxx:621-624`).

---

## 7. Key parameters

### 7.1 Configurable (with header defaults — note the SBND jsonnet overrides several)

| Parameter | Header default | SBND jsonnet value | Meaning | Ref |
|---|---|---|---|---|
| `nchan` | 312 | 312 | number of OpDet channels | `QLMatching.h:51` |
| `pmts` | true | (default) | use SBND PMT set | `:55` |
| `data` | true | `reality == 'data'` | data vs MC (MC masks saturated PMTs) | `:56` |
| `beamonly` | false | false | bypass strength pruning, beam window only | `:58` |
| `flash_minPE` | 500 | **50** | min total PE to keep a flash | `:59` |
| `flash_mintime` | −1500e3 ns | (default) | flash time window low | `:60` |
| `flash_maxtime` | +1500e3 ns | (default) | flash time window high | `:61` |
| `beam_mintime` | −5e3 ns | (default) | beam window low (QA) | `:62` |
| `beam_maxtime` | +5e3 ns | (default) | beam window high (QA) | `:63` |
| `QtoL` | 0.5 | **1.0** | charge→light scale | `:64` |
| `drift_speed` | 1.563e-3 mm/ns | `params.lar.drift_speed` | x drift correction | `:69` |
| `strength_cutoff` | 0.05 | (default) | LASSO strength keep‑threshold | `:73` |
| `ch_mask` | (none) | 19 channels | statically disabled OpDets | `:57`, jsonnet `:26` |
| `VUVEfficiency` / `VISEfficiency` | 312‑entry arrays | (default) | per‑OpDet efficiencies | `:77-82` |
| `semimodel_file` | `sbnd/photodet/semi-analytical-sbnd.json` | (passed in) | optical model JSON | `:93` |

`drift_speed` and `strength_cutoff` were already pulled out of inline hard‑coding into config
(see the header comments at `:65-73`); the rest below have not been.

### 7.2 Hard‑coded inside the algorithm

> **Update (Pass 1):** most of the constants below — and the bare-unit literals — have since
> been given proper Wire-Cell units and pulled up to `QLMatching` jsonnet config keys (defaults
> bit-identical to these values). See `match/docs/improve_progress.md` for the per-constant
> status and the config key names. The table is kept as the historical map of where each
> literal lived.

| Constant | Value | Where | Meaning |
|---|---|---|---|
| `lambda` (λ) | 0.1 | `QLMatching.cxx:419` | LASSO L1 strength |
| `delta_charge` | 0.01 | `:420` | charge‑constraint coefficient (1/δ) |
| `delta_light` | 0.025 | `:421` | flash‑background coefficient (1/δ) |
| `delta_shape` | 0.01 | `:422` | round‑2 KS penalty scale |
| `lo_x_bound`/`hi_x_bound` | ±2000 / 0 mm | `:234-235` | per‑TPC drift bounds |
| Y bound | ±2000 mm | `:287` | y acceptance |
| Z bound | 0–5000 mm | `:287` | z acceptance |
| close‑to‑PMT | 1950 mm | `:291` | flag if `|x| > 1950` |
| drift‑out fraction | 0.25 | `:293` | kill bundle if >25 % points drift out |
| pred‑PE floor | 10 | `:321` | drop bundle below this total pred PE |
| χ²/ndf preselect | 1e4 | `:327` | drop extreme misfits |
| PE‑mismatch knee | 0.3 (30 %) | `:477`, `:562` | weight knee |
| background weight | 0.5 | `:485` | round‑1 background column weight |
| MC saturation total PE | 5000 | `:251` | saturated‑PMT mask trigger |
| PE‑err rule | 0.3 | `Opflash.cxx:24` | 30 % error, floor 0.3 |
| KS preselect | 1.0 | `:323` | drop if KS uncomputable |
| high‑consistent KS/ndf | 0.06 / 3 | `TimingTPCBundle.cxx:193` | quality flag |
| merge KS/χ² gates | 0.2 / 20 | `TimingTPCBundle.cxx:109` | accept cluster merge |
| out‑of‑beam KS/χ²/PE | 0.2 / 20 / 0.5 | `QLMatching.cxx:812-817` | out‑of‑beam QA |
| gid stride | 1,000,000 | `:30` | global flash id stride |

---

## 8. Why matching is per‑APA, and limits to a joint two‑TPC match

The matcher is invoked **once per APA / per TPC**, and several design choices bake that in.
For the planned improvement (matching the two TPCs of an SBND drift region *together*), these
are the structural obstacles, roughly in order of how deep they cut.

1. **Single‑anode component.** `QLMatching` is configured with exactly one `IAnodePlane`
   (`anode`) and one `IDetectorVolumes` (`QLMatching.cxx:37-38`). The drift sign and x bounds
   come from `m_anode->ident()` (`:232-235`). There is no anode *vector* and no loop over TPCs;
   running on two TPCs today means two independent component instances that never share state.

2. **Charge is confined to one drift volume.** Each cluster's points are required to fall in
   that TPC's `[lo_x_bound, hi_x_bound]`; points outside are dropped, and a bundle losing >25 %
   of its points is killed (`:286-293`). A particle crossing the central cathode deposits charge
   in *both* drift volumes; such an event cannot be represented as one cluster here. Any joint
   match presupposes either clusters that may legitimately span `x = 0`, or an explicit pairing
   of a TPC0 cluster with its TPC1 counterpart.

3. **The OpDet space encodes TPC identity.** The per‑TPC split (now by `center.x` vs
   `cathode_x`) makes "which PMTs are visible" a function of which TPC is being processed,
   and the mask is applied **per flash**, not per bundle. A joint fit would need a *unified*
   OpDet/PE space (drop the per‑TPC split, use the full PMT set) — or per‑bundle masks, which
   would mean restructuring the LASSO design matrix so different bundles can illuminate
   different PMT subsets.

4. **Flashes are physically per‑APA.** Separate `opflash_apa<n>.tar.gz` archives, per‑APA row
   indexing, and the `gid = anode*1e6 + idx` workaround (`:26-30`) all exist precisely because
   each APA's optical reconstruction is isolated. A joint match needs both TPCs' flashes (and a
   decision about whether the two TPCs' PMTs see the *same* physical flash) inside one call.

5. **A single drift sign per call.** `sign_offset` is `−1` for TPC0 and `+1` for TPC1
   (`:233`), and the x correction `x += sign_offset * time * drift_speed` (`:245`) assumes all
   charge drifts the same way. Charge on the two sides of the cathode drifts in opposite
   directions; a joint correction must apply the correct sign per point/sub‑cluster.

**Implication.** A genuine two‑TPC joint match is not a small tweak to `QLMatching`. It needs,
upstream, clusters that can be associated across the cathode (or a dual‑anode matching
component fed both APAs), and, inside the matcher, a unified PMT/PE space, a per‑side drift
correction, and a design matrix in which one cluster's charge can predict light on both PMT
sets simultaneously. The physical pay‑off is correctly handling cathode‑crossing tracks, whose
light is shared between the two PMT sets but whose charge is currently split and matched twice.

---

## 9. Hard‑coded numbers to make configurable (for later code work)

> **Update (Pass 1 — DONE):** the tuning-relevant numbers here are now Wire-Cell-unit literals
> exposed as `QLMatching` jsonnet config keys (cross-class ones forwarded via `Match::PEErr` and
> `Match::BundleQualityParams`), with defaults bit-identical to the values below. The remaining
> items (the §D6 `std::abs(x)` logic bug, sentinels, and the SBND layout tables / even-odd mask)
> are intentionally deferred. See `match/docs/improve_progress.md` for the authoritative status.

These should later be exposed through jsonnet so they can be tuned without rebuilding. Per the
project convention, any new toggle must **default to the current behaviour** so existing SBND
production configs stay bit‑identical.

- LASSO weights: `lambda = 0.1`, `delta_charge = 0.01`, `delta_light = 0.025`,
  `delta_shape = 0.01` (`QLMatching.cxx:419-422`).
- Geometry acceptance: x bounds (currently literal ±2000/0), `|y| ≤ 2000`, `0 ≤ z ≤ 5000`,
  close‑to‑PMT `1950`, drift‑out fraction `0.25` (`:234-235`, `:287-293`). These should come
  from `IDetectorVolumes` rather than literals.
- Pre‑selection gates: pred‑PE floor `10`, `chi2/ndf > 1e4` (`:321`, `:327`).
- Weight knee `0.3`, round‑1 background weight `0.5` (`:477-485`, `:562`).
- MC saturation trigger `5000` and the PE‑err `0.3` rule (`:251`, `Opflash.cxx:24`).
- Quality thresholds: high‑consistent `0.06`/`3`, merge gates `0.2`/`20`, out‑of‑beam
  `0.2`/`20`/`0.5` (`TimingTPCBundle.cxx:109,193`, `QLMatching.cxx:812-817`).
- The even/odd OpDet split (`:238-240`) — this is the single biggest obstacle to a joint
  two‑TPC fit and should become a configurable per‑TPC OpDet mapping.

---

## 10. Comparison with the MicroBooNE prototype (`prototype_algorithm.md`)

The current SBND matcher is a **port of the MicroBooNE Wire-Cell prototype**
(`prototype_base/2dtoy/src/ToyMatching.cxx`, written up in
[`prototype_algorithm.md`](prototype_algorithm.md)). The core skeleton survived the port —
predict light per (flash, cluster), score with a *shape* (KS) and a *magnitude* (χ²) metric,
prune, then a **two-round non-negative L1 (Lasso) global solve** with a "track used at most
once" charge constraint and a `ks·(chi2/ndf)^0.8` tie-break. The Lasso constants are even
identical (`λ = 0.1`, `delta_track/charge = 0.01`, boundary seed `0.5`).

But the port **dropped most of the prototype's physics-aware handling of the hard cases** — and
this is exactly the part the planned improvements need to restore. Concretely:

| Prototype mechanism | Where (prototype) | Status in current SBND port |
|---|---|---|
| End-trimming walk (clip cluster ends that poke past anode/cathode) | `ToyMatching.cxx:176-264` | **Ported, for flags only.** `compute_endpoint_flags()` (`QLMatching.cxx:889-935`) does the inward walk to set the boundary flags, but the PE-prediction path still drops out-of-bounds points and >25% out kills the bundle (`QLMatching.cxx:418`) — the trim does not yet rescue charge. |
| `flag_spec_end` (demote truncated-end bundles when a clean match exists) | `ToyMatching.cxx:207-251`, `790-804` | **Filled, inert.** Now set by `compute_endpoint_flags()` (`QLMatching.cxx:904,928,942`); still never read (the prototype `:790-804` demotion is not ported). See `prototype_algorithm.md` §5.4. |
| `flag_at_x_boundary` down-weight to 0.2 in Lasso + suspend auto-reject | `ToyMatching.cxx:1437,1603`; `FlashTPCBundle.cxx:577-600` | **Filled, inert.** Now set by `compute_endpoint_flags()` (`QLMatching.cxx:947,951`) and propagated on merge (`TimingTPCBundle.cxx:144`); never used in a gate or weight. |
| `flag_close_to_PMT` → χ² error inflation + loosest cut ladders | `FlashTPCBundle.cxx:480-485`; `ToyMatching.cxx:518,597,702` | **Filled, inert.** The old buggy block is replaced; now set correctly by `compute_endpoint_flags()` (`QLMatching.cxx:946`) and propagated on merge (`TimingTPCBundle.cxx:143`); still never read downstream. |
| Fired-PMT test (`nfired/ntot`) fallback when not high-consistent | `FlashTPCBundle.cxx:550-600` | **Absent.** `Opflash` tracks `fired_channels` (`Opflash.cxx:26`) but the matcher never uses the fraction. |
| Cosmic-discriminator masking (`cos_pe_low/mid`) | `FlashTPCBundle.cxx:461-464` | **Absent** (different trigger/electronics; would need an SBND equivalent). |
| Long staged consistency/competition ladders (7 passes) | `ToyMatching.cxx:470-1358` | **Collapsed.** The port keeps a single pre-selection + the two Lasso rounds + a light `organize_bundles` merge. |
| Run/timestamp light-yield scaling, single-PMT veto | `ToyMatching.cxx:341-354` | **Absent** (MicroBooNE-specific). |

The net effect: the SBND matcher is **less tolerant precisely where the physics is hardest** —
at the anode, near the PMTs, and at the readout-window edges — because the prototype's leniency
machinery is present only as inert flags. The improvement program below is largely about
*re-instating that machinery in a two-TPC-aware form*, plus exploiting the second TPC, which the
single-TPC prototype never had.

---

## 11. Improvement opportunities

Organized around the four considerations raised for the planned work, with the dropped-prototype
features folded in. Per the project convention, anything that changes matching numerics should
be **jsonnet-togglable and default to current behaviour** so existing SBND production stays
bit-identical until validated.

### 11.1 Readout-window truncation (charge/light that predate the window)

**Problem.** `x` is reconstructed from drift time via `flash_x_offset = sign·t·v_drift`. A flash
early or late in the (~±1.5 ms) window implies charge that partly drifted in *before* readout
began (or after it ended): the cluster is **clipped**, so its predicted light is systematically
*low*. The current code treats this as a *defect* — it drops the out-of-window points and, if
>25% are lost, **kills the bundle** (`QLMatching.cxx:286-293`). That throws away exactly the
genuine, physically-truncated matches the prototype tries to rescue.

**Improvements (restore prototype behaviour, two-TPC aware):**
- Re-introduce the **end-trimming walk** (`ToyMatching.cxx:176-264`): when a cluster end pokes
  just past a drift edge, snap it to the edge and set `flag_spec_end` rather than discarding
  points wholesale.
- **Detect window-edge truncation explicitly** (cluster's leading/trailing slice within N ticks
  of the readout window edge) and, when present, **suspend the under-prediction penalties**:
  truncated charge legitimately under-predicts light, so the magnitude term (χ², total-PE
  mismatch gate at `:816`, round-1/2 weight knees) should be relaxed or one-sided (penalize
  *over*-prediction, forgive *under*-prediction).
- Make the flash time window, and its relationship to the **TPC readout/drift window**,
  explicit and configurable (today `flash_mintime/maxtime` are decoupled from the drift window).

### 11.2 Track at the anode / active-volume boundary (missing charge)

**Problem.** At the anode (`|x| ≈ 2000 mm`) a track may extend past the instrumented volume
(into the field cage / un-wired region), so charge is missing and the predicted light is again
too low — and this happens **right next to the PMTs**, where the semi-analytical model is least
accurate. The prototype handles this with `flag_close_to_PMT` (error inflation + loosest cuts)
and `flag_at_x_boundary` (Lasso down-weight + no auto-reject). In the port these flags are
inert (§10).

**Improvements:**
- **Actually set and use `flag_at_x_boundary`** (cluster touching either drift edge) and
  **down-weight those bundles in the Lasso** (prototype uses `0.2`), so a truncated track is not
  out-competed by a fully-contained one.
- **Re-instate the χ² error-inflation** near the PMTs (the prototype's `pe-pred>350 && pe>1.3·pred`
  rule) so a few saturated/over-bright near-PMT channels don't dominate χ².
- **Fix the dead/buggy `close_to_PMT` logic** at `QLMatching.cxx:289` (the `std::abs(x)` truthy
  test) and wire the flag into the gates/weights, or remove it if not used.
- Distinguish the two boundaries physically: in SBND **`x = 0` is the cathode** (shared between
  the two TPCs, max drift, dim/reflected light) and **`|x| = 2000` is the anode** (near PMTs,
  min drift). They deserve *different* treatment — the current single `flag` conflates them.

### 11.3 Cross-TPC connectivity near the cathode (a new constraint)

**Problem / opportunity.** SBND has two drift volumes sharing a **central cathode at `x = 0`**.
A cosmic crossing the cathode deposits charge in *both* TPCs, and the two charge pieces are
spatially continuous at `x ≈ 0` (same `y, z`, consistent drift-corrected time). The single-TPC
prototype had no analogue of this; the current per-APA matcher matches the two halves
*independently* and cannot use the link.

**Improvements:**
- **Pre-link cross-cathode cluster pairs**: a TPC0 cluster and a TPC1 cluster whose cathode-side
  ends are proximate in `(y, z)` and consistent in drift-corrected `t0` are almost certainly one
  physical track. (Reuse the clustering/proximity facilities in `clus/` for the `(y,z)`
  adjacency test.)
- **Tie linked pairs to a single flash** in the global solve: one shared "track used once"
  charge-constraint row spanning both TPCs, so the combined charge predicts light on **both**
  PMT sets simultaneously. This gives a far stronger lever arm than either half alone and
  resolves ambiguities per-TPC matching cannot — it is the concrete physics payoff that
  motivates the two-TPC redesign of §8.
- Even short of full cross-TPC clusters, **connectivity proximity at the cathode can act as a
  prior/penalty** linking the two per-APA matches to a common `t0`.

### 11.4 Are the two TPCs' light signals separated or bound? (and how to bind them)

**Current state — separated.** Light is delivered **per APA/TPC**: separate `opflash_apa0.tar.gz`
/ `opflash_apa1.tar.gz` archives (`qlmatching.jsonnet:30-37`), and inside the matcher an
**even/odd PMT split** restricts each TPC to half the OpDets (`QLMatching.cxx:238-240`). So the
two TPCs' flashes are two *independent* objects with **no link** — they are neither naturally
bound nor currently bindable.

**Physics.** A scintillation flash is *one* event in time. In SBND each PMT mostly sees its own
side's light (the cathode sits between the two arrays), so the two halves are "naturally
separated *in PMT space*" but "naturally coincident *in time*". The right way to bind them is
therefore **temporal coincidence + cathode connectivity**, not merging the PMT patterns.

**Open question to confirm upstream:** does the SBND opflash producer write a near-cathode flash
**once** (to one APA) or **duplicated** into both archives, and does a per-APA flash carry PE on
all 312 channels or only its even/odd half? This determines whether binding means *merging two
flash objects* or *adding a coincidence constraint between them*. This is upstream of `match/`
(larsoft/opreco) and should be checked before coding.

**Improvements (pick per the answer above):**
- If flashes are per-TPC: introduce a **time-coincidence association** that pairs an APA0 flash
  with an APA1 flash within a small `Δt`, forming a *global flash hypothesis* that the joint fit
  can score across the **full** PMT set (drop the even/odd mask for linked, cathode-crossing
  clusters).
- Make the **even/odd OpDet split a configurable per-TPC OpDet mapping** rather than a hard-coded
  parity (`:238-240`) — the single biggest structural blocker to a joint fit (§8).
- Feed both APAs into one matching call (dual-anode component, or an upstream fan-in) so the
  design matrix can carry both PMT sets and the shared cross-cathode charge constraint.

### 11.5 Lower-hanging, detector-agnostic cleanups

- Promote the §9 hard-coded constants to config (geometry bounds should come from
  `IDetectorVolumes`, not literal `±2000/0/5000/1950`).
- The pre-selection `chi2/ndf > 1e4` and pred-PE `< 10` gates, and the `0.3` weight knee, are
  unmotivated round numbers — tune on truth-matched MC and expose them.
- Add **per-bundle quality outputs** (KS, χ²/ndf, total-PE ratio, boundary flags) to the pctree
  so matching performance can be monitored and the cuts re-tuned without a rebuild.

---

## 12. Machine-learning approaches (research directions)

The hand-tuned cut ladders, the semi-analytical light model, and the greedy→Lasso assignment are
all replaceable or augmentable with learned components. These are exploratory; all need
**truth-matched MC** for training (the optical sim already provides per-flash truth), and all
must respect the toolkit's **run-to-run determinism** requirement (fixed seeds, deterministic
inference, no nondeterministic GPU reductions) and bounded inference cost. Roughly in increasing
order of ambition:

1. **Learned pair classifier / ranker (lowest risk).** Replace the `(ks, chi2/ndf, …)` cut
   ladders with a gradient-boosted-tree or small MLP that scores each candidate (flash, cluster)
   pair from engineered features — KS, χ²/ndf, total-PE ratio, per-PMT residuals, `Δt`,
   boundary/window flags, drift fraction, cross-cathode-connectivity score. Trained on MC
   truth, it subsumes the heuristic passes and is trivially deterministic at inference. The Lasso
   global solve can stay, fed by learned per-pair scores instead of hand cuts.

2. **Neural light-prediction surrogate (replaces/augments the semi-analytical model).** Train a
   coordinate network (e.g. a SIREN/MLP "PhotonLib surrogate") to predict per-OpDet visibility
   from `(x, y, z)` — and, importantly, a *calibrated uncertainty* — from full optical-sim or
   data. This directly attacks the near-PMT and near-cathode regions where the analytical model
   is weakest (§11.2), and a differentiable model enables gradient-based `t0`/assignment fits.

3. **Graph / bipartite neural network for assignment.** Represent clusters and flashes as the two
   sides of a bipartite graph; edges carry the pair features above plus the cross-TPC link.
   A GNN (or an edge-classifier + differentiable matching layer) predicts match probabilities
   with full global context and naturally handles many-to-one and the two-TPC topology — a clean
   home for the §11.3 cathode-connectivity constraint as an edge/message feature.

4. **Differentiable optimal-transport / Sinkhorn matching layer.** Replace (or initialize) the
   Lasso assignment with a differentiable soft-permutation that encodes "one flash per cluster"
   as a transport constraint and can be trained end-to-end against truth, with the cross-cathode
   tie as a coupling between rows. Keeps the convex-assignment spirit of the current solve while
   making it learnable.

5. **End-to-end / amortized inference.** A network that ingests a cluster's points and the
   candidate flashes and directly regresses `t0` (or the matched-flash posterior), trained
   against truth — amortizing the per-event fit. Useful as a fast prior to seed the global solve.

6. **Probabilistic (normalizing-flow / likelihood) matching.** Learn the full `p(PE | charge,
   geometry)` so the match score is a proper likelihood with uncertainties, giving principled
   handling of truncated charge (down-weight via the learned variance instead of ad-hoc rules).

**Pragmatic staging.** (1) and (2) are the highest value-to-risk: a learned pair score and a
learned light model with uncertainties would address the §11.1–11.2 boundary/window problems
*and* improve the existing Lasso pipeline without discarding it. (3)/(4) become attractive once
the two-TPC joint formulation (§8, §11.3) is in place, since their natural strength is global,
topology-aware assignment.

---

## 13. Performance: memory and runtime (measured)

Measured on the standalone SBND chain (`sbnd_xin/run_ql_evt.sh mc <idx>`,
`wct-clus-matching-perevt.jsonnet`; `DL=6.2, DT=9.8, lifetime=6, driftSpeed=1.563`) on a
64-core / 251 GB host, single-threaded per component. Numbers come from `QLMatching`'s built-in
`ExecMon` (`QLMatching.cxx:150,186`) and the debug-log timestamps; "matching time" brackets
`got live pctree` → `done with matching` (i.e. the component's own work, excluding upstream
clustering/imaging and downstream MABC/BEE).

### 13.1 What was measured (3 MC events, both APAs)

| Event | APA | flashes | preselected bundles | Σ cluster-points (logged) | matching time |
|---|---|---|---|---|---|
| EVT2  | 1 | 14 | 53  | —       | 0.59 s |
| EVT2  | 0 | 19 | 139 | 102,898 | 3.26 s |
| EVT9  | 1 | 13 | 41  | —       | 0.17 s |
| EVT9  | 0 | 15 | 230 | 6,557   | 0.19 s |
| EVT11 | 1 | 25 | 324 | —       | 2.18 s |
| EVT11 | 0 | 27 | 490 | 71,422  | 2.37 s |

Per-event matching total (both APAs): **EVT9 ≈ 0.36 s, EVT2 ≈ 3.8 s, EVT11 ≈ 4.6 s**. For
context, the whole per-event job wall (config load + clustering + matching + MABC + BEE zip) was
**1.3 s (EVT9)** and **6.6 s (EVT11)** — so on dense events the **matcher dominates the job**.

Memory: process RSS at the matcher's `ExecMon` checkpoints was **150–192 MB**, peaking at
**~215 MB** across the whole job; `QLMatching`'s own RSS *increment* is **≈ 0** (`Memory: MEM …
increment: res=0K`). The pctree load (`got live pctree` TICK) is **11–30 ms**.

### 13.2 The runtime is dominated by the per-point light prediction, not the fit

The matching time does **not** track the number of bundles — EVT9-APA0 has 230 bundles yet runs
in 0.19 s, while EVT2-APA0 has 139 bundles yet takes 3.26 s. It tracks the **total number of
cluster 3-D points** pushed through the predicted-light loop, at a near-constant
**≈ 30 µs per (flash × cluster-point)**:

```
EVT2-APA0:  3257 ms / 102,898 pts = 31.7 µs/pt
EVT9-APA0:   192 ms /   6,557 pts = 29.3 µs/pt
EVT11-APA0: 2372 ms /  71,422 pts = 33.2 µs/pt
```

So, to first order:

```
matching_time  ≈  30 µs  ×  Σ_(flash × cluster) N_points
              ≈  C · N_flash · N_points_total · (cost of the SemiAnalyticalModel call)
```

The cost is the **`SemiAnalyticalModel` visibility evaluation** (`QLMatching.cxx:299-301`):
for every cluster point, under every flash hypothesis, it computes the direct (VUV) and
reflected (VIS) visibility for all 312 OpDets. By contrast the **two LASSO solves are negligible
— 1–4 ms** even at 490 bundles (`solving (round 1)` → `done with matching` spans ≤ 4 ms in every
case). Loading the pctree is ~11–30 ms. The matcher is firmly **CPU-bound in the light model**,
not memory-bound and not solver-bound.

### 13.3 Two avoidable costs in the current loop

1. **The masked OpDets are still computed.** `detectedDirectVisibilities` /
   `detectedReflectedVisibilities` fill all 312 channels (`:299-301`), and the even/odd TPC mask
   is applied only at accumulation (`:304`). A per-TPC run therefore spends ~half its visibility
   math on channels it discards.
2. **The same point's visibility is recomputed for every flash.** The loop nests
   flash → cluster → point → model-call (`:243,260,280,299`), so a point is evaluated `N_flash`
   times. Only the drift offset `flash_x_offset` differs per flash; the model is re-run from
   scratch each time. This is the `N_flash` factor in the scaling above.

### 13.4 Optimization opportunities (no behaviour change intended)

- **Replace per-point analytic evaluation with a voxel/library lookup.** The MicroBooNE
  prototype used a voxelized photon library (`convert_xyz_voxel_id`, O(1) per point); the SBND
  port swapped that for an on-the-fly semi-analytical computation, which is the root cost.
  A precomputed per-voxel visibility table (built once from the same model) would turn the
  ~30 µs/point analytic call into a table lookup.
- **Compute only the masked-in OpDets** (push the even/odd mask into the model call) — roughly a
  2× win per TPC for free.
- **Downsample / aggregate charge points** before prediction (e.g. predict per blob centroid or
  per voxel rather than per point) — the predicted PE is a linear sum, so coarse-graining trades
  a controlled accuracy loss for a large speed-up on dense clusters (the EVT2/EVT11 case).
- **Hoist the flash loop**: compute each point's per-OpDet visibility once and reuse it across
  flashes (only the drift x-shift changes), instead of recomputing `N_flash` times.
- **Parallelize the flash (or APA) loop** — the bundle predictions are independent.

### 13.5 Implication for the two-TPC improvement

A joint two-TPC fit (§8, §11.3) does *not* materially change the dominant cost: the light
prediction is already per-point-per-flash, and a joint fit predicts the same points against a
similar set of flashes — the LASSO that grows (more columns) is the cheap part. The main new cost
is dropping the even/odd mask so a cathode-crossing cluster predicts onto the **full** 312-OpDet
set, which roughly doubles the inner OpDet loop for the linked bundles. Given §13.3/§13.4, the
right sequencing is to **fix the light-prediction cost first** (voxel lookup + mask-aware
evaluation), which both speeds up today's matcher and removes the only real performance objection
to the joint formulation.

---

## 14. Source map

| Component | Files |
|---|---|
| Main matcher | `match/src/QLMatching.cxx`, `match/inc/WireCellMatch/QLMatching.h` |
| Bundle (χ², KS, merge) | `match/src/TimingTPCBundle.cxx`, `.../TimingTPCBundle.h` |
| Flash adapter (PE, PE_err) | `match/src/Opflash.cxx`, `.../Opflash.h` |
| Optical model | `match/src/SemiAnalyticalModel.cxx`, `.../SemiAnalyticalModel.h` |
| LASSO solver | `util/src/Ress.cxx`, `util/inc/WireCellUtil/Ress.h` |
| Flash tensor → PCs | `aux/src/FlashTensorToOpticalPCs.cxx` |
| SBND graph config | `cfg/pgrapher/experiment/sbnd/qlmatching.jsonnet` |
