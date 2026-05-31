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

**OpDet masking:**
- a static channel mask `ch_mask` disables known‑bad channels (config; 19 channels for SBND);
- a **per‑TPC even/odd PMT split** — TPC0 keeps even‑indexed OpDets, TPC1 keeps odd
  (`QLMatching.cxx:238-240`);
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

3. **The OpDet space encodes TPC identity.** The even/odd PMT split (`:238-240`) makes "which
   PMTs are visible" a function of which TPC is being processed, and the mask is applied
   **per flash**, not per bundle. A joint fit would need a *unified* OpDet/PE space (drop the
   even/odd split, use the full PMT set) — or per‑bundle masks, which would mean restructuring
   the LASSO design matrix so different bundles can illuminate different PMT subsets.

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

## 10. Source map

| Component | Files |
|---|---|
| Main matcher | `match/src/QLMatching.cxx`, `match/inc/WireCellMatch/QLMatching.h` |
| Bundle (χ², KS, merge) | `match/src/TimingTPCBundle.cxx`, `.../TimingTPCBundle.h` |
| Flash adapter (PE, PE_err) | `match/src/Opflash.cxx`, `.../Opflash.h` |
| Optical model | `match/src/SemiAnalyticalModel.cxx`, `.../SemiAnalyticalModel.h` |
| LASSO solver | `util/src/Ress.cxx`, `util/inc/WireCellUtil/Ress.h` |
| Flash tensor → PCs | `aux/src/FlashTensorToOpticalPCs.cxx` |
| SBND graph config | `cfg/pgrapher/experiment/sbnd/qlmatching.jsonnet` |
