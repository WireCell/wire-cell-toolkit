# Prototype charge–light (Q–L) matching algorithm

Source-faithful writeup of the **original MicroBooNE Wire-Cell prototype**
matching code, to anchor the standalone `match/` port (see
[`porting-summary.md`](porting-summary.md)). Every claim carries a `file:line`
reference so it can be checked against source.

Primary source: `prototype_base/2dtoy/src/ToyMatching.cxx`. The metric and
"examine" helpers live in `prototype_base/data/src/FlashTPCBundle.cxx`.

> **Ignore the commented-out regions.** In `ToyMatching.cxx`, the
> `Photon_Library` constructor (lines 43–143) and the entire
> `tpc_light_match_ana` (lines 1965–3607) are disabled. The active code is
> lines 1–1962: `calculate_pred_pe`, `tpc_light_match`, and
> `organize_matched_bundles`.

---

## 1. Overview & I/O

The goal is to associate each TPC charge cluster (**Q**) with at most one
optical flash (**L**), and vice versa. The entry point is

```
FlashTPCBundle tpc_light_match(eventTime, time_offset, nrebin, pl,
                               group_clusters, flashes, run_no, ...)
```
(`ToyMatching.cxx:368`).

- **Inputs:** a voxelized photon library `pl` (per-voxel, per-PMT visibility),
  `group_clusters` (each "main" cluster plus associated additional clusters),
  and `flashes` (optical flashes, each with 32 PMT PE values + errors + a type).
- **Output:** a `FlashTPCBundleSelection` — one bundle per cluster, carrying the
  matched flash (or `flash==0` if unmatched) and a match strength.
- Only the `solv_type==1` branch is live (`ToyMatching.cxx:412`); it is the
  "new matching code".

The matcher works on **bundles**. A `FlashTPCBundle` is a candidate
(flash, main_cluster) pair plus its predicted per-PMT light, its consistency
metrics (`ks_dis`, `chi2`, `ndf`), and a set of boolean flags
(`FlashTPCBundle.cxx:5`).

The high-level flow is:

1. Build a bundle for **every** (flash, cluster) pair and predict its light
   (`calculate_pred_pe`).
2. Score each bundle with two complementary metrics (`examine_bundle`).
3. Prune the bundle set through a long sequence of staged consistency / competition passes.
4. Resolve the remaining ambiguity with a **two-round L1 (Lasso) global solve**.
5. Tidy up the surviving matches (`organize_matched_bundles`).

---

## 2. Predicting the light for a (flash, cluster) pair — `calculate_pred_pe`

`ToyMatching.cxx:145`. For a candidate bundle:

### 2.1 The drift / `offset_x` hypothesis

A flash has a time but no x; a cluster has a drift coordinate (its `x`) but no
absolute t0. Assuming *this* flash gave *this* cluster fixes the t0, which
converts to an x offset:

```
offset_x = (flash->get_time() - time_offset) * 2./nrebin * time_slice_width
```
(`ToyMatching.cxx:167`). The cluster's *true* x under this hypothesis is then
`x - offset_x`. The drift window runs from the **anode** at
`low_x_cut = 0` to the **cathode** at `high_x_cut = 256 cm`
(`ToyMatching.cxx:158-162`). Tolerance bands around the nominal edges are the
`*_ext1/ext2` constants (`ToyMatching.cxx:158-164`).

### 2.2 Predicted PE from the photon library

For each merged cell (mcell) with positive charge, the charge is spread evenly
over its sampling points (`charge = q/pts.size()`, `ToyMatching.cxx:315`). Each
point is shifted by `offset_x`, mapped to a library voxel
(`convert_xyz_voxel_id`, `ToyMatching.cxx:20`), and the library's per-PMT
visibility is accumulated into `pred_pmt_light`, with an **electron-lifetime
attenuation correction** that depends on the drift distance
(`ToyMatching.cxx:322-331`). Channel indices are remapped library↔PE via
`map_lib_pmt` (`ToyMatching.cxx:330`).

Finally the prediction is scaled by `scaling_light_mag`, and **PMT 18
(index 17) is vetoed** (`norm_factor[17]=0`) for certain run numbers /
timestamps (`ToyMatching.cxx:341-354`).

### 2.3 Keep / drop

The bundle is only kept (`flag_good_bundle = true`, `ToyMatching.cxx:363`) if
the (trimmed — see §5) cluster physically fits inside the drift window under
this flash hypothesis (`ToyMatching.cxx:272-275`). Bundles that don't fit are
deleted immediately (`ToyMatching.cxx:459-461`).

---

## 3. The two consistency metrics — `examine_bundle`

`FlashTPCBundle::examine_bundle` (`FlashTPCBundle.cxx:446`) builds two 32-bin
histograms — measured PE (`h1`) and predicted PE (`h2`) — and computes two
complementary numbers. **The whole filtering philosophy rests on requiring both
a good shape *and* a good magnitude.**

- **`ks_dis`** — `h1->KolmogorovTest(h2,"M")` (`FlashTPCBundle.cxx:469`). The
  `"M"` option returns the maximum-distance KS statistic. This compares the
  *shape* of the light pattern across PMTs and is essentially
  normalization-independent.
- **`chi2` / `ndf`** — `sum (pred-meas)²/err²` over the 32 PMTs; `ndf` counts
  PMTs where pred or meas is non-zero (`FlashTPCBundle.cxx:477-500`). This is the
  *magnitude* comparison.

Important refinements:

- **Cosmic-discriminator masking** (flash `type==1`, i.e. non-beam): a predicted
  PE below the per-PMT threshold `cos_pe_low[j]` (or below `cos_pe_mid[j]*1.1`
  when the measurement is 0) is zeroed before comparison
  (`FlashTPCBundle.cxx:461-464`) — light too dim to fire the discriminator
  shouldn't count.
- **`close_to_PMT` error inflation:** when a PMT measures far more than
  predicted (`pe-pred>350 && pe>1.3·pred`), its error is inflated by
  `(0.5·pe)²` (`FlashTPCBundle.cxx:480-485`). Near the PMTs the model
  under-predicts and cosmic contamination is large, so the chi2 there is
  de-weighted.
- **One inefficient PMT allowed:** if the single worst-chi2 PMT measured 0 while
  predicting >0, that contribution is subtracted off
  (`FlashTPCBundle.cxx:501-502`).

The function then sets `flag_high_consistent` from a ladder of
`(ks_dis, ndf, chi2)` cuts (`FlashTPCBundle.cxx:504-522`). If a bundle isn't
high-consistent, it falls back to a **fired-PMT** test: `nfired/ntot` is the
fraction of predicted-lit PMTs that actually saw signal; too few fired sets
`flag_potential_bad_match` and (for non-boundary clusters) rejects the bundle
(`FlashTPCBundle.cxx:550-600`).

There are three related variants: `examine_bundle()` (self, above);
`examine_bundle(other, ...)` and `examine_bundle_rank(other, ...)` which test
the *combined* prediction of two bundles (`FlashTPCBundle.cxx:38, 153`); and
`examine_beam_bundle()` for beam flashes (`FlashTPCBundle.cxx:383`).

---

## 4. General Q–L filtering strategy

`tpc_light_match` does **not** try to score everything once and threshold. It
runs a long sequence of staged passes that progressively converge toward a
one-to-one(ish) cluster↔flash assignment. In order:

1. **Basic cull** — run `examine_bundle` on every bundle; drop those that fail
   the consistency/fired-PMT test (`ToyMatching.cxx:470-503`).
2. **Bad-match removal** — drop `potential_bad_match` bundles for a cluster, but
   only if that cluster already has another consistent bundle
   (`ToyMatching.cxx:506-553`).
3. **Per-cluster consistency ladders** — set `consistent_flag` via many
   `(ks_dis, ndf, chi2)` cut combinations, with relaxed cuts for
   `flag_close_to_PMT` / `flag_at_x_boundary`; consider the best (`min_bundle`)
   and second-best (`min_bundle1`) bundle by `ks_dis`; tighten "highly
   consistent" ones; remove `ndf==1` bundles when a flash has a unique
   consistent match elsewhere (`ToyMatching.cxx:563-780`).
4. **spec_end & over-tight-boundary cleanup** — if a cluster has both
   `spec_end` and non-`spec_end` consistent bundles, drop the non-`spec_end`
   ones; also drop matches that are *suspiciously* tight at a boundary
   (`ks<0.06 && chi2<3·ndf && at_x_boundary`), which are ambiguous
   (`ToyMatching.cxx:783-826`).
5. **Cross-cluster / cross-flash competition** — when a cluster could match
   several flashes (or vice versa), keep the best and remove the rest using
   `ks_dis` (with a small `+0.003·#bundles` penalty) and `chi2/ndf` comparison
   rules (`ToyMatching.cxx:896-1183`).
6. **Beam-flash rescue** — for beam flashes (`type==2`) left with no consistent
   bundle, give the best one a second chance via `examine_beam_bundle`
   (`ToyMatching.cxx:1246-1282`).
7. **Joint cross-flash check** — a flash with a consistent bundle keeps a
   remaining bundle only if it is *jointly* consistent with a consistent one
   (`examine_bundle(other)`) (`ToyMatching.cxx:1307-1358`).

### 4.1 The global L1 (Lasso) solve

After the heuristic passes, the residual ambiguity is resolved with a sparse
linear solve (`ToyMatching.cxx:1377-1640`):

- Measurement vector `M` = `PE/err` per PMT per flash; the response matrix `R`
  maps each (flash, cluster) unknown to its `pred/err` contribution. Errors fold
  in fudge factors and the relative light-yield error
  (`ToyMatching.cxx:1417, 1585`).
- A penalty matrix `RF` enforces **"a track is used at most once"**
  (`delta_track = 0.01`, `ToyMatching.cxx:1397, 1469-1477`) and flash
  normalization (`delta_flash`, first round only).
- Solved with `WCP::LassoModel` (L1, `lambda=0.1`), so most coefficients are
  driven to zero — i.e. each cluster keeps ≤1 flash
  (`ToyMatching.cxx:1492, 1634`).
- **`at_x_boundary` bundles are down-weighted to 0.2** (`ToyMatching.cxx:1437,
  1603`); init values are seeded from `consistent_flag` (1.0) /
  `at_x_boundary` (0.5) (`ToyMatching.cxx:1498-1506`).
- The first round (with a flash-alone term) drops bundles with `beta<0.05`
  (`ToyMatching.cxx:1518`); the second round re-solves without it. The final
  `beta` gives the per-(flash,cluster) match strength; the strongest flash per
  cluster wins (`ToyMatching.cxx:1649-1663`).

### 4.2 Final organization — `organize_matched_bundles`

`ToyMatching.cxx:1744`. For each flash matched to several clusters, the
`best_bundle` is chosen by the combined score `ks·(chi2/ndf)^0.8`, discounted
×0.8 for `at_x_boundary` and again ×0.8 for `close_to_PMT`
(`ToyMatching.cxx:1807-1819`). Other compatible clusters are merged into it
(`examine_bundle_rank` + `add_bundle`); for beam flashes, geometrically nearby
clusters are merged via `examine_merge_clusters` (`ToyMatching.cxx:1863-1866`,
`FlashTPCBundle.cxx:238`). Second and third rounds re-home or release leftover
bundles so that every cluster ends with at most one flash.

---

## 5. The three hard cases: anode and the readout-window edges

The user specifically asked how the algorithm handles a track **at the anode**
and **at the beginning / end of the readout window**. These map to *distinct*
mechanisms, even though the anode is itself one of the window edges. The
one-line mnemonic:

> **anode = boundary + PMT-proximity; window edges = boundary + trimming + spec_end.**

### 5.1 At the anode (x ≈ 0)

In MicroBooNE the PMTs sit on the anode (low-x) side. A cluster reaching the
anode is therefore both an x-boundary case **and** physically close to the
PMTs. The low-x block sets **both** flags together
(`ToyMatching.cxx:281-287`):

- `flag_close_to_PMT` — drives the chi2 **error-inflation** branch
  (`FlashTPCBundle.cxx:480-485`, the `pe-pred>350 && pe>1.3·pred` rule) and the
  *loosest* consistency ladders throughout `tpc_light_match` (e.g.
  `ToyMatching.cxx:518, 597, 702`). Near the PMTs the light model is least
  reliable, so these matches are judged more leniently.
- `flag_at_x_boundary` — see below.

### 5.2 Beginning / end of the readout window

Because x is reconstructed from drift time via `offset_x`, the **window/drift
edges *are* the x ∈ [0, 256 cm] limits** after the t0 shift. A track that
entered early or late has part of its drift clipped, so its charge — and hence
its predicted light — is incomplete. Three mechanisms cope with this:

1. **End-trimming walk** (`ToyMatching.cxx:176-264`). When the cluster's first
   (or last) time slice pokes just past the anode (or cathode), the code walks
   inward over the leading/trailing time slices counting how many mcells lie
   outside. If only a small fraction are outside (e.g. `≤36` slices and
   `<5%` of mcells, with looser tiers up to `≤60` slices), it accepts snapping
   `first_pos_x` / `last_pos_x` to the boundary — i.e. treats the small
   protrusion as a spurious tail and considers the cluster to *start/end at the
   edge*.
2. **`flag_spec_end`** (`ToyMatching.cxx:207-208, 250-251`). Marks clusters
   whose trimmed end is a "special end". Later, if a cluster has both
   `spec_end` and non-`spec_end` consistent bundles, the algorithm **prefers the
   `spec_end` ones** (`ToyMatching.cxx:790-804`).
3. **`flag_at_x_boundary`** — for clusters touching either boundary. Two
   consequences: (a) the bundle is **down-weighted to 0.2 in the Lasso solve**
   (`ToyMatching.cxx:1437, 1603`); and (b) `examine_bundle` **suspends
   auto-rejection** for a low fired-PMT count — `flag_potential_bad_match` is
   still set, but the bundle is *not* rejected outright the way a non-boundary
   cluster would be (`FlashTPCBundle.cxx:577-600`). Both reflect that truncated
   charge legitimately under-predicts the light, so such matches deserve the
   benefit of the doubt.

---

## 6. Other things worth noting

### 6.1 MicroBooNE-specific hardcoding (most relevant for the port)

- **32 PMTs** everywhere (histogram size, loops, `cos_pe_*` arrays).
- **256 cm** drift / `high_x_cut` (`ToyMatching.cxx:158`).
- **Voxel-id geometry constants** baked into `convert_xyz_voxel_id`
  (`ToyMatching.cxx:25-39`): offsets/pitches and the `75×75×400` voxel grid.
- **`uboone_photon_library.root`** as the visibility source (loaded in the
  disabled `Photon_Library` ctor, `ToyMatching.cxx:99-100`).
- **Library↔PE channel remap** `map_lib_pmt` / `map_pmt_lib` — the library and
  the readout do not share channel ordering (`ToyMatching.cxx:330`).
- **Run/timestamp-dependent light-yield scaling** and the **PMT-17 veto**
  (`ToyMatching.cxx:345-354`; yield table in the disabled ctor,
  `ToyMatching.cxx:124-140`).

### 6.2 Design philosophy

- **Two metrics, both required.** A shape metric (`ks_dis`, normalization-free)
  guards against accepting a match that has the right total light but the wrong
  pattern; a magnitude metric (`chi2/ndf`) guards the converse. Almost every cut
  in the code is a joint `(ks_dis, ndf, chi2)` condition.
- **Staged, greedy → global.** Local consistency and pairwise competition prune
  the obvious cases first; the L1 solve then resolves the genuinely contested
  assignments under the "track used once" constraint. The combined ranking
  score `ks·(chi2/ndf)^0.8` recurs as the tie-breaker
  (`ToyMatching.cxx:1258, 1807`; `FlashTPCBundle.cxx:356`).
- **Lenient where the physics is hard.** Near the PMTs and at the window edges
  the cuts are deliberately relaxed and the global weights reduced, because the
  light model and the charge are both least trustworthy there.
