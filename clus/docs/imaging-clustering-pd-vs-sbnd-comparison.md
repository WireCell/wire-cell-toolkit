# Imaging & Clustering Chain: ProtoDUNE-HD / ProtoDUNE-VD vs SBND

**Status:** analysis / reference, no code change.
**Date:** 2026-06-06.
**Author:** generated comparison, verified against source `file:line`.

## Motivation

We are porting the SBND-style Wire-Cell imaging + clustering chain to ProtoDUNE-HD
(`./pdhd`) and ProtoDUNE-VD (`./pdvd`). ProtoDUNE does **not yet have Charge–Light (Q/L)
matching**; SBND (`./sbnd_xin`) is the most mature reference chain and *does*. The goal is
to **improve and validate the PD imaging/clustering chain before adding matching**, so this
document maps the three chains step-by-step and isolates:

1. what is identical (and can be trusted),
2. what differs only parametrically/geometrically,
3. the structural differences in clustering, and
4. which SBND improvements are portable to PD **now** (matching-independent) vs which are
   fundamentally **gated on Q/L T0** (and should wait).

### Sources verified

| Role | SBND | PDHD | PDVD |
|---|---|---|---|
| Imaging | `cfg/pgrapher/experiment/sbnd/img.jsonnet` | `pdhd/img.jsonnet` | `pdvd/img.jsonnet` |
| Clustering | `cfg/pgrapher/experiment/sbnd/clus.jsonnet` | `pdhd/clus.jsonnet` | `pdvd/clus.jsonnet` |
| Standalone entry | `sbnd_xin/wct-clustering.jsonnet` | `pdhd/wct-clustering.jsonnet` | `pdvd/wct-clustering.jsonnet` |

`sbnd_xin/clus.jsonnet` is a thin re-export of the canonical
`cfg/pgrapher/experiment/sbnd/clus.jsonnet`. The clustering-algorithm factory shared by all
three detectors is `cfg/pgrapher/common/clus.jsonnet` (`clustering_methods()`); the per-step
defaults (e.g. `use_flash_t0=false`, `iso_max_dis=null`) quoted below come from there.
`pdhd/pdvd/sbnd_xin` are symlinks into the separate `wcp-porting-img` repo; this doc is a
toolkit-repo artifact.

---

## 1. Imaging

### 1.1 Shared skeleton (identical across all three)

All three imaging graphs are the same WCT subgraph, per anode:

```
(per-anode gauss%d / wiener%d frames)
  → ChannelSelector / CMMModifier / FrameMasking          (per-anode bad-channel prep)
  → ChargeErrorFrameEstimator                              (charge uncertainty)
  → MaskSlices                                             (traces → slices)
  → SliceFanout → GridTiling(face0) + GridTiling(face1) → BlobSetSync   (per-face tiling)
  → { ACTIVE fork:  multi_active_slicing_tiling  → solving(full_deghost) }
  + { MASKED fork:  multi_masked_2view_slicing_tiling → BlobClustering   }
  → ClusterFileSink   (active stream + masked stream)
```

The **charge-solving** sub-pipeline is identical in all three (`solving()` with
`full_deghost=true`):

```
BlobClustering(policy="uboone")
  → ProjectionDeghosting (×2, before 1st and 2nd ChargeSolving)
  → [ BlobGrouping → ChargeSolving(uniform) → LocalGeomClustering → ChargeSolving(uboone) ]  ×3 rounds
  → InSliceDeghosting (×3, one after each ChargeSolving round; good_blob_charge_th=300)
  → GlobalGeomClustering(policy="uboone")
```

Parameters identical across all three:

- `MaskSlices`: active `tick_span=4`, masked `span=500`; `nthreshold=[3.6, 3.6, 3.6]`
  (per-plane SNR threshold).
- Active multi-pass plane sets: `[0,1,2], [0,1], [1,2], [0,2]` (complementary masked
  `[], [2], [0], [1]`), merged with `BlobSetMerge`.
- Masked 2-view multi-pass: dead planes `[0,1], [1,2], [0,2]`, dummy `[2], [0], [1]`.
- `ChargeErrorFrameEstimator`: `rebin=4`, `fudge_factors=[2.31, 2.31, 1.1]` (U, V, W),
  `time_limits=[12, 800]` ticks.
- `GridTiling`: `nudge=1e-2`.
- `ChargeSolving`: `solve_config="uboone"`, `whiten=true`.

### 1.2 Differences (parametric / geometric, not structural)

| Aspect | SBND | PDHD | PDVD |
|---|---|---|---|
| Anodes (APAs) | **2** (1 TPC each) | **4** | **8** (4 bottom-drift + 4 top-drift) |
| Faces per anode | **1** (face 0 only) | **2** | **2** |
| Drift geometry | horizontal, ±X | horizontal, ±X | **vertical-drift CRP**, X-drift |
| Channels / anode | 5638 | (HD geom) | (VD geom) |
| Charge-error file | `sbnd-charge-error.json.bz2` | `microboone-charge-error.json.bz2` | `microboone-charge-error.json.bz2` |
| `MaskSlices max_tbin` | 3427 | 8500 | 8000 |
| Drift speed | 1.563 mm/µs | 1.6 mm/µs | 1.6 mm/µs |
| Tick | 0.5 µs | 0.5 µs (sim) / 512 ns (data, resampled) | 0.5 µs |
| Default `multi_slicing` | `multi-3view` | `multi` (also `single`, `pdhd1` variants) | `multi` (dual active+masked fork) |
| `full_deghost` | true | true | true |
| Dead-region registry | **yes** (see 1.3) | **no** | **no** |

The number of anodes/faces is the only difference that propagates into the *structure* of
clustering (Section 2). The charge-error file and `max_tbin` are calibration/readout
specifics. Everything else in imaging is shared.

### 1.3 SBND-only: the dead-gap registry

SBND has `cfg/pgrapher/experiment/sbnd/dead_regions.jsonnet`, describing one fixed defect:
the middle 6-channel W-plane dead column (TPC0 W channels `[4800, 4806)`, TPC1 mirror
`[10438, 10444)` at the +5638 per-APA offset), declared dead for the full drift extent with
`gap: true`. It is injected into per-APA clustering via `inject_dead_winds` in
`PointTreeBuilding` (`cfg/pgrapher/experiment/sbnd/clus.jsonnet:167`), routing the W winds
into the `Grouping` dead-gap registry so the relaxed connectivity graph bridges the gap
instead of fragmenting a cosmic that crosses it.

**PD has no analog.** PD relies only on the generic masked (dead) imaging fork. If PD has a
comparable systematic dead column, a `dead_regions.jsonnet`-style registry would be the
SBND-proven way to prevent fragmentation — but this is a per-detector geometry question to
confirm before porting.

---

## 2. Per-APA clustering — the biggest structural difference

### 2.1 Tiering: 2-tier (SBND) vs 3-tier (PD)

The tiering is driven by **face count**:

- **SBND** anode = a single drift face, so clustering is **2-tier**:
  per-APA (`clus_per_face` invoked with `face=0`) → all-APA.
- **ProtoDUNE** anode = **two** drift faces (wires read from both sides), so clustering is
  **3-tier**: per-**face** → per-APA **merge** (combine the two faces) → all-APA.

```
SBND:   per-APA (face 0)  ───────────────►  all-APA
PD:     per-face0 ┐
                  ├─ PointTreeMerging ► per-APA-merge ► (×N APAs) ► PointTreeMerging ► all-APA
        per-face1 ┘
```

Source: `pdhd/clus.jsonnet:123` `clus_per_face`, `:258` `clus_per_apa` (with
`clus_per_face(face=0)` + `clus_per_face(face=1)` at `:282-283`), `:372` `clus_all_apa`.
PDVD mirrors this (`pdvd/clus.jsonnet:133/268/382`). SBND only has `clus_per_face`
(`cfg/.../sbnd/clus.jsonnet:149`) and `clus_all_apa` (`:279`); its `per_apa` method just
calls `clus_per_face(face=0)`.

### 2.2 Per-face / per-APA pipeline (`cm_pipeline`)

This is the single-anode (SBND: single-TPC; PD: single drift face) algorithm list. Steps
1–10 are identical; SBND adds six more (11–16) that PD pushes downstream into later tiers.

| # | Step | SBND per-APA (16) | PDHD per-face (11) | PDVD per-face (11) |
|---|---|---|---|---|
| 1 | `pointed` | ✓ | ✓ | ✓ |
| 2 | `live_dead` (offset=2) | ✓ | ✓ | ✓ |
| 3 | `extend` (60/15 cm) | ✓ | ✓ | ✓ |
| 4 | `regular-one` (60 cm, no-extend) | ✓ | ✓ | ✓ |
| 5 | `regular-two` (30 cm, extend) | ✓ | ✓ | ✓ |
| 6 | `parallel_prolong` (35 cm) | ✓ | ✓ | ✓ |
| 7 | `close` (1.2 cm) | ✓ | ✓ | ✓ |
| 8 | `extend_loop` (num_try=3) | ✓ | ✓ | ✓ |
| 9 | `separate` | `use_ctpc, max_hull_points=100000, sbnd_boundary_tag=true` | `use_ctpc` only | `use_ctpc` only |
| 10 | `connect1` | `iso_max_dis=5cm` | (defaults) | (defaults) |
| 11 | `deghost` | ✓ | ✗ → per-APA merge | ✗ → per-APA merge |
| 12 | `examine_x_boundary` | ✓ | ✗ (commented in all-APA) | ✗ (commented in all-APA) |
| 13 | `protect_overclustering` | ✓ | ✗ → per-APA merge | ✗ → per-APA merge |
| 14 | `neutrino` | ✓ | ✗ → all-APA | ✗ → all-APA |
| 15 | `isolated` | `length_cut=15cm` | ✗ → all-APA (default 20 cm) | ✗ → all-APA (default 20 cm) |
| 16 | `examine_bundles` | `relaxed` | ✗ → all-APA | ✗ → all-APA |

Source: SBND `cfg/.../sbnd/clus.jsonnet:177-211`; PDHD `pdhd/clus.jsonnet:190-204`; PDVD
`pdvd/clus.jsonnet:200-214`.

### 2.3 PD-only: the per-APA "merge" tier (steps that combine two faces)

After per-face clustering, PD merges the two faces with `PointTreeMerging(multiplicity=2)`
and runs a short pipeline:

- **PDHD** (`pdhd/clus.jsonnet:303-306`): `deghost()` + `protect_overclustering()`.
- **PDVD** (`pdvd/clus.jsonnet:313-316`): `protect_overclustering()` only —
  **`deghost()` is commented out** (`pdvd/clus.jsonnet:314`).

So PDVD currently performs **no per-APA deghosting at all** before the all-APA stage. This
is worth flagging as a likely regression/oversight relative to PDHD and SBND.

### 2.4 Takeaways (per-APA)

- **SBND front-loads semantic cleanup into the single-TPC pass.** `deghost`,
  `examine_x_boundary`, `protect_overclustering`, `neutrino`, `isolated`, and
  `examine_bundles` all run where the geometry is unambiguous (one TPC, one drift solution).
  PD defers `deghost`/`protect` to the face-merge and pushes
  `neutrino`/`isolated`/`examine_bundles` all the way to all-APA — where, in PD, the drift
  coordinate is only approximately known (Section 3).
- **SBND-specific guards absent in PD and portable now (matching-independent):**
  - `sbnd_boundary_tag=true` in `separate` — topology tag to split two beam-inclined
    crossing cosmics (SBND-tuned; the band geometry would need a PD analog before enabling).
  - `iso_max_dis=5cm` in `connect1` — closest-point guard against isochronous over-merge
    (merging tracks that are ~7 cm apart but collinear in the infinite-line metric).
  - `max_hull_points=100000` in `separate` — raise the hull-point cap so full-detector
    overclusters are not silently skipped (C++ default 10000 → empty hull → no separation).
  - `isolated length_cut=15cm` (vs default 20 cm) — prevent ~16 cm EM blobs from being
    classed "small" and absorbed by the 80 cm small→big merge.

  All of these are pure clustering-geometry knobs (no flash/T0 dependence) and are the
  highest-value matching-independent improvements to evaluate on PD.

---

## 3. Whole-detector (all-APA) clustering

### 3.1 Pipeline comparison

| # | Step | SBND all-APA | PDHD all-APA | PDVD all-APA |
|---|---|---|---|---|
| 1 | `switch_scope` (T0Correction) | ✓ | ✓ | ✓ |
| 2 | `extend` | `use_flash_t0=true` | default (false) | default (false) |
| 3 | `regular-1` (60 cm) | `use_flash_t0=true` | false | false |
| 4 | `regular-2` (30 cm) | `use_flash_t0=true` | false | false |
| 5 | `parallel_prolong` (35 cm) | `use_flash_t0=true` | false | false |
| 6 | `close` (1.2 cm) | `use_flash_t0=true` | false | false |
| 7 | `extend_loop` (num_try=3) | `use_flash_t0=true` | false | false |
| 8 | `separate` | ✗ | ✓ `use_ctpc` | ✓ `use_ctpc` |
| 9 | `neutrino` | ✗ | ✓ | ✓ |
| 10 | `isolated` | ✗ | ✓ | ✓ |
| 11 | `cathode_connect` | ✓ (cross-TPC, all-APA only) | ✓ (`use_flash_t0=false`, 2026-06-09) | ✓ (`use_flash_t0=false`, 2026-06-09) |
| 12 | `examine_bundles` | `use_flash_t0=true` | default | default |

Source: SBND `cfg/.../sbnd/clus.jsonnet:309-323`; PDHD `pdhd/clus.jsonnet:407-420`; PDVD
`pdvd/clus.jsonnet:416-429`.

Coordinate frame:

- **SBND** all-APA: `['x_t0cor', 'y_cor', 'z_cor']` (`pos_offset_on=true`, symmetric per-TPC
  transverse shift: TPC0 `[0, -0.11, +0.67] cm`, TPC1 `[0, +0.11, -0.67] cm`;
  `cfg/.../sbnd/clus.jsonnet:35-40`).
- **PD** all-APA: `['x_t0cor', 'y', 'z']` (`common_corr_coords`, `pdhd/clus.jsonnet:17`).
  Per-face/per-APA tiers use uncorrected `['x','y','z']` in all three detectors.

### 3.2 THE CENTRAL FINDING — T0 / drift coordinate

This is the difference that matters most for the user's "validate before matching" goal.

- **SBND** has a real **per-cluster T0 from Q/L matching** (the matched optical flash).
  `switch_scope` reads it via `T0Correction`, and **every** merge step sets
  `use_flash_t0=true` (`cfg/.../sbnd/clus.jsonnet:311-323`), so a merge only joins clusters
  that are coincident in flash time (`flash_t0_window`, default 80 ns). `x_t0cor` is
  therefore physically meaningful, which is exactly what makes cross-TPC merging,
  `cathode_connect`, and the `pos_offset` transverse correction valid.

- **PDHD / PDVD** run `switch_scope → T0Correction` **with no flash association** (no Q/L
  yet). All all-APA merge steps keep the factory default `use_flash_t0=false`
  (`cfg/pgrapher/common/clus.jsonnet`: `extend/regular/parallel_prolong/close/extend_loop`
  all default `use_flash_t0=false`). PDHD `clus.jsonnet:469` says so explicitly:

  > `coords: ["x", "y", "z"],  // Coordinates to use (uncorrected; x_t0cor needs flash-associated t0)`

  So in PD, `x_t0cor ≈ x` under an **assumed-prompt (t0 ≈ 0) drift**. Out-of-time cosmics
  are misplaced along X, and PD's all-APA cross-anode merging therefore happens at a
  **fixed assumed T0**. PD compensates by running the heavier `separate`/`neutrino`/
  `isolated` here (which SBND already did per-TPC), because the per-face passes were lighter.

### 3.3 Implication for the user's goal

- **Per-APA / per-face imaging and clustering quality is improvable & validatable
  independently of matching.** Porting the SBND per-APA enrichment (move `deghost` etc.
  earlier where geometry is safe) and the SBND guards (`iso_max_dis`,
  `separate max_hull_points`, `isolated length_cut`, and a PD analog of `sbnd_boundary_tag`)
  is the high-value, matching-independent work to do tomorrow.
- **All-APA cross-anode merge correctness for out-of-time activity is fundamentally gated
  on T0** — which is precisely what Q/L matching will supply. Validate it *qualitatively*
  now (it works for in-time/prompt activity), but defer the flash-gated machinery
  (`use_flash_t0`, `cathode_connect`, `pos_offset`) until matching lands.

---

## 4. Recommendations for tomorrow (ordered, matching-independent first)

1. **Re-enable per-APA `deghost` in PDVD** — currently commented (`pdvd/clus.jsonnet:314`);
   PDVD does no per-APA deghosting before all-APA. Lowest-risk, likely a real gap vs PDHD.
2. **Evaluate porting SBND per-APA guards to PD** (per face/per-APA tier), each as a
   jsonnet toggle defaulting OFF so existing configs stay bit-identical (repo convention):
   - `connect1 iso_max_dis` (isochronous over-merge guard),
   - `separate max_hull_points=100000` (don't skip large overclusters),
   - `isolated length_cut` (tighten the small→big absorption),
   - `separate sbnd_boundary_tag` — only after defining the PD upstream-boundary geometry.
3. **Consider moving `deghost`/`examine_x_boundary`/`protect_overclustering` earlier in PD**
   (into the per-face or per-APA-merge tier) to match SBND's "clean per-TPC first" ordering.
   This is a behavior change → gate behind a toggle and A/B on a regression event set.
4. **Consider a PD dead-gap registry** analog to SBND `dead_regions.jsonnet` if PD has a
   systematic dead wire column worth bridging.
5. **Defer** the flash-gated all-APA machinery (`use_flash_t0=true` on merges,
   `pos_offset`) until Q/L matching provides per-cluster T0. Until then,
   PD all-APA is valid for prompt activity and approximate for out-of-time cosmics.
   *(Update 2026-06-09: `cathode_connect` is now enabled on PDHD/PDVD with
   `use_flash_t0=false` — the geometric conjunction alone gates pairing, which only
   admits near-trigger-time crossers on a no-T0 detector; see
   `clus/docs/cathode-crossing-clustering.md` §6.)*

### Validation harness (already present per detector)

Each PD dir has `run_img_evt.sh`, `run_clus_evt.sh`, `run_bee_img_evt.sh`,
`wct-img-all.jsonnet`, `wct-clustering.jsonnet` mirroring the SBND scripts — so any toggle
above can be A/B'd per event and uploaded to Bee for visual comparison, the same workflow
already used for the SBND tuning work.
