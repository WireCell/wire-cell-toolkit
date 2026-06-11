# Separation trigger: fiducial-volume insets and drift-side x (PDHD/PDVD)

Why `ClusteringSeparate` never split obvious multi-track over-clusters on
PDHD (run 27409 evts 40900/40904/40908/40912/40924, drift group 0/2), and the
fix: FV insets in the detector configs plus three knob-gated mechanisms in
`clus/src/clustering_separate.cxx` (all default OFF/prototype-exact; PDHD and
PDVD enable them).

## The blindness

`JudgeSeparateDec_2` decides "this cluster has several track endpoints on the
detector surface, separate it" by counting dense hull points **outside** the
configured FV (`num_outside_points`, `independent_surfaces`,
`independent_points`) plus an out-of-time drift-x channel
(`flag_outx`/`num_outx_points`).  The strict angle ladder in the caller
(`flag_top`) likewise compares independent points against `FV_ymax`.

PDHD configured the *exact active envelope* as its FV (y 7.61–606.0 cm,
z 0.234–462.297 cm, x ±357.985 = 4.8 cm beyond the anodes).  Sampled points
reach exactly these values and can never exceed them, so with instrumented
gates every cluster in five diagnosed events showed
`nout=0 noutx=0 nsurf=0` while their hull extremes sat ON the walls
(hy=605.97 vs top 606.0, ly=7.86 vs floor 7.61, …).  Dec_2 was unsatisfiable,
`flag_top` unreachable, and only the narrow no-top angle ladder (near-⊥-beam
clusters) could ever fire.  The same blindness starves `Separate_1` of path
endpoints — it carves between independent points — so even clusters that DID
proceed were carved with only the dense-extreme fallback points.

The reference behavior is uboone's prototype
(`ToyClustering_separate.h`): its hardcoded boundaries sit **~15 cm inside**
the active volume (y>101.5 vs top 116.25, z<15/z>1022 vs 0–1036.8), i.e.
"on the surface" was designed to mean *within ~15 cm of the wall*.  Every
Dec_2 threshold was tuned at that operating point.

## Fix 1 — 15 cm y/z FV insets (config-only)

`cfg/pgrapher/experiment/pdhd/clus.jsonnet` overall block:

| field | active | new FV |
|---|---|---|
| FV_ymin | 7.61 cm | 22.61 cm |
| FV_ymax | 606.0 cm | 591.0 cm |
| FV_zmin | 0.234 cm | 15.23 cm |
| FV_zmax | 462.30 cm | 447.30 cm |

PDVD (`protodunevd/clus.jsonnet`): y ±336.4 → ±321.4 cm, z 0.05–299.25 →
15.05–284.25 cm.  FV_x and all margins unchanged (margins extend outward;
x is the out-of-time channel below).  Per-face blocks define only x and
inherit the overall y/z.  Note `select_scope_fv` is also consumed by
`clustering_neutrino`, which deliberately moves to the same uboone-consistent
operating point.  SBND is untouched (its FV was already inset, plus the
`sbnd_boundary_tag`).

## Fix 2 — `drift_side_fv_x` (knob, default OFF)

`select_scope_fv` falls back to the cryostat envelope whenever the dv spans
more than one APA, discarding the drift-side x-range even when all faces share
it (PDHD group02: x ∈ [−357.985, −2.54] cm).  Out-of-time tracks with
apparent x past the cathode (observed up to +126 cm in group02) should trigger
the uboone no-T0 x channel but the overall ±357.985 hid them.  With
`drift_side_fv_x: true`, a multi-APA scope whose faces all carry identical
`FV_xmin/FV_xmax` metadata (a drift group) keeps that common x-range (and its
margins).  Mixed-face scopes (e.g. SBND all-APA spanning both drift sides)
never agree and keep the overall x even when the knob is on.

## Fix 3 — far-point knobs (default = prototype-exact)

In Dec_2's two-endpoint test (`independent_points.size()==2`), a boundary
point deviating from the endpoint-connecting line promotes separation
(`num_far_points`).  Two prototype details kept real two-track clusters
blocked:

- `far_point_x_cut` (default 140 cm): the prototype expression
  `fabs(dir_3.x()/units::cm) > 14*units::cm` compares a cm-number against
  14·cm = 140, i.e. demands a 140 cm drift-x deviation — effectively dead.
  This is faithful prototype behavior (WCP uses the same mm=1/cm=10 units),
  NOT a toolkit porting bug, so the default reproduces it bit-for-bit and
  PDHD/PDVD set the evidently intended `14*wc.cm`.
- `far_point_mid_dis` (default 25 cm): far points are discarded when the
  midpoint of the two endpoints lies further than this from the cluster.
  Two diverging/forking near-parallel tracks (27409 evt 40904: midpoint
  47.6 cm out, in the empty space between the prongs) lose all their far-point
  evidence.  PDHD/PDVD raise it to `60*wc.cm`.

## Fix 4 — `track_recarve` (knob, default OFF)

`Separate_1`'s path carve can leave both arms of two *crossing* cosmics in one
piece, and the closing `Separate_2` (pure 5 cm connectivity) re-links them at
the crossing — no connectivity step can ever hold an "X" apart.  New
post-separation step (after `collinear_recover`/`band_recarve`, sharing the
family bookkeeping): for each family member >100 cm with
`0.01 < eval1/eval0 < 0.5`, fit the blob centers with a k=2 **3D** line model
(seeded by the 2nd-PCA-axis sign split; 10 fixed assign/refit iterations with
a fixed-count power iteration for the line directions).  Split only when ALL
hold:

- each arm ≥ 20% of npoints and ≥ 60 cm axial extent,
- each arm thin: npoints-weighted rms residual ≤ 10 cm (blob centers scatter
  several cm on genuine arms; a wide isochronous band "fits" two parallel
  lines instead and is vetoed by the crossing test, not by this),
- the two fitted lines genuinely cross: closest approach ≤ 10 cm, located
  interior ([0.15, 0.85] of the axial extent) to **at least one** arm and
  within [−0.1, 1.1] of both.  A kinked single track (a "V") meets at the END
  of both arms and is rejected; an "X" (both interior) or a "T" (one track
  ends on the other, which it crosses mid-length) splits.

On 27409 evt 40908 this fires exactly once (the t1/t3 X: arms 2584/8352
npoints, cross fractions 0.93/0.42, residuals 5.3/8.7 cm).  It fires nowhere
else in the PDHD/PDVD regression set.

## Fix 5 — `dec1_guard_main_angle` (knob, default = legacy)

`JudgeSeparateDec_1` carries a toolkit-added protection guard (e28db401, not
in the prototype): `angle1 < 5 && ratio2 < 0.05 && length > 300 cm` returns
false to protect long thin drift-aligned through-going cosmics.  But `angle1`
tests the **2nd** PCA axis against perpendicular-to-drift, which a wide
isochronous complex satisfies trivially — a signal-processing rerun nudged the
PDHD floor-to-top giant to r2 = 0.031 and the PDVD band complex (which lost
its out-of-time appendage to slightly different imaging) into the guard's
window, silently vetoing exactly the multi-track over-clusters separation
exists for.  With `dec1_guard_main_angle` set (degrees), the guard
additionally requires the cluster MAIN axis within that angle of the drift
axis (the actual through-goer topology).  Default <0 keeps the unconditional
legacy guard bit-for-bit; PDHD/PDVD set 45.

## track_recarve standalone trigger

Two crossing tracks whose arms end INSIDE the volume expose at most one
surface contact, so neither Dec_2 nor the angle ladders can ever fire on them
(27409 evt 40908 after the SP rerun: the X host had one upstream contact,
r1 = 0.556, 46° from ⊥-beam).  When `track_recarve` is on, non-proceeding
clusters >100 cm are offered directly to the k=2 3D-line splitter — its own
validation (two long thin arms, fair shares, interior crossing) does the
actual discrimination, so the loose trigger is safe: on the full PDHD+PDVD
regression it fires only on the genuine X/T topologies.  RATIO_MAX is 0.8
(a wide-opening X reaches ~0.56).

## band_recarve seed gates (revision)

With separation now firing on track-topology events, `recarve_two_bands`
needed gates so its x-blind (y,z) 2-line re-carve never pools thin crossing
line-tracks (it had re-mixed the two crossing tracks of evt 40900):

- `SEED_WIDTH_MIN = 6 cm`: both seed members must be ribbon-like in rms width
  across the main axis (39324 seed pair: 17.1/6.2 cm; the PDHD track slices
  that must not seed: ≤ 4.6 cm at the qualifying angles).
- `SEED_ANG_MIN` raised 8° → 15°: surviving crossing-track pairs pair up at
  ~11° (PDHD evt 40900: 16.0/10.0 cm widths at 10.9°), genuine two-band seeds
  at 25–39° (39324 evt 0: 25.3°, evt 8: 38.8°).

The 39324 evt 0 invariants (three user pairs in one cluster each, band
complex = exactly 2 clusters, vertical track whole) still pass.

## Configuration

```jsonnet
cm.separate(use_ctpc=true, max_hull_points=1000000,
            collinear_recover=true, collinear_interior=true,
            collinear_member_merge=true,
            track_repartition=true, band_merge_back=true,
            band_recarve=true, drift_side_fv_x=true,
            far_point_x_cut=14*wc.cm, far_point_mid_dis=60*wc.cm,
            track_recarve=true, dec1_guard_main_angle=45,
            iso_slab_split=true),
```

(`collinear_interior` — interior fragment reclaim, added when the DNN-SP
imaging inputs exposed a mid-track bite on 27409 evt 40900;
`track_repartition` / `band_merge_back` — crossing-track segment swap fix and
hatched-band re-assembly, added for the PDVD 39324 over-separation round;
`collinear_member_merge` / `iso_slab_split` — rejoin a carve-cut straight
cosmic and x-slab-split mixed band+drift-track clusters, added in round 5;
see `clustering-separate-refine.md`.  `max_hull_points` was raised 100k → 1M
after a 102,129-pt giant in 39324 evt 339890 slipped past the 100k cap and
silently skipped separation — the same silent-skip failure the knob was
introduced for.)

All new keys use the `[if flag/value]` omit-when-default pattern: existing
configs are bit-identical (verified: new binary + old configs reproduce
027409 evt 40908's three zips md5-identically).

## Diagnosed events (PDHD run 27409, group 0/2; ids from the pre-fix global)

| evt | old cluster | topology | blocked by | after fix |
|---|---|---|---|---|
| 40900 | 23, 26 | iso band + 3 tracks (one 5795-blob cluster) | carve starved (6 fallback endpoints); band_recarve mis-pooled the tracks | 3 tracks in 3 clusters; recarve correctly silent |
| 40904 | 31 | two forking near-parallel tracks | Dec_2 `nindep=2, nfar=0` (middle_dis 47.6 > 25) | split (far_point_mid_dis) |
| 40908 | 27 | 3 tracks, two crossing as an X | ladder 43°; Dec_2 blind; X re-linked by Separate_2 | 3 tracks in 3 clusters (track_recarve) |
| 40912 | 16 | 2 tracks crossing near top (y_max 8.2 cm below wall) | `flag_top` unreachable | split (top ladder via inset) |
| 40924 | 22 | 3 tracks | ladder 36° > 35; Dec_2 blind | 3 tracks in 3 clusters |

## Regression

- OFF-check: new binary + old configs → byte-identical (027409 evt 40908, all
  three mabc zips md5-equal to the archived baselines).
- Knobs-ON sweep, clustering-only rerun against archived pre-change outputs
  (PDHD 027409 evts 0–7,12 + 027380 evts 0–7; PDVD 39252/39253 evts 0–4 +
  39324 evts 0–10): nearly every event changes, as expected for an
  operating-point change; flow analysis shows modest reorganizations only
  (typically +1–2 clusters/event from 1–2 new splits, the largest the
  five diagnosed events and 39324 evt 4 at 5 splits), no shredding.
- A handful of display points (≤138 on a 60k-point cluster, ≤0.3%) can leave
  or enter the Bee dump entirely when a small carve fragment crosses the
  isolated-cluster trash threshold; the underlying blobs are unchanged.
- 39324 evt 0 refinement invariants re-verified after every tuning step.
