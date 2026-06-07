# clustering_connect1 isochronous over-merge guard (`iso_max_dis`)

## Problem (SBND evt 183888)

Two genuinely-separate isochronous cosmics (tracks running ⊥ to drift, i.e. at
nearly constant drift time) were merged into one full-height overcluster during
per-APA clustering. User-reported representative points (T0-corrected frame):
`A=(22.0,144.9,244.1)` and `B=(28.5,30.8,291.8)`.

## Diagnosis (where it happens)

The two tracks stay in **separate** clusters through every parallel pass
(`clustering_regular`, `clustering_extend`, `clustering_parallel_prolong`); those
passes only grow each track individually. They are joined in
**`clustering_connect1`** (clustering_connect.cxx), in the PCA-comparison block's
**isochronous-relaxed branch** (the `if` near line 575).

That branch, for two clusters whose PCA axes are both within ~5° of perpendicular
to drift, relaxes the collinearity tolerance from the general `angle_diff < 5°`
to **`angle_diff < 30°`** and accepts when the **infinite PCA-axis line distance**
`dis < 1.5 cm`. For slightly-tilted isochronous axes (~1–3° off perpendicular) the
two extrapolated lines near-intersect, so `dis` is ~0 even when the two point sets
are far apart. The branch has **no actual closest-point requirement**.

Measured for the 183888 merge: `angle_diff=15.6°` (not collinear), `dis(line)=0.22 cm`,
`dis1(center)=117 cm`, both axes ~88–89° from drift, skeleton overlap 10100 pts
(the 150 cm parallel extrapolation bridges the gap). The **real closest-point
distance between the two clusters is 7.29 cm** (almost entirely a drift-x gap).

This is a MicroBooNE-tuned heuristic: when isochronous tracks are broken by
coherent-noise removal it pays to rejoin them generously. SBND's noise removal
damages tracks far less, so the generous bridge over-merges separate cosmics.

## Fix

New config `iso_max_dis` on `ClusteringConnect1` (default `-1` = OFF, byte-identical).
When `> 0`, the isochronous-relaxed branch additionally requires the **real
cluster-to-cluster closest-point distance** (`Find_Closest_Points`) `≤ iso_max_dis`.
The general collinear branch (`angle_diff < 5°`) and all other branches are
untouched, so truly-broken collinear tracks still rejoin.

- `clus/src/clustering_connect.cxx`: factor the inline isochronous sub-clause into
  `bool iso_relax`; gate it on the closest-point distance only when `iso_max_dis > 0`
  (so OFF makes no extra call and is bit-identical).
- `cfg/pgrapher/common/clus.jsonnet`: `connect1(..., iso_max_dis=null)`; key omitted
  when null.
- `cfg/pgrapher/experiment/sbnd/clus.jsonnet`: `cm.connect1(iso_max_dis=5*wc.cm)`
  (per-APA; SBND's only connect1 usage — the all-APA pipeline has no connect1).
  5 cm < the 7.29 cm real gap, above SBND broken-track gaps.

## Verification

- evt 183888 ON (5 cm): A and B land in **different** final clusters
  (A→7839-pt clean track, B→11834-pt). OFF (-1): reproduces the original over-merge
  **exactly** (one 19676-pt cluster).
- 40-event regression (ON vs OFF vs a second OFF run for the known `clus_all_apa`
  nondeterminism baseline): **only 183888 changed**; the other 39 events are
  identical ON≡OFF≡OFF2. The 5 cm value is surgical — no collateral fragmentation.
- Full hand-scan + 3-file data sample (10 data + 10 MC + 3×50 data = 170 events)
  was regenerated with the fix for Bee review; all 170 produced output.
