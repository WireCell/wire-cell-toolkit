# clustering_separate hull cap and big over-clusters (`max_hull_points`)

## Problem (PDVD run 39324 evt 0, drift group 4567)

A steep "vertical" track (large drift-x extent) and two broad isochronous bands
ended up in one 32268-point cluster (id 51 in the per-drift-group bee output).
User-reported representative points: `(-78.3,-187.7,133.0)` and `(13.6,-203.0,95.9)`
(both on the vertical track) and `(-150.6,-96.9,109.6)` (on one iso band).

## Diagnosis

Trace of the per-anode (stage-2) components and the stage-3 (per-drift-group)
merge decisions (temporary instrumentation of every `return true` in the merge
predicates):

- The two isochronous bands **physically cross inside CRP 4**: anode4 cluster 6
  (11137 pts, PCA planarity ratio 0.63 — a 2-D X-shaped blob, not a track) is
  already a spatial merge of both bands at the per-face stage.  No angle or
  distance cut can prevent merging clusters that overlap in space; the
  MicroBooNE design intent is *merge crossing tracks, then `separate` splits
  them*.
- The remaining stage-3 iso↔iso merges are legitimate band continuations across
  CRP boundaries (closest distances 1.6–4.6 cm at the shared boundary):
  - `clustering_regular`'s cross-APA/face long-track branch
    (`length_1,2 > 100 cm && dis < 5 cm`, no angle requirement — the "special
    fix for very long tracks PDHD") joined the band pieces at dis = 1.57 / 3.81 cm.
  - `clustering_regular`'s parallel branch with `para_angle_cut = 60°` at
    `dis < 5 cm` joined further pieces (measured angle_diff 28–40°).
- The **vertical↔iso merge is the one genuinely wrong link**:
  `clustering_extend` flag=4 (dead-region pass), non-parallel
  `is_angle_consistent(dir1,dir2) && is_angle_consistent(dir3,dir2)` branch,
  at closest distance 6.07 cm.  For an isochronous target the wire-view
  projections are degenerate, so the 10° consistency test passes spuriously.

All of this is supposed to be undone by **`clustering_separate`** (stage-3 pass
right after the merges).  It was not, because:

- `JudgeSeparateDec_1` starts from `cluster->get_hull(max_hull_points)`.
- `Cluster::get_hull` (Facade_Cluster.cxx) **returns an empty vector when
  `npoints() > cap`** (default cap `Constants::MaxHullPoints = 10000`).
- An empty hull makes `JudgeSeparateDec_1` return false, so `separate`
  **silently skips exactly the large over-merged clusters that need it most**
  (32268 pts here).

## Fix

`ClusteringSeparate` already has the `max_hull_points` config knob (default `-1`
= the 10000 constant, byte-identical).  SBND has been running with
`max_hull_points=100000` since the full-detector overcluster work; the same
value is now set for **PDVD and PDHD**:

- `cfg/pgrapher/experiment/protodunevd/clus.jsonnet`:
  `cm.separate(use_ctpc=true, max_hull_points=100000)`
- `cfg/pgrapher/experiment/pdhd/clus.jsonnet`: same one-line change.

Result on run 39324 evt 0: the vertical track separates cleanly into its own
cluster (all 1409 points); the crossing iso bands are carved into band-aligned
pieces.  The user's three representative points go from one cluster (51) to two
(53 / 51 / 53 — the first and third points are both on the vertical track).

## Residual limitation

Two *wide* isochronous bands that cross each other are merged at the imaging /
per-face level (they share volume) and `separate`'s path-based split does not
reproduce the two bands perfectly — pieces of both orientations can remain in
one output cluster.  This is a `separate`-algorithm limitation, not a merge-cut
tuning issue: tightening the iso merge cuts cannot help because the bands
genuinely touch.  See `clustering-connect-isochronous.md` for the related (but
different) case of two *separated* isochronous tracks joined by `connect1`'s
infinite-line branch (`iso_max_dis` guard).

## Side effect to be aware of

`Cluster::get_hull` caches its result.  Before the fix, big clusters (>10k pts)
never had hull points cached; after `separate` runs with the raised cap, later
passes that use hull vertices as closest-point candidates
(`Find_Closest_Points`'s candidate builder in clustering_extend.cxx) see the
cached hull for big clusters.  This can slightly change closest-point picks in
passes downstream of `separate` even for clusters that `separate` did not split.
SBND has run with this behavior since `max_hull_points=100000` was introduced
there.
