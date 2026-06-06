# Connecting cathode-crossing tracks across the two TPCs (post-QLMatching clustering)

## The question

SBND has two TPCs split by a central cathode. After Q/L matching assigns a flash T0,
a cosmic that crosses the cathode is reconstructed as two halves — one in TPC0 (x<0),
one in TPC1 (x>0). The user observes some of these as **two separate clusters** and asked:

1. How does the existing post-QLMatching cross-TPC clustering handle the inter-TPC
   **gap** (a physical drift gap near the cathode, *not* dead channels), and how to improve it?
2. Should a **separate, retireable** algorithm be added for the **misalignment** between
   the two halves (a cathode-region calibration distortion), so it retires when the
   calibration improves?
3. (Follow-up) What does `examine_bundles` do after QLMatching?

Scope: only tracks that **end at the cathode** (halves separated by a readout-window gap
are out of scope). The user asked for **evidence from running the code** on whether the
halves are (A) never connected by the merge passes, or (B) connected and then disconnected
by `examine_bundles`.

**Bottom line:** the run evidence says **(A)** — for the failing cases the geometric merge
passes never connect the two halves; `examine_bundles` is *not* the disconnector (it
flash-*merges* them into one bucket, §3). The single binding gate is a **3 cm cap on the
"lenient" collinear-merge path** in `Clustering_1st_round`: a cathode-crosser's closest-point
distance (physical gap ~1.25–1.5 cm + transverse-offset residual) is pushed just past 3 cm
for the worst cases, and above 3 cm the remaining paths demand a connection-vector alignment
the gap+offset destroys. Fix: a separate, default-OFF, retireable cathode-crossing connector
(or, minimally, a larger distance budget on that lenient path for cathode-end cross-TPC pairs).

## 1. How the cross-TPC clustering is wired

The all-APA clustering runs *after* QLMatching, on a single "live" grouping that holds the
clusters of **both** TPCs (`cfg/pgrapher/experiment/sbnd/clus.jsonnet`, `clus_all_apa`,
`cm_pipeline`):

```
switch_scope            // T0-correct every blob into x_t0cor[,y_cor,z_cor]; cathode ≈ x=0
extend       (use_flash_t0)   length_cut 60
regular '1'  (use_flash_t0)   length_cut 60, no extend
regular '2'  (use_flash_t0)   length_cut 30, extend
parallel_prolong (use_flash_t0) length_cut 35
close        (use_flash_t0)   length_cut 1.2
extend_loop  (use_flash_t0)
examine_bundles (use_flash_t0)
```

- `switch_scope` (`clustering_switch_scope.cxx`) rewrites both TPCs into one common
  T0-corrected frame, so a cathode-crosser's two halves meet near x≈0. The frame is
  `common_corr_coords` = `{x_t0cor, y_cor, z_cor}` with `pos_offset_on=true` (SBND default):
  the rigid per-TPC transverse correction (±(−0.11,+0.67) cm, from
  `sbnd_xin/docs/cathode-crossing-diagnostic.md`) is **already applied** in this frame.
- The merge passes (`extend`/`regular`/`parallel_prolong`/`close`/`extend_loop`) each build a
  graph: for every pair of clusters **in the same flash-T0 group** (`assign_flash_t0_groups`,
  ±80 ns) and above an internal length cut (15 cm), a geometric test decides an edge, then
  `merge_clusters` fuses connected clusters into one. The core test is
  `Clustering_1st_round` (`clustering_regular.cxx:63`): closest-point distance
  (`Find_Closest_Points`) + Hough/PCA direction angles. **It does no path / good-point
  validation** — so the bare inter-TPC gap does not block it.

### How the gap itself is treated

- In the **merge passes**: not at all — `Clustering_1st_round` never walks the path between
  the two clusters, so the physical ~1.25–1.5 cm cathode gap is invisible to them. They
  merge across it whenever the endpoints are close and aligned. (The gap matters only through
  the *closest-point distance* it contributes — see §4.)
- In `examine_bundles`' **relaxed graph** (`connect_graph_relaxed.cxx`): the path-quality
  check steps along the line at 1 cm and, for any step that lands *between* the TPCs
  (`get_wireplaneid → apa()==-1`), counts it **fully bad on all planes** (the `else` branch:
  "no signal from any plane can validate this gap") — the opposite of `connect_graph_ctpc`,
  which *skips* such steps. But a ~1.5 cm gap is ~1–2 steps, far below the rejection threshold
  (>7 bad), and this only shapes the post-merge main/associated sub-structure (§3), not the
  merge. So it is not what produces the two-clusters symptom.

## 2. The A/B evidence (run on the 5 hand-scan data crossers)

Data hand-scan cathode-crossers (`sbnd_xin/docs/cathode-crossing-diagnostic.md`):
evt **686, 1302, 1346, 1852, 2028**. Each was run through the SBND all-APA pipeline two ways:
the **full** pipeline and a variant with the all-APA `examine_bundles` **removed** (in the
truncated run `real_cluster_id == cluster_id`, so `cluster_id` *is* the geometric-pass
result). The QLMatching halves were spatially matched (≤0.6 cm) to the all-APA clusters and
their ids compared. Runs + parsers: `/home/xqian/tmp/cathode_ab/`.

| evt | transverse offset¹ | FULL `real_cluster_id` (TPC0/TPC1) | TRUNC `cluster_id` (TPC0/TPC1) | verdict |
|---|---|---|---|---|
| 686  | 2.56 cm | 20 / **34** | 20 / **34** | **NOT connected (A)** |
| 1302 | small  | 65 / 65   | 65 / 65   | merged by geom passes |
| 1346 | small  | 30 / 30   | 30 / 30   | merged by geom passes |
| 1852 | 2.18 cm | 4 / **20** | 4 / **20** | **NOT connected (A)** |
| 2028 | 1.29 cm | 39 / 39   | 39 / 39   | merged by geom passes |

¹ artifact-immune transverse `|perp_yz|` from the QLMatching diagnostic.

Full and truncated agree exactly. **The merge passes already connect 3 of 5** (1302, 1346,
2028) into a single cluster across the gap; they **fail on the two largest-offset cases**
(686, 1852), which stay two distinct geometric clusters. So the answer is **(A): the merge
passes never connect them.** `examine_bundles` does not disconnect — it is downstream of the
failure (§3). (Hypothesis B is structurally impossible: `examine_bundles` only *merges* a
flash group and relabels; it never splits a geometrically-merged cluster.)

## 3. What `examine_bundles` does after QLMatching

In the all-APA scope (`use_flash_t0=true`) `examine_bundles`
(`clustering_examine_bundles.cxx:80-260`) runs **last** and does two things, both *additive*
— it never breaks an existing geometric connection:

1. **Flash-time "bundle" collapse** (`:99-144`). It groups every cluster by its matched
   flash-T0 group (±80 ns) and `merge_clusters` fuses all clusters in a group into **one**
   cluster — so each matched flash becomes a single bundle, *regardless of geometry* (clusters
   that never touch but share a flash are bucketed together). Before merging it records, per
   blob, the **pre-merge** cluster ident into a `real_cluster_id`/`perblob` array so the Bee
   can still colour the far-apart members distinctly.
2. **Per-bundle main + associated labelling** (`:146-240`). For each bundle it runs the
   relaxed connectivity graph (`connected_blobs(dv, pcts, "relaxed")`), partitions the blobs
   into connected sub-components, and tags **one** as the "main" track (`-1`, by overlap with
   the prior main, else the longest) and the rest as "associated", into the `isolated`/`perblob`
   array. This is `clustering_isolated`'s MicroBooNE-style main+associated grouping, but keyed
   by flash time instead of geometry — it is what the downstream group-aware reconstruction
   consumes.

Bee field meaning (`MultiAlgBlobClustering.cxx:1463-1473`): **`cluster_id`** = the bundle
(post-flash-merge) id; **`real_cluster_id`** = the pre-merge geometric id. So two halves
sharing one `cluster_id` only means "same flash bundle," **not** "geometrically connected."

**Consequence for cathode crossers:** for 686/1852 the geometric passes leave two clusters;
`examine_bundles` then flash-merges them into one `cluster_id` *bucket* but keeps their
distinct `real_cluster_id`s — which is exactly why a single physical track paints as two
coloured objects in the Bee. `examine_bundles` is not the cause; the missing geometric merge
upstream is.

## 4. Why 686/1852 fail — the binding gate (instrumented)

`Clustering_1st_round` was instrumented (temporary, reverted) to log, for cross-TPC pairs,
the closest-point distance `dis`, the connection-to-track angles (`angle_diff1/2`), the
direction–direction collinearity (`angle_diff3`), and the per-pass decision. `parallel_prolong`
was instrumented the same way. Logs: `/home/xqian/tmp/cathode_ab/{dbg,pp,pin}_<evt>.log`.

| evt | code `dis` | merges? | how |
|---|---|---|---|
| 1302 | 2.02 cm | yes | `regular` **`dis≤3 cm` lenient path** (`clustering_regular.cxx:248-252`) |
| 1346 | 2.58 cm | yes | same lenient path |
| 2028 | 2.39 cm | yes | same lenient path |
| 686  | **3.32 cm** | no | misses the 3 cm cap → strict paths require conn-angle <12°; it is 52–90° |
| 1852 | **3.33 cm** | no | misses the 3 cm cap → strict paths require conn-angle <12°; it is 22–90° |

The single binding gate is the **`dis ≤ 3 cm` cap on the lenient collinear-merge path**
(`clustering_regular.cxx:248`). That path accepts a pair on **direction collinearity**
(`angle_diff3` < ~18°) plus a *loose* connection-angle tolerance (`angle_diff1/2` < ~54°) —
exactly the regime a cathode-crosser sits in. The three mergers have `dis` 2.0–2.6 cm and
take it. 686/1852 have `dis` 3.32–3.33 cm — over the cap by ~0.3 cm — so they fall through to
the **strict** paths, which require the closest-point *connection vector* to lie within ~12°
of each track. For a cathode crosser the connection is **oblique** (52–90° for 686, 22–90°
for 1852) because the gap+offset displace the two collinear halves laterally, so the strict
paths reject them in every pass.

Two facts this pins down:
- **It is not the gap alone, nor the offset alone — it is the closest-point distance they
  jointly produce** crossing 3 cm. The physical gap (~1.5 cm) already spends half the 3 cm
  budget; the transverse-offset residual spends the rest. (686's loose-path connection angles,
  51.96°/48.97°, sit *just under* the 54° tolerance — so it is *only* the distance cap that
  stops it; a small cap increase would recover it.)
- **The existing code already uses direction collinearity** (`angle_diff3` is the dir–dir
  angle, gated on `clustering_regular.cxx:249/262/278`). The crossers all have excellent
  dir–dir (≤4.5°). What they lack is a *distance* budget large enough for the inter-TPC gap on
  the one path that leans on collinearity. The strict paths additionally demand the
  offset-sensitive connection angle, which is the real obstacle above 3 cm.

The transverse offset is a **data calibration artifact**: the rigid `pos_offset` (already ON)
removes the coherent ~1.34 cm piece, but a position-dependent (SCE-like) residual ~1 cm
remains and is largest for 686/1852. MC's transverse offset is ~0.48 cm (≈⅓ of data, per the
diagnostic), so MC crossers keep `dis` under 3 cm and merge — i.e. this is a data-specific tail.

## 5. Recommendations

### The gap — no gap-specific change

The merge passes ignore the gap (no path validation) and bridge it for 3/5 crossers;
`examine_bundles`' fully-bad-gap penalty is post-merge sub-structure only and below threshold
for the ~1.5 cm cathode gap. **Do not** close the physical x-gap (it is real, ~1.25 cm,
MC≈data) and do not add gap-special path handling. The gap matters only as its contribution to
the closest-point distance (§4), which the distance budget below already covers.

### The misalignment — a separate, default-OFF, retireable cathode-crossing connector

This is the user's preferred form and the right lever. Add a dedicated all-APA pass (own
`clustering_*.cxx` + jsonnet method, fully duplicated, not woven into `Clustering_1st_round` —
`feedback_fork_by_duplication`), default-OFF / bit-identical when off
(`feedback_toggleable_behavior_changes`), placed among the merge passes before
`examine_bundles`:

- **Candidate pairs:** opposite-TPC clusters in the same flash-T0 group, **both** with an
  extreme point near the cathode `|x_t0cor| < ~5 cm` (reuse `get_extreme_wcps` /
  `get_two_extreme_points`, `Facade_Cluster.h:402,431`). The both-ends-at-cathode requirement
  is what excludes the out-of-scope readout-gap pairs.
- **Accept on direction collinearity at a cathode-sized distance budget.** Require the two
  half-directions collinear (the `dis≤3 cm` path's dir–dir test, all crossers ≤4.5°) but allow
  the closest-point distance up to a tuned **~6–10 cm** (covers gap + offset; the failures sit
  at 3.3 cm). Mechanically this is the existing lenient path with the 3 cm cap raised, scoped
  to cathode-end cross-TPC pairs so the generic 3 cm cut is untouched elsewhere.
- **Purity comes from the cathode-end gate + the bounded distance, not from collinearity.**
  The diagnostic warns collinearity alone admits false pairs (an MC false at 3.6° is *more*
  collinear than a true crosser at 2.8°), which is why QLMatching's `flag_xtpc_consistent`
  needed a distance ceiling for 100% purity. Tune the dir–dir and distance cuts on the 10 data
  + 10 MC hand-scan crossers against random same-flash pairs, on the real C++ `vhough` values.
- **Cheaper alternative to weigh:** QLMatching *already* pairs these as
  `flag_xtpc_consistent` (`match/src/QLMatching.cxx`, `cull_cross_tpc`; `d<5 cm`, or collinear
  &`d<300 cm`). Propagating that confirmed-pair flag into the all-APA clustering as a direct
  merge hint would connect the halves with **no new geometry** — the implementer should weigh
  this against a fresh pass.
- **Why separate / retireable:** the gap is permanent geometry (already handled), but the
  misalignment is a *calibration artifact*. As `pos_offset`/SCE transverse calibration
  tightens, the residual offset shrinks, `dis` falls back under 3 cm, the existing lenient path
  merges these crossers, and this connector can be switched off and retired without touching
  production logic.

Expected effect on the hand-scan set: recover 686 and 1852 (3/5 already merge), at no cost to
the others when off.

## Artifacts

- Diagnostic runs + parsers: `/home/xqian/tmp/cathode_ab/` (`full_<evt>.zip`, `trunc_<evt>.zip`,
  `dbg_*/pp_*/pin_*` instrumented logs, `run_*` logs).
- Truncated run = comment `cm.examine_bundles(use_flash_t0=true)` (`clus.jsonnet:302`); gate
  trace = temporary `std::cerr` in `clustering_regular.cxx` / `clustering_parallel_prolong.cxx`.
  Both reverted after measurement.
- Offsets + hand-scan crosser list: `sbnd_xin/docs/cathode-crossing-diagnostic.md`.
- Code: `clustering_regular.cxx:248-287` (the 3 cm lenient path + strict paths),
  `clustering_examine_bundles.cxx:99-240`, `connect_graph_relaxed.cxx` (inter-TPC fully-bad
  branch), `MultiAlgBlobClustering.cxx:1463-1473` (Bee id semantics).
