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
the gap+offset destroys. **Fix implemented and SBND-on:** a separate, retireable
`ClusteringCathodeConnect` pass (§5) — verified to merge exactly the two failing crossers
(686, 1852) and to leave the other 8 of 10 data events byte-identical (no false merges).

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

## 5. Implementation: `ClusteringCathodeConnect` (SBND-on)

### The gap — no gap-specific change

The merge passes ignore the gap (no path validation) and bridge it for 3/5 crossers;
`examine_bundles`' fully-bad-gap penalty is post-merge sub-structure only and below threshold
for the ~1.5 cm cathode gap. **The physical x-gap is real (~1.25 cm, MC≈data) and is not
closed.** The gap matters only as its contribution to the closest-point distance (§4), which
the connector's distance budget covers.

### The misalignment — a separate, retireable cathode-crossing connector

Implemented as a dedicated all-APA pass `ClusteringCathodeConnect`
(`clus/src/clustering_cathode_connect.cxx`, modelled on `clustering_close.cxx` — own file, not
woven into `Clustering_1st_round`, per `feedback_fork_by_duplication`), with jsonnet method
`cathode_connect()` (`pgrapher/common/clus.jsonnet`). It is placed in the SBND all-APA pipeline
**after** the generic merge passes (so it can only ADD merges they missed) and **before**
`examine_bundles` (so a connected crosser is one cluster before the flash-bundle collapse). It
is driven by the toggle `cathode_connect_on` (`cfg/pgrapher/experiment/sbnd/clus.jsonnet`),
**committed ON for SBND**; set it `false` to recover the pre-connector all-APA pipeline
(bit-identical — verified below).

For every pair of clusters in the same flash-T0 group, above `min_length` (10 cm), it accepts
the pair when the following hold (defaults in the jsonnet method). **Always required:**

1. **collinear** — the two half-directions (`vhough_transform` at the closest points, radius
   20 cm) within `angle_cut` (10°, unsigned axis → near 0° or 180°);
2. **separate TPCs** — the two closest points' `wpid().apa()` differ (this is what makes the
   pass incapable of acting within a single TPC);
3. **both ends at the cathode** — `|x − cathode_x| < cathode_x_cut` (3.5 cm) for each closest
   point, in the T0-corrected frame where the cathode is `x≈0`;
4. **same drift depth** — `|x₁ − x₂| < drift_cut` (4 cm). This is the binding cut: between the
   two halves there is only the ~1.5 cm cathode gap plus a drift-x calibration residual
   (measured up to **3.22 cm** for 686), so the closest points sit at nearly the same drift
   depth even when they are offset within the cathode plane.

The **3D closest-point distance** is then handled in two regimes — this is the improvement over
a single 3 cm-style distance cap:

- **close** (`dis < dis_cut`, 5 cm): the p1→p2 connection vector here is **dominated by the
  drift-x offset** (the ~3 cm calibration artifact, nearly along the drift axis), *not* by the
  track. This is precisely why the generic passes miss these crossers: their merge logic
  (`clustering_regular.cxx`) requires the connection vector to align with the track, which the
  drift-dominated connection fails — only the `dis ≤ 3 cm` lenient path (which drops the
  connection test) can save them, and 686/1852 fall just past it. The connector fills that hole
  by accepting on the half-track collinearity (cut 1) alone. *All four sample crossers live here*
  (dis 3.17–3.33 cm).
- **far** (`dis_cut ≤ dis < max_dis`, 5–25 cm): a shallow-angle crosser travels *inside* the
  cathode plane, so the two halves are offset transversely (large `dis`) while still being one
  track. At this longer baseline the connection vector becomes a *reliable* direction, so the
  pass **borrows the generic passes' distance-graded alignment** (`clustering_regular.cxx:209-215`
  — `7.5°` out to 15 cm, `5°` beyond) and requires **both** the (tightened) track collinearity
  **and** the p1→p2 connection vector to align with the track. Two unrelated cosmics that merely
  graze the cathode at the same drift depth cannot fake both at a long baseline, so the relaxed
  distance does not open a false-merge hole. (`max_dis = 25 cm` matches the generic passes'
  long-track regular-merge ceiling.)

Purity comes from cuts 2 + 3 + 4 + the two-regime distance, **not** from collinearity alone (the
diagnostic warns an MC false pair at 3.6° is *more* collinear than a true crosser at 2.8°).
This is the same geometry QLMatching uses to flag cross-TPC pairs (`flag_xtpc_consistent` /
`cull_cross_tpc`), here used to **connect** rather than flag.

**Why retireable:** the gap is permanent geometry, but the misalignment is a *calibration
artifact*. As `pos_offset`/SCE transverse calibration tightens, the drift residual shrinks,
`|x₁ − x₂|` falls back under 3 cm, the generic lenient path catches these crossers, and the
connector can be flipped off and retired without touching production logic.

### Blast radius

The connector lives in the **all-APA** pipeline, which runs only in the standalone
post-QLMatching dev chain (`sbnd_xin/wct-clus-matching-perevt.jsonnet` via `all_apa()`).
**LArSoft production (`wcls-img-clus.jsonnet`) uses `per_volume()` only — it does not run the
all-APA pipeline** — so the SBND-on default does not touch LArSoft production clustering.

### Verification (10 data + 10 MC events, connector ON vs OFF)

Runs in `/home/xqian/tmp/cc_verify/` (ON vs OFF) and `/home/xqian/tmp/cc_geom_*.log` (the
per-pair geometry, from a temporary `std::cerr` that was reverted). The connector's closest-point
geometry for the four crossers — the numbers that set the cuts:

| crosser | dis | **\|x₁−x₂\| (drift)** | transverse | collinear | x₁ / x₂ |
|---|---|---|---|---|---|
| data 686  | 3.32 | **3.22** | 0.77 | 5.9° | −0.68 / 2.55 |
| data 1852 | 3.33 | 2.77 | 1.85 | 2.2° | 2.17 / −0.60 |
| MC 18     | 3.18 | 1.14 | 2.96 | 3.6° | 0.57 / −0.57 |
| MC 42     | 3.17 | 1.34 | 2.88 | 3.2° | 0.67 / −0.67 |

- **`drift_cut` is set from the data:** 686's drift separation is 3.22 cm, so `drift_cut = 4 cm`
  keeps it with margin; `cathode_x_cut = 3.5 cm` keeps 686's `x₂ = 2.55 cm` off the edge. (A
  naive `drift_cut = 3 cm`, as first sketched, would have *rejected* 686.) All four crossers
  have `dis < 5 cm`, so they are accepted in the **close** regime; the far regime adds capability
  but is not exercised by this sample.
- **Fires exactly 4× across the 20 events, all genuine crossers; no false merges.** Comparing the
  output zips member-CRC ON vs OFF, **only** data 686, data 1852, MC 18, MC 42 differ; the other
  8/10 data and 8/10 MC events are byte-identical. The connector only adds edges, so a differing
  event *is* a merge — and only the four true crossers merge.
- **The dangerous case is rejected.** In MC evt42 a *non*-collinear cross-TPC cathode-region pair
  (local dir–dir 27.5°) is **rejected by cut 1** — two unrelated tracks meeting near the cathode
  are not merged.
- **Additive / production-safe.** The toggle wraps the whole pass (`if cathode_connect_on then
  [cm.cathode_connect()] else []`); OFF the all-APA pipeline is structurally identical to the
  pre-connector one, and the 16 byte-identical ON≡OFF events confirm the pass is a no-op where it
  must not fire.

**Scope caveats:** purity is established on the 20 hand-sample events (the population the close-
regime cuts were tuned on); a larger-sample firing-rate scan is the next purity probe. The
**far-regime** distance relaxation (`dis ≥ dis_cut`, with the graded connection-alignment
borrowed from `clustering_regular.cxx`) is **implemented but unexercised** by this sample — no
crosser here has a transverse offset above ~3 cm, so all four merge in the close regime and the
graded far cut never fires. The graded values (`7.5°`/`5°`, `max_dis = 25 cm`) inherit the
generic passes' tuning rather than being fit here; retune when a longer in-plane crosser is
hand-scanned. **Efficiency caveat:** `cathode_x_cut = 3.5 cm` and `drift_cut = 4 cm` catch 686
(x₂ = 2.55 cm, drift = 3.22 cm) with ~0.8–1 cm margin — a crosser with a substantially larger
drift residual would still be missed; this is an efficiency limit, not a purity risk.

## Artifacts

- Implementation: `clus/src/clustering_cathode_connect.cxx`, `cathode_connect()` in
  `cfg/pgrapher/common/clus.jsonnet`, toggle `cathode_connect_on` + pipeline wiring in
  `cfg/pgrapher/experiment/sbnd/clus.jsonnet`.
- Connector verification (ON vs OFF, 10 data events): `/home/xqian/tmp/cathode_connect/`.
- Diagnostic runs + parsers: `/home/xqian/tmp/cathode_ab/` (`full_<evt>.zip`, `trunc_<evt>.zip`,
  `dbg_*/pp_*/pin_*` instrumented logs, `run_*` logs).
- Truncated run = comment `cm.examine_bundles(use_flash_t0=true)` (`clus.jsonnet:302`); gate
  trace = temporary `std::cerr` in `clustering_regular.cxx` / `clustering_parallel_prolong.cxx`.
  Both reverted after measurement.
- Offsets + hand-scan crosser list: `sbnd_xin/docs/cathode-crossing-diagnostic.md`.
- Code: `clustering_regular.cxx:248-287` (the 3 cm lenient path + strict paths),
  `clustering_examine_bundles.cxx:99-240`, `connect_graph_relaxed.cxx` (inter-TPC fully-bad
  branch), `MultiAlgBlobClustering.cxx:1463-1473` (Bee id semantics).
