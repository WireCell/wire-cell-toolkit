# SBND clustering timing profile (lan-reco2 reprocess) — optimization targets

**Purpose.** Per-event, per-stage wall-time breakdown of the SBND all-APA clustering
chain, taken from the `input-3files-lan-reco2` reprocess (150 real BNB+cosmic data
events). Identifies *which clustering stage* is slow and *which events* dominate, so
the optimization work can be targeted. **Numbers below are LOCKED from the final
150-event auto_mask-on reprocess (2026-06-04, `f{1,2,3}on/match.log`).**

> **⚠ Provenance (2026-06-05): the §"Headline" clustering numbers below were measured on
> yuhw's pre-imaged `*-active.npz` (foreign live input), which we have since stopped using
> ([[feedback_own_imaging_only]]). Our own local multi-3view active imaging changes the
> clustering input, so those per-stage/per-event clustering shares are expected to shift and
> should be treated as FOREIGN-INPUT BASELINE pending a re-pull from a crash-clean
> local-imaging run. The new §0 imaging numbers below ARE from our own local imaging.**

## 0. Imaging time & memory (our own local imaging, 2026-06-05)

Measured on file 1 of `input-3files-lan-reco2` (50 real data events), imaging active+dead
blobs ourselves from the SP frames — no foreign npz. Two separate wire-cell processes,
each wrapped by `profrun.sh` (polls `/proc/PID/status` `VmHWM` for peak RSS); wall from the
process log timestamp span.

| stage | what | wall (50 evt) | per-event | peak RSS |
|---|---|--:|--:|--:|
| **active imaging** | `wct-img-all.jsonnet` (multi-3view + full_deghost, both APAs) | **~382 s** | ~7.6 s | **1131 MB** |
| **matching process** | dead imaging + per-APA + all-APA clustering + QLMatching | ~658 s | ~13 s | 985 MB |

The active-imaging job (multi-3view + `full_deghost` = ProjectionDeghosting ×2 + 3×
ChargeSolving + 3× InSliceDeghosting per APA) is the heavier standalone imaging step and is
where W-plane dead-channel charge recovery happens. Dead/masked imaging is folded into the
matching process (in-graph) and is not separately timed here. Within the matching process,
**QLMatching itself is ~369 s (~56 %)** — see `match/docs/chisquare_flags_comparison.md` §17
for the per-event matcher breakdown.

**Per-stage imaging timing is still NOT instrumented** (the imaging nodes emit no `took`
lines — see the caveat below); the §0 numbers are whole-process only. Adding `took` timing
to `MaskSlices`/`GridTiling`/`BlobClustering`/`ChargeSolving` remains the prerequisite for a
per-stage imaging breakdown.

> **Open blocker (2026-06-05):** the full local-imaging reprocess of files 1–3 is blocked by
> an **intermittent clustering heap-corruption** in `clus_all_apa`, newly exposed by the
> richer local active imaging (segfaults on a layout-dependent ~50–70 % of runs; clean under
> gdb and in single-event isolation → a memory-corruption heisenbug, distinct from the
> already-fixed `fill_wrap_points` off-by-one). Root-causing via valgrind memcheck (Go/
> gojsonnet startup noise suppressed). Until fixed, the §0 wall/RSS come from runs that
> happened to complete all 50 events.

## How it was measured

Each `MultiAlgBlobClustering` node logs `MABC timing: <Stage>:<scope> took <ms> ms`
per clustering method per event. The parser (`sbnd_xin/clustering_timing_analysis.py`)
attributes every timing line to its MABC node instance (`apa0-0`, `apa1-0` = per-TPC;
`clus_all_apa` = all-APA merge) and to the event ident currently loaded in that node,
then aggregates. Run:

```
python3 timing_analysis.py f1/match_f1_fixed.log f2on/match.log f3on/match.log
```

(one non-overlapping log per file: file 1 from the complete `f1off` run, files 2–3
from `f2on`/`f3on`). Re-run against the final `f{1,2,3}on/match.log` once the
auto_mask-on reprocess finishes for the locked numbers.

**Caveat — imaging is NOT in this profile.** The live (active) clusters are read from
pre-imaged `icluster-apa*-active.npz`; the dead (masked) view is imaged in-graph but
emits no `took` lines. So this is a *clustering* profile only. To profile imaging
per-event, the imaging stages (`MaskSlices`/`GridTiling`/`BlobClustering`/
`PointTreeBuilding`) need the same `took`-style instrumentation MABC has — **a small
follow-up worth doing** before optimizing imaging, since right now we are blind to it.

## Headline: one stage and a handful of events dominate

Total clustering wall summed over all 150 events: **~2589 s**.

| stage | total | share | worst single event |
|---|--:|--:|---|
| **ClusteringExamineBundles** | **1536 s** | **59 %** | **310 s** (ev 141530, apa0-0) |
| ClusteringProtectOverclustering | 432 s | 17 % | 171 s (ev 141530) |
| ClusteringDeghost | 366 s | 14 % | 164 s (ev 141530) |
| ClusteringSeparate | 123 s | 5 % | 73 s (ev 59789) |
| ClusteringExtendLoop | 41 s | 2 % | 1.1 s |
| ClusteringConnect1 | 35 s | 1 % | 1.0 s |
| (all others) | ~58 s | 2 % | <0.5 s each |

**Top 4 stages = ~95 % of all clustering time.** Everything below ExtendLoop is noise.

### The work is concentrated in a few pathological events

| rank | event ident | file | total clustering | dominant stage |
|--:|--:|:--:|--:|---|
| 1 | **141530** | 2 | **952 s** | ExamineBundles 613 s |
| 2 | **59789** | 1 | **466 s** | ExamineBundles 239 s |
| 3 | **138824** | 2 | **236 s** | ExamineBundles 150 s |
| 4 | 187650 | 3 | ~92 s | ExamineBundles 46 s |
| 5 | 61637 | 1 | ~72 s | ExamineBundles 48 s |
| … | | | | |

**One event (141530) = ~37 % of the entire clustering wall; the top 3 events = ~64 %.**
The median event is a few seconds; the mean is dragged up entirely by this tail. Every
slow event is ExamineBundles-dominated — the busy cosmic pile-up events where the
per-TPC cluster/blob count blows up. (Note the worst event shifted from 59789 to 141530
vs the mid-run snapshot: the tail is driven by individual dense events, so the exact
ranking moves but the shape — one stage, a few events — is invariant.)

### Per-node (per-MABC-instance) share

| node | share | note |
|---|--:|---|
| `apa0-0` (TPC 0) | ~46 % | consistently slower than TPC 1 — and it owns the worst event 59789 |
| `clus_all_apa` (all-APA merge) | ~28 % | second full clustering pass over the merged both-TPC tree |
| `apa1-0` (TPC 1) | ~26 % | |

## Why ClusteringExamineBundles is the hotspot (starting point for tomorrow)

Source: `clus/src/clustering_examine_bundles.cxx` (see also the algorithm review
`clus/docs/clustering/review_neutrino_isolated_examine_bundles.md`). Two structural
costs that scale badly with cluster/blob count:

1. **O(n²) pairwise cluster loop** (`:113–129`): `for i … for j=i+1 …` over
   `pre_clusters`, building the connectivity graph. Quadratic in the number of
   pre-clusters — on a dense event with many clusters this is the explosion.
2. **Per-cluster blob-graph build** (`:148–162`): `live_clusters[i]->connected_blobs(dv, pcts, graph_name)`
   for every cluster — builds/walks a per-cluster blob graph; cost grows with blobs
   per cluster.

Likely-fruitful directions (to evaluate tomorrow, none done yet):
- **Spatial pre-filter on the O(n²) pair loop** — skip cluster pairs whose bounding
  boxes / drift-time spans cannot connect, so `j` only ranges over near neighbours
  (turn O(n²) into ~O(n·k)). Biggest expected win on the tail events.
- **Profile `connected_blobs`** itself (it recurs in Deghost and Separate too — see
  below) — it is the shared primitive under the top-4 stages, so speeding it helps
  broadly.
- **Cache/reuse** geometry the pairwise loop recomputes per (i,j).

The next tier — **ProtectOverclustering, Deghost, Separate** (together ~37 %) — also
peak on the same events (59789), which strongly suggests a **shared underlying
primitive** (blob connectivity / point-cloud distance queries via `connected_blobs` /
`DynamicPointCloud`) rather than four independent problems. Optimizing that primitive
is plausibly a single change that moves all four. (Note: the crash fix committed today
was in exactly this primitive — `DynamicPointCloud::fill_wrap_points` — so it is
already on the radar.)

## Suggested order for tomorrow

1. **Instrument imaging** with `took` timing (we are currently blind to it).
2. **Re-pull locked numbers** from the finished `f{1,2,3}on` logs.
3. **Attack ExamineBundles** first (56 %): spatial pre-filter on the `:113–129` pair
   loop; profile `connected_blobs`.
4. Check whether the shared connectivity primitive explains the ProtectOverclustering
   / Deghost / Separate tail — one fix may cover all four.
5. Validate any change is **clustering-output-identical** on the 10 hand-scan events
   before trusting it on the tail events.
6. **Follow-up on today's `fill_wrap_points` crash fix** (`DynamicPointCloud.cxx:964`):
   the fix restores the documented `wires_all.size() - 1` clamp (a strict improvement —
   it converts an out-of-bounds read into the long-standing pre-2025-07-28 behavior). But
   the *correctness* of clamping when `point2wind` overshoots is a pre-existing question:
   log `wind - wires_all.size()` where the clamp fires. If it is always ~+1 these are
   genuine edge-of-plane points and clamping to the last wire is right; if it is *large*,
   `point2wind` may have a coordinate/pitch mismatch and the clamped wire is arbitrary on
   exactly the busy events — i.e. non-crashing but subtly wrong. This session verified
   non-crash on all 150 events and byte-identical matched output on the quiet hand-scans
   (where the clamp does not fire), but did NOT separately verify the clamp value on the
   busy tail.
