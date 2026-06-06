# SBND imaging timing profile (lan-reco2) — hotspots & optimization targets

**Purpose.** CPU profile of the SBND **imaging step** (SP frames → blobs/clusters npz,
via the `img/` tiling + charge-solving + geometry-clustering + deghosting chain) on the
10 busiest events of the `input-3files-lan-reco2` reprocess (150 real BNB+cosmic data
events). Identifies *which imaging stage / leaf function* dominates and how it scales
with event busyness, so optimization can be targeted. Companion to
`clus/docs/clustering-timing-profile.md` (which covers the downstream MABC clustering).

**This is a measurement-only study — no source was changed, output is unaffected.**

## Method

- **What is profiled:** the standalone single-event imaging graph
  `sbnd_xin/wct-img-all.jsonnet` (multi-3view, `full_deghost=true` default), run on each
  event's `sp-frames.tar.bz2`. This run produces only `icluster-apa{0,1}-{active,masked}.npz`;
  the downstream MABC clustering happens later in the matching process, so **the whole
  profiled process *is* imaging** — clean attribution with no instrumentation.
- **Engine:** `Pgrapher` (`wct-img-all.jsonnet:87`) = **single-threaded sequential**. Each
  imaging stage is a separate INode dispatched through a virtual `operator()`, so
  `google-pprof --cum` separates the stages natively at the call boundary — this is why
  imaging needs no `took` instrumentation where MABC's single-node sub-stages did.
- **Profiler:** gperftools CPU sampling —
  `LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libprofiler.so.0 CPUPROFILE=… CPUPROFILE_FREQUENCY=250`,
  analysed with `google-pprof --text|--cum build/apps/wire-cell <prof>` (binary + img/clus
  libs are unstripped with debug info). Driver: `/home/xqian/tmp/imgprof/profimg.sh`.
- **Profiled path is the HEAVY chain:** `full_deghost=true` ⇒ **ProjectionDeghosting ×2 +
  ChargeSolving ×3 + InSliceDeghosting ×3** per APA, plus the 4-way active + 3-way masked
  tiling fanout. `run_local.sh`/`run150.sh` use this default, so it is production-equivalent.
- **Single-thread ⇒ CPU ≈ wall:** captured CPU-seconds match wall to <1 % on every event
  (e.g. 141530: 7413 samples / 250 Hz = 29.7 CPU-s vs 29.8 s wall) → imaging is **fully
  CPU-bound**; bz2 decompress (~0.2 %) and npz write are negligible. The hotspot is compute.

### Gotchas (so this reproduces)
- **Pre-compile the jsonnet first.** Profiling `wire-cell <jsonnet>` directly aborts ~70 % of
  the time with a **Go-runtime GC backtrace** — libprofiler's `SIGPROF` colliding with
  gojsonnet's Go runtime during config compilation. Compile once with
  `wcsonnet -o cfg.json -A input=… -S 'anode_indices=[0,1]' -A output_dir=… wct-img-all.jsonnet`
  (unprofiled), then profile `wire-cell -c cfg.json`. With no Go in the profiled process all
  11 runs were clean, single-try. (This is **not** the `clus_all_apa` heap heisenbug of
  `clus/docs/clustering-timing-profile.md` §0 — that lives in the full-file *matching*
  process; single-event imaging is clean natively.)
- **Isolated walls ≈ 0.68× the `timings150.csv` `img_wall_s`** — the CSV was measured under
  `run150` parallel contention (many events at once). Ranking is preserved; absolute CPU
  attribution (the point of this study) is load-independent.
- **Do not merge profiles across events** for leaf attribution — separate processes have
  different ASLR bases and `google-pprof` mis-resolves symbols in a merged set. Use
  per-event profiles (they resolve cleanly).

## The 10 busiest events

Ranked by `timings150.csv` `img_wall_s`; isolated profiling wall ≈ CPU-seconds.

| rank | evt | img_wall_s (csv) | isolated CPU-s | npts | source file |
|---|---|--:|--:|--:|---|
| 1 | 141530 | 43.2 | 29.7 | 608,648 | f1 |
| 2 | 59789  | 26.6 | 17.4 | 157,817 | f1 |
| 3 | 62209  | 22.5 | 14.4 | 452,982 | f1 |
| 4 | 58821  | 20.5 | 12.5 | 448,534 | f1 |
| 5 | 60669  | 20.0 | 12.9 | 336,051 | f1 |
| 6 | 60933  | 19.5 | 12.8 | 370,738 | f1 |
| 7 | 138824 | 19.2 | 13.3 | 328,849 | f1 |
| 8 | 188706 | 17.6 | 14.3 | 321,535 | f3 |
| 9 | 185428 | 17.3 | 11.8 | 378,193 | f3 |
| 10 | 184262 | 16.5 | 11.0 | 244,492 | f3 |
| (10t) | 65025 | 16.5 | 10.9 | 398,173 | f1 |

## Headline: imaging is graph-lifecycle bound, not numeric

Per-event leaf buckets (sample %, gperftools `--text`), representative events:

| evt | allocator+destroy | Rb_tree (setS edge tree) | LASSO Fit | Eigen sparse | MaskSlice thresh |
|---|--:|--:|--:|--:|--:|
| 141530 | 30.0 % | 26.3 % | 5.7 % | 12.2 % | 2.8 % |
| 59789  | 32.0 % | 20.8 % | 7.6 % | 9.2 %  | 4.8 % |
| 60933  | 34.2 % | 26.8 % | 0.8 % | 3.5 %  | 6.2 % |
| 184262 | 33.7 % | 24.5 % | 0.5 % | 3.5 %  | 7.4 % |

**~52–58 % of all CPU is memory management + `std::_Rb_tree` (red-black tree) operations**,
stable across the whole size range. The *actual numeric content* — the LASSO solve
(`LassoModel::Fit`) plus Eigen sparse-matrix ops — is only **~10–15 %**.

**The allocator churn and the Rb_tree churn are one root cause: the cluster-graph type.**
`cluster_graph_t = boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
cluster_node_t>` (`iface/inc/WireCellIface/ICluster.h:167`) uses a **`setS` out-edge container
— a `std::set` (red-black tree) per vertex**. So every `add_edge` is a tree insert and every
graph teardown is a tree erase. A `--focus=_Rb_tree` on 141530 confirms the Rb_tree callers are
*entirely inside boost*: `boost::add_edge`/`push_dispatch` 38 %, the `stored_vertex`/
`vec_adj_list_impl` destructors 30 %, `copy_graph` 15 % — **not** any application-level
`std::map`. The `--cum` view shows the same lifecycle dominating raw CPU: on 141530
`~vec_adj_list_impl` **20.9 %**, `~adjacency_list` 16.6 %, `add_edge` 15.3 %, `copy_graph` 14.6 %.

The cluster graph is **built, copied, and torn down repeatedly**: every `ChargeSolving` pass
unpacks it into per-slice / per-connected-component subgraphs and repacks it (×3),
`ProjectionDeghosting` copies clusters and rebuilds per-cluster projections (×2), and
geom-clustering re-derives blob-blob edges — each `add_edge`/destroy paying the `setS` tree
cost, and `full_deghost=true` multiplying the passes by 2–3. **That churn — not the tiling and
not the LASSO math — is the imaging hotspot.**

## Per-stage breakdown (cumulative %, all 11 events)

`google-pprof --cum` share of each imaging INode `operator()` (sums include each stage's own
allocator/graph children). All instances of a class are aggregated (e.g. ChargeSolving = the 3
solve passes; ProjectionDeghosting = both passes).

```
evt        CPU-s | PrjDg ChSlv InSlD LoGCl BlbGr MskSl BlbCl GlGCl GrdTl
141530     29.7  |  27.2  21.0   9.2   7.2   9.2   3.9   3.7   1.9   0.9
59789      17.4  |  20.8  20.8   9.4   8.3   6.8   6.7   3.4   2.9   0.8
62209      14.4  |  11.9  17.4   9.3   8.7   7.9   8.1   2.3   2.7   1.2
58821      12.5  |   9.1  16.9   9.9   9.2   6.8   9.1   2.2   2.8   1.8
60669      12.9  |  11.9  18.5   9.2   8.6   7.1   8.9   2.6   2.9   1.6
60933      12.8  |  10.3  17.1   9.8   9.2   7.5   9.1   2.3   2.9   1.5
138824     13.3  |  16.3  16.8   9.5   8.6   7.0   8.5   3.0   2.8   1.0
188706     14.3  |  11.5  17.0  10.1   9.4   8.7   8.0   2.3   2.9   1.6
185428     11.8  |  10.3  17.2   9.3   8.8   7.0   9.6   2.3   2.9   1.7
184262     11.0  |  10.3  15.6   8.9   8.5   6.6  10.3   2.4   2.6   1.5
65025      10.9  |  10.6  15.9   9.2   8.3   7.0  10.5   2.3   2.7   1.4
```
PrjDg=ProjectionDeghosting, ChSlv=ChargeSolving, InSlD=InSliceDeghosting,
LoGCl=LocalGeomClustering, BlbGr=BlobGrouping, MskSl=MaskSlices, BlbCl=BlobClustering,
GlGCl=GlobalGeomClustering, GrdTl=GridTiling.

### Scaling — what grows with busyness
- **ProjectionDeghosting is the tail-driving stage** — its share climbs from ~9–12 % on the
  lighter events to **27 % on the worst event (141530)** — but the growth is **graph-lifecycle
  churn, not the O(n²) compare.** On 141530 its 27 % splits as `Projection2D::judge_coverage`
  **5.9 %** (the documented O(n_clusters²)-per-plane pairwise check, `efficiency-concerns.md`
  #2), `Projection2D::get_projection` **5.8 %** (per-cluster Eigen sparse build),
  `connected_components` ~0.8 %, and the remaining **~14 % is cluster-graph copy/destroy**
  (each cluster is copied to its own graph, projected, then discarded). So judge_coverage stays
  ~6 % even on the worst event; what *grows* is the per-cluster graph build/copy/destroy — which
  reinforces the graph-lifecycle headline rather than an O(n²) blow-up.
- **ChargeSolving is the steady #2 (~16–21 %)**, growing modestly with size. `CS::solve`
  (the LASSO) is ~11 % of it; the rest is `unpack`/`repack`/`prune` graph surgery (×3 passes).
- **InSliceDeghosting (~9.4 %), LocalGeomClustering (~8.6 %), BlobGrouping (~7.2 %)** are
  flat-share — they track event size roughly linearly.
- **MaskSlices' share *shrinks* with busyness** (10.5 % → 3.9 %): per-channel-per-slice
  thresholding is ~linear, so it is *diluted* as the O(n²) deghosting grows. Not a target.
- **GridTiling is negligible everywhere (0.8–1.8 %).** The `RayGrid::make_blobs` tiling that
  is intuitively "the imaging" is **not** a hotspot — important negative result.

### The 59789 discriminator (resolved)
59789 is #2 in wall but has the *fewest* points of the ten (158 k vs 141530's 608 k) — ~2.2×
more CPU *per point*. Its stage profile is the same shape as the giant's (deghost+solve
dominated, tiling negligible), confirming imaging cost is driven by **blob/cluster topology
(number of clusters, slices, connected components, cluster pairs)**, not raw point count.
`npts` is the *post-imaging* point-cloud size and does not capture the intermediate blob
fragmentation that feeds ProjectionDeghosting's O(n²) and ChargeSolving's per-component solves.

### Cross-checks
- **CPU ≈ wall** on every event (no I/O confound) — see Method.
- **InSliceDeghosting's TimeKeeper** lines (it is the one already-instrumented stage) sum to
  ~1.82 s on 59789 = 10.4 % of 17.5 s, vs the profile's 9.4 % cum — agree within sampling
  noise.

## Recommendations

Ordered by leverage. All are follow-up proposals; if pursued, each must keep outputs
bit-identical or be jsonnet-togglable default-OFF per the project's change discipline.

**R1 — Cut the cluster-graph build/copy/destroy churn (biggest lever, the ~52–58 % allocator+
Rb_tree).** The `setS`-edged boost `adjacency_list` is repeatedly unpacked into subgraphs,
copied (`CS::prune`/`repack`, per-cluster projection copies), and destroyed across
ProjectionDeghosting ×2, ChargeSolving ×3, and geom-clustering — each `add_edge`/destroy paying
the `std::set` tree cost. Options, roughly increasing risk:
  - *Avoid full-graph copies* in `CS::prune`/`repack`, the per-cluster projection loop, and the
    deghosting filter — operate in place or via `boost::filtered_graph` vertex masks instead of
    materialising and tearing down new graphs.
  - *Reuse scratch graphs/containers* across the 3 ChargeSolving strategy passes and the 2
    ProjectionDeghosting passes instead of fresh-allocating each pass.
  - *Config lever (large saving, but a probable physics non-starter — list, don't chase):*
    `full_deghost=false` drops ProjectionDeghosting entirely (−~10–27 %) and trims the
    solve/deghost passes. **But `full_deghost=true` is where SBND's W-plane dead-channel charge
    recovery happens** (`clus/docs/clustering-timing-profile.md` §0,
    [[project_sbnd_dead_gap_column]] territory), so turning it off likely regresses the very
    feature the local active imaging was built for. It is a clean toggle
    (`wct-img-all.jsonnet:30`), so it is worth *quantifying*, but it is not a free win.

**R2 — Spatial pre-filter on ProjectionDeghosting's `judge_coverage` (the explicitly-O(n²)
part, ~6 %).** `judge_coverage` compares every cluster pair per plane; pre-filter pairs by their
2D (channel × slice) projection bounding boxes and skip the sparse comparison for disjoint pairs
(O(n²)→≈O(n·k); mirrors the ExamineBundles pair-loop prefilter). Scope it honestly: this is the
~6 % named compare, **not** the 27 % stage — most of ProjectionDeghosting is the graph copy/
projection-build that R1 targets. Separately, `reserve()`/reuse the Eigen triplet buffers in
`get_projection` (~6 %, per-cluster), which currently rebuilds a sparse matrix per cluster.

**R3 — Reduce the `setS` edge-tree cost (~21–27 % Rb_tree) — same root cause as R1.** The Rb_tree
is the cluster graph's own `setS` out-edge container (`ICluster.h:167`), not any application
`std::map`, so the lever is structural, not a local container swap: (a) the high-value, lower-risk
path is **fewer graph builds/copies** (R1) — every avoided `add_edge`/destroy is an avoided tree
op; (b) a deeper option is changing the out-edge selector (`setS`→`hashS`/`hash_setS` for O(1)
edge insert, or `vecS` where parallel-edge dedup isn't needed), but `cluster_graph_t` is a core
`iface` type shared by every consumer, so this is a broad, carefully-validated change, not a
surgical one.

**R4 — Lower priority.** InSliceDeghosting's nested 2-view×3-view loops
(`efficiency-concerns.md` #1, stable ~9.4 %) would need grid/spatial indexing — real but not
tail-driving. MaskSlice thresholding and GridTiling are **not** worth touching (small and/or
shrinking share).

## Reproduce

```
# 1. pre-compile config (Go runtime here, unprofiled)
wcsonnet -o cfg.json -A input=<work>/sp-frames.tar.bz2 \
  -S 'anode_indices=[0,1]' -A output_dir=<scratch> \
  sbnd_xin/wct-img-all.jsonnet
# 2. profile (no Go in process)
LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libprofiler.so.0 \
  CPUPROFILE=img.prof CPUPROFILE_FREQUENCY=250 \
  wire-cell -l img.log:debug -L info -c cfg.json
# 3. analyse
google-pprof --cum  --text build/apps/wire-cell img.prof   # per-stage
google-pprof        --text build/apps/wire-cell img.prof   # per-leaf
```
Driver + per-event profiles: `/home/xqian/tmp/imgprof/` (`profimg.sh`, `aggregate.sh`,
`img_<evt>.prof`). Inputs: `sbnd_xin/work/evt<ID>/sp-frames.tar.bz2`.
