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

**R2 — Cut `judge_coverage`'s per-call sparse-matrix allocations (6.1 % → 1.1 %, SHIPPED).**
*Re-scoped after measurement — see "R2/R3 investigation" below.* The original idea (bbox pre-filter
+ triplet `reserve()`) does **not** pan out: `judge_coverage` is already gated by `cs_graph` edges
(not an all-pairs O(n²) compare), bbox-disjoint pairs are only 0.16–2.9 % of calls, and
`get_projection`'s cost is `setFromTriplets` (inherent), which `reserve()` cannot touch. The real
lever is that each `judge_coverage` call allocated **three** sparse matrices — `mask(ref)` +
`mask(tar)` + `(ref_mask − tar_mask)` (the subtraction alone is 54 % of the function) — and ~99 % of
calls hit that path. Replaced with a single union-pass over the two matrices' non-zeros (tracking
ref_nz/tar_nz/has_pos/has_neg with early-exit) — no allocations. Bit-identical and shipped; see "R2/R3
investigation" for the equivalence proof and result.

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

**R5 — `BlobShadow::shadow` (NEW, surfaced during the R2 investigation; ~6.3 % → 4.3 % SHIPPED).**
On 141530 `Aux::BlobShadow::shadow` (`aux/src/BlobShadow.cxx`) is the **single largest
ProjectionDeghosting sub-cost — larger than `judge_coverage` (6.1 %) and `get_projection`
(5.8 %)** — yet it was absent from R1–R4. It builds the per-wire blob-overlap graph that gates the
whole deghosting pair loop. Investigated and a bit-identical fix shipped — see "R5 investigation &
fix" below.

## R1 investigation & Tier-A implementation

**Where the copies come from.** With `full_deghost=true` the per-APA chain runs ~18 full
`cluster_graph_t` (re)builds + 6 per-slice unpack restructurings. The single largest is
**`ChargeSolving::repack`** (`CSGraph.cxx:438-490`), which rebuilds the *entire* graph (re-adds
every surviving node + edge, each a `setS` Rb_tree insert) just to overwrite blob values and
drop pruned b/m — and it runs **6×** (2 ChargeSolving per solving block × 3 blocks). The code
comment concedes it: *"removal of vertices… is dreadfully slow due to using vecS… we use a
somewhat verbose construction instead of copy+remove."* `CS::unpack` additionally restructures
into per-slice b-m subgraphs (inherent to the per-slice LASSO — survives any dataflow change).

**Why the big levers are hard.**
- *Selector change (`setS`→`vecS`/`hashS`, Tier B)* would cut both Rb_tree and most per-edge
  allocs, but `ClusterArrays.cxx:482/511` writes the node/edge npz arrays in **unsorted
  `boost::edges()` order**, so the change reorders the output — **not bit-identical** without a
  canonical output sort + an order-dependence audit. Core `iface` type shared by clus/img/aux.
- *Mutable/move dataflow (Tier C)* addresses the ~14 full copies but **not** `unpack`'s
  restructuring nor `ProjectionDeghosting`'s projections, so it recovers a large fraction, not
  "most," of the 55 %. Contract change to every `IClusterFilter`.

**Validation gotcha — imaging npz are NOT byte-deterministic.** Imaging builds graphs/arrays by
iterating `std::unordered_map` keyed by pointers/descriptors, so the npz **row order varies
run-to-run** (heap-address dependent) even though the physics is identical — two unmodified-binary
runs give different md5. So bit-identical here means **canonical-content-identical**: compare each
npz up to row permutation (node rows as a multiset minus the descriptor column; edges remapped
through their endpoints' content-signatures). Tool: `/home/xqian/tmp/imgprof/cmpnpz.py`, validated
to call two unmodified runs IDENTICAL.

**Tier A — implemented, committed, output bit-identical, ~3 % CPU.** Three safe, surgical edits:
- `InSliceDeghosting` round-1: the two back-to-back `copy_graph` passes (remove bad blobs →
  `cg_old_bb`, then strip b-b edges → `cg_new_bb`) collapse into **one** combined vertex+edge
  filtered copy (`InSliceDeghosting.cxx`).
- `Projection2D::get_geom_clusters`: the materialized blob subgraph (`cg_blob` via
  add_vertex/add_edge) → a lazy `filtered_graph` + `connected_components` (mirrors the existing
  `LocalGeomClustering:60-64` pattern); no graph allocation.
- Guard eager debug: `LocalGeomClustering`'s unconditional `connected_components` (debug-only,
  the bug-class `efficiency-concerns.md` #12 already fixed in *GlobalGeomClustering*) and the
  whole-graph `dumps()` calls (evaluated as `log->debug` arguments even at `-L info`, since
  `img/` uses `log->debug(...)` not the short-circuiting `SPDLOG_LOGGER_DEBUG` macro).

Measured (single-event imaging, CPU-seconds = samples/250 Hz):

| evt | baseline | Tier A | Δ |
|---|--:|--:|--:|
| 141530 | 29.7 | 28.8 | −2.9 % |
| 59789  | 17.4 | 16.8 | −3.8 % |
| 60933  | 12.8 | 12.5 | −2.6 % |

Output **canonically identical on all 44 npz** — every one of the 11 profiled events ×
{apa0,apa1} × {active,masked} — via `cmpnpz.py` (old binary vs final modified binary; the
deghosting topology varies event-to-event, so the 11-event breadth is the real correctness check,
not just the 608k giant). This is the safe ceiling for R1; the larger wins (Tier B/C) remain open
and are the user's call.

## R2/R3 investigation

**Where ProjectionDeghosting's time actually goes (141530, `--cum`).** The 27.8 % stage splits
into: `BlobShadow::shadow` 6.3 %, `judge_coverage` 6.1 %, `get_projection` 5.8 %,
`ClusterShadow::shadow` 1.1 %, `remove_blobs` 0.6 %, `judge_coverage_alt` 0.4 % — the rest is the
per-cluster graph copy/destroy R1 targets. So R2's `judge_coverage` is **not** the biggest
sub-cost (R5's `BlobShadow::shadow` is).

**R2 — `judge_coverage` (line-level, 441 samples on 141530).** Costs concentrate in three sparse
allocations: `mask(ref)` line 303 = 87, `mask(tar)` line 304 = 44, and `ref_mask − tar_mask`
line 326 = **239 (54 %)**; the `loop_exist` scans are 56+7+7. The pair loop is **gated by
`boost::edge_range(cs_cluster, clust3D, cs_graph)`** (Projection2D.cxx caller), so only
shadow-overlapping cluster pairs are ever judged — it is *not* an all-pairs O(n²) compare.

Outcome distribution, measured with throwaway counters (since reverted):

| evt | judge_coverage calls | both-non-zero (line-326 path) | bbox-disjoint |
|---|--:|--:|--:|
| 141530 | 24,157 | 24,017 (99.4 %) | 38 (0.16 %) |
| 59789  | 14,392 | 14,383 (99.9 %) | 417 (2.9 %) |

Conclusions:
- **bbox pre-filter is marginal, not ineffective-by-construction.** A `'w'` shadow edge shares a
  *wire* (channel axis) regardless of *slice*, so disjoint (channel × slice) live-bboxes are
  possible and deterministically return `OTHER` — a valid bit-identical short-circuit. But it fires
  on only 0.16–2.9 % of calls, so caching a bbox on `Projection2D` is not worth it.
- **triplet `reserve()` is worthless.** `get_projection`'s 419 samples are `setFromTriplets`
  line 215 = **220 (52 %)**, the 3 empty `Projection2D(nchan,nslice)` ctors = 66, and
  `SimpleSlice::activity()`'s by-value map copy line 150 = 64. `reserve()` touches none of these;
  the triplet `push_back`s do not even appear.
- **The real lever** (SHIPPED): replace `mask(ref)+mask(tar)+(ref_mask−tar_mask)` with a single
  per-column merge over the two matrices' (row-sorted) non-zeros computing the same four booleans
  (`ref_nz`, `tar_nz`, `has_pos`, `has_neg`) with early-exit, then the identical branch logic.
  **Equivalence:** `mask` sets every *stored* cell to exactly `0`/`1` and keeps the sparsity, so
  `ref_mask − tar_mask ∈ {−1,0,1}` over the union of stored cells; `diff=+1 ⟺ (ref live ∧ tar
  not-live)` and `diff=−1 ⟺ (tar live ∧ ref not-live)`, where "live" = stored ∧ `value > −uncer_cut`
  and "not-live" subsumes not-stored — exactly what `ref_m_tar_pos`/`ref_m_tar_neg` test, and the
  ±0.01 thresholds are exact on integer diffs. The early-exit to `OTHER` is sound because the four
  flags are monotone. Removed the now-orphaned `loop_exist` helper (kept `mask`, still used by
  `judge_coverage_alt`). `judge_coverage` has a **single caller** (`ProjectionDeghosting.cxx:264`),
  so the imaging validation covers it fully. **Measured 6.1 % → 1.1 %** on 141530 (~83 % of the
  function removed; ProjectionDeghosting 27.8 % → 22.3 % with R5+this); output canonically identical
  on all 44 npz via `cmpnpz.py`.

**R3 — confirmed scoping.** (a) "fewer graph builds/copies" = R1 (Tier A shipped).
(b) Out-edge selector `setS`→`hashS` (O(1) edge insert, preserves `setS` parallel-edge dedup
semantics — *not* `vecS`, which would admit duplicate edges) is the structural lever for the
~21–27 % Rb_tree, but `cluster_graph_t` (`ICluster.h:167`) is a **core `iface` type shared by every
clus/img/aux consumer**. **Spike attempted and rejected — see "R3b spike" below.**

## R3b spike — `setS`→`hash_setS` selector change (ATTEMPTED, REJECTED)

Ran the one-line change as a throwaway spike to settle the effort/payoff question empirically.
Result: **not a drop-in, and it changes physics output.** Three concrete findings, in order hit:

1. **It cascades past one line.** `boost::hashS` does not exist — the real selector is
   `boost::hash_setS`. And `cluster_indexed_graph_t = IndexedGraph<cluster_node_t>`
   (`util/inc/WireCellUtil/IndexedGraph.h:37`) hard-codes its *own* `setS` graph and feeds
   `.graph()` straight into `SimpleCluster(cluster_graph_t)` at `BlobClustering.cxx:123`, so
   `IndexedGraph` must flip too. (It is the only non-test instantiation, so for the spike both
   were flipped.) ~11 img/aux/sio files mix graph mutation with adjacency walking — the audit
   surface for the change is broad.

2. **It segfaults — the code depends on `setS` iterator stability.**
   `ClusterArrays::bodge_channel_slice` (`aux/src/ClusterArrays.cxx:131`) walks
   `adjacent_vertices(mvtx)` while calling `add_edge(active_cvtx, mvtx)` **inside** the loop,
   mutating `mvtx`'s out-edge set mid-iteration. A node-based `setS` (Rb_tree) keeps the live
   iterator valid across insertion; `hash_setS` (boost `unordered_set`) rehashes and invalidates
   it → SIGSEGV on `++`. A local collect-then-add fix (output-equivalent: `add_edge` dedups, so
   deferring is edge-set identical) cleared this one site, but every mutate-while-iterating site
   would need the same treatment, each a latent UB / crash.

3. **Even patched, the output drifts.** After the fix, 141530 imaging *completed* (the only crash
   site for that event), but `cmpnpz.py` vs the `setS` baseline split:
   **masked APAs IDENTICAL, active APAs DIFFERENT.** The masked-identical result is a control: it
   proves the writer fix is edge-set-correct *and* that `cmpnpz` canonicalization factors ordering
   out across selectors — so the active difference is **genuine content drift**, the hash edge
   iteration order perturbing the active-imaging deghost/solve chain (FP-accumulation order /
   tie-breaks / component labeling). A compile-time typedef cannot be runtime-toggled, so a change
   that alters physics output is **blocked by the default-OFF togglability rule**.

**Timing is moot.** A single uncontrolled run hinted ~24.5 s vs `setS` ~29 s, but a result that
computes *different* physics is not shippable, and the predicted small-degree-out-edge regime (the
same one where the BlobShadow per-pair hash cache lost to a tiny tree, see R5) makes even a clean
win unlikely. **Net: `setS`→`hash_setS` converts an uncertain payoff into a high, latent-UB cost
for a likely-modest and output-changing payoff — not pursued.** All spike edits reverted; tree is
back to the original `setS`.

## R5 investigation & fix (BlobShadow::shadow)

**Where the time goes (141530, 457 samples).** `BlobShadow::shadow` builds the per-wire
blob-overlap graph by looping, per slice, over every *(shared-leaf, blob-pair)* and ensuring one
edge per layer between the pair (extending its `[beg,end]` index range). The hot line was
`existing_layer_edges(bsgraph, v1, v2)` (line 185, **212 samples = 46 %**) + `add_edge` (line 197,
83). Under `--focus=BlobShadow::shadow`, that 46 % splits into the `edge_range` multiset scan
(`multiset::equal_range` **27.4 %**) and a **fresh per-call `unordered_map`** (`insert` 16.6 % +
rehash `_M_need_rehash` 11.2 % / `_M_next_bkt` 6.6 % + `find` 7.2 %): `existing_layer_edges` built a
map of *all* layers' edges between the pair on **every call**, then read back a **single** layer
(`wpid.layer()`, which is constant per leaf) — the map grew from empty each time, forcing rehashes,
and was discarded immediately.

**Fix (shipped, bit-identical).** Replace `existing_layer_edges` (returns a layer→edge map) with
`existing_layer_edge(..., layer)` — a direct `edge_range` scan returning the single matching-layer
edge (first match in `edge_range` order, exactly the previous insert-then-`find` semantics). No
per-call map, no rehash. The `add_edge` order, edge set, and `[beg,end]` values are unchanged, so
the BlobShadow graph — and everything downstream — is identical. Removed the now-orphaned
`bs_layer_edge_t` typedef; left the pre-existing-unused `pair_hash` struct in place (it hints at the
*next* lever — a persistent per-pair edge cache to also kill the 27 % `edge_range` scan).

| evt | BlobShadow share (HEAD → fix) | samples cut | total CPU-s (HEAD → fix) |
|---|---|--:|---|
| 141530 | 6.3 % → **4.3 %** | −147 (−33 % of fn) | 28.8 → 28.1 |
| 59789  | 3.0 % → **2.0 %** | −43 (−35 % of fn)  | 16.8 → 16.5 |
| 60933  | 1.2 % → **1.0 %** | −7                 | 12.5 → 12.2 |

BlobShadow's share is event-dependent (biggest on the busiest, most-overlapping events). The
within-profile drop (~⅓ of the function removed) is the robust signal and matches the predicted
map-removal. Output **canonically identical on all 44 npz** (11 events × {apa0,apa1} ×
{active,masked}) via `cmpnpz.py`.

**Per-pair edge cache — TRIED, rejected (it is slower).** The remaining ~27 % of the function is
the `edge_range` multiset scan, and the dormant `pair_hash` struct hints the author once meant to
cache it (pair → layer → edge), turning each scan into a hash lookup. Implemented and measured
(output still bit-identical, 44/44): it is **substantially slower**, not faster. With the stock
`pair_hash` (`h1 ^ h2`) it was catastrophic — symmetric XOR collides for `(min,max)` index pairs
(e.g. `(0,3)` and `(1,2)` both hash to 3), degenerating the map into linear bucket scans (141530
BlobShadow 4.3 % → 63 %). Even with a proper `hash_combine`, the per-slice `unordered_map` +
per-pair `vector` overhead beats the savings: 141530 BlobShadow **4.3 % → 10.9 %** (total 7014 →
8239 samples), 59789 **2.0 % → 4.0 %**. Boost's `edge_range` is a sorted O(log·deg) multiset
lookup over the ≤3 edges of a pair — already cheap — so a hash cache cannot win here. That is why
the stub stayed dormant; reverted to the direct scan. **No further lever inside `BlobShadow`.**

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
