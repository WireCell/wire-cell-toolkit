# PDHD / PDVD imaging + clustering optimization log

Running log of the memory-footprint and wall-time optimization work that
follows the baseline measured in [imgclus-resource-profile.md](imgclus-resource-profile.md).

**Ground rule:** every change applied here must produce *functionally
identical* outputs, verified by content-hash A/B comparison on a fixed event
set. Ideas that would change results (even benignly) are **not** applied;
they are collected in the "Result-changing ideas" section at the end for
review.

## A/B methodology

- Harness: `wcp-porting-img/abtest/` (`run_events.sh`, `ab_compare.sh`,
  `hash_archive.py`, `timecmd.py`). Snapshots of all per-event output
  archives + wall/peak-RSS metadata land in `abtest/snap/<label>/`.
- Output archives embed wall-clock timestamps (tar member `mtime=time(0)` in
  `util/inc/WireCellUtil/custard/custard.hpp`; zip likewise), so raw archive
  bytes are never reproducible. Comparison hashes **member payloads**
  (sha256 over name+payload, members sorted by name).
- Peak RSS via `getrusage(RUSAGE_CHILDREN)` (exact high-water mark over the
  full child tree); wall via wrapper timer. Events run 6-way concurrent on a
  64-core host, `OMP_NUM_THREADS=MKL_NUM_THREADS=1`.

### Fixed event set

| Tag | Detector | Run / Evt | Why |
|---|---|---|---|
| hd-typ  | PDHD | 027409 / 0  | typical |
| hd-typ2 | PDHD | 027980 / 3  | second typical |
| hd-busy | PDHD | 028084 / 18 | busy tail (981 s / 7.2 GB in profile) |
| hd-max  | PDHD | 027305 / 0  | extreme (1501 s / 10.1 GB in profile) |
| vd-typ  | PDVD | 039349 / 0  | typical |
| vd-busy | PDVD | 039252 / 5  | PDVD tail |

A 22-event spot-check list spanning every profiled run is in
`abtest/spot_events.txt`; content hashes of the pre-optimization outputs are
stashed in `abtest/snap/preexist/hashes.txt`.

## Baseline (this work, current master add080f5, snapshot `base1`)

6-way concurrent, wall = stage wall-clock (s), RSS = exact peak over the
process tree (MB).

| Event | img wall | img RSS | clus wall | clus RSS |
|---|---|---|---|---|
| hd-typ  027409/0  | 66   | 703   | 28  | 704 |
| hd-typ2 027980/3  | 70   | 691   | 31  | 769 |
| hd-busy 028084/18 | 656  | 7180  | 338 | 2833 |
| hd-max  027305/0  | 1028 | 10113 | 530 | 3284 |
| vd-typ  039349/0  | 140  | 441   | 26  | 745 |
| vd-busy 039252/5  | 327  | 1918  | 145 | 1879 |

## Determinism findings

- `base1` vs `base2` (two back-to-back full reruns, same binary, same
  inputs): **all 138 output archives byte-identical at member level** for
  all 6 events, both stages, including the all-APA clustering outputs. The
  imaging+clustering chain is fully deterministic on PDHD/PDVD, so the
  byte-identity gate applies to every output file with no exclusions.
- Content hashes of outputs that pre-dated this work (stash
  `snap/preexist`) partially disagree with `base1` (notably masked-fork
  imaging archives). Since run-to-run determinism is proven, those
  disagreements reflect input/config state during the original profile
  batch (e.g. PDVD `sp-frames` symlinks re-pointed mid-batch), not
  nondeterminism. The preexist stash is therefore NOT used as a comparison
  baseline; the spot-check baseline is regenerated with current code.

## Optimization entries

### 1. Per-anode sequential imaging (`run_img_evt.sh -P`) — script-level, no C++ change

One `wire-cell` process per anode instead of all anodes in one process.
Per-anode configs are pre-compiled in parallel with `wcsonnet` (also removes
the single large gojsonnet compile from the critical path); the per-anode
processes then run sequentially. Implemented in
`wcp-porting-img/{pdhd,pdvd}/run_img_evt.sh` (wcp-porting-validation repo).

**A/B verdict: PASS** — all imaging archives and all downstream clustering
outputs byte-identical to `base1` on all 6 events (snapshot `peranode`).

| Event | img wall (s) | img RSS (MB) |
|---|---|---|
| hd-typ  | 66 → **50**  | 702 → 657 |
| hd-typ2 | 70 → **53**  | 691 → 663 |
| hd-busy | 656 → **595** | 7180 → 7152 |
| hd-max  | 1028 → **996** | 10113 → 10090 |
| vd-typ  | 140 → **50** | 440 → 396 |
| vd-busy | 327 → **183** | 1918 → 1876 |

Take-aways:
- **Wall improves everywhere** (PDVD typical 2.8x: the 8-anode jsonnet
  compile + geometry init dominated). Clustering unchanged.
- **The RSS tail is NOT a sum-of-anodes effect.** On hd-max the peak barely
  moved: one anode alone (anode0, 809 s of the 996 s) holds the ~10 GB.
  Per-anode wall split on hd-max: 809/135/18/31 s; on hd-busy: 106/12/465/8 s.
  The tail memory therefore lives *inside* a single anode's imaging pipeline
  (cluster graph + per-pass rebuilds), and is attacked by the in-process
  work below, not by process splitting.
- `-P` is opt-in; default behavior of the scripts is unchanged.

### 2. BlobGrouping: per-wire channel cache (imaging)

`img/src/BlobGrouping.cxx` `doit()`: the blob→wire→channel traversal
iterated each wire's full `setS` adjacency set — dominated by the *other
blobs* sharing the wire — to find the wire's channel(s), for every
(blob, wire) pair: cost ∝ Σ deg(wire) per blob ⇒ quadratic in blob count
(65% of total CPU on the worst PDHD anode, see Phase-2 findings). Now the
wire's adjacent channels are cached per wire on first visit, in the same
adjacency iteration order, and replayed afterwards — the visit sequence is
provably identical (wires never gain edges during `doit`; only measure
edges to blobs/channels are added).

**A/B verdict: PASS** (snapshot `opt1`, all archives both stages).

| Event | img wall (s) | clus wall (s) |
|---|---|---|
| hd-max  | 1028 → **439** | 530 → 503 |
| hd-busy | 656 → **247**  | 338 → 318 |
| vd-busy | 327 → **281**  | 145 → 136 |
| hd-typ / hd-typ2 / vd-typ | unchanged | unchanged |

(The clustering deltas in this table come from entry 3, run in the same
A/B round; imaging RSS unchanged as predicted — the live peak is the
projection cache + graph, not BlobGrouping.)

### 3. Grouping::fastgeom — per-(apa,face) conversion-constant memo (clustering)

`clus/src/Facade_Grouping.cxx`: `convert_3Dpoint_time_ch()` did, per call:
a heap `std::vector` for the 3 wire angles, a by-value
`IAnodePlane::faces()` vector copy, and ~10 nested `.at()` map-chain
lookups; `has_closest_point()`/`get_closest_points()` repeated the
wire-angle chains. These run per point per plane inside the good-point
tests under Separate/Deghost/connect passes. New `fastgeom(apa,face)` memo
(pattern follows the existing `m_kd2d_scope_cache`) holds copies of the
constants + the `IAnodeFace` pointer; `convert_time_wire_2Dpoint()` also
converted. Values are identical by construction.

**A/B verdict: PASS** (same `opt1` round as entry 2). Clustering wall:
hd-max 530→503 s, hd-busy 338→318 s, vd-busy 145→136 s (~5-6% on busy
events; the remaining Separate cost is kd-tree geometry and
DynamicPointCloud construction).

### 4. ISlice::activity() by const-ref + ProjectionDeghosting set-copy removal (imaging)

Post-entry-2 re-profile made ProjectionDeghosting the imaging leader
(39.5% of the busy anode); inside it, `std::set` copy-construct/destroy was
~36% (the eight `auto b_cluster = c2b[...]` statements each copied the
cluster's blob-descriptor set) and `Slice::activity()` map copies were
~12% (by-value return, called per slice per cluster in
`Projection2D::get_projection`; GridTiling called it 3x per slice).

- `ISlice::activity()` now returns `const map_t&` (implementers SimpleSlice
  and ImgData::Slice store the map; all callers compile unchanged or were
  bound by const-ref). Call sites that relied on `operator[]`
  default-insert on a local copy (Projection2D, BlobGrouping) use
  find-with-zero-default — exactly equivalent.
- ProjectionDeghosting: `c2b[...]` results bound by const-ref;
  `remove_blobs` takes the tagged set by const-ref.

**A/B verdict: PASS** (snapshot `opt2`). Cumulative vs baseline:
hd-max img 1028→**408 s**, hd-busy img 656→**216 s** (3.0x), vd-busy img
327→275 s; typical events unchanged.

### 5. tcmalloc for the imaging stage (deployment, run scripts)

With the algorithmic fixes in, imaging is allocator-bound (malloc/free
~50% + `setS` Rb_tree ops ~25% of the busy anode, spread over all passes
— the graph build/copy/teardown lifecycle). `LD_PRELOAD`ing
`libtcmalloc_minimal` was A/B'd over all 6 events, both stages
(snapshot `opt2tc`, vs `base1`):

- **Imaging: all archives byte-identical**, large wins:
  hd-max 1028→**308 s**, hd-busy 656→**163 s** (4.0x vs baseline), and
  busy-event RSS down (10113→8862 MB, 7180→6857 MB).
- **Clustering: FAILED the gate on hd-typ** → tcmalloc stays OFF for
  clustering (see finding below).

Adopted in `run_img_evt.sh` (both detectors, wcp-porting-validation):
imaging runs under tcmalloc by default; `WCT_TCMALLOC=off` reverts.

### 6. kd existence query early-exit + DynamicPointCloud real move (clustering)

- `NFKDVec::Tree::exists_within(r2, q)`: terminates the k-d search at the
  first point strictly inside the squared radius (nanoflann `addPoint`
  returning false). Used by `Grouping::has_closest_point` — the hottest
  clustering leaf (21.6% cumulative on hd-max; previously knn(1) resolved
  the true nearest neighbor only to compare against the radius). Boolean-
  identical by construction.
- `DynamicPointCloud::add_points` "moved" its input through
  `make_move_iterator` over a **const** range — a silent deep copy of
  every DPCPoint (5 nested vectors each; DPCPoint ctor was 7.2% flat).
  Added a genuine rvalue overload (all hot callers pass temporaries) and
  re-pointed the kd-indexing tail to read from `m_points`, which also
  removes a latent read-after-move had the move ever become real.

**A/B verdict: PASS** (snapshot `opt3`, includes entry-5 tcmalloc imaging).
Clustering wall vs baseline: hd-max 530→**454 s**, hd-busy 338→**298 s**,
vd-busy 145→134 s. Imaging vs baseline (with tcmalloc): hd-max
1028→**304 s**, hd-busy 656→**161 s**, vd-busy 327→247 s.

## Closing validation

- **22-event spot check** (`spot_events.txt`, spans every profiled run,
  both detectors): optimized code + tcmalloc imaging vs pre-optimization
  baseline — **all 706 output archives byte-identical** (snapshots
  `spotbase` vs `spotopt`).
- **Determinism re-run** (`opt3` vs `opt3b`, full 6-event set twice with
  final code): all archives identical (hd-max clustering needed a re-run
  in the second pass — see crash finding below).

### FINDING: intermittent clustering crash is the gojsonnet Go runtime

During the determinism re-run, hd-max (027305/0) clustering aborted
(rc=134) with `fatal error: traceback did not unwind completely` raised
by the **embedded Go runtime** (gojsonnet's GC threads, which keep
running after config parsing). The identical binary+input passed
byte-identically an hour earlier, and the pre-optimization 283-event
batch had a similar 027305 segfault — this is the long-standing
intermittent "clus heisenbug", now with a concrete signature: it is the
Go GC crashing, not (necessarily) the clustering C++.

**Mitigation to try:** pre-compile the clustering config with `wcsonnet`
and feed `wire-cell` pure JSON (as `-P` imaging mode already does). If
the Go runtime is only initialized when jsonnet is actually parsed, this
removes the Go GC from the production process entirely; it also takes
the jsonnet compile off the clustering critical path. Needs a check that
`wire-cell -c cfg.json` does not still spin up Go threads (inspect
/proc/<pid>/task count or the crash disappearing over many runs).

**RESOLVED — see round-2 entry 7** (pure-JSON config + `GOGC=off`).

### FINDING: clustering has a pointer-order-dependent merge decision

Under tcmalloc, hd-typ (027409/0) per-face clustering of apa2-face0
produced **119 clusters instead of 120** — same total points, same
charges, but a different partition (one merge decision flipped), which
then cascades through group02 / all-apa outputs. So the chain is NOT
allocator-independent: some container keyed/ordered by heap pointers
(e.g. `Cluster*`/`Blob*` sets or maps, or a pointer tie-break in a sort)
steers a merge. Under glibc the outcome is reproducible run-to-run only
because the allocation sequence happens to be identical. This is very
likely the same family as the long-standing SBND `clus_all_apa`
nondeterminism. Worth hunting and fixing for robustness (and it would
unlock tcmalloc for clustering too): instrument with
`LD_PRELOAD=libtcmalloc` A/B on a small event and bisect by stage
(per-face → per-APA → group) to find the first pass whose output flips.

**RESOLVED — see round-2 entry 8** (`get_closest_blob` pointer-keyed map).

## Round 2 (2026-06-11): follow-ups

Continuation of the campaign on the three follow-up items picked from the
findings/ideas above: the Go-runtime crash mitigation, the pointer-order
merge hunt, and the ProjectionDeghosting projection-cache bound. Same A/B
methodology and event set. One deliberate exception to the byte-identity
ground rule is entry 8 (a nondeterminism *fix* necessarily changes the
arbitrary tie-break it repairs); its output impact is quantified there.

### 7. Clustering config pre-compiled to JSON + `GOGC=off` (crash mitigation, script-level)

`run_clus_evt.sh` (both detectors, wcp-porting-validation) now compiles
`wct-clustering.jsonnet` with `wcsonnet` and feeds `wire-cell` the pure
JSON, running the process with `GOGC=off` (the imaging `-P` path got the
same `GOGC=off` treatment).

What was learned while verifying the original mitigation hypothesis:

- `libgojsonnet.so` is **hard-linked** into `wire-cell` (via
  WireCellUtil), and merely loading it starts the Go runtime: 5 runtime
  threads appear at `dlopen` time with no jsonnet evaluated. So a pure
  JSON config does *not* remove the Go threads.
- It does remove the Go *heap*: `Persist` only enters the jsonnet VM for
  `.jsonnet` files, and `GODEBUG=gctrace=1` shows the old mode runs 35 GC
  cycles during config evaluation inside the production process (64-P Go
  scheduler, 56 threads total), with the 2-minute *forced* periodic GC
  continuing for the lifetime of long jobs — which is exactly why only
  long busy events ever crashed.
- With pure JSON + `GOGC=off`, the production wire-cell process logs
  **zero GC cycles** end-to-end and runs 7 threads instead of 56.
  `GOGC=off` disables the Go collector *including* the periodic forced
  GC, so the crash vector (GC traceback walking) is removed
  deterministically rather than made statistically rarer. The jsonnet
  evaluation (and its Go GC) lives in the short-lived `wcsonnet` process.

**A/B verdict: PASS** — full 6-event set, all clustering archives
byte-identical vs `opt3` (snapshot `gogc1`); wall/RSS unchanged. A
repetition run of hd-max clustering exercises the formerly crash-prone
path (see closing notes of this round).

### 8. Pointer-order-dependent merge: found and fixed → tcmalloc ON for clustering

The hunt (hd-typ 027409/0, glibc vs `LD_PRELOAD=libtcmalloc`, per-pass
TRACE counts + temporary per-decision instrumentation in
`clustering_connect1`):

- Per-pass cluster counts bisect the flip to **`ClusteringConnect1`**
  (apa2-face0: identical counts through ExtendLoop, then 120 vs 119).
- Per-decision dumps showed identical merge-candidate inputs but a
  flipped angle gate: the *stored* `dir1`/`dir2` of clusters differed
  between allocators — two mirror hough bins (±0.5° off the drift axis)
  swapping winners.
- Full-precision printing then exposed the root: `get_two_extreme_points()`
  returned coordinates differing in the **last FP bit**. Its
  `calc_ave_pos()` (and `calc_dir()`) iterate
  `Cluster::get_closest_blob()` and accumulate `center_pos()*q` —
  and that map was `std::map<const Blob*, geo_point_t>`, i.e. keyed by
  **raw heap pointer**, so the floating-point summation order followed
  the allocator's address assignment. The last-bit difference shifts a
  near-boundary point into the adjacent hough bin, the bin argmax flips,
  one connect1 merge flips, and the difference cascades to the
  group/all-APA outputs.

**Fix:** key `const_blob_point_map_t` with the existing content-based
`BlobLess` comparator (`clus/inc/WireCellClus/Facade_Cluster.h`). All
three consumers (`calc_ave_pos` x2, `calc_dir`) now sum in
content-deterministic order. `blob_less` compares wpid, npoints, charge,
slice and wire ranges before its documented pointer fallback, so
distinct-blob ties are practically impossible.

**Verification (snapshots `pofix` glibc / `pofixtc` tcmalloc):**

- **Determinism gate: PASS** — with the fix, glibc and tcmalloc produce
  byte-identical archives for the full 6-event set (114 archives).
- **This is the one deliberate not-byte-identical change of the
  campaign**: repairing the tie-break necessarily re-resolves it.
  Vs the old baseline, 17/114 clustering archives change, confined to 3
  PDHD events (hd-max 9, hd-typ 4, hd-typ2 4 — per-face/per-APA zips on
  the flipped faces plus their group/all-APA descendants); hd-busy and
  both PDVD events are byte-identical to the old baseline. The changes
  are single merge-decision flips of the same kind tcmalloc used to
  induce — arbitrary tie-break resolutions, now stable.

**tcmalloc adopted for clustering** (`run_clus_evt.sh`, both detectors,
`WCT_TCMALLOC=off` reverts), measured on the gate run (6-way load):
hd-max 458→**372 s**, hd-busy 303→**238 s**, vd-busy 137→**112 s**,
typical events ~20% faster; RSS slightly lower on PDHD (vd-busy RSS
rose 1.9→2.4 GB — tcmalloc trades fragmentation for cache retention;
acceptable against the 3.3 GB clustering tail).

This very likely also resolves (or at least stabilizes) the
long-standing SBND `clus_all_apa` run-to-run nondeterminism family —
worth a dedicated SBND A/B before relying on it there.

### 9. ProjectionDeghosting: evict projection matrices after last use (imaging memory)

Implements option (b) from the memory ideas: the 3.7 GB live peak on the
worst anode was the per-cluster Eigen sparse projection cache
(`id2lproj`) holding every cluster's 3 layer matrices for the whole
pass. The exact liveness rule exploited
(`img/src/ProjectionDeghosting.cxx`):

- Every cache entry is created when its cluster is the *current* vertex
  of the main loop (later accesses at the judge sites only ever touch
  clusters that are still **keys** of a per-layer `clus_2D_3D_map`, and
  keys are only ever added when current).
- Everything after the secondary judge loop (`m_saved_flag` summing, the
  final per-cluster cut, `remove_blobs`) reads only the scalar metadata,
  which survives eviction.

So: clear `m_layer_proj` (keep the metadata struct) as soon as a cluster
has been processed and is not a key of any layer map — ghosts/duplicates
die during the main loop, erased map entries die at erasure, and the
whole cache is cleared after the secondary loop, *before* `remove_blobs`
allocates the output graph it used to coexist with.

**A/B verdict: PASS** (snapshot `imgr2`, all 64 imaging archives
byte-identical vs `opt3`; downstream clustering unaffected by
construction). Peak RSS: hd-max 8862→**6078 MB**, hd-busy
6855→**3339 MB** (−51%), vd-busy 1768→1569 MB; wall unchanged within
noise.

### 10. remove_blobs / CS::prune: old2new remap unordered_map → vector (imaging)

The graph-rebuild remap in `ProjectionDeghosting::remove_blobs` and
`CS::prune` keyed an `unordered_map` by vecS vertex descriptors — which
are dense indices. Replaced with a flat `std::vector` (sentinel =
`max`), removing a hash insert+find per vertex/edge of every rebuilt
graph. Same insertion order, byte-identical by construction; gated in
the same `imgr2` snapshot as entry 9. Wall effect small (a few % of the
rebuild sites), bundled here as the safest of the graph-rebuild
reductions; the `copy_graph` sites remain on the target list.

### 11. NFKDVec::knn1 — allocation-free single-NN query (clustering)

`get_closest_point_blob` (12.5% of hd-max clustering, the inner call of
`Find_Closest_Points`' alternating-projection loops) heap-allocated a
one-element result vector per query via `knn(1)` and resolved the scoped
view twice (once for the query, again inside `point3d`). Added
`NFKDVec::Tree::knn1(query, index, metric)` (stack scalars, identical
nanoflann search) and re-pointed the four single-NN `Cluster` wrappers
(`get_closest_point_blob`, `get_closest_wcpoint`,
`get_closest_point_index`, `get_closest_dis`) at it, resolving the
scoped view once.

**A/B verdict: PASS** (snapshot `knn1` vs `pofixtc`, all archives
byte-identical). Wall effect is noise-level on the gate run (hd-max
372→366 s) — with tcmalloc already adopted the malloc share this
removes had become cheap; kept as the right primitive for these hot
paths.

## Phase-2 profiling findings (PDHD/PDVD-specific)

CPU profile of the pathological anode (hd-busy 028084/18 anode2, 465 s solo;
gperftools, 154835 samples):

- **`BlobGrouping::doit` is 77.8% of the whole imaging job.** 65% of ALL CPU
  is `std::_Rb_tree_increment`, and 99.3% of those samples come from the
  blob→wire→channel traversal in `img/src/BlobGrouping.cxx`: for every
  (blob, wire) pair it iterated the wire's full `setS` adjacency set — which
  contains every *other blob* sharing that wire — to find the wire's single
  channel. Cost ∝ Σ_(b,w) deg(w) → quadratic in blob count. This, not
  ProjectionDeghosting (9.5% here), is the PDHD/PDVD busy-anode hotspot —
  a different balance than the SBND profile in
  `img/docs/imaging-timing-profile.md`.
- Remaining imaging stages on this anode: ProjectionDeghosting 9.5%,
  BlobClustering 2.4%, ChargeSolving 2.2%, InSliceDeghosting 1.2%,
  LocalGeomClustering 1.0%, GridTiling 0.3%.
- **Imaging RSS is an early plateau, not a leak**: the profiled anode jumps
  to 4.1 GB at ~30 s, 7.16 GB at ~60 s (slicing/tiling/graph build), then
  stays flat to the end.
- **Heap profile (tcmalloc) of the same anode**: peak *live* heap is
  5.9 GB; cumulative allocation churn is ~140 GB. So roughly 1.2 GB of the
  glibc 7.16 GB plateau is allocator retention/fragmentation of the churn
  (glibc never returns it; live drops to ~1 GB later in the job while RSS
  stays at the high-water mark). The live peak splits as:
  - **3.7 GB (63%): `Eigen::SparseMatrix<float>::resize`** — the
    per-shadow-cluster 2-D projection matrices that `ProjectionDeghosting`
    memoizes in `id2lproj` for *every* cluster for the duration of each
    pass (`img/src/ProjectionDeghosting.cxx:186`). Each cluster holds 3
    layer matrices dimensioned (nchan × nslice); the dense outer-index
    array (~nslice × 4 B each) dominates when busy events have O(10⁴⁻⁵)
    shadow clusters.
  - ~1.3 GB: boost `setS` adjacency edge trees + `stored_vertex` vectors
    of the cluster graph(s).

Clustering profile (hd-busy full clustering, 86296 samples), cumulative:

- MABC pass split (from the always-on `perf` logs, busy events):
  **Separate 32-35%**, Deghost 9-18%, ExamineBundles 10-15%, then
  Close/Connect1/ProtectOverclustering ~5% each. Typical events: Deghost
  ~20%, Connect1 ~11%, ExamineBundles ~9%.
- Inside Separate: `connect_graph_ctpc` 22.5% → `Grouping::is_good_point`
  20.4% → `has_closest_point` 14.4% (kd2d knn) + `get_closest_dead_chs`
  9.4% + `convert_3Dpoint_time_ch` 8.3%.
- `convert_3Dpoint_time_ch` did, per call: a heap `std::vector` for 3
  angles, a by-value `IAnodePlane::faces()` vector copy, and ~10 nested
  `std::map`/`unordered_map` `.at()` chain lookups — visible in the profile
  as `_Vector_impl_data` 5.3%, `_M_lower_bound` 4.8%, malloc/free ~20%.
- nanoflann `searchLevel`+`evalMetric` ~19% (genuine k-d geometry,
  irreducible), `DynamicPointCloud` point construction ~12%.

## Next runtime targets (profiled, not yet attempted)

From the post-fix profiles (imaging: hd-busy anode2 now ~95 s CPU;
clustering: hd-max ~510 s):

- **Imaging** (now allocator/graph-lifecycle bound, spread across passes:
  PD 25%, ChargeSolving 22%, BlobClustering 9%, InSliceDeghosting 8%,
  LocalGeomClustering 8%): remaining structural lever is reducing
  cluster-graph rebuilds (`remove_blobs` full rebuild, `CS::prune`
  copy-then-filter ×3 passes, the `copy_graph` sites). Each is a few %;
  the `old2new` unordered_map→vector remap in `remove_blobs`/`CS::prune`
  is the safest first step.
- **Clustering** (hd-max): `Find_Closest_Points` 14.3%,
  `Cluster::sv3d`/scoped-view lookups under ExamineBundles 13.2%,
  `get_closest_point_blob` 12.5%, vhough 10.6%, `DynamicPointCloud`
  construction ~9% (DPCPoint is AoS with 5 nested vectors per point — an
  SoA restructure is the big-ticket item, byte-identity achievable but a
  larger surgery), nanoflann searchLevel ~20% (genuine geometry).
- **Clustering pointer-order dependence** (see FINDING below): fixing it
  is both a robustness win and unlocks tcmalloc for clustering (~15-20%
  more, judging by the opt2tc walls).

## Ideas not yet applied — for review

### Memory (imaging tail), in suggested order of attack

1. **Bound or restructure the `ProjectionDeghosting` projection cache**
   (3.7 GB live on the worst anode). Options, in increasing risk:
   (a) shrink each matrix's dense outer dimension by storing projections
   with a global slice-offset (helps only if activity does not span the
   readout — busy events likely do); (b) evict cache entries once a
   cluster's last cs_graph neighbor has been visited (exact, byte-identical,
   needs a pre-pass to count uses; recompute-on-miss is also byte-identical
   since projections are pure functions of the graph); (c) replace the
   Eigen sparse layer matrices with hash/CSR keyed only on occupied
   (channel, slice) — touches `judge_coverage`, byte-identity risk, last
   resort.
   **(b) DONE in round-2 entry 9** (map-key-membership liveness instead of
   a use-count pre-pass; hd-busy imaging RSS −51%). (a)/(c) remain
   unattempted and likely unnecessary now.
2. **tcmalloc deployment** (`LD_PRELOAD` in run scripts): with ~140 GB of
   churn, glibc retains the high-water mark forever while live memory drops
   to ~1 GB; tcmalloc returns freed pages. Helps the *sustained* footprint
   under 16-way concurrency (not VmHWM of the early peak). Must pass the
   full A/B gate first (allocator changes can flip pointer-keyed iteration
   order if any exists).
   **DONE for both stages** (imaging in round-1 entry 5; clustering in
   round-2 entry 8 after the pointer-order fix).
3. `malloc_trim(0)` after each `ClusterFileSink` write (glibc-only, helps
   between anodes in all-anode mode; mostly superseded by `-P`).

### Result-changing ideas (NOT applied)

- Reduce stepped sampling density of blob point clouds (clustering input
  size) — changes point clouds, needs physics review.
- `full_deghost=false` or `pipe_type="single"` imaging variants — large
  CPU/memory savings but changes blob content.
- Raising the slicing threshold (`nthreshold`) on busy/incomplete-readout
  events — directly reduces blob combinatorics; changes imaging output.
- `setS`→`hash_setS` edge container for `cluster_graph_t` — rejected
  previously (changes physics output, iterator invalidation).
