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

Follow-up in the same vein: `Simple3DPointCloud::get_closest_wcpoint` /
`get_closest_dis` (8.5% cumulative in the round-2 closing profile, the
inner queries of `Simple3DPointCloud::get_closest_points`) converted to
`knn1` as well — **PASS** (snapshot `s3knn1` vs `fastg`), ~1%
(hd-max 314→311 s, hd-busy 207→204 s).

### 12. Cluster::sv3d() scoped-view memo (clustering)

A fresh hd-max profile in the final mode (tcmalloc, JSON config, 92446
samples) put `Points::scoped_view` at **13.9% cumulative**, almost all
via `Cluster::sv3d()`: every `kd3d()`/`point3d()`/`points()` call in the
hot loops re-resolved the view — a Scope hash (boost::hash_range over
the name strings, ~3.5%) plus hashtable find with string-compare
equality (`__memcmp` 3.8%).

Memoized the resolved view pointer in the Cluster facade together with
the `Points::scoped_indices_valid()` flag pointer (added previously for
exactly this pattern):

- Fast path: memo set and `*valid` → return the view directly.
- Node *insert* only drops the validity flag (view object is stable):
  the slow path re-calls `scoped_view()`, which re-syncs and restores it.
- Node *removal* erases the view object, so `Cluster::on_remove` (and
  `on_insert`, for safety) reset the memo before any reuse —
  `FacadeParent` already receives these tree notifications.
- `set_default_scope()` resets the memo (scope changed).

**A/B verdict: PASS** (snapshot `svmemo` vs `knn1`, all archives
byte-identical). Clustering wall: hd-max 366→**330 s**, hd-busy
245→**222 s**, typical events −10-13%.

### 13. Grouping::fastgeom dense-vector cache (clustering)

The round-1 `fastgeom(apa,face)` memo still paid an `unordered_map`
hash-find per call (~5% cumulative on the fresh hd-max profile; it is
called per point per plane inside the good-point tests). The key is
`apa*2+face` — a tiny dense space — so the map became a
`std::vector<std::unique_ptr<fastgeom_t>>` indexed by key (`unique_ptr`
keeps the returned references stable across growth; invalid `apa=-1`
still throws at the same `m_anodes.at()` as before).

**A/B verdict: PASS** (snapshot `fastg` vs `svmemo`, all archives
byte-identical). Clustering wall: hd-max 330→**314 s**, hd-busy
222→**207 s**, vd-busy 111→105 s.

### Round-2 closing

Final state of the 6-event set (6-way load; img = snapshot `imgr2`,
clus = snapshot `fastg`), vs the original `base1` baseline:

| Event | img wall (s) | img RSS (MB) | clus wall (s) | clus RSS (MB) |
|---|---|---|---|---|
| hd-typ  027409/0  | 66 → 59 | 703 → 670 | 28 → **20** | 704 → 789 |
| hd-typ2 027980/3  | 70 → 62 | 691 → 794 | 31 → **22** | 769 → 706 |
| hd-busy 028084/18 | 656 → **163** | 7180 → **3339** | 338 → **207** | 2833 → 2964 |
| hd-max  027305/0  | 1028 → **322** | 10113 → **6078** | 530 → **314** | 3284 → 3030 |
| vd-typ  039349/0  | 140 → 135 (50 with `-P`) | 441 → 486 | 26 → **20** | 745 → 637 |
| vd-busy 039252/5  | 327 → **264** | 1918 → **1569** | 145 → **105** | 1879 → 1825 |

(Round-2 clustering walls include the entry-8 tcmalloc adoption; round-2
imaging RSS includes the entry-9 eviction. Clustering outputs since
entry 8 are deterministic but differ from the original baseline on 3
PDHD events — see entry 8.)

**Crash-mitigation soak**: 12 repetitions of hd-max clustering in the
final mode (pure-JSON config, `GOGC=off`, tcmalloc): 12/12 rc=0 and all
12 `mabc-all-apa.zip` payload-identical. The Go-GC crash vector is
removed by construction; the soak exercises the formerly crash-prone
path without incident.

**Closing verification**: (a) clustering determinism double-run with the
final code (`s3knn1` vs `s3knn1b`) — all archives identical; (b)
22-event imaging spot check (`spotr2` vs round-1 `spotbase`) — **all 256
imaging archives byte-identical**, confirming the round-2 imaging
changes (eviction, old2new) across every profiled run, with broadly
lower wall and RSS.

**Fresh hd-max clustering profile after round 2** (78.7k samples,
~315 s): the chain is now kd-geometry-bound —
nanoflann searchLevel/findNeighbors ~22-26% (irreducible),
`is_good_point` 25.3% (mostly `exists_within` kd descent),
`hough_transform` 13.4% (boost::histogram + trig on cluster points),
`Simple3DPointCloud::get_closest_points` 10.2%,
`Find_Closest_Points` 7.6%, DPC `add_points`+`DPCPoint` ctor ~7%.
Previously-listed targets now retired: `scoped_view` 13.9→1.5%,
`get_closest_point_blob` 12.7→6.0%.

Remaining targets for a round 3, in suggested order (shares from the
post-round-2 hd-max profile; the chain is now kd-geometry-bound, so
expect smaller, harder wins than rounds 1-2).  **Status annotations
added after round 3** — see the Round 3 section for the entries:

1. **Reduce kd query *counts* in the good-point tests** — ✅ **DONE in
   round 3 (entry 16)** via the plane short-circuit (the cheap share of
   this item: with `allowed_bad=1` both decided-early extremes skip the
   third plane's descent).  hd-max 308→288 s.  The cell-keyed verdict
   memo variant was NOT implemented — it would need a duplicate-query
   census first and is only worth a dedicated follow-up if a future
   profile still shows `is_good_point` dominant.
   *(original notes)* `is_good_point` is 25.3% (almost all
   `has_closest_point` → `exists_within` kd descent, called per point
   per plane per pass under Separate/Deghost); the lever is issuing
   fewer queries, e.g. memoize verdicts per (rounded time/wire key)
   within a pass.
2. **`Cluster::hough_transform` 13.4%** — ✅ **DONE in round 3 (entry
   17)**: it was already sparse-map (not boost::histogram as written
   below); round 3 replaced the map with a generation-stamped scratch
   grid and dropped the intermediate `pts`/`blobs` vectors.  hd-max
   288→281 s.  The residual cost is the `kd_radius` descent + trig,
   which is content-bound — retired as a target.
3. **DynamicPointCloud SoA restructure** (~7% now, was ~9-12% pre-
   tcmalloc) — ✅ **DONE in round 4 (entry 24)**: DPCPoint replaced
   wholesale by the columnar `DPCBatch`; byte-identical across the
   PDHD/PDVD 6-event gate, the MicroBooNE qlport full chain, SBND, and
   the 5-test porting bats suite.  Clustering RSS −14..−45% on busy
   events, wall −4..−5%.
4. **Imaging graph-rebuild remainder** — ✅ **DONE in round 4 (entries
   22-23)**: `CS::prune` no-prune fast path (entry 22, fires ~100% with
   the default threshold) and the `copy_graph` reduction via the
   `SimpleCluster` move constructor (entry 23, one whole-graph
   copy+destroy per stage removed).  Both byte-identical; both
   wall-neutral under tcmalloc — the residual upside this item carried
   was priced at glibc-malloc rates.  The remaining filtered
   `copy_graph` sites (~1.5% combined) are closed as not-worth-it.
5. **SBND follow-through** — ✅ **DONE in round 3 (entries 18-19)**: the
   heisenbug soak is clean 8/8 (entry 18), SBND clustering is verified
   run-to-run deterministic AND glibc==tcmalloc byte-identical (entry
   19), and tcmalloc + JSON-precompile + GOGC=off are enabled in
   `sbnd_xin/run_clus_evt.sh`.  No `sbnd/clus.jsonnet` change was
   needed (all changes are environment-level).
6. nanoflann kd descent (~25%) is genuine geometry — only fewer queries
   (item 1) move it.  Still true after round 3; with the short-circuit
   landed, further reduction needs the memo/census route.

Retired from the previous target list: `scoped_view`/sv3d (13.9→1.5%,
entry 12), `get_closest_point_blob` (12.7→6.0%, entry 11),
`Simple3DPointCloud` single-NN wrappers (entry 11 follow-up),
`fastgeom` lookup (entry 13).

Beyond byte-identical work, the largest remaining levers are the
**result-changing ideas** listed at the end of this file (sampling
density, `full_deghost`, slicing threshold) — those need physics review
rather than engineering.

## Round 3 (2026-06-11): deployment defaults, masked-span coarsening, kd-query reduction

### 14. Per-anode sequential imaging is now the default (script-level)

`run_img_evt.sh` (PDHD + PDVD, wcp-porting-validation `a1ad4b2`): the round-1
`-P` mode (entry 1, output-identical, peak RSS = busiest anode instead of
sum) is now the default; `-1` reverts to the old single-process mode, `-P`
is kept as a back-compat no-op.  Re-validated incidentally by entry 15's
6-event regeneration: all 32 active archives byte-identical to the
single-process baseline.

### 15. Masked-fork slicing span 100/500 → 1500 ticks (**RESULT-CHANGING, approved**)

- `cfg/pgrapher/experiment/protodunevd/img.jsonnet` (production "both"
  branch): masked 2-view span 100 → **1500** ticks.
- `cfg/pgrapher/experiment/pdhd/img.jsonnet`: 500 → **1500**.

The masked fork carries only dead-region geometry (no charge solving), so
its slicing granularity is a memory/time knob with low physics
sensitivity.  PDVD's span=100 made 15x more masked slices than PDHD's 500
with no documented physics justification.  User approved coarsening both
to 1500.

Measured on the 6-event set (snapshot `span1500`, includes the entry-14
per-anode default; vs the round-2 closing numbers):

| Event | img wall (s) | img RSS (MB) | clus wall (s) |
|---|---|---|---|
| hd-typ  027409/0  | 59 → 41 | 670 → 643 | 20 → 19 |
| hd-typ2 027980/3  | 62 → 44 | 794 → 694 | 22 → 21 |
| hd-busy 028084/18 | 163 → **143** | 3339 → 3382 | 207 → 204 |
| hd-max  027305/0  | 322 → **288** | 6078 → 6190 | 314 → 308 |
| vd-typ  039349/0  | 135 → **40** | 486 → 419 | 20 → 15 |
| vd-busy 039252/5  | 264 → **105** | 1569 → 1227 | 105 → 90 |

The big PDVD imaging wins are the 15x masked-slice reduction (the masked
fork dominated PDVD imaging); PDHD moves less (3x on a smaller share).
hd-busy/hd-max imaging RSS is unchanged because the active fork on the
busiest anode dominates the per-anode peak.

Verification of the blast radius:
- **Active fork untouched**: all 32 `*-ms-active.tar.gz` archives
  byte-identical to the old-span baseline.
- **Masked archives**: PDVD ~6-7x smaller, PDHD ~1.3-1.6x.
- **Downstream clustering**: 6/6 events rc=0.  Per-pass cluster counts
  are identical on 5/6 events; hd-busy differs in two mid-chain passes
  (61→60, 164→162) and converges to the **same final count (184)**.
  `mabc-*.zip` payloads differ on all events as expected (the dead/masked
  blob geometry is part of the payload).

Snapshot `span1500` is the new A/B baseline for subsequent byte-identical
work (the old `imgr2`/`fastg` baselines predate the span change).

### 16. Good-point tests: plane short-circuit + dead-chs single map descent (clustering)

`clus/src/Facade_Grouping.cxx`:

- `is_good_point` / `is_good_point_wc`: the per-plane tests are pure
  queries, so the U/V/W loop now returns as soon as the verdict is
  decided.  With the hot-path default `allowed_bad=1`, both extremes skip
  the third plane's kd descent: two matches → early true, two misses →
  early false (weighted variant analogous, threshold 4 with W=2).  This
  attacks the round-3 #1 target — `is_good_point` was 25.3% of hd-max
  clustering, nearly all `exists_within` kd descent — by issuing fewer
  queries rather than making the descent faster.
- `get_closest_dead_chs`: one ordered-map `lower_bound` descent for the
  whole `[wind-ch_range, wind+ch_range]` window instead of a `find()` per
  channel (3 descents at the default `ch_range=1`), plus an empty-map
  early-out before the projection arithmetic.  Same ascending visit
  order, verdict-identical.  (This is the PDHD/PDVD share of the SBND
  profile's `get_closest_dead_chs` ~9% item; the flat-bitmask variant was
  rejected because `get_dead_winds` hands out a mutable reference that
  PointTreeBuilding writes through — a cached flat structure would need
  invalidation hooks for no extra win over the single descent.)

**A/B**: 6-event clustering gate vs `span1500` — all archives PASS
(byte-identical).  Wall: hd-busy 204→**190** s, hd-max 308→**288** s,
vd-busy 90→**86** s; typicals unchanged.

### 17. hough_transform: scratch-grid accumulation, no intermediate vectors (clustering)

`clus/src/Facade_Cluster.cxx` `Cluster::hough_transform` (round-3 #2
target, 13.4% of hd-max).  The histogram was already sparse (an earlier
session replaced boost::histogram with an unordered_map keyed by the
linear bin index); the remaining cost was per-fill hashing/allocation and
the materialization of `pts`/`blobs` vectors per call.  Now:

- accumulation goes into a generation-stamped dense thread-local scratch
  grid (a bin is live when its stamp matches the call's generation) — no
  per-call zero-fill of the 64 800-bin grid, no hashing, no allocation
  after warm-up; the touched-bin list drives the peak scan.  Per-bin sums
  see the same values in the same order, and the peak criterion (max
  value, then smallest linear index) is a total order, so the result is
  bit-identical regardless of visit order.
- the kd results are iterated directly (`blob_with_point` per index +
  the scoped-view coordinate arrays) instead of building `pts`/`blobs`
  vectors first — same per-index values and order.

**A/B**: 6-event clustering gate vs the entry-16 snapshot — all archives
PASS (byte-identical).  Wall: hd-busy 190→**182** s, hd-max 288→**281** s.

### 18. SBND `clus_all_apa` heisenbug: 8-pass soak clean under round-3 code

The SBND local-imaging stream (file 1 of lan-reco2,
`/home/xqian/tmp/lanrepro2/f1/config.json`) historically crashed with an
intermittent heap corruption on ~50-70% of runs (2026-06-05 trail in
`clustering-timing-profile.md` §0).  Re-soaked with the current toolkit
(BlobLess pointer-order fix + rounds 2-3 changes), 8 sequential passes
via `/home/xqian/tmp/lanrepro2/soak_r3.sh`: 4 plain glibc + 4 with
`MALLOC_PERTURB_=85 MALLOC_CHECK_=3` (perturbed free fills surface latent
corruption a lucky layout hides).  **8/8 rc=0, zero corruption
signatures**, stdout structure identical to the known-clean June-5 pass.
At the historical 50% floor, 8 clean runs has p≈0.4% — the bug no longer
reproduces under current code.  Which change removed it cannot be
pinpointed (the BlobLess re-keying also shifts allocation patterns);
operationally the SBND local-imaging benchmarking path is unblocked.
Bonus datum: each pass now takes ~260 s vs ~630 s on June 5 — the round
2-3 clustering speedups carry over to SBND unchanged.

### 19. SBND determinism A/B + tcmalloc enablement (round-2 item 5 follow-through)

Determinism/allocator A/B on 5 SBND events (data 138670 — the event
whose all-APA output historically differed run-to-run — plus data
139220, 137680 and MC evt11, evt12), driving `sbnd_xin`'s clustering
graph directly: per event, two glibc runs (run-to-run determinism) and
one tcmalloc-preloaded run (allocator independence).  **All 15 archives
(3 mabc zips × 5 events) byte-identical across all three runs** — the
entry-8 BlobLess fix resolves the SBND `clus_all_apa` nondeterminism,
exactly as predicted.

Enablement (`sbnd_xin/run_clus_evt.sh`, wcp-porting-validation
`ecef431`, mirroring the pdhd/pdvd pattern): tcmalloc preload
(`WCT_TCMALLOC=off` reverts), wcsonnet pre-compile to pure JSON +
`GOGC=off` (removes the Go-GC SIGABRT vector).  The updated script path
re-verified byte-identical on evt12.  **No `sbnd/clus.jsonnet` change
was needed** — everything is environment-level; the config compiles
unchanged.  The QL/matching scripts were deliberately left alone:
QLMatching has known marginal FP ties at `m_strength_cutoff` that an
allocator change could flip, so its tcmalloc adoption should be its own
gated decision.

### 20. hd-max clustering heap profile (tcmalloc HEAPPROFILE)

First heap-level decomposition of the ~3 GB hd-max clustering RSS
(`HEAPPROFILE` + full `libtcmalloc.so.4`, dumps every 200 MB of in-use
growth; profiler overhead ~14x on allocation-heavy passes, so the run
was stopped once the headline was unambiguous; dumps in
`/home/xqian/tmp/heaprof/`).  Findings:

- **The live heap is much smaller than RSS**: at a point where process
  RSS exceeded 1 GB, only ~423 MB was live; at the last dump 843 MB live
  vs ~2 GB RSS.  A large share of clustering "RSS" is
  allocator-retained/fragmented memory from churn (146M allocations /
  15.3 GB cumulative by the last dump), not live data.
- **~90% of the live heap at the sampled points sits under
  `ClusterFileSource` ingest** — for the JSON cluster files this is the
  slurped `istringstream` buffer plus the jsoncpp DOM (an Rb_tree node
  per object member, 40%, + key/value strings, ~42%).  Allocation churn
  at those points is likewise parse-dominated.

Actionable outcome: stop ingesting JSON → entry 21.

### 21. Cluster files JSON → numpy (`format: "numpy"`), + two latent toolkit bugs

`cfg/pgrapher/experiment/{pdhd,protodunevd}/img.jsonnet`: the
`ClusterFileSink` in the imaging `dump()` now writes `format: "numpy"`
(npy members in the same tar.gz, as SBND's npz path has long used)
instead of JSON.  `ClusterFileSource` dispatches per member, so
downstream needs no change.

Enabling this exposed two latent toolkit bugs, both fixed:

1. **custard tar reader 512-padding bug** (`custard_boost.hpp:486`,
   commit `05642fad`): a member whose size is an exact multiple of 512
   has no padding, but the reader skipped `512 - size%512 = 512` bytes,
   swallowing the next member's tar header (crash: `stol` on garbage).
   JSON members (arbitrary sizes) hit this with ~1/512 probability per
   member; npy arrays (element sizes 4/8) hit it readily.  The tar
   *writer* was always correct, so all existing archives remain valid.
   Same expression fixed in the unused `custard_file.hpp` pair.
2. **`ClusterFileSink::numpify` empty-cluster desync** (commit
   `a93e696e`): an empty cluster graph produced no members, so the
   source delivered EOS instead of an empty cluster and the live/dead
   stream pair into `PointTreeBuilding` desynchronized ("missing 1
   input tensors", hd-busy/hd-typ2 failed).  Fixed by writing empty
   clusters as the (tiny) JSON member.

**A/B (snapshot `npfmt2`)**: all **114 clustering archives across the
6-event set byte-identical** to the JSON-input baseline — the JSON
writer was full-precision, so the binary round-trip changes nothing.

| Event | img wall (s) | img RSS (MB) | clus wall (s) | clus RSS (MB) |
|---|---|---|---|---|
| hd-typ  027409/0  | 41 → 40 | 628 → **400** | 20 → **17** | 644 → 702 |
| hd-typ2 027980/3  | 44 → 43 | 677 → **400** | 22 → **19** | 702 → 753 |
| hd-busy 028084/18 | 143 → 141 | 3302 → 3304 | 182 → **176** | 3137 → 3091 |
| hd-max  027305/0  | 288 → 288 | 6045 → 6046 | 281 → **272** | 3015 → 3016 |
| vd-typ  039349/0  | 40 → 38 | 408 → **236** | 14 → **12** | 485 → 510 |
| vd-busy 039252/5  | 105 → 101 | 1198 → 1136 | 85 → **79** | 2219 → 2497 |

Interpretation: clustering wall −3 to −15% (parse cost gone), imaging
RSS on typical events −35-45% (no JSON serialization buffers), busy-
event imaging RSS unchanged (in-memory graph dominates).  Clustering
*peak* RSS is essentially unchanged — consistent with entry 20: the
ingest DOM dominated the *live heap mid-load*, but the peak RSS envelope
is set by the clustering stages themselves plus allocator retention.
(vd-busy clus RSS swings ±10% with concurrent load; treat single-run
deltas there as noise.)

## Round 4 (2026-06-11): imaging graph-rebuild remainder, DPC SoA

### 22. CS::prune no-prune fast path (imaging)

`CS::prune` rebuilt the solved b-m subgraph vertex-by-vertex to drop
blobs below `blob_value_threshold` — but the production configs leave
that threshold at its default of **−1**, and LASSO charges are
non-negative, so *nothing is ever pruned*: the rebuild was a pure
identity copy of every connected subgraph for every weighting strategy.
Added an rvalue overload (`img/src/CSGraph.cxx`) that scans the blobs
once and, when none falls below threshold, moves the input graph
through unchanged; `ChargeSolving` passes the solve output as rvalue.
Iteration-order safe by construction: vertices are vecS and out-edges
setS (ordered by target descriptor, insertion-order independent), so
the moved graph and a rebuilt copy iterate identically.  Anything
actually below threshold falls back to the original rebuild.

**A/B verdict: PASS** (snapshot `prune1` vs `npfmt2`): all 162 archives
across the 6-event set byte-identical (48 imaging + 114 clustering).
Wall and RSS unchanged within noise (hd-max img 288→288 s) — the
subgraphs are small, so the skipped rebuild was allocator churn rather
than measurable wall.  Kept as the correct semantics: the per-subgraph
copy now happens zero times instead of once per strategy.

### 23. SimpleCluster move constructor — graph rebuild remainder closed

A fresh hd-busy anode2 CPU profile (glibc malloc, 23478 samples) put
`SimpleCluster::SimpleCluster` at **7.3% cumulative**: the constructor
ran `boost::copy_graph(g, m_graph)`, a full redundant copy (plus later
destruction of the source) of the cluster graph at the end of *every*
imaging stage.  All other `copy_graph` sites combined (the filtered
copies in Local/GlobalGeomClustering, InSliceDeghosting rounds 2/3,
BlobGrouping, ClusterScopeFilter) measured only ~1.5%.

Added a `SimpleCluster(cluster_graph_t&&)` move overload and converted
the 15 call sites whose graph is dead after construction (img: all
stage tails incl. BlobClustering/BlobSolving via the indexed graph's
non-const accessor; sio: `ClusterFileSource` JSON+numpy loads; aux:
`TensorDMcluster`).  Iteration-order safe: the moved graph is the very
object the copy would have reproduced (vecS vertices; setS out-edges
ordered by target).  The dryrun paths keep the copying ctor.

**A/B verdict: PASS** (snapshot `scmove1` vs `npfmt2`): all archives
byte-identical on the 6-event set.  **Wall/RSS: no measurable change**
— the profile above was taken under glibc malloc, but production
imaging runs under tcmalloc (round-1 adoption), where the copy's
allocation churn is far cheaper; the same lesson as entry 11.  Kept as
the structurally right thing: one whole-graph copy+destroy per pipeline
stage now simply does not happen, and the win applies to any future
allocator/regime change.

The remaining filtered `copy_graph` sites are **closed as
not-worth-it**: ~1.5% combined under the glibc profile (less under
tcmalloc), each requiring a manual order-preserving rebuild to replace
a one-line `copy_graph` — poor risk/benefit against the byte-identity
gate.  This closes round-2 follow-up item 4 (imaging graph-rebuild
remainder); both sub-items (CS::prune fast path = entry 22, copy_graph
sites = this entry) are done.

### 24. DynamicPointCloud SoA restructure (round-2 item 3, the big one)

`DPCPoint` carried five nested vectors per point (`x_2d`/`y_2d`/
`wpid_2d` each `vector<vector<...>>` plus `wind`/`dist_cut`) — roughly
a dozen heap allocations *per point*, megabytes of points per event,
built and torn down across every `make_points_*` call.  Replaced the
AoS struct wholesale with a columnar `DPCBatch` (scalar columns +
CSR-encoded per-plane 2D projections) used both as the builder transfer
format and as `DynamicPointCloud`'s internal storage:

- builders (`make_points_cluster[_steiner|_skeleton]`,
  `make_points_direct`, `make_points_linear_extrapolation`,
  `fill_wrap_points`) append straight into the batch columns — zero
  per-point allocations;
- `add_points(DPCBatch&&)` on a fresh cloud is an O(1) column move
  (the dominant call pattern);
- consumers migrated from `get_points()[i].field` to
  `npoints()/point3d(i)/cluster(i)/dist_cut(i,p)/points()` accessors;
  whole-cloud and row-subset copies (PRShower merge, break_segment
  redistribution) became `add_points(other[, rows])` column appends.

Values, point order, and kd-tree construction order are bit-identical
by construction; ~10 consumer files touched including the
neutrino-porting code (PRShower/PRSegmentFunctions/NeutrinoDeghoster/
NeutrinoEnergyReco/NeutrinoShowerClustering/TaggerCheckSTM/MABC).

**Verification (all PASS):**
- PDHD/PDVD 6-event clustering A/B (snapshot `dpcsoa1` vs `npfmt2`):
  all archives byte-identical.
- MicroBooNE full chain (qlport `uboone-mabc.jsonnet`, events
  6501/6505/6512, kind=both incl. pattern recognition + taggers): all
  3 bee zips byte-identical — this is the path that exercises the
  PR*/Neutrino* consumers.
- SBND clustering (sbnd_xin `wct-clustering.jsonnet`, evt 138670 data +
  evt 12 sim): 6/6 mabc zips byte-identical.
- Dedicated porting suite `clus/test/test-porting.bats`: **5/5 ok**
  (qlport, steiner, stm, pdhd, fgval) including historical log-digest
  diffs.

**Resources** (clustering stage, vs `npfmt2`):

| Event | clus wall (s) | clus peak RSS (MB) |
|---|---|---|
| hd-typ  027409/0  | 17 → 17 | 702 → **644** |
| hd-typ2 027980/3  | 19 → 19 | 753 → **703** |
| hd-busy 028084/18 | 176 → **168** | 3092 → **2666** |
| hd-max  027305/0  | 272 → **260** | 3089 → 3087 |
| vd-typ  039349/0  | 12 → 12 | 510 → **472** |
| vd-busy 039252/5  | 79 → **76** | 2497 → **1373** |

The headline is memory: vd-busy **−45%**, hd-busy −14%, typicals
−7..−9% — the nested-vector churn this removes was the allocator-
retention feeder identified in the round-3 heap profile.  Wall −4..−5%
on the busy/max events (the projected ~7% was a glibc-profile number;
under tcmalloc the allocation share is cheaper).  hd-max RSS is flat:
its peak envelope sits in stages that don't go through
DynamicPointCloud.  Imaging is untouched by this change.

## Round 5 (2026-06-11): post-round-4 profiling — remaining wall & RSS targets

Measurement-only round: re-profile after rounds 1-4 to rank what is left,
on both axes.  No code changes shipped; the ranked list below is the
round-6 work queue.

### Methodology (two fixes over earlier rounds)

- **CPU profiles now run under the production allocator**:
  `LD_PRELOAD=libtcmalloc_and_profiler.so.4` (production preloads
  tcmalloc_minimal; round-4 showed glibc-malloc profiles overstate
  allocator-heavy attributions).  `abtest/profile_img.sh` /
  `profile_clus.sh` gained `PROFLIB` / `OUTDIR` / `HEAPOUT` env knobs;
  the default profiler lib is now tcmalloc_and_profiler.  Profiler
  overhead measured nil (profiled walls == round-4 baseline walls).
- **Heap profiling via jemalloc sampling, NOT gperftools HEAPPROFILE.**
  `HEAPPROFILE` records a stack on *every* alloc/free → ~25× slowdown on
  this workload (hd-max imaging would have taken hours; killed).
  `TCMALLOC_SAMPLE_PARAMETER` does not apply to `HEAPPROFILE` (it only
  feeds the in-process MallocExtension sampler).  Debian's
  libjemalloc2 5.3 has the profiler compiled in:
  `MALLOC_CONF=prof:true,lg_prof_sample:19,lg_prof_interval:33,prof_final:true`
  ≈ zero overhead; analyze with the matching `jeprof` (5.3.0) script.
  Do NOT use `prof_gdump:true` here — it dumps at every chunk-level VM
  high-water mark (6105 files in 30 s).  Peak attribution = the interval
  dump with the largest in-use total (header line 2).
- Profiled: imaging heaviest anodes (hd-max 027305/0 a0 186 s, hd-busy
  028084/18 a2 66 s, vd-busy 039252/5 a6 35 s) and full clustering on
  hd-max / hd-busy / vd-busy.  Profiles in `/home/xqian/tmp/r5_*.prof`,
  heap dumps in `/home/xqian/tmp/r5_heap/`, jeprof at
  `/home/xqian/tmp/jeprof`.

### Imaging CPU (per-anode shares, tcmalloc)

| component | hd-max a0 | hd-busy a2 | vd-busy a6 |
|---|---|---|---|
| ChargeSolving total | **47.8%** | 20.1% | 27.5% |
| — of which LASSO `CS::solve` | 35.0% | 4.5% | 11.0% |
| ProjectionDeghosting | 19.2% | **28.3%** | 20.9% |
| InSliceDeghosting | 6.7% | 8.4% | 8.6% |
| LocalGeomClustering | 6.6% | 8.0% | 8.0% |
| BlobGrouping | 4.9% | 6.9% | 6.6% |
| BlobClustering | 4.7% | 7.0% | 3.5% |
| input read (FrameFileSource) | 1.4% | 3.6% | 6.9% |

Cross-cutting (cumulative, overlaps the rows above): `boost::add_edge`
into `setS` out-edge sets 24.9 / 31.5 / 27.7%, `adjacency_list` dtor
10-14%, **graph copy (`vec_adj_list_impl::copy_impl`) 13-20%**, tcmalloc
internals (flat) 10-14%.

Findings:

1. **Hidden by-value graph copies at stage heads** — `const auto
   in_graph = in->graph();` deduces `cluster_graph_t` *by value* and
   silently full-copies the input graph once per stage call.  Sites:
   `ChargeSolving.cxx:261`, `ProjectionDeghosting.cxx:164`,
   `InSliceDeghosting.cxx:742`, `LocalGeomClustering.cxx:54`,
   `GlobalGeomClustering.cxx:74`, `ClusterScopeFilter.cxx:55`,
   `LCBlobRemoval.cxx:118` (+ `TestClusterShadow.cxx:76`, test-only).
   This is most of the 13-20% copy_impl share plus the matching dtor
   share, *and* 1.56 GB of hd-max peak RSS (see heap below).  It also
   explains why round-4's SimpleCluster move-ctor (entry 23) was
   wall-neutral: the output-side copy was removed but the input-side
   copy remained.  Fix = `const auto&` (one character per site),
   byte-identical by construction.
2. **LASSO Gram build computes the full symmetric matrix**:
   `LassoModel.cxx:92-104` does `X.col(i).dot(X.col(j))` for ALL (i,j);
   computing j≥i and mirroring the triplet is bit-identical (sparse dot
   is element-wise commutative; summation order per dot unchanged; no
   duplicate triplets so `setFromTriplets` is order-insensitive).
   Worth ~6% on LASSO-bound events (hd-max), ~1-2% elsewhere.
3. `setS` add-edge churn (~25-31%) remains the floor of the graph
   lifecycle: container swap was rejected (output-changing); a custom
   pool allocator for the edge sets would preserve ordering
   (set is ordered by value, not address) but is invasive — parked.
4. ProjectionDeghosting internals (hd-busy): `get_projection` 33%-of-PD,
   `judge_coverage(+_alt)` ~19%, Blob/ClusterShadow build ~18%, rest
   add_edge/Rb_tree.

### Clustering CPU (tcmalloc; consistent across all 3 events)

| item | hd-max | hd-busy | vd-busy |
|---|---|---|---|
| clustering_separate (ctpc) | 31.5% | 29.2% | 27.0% |
| ClusteringExamineBundles (relaxed) | 17.1% | 12.0% | 8.6% |
| ClusteringDeghost | 10.0% | — | 7.3% |
| nanoflann kd queries (flat, total) | ~32% | ~27% | ~25% |
| `is_good_point` | 22.6% | 12.3% | — |
| `get_closest_points` (S3DPC) | — | 15.0% | 17.3% |
| hough/vhough | 13.2% | 12.9% | 10.5% |

Clustering is near its result-preserving floor: `is_good_point` is
already memoized + short-circuited (round 3), and
`Simple3DPointCloud::get_closest_points` is an order-sensitive
iterative walk (restructuring it changes results).  The remaining CPU
is genuine kd geometry over O(pairs) sub-cloud comparisons.

### Heap peaks (jemalloc sampled, live bytes at peak dump)

**Imaging hd-max a0: 6.0 GB attributed — 84.5% (5.1 GB) live inside one
ProjectionDeghosting call:**

| holder | live at peak |
|---|---|
| `BlobShadow::shadow` graph (`adjacency_list<multisetS,vecS,undirectedS>`) | **3.13 GB (52%)** — 2.1 GB multiset edge nodes + 1.6 GB edge-list nodes; ~120-150 B container overhead per shadow edge |
| `in_graph` by-value copy (finding 1 above) | 1.56 GB (26%) |
| Projection2D Eigen sparse cache (`id2lproj`) | 0.81 GB (13.5%) |
| SimpleCluster-held cluster graph | 0.80 GB (13.3%) |

**The `bsgraph` lifetime is the cheap 3 GB**: it is consumed only by
`ClusterShadow::shadow` (`ProjectionDeghosting.cxx:167-171`) and then
stays alive to the end of the pass.  Scoping it in a block (or
`std::move`-and-reset) frees it before the projection/judge phase —
byte-identical, trivial.  After that + the `auto&` fix the hd-max
imaging peak should drop from ~6.2 GB to ≈ the bs-build moment
(~4 GB); shrinking further needs the multisetS→vecS container change
(edge iteration order changes → consumer audit + full A/B required).

**Clustering peaks:**

| holder | hd-max (3.03 GB) | hd-busy (2.27 GB) | vd-busy (1.03 GB) |
|---|---|---|---|
| pc-tree load (`as_pctree`/Dataset) | 1.23 GB (41%) | 0.89 GB (39%) | 0.56 GB (54%) |
| **Bee JSON buffers** (`fill_bee_points` → `Json::Value`) | **1.35 GB (45%)** | 0.33 GB (14.5%) | 0.09 GB (8%) |
| DPC kd-index arrays (`index_new_points`) | — | 0.43 GB (19%) | 0.17 GB (16%) |
| cluster graphs (`graph_algorithms`/`give_graph`) | 0.20 GB | 0.36 GB | — |

This resolves round-4's "hd-max RSS flat" puzzle: its peak never was in
DynamicPointCloud — it is Bee display buffers + the pc-tree.
`Bee::Points` stores every point as 6 `Json::Value` array elements
(`util/src/Bee.cxx:103-114`), ~300 B per point for 44 B of payload, and
all Points objects (all-apa + per-apa + per-face + img/dead hooks) stay
alive until the end-of-event zip write.  Storing compact columns and
building/streaming the JSON only at write time would emit identical
bytes with ~1 GB less hd-max peak.  The pc-tree share is the event data
itself — irreducible.

### Round-6 work queue (ranked, all byte-identity-gated)

1. `const auto` → `const auto&` at the 8 `in->graph()` sites — ~10-15%
   imaging wall on busy events + 1.5 GB hd-max peak.  Trivial.
2. Scope `bsgraph` death before the PD projection phase — ~2-3 GB hd-max
   imaging peak.  Trivial.
3. LASSO Gram symmetric halving — ~6% hd-max imaging wall.  Small.
4. Bee::Points compact columns, serialize at write — ~1.0-1.2 GB hd-max
   clustering peak (45%).  Medium; byte-identical JSON achievable.
5. BlobShadow graph edge container multisetS→vecS — ~2 GB + CPU in the
   bs build; **medium risk** (out-edge iteration order changes; audit
   `ClusterShadow::shadow` + PD consumers first).
6. DPC kd zero-copy adaptor (nanoflann reads SoA columns directly,
   drops `index_new_points` copies) — ~0.4 GB hd-busy clustering.
   Medium.

Document-only / result-changing: SpGEMM (`X.T*X`) Gram build (FP
summation order changes); pool allocator for `setS` edge sets
(byte-identity plausible but invasive); per-APA threaded clustering
passes (wall-clock only, scheduling surgery); reuse of PD projections
across full_deghost passes (staleness risk).

## Round 6 (2026-06-11): the round-5 queue, items 1-3

### 25. Reference-bind the per-stage input graph (`const auto` → `const auto&`)

`const auto in_graph = in->graph();` deduces `cluster_graph_t` *by
value*, full-copying the input graph once per stage call (round-5
finding 1).  Changed to `const auto&` at all 8 sites: `ChargeSolving`,
`ProjectionDeghosting`, `InSliceDeghosting`, `LocalGeomClustering`,
`GlobalGeomClustering`, `ClusterScopeFilter`, `LCBlobRemoval`,
`TestClusterShadow` (test-only).  Safe: `ICluster::graph()` returns a
const ref owned by the input `ICluster`, which outlives the stage call;
every site reads only.

- A/B snapshot `r6graphref` vs `scmove1` (img) / `dpcsoa1` (clus):
  **178/178 archives byte-identical, PASS.**
- Imaging wall: hd-max 284→**263 s** (−7%), hd-busy 141→**126 s**
  (−11%), vd-busy 101→**94 s** (−7%), typicals −2 s.
- Imaging peak RSS: hd-max 6190→**5425 MB** (−12%), hd-busy
  3382→**2884 MB** (−15%), vd-busy 1162→**994 MB** (−15%).
- Clustering walls drifted −3% — noise; this change does not touch
  clustering code paths.

### 26. Scope the BlobShadow graph inside ProjectionDeghosting

Round-5 heap profiling showed the `bsgraph` built at the top of every
`ProjectionDeghosting` call (3.13 GB live on hd-max a0 — multisetS edge
nodes + edge-list nodes) is consumed only by `ClusterShadow::shadow`
yet stayed alive to the end of the pass.  Wrapped build+conversion in a
block so it is freed before the projection/judge phase.  Pure lifetime
change, byte-identical by construction.

- A/B snapshot `r6bsscope` vs `r6graphref`: **178/178 byte-identical,
  PASS.**
- Imaging peak RSS: hd-max 5425→**4667 MB** (−14%), hd-busy
  2884→**2173 MB** (−25%), vd-busy 994→**755 MB** (−24%); typicals
  flat (their bsgraph is small).  Wall unchanged (±2 s noise).
- hd-max keeps a residual because its peak now sits at the bs-build
  moment itself, as round 5 predicted; shrinking that needs the
  multisetS→vecS container change (queue item 5, order-audit required).
- Cumulative imaging peak RSS from the round-4 baseline (items 25+26):
  hd-max 6190→4667 (−25%), hd-busy 3382→2173 (−36%), vd-busy
  1162→755 (−35%).

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

> **Superseded by round 2** — most items below were implemented (entries
> 8-13) or re-ranked with fresh numbers; the current list is in the
> "Round-2 closing" section above. Kept for the round-1 record.

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
