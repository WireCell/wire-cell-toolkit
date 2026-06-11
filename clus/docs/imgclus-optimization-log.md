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
2. **tcmalloc deployment** (`LD_PRELOAD` in run scripts): with ~140 GB of
   churn, glibc retains the high-water mark forever while live memory drops
   to ~1 GB; tcmalloc returns freed pages. Helps the *sustained* footprint
   under 16-way concurrency (not VmHWM of the early peak). Must pass the
   full A/B gate first (allocator changes can flip pointer-keyed iteration
   order if any exists).
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
