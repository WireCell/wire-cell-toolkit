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
> already-fixed `fill_wrap_points` off-by-one). Until fixed, the §0 wall/RSS come from runs
> that happened to complete all 50 events. See the debugging trail below.
>
> **Debugging trail (what has and hasn't reproduced it) — for future work:**
> The bug is rare and layout-sensitive; every "definitive" tool so far has come up empty,
> so do **not** read these clean results as "no bug" — read them as "this probe didn't
> trip it." Catalogue of attempts on file 1 (the only file that has ever crashed; files 2–3
> = 50 events each, always clean):
> - **Native, single event 60119 in isolation:** 25/25 clean → the crash is cross-event /
>   state-dependent, not a property of one event's data.
> - **gdb:** never reproduces (layout shift hides it) — classic heisenbug signature.
> - **valgrind memcheck** (Go/gojsonnet startup noise suppressed via `go.supp`, JSON config
>   to skip the Go eval path): ~35 min run, **0 errors**. Consistent with a *within-capacity*
>   overrun (write stays inside an allocation's real heap block → invisible to memcheck).
> - **`-D_GLIBCXX_ASSERTIONS`** (libstdc++ bounds-checked `vector::operator[]` etc.): 4/4
>   file-1-alone passes **silent** → not a checked-container overrun; rules out the obvious
>   `std::vector[]` OOB family.
> - **AddressSanitizer** (`-fsanitize=address` on `WireCellClus` only; red zones catch
>   heap-buffer-overflow + use-after-free even within-capacity-tolerant native): **12 parallel
>   workers × 12 independent ASLR layouts, each ~43+ of file 1's 50 events, 0 catches.** No
>   red-zone hit, no captured SEGV. (~2.2 min/event under ASan, ~7.5 GB RSS/worker.)
> - **Net:** not reproduced under any instrumented build this session. The crash is rarer
>   than ~1-in-12-layouts-per-pass under ASan, or its bad access lands inside a live
>   allocation on the layouts probed. Next angles untried: longer ASan soak (many passes,
>   not one wave), `MALLOC_PERTURB_`/`MALLOC_CHECK_` on the native build, the busy-data
>   correlation from the [[fill_wrap_points OOB]] memory (crash favors busy events →
>   instrument the wire/channel-projection writes on W-dead-recovery blobs specifically),
>   or a TSan pass if any clustering state is touched off-thread. All temporary diagnostic
>   build flags were reverted after the session (`waft/smplpkgs.py` back to HEAD; clus
>   rebuilt with 0 `__asan` symbols).

## 1. Local-imaging per-event re-measurement (2026-06-05) — the foreign hotspots dissolve

The five "pathological" events from the foreign-input headline (below) were re-measured
under **our own local imaging**. Four already had complete timing on disk from the clean
full-file `run_local.sh` reprocesses (`f2`/`f3` rc=0, 50 evt each; `f1` 59789 from the
partial run); **61637** (never reached before the f1 heisenbug crash) and **59789** (its
on-disk number came from the crashed f1 run) were re-run **single-event in isolation**
(crash-safe: the `clus_all_apa` heisenbug is cross-event) via
`/home/xqian/tmp/lanrepro2/{mkev.sh,runev.sh}` → `f1/ev<ID>/ev.log`. 59789's single-event
total (55.2 s) matches its full-file total (54.7 s), confirming the single-event method
and validating the crashed-run number.

| rank (old) | event | file | **foreign baseline** (superseded) | **local imaging** (total clustering / ExamineBundles) |
|--:|--:|:--:|--:|--:|
| 1 | **141530** | 2 | 952 s (ExBundles 613 s) | **6.1 s** / 2.9 s |
| 2 | **59789**  | 1 | 466 s (ExBundles 239 s) | **55.2 s** / 29.7 s  *(single-event re-run; full-file 54.7 s)* |
| 3 | **138824** | 2 | 236 s (ExBundles 150 s) | **47.5 s** / 29.7 s |
| 4 | 187650 | 3 | ~92 s (ExBundles 46 s) | **56.7 s** / 39.1 s |
| 5 | 61637  | 1 | ~72 s | **4.0 s** / 2.3 s  *(single-event)* |

Per-file local-imaging totals: **f2 = 144.5 s / 50 evt, f3 = 173 s / 50 evt** (≈ 3 s/event
mean vs ≈ 17 s/event on the foreign baseline). ExamineBundles is **still the #1 stage
(50–54 %)** and still long-tailed, but the worst *single* node-pass is now ~15–20 s, not
~310 s.

**Two premises overturned.**

1. **141530 — the isochronous event in the screenshot
   (`sbnd_xin/pics/Screenshot 2026-06-05 at 6.26.44 PM.png`) — is NOT a hotspot under our
   pipeline** (6.1 s; per-node ExamineBundles 228 + 1213 + 1422 ms). The "isochronous look"
   is not what drives clustering time. The genuinely slow events are 59789 / 138824 / 187650.

2. **Neither cluster count nor blob/point count predicts ExamineBundles time.** Final live
   `nclusters` is tiny for *every* event (≤ 30; 141530 has the most at 30, yet is fastest),
   so the O(n²) cluster *pair* loop (`clustering_examine_bundles.cxx:113-129`) cannot be the
   explosion — it ranges over ≤ 30 items. And raw blob/point counts (from the active npz)
   *anti*-correlate if anything:

   | event | clustering | total points | total blobs |
   |--:|--:|--:|--:|
   | 141530 | **6.1 s** | 8853 | **9187** (most) |
   | 138824 | 47.5 s | 8178 | 6689 |
   | 187650 | 56.7 s | 7436 | 6014 |
   | 59789  | 55.2 s | 6964 | 6567 |
   | 61637  | **4.0 s** | 4939 | 2873 |

   141530 has the most blobs and points but the fastest clustering. Even within a TPC:
   141530/apa0 (5666 blobs) = 1.2 s vs 59789/apa0 (4886 blobs) = 14.6 s — comparable blob
   count, **12× the time**. So the cost is **topological** (the per-cluster
   `connected_blobs` blob-graph: edge density / how tangled the blobs are in a few large
   clusters), **not** count-driven.

**Attribution caveat (do not over-claim the 5–150× drop).** The foreign baseline and these
local runs differ in *three* uncontrolled ways — imaging input (foreign larsoft
`*-active.npz` vs our multi-3view active), code version (the `fill_wrap_points` fix landed
between), and config (foreign was auto_mask-on; `run_local.sh` is not). So treat the local
numbers as the **live-pipeline measurement that supersedes** the foreign baseline; the gap
is not attributable to any single cause. (An isolation control — current binary on foreign
141530 input — would pin imaging-vs-code, but is not needed for the optimization decision.)

### Revised optimization strategy (replaces the foreign-baseline "spatial pre-filter" plan)

The original suggestion (spatial pre-filter on the O(n²) cluster pair loop) is **misdirected**
under local imaging: with ≤ 30 final clusters that loop is negligible. The real cost is the
per-cluster `connected_blobs` (`clustering_examine_bundles.cxx:148-162` →
`Facade_Cluster.cxx:2949`) on a few large, densely-connected clusters, and it scales with
blob-graph **edge count**, not blob count.

1. **Split the ExamineBundles timer** into *pair-loop* vs *connected_blobs* (and, inside
   `connected_blobs`, graph-build vs connected-components) on one slow event (138824 —
   ExamineBundles concentrated in `apa1-0` + `clus_all_apa` ≈ 30 s). A small debug-only
   `took` timer in `MultiAlgBlobClustering.cxx` / `Facade_Cluster.cxx`; no behavior change.
   This is the prerequisite measurement — it tells us *which half* to optimize.
2. **Profile `connected_blobs` per cluster**: log per-cluster blob count + graph edge count
   + time, to confirm a handful of large tangled clusters dominate and find the edge-count
   blow-up. The blob-graph build (proximity search via the `DynamicPointCloud` kd-tree) and
   the edge set are the candidate hot loops.
3. `connected_blobs` is the **shared primitive** under Deghost / ProtectOverclustering /
   Separate (which peak on the same events), so one fix there plausibly moves all four.
4. Validate any change is **clustering-output-identical** on the quiet hand-scan events
   before trusting the tail events.

> **Resolved by §2 below (2026-06-05):** the gperftools profile of event 187650 answers
> steps 1–2 directly. The cost is **not** in `connected_components` and **not** edge-count
> per se — it is a `boost::format` scope-key string rebuilt millions of times inside
> `connect_graph_relaxed`. See §2.

## 2. Root-cause CPU profile of event 187650 (gperftools, 2026-06-05) — the real hotspot is a format() in a hot loop

**Method (no rebuild).** Ran the single-event 187650 matching graph (local active +
in-graph dead) under the gperftools CPU profiler — `LD_PRELOAD=libprofiler.so.0`,
`CPUPROFILE=cpu.prof`, `CPUPROFILE_FREQUENCY=250` — then analyzed with `google-pprof`
against `build/apps/wire-cell` (full symbols; 17434 samples ≈ 70 s CPU). Repro in
`/home/xqian/tmp/lanrepro2/f3/ev187650/` (`cpu.prof`, `evprof.log`).

**Call tree (cumulative % of whole-process CPU):**

```
MultiAlgBlobClustering              85.7%   (clustering dominates; imaging/IO are noise)
  ClusteringExamineBundles          58.9%   -> connect_graph_relaxed -> get_closest_points
      make_graph_relaxed            58.8%
        connect_graph_relaxed       57.9%   <-- the long-range "relaxed" bridging
          get_closest_points        57.6%
            kd2d                     47.3%   <-- 2D kd-tree accessor by (apa,face,plane)
              String::format        40.8%   <-- !!! boost::format rebuilding a scope key
            test_good_point         43.4% / is_good_point 19.3%
          vhough_transform          14%     <-- ALSO under connect_graph_relaxed (Hough direction)
  ProtectOverclustering             14.8%   -> Separate_overclustering -> kd2d (66% of Protect)
  ClusteringDeghost                  8.0%   -> connect_graph_ctpc -> kd2d (62% of Deghost)
```

**`kd2d` is the shared hotspot** (verified by the upward caller chain): its 47.3% splits
**ExamineBundles 71% / ProtectOverclustering 21% / Deghost 8%** (Examine & Protect reach it
via `connect_graph_relaxed`, Deghost via `connect_graph_ctpc`). So fixing `kd2d` speeds all
three at once. The other big item, `vhough_transform` (~14% here), is **also inside
`connect_graph_relaxed`** (the per-candidate Hough direction) — not a separate stage — and is
the next hotspot once the format cost is removed (see §3).

**The hotspot is `Grouping::kd2d` (`clus/src/Facade_Grouping.cxx:297-306`):**

```cpp
const Grouping::kd2d_t& Grouping::kd2d(const int apa, const int face, const int pind) const {
    std::vector<std::string> plane_names = {"U", "V", "W"};                       // heap alloc / call
    const auto sname = String::format("ctpc_a%df%dp%d", apa, face, plane_names[pind]); // line 300: 7171 samples = 41%
    Tree::Scope scope = {sname, {"x", "y"}, 1};
    const auto& sv = m_node->value.scoped_view(scope);                            // line 303: 781 samples (4.5%)
    return sv.kd();
}
```

`kd2d` is called per-point-per-plane inside `connect_graph_relaxed`'s good-point tests
(`get_closest_points` → `test_good_point`/`is_good_point`), i.e. **millions of times**, and
every call rebuilds — via slow `boost::format` — one of only **≤ 12 distinct, stable**
`(apa, face, pind)` scope keys, plus a fresh `vector<string>{"U","V","W"}`. Line 300 alone
is **41% of whole-process CPU (~48% of clustering)**; `kd2d` total is **47.3%** (the
`__dynamic_cast` 10%-self in the flat profile is *also* under `kd2d → scoped_view`).

**This explains the §1 paradox** (blob count doesn't predict time). The overhead is a flat
per-good-point-test cost, so an event's clustering time tracks **how many long-range
good-point tests `connect_graph_relaxed` performs**, which depends on cluster geometry/extent
(how much candidate bridging it probes), *not* raw blob count. The isochronous 141530 — many
blobs but geometry that triggers far fewer relaxed probes — is fast; the tangled 59789 /
138824 / 187650 probe far more.

### Fix (IMPLEMENTED 2026-06-05) — memoize the kd2d scope key

The per-call `boost::format` + `vector<string>` are replaced by a memoized, data-independent
`(apa,face,pind) → Tree::Scope` lookup (`mutable m_kd2d_scope_cache` on `Grouping`); the
`scoped_view(scope)` call is unchanged, so the tree's own kd-tree caching/invalidation — and
thus the clustering result — is **bit-identical**. Changed files:
`clus/src/Facade_Grouping.cxx` (kd2d) + `clus/inc/WireCellClus/Facade_Grouping.h` (cache member).

**Measured result** (single-event, local imaging; total = sum over all 3 MABC nodes):

| event | pre-fix | post-fix | speedup | ExamineBundles pre→post |
|--:|--:|--:|--:|--:|
| 187650 | 56.7 s | 26.5 s | **2.14×** | 39.1 → 17.8 s |
| 59789  | 55.2 s | 27.6 s | **2.00×** | 29.7 → 17.9 s |
| 138824 | 47.5 s | 25.1 s | **1.89×** | 29.7 → 18.0 s |
| 61637  | 4.0 s  | 2.9 s  | 1.38×    | 2.3 → 1.6 s |
| 141530 | 6.1 s  | 5.3 s  | 1.15×    | 2.9 → 2.2 s |

**Output-identical: verified** — pre-fix vs post-fix `mabc.zip` is **byte-for-byte identical**
across all Bee members on 59789, 61637 and 187650 (including the worst event). The tail events
(the ones that matter) get ~2×; the already-fast events less, because they spent proportionally
little time in the format.

**Post-fix profile** (re-run of 187650 under gperftools): `String::format` is **gone** from the
top; `kd2d` dropped **47.3% → 8.9%** (now just the legitimate `scoped_view` lookup); total CPU
samples fell ~44% (17434 → ~9840), consistent with the ~2× wall speedup. The new top clustering
costs are genuine work — `vhough_transform` (now ~31%, the boost::histogram Hough in
Separate/Protect) and `nanoflann` kd-tree neighbour searches (~14%) — which are the next targets.

## 3. Next hotspot after the kd2d fix: `vhough_transform` (post-fix profile of 187650, 2026-06-05)

Re-profiling 187650 **after** the kd2d fix (`cpu_postfix.prof`; 9844 samples) shows the new
clustering shape — ExamineBundles is still #1 (61%), and inside its `connect_graph_relaxed`
the cost is now genuine geometry, dominated by one item:

```
ClusteringExamineBundles               61.0%
  connect_graph_relaxed                58.6%
    vhough_transform / hough_transform 30.9%   <-- #1 remaining: per-candidate Hough direction
        std::max_element               20.1%   (65% of vhough)  <-- scans ALL histogram bins
        boost::histogram::fill         19.5% (of vhough)
    get_closest_points (good-point)    27.5%
        nanoflann kd-tree search       14.3%   <-- #2 remaining: genuine radius queries
ProtectOverclustering                   9.5%
ClusteringDeghost                       4.6%
```

**Root cause of `vhough_transform`** (`clus/src/Facade_Cluster.cxx:1609-1680`,
`hough_transform`): it builds a **180 × 360 = 64 800-bin** `boost::histogram` with the
**default `unlimited_storage`** (an adaptive type-erased *variant* per bin), fills only the
handful of points within the search radius, then finds the peak with

```cpp
auto hist = bh::make_histogram(bh::axis::regular<>(180,...), bh::axis::regular<>(360,...));
... // fill ~few hundred points
auto indexed = bh::indexed(hist);
auto it = std::max_element(indexed.begin(), indexed.end());   // scans ALL 64 800 bins
```

So the histogram is **sparse to fill but the peak-scan is full-grid**, and every one of the
64 800 `max_element` comparisons dispatches through `unlimited_storage::buffer_type::visit` (a
std::variant visit — 31% self time). `max_element` is **65% of `vhough_transform`**. This is
called once per relaxed-bridging candidate inside `connect_graph_relaxed`.

### Fix (IMPLEMENTED 2026-06-05) — dense scalar Hough-histogram storage

`hough_transform` now builds the histogram with **`bh::dense_storage<double>`**
(`bh::make_histogram_with(bh::dense_storage<double>(), …)`) instead of the default
`unlimited_storage`.  Each bin is a plain `double`, so the `std::max_element` comparisons are
scalar compares instead of `std::variant` visits.  The axes, fill order, `bh::indexed` +
`std::max_element` peak search and its tie-break are **unchanged**, so the selected peak bin is
bit-identical.  Changed: `clus/src/Facade_Cluster.cxx` (`hough_transform`).

**Measured result** (single-event, on top of the kd2d fix; total = sum over 3 MABC nodes):

| event | original | after kd2d fix | **after vhough fix** | total vs original |
|--:|--:|--:|--:|--:|
| 187650 | 56.7 s | 26.5 s | **22.1 s** | **2.57×** |
| 59789  | 55.2 s | 27.6 s | **23.2 s** | **2.38×** |
| 138824 | 47.5 s | 25.1 s | **20.1 s** | **2.36×** |
| 61637  | 4.0 s  | 2.9 s  | **2.7 s**  | 1.48× |
| 141530 | 6.1 s  | 5.3 s  | **4.9 s**  | 1.24× |

**Output-identical: verified** — post-fix `mabc.zip` is **byte-identical (CRC-32) to the
original pre-any-fix baseline** on 59789, 61637 and 187650 (so the *combined* kd2d + vhough
chain preserves clustering output exactly).

**Post-fix profile** (187650, `cpu_vhough.prof`, 6904 samples): `buffer_type::visit` (the
variant dispatch) is gone (612 → 7 samples); `vhough_transform` dropped 30.9% → 20.9%.

### Remaining hotspots / next targets

1. **(Now #1, ~18%) `nanoflann` kd-tree radius searches** in the good-point tests
   (`get_closest_points`/`test_good_point`/`is_good_point` under `connect_graph_relaxed`) —
   genuine geometry. Algorithmic only (fewer/cheaper queries); no pure-overhead win left here.
2. **(~17% combined) deeper `vhough` rewrite** — two residual costs remain because the grid is
   still dense-but-sparsely-filled: a **64 800-double zero-fill on every call**
   (`std::__fill_a1`, ~8.5%) and the **full-grid `max_element` scan** (~8%). Both vanish if we
   **track the touched bins during fill and take the max over only those** (O(#points) instead
   of O(64 800)) — but this must reproduce `std::max_element`'s tie-break (smallest grid index
   among equal maxima) to stay byte-identical, so it is a more careful change than the storage
   swap. Largest remaining vhough lever.
3. **(Optional ~5–10%) Cache the `kd2d_t&` reference** by `(apa,face,pind)` to also skip the
   residual `scoped_view` hash-lookup — only if verified the scoped view is not invalidated
   between calls within a pass (otherwise keep the current safe Scope-only memo).

> Validate any of these the same way: **byte-identical `mabc.zip`** (CRC-32) on
> 59789 / 61637 / 187650 vs the original baseline, then re-pull per-event timing.

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

> **⚠ SUPERSEDED (2026-06-05) — foreign-input baseline.** The numbers in this section were
> measured on yuhw's foreign larsoft `*-active.npz` (pre-`fill_wrap_points`-fix, auto_mask-on)
> and are kept only for history/comparison. Under our own local imaging the per-event times
> are 5–150× smaller and the ranking changes — see **§1** above. In particular the
> "worst event, 141530 = 952 s" is **6.1 s** under local imaging, and ExamineBundles cost is
> topological, not count-driven. Use §1, not this section, for optimization decisions.

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

> **⚠ Partially revised by §1 (2026-06-05).** Direction #1 below (spatial pre-filter on the
> O(n²) pair loop) is the wrong target under local imaging — final `nclusters` ≤ 30, so that
> loop is negligible. The cost is the per-cluster `connected_blobs` (#2), and it scales with
> blob-graph *edge density*, not raw blob count (see §1's count tables). Use §1's strategy.

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
