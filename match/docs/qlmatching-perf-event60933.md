# QLMatching per-phase profile — where the 31 s on event 60933 goes

**Goal.** The §17 timing study (`chisquare_flags_comparison.md`) measured QLMatching at
~5 s median / 31 s max per event but only logged the **whole-event total**. This drills
into *which phase* of `QLMatching::operator()` spends the time, using the busiest event
(60933, lan-reco2 file 1) and a median event (58667) for the scaling comparison.

**Verdict (surprising):** the cost is **not** `build_bundles`/SemiAnalyticalModel and
**not** the two-round LASSO solve (both small). **91 % of the wall is `cull_cross_tpc`**
— the SBND cross-TPC cathode-crossing pre-fit cull (`xtpc_flag=true`) — and inside it,
the brute-force closest-approach loop whose per-point accessor `Cluster::point3d()`
re-resolves a scoped k-d-tree view (hash + hashtable lookup) on **every** point access.

## How it was measured

- **Per-phase markers** added to `match/src/QLMatching.cxx` (logging only, outputs
  bit-identical): `std::chrono` splits emitted at debug under the `QLtiming` tag in
  `operator()`, `run_one_apa_prefit`, `build_bundles`, `run_one_apa_fit`,
  `fit_round1/2`, and `cull_cross_tpc` (with pair / point-pair counters).
- **gperftools** CPU sampling (`libprofiler.so.0`, `CPUPROFILE_FREQUENCY=1000`) for
  leaf-level confirmation; `libWireCellMatch.so` is unstripped with debug info.
- **Input:** the user's own multi-3view-imaged single-event clusters staged by
  `/home/xqian/tmp/lanrepro2/f1/mkev.sh <ident>`, matched with
  `wct-clus-matching-standalone.jsonnet` (`reality=data joint=true pmt_nl=true`).
  No re-imaging (avoids the `clus_all_apa` heisenbug); deterministic ~30 s reproduces
  the documented 31 s (single-event vs the in-file run).

## Per-phase breakdown

`QLtiming operator:` line (ms), joint two-TPC matching (two `ApaRun`s → the cross-TPC
cull fires):

| phase | 60933 (busy) | 58667 (median) | ratio | note |
|---|--:|--:|--:|---|
| prefit (both APAs) | 2537 | 1013 | 2.5× | of which `build_bundles` 2448 / 953 |
| └ build_bundles `vis_loop` | 2380 | 910 | 2.6× | SemiAnalyticalModel — the *assumed* hotspot, only ~8 % |
| **xtpc_cull** | **26949** | **5713** | **4.7×** | `cull_cross_tpc` |
| fit (both LASSO rounds) | 1.4 | 1.5 | — | **negligible** (X is ~60×60) |
| output build | 125 | 94 | — | tensor/merge tail |
| **total (`took`)** | **29613** | **6822** | **4.3×** | |

`QLtiming cull_cross_tpc:` detail:

| | 60933 | 58667 |
|---|--:|--:|
| candidate bundles (both APAs) | 70 | 72 |
| APA0×APA1 pairs evaluated | 1104 | 1295 |
| time-coincident pairs (→ `xtpc_pair_consistent`) | 70 | 121 |
| **point-pairs (Σ npts0·npts1)** | **2.7×10⁸** | **5.4×10⁷** |
| pairing-loop wall | 26949 ms | 5713 ms |
| **wall / point-pair** | **~100 ns** | **~106 ns** |

**Reading it:** the cull wall is `point_pairs × ~100 ns`. The per-event time variance is
almost entirely the cross-TPC cull's point-pair count (5× more point-pairs on 60933 →
4.7× more cull time); `build_bundles` scales gentler (2.5×) and LASSO is flat. So the
cull dominates at *every* busyness level (84–91 %), not just the tail.

## Why ~100 ns per point-pair (not ~1 ns)

`xtpc_pair_consistent` (`QLMatching.cxx:1925`) finds the closest approach by a
brute-force double loop over **all** 3D points of both clusters:

```cpp
for (int a = 0; a < m0.c->npoints(); ++a) {
    const geo_point_t ra = m0.c->point3d(a);          // <-- scope re-resolution
    for (int b = 0; b < m1.c->npoints(); ++b) {
        const geo_point_t rb = m1.c->point3d(b);      // <-- scope re-resolution
        ... d2 = dx*dx+dy*dy+dz*dz ...
    }
}
```

`Cluster::point3d(i)` = `kd3d().point3d(i)` (`Facade_Cluster.cxx:832`), and `kd3d()`
resolves the scoped k-d view by **hashing a `Scope` (string keys `{"x","y","z"}`) and a
hashtable lookup every call**. The arithmetic is ~6 subtractions; the accessor is ~100×
that. gperftools confirms — focused on `xtpc_pair_consistent` (70.6 % of all samples):

```
Scope::hash                     51.4% cum   boost::hash_mix      28.1%
_Map_base::operator[]           58.9% cum   __memcmp_evex_movbe  11.8%
_Hashtable::_M_find_*           ~32% cum    boost::hash_range    18.3%
```

i.e. essentially all of the cull time is `Scope` hashing + hashtable lookups behind
`point3d()`, **not** distance math, PCA, `vhough_transform`, or the LASSO solver. (Note:
the inner `rb = point3d(b)` is also fully re-resolved for every `a`, so the inner cluster
is re-hashed npts0 times.)

## Does the closest-point search use the k-d tree? No.

The hand-rolled double loop **bypasses** the k-d tree for the *search* while paying its
*per-access* cost (the `point3d` → `kd3d()` Scope re-resolution above). `Facade::Cluster`
already provides the correct primitive:

- **`Cluster::get_closest_points(const Cluster& other)`** (`Facade_Cluster.cxx:1357`)
  returns `(idx1, idx2, distance)` and **is k-d-backed**: ~20 probe points (stride
  `npoints()/20`) → `other.kd3d().knn(5, p1)` → alternating closest-point refinement
  (≤3 iters, early-out at 0.5 cm). ~O(log N) per cluster-pair, not O(npts0·npts1).
- Related: `get_closest_point_index`, `get_closest_dis`, `kd_knn`, `kd_radius`.

The cull reinvents `get_closest_points`, the slow way.

## Result — fix (1) implemented and measured (bit-identical)

The hoist below was applied to `xtpc_pair_consistent` (`QLMatching.cxx:1924`). Same
single-event reruns, `point_pairs` unchanged (same points scanned), cross-TPC decisions
identical (same consistent pairs, same `d`/angles/scenario on both events):

| event | `xtpc_cull` before → after | total before → after | cull / point-pair |
|---|---|---|---|
| 60933 (busy)   | 26949 ms → **300 ms** (~90×) | 29613 ms → **2956 ms** (~10×) | 100 ns → ~1.1 ns |
| 58667 (median) | 5713 ms → **83 ms** (~69×)   | 6822 ms → **1171 ms** (~5.8×) | 106 ns → ~1.5 ns |

After the fix the cull is no longer the bottleneck; `build_bundles` (SemiAnalyticalModel,
~2.5 s on 60933) becomes the dominant phase — the next target if more is wanted.

## Optimization headroom

Two distinct fixes — keep them separate; the first is the headline (now done):

1. **Hoist the scope resolution out of the loop — dominant win, bit-identical, no toggle.**
   *(Implemented — see Result above.)*
   The arithmetic and the T0 offset are already cheap; the slow part is fetching the *raw*
   coordinates. `m1.c->point3d(b)` = `m1.c->kd3d().point3d(b)`, and `kd3d()` re-resolves the
   scoped view (Scope hash + hashtable lookup) on **every** call — and `rb = point3d(b)` is
   re-fetched for every `(a,b)`, so the scope is hashed `npts0·npts1` times per pair
   (2.7×10⁸ on 60933). The raw coordinates are **T0-independent**, so grab the flat raw
   arrays **once** per cluster and index them; the per-point T0 shift (`+off/+dy/+dz`) stays
   a runtime **scalar** add exactly as today (the distance only depends on the *relative*
   offset, which those scalar adds already encode):

   ```cpp
   const auto& P0 = m0.c->points();   // resolve scope ONCE -> {xs,ys,zs}
   const auto& P1 = m1.c->points();   // resolve scope ONCE
   for (int a=0; a<n0; ++a) {
       const double ax = P0[0][a] + m0.off, ay = P0[1][a] + m0.dy, az = P0[2][a] + m0.dz;
       for (int b=0; b<n1; ++b) {
           const double dx = ax - (P1[0][b] + m1.off); /* + dy,dz */ ...
       }
   }
   ```

   Same point indices, iteration order, and first-min-wins tie-break → **bit-identical**.
   ~100 ns → ~1–2 ns/pair: 60933's 26949 ms → **~0.3–0.5 s**. **Nothing T0-dependent is
   precomputed** — only the scope resolution is hoisted; the shift remains per-pair scalar.
   *Caveat:* keep the min-reduction scalar — don't `#pragma`/auto-vectorize it, or
   equal-distance ties could pick a different index (physically negligible but not
   bit-identical).

2. **k-d nearest-neighbour — further speedup, but CHANGES results (needs toggle +
   validation).** `Cluster::get_closest_points(other)` (k-d `knn` + alternating refine)
   is ~O(log N) and would cut the ~0.5 s above to <0.01 s, **but it is a heuristic**: ~20
   probes + `knn(5)` + ≤3 refinement iters find a *near*-global closest approach, not the
   guaranteed global minimum the brute-force loop returns. A different closest pair →
   different `d` and `vhough` angles → a scenario1/2 verdict can flip → matching output
   changes. Also it works on **raw** `point3d` (its refinement queries the raw k-d trees),
   so the cull's offset frame (`m.off/dy/dz`) is not a clean drop-in — making it
   offset-correct means threading offsets through refinement or building a k-d on shifted
   coords. Treat as an opt-in variant per [[feedback_toggleable_behavior_changes]]
   (default OFF, MC/data validate), not a free swap. After fix (1) the cull is no longer
   the bottleneck, so this is optional.
3. **Memoize / pre-filter (cheap, exact).** Multiple bundles share a main cluster; the
   same (cluster0, cluster1) closest-approach is recomputed per bundle pair — cache it
   keyed on the two clusters + their flash-T0 offsets. And reject obviously-far pairs by
   AABB / centroid distance before the per-point scan (a real crosser's two halves sit at
   d ≤ ~264 cm). Both preserve the exact result.

`build_bundles`' `vis_loop` (SemiAnalyticalModel, ~2.4 s on 60933) is the distant second
and a separate target if the cull is fixed.

## Correction to §17.1

§17.1 attributed the per-event cost to "per-cluster geometry (PCA, `vhough`,
`get_extreme_wcps`) and the two-round LASSO." That is **incorrect**: the LASSO solve is
~1.4 ms and `vhough`/PCA are not the hot leaf. The dominant cost is the cross-TPC cull's
`point3d()` closest-approach loop (this doc). The 5 s→31 s tail is driven by the cull's
point-pair count, not by PCA/vhough/LASSO scaling.

## Full 150-event reprocess (post-fix) and the next hotspot

Re-ran **all 150 lan-reco2 data events** (3 files × 50) with the cull fix: own per-event
imaging (stored in `sbnd_xin/work/evt<ID>/` for reuse) + per-event clustering/matching
(`work/ql_evt<ID>/`), one event per process, 32-way parallel. **149 matched; evt 62561 is
empty (no clusters/flashes). Zero crashes** — the per-event isolation avoided the
`clus_all_apa` heisenbug ([[project_local_imaging_clus_heisenbug]]). Timings recovered
from the per-event logs (`recover_timings.py`).

**QLMatching per-event (ms), 149 events — the cull fix holds across the whole sample:**

| phase | median | p90 | p99 | max |
|---|--:|--:|--:|--:|
| `ql_total` | 1146 | 2952 | 4363 | **5630** |
| prefit | 996 | 2581 | 4056 | 5370 |
| └ `build_bundles` vis | 887 | 2362 | 3761 | 5060 |
| **`xtpc_cull` (post-fix)** | **56** | 188 | 346 | **406** |
| fit (LASSO) | 1.4 | — | — | 5.5 |
| output | 88 | — | — | 206 |

QLMatching max dropped from **31 s → 5.6 s**; `xtpc_cull` max from **27 s → 0.4 s**.

**Next hotspot — two levels:**

1. **Within QLMatching:** `build_bundles` / `SemiAnalyticalModel` is now dominant —
   median 887 ms, max 5060 ms (~77 % of `ql_total`). It's the per-point
   `detectedDirect/ReflectedVisibilities` over all opdets (`build_bundles:927-964`). Headroom:
   the visibility loop redundantly re-evaluates per-point solid-angle/Gaisser-Hillas; a
   first cut is to skip masked opdets earlier and cache per-opdet geometry, or coarsen the
   per-point sampling for large blobs.

2. **In the whole per-event chain (QLMatching is no longer significant):** per-event wall
   (timestamps) is now dominated by

   | stage | median | max | total (149 evt) |
   |---|--:|--:|--:|
   | imaging | 9.9 s | 43 s | **1589 s** |
   | match-graph (clustering + QLMatching + IO) | 5.2 s | 40 s | 971 s |
   | └ clustering (≈ match-graph − QLMatching) | 4.1 s | 39 s | ~754 s |
   | └ QLMatching | 1.1 s | 5.6 s | 217 s |

   The clustering tail is **`ClusteringExamineBundles`** (the O(n²) pair loop,
   `clus/src/clustering_examine_bundles.cxx`): on evt 59789 it was 10.7 s (`:all`) + 11.8 s
   (`:apa0-0`), ~95 % of that event's clustering. So the real next targets are **imaging**
   (largest aggregate) and **`ClusteringExamineBundles`** — see
   `clus/docs/clustering-timing-profile.md` and [[project_clustering_examinebundles_hotspot]],
   not QLMatching.
