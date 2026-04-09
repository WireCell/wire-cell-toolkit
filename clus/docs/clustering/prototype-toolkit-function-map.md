# Prototype → Toolkit Clustering Function Map

This document maps clustering functions from the prototype (`prototype_base/2dtoy/src/`) to
their corresponding toolkit components used in the MABC (MultiAlgBlobClustering) pipeline,
as configured in `pdhd/wct-clustering.jsonnet` via `clus/test/test-porting/pdhd/clus.jsonnet`
and the common factory in `cfg/pgrapher/common/clus.jsonnet`.

## Pipeline Structure Overview

The PDHD clustering runs in three nested stages:

1. **Per-Face** — one MABC instance per APA face (2 faces × N APAs)
2. **Per-APA** — one MABC instance per APA, runs after faces are merged
3. **All-APA** — one MABC instance for the full detector, after all APAs are merged

---

## Per-Face Pipeline (`cm_pipeline` in `clus_per_face`)

| `clus.jsonnet` method | Toolkit Component | Prototype Function | Prototype File |
|---|---|---|---|
| `cm.pointed()` | `ClusteringPointed` | *(no direct equivalent — new in toolkit)* | `clus/src/clustering_pointed.cxx` |
| `cm.live_dead(dead_live_overlap_offset=2)` | `ClusteringLiveDead` | `Clustering_live_dead()` | `ToyClustering_dead_live.h` |
| `cm.extend(flag=4, length_cut=60cm, num_try=0, length_2_cut=15cm, num_dead_try=1)` | `ClusteringExtend` | `Clustering_extend(flag=4, 60cm, 0, 15cm, 1)` | `ToyClustering_extend.h` |
| `cm.regular(name="-one", length_cut=60cm, flag_enable_extend=false)` | `ClusteringRegular` | `Clustering_regular(60cm, false)` | `ToyClustering_reg.h` |
| `cm.regular(name="_two", length_cut=30cm, flag_enable_extend=true)` | `ClusteringRegular` | `Clustering_regular(30cm, true)` | `ToyClustering_reg.h` |
| `cm.parallel_prolong(length_cut=35cm)` | `ClusteringParallelProlong` | `Clustering_parallel_prolong(35cm)` | `ToyClustering_para_prol.h` |
| `cm.close(length_cut=1.2cm)` | `ClusteringClose` | `Clustering_close(1.2cm)` | `ToyClustering_close.h` |
| `cm.extend_loop(num_try=3)` | `ClusteringExtendLoop` | Inline `for` loop (num_try=3) calling `Clustering_extend()` variants | `ToyClustering.cxx` (lines ~293–330) |
| `cm.separate(use_ctpc=true)` | `ClusteringSeparate` | `Clustering_separate()` | `ToyClustering_separate.h` |
| `cm.connect1()` | `ClusteringConnect1` | `Clustering_connect1()` | `ToyClustering_connect.h` |

**Note:** `cm.ctpointcloud()` (`ClusteringCTPointcloud`) maps to `Clustering_CTPointCloud()` in
`ToyClustering_CTPointCloud.h`, but is **commented out** in both the toolkit pipeline and
`ToyClustering.cxx` (line 252).

---

## Per-APA Pipeline (`cm_pipeline` in `clus_per_apa`)

| `clus.jsonnet` method | Toolkit Component | Prototype Function | Prototype File |
|---|---|---|---|
| `cm.deghost()` | `ClusteringDeghost` | `Clustering_deghost()` | `ToyClustering_deghost.h` |
| `cm.protect_overclustering()` | `ClusteringProtectOverclustering` | `Clustering_protect_overclustering()` | `ToyClustering_protect_overclustering.h` |

---

## All-APA Pipeline (`cm_pipeline` in `clus_all_apa`)

| `clus.jsonnet` method | Toolkit Component | Prototype Function | Prototype File |
|---|---|---|---|
| `cm_old.switch_scope()` | `ClusteringSwitchScope` | *(no direct equivalent — handles T0 coordinate correction)* | `clus/src/clustering_switch_scope.cxx` |
| `cm.extend(flag=4, length_cut=60cm, ...)` | `ClusteringExtend` | `Clustering_extend()` | `ToyClustering_extend.h` |
| `cm.regular(name="1", length_cut=60cm, flag_enable_extend=false)` | `ClusteringRegular` | `Clustering_regular(60cm, false)` | `ToyClustering_reg.h` |
| `cm.regular(name="2", length_cut=30cm, flag_enable_extend=true)` | `ClusteringRegular` | `Clustering_regular(30cm, true)` | `ToyClustering_reg.h` |
| `cm.parallel_prolong(length_cut=35cm)` | `ClusteringParallelProlong` | `Clustering_parallel_prolong(35cm)` | `ToyClustering_para_prol.h` |
| `cm.close(length_cut=1.2cm)` | `ClusteringClose` | `Clustering_close(1.2cm)` | `ToyClustering_close.h` |
| `cm.extend_loop(num_try=3)` | `ClusteringExtendLoop` | Inline `for` loop (num_try=3) calling `Clustering_extend()` variants | `ToyClustering.cxx` (lines ~506–545) |
| `cm.separate(use_ctpc=true)` | `ClusteringSeparate` | `Clustering_separate()` | `ToyClustering_separate.h` |
| `cm.neutrino()` | `ClusteringNeutrino` | `Clustering_neutrino()` | `ToyClustering_neutrino.h` |
| `cm.isolated()` | `ClusteringIsolated` | `Clustering_isolated()` | `ToyClustering_isolated.h` |
| `cm.examine_bundles()` | `ClusteringExamineBundles` | `ExamineBundles()` / `ExamineBundle()` | `ExamineBundles.cxx` |
| `cm.retile(retiler=retiler)` | `ClusteringRetile` | *(no equivalent — new in toolkit, re-tiles blobs across APA boundaries)* | `clus/src/clustering_retile.cxx` |

**Note:** `Clustering_examine_x_boundary()` (defined in `ToyClustering_neutrino.h`, used in
prototype `Clustering_jump_gap_cosmics`) is exposed in the toolkit as `cm.examine_x_boundary()`
→ `ClusteringExamineXBoundary`, but is **commented out** in the all-APA `cm_pipeline`.

---

## Complete Prototype Function List (for reference)

The following `Clustering_*` functions appear in `ToyClustering.cxx` and the included headers:

| Prototype Function | Header File |
|---|---|
| `Clustering_live_dead()` | `ToyClustering_dead_live.h` |
| `Clustering_regular()` | `ToyClustering_reg.h` |
| `Clustering_parallel_prolong()` | `ToyClustering_para_prol.h` |
| `Clustering_close()` | `ToyClustering_close.h` |
| `Clustering_extend()` | `ToyClustering_extend.h` |
| `Clustering_separate()` | `ToyClustering_separate.h` |
| `Clustering_deghost()` | `ToyClustering_deghost.h` |
| `Clustering_connect1()` | `ToyClustering_connect.h` |
| `Clustering_protect_overclustering()` | `ToyClustering_protect_overclustering.h` |
| `Clustering_neutrino()` | `ToyClustering_neutrino.h` |
| `Clustering_examine_x_boundary()` | `ToyClustering_neutrino.h` |
| `Clustering_isolated()` | `ToyClustering_isolated.h` |
| `Clustering_CTPointCloud()` | `ToyClustering_CTPointCloud.h` |
| `ExamineBundles()` / `ExamineBundle()` | `ExamineBundles.cxx` |

---

## Functions to Examine

The following toolkit–prototype pairs are the primary targets for detailed comparison:

1. **`ClusteringLiveDead` ↔ `Clustering_live_dead`** — live/dead channel overlap merging
2. **`ClusteringExtend` ↔ `Clustering_extend`** — cluster extension by flag mode
3. **`ClusteringRegular` ↔ `Clustering_regular`** — regular length-based merging
4. **`ClusteringParallelProlong` ↔ `Clustering_parallel_prolong`** — parallel/prolongation merging
5. **`ClusteringClose` ↔ `Clustering_close`** — close-proximity merging
6. **`ClusteringExtendLoop` ↔ inline loop** — looped extension (wraps `Clustering_extend`)
7. **`ClusteringSeparate` ↔ `Clustering_separate`** — cluster separation
8. **`ClusteringConnect1` ↔ `Clustering_connect1`** — connection via global point cloud
9. **`ClusteringDeghost` ↔ `Clustering_deghost`** — ghost removal
10. **`ClusteringProtectOverclustering` ↔ `Clustering_protect_overclustering`** — overclustering protection
11. **`ClusteringNeutrino` ↔ `Clustering_neutrino`** — neutrino-specific clustering
12. **`ClusteringIsolated` ↔ `Clustering_isolated`** — isolated cluster handling
13. **`ClusteringExamineBundles` ↔ `ExamineBundles`** — bundle examination
14. **`ClusteringPointed`** — new in toolkit (no prototype counterpart)
15. **`ClusteringRetile`** — new in toolkit (no prototype counterpart)
16. **`ClusteringSwitchScope`** — new in toolkit (T0 coordinate correction, no prototype counterpart)
