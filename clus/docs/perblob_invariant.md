# The "perblob" row-order invariant

**Invariant.** A cluster's node-local `"perblob"` point-cloud Dataset (arrays
`isolated`, `assoc_cluster_id`, `assoc_cluster_main`, `real_cluster_id`,
`real_cluster_main`, ...) has exactly `nchildren()` rows and **row *i*
describes blob child *i* in children order**.  Every consumer assumes this:
`TaggerCheckTGM` (`main_component_mode="real"`) maps KD-nearest blob indices
into `real_cluster_main` rows, `ClusteringUnmergeBundle` splits clusters on
the `assoc_*`/`isolated` rows, the Bee dump colors by `real_cluster_id`, and
TensorDM serialization round-trips rows positionally.

**Repro / audit sweep** (how the table below was produced, 2026-07-25):

```
grep -rn "separate(\|take_children\|->merge(\|sort_children" clus/src match/src img/src
grep -rn "perblob" clus/src match/src
```

## Why the invariant is fragile

The tree primitives do **not** maintain it:

- `Grouping::separate(cluster, cc)` — children with `cc[i] < 0` stay (original
  relative order); `cc[i] >= 0` children move to new clusters keyed by group
  id.  The survivor's `local_pcs` are untouched (stale length AND order); new
  split clusters get no perblob arrays; with `remove=true` the original's
  arrays are destroyed with it.  `from()` copies scalars/flags only.
- `FacadeParent::merge()` / `take_children()` — appends the other clusters'
  children at the END of the target's list (`std::map` iteration = ascending
  group id).  The target's `local_pcs` are untouched.

So a **separate-then-merge-back round trip preserves the array length while
permuting the children**: kept (`cc<0`) blobs stay in place, split-off groups
are re-appended at the end in ascending-gid order.  Every length check passes;
row *i* silently stops describing child *i*.  This was the doc 52 §12
QLMatching decompose/recompose corruption (wcp-porting-validation
`sbnd/sbnd_xin/docs/52_isolated-grouping-fix-design.md`), and the 2026-07-25
audit (§13 there) found the identical latent pattern in `RetileCluster::mutate`
and twice in `clustering_neutrino`.

## The primitives enforce the invariant (doc 52 §13, option 2)

Since 2026-07-25 the Clus-level primitives maintain the invariant
themselves:

- **`Grouping::separate()`** snapshots the cluster's `"perblob"` Dataset and
  **carves** it across the survivor (kept `cc<0` rows, original order; the
  entry is erased when no rows remain) and every split (its own rows).  A
  Dataset whose rows already disagree with the children is **dropped with a
  warning** — a loud absence instead of a silent lie.
- **`Grouping::merge()`** (new overloads that deliberately HIDE the
  `NaryTree::FacadeParent` merges) **concatenates** the parts' Datasets onto
  the target's in the same order the children are adopted.  When a truthful
  concatenation is impossible (stale rows, a child-bearing part without a
  Dataset while others carry one, mismatched key sets) the target's Dataset
  is dropped with a warning.  The returned cc array is exactly the base
  merge's.

So any separate/merge sequence through the **Grouping** methods keeps
row *i* == child *i* automatically.  Clusters without the Dataset (the whole
pre-`isolated` clustering stage) pay one `local_pcs` map lookup per
participant — measured zero wall/RSS drift on the SBND 30-event, PDHD/PDVD
abtest, and uBooNE 35-event manifests.

## Rules for any code that mutates a cluster's child set

1. **Prefer `merge_clusters()`** (`ClusteringFuncs.cxx`) for graph-driven
   merging — it carries the perblob provenance keyed by **blob pointer**,
   immune to member-order effects (doc 52 Stage 2), and REBASES the id
   arrays across members (raw concatenation must not be used there: member
   ids collide).
2. **Split / merge-back through the Grouping methods** —
   `Grouping::separate()` / `Grouping::merge()` — and the Dataset follows
   the blobs automatically.  Do NOT re-subset the Dataset afterwards: it is
   already aligned, and a second permutation would scramble it.  (This
   superseded the inline realigns that briefly lived in QLMatching
   recompose, `retile_cluster.cxx` and `clustering_neutrino.cxx`.)
3. **Raw NaryTree primitives** (`node()->separate`, `take_children`,
   `adopt_children`) bypass the enforcement.  If you must use them on a
   perblob-carrying cluster, restore the invariant with
   `realign_perblob_after_regroup(cluster, cc)` (`ClusteringFuncs.h`,
   implemented in `prov_check.cxx`) for a length-preserving round trip, or
   carve by hand (`clustering_switch_scope.cxx`'s carve lambda,
   `ClusteringUnmergeBundle.cxx`).
4. **Intentional drops stay explicit**: `clustering_switch_scope.cxx` erases
   the auto-carved Dataset on each part and re-attaches only its carry list
   (it has always deliberately dropped `isolated`); `ClusteringUnmergeBundle`
   and `ClusteringRecoveringBundle` keep their own carves as defense in
   depth (same values as the auto-carve).
5. **Never** iterate/emit these arrays through a child-set mutation you did
   not account for — and remember `put_pcarray()` on an existing key
   `assign()`s with a new shape, bypassing `Dataset::add`'s major-size
   check, so a wrong write is silent.

## Regression tripwire

`WCT_PROV_CHECK=1` enables `check_perblob_provenance()` at every
MultiAlgBlobClustering boundary (load / after each visitor / save /
serialize-deserialize round trip): per-key row counts vs `nchildren`, plus
spatial compactness (rmax > 25 cm) of each `assoc_cluster_main == 0` group.
Keep it on in A/B harness runs; it is a fast no-op otherwise.

## Audit table (2026-07-25, doc 52 §13)

Classes: (a) silent length-preserving permutation; (b) stale-length arrays;
(c) arrays lost with a destroyed/emptied cluster.

| Site | Pattern | Class | Status |
|---|---|---|---|
| `match/src/QLMatching.cxx` decompose/recompose | separate + merge back | (a) | ENFORCED by the primitives (§13.3); the interim inline realign (doc 52 §12.4) removed, byte-identically |
| `clus/src/retile_cluster.cxx` `RetileCluster::mutate` | separate + merge back, fresh `isolated` only | (a) | ENFORCED by the primitives; reachable only via `cm.retile`, commented out in all configs |
| `clus/src/clustering_neutrino.cxx` cluster1/cluster2 round trips | separate + total take-back | (a) | ENFORCED — loops replaced by `Grouping::merge`; dormant anyway (runs before `isolated` writes perblob in every config) |
| `clus/src/ClusteringRecoveringBundle.cxx` | separate, pieces leave | (b)+(c) | FIXED (carve) + primitive auto-carve. CORRECTION: it IS instantiated — qlport `uboone-mabc.jsonnet:1247` (uBooNE chain); the original "not instantiated" claim only grepped `cfg/pgrapher/experiment/`. uBooNE gate: doc 52 §13.3 |
| `clus/src/clustering_switch_scope.cxx` | separate remove=true | safe | carves `assoc_*`/`real_*` per part (drops `isolated`, rewritten downstream) |
| `clus/src/ClusteringUnmergeBundle.cxx` | separate remove=false | safe | subsets whole Dataset per part |
| 12 graph-merge visitors via `merge_clusters()` | merge | safe | pointer-keyed carry (doc 52 Stage 2) |
| `clus/src/CreateSteinerGraph.cxx` + `ImproveCluster_1/2` | retile a scratch COPY | safe | source cluster's child set never mutated |
| `clus/src/clustering_separate.cxx` (many), `examine_x_boundary`, `protect_overclustering`, `Separate_1` | absorb/pool/split | (b)/(c) if arrays present | dormant: all run before `isolated` writes perblob; fragments intentionally array-less |
| `clus/src/GroupingHelper.cxx` `process_groupings_helper` | separate both groupings | (b)+(c) | dead code (only caller commented out) |
| `clus/src/clustering_test.cxx` | separate remove=true | (c) pattern | test-only, config refs commented out |
| `util` `NaryTree::sort_children` | reorder | would be (a) | zero callers |

Side finding (not fixed here): `clustering_protect_overclustering.cxx` — a
blob with zero points in the default-scope kd3d keeps groupid −1 and is
destroyed with the original cluster (silent blob loss, independent of perblob).
