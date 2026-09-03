# `connect_graph_relaxed_strict` perf — round 1

**Scope.** Four changes to the strict connect-graph family and its two helpers.
Three are byte-identical waste removal with no knob; the fourth is a new graph
flavor carrying doc 78's busy-gated lazy walk, which no job selects by default.
No production default is flipped by this round.

**Status.** Shipped, gated. The lazy flavor is available and measured but is
**not** turned on anywhere — see *The flip is the owner's call* below.

## Repro

```bash
wcbuild
./build/clus/wcdoctest-clus                       # 263/263, 2863 assertions

cd /home/xqian/toolkit-dev/wcp-porting-img/pdvd
# the census that sized the problem (kill once 'OC53CENSUS-S cluster ... ncomp=' prints)
WCT_RELAXED_EDGE_CENSUS=1 ./run_pr_evt.sh -s cens8 \
    -pipe switch_scope,flag_mains,steiner,fiducialutils,tagger_check_tgm,tagger_check_stm,tagger_check_fc,protect_bundle \
    039252 8

# item 4, lazy vs eager on the pathological cluster (both arms open every bundle)
PDVD_PR_TLA="-S protect_stm_only_bundles=false -S protect_graph_name='relaxed_strict_img_2d_rescue_long_wtrack_fast'" \
    ./run_pr_evt.sh -s r10fast  -stm 039252 8
PDVD_PR_TLA="-S protect_stm_only_bundles=false" \
    ./run_pr_evt.sh -s r10eager -stm 039252 8

# items 1+2+3 alone, controlled (same box, concurrent, eager flavor, libs pinned)
LD_LIBRARY_PATH=/home/xqian/tmp/doc25gate/lib_r10old PDVD_PR_TLA="-S protect_stm_only_bundles=false" \
    ./run_pr_evt.sh -s b12old -stm 039252 16
LD_LIBRARY_PATH=/home/xqian/tmp/doc25gate/lib_r10new PDVD_PR_TLA="-S protect_stm_only_bundles=false" \
    ./run_pr_evt.sh -s b12new -stm 039252 16

# shared-component gate (round 10 of the doc 25 series)
NEW_ARM=r10new NEW_LIB=/home/xqian/tmp/doc25gate/lib_r10new \
OLD_ARM=r10old OLD_LIB=/home/xqian/tmp/doc25gate/lib_r10old \
    QL_SUFFIX=d99fix PR_JOBS=6 ./stm/gates/shared_gate.sh arms
NEW_ARM=r10new OLD_ARM=r10old ./stm/gates/shared_gate.sh compare
```

## 1. The problem, sized

From `wcp-porting-img/pdvd/docs/25_pdvd-stopping-muon-michel-chain.md` §13.11.
`ClusteringProtectBundle` is the only consumer of the `relaxed_strict*` graph
flavors. On PDVD 039252/8 the stage took **1621 s**, and **1602 s of that
(98.8 %) was one `connected_blobs` call on one cluster**: 18876 blobs, 93584
points, and — from `WCT_RELAXED_EDGE_CENSUS=1` —

```
OC53CENSUS-S cluster nblobs=18876 npoints=93584 ncomp=509
```

509 components ⇒ 509·508/2 = **129 286 component pairs at 12.4 ms each**. RSS
is flat at 3.58 GB throughout, so it is CPU, not memory: the five `num × num`
tuple arrays are ~15 MB at num = 509.

Per pair the connector does ~220 kd queries in `get_closest_points`, two 30 cm
`vhough_transform` calls when the clouds are big and closer than 80 cm, and up
to **three** 1 cm-step CTPC path walks (closest pair, `dir1`, `dir2`) with a
three-plane `Grouping::test_good_point` per step, plus the S6 gap and S7
corridor kills.

**None of it had a fast path.** The doc 78 round-2/3 busy-gated lazy walk
(`RelaxedFastCfg`, `busy_num_threshold = 200`) — SBND production as
`eb_fast`/`po_fast` — was only ever threaded into `connect_graph_relaxed`.
The strict family never received the cfg and contained no `lazy_walk`.

For contrast, SBND's PR-stage protect runs the **identical** flavor for
31–257 ms per event; the whole difference is scope (a 200–2200 ns beam window
opens ≈ 1 bundle, PDVD's readout-wide window opens every matched bundle by
design). That scope lever is doc 25 §13.11's `stm_only_bundles`. This round is
the orthogonal one: make the work itself cheaper.

## 2. What was changed

| # | change | files | identity |
|---|---|---|---|
| 1 | hoist the closest-pair walk out of the `(j,k)` loop into `eval_pair_verdict` | `connect_graph_relaxed_strict.cxx` | pure refactor |
| 2 | port doc 78's busy-gated lazy Kruskal; new flavor `relaxed_strict_img_2d_rescue_long_wtrack_fast` | `connect_graph_relaxed_strict.cxx`, `connect_graphs.h`, `make_graphs.{h,cxx}`, `Facade_Cluster.cxx` | new name; unselected ⇒ identical by construction |
| 3 | drop the dead `knn(5)` seed query in `get_closest_points` | `Facade_Util.h` | byte-identical |
| 4 | allocation-free `Grouping::test_good_point` for the per-step walks | `Facade_Grouping.{h,cxx}`, `connect_graph_relaxed{,_strict}.cxx` | byte-identical |

### Item 1 — the refactor, kept separate on purpose

`connect_graph_relaxed.cxx` had this done in doc 78 round 1, *before* round 2
could gate anything. The strict file's walk was inline across ~215 lines,
interleaved with `max_ghost_run`, `two_d_gap_kill`, `long_gap_kill` and three
census emissions, so there was no handle to defer. The hoist is its own commit
so a gate failure could be attributed to the move rather than to the lazy mode
built on it. Doc 78 round 1's `forced_kill` early exit was deliberately **not**
ported — that is a separate optimization with its own identity question.

### Item 2 — the lazy walk, and why it is a flavor and not a config key

`connect_graph_relaxed_strict` gained an optional trailing
`const RelaxedFastCfg* fast = nullptr` (the same struct, the same 200
threshold). `nullptr` is the legacy eager path, bit-for-bit. Non-null **and**
`num > threshold` defers the closest-pair verdict and evaluates it by a
single-pass Kruskal over pairs in ascending `(distance, j, k)`: a pair is
walked only when it would bridge two union-find components. Pairs under the
3 cm direct-edge rule and pairs carrying a directional candidate are still
walked eagerly, because the direct-edge rule and MST #2's MAIN-entry recording
need their verdicts regardless of the forest.

Selected by a distinct **flavor name** rather than a component config key, the
same wiring choice doc 78 round 2 made for `relaxed_fast`. A job that does not
ask for `..._long_wtrack_fast` is byte-identical by construction, which is
what makes the SBND/uBooNE gate below meaningful without a knob-off arm.

Four registrations are needed (`find_graph` ×2, `graph_algorithms` ×2). Miss
one and the flavor silently falls back to another graph, presenting as "the
port didn't help" rather than as a wiring bug — so a debug line reports the
gate decision whenever a cfg is supplied:

```
connect_graph_relaxed_strict: nblobs=18876 npoints=93584 ncomp=509 threshold=200 lazy=true
```

### Items 3 and 4 — provably dead work

`Clus::Facade::get_closest_points` sampled ~20 seeds per pair and called
`two.kd().knn(5, p1)` per seed, then discarded **both** returned values: the
first two statements of the loop body overwrite `idx2`, and `dist` is never
read. Its only surviving effect was to run the body `min(5, npoints)` times.
That count is load-bearing — the alternating walk is stateful in `p1` — so it
is now taken directly from `kd().npoints()` and the search is skipped.
nanoflann's `KNNResultSet` truncates rather than pads, and `knn()` returns
empty on an empty cloud, so the count is reproduced exactly.

`Grouping::test_good_point` returned a freshly heap-allocated 6-element
vector, once per 1 cm step of every walked pair. **Note what this is not:** the
doc-78 entry-16 plane short-circuit that `is_good_point`/`is_good_point_wc`
received does *not* apply — `test_good_point` returns all six per-plane scores
and its callers consume every one in the `num_bad`/`num_bad1` accumulators, so
it cannot decide early. The available win is the allocation, not the third
plane's kd descent. An `int (&)[6]` overload now holds the body and the vector
form delegates to it, so the two cannot drift.

## 3. Measured

### Item 2 (PDVD 039252/8, same binary, `-stm`, every bundle opened)

| | eager (`r10eager`) | lazy (`r10fast`) |
|---|---|---|
| `ClusteringProtectBundle:pr` | 1710.99 s | **797.65 s** (**2.15×**) |
| event wall | 1838 s | **927 s** (1.98×) |
| splits | 25 → 302 | 25 → 302, all 25 lines identical |
| `mabc-pr.zip` member hash | `20dc1d01…` | `20dc1d01…` — **byte-identical** |

The gate fired exactly where intended, and the Kruskal census explains the
size of the win:

```
connect_graph_relaxed_strict lazy: num=509 pairs=129127 walked=57227 accepted=413
```

**44.3 % of pairs are still walked**, against 3.2 % (15 941 / 503 506) in doc
78's ExamineBundles monster. The reason is structural and worth recording: the
lazy gate defers the **closest-pair** verdict only. `dir1` and `dir2` keep
their own eager walks, and every pair carrying a directional candidate is also
pre-walked eagerly for MST #2. So the ceiling here is not `1/(pairs walked)` —
it is bounded by the un-gated directional term. 2.15× is close to that ceiling,
not a sign of a broken port.

### Items 1+3+4 alone (PDVD 039252/16, libs pinned, concurrent arms)

| | pre-round lib (`b12old`) | this round (`b12new`) |
|---|---|---|
| `ClusteringProtectBundle:pr` | 250.88 s | 256.39 s |
| `CreateSteinerGraph:pr` | 27.17 s | 26.71 s |
| splits | 47 → 509 | 47 → 509 |
| `mabc-pr.zip` member hash | `d9705233…` | `d9705233…` — **byte-identical** |

**No measurable wall-clock win, and the arithmetic says there should not be
one.** Item 4 removes ~19 M allocations on the heavy event (3 walks × ~50
steps × 129 k pairs); at ~20 ns per tcmalloc round trip that is ~0.4 s out of
1700 s. Item 3 removes 1 kd search of ~11 per seed inside `get_closest_points`,
which is itself a small share of the 12.4 ms per pair. Both are below the
±2–3 % arm-to-arm noise of a loaded 64-core box. They are kept because they
delete work that provably cannot affect the result — not because they were
shown to be fast. **Do not quote them as a speed-up.**

The SBND gate arms tell the same story (48 nuecc48 events, fast flavor not
selected): `ClusteringProtectBundle:pr` 9.91 s → 10.17 s, `CreateSteinerGraph`
31.99 s → 33.84 s — i.e. arm-to-arm load noise on ~200 ms/event stages, with
byte-identical outputs.

## 4. Gates

- **`./build/clus/wcdoctest-clus`**: 263/263 test cases, 2863 assertions.
- **Shared-component gate round 10**
  (`wcp-porting-img/pdvd/stm/gates/shared_gate_round10.txt`). OLD = `6365aa00`
  for this round's 9 files with the concurrent session's state untouched;
  NEW = the same tree plus this round. Both libraries built back to back
  (`lib_r10old` clus `c018c30dc71d`, no fast flavor; `lib_r10new`
  `83055ddaa8f6`, fast flavor present), peer fingerprint `9dcc4029` identical
  before and after, restore build reproduced the NEW md5.
  - **SBND: 201 SAME / 0 DIFF / 0 MISSING** over the 48 + 19 event `d99fix`
    manifest (Bee member hashes, `calib-pr` JSON minus the `dual_chain` timer,
    `nusel-evt` TSV).
  - **uBooNE: ZIPS 35/35 content-identical; TAGGER 34/35 identical.** The one
    diff is idx 22 ev 6805 — the documented bistable event 5384-136-6805 on
    its two documented `kine_pio_angle` states (109.512 ↔ 14.806, the arms
    swapped relative to round 9). It differed in gate rounds 1, 2, 3, 4, 6 and
    9 and was identical in 5, 7 and 8; round 9's control ran the OLD library
    twice and reproduced itself 35/35. Not attributable to this round.
- **PDVD**, both A/B pairs above: `mabc-pr.zip` member hashes byte-identical
  across every arm, including lazy-vs-eager on the 509-component cluster.

## 5. The flip is the owner's call

`relaxed_strict_img_2d_rescue_long_wtrack_fast` is **not selected by any job**.
Two reasons to leave it that way for now:

1. Kruskal and the legacy per-component Prim can differ on exact distance ties,
   and **this port's lazy path has been exercised on exactly one cluster.**
   Be precise about the two separate bodies of evidence:
   - Doc 78's 186 byte-identical SBND gate events exercised
     **`connect_graph_relaxed`**'s lazy Kruskal — the same *scheme*, a
     different function.
   - Gate round 10 below never ran the strict lazy path at all: SBND selects
     `..._long_wtrack`, not `..._long_wtrack_fast`. That is exactly why its
     201/201 is a meaningful identity gate for items 1, 3 and 4 — and exactly
     why it is **zero** witnesses for item 2's tie behaviour.
   - The only witness for **this** code is `r10fast` on 039252/8: one
     `lazy=true` line (cluster 84, ncomp=509), byte-identical `mabc-pr.zip`,
     all 25 split lines identical. Every other cluster in that event came in
     at ncomp 3–10 and took the legacy path.

   So: n = 1 cluster, not 187. Do not read the union.

   **What would validate it** (the owner's next question): run any event with
   the `_fast` flavor and census the `ncomp=` / `lazy=` debug line to find
   which PDVD `_keep` events hold a cluster above the 200 threshold at all,
   then run lazy-vs-eager and `abtest/hash_archive.py` on that subset only.
   Note 039252/16 may contribute zero witnesses — its `b12` arms ran the eager
   flavor, so there is no `ncomp` reading for it.
2. On PDVD the scope knob already dominates it. With
   `protect_stm_only_bundles=true` (doc 25 §13.11, the working mode) the stage
   on 039252/8 is **146 ms**; the 2.15× here applies to the 798 s / 1711 s
   regime you only see with that knob off.

To turn it on for PDVD:
`PDVD_PR_TLA="-S protect_graph_name='relaxed_strict_img_2d_rescue_long_wtrack_fast'"`,
or change the `protect_graph_name` default in
`wcp-porting-img/pdvd/wct-pr-perevt.jsonnet`. The value is for anyone running
with the scope knob off, and for a future detector whose busy clusters land in
the same regime.

## 6. What is left

- **The directional term is the remaining half.** The lazy gate covers the
  closest-pair walk; `dir1`/`dir2` still walk eagerly, which is why 44.3 % of
  pairs are walked. Gating them needs a separate argument, because MST #2
  consumes their verdicts directly rather than through a spanning forest.
- **Doc 78 round 1's `forced_kill` early exit** was not ported. It is an exact
  early exit (the monotone counters can satisfy every kill branch before the
  walk ends) and applies to the strict file's walk too.
- **Widen the port.** Only `..._long_wtrack` got a `_fast` sibling; the other
  six strict flavors would each need one, or the family needs a different
  selection mechanism than one name per combination.
- **Nothing here touches the root cause of the big clusters.** 039252/8's 18876
  blobs and 509 components arrive fixed in the input pctree; see doc 25 §13.11
  for why the 500 e `steiner_terminal_charge` floor is *not* that root cause.
