# Latent-pattern audit: the doc pr/11 defect classes, hunted across the toolkit

**Status:** fixes shipped in two commits (bucket 1 gated byte-identical;
bucket 2 flagged **NOT bit-identical — needs revalidation** on the pr/11
F/G precedent, though no output change was observed on the 78 gated events —
§4).  Gate labels and binary lineage in §4.  Follow-up commit: the §5.10
run-to-run bimodal `nue_score` was root-caused (pointer-ordered
`boost::out_edges` in `broken_muon_id`) and fixed — evidence, gates, and the
residual open layer in §5.10 and §4b.

**What this is.** The 1071-event SBND PR-chain census
(`sbnd_xin/docs/pr/11_pr-chain-population-scores-and-perf.md`, toolkit commit
`289d78e4`) fixed seven defect classes, each traced to a specific line:

| class | shape |
|---|---|
| A | unbounded `flag_continue` loop whose body can silently no-op |
| B | consumer ignores a producer's "may legally be absent/empty" contract |
| C | out-of-volume point → `contained_by()` returns apa/face −1 → keyed `.at(-1)` abort |
| D | self-append of a container to itself → heap blow-up (224 GB) |
| E | fatal throw on a legitimately-degenerate state |
| F | shared_ptr stored where a copy was intended (aliasing) |
| G | the same merge performed twice |

Those seven fixes cured the *sites the data happened to hit*.  This audit asks
the owner's follow-up question: do the same patterns exist elsewhere in the
toolkit, latent only because the 1071 events never reached them?  Answer: yes
— **24 fixed sites** across 14 files, plus two defects **observed live** on an
event the census marked healthy, plus a set of documented-only contract notes.

## Repro

```
# The searches (each class's generalized pattern, from toolkit root):
grep -n "while *(flag_continue" clus/src/*.cxx                 # class A
grep -n "\[0\]\.first\|\.front()\.first" clus/src/*.cxx        # class B/E2
grep -rn "contained_by(" clus/src | grep -v "apa() *== *-1"    # class C
grep -n "dpcloud([a-z_]*name[a-z_]*, " clus/src/*.cxx          # class F
grep -n "add_points(\|append(" clus/src/DynamicPointCloud.cxx  # class D

# The instrumented confirmation (§2, TEMP-AUDITDIAG probes since removed):
cd sbnd_xin && PR_JOBS=6 ./run_pr_chain_batch.sh work-nuecc48-nuf \
    work-auditdiag-nuecc48 data           # probe binary 74e7422d, 48 evts

# The gates (§4):
cd sbnd_xin && PR_JOBS=6 ./run_pr_chain_batch.sh work-mcp1kall-d59k \
    work-mcp1kall-ab30audit1 data <30 evts of the pr/11 ab30 manifest>
# ... same for work-mcp1kall-ab30audit2 with the bucket-2 binary, then
# hash_archive.py member-hash comparison vs work-mcp1kall-pr11 (M2).

# The §5.10 determinism repro (vertex pinned to geometric so the SCN layer
# is out of the picture; inject dl_weights= into a COPY of the driver):
cd sbnd_xin && sed 's|--tla-str  "reality=\$REALITY" \\|&\n            --tla-str  "dl_weights=" \\|' \
    run_pr_chain_batch.sh > tmp_run_pr_chain_geovtx.sh && chmod +x tmp_run_pr_chain_geovtx.sh
for i in 1 2 3 4 5 6; do ./tmp_run_pr_chain_geovtx.sh work-nuecc48-nuf \
    work-detfixgeo$i-469665 sim 469665; done   # pre-fix: 5x brm=4 / 1x brm=5
```

## 1. Method

Three parallel sweeps of `clus/` (loops & throws; apa/face-keyed lookups;
aliasing & duplicate merges), each seeded with the pr/11 root causes as search
patterns.  Every candidate was then re-verified against source before being
classified:

- **FIX bucket 1** — the new branch is entered *only* where the old code hung,
  threw, or was undefined behavior.  Byte-identical on every event that
  completed before (the pr/11 §6.6 argument), gated as such.
- **FIX bucket 2** — changes output on events that pass today (the F/G class).
  Fixed default-ON with no knob per the owner's decision (bug, not behavior),
  gated by attribution.
- **DOCUMENT-ONLY** — real contract holes with no data-reachable trigger
  found, dead-code interactions, or divergences that need an owner decision
  (§5).

A useful meta-observation: most bucket-1 findings are **fork-lag** — a guarded
and an unguarded copy of the *same code* coexist, because a file was forked
(M10 duplication) before the guard landed in the original.  The fix is almost
always "port the guard the sibling already has", cited per site below.

## 2. The two defects observed live (bucket 2)

Temporary WARN probes (`TEMP-AUDITDIAG`, removed after the run) were placed at
every shower-cloud seed/merge site, and the full 13-stage PR chain was run on
51 events (48 nueCC + 172230/166870/350121).  On **nueCC evt 163543** — an
event the pr/11 census completed with rc=0 — both duplicate-merge mechanisms
fired:

```
TEMP-AUDITDIAG set_start_segment RE-MERGE shower=2 seg=-1 cloud=fit             seg_npts=9   shower_npts=147
TEMP-AUDITDIAG set_start_segment RE-MERGE shower=2 seg=-1 cloud=associate_points seg_npts=85  shower_npts=1456
TEMP-AUDITDIAG add_segment       RE-MERGE shower=2 seg=-1 cloud=fit             seg_npts=11  shower_npts=156
TEMP-AUDITDIAG add_segment       RE-MERGE shower=2 seg=-1 cloud=associate_points seg_npts=140 shower_npts=1541
```

**2a. `examine_showers` merges the same segment twice** (defect-G twin).
`NeutrinoShowerClustering.cxx:2342-2344` runs `shower->add_segment(sg)` and
then `shower->set_start_segment(sg)`; both merge `sg`'s clouds into the
shower, so `sg`'s points landed twice — the first two probe lines (9 "fit" /
85 "associate_points" duplicate points on this event).

**2b. Re-added segments re-merge their clouds** (the general mechanism behind
2a).  `Shower::add_segment` / `set_start_segment` merged the segment's DPC
*unconditionally*, while graph membership was already idempotent
(`std::set::insert`).  Any caller that re-adds a member segment duplicates its
points; the audit found the fresh `used_segments` set handed to
`complete_structure_with_start_segment` at `NeutrinoShowerClustering.cxx:2346`
does exactly that for every pre-existing member the worklist reaches — the
last two probe lines (11 / 140 duplicates).

*Fix (one mechanism for both):* a **membership gate** in
`Shower::add_segment` and `Shower::set_start_segment`
(`clus/src/PRShower.cxx`) — capture `has_edge(descriptor)` before the graph
add and skip the cloud merge for an already-member segment.  Seeding an
as-yet-cloudless shower stays allowed (nothing to duplicate).  This makes the
cloud merge idempotent, matching the graph semantics, and subsumes deleting
the redundant call in 2a.

Duplicate points bias every k-d proximity/density query against the shower
cloud and inflate `shower_track` Bee dumps; by the pr/11 §6.7 measurement the
charge sums are distance-based and *mostly* insensitive, but defect F moved a
`nue_score` — hence bucket 2, needs revalidation.

**2c. `add_shower` still seeds by aliasing** (defect-F twin, structural).
`PRShower.cxx:260/:274` seeded the target shower's cloud with the source
segment's own shared_ptr — byte-for-byte the pr/11 defect F, cured at the
file's four *other* seeding branches (`clone_dpc`) but missed here because
`add_shower` batches manually instead of calling `add_segment`.  Once aliased,
a later pass of the same absorb loop appends the accumulated cloud to itself
(the 224 GB defect-D shape), and every consumer reading that *segment's* cloud
answers for the whole shower.  The probes did **not** observe the seeding
branch firing in 51 events (`associate_points`-cloudless showers are needed —
structural reachability only), but the code path is unconditional when
reached.  *Fix:* `clone_dpc` on both seeding branches + `!=` self-batch test +
the same membership skip, mirroring `add_segment`.

## 3. Bucket-1 fixes (crash/hang/UB guards, byte-identical on passing events)

### 3.1 Class A — loops that can spin forever

| site | defect | fix |
|---|---|---|
| `PRShower.cxx` `get_last_segment_vertex_long_muon` (~:545) | `s_vtx == m_start_vertex` branch sets `flag_continue=true` with **no state change**; sole escape is `find_vertices()` moving `s_vtx`, which silently no-ops for an invalid segment descriptor (`PRGraph.cxx` returns `{nullptr,nullptr}`) or a self-loop — the exact class-A mechanism | end-of-pass progress guard: WARN + break when neither `s_seg` nor `s_vtx` moved |
| `NeutrinoPatternBase.cxx` `proto_extend_point` (~:522) | 1-cm stepping walk with **no iteration cap**; its relative `Facade_Cluster.cxx:343` caps at 400 (the authors knew the walk needs a bound); `walk_history` grows unboundedly too | `counter < 4000` cycle backstop, WARN when hit — **not** 400, see the cap-sizing lesson in §4 |
| `TaggerCheckSTM.cxx` crawl loop (~:1096) | third copy of the same walk, also uncapped (forked before the cap) | same 4000 cap + WARN |
| `NeutrinoPatternBase.cxx` `break_segments` kink loop (~:910) | visited-set progress guard is **bypassed when `break_idx == INVALID_STEINER_INDEX`** (empty steiner kNN); the pass then re-runs from the same `test_start_p` — a deterministic period-1 fixed point | exact stationarity detector (two consecutive passes with identical input *and* outcome ⇒ provably periodic ⇒ WARN + break) plus a 1000-pass backstop cap |

### 3.2 Class B/E2 — `[0]`/`front()` on legally-empty kNN results

The producer contract (pr/11 class B): a steiner cloud may be absent *or
present with zero points* — `Dataset::size()` counts **arrays**, not points
(`Clustering_Util.cxx:81`), so `size_major()` is the only correct emptiness
predicate.  `kd_steiner_knn` on a zero-point cloud returns an empty vector and
`[0]` reads past the end.

| site | note | fix |
|---|---|---|
| `NeutrinoPatternBase.cxx` `do_rough_path` (:83) | `has_pc` was checked only *after* the kNN | `has_pc`+`size_major` pre-check + empty-kNN check → return empty path (already a legal outcome) |
| `NeutrinoPatternBase.cxx` `do_rough_path_reg_pc` (:135) | same shape on `kd_knn` | empty-kNN check |
| `TaggerCheckSTM.cxx` `do_rough_path` (~:870) | verbatim fork of the NPB function, same hole | same guard |
| `TaggerCheckSTM.cxx` crawl block (~:1054-1156) | five bare `[0].first` (`:1073, :1122, :1146, :1152, :1156`); the NPB original guards every one (`:509-512, :534, :592`) — fork-lag | one `size_major() > 0` term on the block's gate covers all five (kNN(1) on a non-empty cloud always returns a result) |
| `TaggerCheckSTM.cxx` `adjust_rough_path` (~:911) | **new interaction with the pr/11 class-E fix**: `do_single_tracking` may now legally clear its fit, and `.at(0)` on the empty path throws the very `out_of_range` class E removed | empty-path early return (caller already falls back to the un-adjusted path) |
| `MyFCN.cxx:346` | guard exists but uses the **superseded predicate** (`.size()`), so a zero-point cloud sails through to `[0]` at :358/:407; the corrected predicate already landed twice elsewhere (`Clustering_Util.cxx:87`, `NPB:212`) | add the `size_major()` term |
| `TrackFitting.cxx` steiner block (~:2254) | gate tested *presence* (`has_pc`) but `.front()` at :2268-2270 needs *non-emptiness* | `size_major() > 0` term on the gate |
| `PRGraph.cxx` `find_vertices` (:127) | dereferenced possibly-null graph vertices for the ordering distance; every caller defends against nulls in the *returned pair*, so null is an anticipated state | null check → return the unordered pair |

### 3.3 Class C — apa/face −1 into keyed lookups

Crash primitive: both `IPCTransform` impls (`PCTransforms.cxx` `T0Correction`,
`SCECorrection.h`) do `.at(apa).at(face)` in `forward`/`backward` — an
uncontained point aborts with the bare `std::out_of_range: map::at` of the
class-C crash.  `PR::Fit::paf` defaults `{-1,-1}` (`PRCommon.h`) and two
producers legally leave it there for out-of-volume fit points
(`TrackFitting.cxx multi_trajectory_fit` ~:3990, `PRSegment.cxx clear_fit`
~:101); the guards belong at the consumers because the sentinel default is
load-bearing elsewhere (`MultiAlgBlobClustering.cxx:915` tests `>= 0`).

| site | defect | fix |
|---|---|---|
| `TrackFitting.cxx` `trajectory_fit` init loop (~:4186) | unguarded `backward()` on each path point; **cannot skip** (Eigen `pos_3D` rows are index-aligned) | identity fallback: use the corrected point as its own raw seed when uncontained (such points have no 2D data; the solve regularizes around the seed) |
| `TrackFitting.cxx` `trajectory_fit` path loop (~:4485) | unguarded `forward()` | `continue`, mirroring the `do_single_tracking` 2nd-pass guard (:8533); the three output vectors append together so the skip stays consistent |
| `TrackFitting.cxx` `trajectory_fit` projection loop (~:4618) | `wpid_offsets.find()` **dereferenced without an end() check** — UB, while every sibling (:6081, :6437, :7184) checks; sat behind the :4485 throw | WARN + `continue`; a resulting size mismatch lands in the now-non-fatal class-E consistency check, which drops the fit cleanly |
| `TrackFitting.cxx` `fit_point` (~:3436) | safe today only via an *unstated* caller contract (both callers pre-check containment) | defense-in-depth early `return init_p`, mirroring the function's own no-2D-associations return |
| `NeutrinoStructureExaminer.cxx` (~:1265) | the same function as the fixed class-C crash, 40 lines above the pr/11 guard (:1310), same `sg->fits()` provenance; unhit because this loop breaks on the first non-dead point | `continue`, mirroring :1310 verbatim |
| `TaggerCheckSTM.cxx` kink sweeps (~:1394, ~:1476) | consume `fit.paf` blindly: `backward(paf.at(i)...)` plus `inside_dead_region` (itself an unguarded `convert_3Dpoint_time_ch` wrapper) | `paf.at(i) < 0` → `continue` (an out-of-volume point cannot be a kink candidate) |
| `FiducialUtils.cxx` `inside_dead_region` (:35) | three unguarded converts; the crash primitive one call away from any −1 | `apa/face < 0` → `return false` (no volume ⇒ no dead region) |
| `NeutrinoTaggerSinglePhoton.cxx` (~:2203) | overwrites its safe `(0,0)` init with a possibly −1 wpid; **both sibling taggers guard this and say why** (`TaggerCheckNeutrino.cxx:683`, `NeutrinoTaggerNuE.cxx:2762`); currently dormant (no `ctx.apa` consumer yet) | restore the sibling guard |

### 3.4 Class E — fatal throws on legitimately-degenerate states

| site | defect | fix |
|---|---|---|
| `Facade_Cluster.cxx` `get_two_extreme_points` (~:2069) / `get_main_axis_points` (~:2031) | throws when **every** point is excluded — a legal state: `connect_graph.cxx:266` marks `should_exclude` for any point ≥ 1 cm from its reference cluster, with no non-empty guarantee.  The closest surviving twin of the fixed class-E crash | all-excluded ⇒ WARN + fall back to the *unfiltered* extremes (real points of this cluster); the throw remains for a truly point-less cluster |
| `clustering_neutrino.cxx` (~:419/:428 and ~:597/:605) | `Cluster* largest_cluster = 0` assigned only when a separated cluster has `npoints() > 0` ⇒ **null deref** (segfault, not even a throw) on an empty separation | null guard: skip the temp refinement (`flag_enable_temp` stays false = "no separation result") but still merge back |
| `clustering_connect.cxx` (~:274) | the guard `assert(npoints() > 0)` names the right hazard but tests it **only under !NDEBUG** | real `continue` beside the assert |

### 3.5 Class D — structural self-append guard

`DynamicPointCloud::add_points` (all four overloads, `DynamicPointCloud.cxx`)
now refuses `&other == this` / `&points == &m_pts` with a **loud WARN**
naming the aliasing-defect suspicion, then no-ops (the points are by
definition already present).  Self-append is UB, not a bounded 2× —
`DPCBatch::append` inserts from iterators into the growing destination, which
reallocation invalidates mid-copy.  The pr/11 round deliberately left
`add_points` unguarded so call-site fixes would not *mask* aliasing; the WARN
honors that intent — a future aliasing defect is reported, not hidden.

## 4. Verification

Binary lineage (`md5sum local/lib/libWireCellClus.so`):

| md5 | contents |
|---|---|
| `e154d50b` | pr/11 baseline (`289d78e4` + `04545d80`, no audit change).  Rebuilt from clean HEAD mid-audit and reproduced **exactly** — waf's content signatures make the baseline identity airtight |
| `74e7422d` | baseline + TEMP-AUDITDIAG probes (log-only; since removed) — the output baseline for the nueCC48 arms |
| `a71d0dfe` | + all bucket-1 fixes (everything except `PRShower.cxx`) |
| `566ad005` | + bucket-2 (`PRShower.cxx`), walk caps at 400 — superseded, see the cap lesson below |
| `8f66490a` | same with walk caps at 4000 — **the shipped binary** |
| `582704e0` | shipped + TEMP-DETDIAG probes in `broken_muon_id` (log-only; since removed) — the §5.10 root-cause evidence |
| `07a78447` | + the deterministic out-edge order fix in `broken_muon_id` (§5.10, follow-up commit) |

`libWireCellRoot.so` `1e0908be` throughout (no `root/` file touched).
`wcdoctest-clus`: 49 cases / 565 assertions pass on every build.  Freshness
proof (M1) done before each arm.  Comparisons: `abtest/hash_archive.py`
member content (M2) for `mabc-pr.zip` + `pctree-pr-evt*.tar.gz` (with a
negative control: cross-event comparison correctly reports DIFF), plus
`pr_scores_table.py` rows with timing/RSS columns excluded.

**Baselines are fresh same-code arms, not the census tag.**  The pr/11 census
tag `work-mcp1kall-pr11` holds *first-pass* outputs for events that passed —
i.e. pre-defect-F binaries — so its evt 166870 still carries the pre-F
`nue_score` 2.6070938 (the documented pr/11 §6.7 1/30 change).  Gate arms are
therefore compared against `work-mcp1kall-ab30base` (30 MCP2025C events,
binary `e154d50b`) and `work-auditdiag-nuecc48` (48 nueCC events, binary
`74e7422d`).

| arm | binary | events | archives | score rows |
|---|---|---|---|---|
| `work-mcp1kall-ab30audit1` (bucket 1 only) | `a71d0dfe` | 30 | **60/60 identical** | **30/30 identical** |
| `work-mcp1kall-ab30audit3` (shipped) | `8f66490a` | 30 | **60/60 identical** | **30/30 identical** |
| `work-auditfinal-nuecc48` (shipped) | `8f66490a` | 48 | **96/96 identical** | 46/48 identical; both diffs proven **pre-existing run-to-run noise**, below |
| `work-auditfix-163543` (smoke, the probe-confirmed event) | `566ad005` | 1 | identical | identical |

All 78 + 4 events rc=0; **none of the new guard WARNs fired on any passing
event** (the only WARN in any arm is the pr/11 class-A guard doing its job on
evt 352365).

**The cap-sizing lesson (why 4000, not 400).**  The first bucket-2 build
capped the two ported stepping walks at 400 like `Facade_Cluster.cxx:343`.
The nueCC48 arm then moved 3/48 score rows; none of the capped sites WARNed
(the STM copy of the cap had no WARN — the audit's own mistake), and
bisection attributed two of the three to the cap: these crawls legitimately
walk the full track length in ~1 cm steps, so 400 passes truncates real SBND
tracks (~6.4 m diagonal = ~640 steps).  At 4000 — beyond any legitimate walk,
still milliseconds for a true cycle — both events returned to
baseline-identical.  Rule extracted: **a loop bound ported from a sibling
must be re-derived from this walk's physical scale, and every cap needs a
WARN so a binding cap is observable.**

**The two residual nueCC48 row diffs are pre-existing nondeterminism, proven
on a single binary.**  (1) evt 469665: five same-driver runs across three
binaries produce `nue_score` −0.443180 (`brm_n_mu_segs`=4) *and* +1.217744
(=3), **including both values from two runs of the shipped binary** — the
§5.10 finding, since root-caused and fixed (see §5.10 and §4b).  (2) evt
422851: `kine_reco_Enu` 1455.9768 vs 1455.9769 across two runs of the
shipped binary (evt 138009 flickers identically at the same magnitude).
Neither is attributable to any audit fix.

**Regression check** — the four named pr/11 crash events (352365 class A,
49951 C, 399118 D, 386354 E) complete rc=0 on the shipped binary
(`work-auditreg2-crash4`; also on `566ad005` in `work-auditreg-crash4`),
evt 352365 in ~6 s wall.

**Bucket-2 measured effect.**  The membership gate provably removes in-memory
duplicate points (probe evidence, §2) but produced **no output change on any
of the 78 gated events** — on evt 163543 the duplicates never reached the
dumps or scores.  The NOT-bit-identical flag is therefore structural rather
than observed: on events where an `add_shower` alias or a re-merge feeds a
shower-cloud consumer (the pr/11 defect-F mechanism that moved evt 166870's
`nue_score`), outputs can move, and revalidation at the next population run
should expect that.

## 4b. Follow-up round: the §5.10 bimodality — root cause, fix, gates

All on the fix binary `07a78447` (= `8f66490a` + the `broken_muon_id`
out-edge-order fix, nothing else; committed together with this doc update)
against baseline `8f66490a` rebuilt from clean committed HEAD.

**Repeat-identity gate (the determinism claim).**  Six pinned-vertex runs
(`dl_weights=` → geometric) `work-detfixgeo{1..6}-469665` and two
production (DL-on) runs `work-detfixdl{1,2}-469665`: **every scalar
`T_tagger` branch equal and `hash_archive.py`-member-identical
`pctree`/`mabc` in all pairs**.  Pre-fix, the same six-run pinned-vertex
protocol on probe binary `582704e0` split 5-vs-1 between the two walk modes
(`work-detgeo{1..6}-469665`) with the §5.10 AMBIGUOUS probe line firing on
the flip vertex.

**A/B arms (the scope-of-change claim), current config both sides:**

| comparison | archives | `T_tagger` scalar rows |
|---|---|---|
| `work-detfix-ab30` vs `work-detbase-ab30` (30 MCP data evts) | **60/60 identical** | 17 identical + 13 absent-in-both, **0 diffs** |
| `work-detfix-nuecc48` vs `work-detbase-nuecc48` (48 nueCC evts) | **96/96 identical** | 45 identical + 1 absent-in-both, **2 diffs** (`brm_*` only): evts 389538, 52672 |

**The two row diffs are pre-existing coin flips being pinned, proven on the
baseline binary alone.**  Four baseline samples each
(`work-detbase-nuecc48` + `work-detbaserep{1,2,3}-<evt>`): evt 389538 is
run-to-run bimodal `nue_score` 2.935744 / `brm_acc_length` 67.79 ↔
2.470848 / 82.09 (2-vs-2); evt 52672 flips `brm_acc_length` 67.05 ↔ 64.84
(1-vs-3, `nue_score` unchanged).  The fix deterministically lands each on
one mode (82.09; 64.84).  So across the 78-event population the fix's only
output effect is pinning two latent coin-flip events — everything else is
row- and archive-identical.  (Evt 469665 itself is unchanged in this arm
pair: under the committed config its walk has no ambiguous vertex; its flip
reproduces under the audit-era config and under the pinned geometric
vertex, where the fix pins it.)

**The cross-era config delta (why the audit-era tags are stale as
baselines).**  `work-detbase-nuecc48` vs `work-auditfinal-nuecc48` — the
*same* binary md5 `8f66490a` on both sides — differs on **47/48** events'
`T_tagger` rows including `nu_x/y/z`, because the audit-era arms ran under
then-present uncommitted WIP in the two SBND cfg jsonnets (since reverted).
Two lessons: (1) clean rebuild of committed HEAD reproduces `8f66490a`
byte-for-byte, so the audit *binaries* carried no uncommitted code — the
era delta is cfg-only; (2) the delta is **invisible in the archives**
(`pctree`/`mabc` were 96/96 identical even cross-era) — archive identity
does NOT cover `T_tagger`; row comparison is mandatory in any PR-chain
gate.  Future gates must re-run a fresh baseline arm rather than compare
against `work-auditdiag-*`/`work-auditfinal-*`/`work-mcp1kall-ab30audit*`.

## 5. Documented-only: surfaced, not fixed

Contract holes and divergences that either lack a data-reachable trigger, sit
in dead code, or need an owner decision.  None of these were changed.

1. **`.at()`-as-assert in `TrackFitting.cxx`** (:6023, :6646, :6845 on
   `segment_point_index_map`; :4158 on `m_3d_to_2d`).  Deliberate
   (`// .at() catches missing-key bugs`), and fill/use guards are consistent
   today.  Data-reachable only through the documented setS descriptor
   aliasing hazard (`NeutrinoStructureExaminer.cxx:279`) — two aliased
   segments with different fit counts would collapse to one map entry and the
   longer one's `.at()` kills the job with an unlocated `out_of_range`.
2. **`Facade_Cluster.cxx` ~:2525**: when `contained_by(p_test)` misses, the
   code re-queries with the nearest cluster point but does **not re-check**
   the result; a second miss silently yields `wire_direction == (0,0,0)` and
   all three plane angles collapse to 0 (degenerate 2D projections).  Fixing
   this changes a wrong-but-finite result on unknown events — owner call.
3. **`TaggerCheckTGM.cxx:706-712`**: the `pick` lambda's fallback
   `m.begin()->second` is guarded only by the `wpid_U_dir.empty()` early
   return at :703 plus the (currently true) symmetry that U/V/W maps fill
   together.  Fragile, not broken.
4. **`Shower::add_segment` returns void** and discards
   `TrajectoryView::add_segment`'s `false` (invalid descriptor) at 8
   single-pass call sites — a silently *lost* segment, not a hang (the two
   loop sites got progress guards in pr/11; the long-muon walk got one here).
   Likewise `merge_vertex_into_another` (`NeutrinoPatternBase.cxx` ~:1418)
   returns `true` even when its internal `remove_vertex` failed.
5. **`dpcloud(name, ptr)` aliases by design** (`PRCommon.h:78`) with nothing
   in the name saying so — the root cause of defects D/F and §2c.  A rename
   (`set_dpcloud_shared`?) or a cloning overload would prevent the next
   accident.  Same trap: `Segment::particle_info(shared_ptr)` — all ~30
   current callers pass fresh objects, but the stored object is mutated in
   place (`set_pdg`), so the first `seg2->particle_info(seg1->particle_info())`
   will silently couple two segments.
6. **`PRSegmentFunctions.cxx:55-65`**: replacing a segment's cloud leaves the
   previous `global_indices` in place when the new call passes none — stale
   indices describing a different point set.  Dormant today only because the
   sole consumer (`transfer_info_from_segment_to_cluster`) has zero call
   sites *and* defaults to the misspelled cloud name `"associated_points"`
   (everything else uses `"associate_points"`) — two latent bugs cancelling.
   Worth deciding before anyone wires it up.
7. **`TaggerCheckSTM.cxx` ~:1115 momentum trick** retains
   `dir1 + dir * 5.0 * units::cm` where the NPB original was deliberately
   corrected to the dimensionless `dir1.norm() + dir * 5`
   (`// both terms are now dimensionless`).  A physics-affecting fork
   divergence (M15): surfaced for the owner, not silently "fixed".
8. **Inconsistent out-of-volume semantics in guarded loops**: several
   `is_good_point` sample loops guard apa==−1 with no `else` (point silently
   *not* counted bad: `NeutrinoStructureExaminer.cxx` ~:803, ~:1352, …;
   `NeutrinoOtherSegments.cxx` ~:1104, ~:1251) while three sibling loops in
   the same file count it bad with the comment "outside all TPCs → bad
   sample" (:92, :246, :604).  A physics question, not a stability one.
9. **Producers of the `paf {-1,-1}` sentinel** (`multi_trajectory_fit`,
   `Segment::clear_fit`) are left as-is on purpose: consumers now guard, and
   the sentinel is the documented "no volume" encoding.
10. **The bimodal `nue_score` (evt 469665) — root-caused and FIXED in the
    follow-up commit.**  Found while attributing the gates: nueCC evt 469665
    flips between `nue_score = −0.443180` (`brm_n_mu_segs = 4`,
    `brm_acc_length = 143.3`) and `+1.217744` (3, 55.4) across repeated runs
    of the *same* binary under the same `setarch -R` driver.  Exactly 6
    scalar `T_tagger` branches move, all `brm_*` — the
    `broken_muon_id` long-muon chain (`NeutrinoTaggerNuE.cxx:1370`,
    prototype `NeutrinoID_nue_tagger.h:1010`).

    **Root cause (probe-confirmed).**  The PR graph stores vertex and edge
    descriptors with `boost::setS` (`PRGraphType.h:91-93`), so
    `boost::out_edges` iterates a vertex's edges in *descriptor pointer
    order*, which varies run to run (tcmalloc heap layout, even under
    `setarch -R`).  `broken_muon_id`'s step A walks `out_edges` and takes the
    **first** segment passing the back-to-back collinearity cut
    (180°−angle < 15°, length > 6 cm) with a `break`.  On evt 469665 (vertex
    pinned to geometric via `dl_weights=`), a TEMP-DETDIAG probe caught two
    segments passing the cut at one vertex, and the out-edge order literally
    flipping between runs:

    ```
    geo1: stepA AMBIGUOUS at vtx 16: 2 candidates [29 27], chosen 29 (12.39 cm)
    geo4: stepA AMBIGUOUS at vtx 16: 2 candidates [27 29], chosen 27 (87.75 cm)
    ```

    Five of six probe runs walked seg 29 (`brm` = 4 / 55.9 cm), one walked
    seg 27 (5 / 137.7 cm) — a physics-visible score resolved by heap
    addresses.  The prototype iterates a pointer-keyed
    `std::set<ProtoSegment*>` at the same spot, i.e. its order is equally
    arbitrary — so no parity is broken by choosing a stable one.

    **Fix (no knob — determinism bug):** collect the vertex's out-edge
    segments into an `IndexedSegmentSet` (stable graph-index = insertion
    order) and run the first-match scan on that.  Verification, all on fix
    binary `07a78447`: 6/6 pinned-vertex repeats and 2/2 production (DL-on)
    repeats **fully identical** — every scalar `T_tagger` branch equal and
    `hash_archive.py` member-identical `pctree`/`mabc` archives
    (`work-detfixgeo{1..6}-469665`, `work-detfixdl{1,2}-469665`); plus the
    §4b arms.

    **The apparent "second layer" resolved — a config delta, not SCN
    instability.**  Mid-investigation the DL vertex seemed to sit in one of
    two positions ~45 cm apart depending on the binary.  §4b proves the
    real cause: the audit-era arms ran while the tree carried
    *uncommitted* WIP in `cfg/pgrapher/common/clus.jsonnet` +
    `cfg/pgrapher/experiment/sbnd/clus.jsonnet` (runtime inputs via
    `WIRECELL_PATH`, not baked into the binary), since reverted by the
    concurrent session.  Re-running the byte-identical binary `8f66490a`
    under the committed config moves 47/48 nueCC rows (including
    `nu_x/y/z`), while within any one config every repeat is fully
    reproducible.  No SCN run-to-run instability was observed anywhere in
    this round; the M4 caution stands on its own history, not on this
    evidence.  The one residual: `kine_reco_Enu` flickers ±0.0001 MeV on
    evts 422851/138009 within a fixed binary+config (FP-level,
    unattributed).

    **The residual class.**  `boost::out_edges` is iterated at 132 sites in
    13 `clus/` files; most are order-insensitive (counting, any-match,
    inserting into Indexed*Sets), but any *first-match-break*, tie-broken
    min/max, or float-accumulation loop over raw `out_edges` shares this
    hazard.  A dedicated sweep is a follow-up candidate; this commit fixes
    the one site with observed physics impact.

## 6. Relation to the pr/11 fixes

Nothing in this round touches the pr/11 fixes themselves; the previously
failing events stay fixed (§4 regression check).  The membership gate (§2b)
generalizes the pr/11 defect-G repair from one call site to the API, and the
`add_shower` clone (§2c) completes the defect-F repair at the one seeding
branch pr/11 missed.  The `DynamicPointCloud` self-append guard (§3.5) closes
the defect-D *class* structurally.
