# APA1 has no 2D charge in the PR tail: every TPC1-only cluster loses its steiner graph

**Status: FINDING ONLY — no code changed.** This needs an owner decision (see
§6): the natural fix touches `PointTreeMerging`, which runs inside the
byte-identical-gated all-APA clustering stage.

**Update 2026-07-23:** §7 adds a second, independent failure mode of the same
gap (the boundary wcps collapse to the ORIGIN, so even cathode-straddling
clusters that *do* have a steiner graph get a fabricated exit point), and §8
proves by direct probe that APA1's `ctpc_a1f0p*` reaches the merge intact and
is dropped there — nothing upstream is broken. Owner ruling recorded in §7:
FC is a fiducial-containment verdict, independent of readout-window
truncation, so evt285999 main 18 *should* be FC = true.

Investigated because SBND MCP2025C evt 285185 bundle grp 3 (main cluster 12,
flash t0 = −1171.4 µs) is not labeled FC even though it displays as a compact
65 cm track fully inside the fiducial volume.

## Repro

```bash
cd /nfs/data/1/xqian/toolkit-dev/toolkit && wcbuild
SB=/nfs/data/1/xqian/toolkit-dev/wcp-porting-img/sbnd/sbnd_xin
cd $SB && SBND_INPUT_DIR=$PWD/input_files_reco1/extracted-mcp2025c-10evt \
  SBND_WORK_ROOT=$PWD/work-mcp10 ./run_nusel_evt.sh data all
grep "produced no steiner_graph" work-mcp10/nusel_evt*/wct_nusel_evt*.log
# ctpc datasets actually present in the persisted tree:
tar tzf work-mcp10/ql_evt285185/pctree-evt285185.tar.gz >/dev/null   # see §3
```

Two temporary env-gated probes (`WCT_STEINER_DEBUG` in
`SteinerGrapher::find_steiner_terminals`, `WCT_ACT_DEBUG` in
`ImproveCluster_1::get_activity_improved`) produced the per-plane numbers
quoted below. Both have been reverted; the tree is clean.

## 1. Why cluster 12 is not FC — the proximate chain

It is **not** a containment mis-read (unlike the two bugs fixed in `d48dc647`).
The containment check never runs:

```
create_steiner_tree: only 0 steiner terminal(s) found (need >=2), returning empty graph
CreateSteinerGraph: create_steiner_tree produced no steiner_graph for main 12, skipping transfer
```

`Facade::cluster_fc_check` (`clus/src/Clustering_Util.cxx`) opens with a
conservative guard: no `steiner_pc` ⇒ return `is_fc = false`. So cluster 12 is
reported not-FC by default, without any boundary point ever being tested.

## 2. Why there are no terminals — no charge in the retiled cloud

Terminals come from `find_peak_point_indices`, which requires
`calc_charge_wcp(...)` to return **both** `charge > 4000` **and** a quality flag
that demands all three planes exceed the cut (or be dead).

`WCT_STEINER_DEBUG` on evt 285185, per cluster, for the two graphs the steiner
stage builds:

| cluster | graph | npts | u>cut | v>cut | w>cut | dead(u,v,w) | terminals |
|---|---|---|---|---|---|---|---|
| 12 | `basic_pid` (original) | 107 | 31 | 50 | 28 | 3,0,0 | **8** |
| 12 | `ctpc_ref_pid` (retiled) | 100 | **0** | **0** | **0** | 95,97,74 | **0** |
| 18 | `basic_pid` (original) | 6515 | 4004 | 3705 | 3984 | 0,53,41 | **377** |
| 18 | `ctpc_ref_pid` (retiled) | 6421 | 3 | 4 | 4 | 6268,6140,6224 | **2** |

The *original* cluster's charge is healthy. The *retiled* cluster
(`improve_cluster_2`, which `CreateSteinerGraph` uses as its retiler) has
essentially zero charge and reads as almost entirely dead.

Cluster 18 is the discriminator that rules out "the track is just faint":
it is the **brightest** cluster in the event (Bee per-point charge median
7539, 79 % of points above 4000; cluster 12's median is 3333) and it still
collapses to 2 terminals after retiling.

## 3. Root cause — the merged tree carries ctpc for APA0 only

`WCT_ACT_DEBUG` in `get_activity_improved`, which is what fills the retiled
blobs' charge, is unambiguous — every APA1 query comes back completely empty,
including the dead-channel ranges:

```
ACT_DEBUG clus=12 apa=1 face=0 t=[240,560)  ... good(u,v,w)=0,0,0     dead(u,v,w)=0,0,0
ACT_DEBUG clus=15 apa=0 face=0 t=[52,1840)  ... good(u,v,w)=8071,7717,3078 dead=2,2,1
ACT_DEBUG clus=15 apa=1 face=0 t=[420,1844) ... good(u,v,w)=0,0,0     dead(u,v,w)=0,0,0
ACT_DEBUG clus=17 apa=0 face=0 t=[964,2304) ... good(u,v,w)=2171,1462,820  dead=0,0,0
ACT_DEBUG clus=17 apa=1 face=0 t=[0,2304)   ... good(u,v,w)=0,0,0     dead(u,v,w)=0,0,0
```

`Grouping::get_overlap_good_ch_charge` looks up a local PC named
`ctpc_a<apa>f<face>p<U|V|W>` and **silently returns an empty map when that
dataset is absent** (`Facade_Grouping.cxx:694`). The persisted tree contains:

```
pointtrees/285185/live/pointclouds/namedpcs/ctpc_a0f0pU   (+ pV, pW)
```

and **no `ctpc_a1f0p*` at all**. The time/wire ranges in the APA1 queries are
sensible, so the query arguments are fine — the data simply is not there.

The naming code is correct: `PointTreeBuilding::add_ctpc` uses
`anode_ident = m_anode->ident()` (`PointTreeBuilding.cxx:335,350`), so each
per-APA imaging→clustering job does write its own `ctpc_a{0,1}f0p*`.
(`Aux::add_ctpc` in `aux/src/SamplingHelpers.cxx:348` *does* hardcode the APA
index to 0, but its only caller is `root/src/UbooneClusterSource.cxx` —
single-APA uBooNE, harmless there. It is not this bug.)

The loss happens at the **per-APA → all-APA merge**. Two components implement
that merge with the same `opflash`-only default, so it matters which one is in
this path:

- `clus/src/PointTreeMerging.cxx` (legacy per-APA path), and
- `match/src/QLMatching.cxx:1050` (joint path, `m_root_pcs_to_merge{"opflash"}`
  at `QLMatching.h:824`).

**The active path is QLMatching.** `run_ql_evt.sh` sets `JOINT=true` by default
(line 80) and the evt 285185 log confirms `joint=true`, i.e.
`clus.per_apa (×2) → FlashTensorToOpticalPCs → QLMatching(joint) → MABC`, with
**no `PointTreeMerging` node at all**. The pctree the PR job loads is
QLMatching's output, so QLMatching is the dropper here; `PointTreeMerging`
would be the equivalent culprit under `-s/--per-apa` (`SBND_JOINT=0`).

```cpp
// QLMatching.cxx, multi-APA branch
auto* root_live = runs.front().root_live.get();          // input[0] keeps ALL its root PCs
for (std::size_t k = 1; k < runs.size(); ++k) {
    merge_pct(root_live, runs[k].root_live.get(), m_root_pcs_to_merge);  // only listed names
}
```

`merge_pct`'s own comment names the casualty outright:

```cpp
// Only names in root_pcs_to_merge are merged; everything else
// (flash/light/flashlight, per-anode ctpc_a*, ...) is dropped from the
// source roots, exactly as the standalone PointTreeMerging did.
```

and SBND configures (`cfg/pgrapher/experiment/sbnd/qlmatching.jsonnet:337`):

```jsonnet
root_pcs_to_merge: ['opflash'],
```

So the merged all-APA grouping inherits APA0's `ctpc_a0f0p*` and
`dead_winds_a0*` from input[0] and **drops APA1's entirely**. Everything
downstream — the persisted pctree and the whole PR tail — has 2D charge and
dead-channel information for APA0 only.

The same `root_pcs_to_merge: ['opflash']` appears in the PDHD
(`qlmatching.jsonnet:535`) and PDVD (`qlmatching.jsonnet:875`) configs and in
SBND's `clus.jsonnet:301`, so any multi-APA detector running a steiner/PR tail
on the merged tree is exposed the same way. Not checked here.

## 4. Scale — it is systematic, not one event

All 10 MCP2025C events, every cluster the steiner stage touched, bucketed by
how much of the cluster lies in TPC0 (x<0). NB the TPC fractions come from the
Bee dump of the *original* cluster while the failure is recorded on the
*retiled* one; the point counts differ, but the separation is total:

| bucket | clusters | no steiner graph | rate |
|---|---|---|---|
| pure TPC1 (100 % x>0) | 32 | **32** | **100 %** |
| mostly TPC1 (<10 % x<0) | 11 | 4 | 36 % |
| mixed (10–90 %) | 16 | 0 | 0 % |
| mostly/pure TPC0 (>90 %) | 58 | 6 | 10 % |

**Every cluster lying entirely in TPC1 fails, without exception.** The only six
TPC0 failures are trivially small clusters (4, 6, 9, 12, 15, 28 points) that
legitimately cannot form a steiner tree. Mixed clusters survive on their TPC0
half alone.

## 5. Impact

- **FC/STM/Neutrino verdicts.** Any main cluster wholly inside TPC1 gets the
  conservative `is_fc=false` with no check performed. evt 285185 main 12 is
  exactly this.
- **Degraded, not only absent.** Clusters that straddle the cathode do build a
  steiner graph, but from TPC0 charge only — e.g. evt 285185 cluster 17 keeps
  238 terminals instead of the 570 its original cloud supports.
- **Not confined to out-of-time cosmics.** This is TPC-dependent, not
  t0-dependent. In evt 285185 the in-beam **nu-candidate** cluster 21
  (t0 = +1.5 µs) drops from 203 terminals to 7. The PR tail is therefore
  reasoning about neutrino candidates on a degraded steiner representation
  whenever they have TPC1 content.
- `FiducialUtils`' 2D queries (`get_closest_points`, `inside_dead_region`) read
  the same per-APA ctpc/dead structures, so the SP-gap and dead-volume checks
  are equally blind in APA1. Worth a separate look.

## 6. Why this is not being fixed here — owner decision needed

The obvious fix is to merge the per-APA `ctpc_*` / `dead_winds_*` root PCs
(they are uniquely named per APA, so there is no key collision). Three reasons
it is not a repeat of the "correctness fix, no knob" pattern used for the two
bugs in `d48dc647`:

1. **It is not confined to the PR tail.** `steiner`/`improve_cluster_2` are
   PR-only (verified: they appear solely in the `cm_by_name` block, not in
   either main `cm_pipeline`). But the merge lives in **QLMatching**, which is
   squarely inside the byte-identical Q/L matching gates, and its output tree
   feeds the all-APA clustering stage too. Making APA1 ctpc appear could change
   matching and all-APA clustering output. Per escalation rule 1 that needs a
   default-OFF knob and a full gate run, not a silent change.
2. **It may be deliberate.** Both merge implementations carry explicit comments
   that per-anode root PCs are dropped on purpose ("everything else
   (flash/light/flashlight, per-anode `ctpc_a*`, ...) is dropped"). ctpc is
   large; this may be a memory decision. Per M15 the divergence should be
   surfaced, not silently reverted.
3. **Cost is unmeasured.** Merging both APAs' ctpc roughly doubles that part of
   the tree in memory and in every persisted tarball.

Options for the owner, in increasing blast radius:

- **(a) Knob on the merge** — extend `root_pcs_to_merge` (or add a
  `merge_per_apa_ctpc` flag) on **QLMatching** (and `PointTreeMerging` for the
  `--per-apa` path), defaulting to today's behavior; enable it only for the
  SBND PR job's tree. Byte-identical when off.
- **(b) Rebuild ctpc in the PR job** from the loaded blobs instead of relying on
  the merged root PCs. Leaves the gated stage untouched entirely, at the cost
  of recomputation.
- **(c) Leave as-is and make the failure loud** — at minimum, `cluster_fc_check`
  returning the conservative `is_fc=false` because the steiner graph is missing
  should be distinguishable in the table from a genuine "boundary exits found"
  verdict. Right now both surface as `fc=0`.

Independent of the choice, `get_overlap_good_ch_charge` silently returning an
empty map for a missing dataset is what let this hide; a warning there would
have surfaced it immediately.

---

## 7. Second manifestation — the boundary wcps collapse to the ORIGIN

Investigated because SBND MCP2025C evt **285999** bundle grp 3 (main cluster
**18**, flash t0 = −891.5 µs) is not labeled FC either. Unlike evt 285185 main
12 this cluster *does* have a steiner graph (it straddles the cathode, so its
TPC0 half still carries charge), and it is not in the
"produced no steiner_graph" warning list. The same missing `ctpc_a1f0p*`
breaks it a second, sharper way.

### Repro

```bash
wcbuild
SB=/nfs/data/1/xqian/toolkit-dev/wcp-porting-img/sbnd/sbnd_xin
cd $SB && SBND_INPUT_DIR=$PWD/input_files_reco1/extracted-mcp2025c-10evt \
  SBND_WORK_ROOT=$PWD/work-mcp10-mainflag ./run_nusel_evt.sh data 4 -chord
# ctpc datasets in this event's tree (APA0 only):
#   tar xzf work-mcp10-mainflag/ql_evt285999/pctree-evt285999.tar.gz '*metadata.json'
#   -> ctpc_a0f0p{U,V,W}, dead_winds_a0f0p*, and NO ctpc_a1f0p*
```

Numbers below come from a temporary `WCT_FC_DEBUG=<ident>` probe in
`cluster_fc_check`, `FiducialUtils::check_signal_processing` and
`Cluster::get_two_boundary_wcps`. Reverted; the tree is clean.

### The cluster

1833 points, 204.7 cm, a clean cathode-crossing cosmic:
TPC0 843 pts x ∈ [−92.2, −4.1], TPC1 990 pts x ∈ [+4.4, +92.6],
y 20.8 → 107.9, z 247.9 → 260.9. Every point is far inside the PR fiducial box
(`sbnd_pr_fv` ± margins: |x| ≤ 199.05, |y| ≤ 196.8, 3.85 ≤ z ≤ 497.15), so no
direct containment test can fail.

### The chain

1. `Blob::estimate_total_charge()` sums `Grouping::get_wire_charge()`, which
   reads the per-APA ctpc through `build_wire_cache` — and that function
   silently leaves the cache empty when the dataset is absent
   (`Facade_Grouping.cxx:733`). Probe:

   ```
   FCDBG b2wcps clus 18 blobs apa0 302 (sum est charge 5064693) apa1 337 (sum est charge 0)
   ```

2. `Cluster::get_two_boundary_wcps()` skips every blob with
   `estimate_total_charge() < 1500` (`Facade_Cluster.cxx:3455`). With APA1
   charge identically zero **all 990 TPC1 points are skipped**, so the
   `(apa,face) = (1,0)` entry — which was pre-seeded from `wpids_blob()` —
   stays `initialized == false` and keeps its default-constructed extreme
   points:

   ```
   FCDBG b2wcps clus 18 apa 0 face 0 initialized true  -> (-4.1,63.8,254.8) - (-92.2,20.8,248.0)
   FCDBG b2wcps clus 18 apa 1 face 0 initialized false -> (0.0,0.0,0.0) - (0.0,0.0,0.0)
   ```

3. The function then returns the **farthest-apart pair** among all per-key
   boundary points. The origin is ~250 cm from a track sitting at z ≈ 250, so
   the origin wins:

   ```
   FCDBG r1 PRE-SNAP boundary wcps (-92.2,20.8,248.0) - (0.0,0.0,0.0)
   ```

   This is a **silent uninitialized-value leak**, not a physics choice: (0,0,0)
   is outside the detector entirely.

4. `get_two_boundary_steiner_graph_idx` snaps each boundary wcp to the nearest
   steiner terminal. The steiner cloud is TPC0-only (271 pts, x ∈ [−91.9,
   −4.1] — §3's degradation), and the terminal nearest the origin is an
   **interior** point 51 cm along the TPC0 half:

   ```
   FCDBG r1 steiner x range [-91.9,-4.1] cm; boundary (-91.9,21.0,247.9) - (-41.0,45.0,251.3)
   ```

5. `cluster_fc_check` tests that interior point as if it were a track endpoint.
   The local Hough direction there points *along* the track, the W-plane
   projected angle is 3.8° (< 5°), so `check_signal_processing` runs and walks
   outward straight through the cluster's own charge until it leaves
   FiducialUtils' volume at the CPA slab (x = +0.2 cm):

   ```
   FCDBG SP walk from (-41.0,45.0,251.3) dir (0.92,0.40,0.06):
         left FV after 45 steps at (0.2,62.9,254.1), dead 42 -> ret false
   ```

   42 > 0.8 × 45 ⇒ "charge continues past this endpoint" ⇒ exit.

6. `exit_wcps = {(-41.0,45.0,251.3)}` ⇒ `is_fc = false`.

Both **real** endpoints pass every test — (+92.6,107.8,260.9): SP not
triggered, dead-volume OK; (−92.2,20.8,248.0): SP OK (early-true after 7
steps), dead-volume OK. So absent the origin leak this cluster would be
**FC = true** — an inference, not an observation: round 2 (`flag_cosmic=false`)
never ran here because round 1 already produced an exit, and no leak-free
variant was executed.

### Scope of the leak

`get_two_boundary_wcps` is not FC-only. Its callers are
`CreateSteinerGraph.cxx:142`, `improvecluster_2.cxx:110,151`,
`TaggerCheckSTM.cxx:1766`, `NeutrinoPatternBase.cxx:1717` and
`cluster_fc_check`. Any cluster with APA1 blobs in the merged tree therefore
hands the same origin point to steiner-tree seeding, retiling, STM and the
neutrino PR — not just FC. (Call sites and the returned origin are verified;
no corrupted STM/PR verdict was observed here — that is unmeasured.) On the 10-event sample that is every cluster with
TPC1 content (§4's "mixed" and "TPC1" buckets).

### The physics answer for this bundle

**Owner ruling (2026-07-23): FC is a fiducial-containment verdict and is
independent of the readout-window question — evt285999 cluster 18 should be
FC = true.** Every one of its reconstructed points is deep inside
`sbnd_pr_fv`, and both real endpoints pass all three tests; the only thing
standing between this bundle and FC = true is the leaked origin point (and,
behind it, the missing APA1 charge that produces it). The fact that the same
track is *also* an out-of-time cosmic truncated by the readout window
(ends at |x| = 92.2 / 92.6, symmetric about the cathode; the sibling
cluster 20 at t0 = −1175.9 µs truncates at |x| ≈ 50.3 and is labeled FC = 1)
is a separate concern tracked elsewhere and must not be folded into the FC
verdict.

## 8. Direct proof of where the APA1 charge goes

A temporary `WCT_MERGE_DEBUG` probe around `merge_pct`
(`match/src/QLMatching.cxx:1052`, joint path) on evt285999, printing each
input root's local PC names before the merge and the merged root's names
after (probe reverted):

```
MERGEDBG input[0] root_live local_pcs: ctpc_a0f0pU ctpc_a0f0pV ctpc_a0f0pW
                                       dead_gap_a0f0pW dead_winds_a0f0pU dead_winds_a0f0pV
                                       dead_winds_a0f0pW flash flashlight light opflash
MERGEDBG input[1] root_live local_pcs: ctpc_a1f0pU ctpc_a1f0pV ctpc_a1f0pW
                                       dead_gap_a1f0pW dead_winds_a1f0pU dead_winds_a1f0pV
                                       dead_winds_a1f0pW flash flashlight light opflash
MERGEDBG merged  root_live local_pcs: ctpc_a0f0pU ctpc_a0f0pV ctpc_a0f0pW
                                       dead_gap_a0f0pW dead_winds_a0f0pU dead_winds_a0f0pV
                                       dead_winds_a0f0pW flash flashlight light opflash
```

So the APA1 2D charge **is** produced upstream, correctly named, and present on
input[1] right up to the merge. `merge_pct` copies only the names listed in
`root_pcs_to_merge` (SBND: `['opflash']`) from every source root and then
`take_children`s the blob/cluster nodes across. The APA1 *blobs* therefore
arrive in the merged tree, but the 2D charge and dead-channel arrays that
describe them stay on the source root and die with it. Nothing upstream is
broken — this is purely the merge policy. The same holds for
`dead_gap_a1f0pW` / `dead_winds_a1f0p*`, so APA1 dead-channel queries
(`is_wire_dead`, `inside_dead_region`, the dead-gap column) are blind too.

**Cost of keeping them** (evt285999, measured on the persisted pctree):
`ctpc_a0*` + `dead_winds_a0*` + `dead_gap_a0*` are 2.4 MB of the 9.3 MB of
array payload (26 %); adding APA1's would grow the tree by roughly the same
amount (tarball 2.6 MB → ~3.3 MB). Modest, not the order-of-magnitude cost
§6.3 worried about.

---

## 9. FIX (2026-07-23) — merge the per-anode PCs; guard the boundary scan

Owner decision: clear bug, fixed unconditionally (no knob), same precedent as
`d48dc647`. Three edits:

| edit | file |
|---|---|
| `merge_pct` also carries `ctpc_a*`, `dead_winds_a*`, `dead_gap_a*` | `match/src/QLMatching.cxx` (joint path) |
| same, for the standalone merger | `clus/src/PointTreeMerging.cxx` |
| skip never-initialized (apa,face) keys in the farthest-pair scan | `clus/src/Facade_Cluster.cxx::get_two_boundary_wcps` |

The names are unique per (apa,face,plane) so the merge is collision-free;
`flash`/`light`/`flashlight` stay dropped as before (per-APA duplicates of the
same physical flashes). The boundary guard stays necessary independently: any
key whose blobs all fall below the 1500 charge cut would still leak (0,0,0).

### Repro

```bash
wcbuild
SB=/nfs/data/1/xqian/toolkit-dev/wcp-porting-img/sbnd/sbnd_xin
cd $SB && SBND_MAX_JOBS=4 SBND_INPUT_DIR=$PWD/input_files_reco1/extracted-mcp2025c-10evt \
  SBND_WORK_ROOT=<fresh root, imaging symlinked from work-mcp10-mainflag/evt*> \
  ./run_nusel_evt.sh data all -chord
```

### Verification (freshness: libs 19:34 > sources 19:27; wcdoctest-clus 41/41 518 asserts, wcdoctest-match 4/36)

**SBND, 10-evt MCP2025C (post-fix scratch vs `work-mcp10-mainflag`):**

- `mabc-all-apa.zip` **byte-identical 10/10** — QL matching and the all-APA
  clustering stage are untouched on SBND.
- pctree now carries `ctpc_a1f0p{U,V,W}`, `dead_winds_a1f0p*`,
  `dead_gap_a1f0pW` (tarball +9%, 2.58→2.81 MB on evt285999).
- Steiner-graph failures 42 → 12; every survivor is a ≤28-point speck
  (the legitimately-too-small class). The 32/32 TPC1-only starvation is gone.
- **evt285999 main 18 → FC=1.** evt285185 main 12 → FC=1, and its 480 cm
  main 16 is now tagged TGM. All verdict flips, by column:
  - FC 0→1: 284657 c13, 285185 c12, 285999 c18, 286527 c17.
  - new tags: 284657 c14 →STM, 285185 c16 →TGM, 286241 c8 →STM, c2 →TGM.
  - one tag removed: 286021 c15 STM 1→0 — itself an origin-leak victim
    (8-point TPC1 speck at x≈116 zeroed its apa1 key pre-fix, so its old
    STM=1 came from a corrupted boundary pair).
  - five -1↔0 cells are torn log lines (verified directly: the FC verdict
    text is interrupted mid-line by an interleaved write), not verdicts.

**PDHD+PDVD (`abtest` clus gate, labels `pre_ctpcmerge`/`post_ctpcmerge`,
5 events — pdhd 028084/18 excluded, its work dir has no SP frames):**

- Every per-anode / per-face zip **PASS** (pre-merge stages untouched).
- Merged-scope zips differ on 3/5 events — the expected consequence:
  group/all-APA `connect1`/`deghost`/`separate` now see the non-primary APAs'
  dead-channel maps (previously silently empty). Magnitude:
  - pdhd_027305_0: −20 pts (of 355k), clusters 116→116
  - pdvd_039349_0: −3 pts (of 44k), clusters 77→77
  - pdvd_039252_5: −279 pts (of 179k), clusters 120→**121**
  - pdhd_027409_0, pdhd_027980_3: fully byte-identical.
  **NOT bit-identical — pdhd/pdvd downstream needs revalidation** (Bee spot
  checks of the changed events recommended).

**uBooNE (`qlport` sweep `pre_ctpcmerge`/`post_ctpcmerge`, 35 events):**

- Gate 1: 34/35 Bee zips content-identical. The 35th (`mabc_21`,
  `shower_track` layer) is a pure permutation (same 11613 points, same total
  charge) and **reruns of the post-fix binary reproduce the PRE hash**
  (labels `detcheck1`/`detcheck2`) — a nondeterministic flip, not the fix.
- Gate 2: 32/35 tagger logs differ — but the SAME binary run twice differs
  identically (idx12: `kine_energy_particle` flips 805.621↔833.06, per-track
  arrays permute). Pre-existing PR-tail run-to-run nondeterminism, NOT this
  change. Worth its own investigation: a 28 MeV kine flip is more than the
  documented "benign FP tie + cosmetic gidx relabel" residual.

### Left open

- The PR-tail nondeterminism above (uBooNE repeat_check-style).
- `work-mcp10-mainflag` / `work-mcp1000-mainflag` (the live scan roots on
  port 5010) still hold pre-fix FC/STM/TGM columns.
- pdhd 028084/18 needs its NF+SP frames regenerated before it can gate again.
