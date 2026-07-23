# APA1 has no 2D charge in the PR tail: every TPC1-only cluster loses its steiner graph

**Status: FINDING ONLY — no code changed.** This needs an owner decision (see
§6): the natural fix touches `PointTreeMerging`, which runs inside the
byte-identical-gated all-APA clustering stage.

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
