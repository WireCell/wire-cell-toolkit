# QLMatching — current code reference (input / algorithm / output / BEE dump)

Orientation for **modifying** the SBND charge–light (Q/L) matching code. This is
a code/dataflow reference for the implementation *as it currently is*. For the
larsoft→WCT porting history, JSON schema and build details see
[`qlmatching-port.md`](qlmatching-port.md) and
[`porting-summary.md`](porting-summary.md). For how this component is wired into
a runnable pipeline and packaged for BEE, see the chain doc
`wcp-porting-img/sbnd/sbnd_xin/docs/ql-chain.md`.

All line references are into `match/src/QLMatching.cxx`, `match/src/Util.cxx`,
etc., as of the `apply-pointcloud` branch.

---

## 1. Component shape

`WireCell::Match::QLMatching` is an `ITensorSetFilter` (1→1) + `IConfigurable`
(`QLMatching.h:30`). It is registered as factory `"QLMatching"`.

```
operator()(const input_pointer& in, output_pointer& out)   // QLMatching.cxx
```

| Port | Direction | Contents |
|------|-----------|----------|
| in 0 | cluster tensorset | a point-cloud tree at `inpath` (`pointtrees/<id>`) with `/live` and `/dead` groupings — the imaging/clustering result for this APA, **with the per-event optical flashes attached as a `flash` point cloud on the live root node** (placed there by `Aux::AttachPointCloudToTree`, §1a) |
| out 0 | cluster tensorset | **the input cluster tensorset passed through**, with only each matched cluster's `cluster_t0` mutated to the matched flash time |

It is a *filter* (1→1): charge+light in, charge-with-t0 out. The light is no
longer a separate input port — it rides on the cluster tree's live root node
(mirroring the MicroBooNE `UbooneClusterSource` design, where optical data is
placed on the root node and the matcher reads it from there). The matching
result for the event display is written as a **side effect** to BEE JSON files
(§5), not into the output tensorset.

## 1a. Light I/O — two standard nodes (read + attach)

The light path is **two reused graph nodes**, inserted between clustering and
QLMatching:

1. **`Sio::TensorFileSource`** (existing, plugin `WireCellSio`) — reads the
   opflash tensor archive (`inname: "opflash_apa<n>.tar.gz"`, `prefix:
   "opflash_"`) → an opflash tensor set on its output port.
2. **`Aux::AttachPointCloudToTree`** (`aux/src/AttachPointCloudToTree.cxx`,
   plugin `WireCellAux`) — a generic `ITensorSetFanin` (2→1): port 0 = the
   cluster pctree tensor set, port 1 = the opflash tensor set. It deserializes
   the live tree (`as_pctree(.../live)`), stores port-1's first tensor verbatim
   as a 2-D `value` array in a point cloud named `pcname` on the **live root
   node** (`root_live->value.local_pcs()[pcname]`), re-serializes live, and
   passes `/dead` through. Config: `pcname` (required), `inpath` (default
   `"pointtrees/%d"`).

`AttachPointCloudToTree` does **not** interpret the tensor — it just parks it on
the root node for a downstream consumer. The SBND chain configures `pcname:
"flash"` so it provides the `[nflash, 1+nchan]` matrix QLMatching expects
(`flash_pcname` must match `pcname`).

(This replaced an earlier single `Sio::TensorFileToPCTree` that *internally*
constructed and pumped a `TensorFileSource` — the two-node split reuses
`TensorFileSource` as-is via standard graph wiring, and the framework syncs the
two input ports / manages the source's lifecycle.)

---

## 2. Configuration (`configure()`, `:36-116`)

| Key | Default | Member | Meaning |
|-----|---------|--------|---------|
| `anode` | (req) | `m_anode` | `IAnodePlane` — TPC id / drift sign |
| `detector_volumes` | (req) | `m_dv` | `IDetectorVolumes` — drift & Y/Z bounds |
| `inpath` / `outpath` | `pointtrees/%d` | `m_inpath`/`m_outpath` | pctree path template |
| `flash_pcname` | `flash` | `m_flash_pcname` | name of the live-root PC holding the flash matrix (must match `AttachPointCloudToTree.pcname`) |
| `bee_dir` | `data` | `m_bee_dir` | BEE-dump output dir (empty ⇒ no dump) |
| `semimodel_file` | `sbnd/photodet/semi-analytical-sbnd.json` | `m_semimodel_file` | photon model JSON (`Persist::load`) |
| `pmts` | `true` | `m_pmts` | use the SBND 312-OpDet PMT mask |
| `data` | `true` | `m_data` | data vs MC; MC masks saturated PMTs (`:258`) |
| `beamonly` | `false` | `m_beamonly` | force beam window + skip strength cut |
| `ch_mask` | `[]` | `m_ch_mask` | OpDet indices to disable |
| `flash_minPE` | `500` | `m_flash_minPE` | min total PE to keep a flash |
| `flash_min/maxtime` | ∓1500 ms | … | flash time window (overridden if `beamonly`) |
| `beam_min/maxtime` | ∓5 ms | … | beam window |
| `QtoL` | `0.5` | `m_QtoL` | charge→light scale in the prediction (`:317`) |
| `strength_cutoff` | `0.05` | `m_strength_cutoff` | LASSO solution threshold to keep a bundle |
| `VUVEfficiency` / `VISEfficiency` | 312-elt arrays | … | per-OpDet QE (direct / reflected) |

**Standalone jsonnet overrides** (`wct-clus-matching-standalone.jsonnet`): sets
`flash_minPE: 50` (not 500), `QtoL: 1.0` (not 0.5), `data` from `reality`,
`bee_dir: "data-sep"`, and a 19-element `ch_mask`. Keep these in mind when
comparing code defaults to observed runs.

The `semimodel_file` JSON top-level keys `VUVHits`, `VISHits`, `Geometry`,
`OpDets` are handed straight to the `SemiAnalyticalModel` ctor (`:78-112`).

---

## 3. Input parsing (`operator()`)

- **EOS** guard first: a null `in` returns immediately (single input now).
- **OpDet mask**: start from the hardcoded SBND PMT mask if `m_pmts`, then zero
  out `m_ch_mask` entries; later reduced per-TPC (even/odd OpDets) and, for MC,
  saturated channels (`total_PE>5000 && PE==0`).
- **Clusters**: the input tensorset is read as a pctree
  (`as_pctree(charge_tens, inpath + "/live")`), wrapped as a `Facade::Grouping`,
  anode + detector-volumes attached, and its `children()` taken as the
  `Cluster*` list, **sorted by length descending**.
- **Flashes**: read the `flash_pcname` (`flash`) point cloud from the live root
  node (`root_live->value.local_pcs()`), take its 2-D `value` array
  `[nflash, 1+nchan]`, rebuild an `Aux::SimpleTensor` over it, and construct
  `Opflash` objects per row exactly as before
  (`std::make_shared<Opflash>(ten, iflash, 0.0, nchan)`); each is dropped unless
  its time is in `[flash_mintime, flash_maxtime]` **and** `total_PE >=
  flash_minPE`. (Empty/absent flash PC ⇒ 0 flashes for that event.)

---

## 4. Algorithm (`operator()`)

### 4.1 Build candidate bundles (`:231-352`)
For every (flash × cluster) pair create a `TimingTPCBundle(flash, cluster,
flash_id, cluster_idx)` (`:269`). For each 3-D point in the cluster:
- drift-correct the X by the flash time: `x += sign * flash_time * 1.563e-3`
  (`:252,289`) — **note the drift speed is hardcoded `1.563e-3`**, independent of
  the jsonnet `driftSpeed`; a refactor target.
- drop points outside the drift / Y / Z bounds; set `close_to_PMT` /
  `at_x_boundary` flags (`:293-299`).
- predict light: convert point to cm, call
  `m_semi_model->detectedDirectVisibilities()` and
  `detectedReflectedVisibilities()` (`:306,308`), accumulate per OpDet
  `pred_flash[idet] += q * QtoL * (dir_vis*VUVeff + ref_vis*VISeff)` (`:310-317`).

Then quality-gate the bundle (`:323-337`): skip if >25% of points drifted out of
bounds; `set_pred_flash`; skip if `get_total_pred_light() < 10`;
`examine_bundle()` (computes KS distance + chi²/ndf); skip if `ks_dis == 1` or
`chi2/ndf > 1e4`. Survivors go into `pre_bundles`.

### 4.2 Deterministic ordering (`:355-410`)
Bundles are bucketed into `flash_bundles_map`, `cluster_bundles_map`, and
`flash_cluster_bundles_map`, then iterated in a **stable order** — flashes by
`flash_id`, clusters by global cluster index, inner bundles by cluster index
(`:384-410`). This is the `6eb0de08` determinism fix (replacing
allocator-dependent `std::map<Pointer*>` iteration that flipped marginal
bundles run-to-run). Preserve this ordering in any refactor.

### 4.3 Two LASSO solves (`Ress::lasso`)
The match is posed as a sparse linear solve `min ‖y − X·s‖ + λ‖s‖₁` with per-row
weights; the solution `s` is the per-bundle "strength".

- **Round 1** (`:442-522`): design matrix has bundle columns **plus per-flash
  background columns**, and rows for both per-(opdet,flash) light and per-cluster
  charge constraints. After solving, bundles with `solution ≤ strength_cutoff`
  are pruned (unless `beamonly`) via `remove_bundle_selection`.
- **Round 2** (`:524-634`): rebuilt over the survivors; **bundle columns only**,
  weights include the KS-test penalty (`:572`). `bundle->set_strength(solution)`
  (`:597`), prune again by `strength_cutoff`.

These are sequential phases (rough cut → final weighted selection), not
independent solves.

### 4.4 Selection & organization (`:607-634`)
`matched_pairs` keeps, per cluster, the flash with the highest strength
(`:607-617`). `results_bundles` is built per cluster (`:618-633`): the matched
bundle, **or** — for an unmatched cluster — a placeholder
`std::make_shared<TimingTPCBundle>(nullptr, cluster, 0, cidx)` (`:629`).

> **Null-flash bundle (fixed):** that `nullptr` flash used to crash the
> `TimingTPCBundle` ctor, which did `flash->get_num_channels()`
> (`TimingTPCBundle.cxx:50`) → null-deref. MC never exercised it (every cluster
> matched a flash); `data` mode does (cosmic clusters with no flash). The ctor
> now guards: `m_nchan = flash ? flash->get_num_channels() : 0;` — a null-flash
> bundle carries no predicted light and its `pred_flash` is never read
> (`organize_bundles` already skips null-flash bundles).

`organize_bundles()` (`:735-806`) then merges compatible bundles per flash and
applies a beam-window quality filter (drop out-of-beam bundles with `ks>0.2`,
`chi2/ndf>20`, or PE mismatch >50%).

### 4.5 t0 + output (`:656-680`)
For every bundle in `flash_bundles_map`, set
`main_cluster->set_cluster_t0(flash_time)` (`:659`). Then rebuild the output
tensorset from the (now t0-stamped) live + dead pctrees and emit it under the
same `charge_ident`.

### Branch effects
- `m_data == false` (MC): masks saturated PMTs (`:255-260`).
- `m_beamonly == true`: forces the flash window to the beam window (`:61-64`) and
  **skips** the `strength_cutoff` pruning in both rounds (`:514,598`).

---

## 5. BEE output — how the event-display JSON is dumped

When `bee_dir` is non-empty, after matching (`:642-652`):

```cpp
sub_dir = "<bee_dir>/<bee_index>";              // e.g. data-sep/0
dump_bee_3d (*root_live, "<sub_dir>/<bee_index>-img-apa<tpc>.json");
dump_light  (flashes, flash_bundles_map, global_cluster_idx_map,
             "<sub_dir>/<bee_index>-op-apa<tpc>.json");
++m_bee_index;                                  // one event per call
```

So each call writes `data-sep/<n>/<n>-img-apa<tpc>.json` and
`…-op-apa<tpc>.json`, where `<tpc> = m_anode->ident()` and `<n>` increments per
event. These two files are the **charge** (img) and **light/matching** (op)
layers of one APA for one event.

### `dump_bee_3d` — the charge/img layer (`Util.cxx:22-88`)
Per-point 3-D cloud of the live clusters. Fields: `runNo/subRunNo/eventNo`
(all 0), `geom:"uboone"`, `type:"cluster"`, `x`,`y`,`z` (cm), `q` (1.0 per
point), `cluster_id` (sequential id per cluster, clusters sorted by length).

### `dump_light` — the light/matching layer (`Util.cxx:188-254`)
One entry per flash (the `f8b91803` change: **all** flashes, matched or not).
Fields: `runNo/subRunNo/eventNo` (0), `geom:"sbnd"`,
`op_nomatching_cluster_ids` (placeholder), and per flash:

| field | meaning |
|-------|---------|
| `op_t` | flash time (µs) |
| `op_pes` | measured PE per channel (`flash->get_PEs()`) |
| `op_pes_pred` | predicted PE per channel from the matched bundle (`bundle->get_pred_flash()`); **`[]` if unmatched** |
| `op_peTotal` | total measured PE |
| `cluster_id` | `[cid]` of the matched cluster; **`[]` if unmatched** |

Matched vs unmatched is therefore visible directly: an empty `cluster_id` /
`op_pes_pred` ⇒ a flash with no charge match (still drawn in the BEE light
display). A bundle only counts as "matched" output if its predicted PE total is
≥ 100 (`:164`-style filter inside the dump).

### Related dump helpers (not on the current path)
- `dump_bee_flash` (`Util.cxx:90-134`) — raw flashes from a tensorset, `op_pes_pred = op_pes` (no matching). Unused by `operator()`.
- `dump_bee_bundle` (`Util.cxx:136-186`) — matched bundles **only**; superseded by `dump_light`. Unused by `operator()`.

### From JSON to a BEE link
The per-APA `…-img/op-apa<tpc>.json` files are merged across APAs by
`merge-apa.py` and unioned with the clustering zips by `bee-upload.sh` into
`combined.zip`, whose upload returns the BEE URL. That packaging half lives in
the chain doc (`sbnd_xin/docs/ql-chain.md`); here the contract is just the file
naming + fields above.

---

## 6. Helper classes

- **`Opflash`** (`Opflash.{h,cxx}`) — one optical flash. Holds per-channel `PE`
  (+ `PE_err`), `total_PE`, `time` (ns), `m_nchan`, `flash_id`, fired-channel
  list. Key accessors: `get_PEs()`, `get_PE(ch)`, `get_total_PE()`,
  `get_time()`, `get_num_channels()`, `get_flash_id()`.
- **`TimingTPCBundle`** (`TimingTPCBundle.{h,cxx}`) — one (flash, cluster)
  candidate. Holds `flash`, `main_cluster`, `pred_flash` (predicted PE/ch),
  `ks_dis`/`chi2`/`ndf`, `strength`, opdet mask, and flags
  (`close_to_PMT`, `at_x_boundary`, `potential_bad_match`, `high_consistent`).
  `examine_bundle()` (`:148-192`) fills the KS/chi² metrics. Ctor
  (`:33-56`) takes `m_nchan` from the flash but guards a null flash
  (unmatched cluster) — see §4.4.
- **`SemiAnalyticalModel`** (`SemiAnalyticalModel.{h,cxx}`) — point → per-OpDet
  direct (VUV) + reflected (VIS) visibilities. SBND-minimal scope (dome PMTs +
  flat (X)Arapucas; no lateral PDs / anode reflections / Xe). Details in
  [`qlmatching-port.md`](qlmatching-port.md).

---

## 7. Refactor-relevant notes

- **Hardcoded drift speed** `1.563e-3` in the per-point X correction
  (`:252,289`) — not tied to the jsonnet `driftSpeed`. Thread it through config.
- **No `tagger_info` PC written.** Downstream all-APA tagging
  (`ClusteringTaggerFlagTransfer`) expects a per-cluster `tagger_info` PC
  (`has_beam_flash`, …) to set the `beam_flash` flag; QLMatching writes none, so
  with neutrino tagging enabled `recover_bundle` finds nothing → null
  `main_cluster` → segfault in `TaggerCheckNeutrino`. This is why the chain runs
  with `nu_tagging=false` (see chain doc). Writing `tagger_info` here is the
  proper fix.
- **Two-solve coupling** — round 2 depends on round-1 pruning; keep the
  deterministic ordering (§4.2) if you touch either.
- **`geom:"uboone"`** is hardcoded in the img dump (`Util.cxx:32`) while the op
  dump uses `"sbnd"`; harmless to BEE but worth normalizing.
- **Output carries only `cluster_t0`** — if downstream needs match metadata
  (matched flash id, strength, beam flag), it must be added to the output
  tensorset / a `tagger_info` PC; today it lives only in the BEE JSON.
