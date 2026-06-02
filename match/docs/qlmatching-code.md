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
| in 0 | cluster tensorset | a point-cloud tree at `inpath` (`pointtrees/<id>`) with `/live` and `/dead` groupings — the imaging/clustering result for this APA, **with the per-event optical flashes attached as the canonical `flash`/`light`/`flashlight` point clouds on the live root node** (placed there by `Aux::FlashTensorToOpticalPCs`, §1a) |
| out 0 | cluster tensorset | **the input cluster tensorset passed through**, with each matched cluster's `cluster_t0` set to the matched flash time **and a per-cluster `flash` scalar set to the matched flash row index** (so `Clus::Facade::Cluster::get_flash()` reflects the match); unmatched clusters get `flash = -1`. **Additionally persisted for the downstream Bee op/flash display** (see §5): a per-cluster `matched_flash_gid` scalar (a globally-unique flash id `anode·1e6 + per-APA-flash-index`, surviving the per-APA→all-APA merge; `-1` if unmatched), a per-cluster `flashpred` PC holding the matched bundle's predicted per-channel PE (`op_pes_pred`), and a self-contained per-root `opflash` PC (one row per `(flash, channel)`: `gid`, `time`, `ch`, `pe`) carrying every flash's measured light |

It is a *filter* (1→1): charge+light in, charge-with-t0 out. The light is no
longer a separate input port — it rides on the cluster tree's live root node in
the **same canonical schema as the MicroBooNE `UbooneClusterSource`** (optical
data on the root node, per-cluster matched-flash scalar). The matching result
for the event display is also written as a side effect to BEE JSON files (§5).

## 1a. Light I/O — two standard nodes (read + expand to canonical PCs)

The light path is **two graph nodes**, inserted between clustering and QLMatching:

1. **`Sio::TensorFileSource`** (existing, plugin `WireCellSio`) — reads the
   opflash tensor archive (`inname: "opflash_apa<n>.tar.gz"`, `prefix:
   "opflash_"`) → an opflash *matrix* tensor set `[nflash, 1+nchan]` (col 0 =
   time, cols 1..nchan = per-channel PE) on its output port.
2. **`Aux::FlashTensorToOpticalPCs`** (`aux/src/FlashTensorToOpticalPCs.cxx`, plugin
   `WireCellAux`) — an `ITensorSetFanin` (2→1): port 0 = the cluster pctree
   tensor set, port 1 = the opflash matrix. It deserializes the live tree
   (`as_pctree(.../live)`) and **expands the matrix into the canonical
   `flash`/`light`/`flashlight` point clouds** on the **live root node**
   (`root_live->value.local_pcs()`), re-serializes live, passes `/dead` through.
   Config: `nchan` (default `312`; raises if `ncol-1 != nchan`), `inpath`
   (default `"pointtrees/%d"`).

This is the **same schema** the MicroBooNE `root/UbooneClusterSource` writes, so
SBND flashes interoperate with all `clus` tooling (`get_flash()`,
`ClusterFlashDump`, `retile`, BEE):
- `flash`: `time`, `tmin`, `tmax`, `value`(=ΣPE), `ident`(=row), `type`(=0)
- `light` (sparse, one row per fired channel): `ident`(=channel), `time`,
  `value`(=PE), `error`(=0)
- `flashlight`: `flash`(=flash row), `light`(=light row) — the join table

**Units:** the stored `time` is in WCT units (ns) — the **same** as Uboone
`load_optical` — but taken from the matrix verbatim (no `units::us` scaling)
because the SBND dump is already ns (`ns == 1`) whereas Uboone's tree is µs and
gets scaled. Same canonical units, different native input unit. Per-channel
`error` is `0` (the SBND dump carries none — the `PE_err` 0.3 rule stays inside
`Match::Opflash`).

`FlashTensorToOpticalPCs` carries no detector physics (it only reshapes the matrix), so
it lives in **`aux/` alongside `Aux::AttachPointCloudToTree`** — the earlier generic
"park a matrix as one PC" node, which remains in `aux/` as a generic primitive but is
no longer used by this chain (`FlashTensorToOpticalPCs` supersedes it for flashes). The WCT
factory type string is `"FlashTensorToOpticalPCs"` (unchanged by the package move).

---

## 2. Configuration (`configure()`, `:36-116`)

| Key | Default | Member | Meaning |
|-----|---------|--------|---------|
| `anode` | (req) | `m_anode` | `IAnodePlane` — TPC id / drift sign |
| `detector_volumes` | (req) | `m_dv` | `IDetectorVolumes` — drift & Y/Z bounds |
| `inpath` / `outpath` | `pointtrees/%d` | `m_inpath`/`m_outpath` | pctree path template |
| `nchan` | `312` | `m_nchan` | optical-detector channel count (per-flash PE vector length); the bit-safe source of nchan now that flashes arrive via canonical PCs |
| `bee_dir` | `data` | `m_bee_dir` | BEE-dump output dir (empty ⇒ no dump) |
| `semimodel_file` | `sbnd/photodet/semi-analytical-sbnd.json` | `m_semimodel_file` | photon model JSON (`Persist::load`) |
| `pmts` | `true` | `m_pmts` | apply the OpDet type mask (gate) |
| `active_opdet_types` | `[1]` | `m_active_opdet_types` | OpDet `type`s kept by the mask (1=PMT, 0=(X)Arapuca); derived from the injected `OpDets` table |
| `data` | `true` | `m_data` | data vs MC; MC masks saturated PMTs (`:258`) |
| `beamonly` | `false` | `m_beamonly` | force beam window + skip strength cut |
| `ch_mask` | `[]` | `m_ch_mask` | OpDet indices to disable |
| `flash_minPE` | `500` | `m_flash_minPE` | min total PE to keep a flash |
| `flash_min/maxtime` | ∓1500 ms | … | flash time window (overridden if `beamonly`) |
| `beam_min/maxtime` | ∓5 ms | … | beam window |
| `QtoL` | `0.5` | `m_QtoL` | charge→light scale in the prediction (`:317`) |
| `strength_cutoff` | `0.05` | `m_strength_cutoff` | LASSO solution threshold to keep a bundle |
| `drift_speed` | `1.563e-3` | `m_drift_speed` | drift speed for the per-flash X correction, in WCT units (mm/ns). Pass `params.lar.drift_speed` from the common config. |
| `VUVEfficiency` / `VISEfficiency` | 312-elt arrays | … | per-OpDet QE (direct / reflected) |
| `anode_ext1` | `-2.0 cm` | `m_anode_ext1` | anode-side inclusion / flag-window edge in `u` (`anode_in`); see §4.1a |
| `anode_ext2` | `4.0 cm` | `m_anode_ext2` | anode flag-window outer edge in `u` (close-to-PMT / x-boundary) |
| `cathode_ext1` | `1.2 cm` | `m_cathode_ext1` | cathode-side inclusion edge (`cathode_in = u_cathode + ext1`) |
| `cathode_ext2` | `-2.0 cm` | `m_cathode_ext2` | cathode flag-window inner edge (x-boundary) |
| `window_edge_ticks` | `4` (SBND: **24**) | `m_window_edge_ticks` | raw-tick edge band for `window_truncated`; SBND jsonnet sets 24 |
| `readout_window_ticks` | `3427` (SBND: **3428**) | `m_readout_window_ticks` | exclusive readout end for the trailing-edge truncation test; SBND jsonnet sets 3428 (rebin-4 `slice_index_max`) |
| `require_containment` | `false` (SBND: **true**) | `m_require_containment` | when true, discard bundles whose cluster fails the TPC-box containment guard (prototype `flag_good_bundle`); see §4.1a. Default OFF keeps non-SBND configs bit-identical |

**Standalone jsonnet overrides** (`wct-clus-matching-standalone.jsonnet`): sets
`flash_minPE: 50` (not 500), `QtoL: 1.0` (not 0.5), `data` from `reality`,
`bee_dir: "data-sep"`, and a 19-element `ch_mask`. Keep these in mind when
comparing code defaults to observed runs.

The `semimodel_file` JSON top-level keys `VUVHits`, `VISHits`, `Geometry`,
`OpDets` are handed straight to the `SemiAnalyticalModel` ctor (`:78-112`).

---

## 3. Input parsing (`operator()`)

- **EOS** guard first: a null `in` returns immediately (single input now).
- **OpDet mask**: if `m_pmts`, build per-channel from the injected `OpDets` table
  (on iff `type ∈ active_opdet_types`, default `[1]` = PMTs), then zero out
  `m_ch_mask` entries; later reduced per-TPC (by OpDet `center.x` vs `cathode_x`)
  and, for MC, saturated channels (`total_PE>5000 && PE==0`).
- **Clusters**: the input tensorset is read as a pctree
  (`as_pctree(charge_tens, inpath + "/live")`), wrapped as a `Facade::Grouping`,
  anode + detector-volumes attached, and its `children()` taken as the
  `Cluster*` list, **sorted by length descending**.
- **Flashes**: enumerated through the shared facade — `grouping->flashes()`
  (`Clus::Facade::Grouping`) returns one `Facade::Flash` per row of the canonical
  `flash`/`light`/`flashlight` PCs, owning the `flashlight`-join walk. Each is
  wrapped in an `Opflash` matching-adapter (`Opflash(const Facade::Flash&, 0.0,
  m_nchan)` pulls `time()` + the dense `pes(m_nchan)` and the flash `ident()`).
  Dropped unless its time is in `[flash_mintime, flash_maxtime]` **and**
  `total_PE >= flash_minPE`. `flash_id == ident() ==` the canonical flash row
  index, used for the per-cluster `flash` scalar writeback. (`Opflash` adds only
  the matching conventions: `PE_err` 0.3 rule, `fired_channels`, `OpFlashCompare`.)

---

## 4. Algorithm (`operator()`)

### 4.1 Build candidate bundles (`:231-352`)
For every (flash × cluster) pair create a `TimingTPCBundle(flash, cluster,
flash_id, cluster_idx)` (`:269`). For each 3-D point in the cluster:
- drift-correct the X by the flash time: `x += sign * flash_time * m_drift_speed`
  (`:252,289`); `m_drift_speed` is the configurable `drift_speed` (WCT units,
  default `1.563e-3` = the common SBND value), wired from `params.lar.drift_speed`.
- drop points outside the drift / Y / Z bounds; the boundary / truncation flags
  (`close_to_PMT`, `at_x_boundary`, `spec_end`, `window_truncated`) are filled by
  `compute_endpoint_flags` — see §4.1a.
- predict light: convert point to cm, call
  `m_semi_model->detectedDirectVisibilities()` and
  `detectedReflectedVisibilities()` (`:306,308`), accumulate per OpDet
  `pred_flash[idet] += q * QtoL * (dir_vis*VUVeff + ref_vis*VISeff)` (`:310-317`).

Then quality-gate the bundle (`:323-337`): skip if >25% of points drifted out of
bounds; `set_pred_flash`; skip if `get_total_pred_light() < 10`;
`examine_bundle()` (computes KS distance + chi²/ndf); skip if `ks_dis == 1` or
`chi2/ndf > 1e4`. Survivors go into `pre_bundles`.

### 4.1a Endpoint / boundary flags (`compute_endpoint_flags`, `:961-1096`)
Called once per candidate bundle from the build loop (`:472`). Fills four bundle
flags. Port of the prototype `ToyMatching.cxx::calculate_pred_pe` (~176-290).

**Drift coordinate `u`.** Everything below is in a per-TPC drift coordinate
`u = s·(x_corrected − anode_x)`, with `x_corrected = blob_x + flash_x_offset`
(`:403-408`, `:1014`). So `u = 0` at the anode plane and increases to
`u = u_cathode (>0)` at the cathode, regardless of the TPC's physical drift sign
(both reversed-drift SBND APAs share one convention). `flash_x_offset =
±flash_time·drift_speed` makes these flags **T0-dependent** (re-evaluated per
flash hypothesis) — unlike `window_truncated` below, which is T0-independent.

**Endpoints `first_u` / `last_u` (`:1006-1077`).** Collect one `(u, nblobs)` per
non-empty slice, sort ascending in `u`; `first_u` = smallest u = **anode end**,
`last_u` = largest u = **cathode end**. Each end is then *trimmed inward*
(`:1031-1053` anode, `:1055-1077` cathode) to discard sparse, low-charge stubs
poking past the boundary, so the flag tests fire on the real track end. If almost
no blobs lie strictly outside, the end **snaps** exactly to `0` (anode) /
`u_cathode` (cathode).

**Guard / TPC-box containment.** Flags are only set if the trimmed endpoints sit
sensibly inside the TPC: `first_u > anode_in − 1cm`, `0 < last_u < cathode_in`,
`first_u < u_cathode`, where `anode_in = m_anode_ext1 = −2cm` and
`cathode_in = u_cathode + m_cathode_ext1 = u_cathode + 1.2cm`. This 4-part test is
exactly the prototype's `flag_good_bundle` / TPC-containment gate
(`ToyMatching.cxx` 272-275); `compute_endpoint_flags` **returns** it. When
`require_containment` is true (SBND default), the build loop **discards** any
bundle for which it is false — i.e. the cluster, after the flash T0 x-offset, is
not contained in the drift box — before the predicted-light loop, so the bundle
never competes in matching. A cluster with no 3d points under this anode also
returns false (it cannot be contained) and is discarded. Default OFF leaves
non-SBND configs bit-identical (no `continue` is ever taken).

- **`flag_close_to_PMT` — anode end only (`:1087-1090`).** Set true when the
  trimmed anode end `first_u ∈ (anode_in − 1cm, m_anode_ext2] = (−3cm, +4cm]`.
  The SBND PMTs sit at the anode plane (`u≈0`), so this is intentionally
  anode-side only. When it fires it also sets `at_x_boundary`.
- **`flag_at_x_boundary` — anode OR cathode (`:1087-1094`).** True if the anode
  end is in the anode window (above) **or** the cathode end
  `last_u ∈ [u_cathode + m_cathode_ext2, cathode_in) = [u_cathode − 2cm,
  u_cathode + 1.2cm)`. Marks "cluster touches an x (drift) boundary" on either
  side; the cathode case does **not** set `close_to_PMT` (no PMTs there).
- **`flag_spec_end` — set during trimming (`:1046`, `:1070`).** Not a position
  test on the final endpoint: raised when an inward trim walked off **> 10
  slices** of sparse charge (`< 5%` of the cluster's blobs) **and** the last step
  in `u` was a small gap (`< 10cm`) — i.e. a long, contiguous, low-charge tail
  (a genuine track extension, not a detached fragment). Committed via
  `set_spec_end_flag` inside the guard (`:1084`).

**`flag_window_truncated` — raw readout-window edge (`:986-1004`).**
T0-independent and APA-agnostic: a property of the *raw* readout window, computed
from the matched anode's slice indices (raw ticks, `slice_index = start/tick`).
With `min_tick` = smallest `slice_index_min` and `max_tick` = largest
`slice_index_max` over the cluster's non-empty slices:
```
truncated = (min_tick ≤ window_edge_ticks)                       // leading edge
         || (readout_window_ticks − max_tick ≤ window_edge_ticks); // trailing edge
```
SBND config sets `window_edge_ticks = 24` and `readout_window_ticks = 3428`
(`cfg/.../sbnd/qlmatching.jsonnet` §H), giving a **leading band `[0, 24]`** and a
**trailing band `[3404, 3427]`**. `readout_window_ticks` is the *exclusive*
window end: `daq.nticks = 3427`, but with rebin 4 the final 4-tick slice's
`slice_index_max` reaches 3428, so 3428 is the correct trailing reference and the
band starts on the slice boundary 3404 (= 3428 − 24). The flag is currently
**inert** (no consumer reads it); it is only emitted on the `bundle_flags:` debug
line (`:890-899`), so changing the cut does not alter production matching output.

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
Clusters are pre-initialized with `set_scalar<int>("flash", -1)` (alongside the
`set_cluster_t0(-1e12)` init pass). For every bundle in `flash_bundles_map`, set
`main_cluster->set_cluster_t0(flash_time)` **and** `set_scalar<int>("flash",
flash->get_flash_id())` (the canonical flash row index), so
`Clus::Facade::Cluster::get_flash()` returns the matched flash downstream. Then
rebuild the output tensorset from the (now t0- and flash-stamped) live + dead
pctrees and emit it under the same `charge_ident`.

### Branch effects
- `m_data == false` (MC): masks saturated PMTs (`:255-260`).
- `m_beamonly == true`: forces the flash window to the beam window (`:61-64`) and
  **skips** the `strength_cutoff` pruning in both rounds (`:514,598`).

---

## 5. BEE output — how the event-display JSON is dumped

> **Relocated to the all-APA MABC.** The Bee output below (QLMatching writing
> per-APA `data-sep/*-{img,op}-apa<n>.json`, later stitched by `bee-upload.sh` /
> `merge-apa.py`) is **superseded** by a dump inside the all-APA
> `Clus::MultiAlgBlobClustering`, which now writes the `op` (light/flash) layer
> directly into the same `mabc-all-apa.zip` as the charge `img`/`clustering` and
> the dead-area patches — no separate combine step. QLMatching's job is now to
> **persist** the match into the pctree (the `matched_flash_gid` / `flashpred` /
> root `opflash` PCs of the out-port table above); the MABC reads those at its
> pre-clustering `img` dump point (where the live clusters are still the per-APA
> matched clusters and their ids equal the `img` enumeration) and emits the `op`
> JSON via the new `WireCell::Bee::Flashes` object (`save_opflash: true`).
> The legacy `dump_*` helpers described here are retired (see porting-summary).
>
> **Known degeneracy (to fix with later matching-algorithm tuning):** the pctree
> match is **one flash per cluster** (the MicroBooNE model: a single `cluster_t0`
> / `flash` scalar). The legacy `dump_light` emitted a row per *bundle*, so it
> could show the same cluster under two flashes, or a duplicate `(flash,cluster)`
> bundle. Those LASSO degeneracies collapse to the single persisted match in the
> MABC dump (last-wins, identical to the existing `cluster_t0`/`flash` scalar).
> Observed on the 10-event mc sample as ~5 display rows across 4 events; the
> measured-light flashes and the one-flash→many-clusters case are unaffected.

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

- **`Opflash`** (`Opflash.{h,cxx}`) — the matcher's per-flash **adapter over
  `Clus::Facade::Flash`**. Holds per-channel `PE` (+ `PE_err`), `total_PE`, `time`
  (ns), `m_nchan`, `flash_id`, fired-channel list. Two ctors: from a `Facade::Flash`
  (pulls `time()` + `pes(nchan)` + `ident()`) and from `(time, pe vector, threshold,
  nchan)` (the facade ctor delegates to it). `PE_err` is synthesized in-class (the
  0.3 rule) — a matching convention kept out of the generic facade/canonical data.
  No tensor/PC knowledge of its own. Key accessors: `get_PEs()`, `get_PE(ch)`,
  `get_total_PE()`, `get_time()`, `get_num_channels()`, `get_flash_id()`.
- **`Aux::FlashTensorToOpticalPCs`** (`aux/.../FlashTensorToOpticalPCs.{h,cxx}`, plugin
  `WireCellAux`) — `ITensorSetFanin` (2→1) that expands the SBND opflash matrix into
  the canonical `flash`/`light`/`flashlight` PCs on the live root node (§1a). The SBND
  analog of `root/UbooneClusterSource`'s optical loading. **Lives in `aux/`** (not
  `match/`): it is detector-agnostic plumbing, sibling to `Aux::AttachPointCloudToTree`.
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
- **Output carries `cluster_t0` + matched-flash scalar + canonical flash PCs** —
  `get_flash()` works downstream. Other match metadata (strength, beam flag) is
  still BEE-JSON-only; add it to a `tagger_info` PC if downstream needs it.
