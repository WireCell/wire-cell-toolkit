# Joint two-TPC (both-APA) charge-light matching — design & status

**Status:** the *infrastructure* for joint both-APA matching is built and merged
(Steps 1–2 below). The actual joint *algorithm* (combine ±80 ns flashes, drop the
per-TPC OpDet mask, one global LASSO) is still deferred.

**Decision (supersedes the earlier "fork" plan):** **convert `QLMatching` in
place**, do not fork a new component. One set of files to maintain. This was an
explicit, authorized reversal of the earlier "fork a new `JointQLMatching` by
duplication" decision (and of the general [[feedback_fork_by_duplication]]
convention) — see that memory's recorded exception. `QLMatching` keeps the
single-APA path exactly (default config), so pdhd/pdvd/uboone and the SBND
per-APA dev chain are byte-for-byte unaffected; SBND opts into the joint path via
config.

---

## 1. Motivation

SBND Q/L matching historically ran **once per APA**: two independent `QLMatching`
nodes, each seeing only its own TPC's clusters and that TPC's optical flashes
(`opflash_apa{0,1}.tar.gz`), joined only downstream by the all-APA
`PointTreeMerging`. Running matching jointly over both APAs will eventually let us
(a) combine the two per-TPC flashes that are one physical scintillation event into
one joint flash (full-PMT PE), and (b) use both TPCs' cluster geometry in one
solve so a cross-cathode interaction's light is not artificially split.

---

## 2. What is built (merged)

### Step 1 — factor the pipeline into per-APA stage functions (bit-identical)
`QLMatching::operator()` was a ~700-line inline loop. It is now split into named
stage functions, each operating on a per-run `ApaRun` state struct that carries
the bundle maps, index maps, flash vector, geometry and charge bookkeeping:

    build_opdet_mask -> read_flashes -> decompose_cluster_groups ->
    compute_geometry -> build_bundles [Stage 1] -> build_bundle_maps ->
    cull_inconsistent [Stage 1] -> fit_round1 [Stage 2] -> fit_round2 [Stage 3,
    which calls organize_bundles for Stage 4 merge + Stage 5 out-of-beam QA] ->
    apply_matched_t0s -> write_opflash_pc

`run_one_apa()` runs the sequence. **Invariant:** every per-run container lives in
`ApaRun`; nothing per-run is a `QLMatching` member or static. That isolation is
what keeps the pointer-keyed map iteration deterministic when one node processes
both APAs (see [[project_ql_matching_pointer_nondeterminism]]).

### Step 2 — `QLMatching` becomes a multi-APA fanin (default off)
`QLMatching` is now an `ITensorSetFanin` (was `ITensorSetFilter`) with a
configurable multiplicity = number of configured anodes.

- **Framework constraint that drove the shape:** a Wire-Cell node returns a single
  `INode::category()`, so it cannot be both a filter and a fanin. A fanin with
  **multiplicity 1** declares exactly one input port and wires/behaves identically
  to the old single-input filter — so the single-anode default is bit-identical.
- Single `anode` (historical config) => multiplicity 1, one input port, exact old
  behavior. A list `anodes` => one input port per APA.
- Each input is matched in its own fresh isolated `ApaRun` (per-APA flash set,
  per-TPC OpDet mask, per-TPC geometry, per-APA `flash_gid` stride), reusing
  `run_one_apa()` unchanged. `compute_endpoint_flags` now takes the per-APA anode
  ident instead of reading `m_anode`.
- With >1 input the per-APA trees are **merged in-node**, reproducing the
  standalone `clus_all_apa` `PointTreeMerging` it replaces: primary = input 0,
  concatenate `root_pcs_to_merge=['opflash']`, `take_children`, normalize local-PC
  keys. The disjoint per-APA `flash_gid` stride (`anode_ident*1e6 + idx`) keeps
  merged ids collision-free, so the union needs no remap.

### Config / drivers
- `cfg/.../sbnd/qlmatching.jsonnet`: shared `data` block factored into
  `match_data()` (single source of truth); `matching()` (per-APA, unchanged
  output) and new `matching_joint(anodes, dv, ...)` both use it. The joint builder
  declares the all-anode DetectorVolumes as a `uses` dependency.
- `cfg/.../sbnd/clus.jsonnet`: `all_apa(..., premerged=false)`; when true it skips
  the `PointTreeMerging` fanin and feeds the single pre-merged tree to MABC.
- `sbnd_xin/wct-clus-matching-perevt.jsonnet` (in the wcp-porting-img repo): a
  `joint=false` TLA toggle. Off = historical per-APA path; on = one
  `matching_joint` node fed by both `flash_attach` outputs, then MABC directly.

### Verification (the bit-identical contract)
On SBND mc evt 2 and evt 11 (both APAs), `work/ql_evt<id>/mabc-all-apa.zip`:
- Step 1 and Step 2-OFF: byte-identical to the pre-change binary.
- Step 2-ON (joint toggle): byte-identical to the OFF per-APA + `PointTreeMerging`
  path. This proves the in-node merge reproduces `PointTreeMerging`, the per-APA
  matching matches the standalone matchers, the all-anode DV yields identical
  per-TPC geometry, and no cross-APA state leaks.

---

## 3. Verified structural facts (grounding for the deferred algorithm)

1. `Blob::wpid()` / `WirePlaneId::apa()` expose each cluster's TPC, so the joint
   solve can read a cluster's TPC directly (no x-sign heuristic).
2. The ±80 ns flash combination is essentially new work. Today's
   `MultiAlgBlobClustering::store_flash_groups` (`flash_group_window: 80*wc.ns`)
   only computes a per-flash group-ID for Bee display, after matching; it does not
   produce combined-PE flash vectors. The greedy rule is reusable as the
   coincidence test; the PE summation across the full PMT set is new.
3. The `flashlight` PC is an index-based join (`FlashTensorToOpticalPCs`), so the
   per-APA optical PCs must NOT be blindly concatenated — flashes must be combined
   in raw opflash-matrix space, then one `FlashTensorToOpticalPCs` builds a
   self-consistent set. (Why the joint *algorithm* combines flashes upstream of the
   matcher, rather than inside it.)

---

## 4. Deferred — the actual joint algorithm

Now that `QLMatching` ingests both APAs in one node and the wiring is proven
bit-identical, the remaining work turns the per-APA loop into a true joint solve:

- **Combine the flash:** a pre-step (matrix space) sums ±80 ns-coincident TPC0/TPC1
  flashes over the full PMT set into one joint flash; one-sided flashes pass
  through. (Open: confirm each per-APA opflash row spans all `nchan` PMTs with the
  other TPC's channels ≈ 0.)
- **Drop the per-TPC OpDet mask** so predicted light is summed over all PMTs and
  compared to the combined flash's full-PMT PE. (Confirm the semi-analytical model
  returns visibilities for all PMTs from a point in either TPC — it takes a global
  xyz, so expected yes. Note the hard cathode optical cut in
  `SemiAnalyticalModel`: opposite-TPC OpDets currently get exactly 0 predicted PE.)
- **One global LASSO:** all (cluster, joint-flash) bundles — clusters from both
  TPCs — in one solve, per-cluster TPC geometry selected by `wpid().apa()`.
  Re-check the chi²/KS gates and the deterministic iteration order over the now
  mixed-TPC bundle set; re-tune thresholds explicitly if needed, don't silently.
- **Cross-cathode cluster joining** (further out): pair a track's two cross-cathode
  halves and predict combined light from the full track. Partly falls out for free
  if both halves match the same joint flash (same `cluster_t0`, aligned at the
  cathode under the downstream `x_t0cor`).

---

## 5. Files

- `match/src/QLMatching.cxx`, `match/inc/WireCellMatch/QLMatching.h`
- `cfg/pgrapher/experiment/sbnd/qlmatching.jsonnet`,
  `cfg/pgrapher/experiment/sbnd/clus.jsonnet`
- `sbnd_xin/wct-clus-matching-perevt.jsonnet` (wcp-porting-img repo)
- Reference: `clus/src/PointTreeMerging.cxx` (the merge reproduced in-node),
  `aux/src/FlashTensorToOpticalPCs.cxx`, `match/src/SemiAnalyticalModel.cxx`
  (cathode optical cut), `iface/inc/WireCellIface/{INode,IFaninNode}.h`.
