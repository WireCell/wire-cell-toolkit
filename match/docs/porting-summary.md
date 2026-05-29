# Porting summary: larwirecell `qlmatch` → standalone `WireCellMatch`

Orientation doc for future work. For the deep dive (larsoft→WCT API mapping,
JSON schema, gotchas) see [`qlmatching-port.md`](qlmatching-port.md).

## Goal

Port SBND charge–light (QL) matching out of `larwirecell/larwirecell/qlmatch/`
(which depends on larsoft/`art`/`fhicl`) into a self-contained
wire-cell-toolkit subpackage that runs under plain `wire-cell`, with **no
larsoft dependency at runtime**. The larsim `phot::SemiAnalyticalModel` and its
`SBNDOpticalPath` tool had to be ported too.

## What was built

New subpackage **`wire-cell-toolkit/match/`** (`WireCellMatch`):

| File | Role |
|------|------|
| `SemiAnalyticalModel.{h,cxx}` | Port of larsim's `phot::SemiAnalyticalModel` (SBND scope: dome PMTs + flat (X)Arapucas at anode/cathode orientation; VUV direct + VIS reflected). No larsoft deps. |
| `Opflash.{h,cxx}` | The matcher's per-flash adapter over `Clus::Facade::Flash`. Two ctors: from a `Facade::Flash` (pulls `time()`+`pes(nchan)`) and from `(time, pe, threshold, nchan)`. Adds the matching-only `PE_err` (0.3 rule), `fired_channels`, and `OpFlashCompare`. No tensor/PC knowledge. |
| `TimingTPCBundle.{h,cxx}` | Moved from larwirecell, interface unchanged. No canonical-schema equivalent — algorithm working object only. |
| `QLMatching.{h,cxx}` | `ITensorSetFilter` + `IConfigurable` component. Reads a JSON model file at `configure()` via `Persist::load`. Flash arrives via the canonical `flash`/`light`/`flashlight` PCs on the live root node (written by `Aux::FlashTensorToOpticalPCs`), not a 2nd input port; writes back a per-cluster matched-flash scalar. |
| `Util.{h,cxx}` | BEE-JSON dump helpers (`dump_bee_3d`, `dump_bee_bundle`, `dump_light`). |
| `wscript_build` | `bld.smplpkg('WireCellMatch', use='WireCellClus WireCellAux WireCellIface WireCellUtil')` |

The opflash-matrix → canonical-optical-PC converter started life as
`Match::OpflashToFlashPCs`, was relocated to `aux/` as `Aux::OpflashToFlashPCs` (no detector
physics — it only reshapes a `[nflash, 1+nchan]` matrix into the canonical PCs, so it belongs
with the generic tensor/pctree plumbing next to `Aux::AttachPointCloudToTree`), and was then
**renamed to `Aux::FlashTensorToOpticalPCs`** (`aux/{inc/WireCellAux,src}/FlashTensorToOpticalPCs.{h,cxx}`,
plugin `WireCellAux`) to drop the confusing `Opflash→Flash` echo — its name now reflects
input (a flash *tensor*) → output (the canonical *optical* PCs). The rename changed the WCT
factory type string to `"FlashTensorToOpticalPCs"`, so the one config that names it
(`cfg/pgrapher/experiment/sbnd/qlmatching.jsonnet`) was updated in lock-step; matching output
is unchanged (mc 40/40 byte-identical).

The standalone matching graph nodes (opflash reader → converter → `QLMatching`) plus the
SBND matching constants (`nchan`, `ch_mask`) are built by a canonical config helper
**`cfg/pgrapher/experiment/sbnd/qlmatching.jsonnet`** (a factory `function(params)` mirroring
`clus.jsonnet`/`img.jsonnet`), re-exported by the standalone via a thin
`sbnd_xin/qlmatching.jsonnet` shim. `drift_speed` flows from the caller's `params` (the common
SBND `1.563 mm/us`).

Key larsoft→WCT substitutions: `geo::Point_t`→`WireCell::Point`;
`fhicl::ParameterSet`→`WireCell::Configuration` (Json) loaded via
`Persist::load`; `cet::exception`→`raise<ValueError>`;
`mf::LogInfo`→`Aux::Logger`; `SBNDOpticalPath` tool → inlined same-TPC (X-sign)
check; geometry/opdet arrays → plain data loaded from a JSON file.

## SBND data file

`QLMatching` reads one JSON (default `sbnd/photodet/semi-analytical-sbnd.json`,
resolved via `WIRECELL_PATH`) with keys `VUVHits`, `VISHits`, `Geometry`,
`OpDets`. It currently lives at
`/exp/sbnd/app/users/yuhw/wire-cell-data/sbnd/photodet/semi-analytical-sbnd.json`.

Generated (one-off) by tools in
`wcp-porting-img/sbnd/standalone-sample/build-semi-analytical-data/`:
- `SBNDOpDetDumper_module.cc` + `dump_sbnd_opdets.fcl` — larsoft analyzer that
  prints the 312 SBND OpDets (`OPDET:` CSV) plus active-volume / cathode /
  `vuv_absorption_length` constants (`GEOM:` lines). Built in larwirecell's
  `qlmatch/CMakeLists.txt`.
- `build_semi_analytical_sbnd_json.py` — parses the `fhicl-dump` of
  `semimodel_sbnd.fcl` + the OpDet CSV into the final JSON.

**Critical constant:** `vuv_absorption_length = 2000 cm` — take it from
`LArPropertiesService::AbsLengthSpectrum()` @ 9.7 eV, *not* a guess. An early
85 cm placeholder silently halved the predicted PE.

## Build / run / verify

Everything builds in the SL7 apptainer with the sbndcode env, installed to
`/exp/sbnd/app/users/yuhw/opt`. The wrapper is
`/exp/sbnd/app/users/yuhw/claude-utilities/in-gpvm-sl7.sh`.

```bash
# build + install
in-gpvm-sl7.sh bash -c '
  source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
  setup sbndcode v10_14_02_03 -q e26:prof
  cd /exp/sbnd/app/users/yuhw/wire-cell-toolkit
  ./wcb -p --notests install'        # configured with --prefix=/exp/sbnd/app/users/yuhw/opt

# run the standalone (reads pre-dumped icluster-*.npz + opflash_*.tar.gz)
source /exp/sbnd/app/users/yuhw/wcp-porting-img/sbnd/setup-local-opt.sh
cd <dir with icluster-apa{0,1}-{active,masked}.npz and opflash_apa{0,1}.tar.gz>
wire-cell -l stdout -L info \
  -V reality=sim -V DL=6.2 -V DT=9.8 -V lifetime=6 \
  -V input=. -V semimodel_file=sbnd/photodet/semi-analytical-sbnd.json \
  -c wct-clus-matching-standalone.jsonnet
```

`setup-local-opt.sh` points the env at `/opt` for both WCT and (if present)
larwirecell, and puts `wire-cell-data` on `WIRECELL_PATH`.

Inputs for SBND are produced once from an artROOT file by
`wcls-img-dump.fcl` (→ `icluster-*.npz`) and `wcls-flash-dump.fcl`
(→ `opflash_*.tar.gz`).

### Three runnable paths (same physics, used for cross-checks)

1. **Standalone** `wire-cell -c wct-clus-matching-standalone.jsonnet` —
   loads `WireCellMatch`; reads npz. The deliverable.
2. **Legacy `lar` standalone** `lar -c standalone-sample/wct-clus-matching.fcl` —
   larwirecell `WireCellQLMatch` (larsim model); also reads npz via WCT's
   `ClusterFileSource`.
3. **Full in-job WCLS** `lar -c wcls-img-clus-matching.fcl -s <art.root>` —
   sig→img→clustering→matching in one job, **no npz round-trip**; the
   reference for the deadarea check.

## Fixes discovered while validating (all on branch `match`)

| Commit | Fix |
|--------|-----|
| `4b196e5e` | New `match` subpackage (the port itself). |
| `ccfee6ff` | Doc: `vuv_absorption_length=2000` + numerical-agreement note. |
| `6eb0de08` | `QLMatching`: configurable `strength_cutoff` (was hard-coded 0.05) **and** deterministic LASSO column/row order (sort by `flash_id` / global cluster index instead of `std::map<Pointer*>` iteration, which was allocator-dependent and flipped marginal bundles run-to-run). |
| `cdecb9e6` | Doc: legacy `WireCellQLMatch` factory-name coexistence + `/opt` install layout. |
| `b2352c39` | `sio/ClusterFileSource`: don't re-`load_filename()` over a header that `load_numpy()` already peeked at the ident boundary — fixed EOS-after-one-event when a single npz holds many events (multi-event dumps). |
| `2fa2e5b3` | `aux/ClusterArrays::to_cluster`: pass `nudge=1e-3` (matching `Img::GridTiling`) to `RayGrid::Blob::add` so deserialized blob corners match imaging. Fixes dead-region polygons collapsing (quad→triangle) and a dropped boundary strip at the z≈501 cm edge. |
| `58dad76b` | `match/SemiAnalyticalModel`: finite-safe `angle_bin()` clamp on the VUV `j`, VIS `k`, and dome-model `j` indices — guards the unchecked `[bin]` table reads against `theta==90°` / NaN points (latent segfault). Defensive only; results unchanged. |
| `f8b91803` | `match/Util`: new `dump_light()` — dump **every** flash to `*-op.json`, not just matched ones. Matched flashes keep `cluster_id`/`op_pes_pred` (same filter as `dump_bee_bundle`); unmatched flashes are emitted with an empty `cluster_id` so they still show in the BEE light display. `QLMatching` now calls this instead of `dump_bee_bundle`. Mirrored in larwirecell `qlmatch` (commit `d634638` on `dev-v10_14_02_02`). |
| `985ec8fb` | **Flash representation consolidated onto the toolkit's canonical schema.** SBND flashes now ride the cluster pctree root node as the same `flash`/`light`/`flashlight` PCs the MicroBooNE `root/UbooneClusterSource` writes (not a bespoke `[nflash,1+nchan]` matrix). New `match/FlashTensorToOpticalPCs` expands the SBND opflash matrix into those PCs; `QLMatching` rebuilds its `Opflash` objects from them (new `Opflash(time,pe,thr,nchan)` ctor; `nchan` config, default 312) and writes back a per-cluster `flash` scalar so `Clus::Facade::Cluster::get_flash()` reflects the match. Deliberate divergences from Uboone `load_optical`: flash/light time stored **raw** (no `units::us`, to keep matching identical); per-channel `error`=0 (the `PE_err` 0.3 rule stays in `Opflash`). `TimingTPCBundle` has no canonical equivalent — unchanged. Output bit-identical: mc 40/40 byte-identical to baseline; data runs all 10 events. |
| `3ba2698a` | `match/QLMatching`: drift speed for the per-flash X correction is now the configurable `drift_speed` (WCT units), defaulting to the historical `1.563e-3`; the standalone jsonnet wires `params.lar.drift_speed` (the common SBND `1.563 mm/us`) instead of the previous hard-coded literal and a separate `driftSpeed` run-script override. Output unchanged (40/40 byte-identical) since the common value equals the old hard-coded one. |
| `5ef34045` | `match/TimingTPCBundle`: null-safe ctor. An unmatched cluster is a bundle with a null flash (`QLMatching:629`); the ctor used to deref `flash->get_num_channels()` and segfault. MC never hit it (all clusters match); `data` mode does (cosmics with no flash). Now `m_nchan = flash ? flash->get_num_channels() : 0;`. MC output unchanged (40/40 byte-identical); `data` runs all 10 events. |
| _relocate_ | **`OpflashToFlashPCs` relocated `match/` → `aux/`** as `Aux::OpflashToFlashPCs` (plugin `WireCellAux`). It has no detector physics — only reshapes the `[nflash,1+nchan]` matrix into the canonical PCs — so it sits with the generic tensor/pctree plumbing next to `Aux::AttachPointCloudToTree`. Factory type string unchanged at this step ⇒ configs untouched. Verified: factory registers once (in `libWireCellAux`, 0 refs left in `libWireCellMatch`); resolved config byte-identical; mc 40/40 byte-identical. |
| _config_ | **SBND matching config hoisted** into `cfg/pgrapher/experiment/sbnd/qlmatching.jsonnet` (factory `function(params)` → `opflash_source`/`flash_attach`/`matching`, mirroring `clus.jsonnet`), holding the matching-only constants `nchan`/`ch_mask` and sourcing `drift_speed` from `params.lar.drift_speed`. The standalone `wct-clus-matching-standalone.jsonnet` now builds its nodes via this helper (re-exported by a thin `sbnd_xin/qlmatching.jsonnet` shim), dropping ~60 inline lines. Pure jsonnet refactor: resolved config byte-identical (`wcsonnet` diff), mc 40/40 byte-identical. |
| _rename_ | **`OpflashToFlashPCs` → `FlashTensorToOpticalPCs`** (class + `WIRECELL_FACTORY` type string + the one config that names it). Drops the confusing `Opflash→Flash` echo; the name now reads input (flash *tensor*) → output (canonical *optical* PCs). Resolved-config type string changes (by design); mc 40/40 byte-identical; factory registers once under the new name. |
| _facade_ | **Canonical flash facade expanded + matcher consolidated onto it.** Promoted the nested `Cluster::Flash` to a standalone `Clus::Facade::Flash` (additive: `idents()`/`pes(nchan)` accessors); added `Grouping::flashes()` enumerating all flashes via the single shared flashlight-join walk (also used by `Cluster::get_flash()`). `QLMatching` now reads flashes via `grouping->flashes()` and `Opflash` became a thin adapter over `Facade::Flash` (new facade ctor; dropped the dead tensor ctor + its `ITensor`/boost deps). uboone & SBND now share one flash representation. Backward-compatible (existing `get_flash()` consumers untouched); mc 40/40 byte-identical, data all 10 events. |

## Verification status (10 SBND events, ids 2,9,11,12,14,18,31,35,41,42)

- **Matching**: standalone vs legacy-`lar` op.json — **173/173 bundles
  identical**, 0 only-legacy / 0 only-standalone, max rel. pred-PE diff
  ~3e-4 (float noise). The determinism fix removed the previous marginal-bundle
  flapping.
- **Deadarea**: after the `to_cluster` nudge, standalone `channel-deadarea`
  matches the full in-job WCLS pipeline within **0.01 cm** (identical polygon
  counts + topology on both APAs; was 10-vs-9 polygons + triangle/quad
  mismatches before).
- All three paths run all 10 events to completion and were uploaded to BEE for
  visual inspection.
- **`dump_light` cross-check** (10 events, both APAs, all three paths): `op_t`,
  `op_pes`, `op_peTotal`, and **`cluster_id` are identical everywhere** — the
  matching outcome (incl. the new empty-`cluster_id` unmatched flashes) agrees
  bit-for-bit. The legacy `lar` and full-pipeline `op_pes_pred` are byte-equal
  (both use larsim's model); the standalone differs only by ≤1.3% on tiny
  per-channel predicted PE (~1–8 PE) — the expected ported-model float drift,
  changing no match.

## Known caveats / not done

- **larwirecell is unchanged** and still uses larsim's `phot::SemiAnalyticalModel`
  (cvmfs, not vendored). The GH-index clamp and the SemiAnalyticalModel port
  live only WCT-side. The legacy `lar` paths still depend on larsoft.
- **Factory-name overlap**: both `WireCell::Match::QLMatching` and the legacy
  `WireCell::QLMatch::QLMatching` register a factory named `"QLMatching"`.
  They never coexist in one process (different plugins/jsonnets), but don't
  load both at once. Rename one `WIRECELL_FACTORY` first arg if you ever must.
- **SemiAnalyticalModel scope is SBND-minimal**: lateral PDs, anode
  reflections, Xe absorption, field-cage transparency, vertical-border mode and
  disk PMTs are **not** ported. Add the corresponding larsim branches if a
  future detector needs them.
- The `match` branch is local only (nothing pushed to any remote).

## Making matching generic for other detectors (consolidation roadmap)

`QLMatching` today runs SBND only, but other detectors (ProtoDUNE-HD/VD) will want the
same matcher with a different light system + photon library. The **idiomatic toolkit
pattern** (confirmed by `FieldResponse`/`Drifter`) is **one configurable C++ class fed
detector-specific data files + a per-detector jsonnet helper**, *not* detector subclasses.
This round created the SBND helper (`cfg/.../sbnd/qlmatching.jsonnet`); a sibling per
detector + a per-detector `semi-analytical-<det>.json` is the extension point.

What is still SBND-locked in the C++, and the path to consolidate — **keeping results
byte-identical first** (defaults equal to today's literals, re-gate, *then* switch the
source):

1. **Fiducial / geometry bounds — top priority; reuse the existing clustering tools.**
   `QLMatching` hardcodes the active-volume cuts (X ±2000 mm, Y ±2000 mm, Z [0,5000] mm).
   It **already holds an `IDetectorVolumes` (`m_dv`)** but uses it only for name/volume
   lookup. The official detector-volume machinery clustering already uses is the
   `DetectorVolumes` pnode + the `dvm` fiducial-volume structure in `clus.jsonnet`
   (`FV_xmin/xmax/ymin/ymax/zmin/zmax` with margins) + `PCTransformSet`. The standalone
   **already passes** `clus_maker.detector_volumes([anode])` into `QLMatching`, so the
   wiring exists — the consolidation is to **source the matching FV cuts from `m_dv`**
   instead of the literals. ⚠ These feed the matching cuts, so this is its own
   incremental, separately bit-identical-gated change — do not bundle with config-only work.
2. **OpDet → TPC assignment.** The even/odd-OpDet → TPC-0/1 (drift-X-sign) rule is hardcoded;
   derive it from the OpDet X positions already present in the `semimodel_file` JSON.
3. **Efficiency-array defaults.** The 312-entry VUV/VIS efficiency defaults are baked into the
   header (already config-overridable); move the SBND defaults out into the per-detector
   semimodel JSON.
4. **Per-detector config.** Add a `qlmatching.jsonnet` per detector (SBND done) + a
   per-detector `semi-analytical-<det>.json` (Geometry + OpDets + VUVHits/VISHits).
5. **`SemiAnalyticalModel` porting-completeness gap (separate from architecture).** Lateral
   PDs, anode reflections, Xe absorption, field-cage transparency, and disk PMTs are **not**
   ported (see *Known caveats* above) — a physics-coverage TODO, not a config decision.

**Where the converter belongs (resolved):** `FlashTensorToOpticalPCs` was relocated to `aux/` this
round because it is detector-agnostic. A future polish (not now) is a config-driven column
mapping (time column / PE columns) so non-SBND opflash dumps with a different matrix layout can
reuse it without a code change.

**`Opflash` is now a thin adapter over `Clus::Facade::Flash` (consolidation done).** `Opflash`
stays in `match/` because it is **algorithm-internal state** (used only by `QLMatching` +
`TimingTPCBundle`) — but it no longer duplicates the flash-reading. The canonical flash facade
was promoted to a standalone `Clus::Facade::Flash` and given a Grouping-level enumerator
`Grouping::flashes()` (one shared flashlight-join walk, also used by `Cluster::get_flash()`).
`QLMatching` now reads its flashes via `grouping->flashes()` and wraps each `Facade::Flash` in an
`Opflash` (new `Opflash(const Facade::Flash&, threshold, nchan)` ctor pulling `time()` +
`pes(nchan)`); `Opflash` keeps only the matching-specific conventions — the `PE_err` 0.3 rule,
`fired_channels`, and the `OpFlashCompare` LASSO ordering — and dropped its old tensor-parsing ctor
(now detector/tensor-agnostic). uboone and SBND share the one `Facade::Flash` representation; the
matcher adds its conventions on top. Output unchanged: mc 40/40 byte-identical; data all 10 events.
`TimingTPCBundle` remains transient and unchanged (see below).

**`TimingTPCBundle` stays transient — matches the uboone pctree pattern (record only).** There
is **no** Facade/clus equivalent of a flash↔cluster *bundle*, and there is no precedent for
persisting one: `root/UbooneClusterSource` writes the optical PCs + a per-cluster `flash` index
scalar (the *association*) and nothing else — no bundle, predicted light, chi²/KS, or strength.
SBND already matches this exactly (`QLMatching` writes `cluster_t0` + the per-cluster `flash`
scalar; the bundle objects, `pred_flash`, and match-quality stay in-memory and go only to the BEE
JSON). So `TimingTPCBundle` correctly stays a transient algorithm working object. Persisting match
quality (`pred_flash`/`chi2`/`strength`) into the pctree would be a **new feature beyond uboone
parity, with no current consumer** — deliberately out of scope. If a downstream consumer ever
needs it, the idiomatic slot is per-cluster scalars/arrays (the `set_scalar`/`tagger_info`
pattern), not a new persistent bundle type.

## Where things live

- WCT source / branch: `/exp/sbnd/app/users/yuhw/wire-cell-toolkit` (branch `match`)
- Install: `/exp/sbnd/app/users/yuhw/opt` (WCT) and `/opt/larwirecell/...`
- larwirecell (legacy, reference): `/exp/sbnd/app/users/yuhw/larwirecell`
- SBND working area, jsonnets, data builders, archived runs:
  `/exp/sbnd/app/users/yuhw/wcp-porting-img/sbnd/standalone-sample/`
- SBND model JSON: `/exp/sbnd/app/users/yuhw/wire-cell-data/sbnd/photodet/`
