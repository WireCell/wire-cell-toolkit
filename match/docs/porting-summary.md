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
| `Opflash.{h,cxx}` | Moved from larwirecell, interface unchanged. |
| `TimingTPCBundle.{h,cxx}` | Moved from larwirecell, interface unchanged. |
| `QLMatching.{h,cxx}` | `ITensorSetFilter` + `IConfigurable` component. Reads a JSON model file at `configure()` via `Persist::load`. Flash arrives via a `flash` point cloud on the cluster tree's live root node (attached upstream by `Aux::AttachPointCloudToTree`), not a 2nd input port. |
| `Util.{h,cxx}` | BEE-JSON dump helpers (`dump_bee_3d`, `dump_bee_bundle`, `dump_light`). |
| `wscript_build` | `bld.smplpkg('WireCellMatch', use='WireCellClus WireCellAux WireCellIface WireCellUtil')` |

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
  -V reality=sim -V DL=6.2 -V DT=9.8 -V lifetime=6 -V driftSpeed=1.565 \
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

## Where things live

- WCT source / branch: `/exp/sbnd/app/users/yuhw/wire-cell-toolkit` (branch `match`)
- Install: `/exp/sbnd/app/users/yuhw/opt` (WCT) and `/opt/larwirecell/...`
- larwirecell (legacy, reference): `/exp/sbnd/app/users/yuhw/larwirecell`
- SBND working area, jsonnets, data builders, archived runs:
  `/exp/sbnd/app/users/yuhw/wcp-porting-img/sbnd/standalone-sample/`
- SBND model JSON: `/exp/sbnd/app/users/yuhw/wire-cell-data/sbnd/photodet/`
