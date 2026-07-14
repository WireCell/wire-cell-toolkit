# `dune10kt-1x2x6` — status of this configuration directory

WireCell config for the **full DUNE Far Detector, Horizontal Drift** — 1×2×6
module (1 cryostat × 2 APA columns × 6 rows = 12 anode planes). HD, not VD.

Last reviewed / updated: **2026-05-18**.

---

## 1. What this is

- DUNE FD-HD, 12 APAs (`params.jsonnet:91`, `std.range(0,11)`).
- Horizontal-drift geometry: `apa_cpa = 3.63075 m`, `apa_w2w = 60.031 mm`,
  3-plane U/V/W.
- Builds on the shared base `pgrapher/common/params.jsonnet`.
- Distinct from `pdhd` (the ProtoDUNE-HD prototype, 4 APAs) and from `dune-vd`
  (vertical drift) — see §6.

## 2. Current status

After the 2026-05-18 changes (§3) **every entry config in this directory
compiles.** Verified with `wcsonnet` (`WIRECELL_PATH=cfg:wire-cell-data`):

| Entry config | Status |
|---|---|
| `wct-sim-check.jsonnet` | ✅ compiles |
| `wct-sim-ideal-sn.jsonnet` | ✅ compiles |
| `wct-sim-ideal-sig.jsonnet` | ✅ compiles |
| `wct-sim-ideal-sn-nf-sp.jsonnet` | ✅ compiles |
| `wcls-nf-sp.jsonnet` (LArSoft NF+SP) | ✅ compiles |
| `wcls-sp.jsonnet` (LArSoft SP) | ✅ compiles |
| `wcls-sim-drift-simchannel.jsonnet` | ✅ compiles |
| `wcls-sim-drift-simchannel-nf-sp.jsonnet` | ✅ compiles |
| `wcls-sim-drift-simchannel-splusn.jsonnet` | ✅ compiles |
| `wcls-blip-sim-drift-simchannel.jsonnet` | ✅ compiles |

"Compiles" means the Jsonnet evaluates. It does **not** prove a job runs
end-to-end, nor that every referenced data file is loadable — the data files
were checked separately (§3) and are present, but the configs are
compile-verified, not run-verified.

## 3. What was changed (2026-05-18)

Two rounds of fixes were applied; together they leave the directory fully
compiling.

### 3a — `params.jsonnet`: inputs aligned to PDHD

The directory did not compile at all before this date: commit `41e02736`
("cfg/pdvd: sync NF+SP configs from dunereco") removed the `elec.gain` and
`adc.resolution` defaults from the shared base `pgrapher/common/params.jsonnet`.
`dune10kt-1x2x6` had relied on those inherited defaults and never set its own,
so every entry config failed with `RUNTIME ERROR: Field does not exist: gain`
(the `uboone` config was hit the same way and was repaired in `dcaa9f99`; this
directory was missed because nothing tests or runs it).

`params.jsonnet` was updated — **field response, noise spectra, FE gain and
shaping now follow the PDHD config**:

| Field | Before | After |
|---|---|---|
| `elec.gain` | *(unset — broke compile)* | `14.0 mV/fC` — PDHD default |
| `elec.shaping` | `2.2 µs` | `2.2 µs` — unchanged, already matches PDHD |
| `adc.resolution` | *(unset — broke compile)* | `14` bits — same as PDHD / DUNE FD-HD |
| `files.fields` | `dune-garfield-1d60563.json.bz2` *(absent from wire-cell-data)* | `dune-garfield-1d565.json.bz2` — the generic DUNE HD response PDHD uses |
| `files.noise` | `protodune-noise-spectra-v1.json.bz2` *(MicroBooNE-derived, self-flagged "bogus")* | `protodunehd-noise-spectra-14mVfC-v1.json.bz2` — PDHD, matched to the 14 mV/fC gain |

`postgain` (`1.1365`) was left untouched — its calibration comment already
assumes 14 mV/fC, consistent with the new `elec.gain`. Both new data files are
present in `wire-cell-data`.

### 3b — entry configs: removed dependence on deleted helpers

The 7 simulation entry configs imported helper files that no longer exist —
`pgrapher/common/fileio.jsonnet` and `pgrapher/ui/cli/nodes.jsonnet` (both
removed repo-wide; see §4) — and three also referenced an undefined
`sim.frame_sink`. They were fixed following the **PDHD pattern**: the working
`pdhd_sim` configs use no such helpers; they inline a `FrameFileSink` and the
`wire-cell` command-line node directly.

- **4 `wcls-sim-drift-*` configs** — their `io` (fileio) use was already fully
  commented out, so the `import` was dead code. The one dead line was removed
  from each.
- **3 `wct-sim-ideal-*` configs** — rewritten to be self-contained like the
  `pdhd_sim` configs: dropped the `fileio.jsonnet` / `cli/nodes.jsonnet`
  imports; the npz frame-saver and the undefined `sim.frame_sink` were replaced
  by an inline `FrameFileSink` (**output is now `.tar.bz2`, not `.npz`**); the
  command-line node is inlined; the depo-file save was dropped (no pdhd/pdvd
  sim config saves depos — re-add an inline depo sink if depo output is ever
  needed). `wct-sim-ideal-sn.jsonnet` additionally now imports
  `simparams.jsonnet` — it had wrongly imported the base `params.jsonnet`,
  which lacks the sim-only `sys_resp`/`sys_status` fields — consistent with its
  two sibling sim configs.

## 4. Known limitations / remaining issues

- **`pgrapher/common/fileio.jsonnet` and `pgrapher/ui/cli/nodes.jsonnet` are
  missing repo-wide.** Commit `6b8ef2e2` ("Remove helpers in favor of layers")
  removed them, yet ~75 configs across nearly every experiment (`pdsp`,
  `uboone`, `protodunevd`, `pdhd`, `sbnd`, `icarus`, `dune-vd`, …) still import
  them. **This directory no longer depends on either** (see §3b), but the
  repo-wide breakage remains and is out of scope here. Note
  `cfg/pgrapher/common/io.jsonnet` is not a usable replacement — it is itself
  incomplete (a syntax error at `io.jsonnet:34`, missing its closing brace).
- **Noise spectra are PDHD's, not FD-tuned.** `protodunehd-noise-spectra-14mVfC-v1.json.bz2`
  is a ProtoDUNE-HD measurement used here as the best available HD stand-in. A
  real DUNE FD-HD electronics-noise model would eventually be preferable.
- **Field response is the generic DUNE HD Garfield 1D file.** PDHD's
  APA0-specific best-fit response (`np04hd-garfield-6paths-mcmc-bestfit.json.bz2`)
  is a ProtoDUNE-HD detector calibration and is deliberately *not* used here.
- **NF/SP algorithm files are stale** (`nf.jsonnet` ~2021, `sp-filters.jsonnet`
  ~2020, `sp.jsonnet` ~2023). They have not been revisited; the 2026-05-18 work
  only addressed the four input settings above. The NF still uses MicroBooNE
  noise components (`mbOneChannelNoise`, flagged in `nf.jsonnet` itself).

## 5. What uses this config

Searched the whole repo. **No test, CI, script or build target runs it.** Its
only live consumer is `dune-vd`, which imports one file —
`chndb-perfect.jsonnet`:

- `cfg/pgrapher/experiment/dune-vd/wcls-sim-drift-simchannel.jsonnet:118`
- `cfg/pgrapher/experiment/dune-vd/wcls-sim-drift-simchannel-nf-sp.jsonnet:158`
- `cfg/pgrapher/experiment/dune-vd/wcls-sim-drift-simchannel-nf-sp-img.jsonnet:158`

Note: `dune-vd/sp.jsonnet:6` carries a stale "BIG FAT FIXME: we are taking from
dune10kt-1x2x6" comment — outdated, that file now imports its own
`dune-vd/sp-filters.jsonnet`. Worth correcting on the `dune-vd` side.

## 6. Relationship to `pdhd`

`pdhd` and `dune10kt-1x2x6` are **independent lineages**, not a fork pair, and
have never been co-maintained (no commit touches both). `pdhd` was created
2025-04 and synced from an external dunereco fork; `dune10kt-1x2x6` dates to
2020. They are also **different detectors** — full FD vs. the 4-APA prototype.

The 2026-05-18 update borrows PDHD's *input settings* (gain, shaping, field,
noise) because PDHD is the actively-tuned HD reference. It does **not** copy
PDHD's prototype-specific machinery (FEMB neg-pulse groups, per-APA APA0/APA1
tuning, hand-listed bad channels, L1SP / DNN-ROI) — none of that applies to the
full FD. A wholesale "sync to pdhd" is neither needed nor correct.

## 7. Recommendation

The directory is now self-consistent and **every entry config compiles**; it is
suitable as an HD reference config. Before it could be used for real DUNE FD-HD
processing, FD-specific electronics-noise spectra would eventually replace the
borrowed PDHD file (§4).

Until there is an active FD-HD WCT workflow, no further investment (NF/SP
modernization, per-APA tuning) is warranted — keep it as a maintained reference.
