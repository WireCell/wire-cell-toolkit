# `dune-vd-coldbox` — response-config sync status

WireCell config for the **DUNE Vertical-Drift CRP1 cold-box** test stand — a
single Charge Readout Plane (`dunevdcb1-3view-wires-v2-splitanode`). Builds on
the shared base `pgrapher/common/params.jsonnet`.

Last reviewed / updated: **2026-05-18**.

---

## 1. What this is

A single-CRP test config. The `elec` block branches on the `active_cru`
external variable: `'tde'` selects **top-drift** electronics, any other value
selects **bottom-drift** electronics. This sync aligns each branch with the
matching `protodunevd` drift.

## 2. Response-config sync vs `protodunevd`

`protodunevd` (ProtoDUNE-VD) is the freshly-tuned reference:

| Item | Before | protodunevd | After | Risk |
|---|---|---|---|---|
| Field response | `dunevd-resp-isoc3views-18d92.json.bz2` | `protodunevd_FR_imbalance3p_260501.json.bz2` | **synced** | MED |
| `response_plane` | `18.92 cm` | `18.1 cm` | **synced** | MED |
| Elec — top (`active_cru='tde'`) | `JsonElecResponse` top file, `postgain 1.0` | `JsonElecResponse` top file, `postgain 1.36` | **synced** | LOW |
| Elec — bottom (else) | `ColdElecResponse`, `postgain 1.1365`, no explicit gain | `ColdElecResponse`, `gain 7.8 mV/fC`, `postgain 1.0` | **synced** | LOW |
| Extra gain (top `postgain`) | `1.0` | `1.36` | **synced** | LOW |
| Noise spectra | `protodune-noise-spectra-v1.json.bz2` (single; MicroBooNE-derived, self-flagged "bogus") | top `pdvd-top-noise-spectra-v3` / bottom `pdvd-bottom-noise-spectra-7d8mVfC-v1` | **synced** (now an `active_cru` conditional, mirroring `elec`) | LOW |

## 3. What was changed (2026-05-18)

`params.jsonnet` only:

- `det.response_plane`: `18.92 → 18.1 cm`.
- `files.fields[0]`: → `protodunevd_FR_imbalance3p_260501.json.bz2`.
- `elec` — top branch: `postgain 1.0 → 1.36`; added explicit
  `gain: 14 mV/fC` placeholder. Bottom branch: `postgain 1.1365 → 1.0`; added
  explicit `gain: 7.8 mV/fC`; `shaping 2.2 µs` unchanged.
- `files.noise`: single legacy file → an `active_cru` conditional
  (top `pdvd-top-noise-spectra-v3.json.bz2`,
  bottom `pdvd-bottom-noise-spectra-7d8mVfC-v1.json.bz2`).
- `adc.resolution: 12` **restored** — see the breakage note below.

### Pre-existing breakage repaired as a side effect

`dune-vd-coldbox` **did not compile at all** before this work. Commit `41e02736`
removed the `elec.gain` and `adc.resolution` defaults from the shared base
`pgrapher/common/params.jsonnet`, and this directory had never restored either
(the same regression that hit `dune10kt-1x2x6` and `uboone`). Every entry config
failed with `RUNTIME ERROR: Field does not exist: gain` / `... resolution`.

The electronics sync (§3) adds an explicit `gain` to both `elec` branches, which
clears the `gain` error; `adc.resolution: 12` was additionally restored (the
value the sibling `dunevd-crp2` and `dune-vd` configs use; pre-`41e02736` the
base supplied 12). `adc.resolution` is **not** one of the four protodunevd sync
items — it is a breakage repair, restoring the pre-`41e02736` default, recorded
here for transparency.

## 4. Assumptions & open questions

- The `active_cru == 'tde'` ⇒ top-drift / else ⇒ bottom-drift convention is
  pre-existing in this file and is taken as correct.
- The `protodunevd_FR_imbalance3p_260501` field response is a ProtoDUNE-VD CRP
  Garfield response, assumed applicable to the CRP1 cold box (same CRP design).

## 5. Out of scope

Not touched: `chndb-*` (including the handwritten `chndb-resp.jsonnet`),
`nf.jsonnet`, `sp-filters.jsonnet`. Only the four response items were synced
(plus the `adc.resolution` breakage repair).

## 6. Verification

`wcsonnet` (`WIRECELL_PATH = wire-cell-data:cfg`), checked with both
`active_cru=tde` and `active_cru=bottom`:

| Entry config | Status |
|---|---|
| `wct-sim-check.jsonnet` | ✅ compiles (both `active_cru` values) |
| `wcls-sim-drift-simchannel.jsonnet` | ✅ compiles (both `active_cru` values) |
| `wcls-nf-sp.jsonnet` | ✅ evaluates past every real error; stops only at the LArSoft wrapper extVar `geo_planeid_labels` (supplied by the dunereco fcl at runtime) |

"Compiles" means the Jsonnet evaluates; it is not a run-verified test. The
referenced data files were confirmed present in `wire-cell-data`.
