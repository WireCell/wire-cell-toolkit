# `dunevd-crp2` — response-config sync status

WireCell config for the **DUNE Vertical-Drift CRP2 cold-box** test stand — a
single Charge Readout Plane (`dunevdcrp2-wires-larsoft-v1`, 4 split-anode
volumes). Builds on the shared base `pgrapher/common/params.jsonnet`.

Last reviewed / updated: **2026-05-18**.

---

## 1. What this is

A single-CRP test config. It runs one electronics response (no top/bottom
split) and is configured with the **top-drift** `JsonElecResponse` — so for the
purpose of this sync it is treated as a **top-drift CRP** (see §4).

## 2. Response-config sync vs `protodunevd`

`protodunevd` (ProtoDUNE-VD) is the freshly-tuned reference. The four
response-related items were synced to its **top-drift** branch:

| Item | Before | protodunevd (top) | After | Risk |
|---|---|---|---|---|
| Field response | `dunevd-resp-isoc3views-18d92.json.bz2` | `protodunevd_FR_imbalance3p_260501.json.bz2` | **synced** | MED |
| `response_plane` | `18.92 cm` | `18.1 cm` | **synced** | MED |
| Electronics response | `JsonElecResponse`, top-CRP file | `JsonElecResponse`, top-CRP file | unchanged (already top) | — |
| Extra gain (`postgain`) | `1.0` | `1.36` | **synced** | LOW |
| Noise spectra | `protodune-noise-spectra-v1.json.bz2` (MicroBooNE-derived, self-flagged "bogus") | `pdvd-top-noise-spectra-v3.json.bz2` | **synced** | LOW |

## 3. What was changed (2026-05-18)

`params.jsonnet` only — four values:

- `det.response_plane`: `18.92 → 18.1 cm` (must track the new field-response file).
- `files.fields[0]`: → `protodunevd_FR_imbalance3p_260501.json.bz2`.
- `elec.postgain`: `1.0 → 1.36` — the "extra gain" (comment: 11 mV/fC, 1.94 → 14 mV/fC).
- `files.noise`: → `pdvd-top-noise-spectra-v3.json.bz2`.

`elec.gain: 14 mV/fC` was **kept** — it is a jsonnet-required placeholder
(`JsonElecResponse` loads its response from file, but `tools.jsonnet` still
reads the `gain` field, and the shared base no longer supplies a default).

## 4. Assumptions & open questions

- **`dunevd-crp2` is treated as a top-drift CRP.** It already used the top-drift
  `JsonElecResponse` file, so it was synced to `protodunevd`'s top branch. If
  CRP2 actually represents a *bottom*-drift CRP, switch the elec to
  `protodunevd`'s bottom branch (`ColdElecResponse`, `gain 7.8 mV/fC`,
  `postgain 1.0`) and the noise to `pdvd-bottom-noise-spectra-7d8mVfC-v1.json.bz2`.
- The `protodunevd_FR_imbalance3p_260501` field response is a ProtoDUNE-VD CRP
  Garfield response; it is assumed applicable to the CRP2 cold box (same CRP
  design). Flag for an owner if the CRP2 geometry differs.

## 5. Out of scope

Not touched: `chndb-*`, `nf.jsonnet`, `sp-filters.jsonnet`, and the handwritten
`chndb-resp.jsonnet`. Only the four response items above were synced.

## 6. Verification

`wcsonnet` (`WIRECELL_PATH = wire-cell-data:cfg`). All three entry configs
evaluate:

| Entry config | Status |
|---|---|
| `wct-sim-check.jsonnet` | ✅ compiles |
| `wcls-nf-sp.jsonnet` | ✅ compiles (with LArSoft extVars supplied) |
| `wcls-sim-drift-simchannel.jsonnet` | ✅ compiles (with LArSoft extVars supplied) |

"Compiles" means the Jsonnet evaluates; it is not a run-verified test. The five
referenced data files were confirmed present in `wire-cell-data`.
