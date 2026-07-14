# `dune-vd` — response-config sync status

WireCell config for the **full DUNE Far Detector, Vertical Drift**. Two params
files:

- `params.jsonnet` — single-CRP debug geometry (36 CRMs, one `upper` drift).
- `params-10kt.jsonnet` — dual-drift FD geometry (bottom CRMs `ident < ncrm/2`,
  top CRMs `ident ≥ ncrm/2`).

Both build on the shared base `pgrapher/common/params.jsonnet`.

Last reviewed / updated: **2026-05-18**.

---

## 1. Response-config sync vs `protodunevd`

`protodunevd` (ProtoDUNE-VD) is the freshly-tuned reference. It has **two**
electronics responses and **two** noise spectra — one per drift volume. The
field-response file is the *same* for both drifts.

| Item | Before | protodunevd | After | Risk |
|---|---|---|---|---|
| Field response | `dunevd-resp-isoc3views-18d92.json.bz2` | `protodunevd_FR_imbalance3p_260501.json.bz2` | **synced** | MED |
| `response_plane` | `18.92 cm` | `18.1 cm` | **synced** | MED |
| Elec — bottom (`params-10kt`) | *(single elec, see below)* | `ColdElecResponse`, `gain 7.8 mV/fC`, `postgain 1.0` | **synced** (`elecs[0]`) | MED |
| Elec — top (`params-10kt`) | *(single elec)* | `JsonElecResponse` top file, `postgain 1.36` | **synced** (`elecs[1]`) | MED |
| Elec (`params.jsonnet`, single drift) | `ColdElecResponse`, `gain 14`, `postgain 1.1365` | top: `JsonElecResponse`, `postgain 1.36` | **synced to top** | MED |
| Extra gain (`postgain`) | `1.1365` | top `1.36` / bottom `1.0` | **synced** | LOW |
| Noise — bottom (`params-10kt`) | *(single spectrum, see below)* | `pdvd-bottom-noise-spectra-7d8mVfC-v1.json.bz2` | **synced** (`noises[0]`) | MED |
| Noise — top (`params-10kt`) | *(single spectrum)* | `pdvd-top-noise-spectra-v3.json.bz2` | **synced** (`noises[1]`) | MED |
| Noise (`params.jsonnet`, single drift) | `dunevd10kt-1x6x6-3view30deg-noise-spectra-v1.json.bz2` | top `pdvd-top-noise-spectra-v3.json.bz2` | **synced to top** | LOW |

Before this work `dune-vd` had a **single-electronics / single-noise /
single-field** pipeline; the full FD has both drifts. The sync gives the
dual-drift `params-10kt.jsonnet` the 2-element form and makes the pipeline
drift-aware, following the `protodunevd` pattern.

## 2. What was changed (2026-05-18)

**`params-10kt.jsonnet`** (dual-drift FD):
- Single `elec` → 2-element `elecs: [bottom, top]`, plus `elec: $.elecs[0]`
  (kept as the nominal for components that read `params.elec`).
- `files.fields` → the protodunevd FR listed **twice** (one per drift, so
  `tools.jsonnet` builds one PIR per drift, paired with `elecs[0]`/`elecs[1]`).
- `files.noise` (single) → `files.noises: [bottom, top]`.
- Restored `adc.resolution: 12` — pre-existing `41e02736` breakage (see §4).

**`params.jsonnet`** (single-CRP debug geometry):
- Single `elec` → protodunevd **top** elec (`JsonElecResponse`, `postgain 1.36`).
- `files.fields` → protodunevd FR; `files.noise` → `pdvd-top-noise-spectra-v3`.
- This geometry's one CRP is `upper_crp_x` and is treated as the top drift (§4).

**`sp.jsonnet`**: `elecresponse: wc.tn(tools.elec_resp)` (single) →
`wc.tn(elec_for(anode))`, where `elec_for` selects `elec_resp_for(0)` for
bottom CRMs (`ident < n_volumes/2`) and `elec_resp_for(1)` for top.
`elec_resp_for` clamps, so single-electronics params are unaffected.

**`sim.jsonnet`**: per-drift PIR selection (`tools.pirs[drift_idx(anode)]`);
the noise model is now built per drift volume with the matching spectrum
(`make_noise_model` takes an explicit spectra file); the single shared
`mega_anode` became one mega-anode per drift (`meganodes0` / `meganodes1`);
the dead `use_shared_model` toggle was removed.

**Entry configs**: `wct-sim-fans`, `wct-depo-sim-fans`, `wct-depo-sim-img-fans`
— the `files` override's `fields`/`noise` literals → protodunevd values, and
`response_plane` `18.92 → 18.1`. `wct-sim-npz` — `response_plane` `18.92 → 18.1`.

### Behavior note

For the single-drift `params.jsonnet` path the drift index is always `0`, so
the `sp.jsonnet` / `sim.jsonnet` changes are behavior-neutral apart from
component renames (`elecresp0`, `meganodes0`). The dual-drift behavior is
exercised only through `params-10kt.jsonnet`.

## 3. Pre-existing breakage repaired as a side effect

`params-10kt.jsonnet` did not compile: commit `41e02736` removed the
`adc.resolution` default from the shared base and `params-10kt` never restored
it (`params.jsonnet` already had). `adc.resolution: 12` was restored — the value
the sibling configs use, restoring the pre-`41e02736` default. This is a
breakage repair, not one of the four protodunevd sync items.

## 4. Assumptions & open questions

- **`params.jsonnet` single-CRP geometry is treated as the top drift.** Its one
  CRP is `upper_crp_x` (the upper CRP), so its single `elec`/`noise` were synced
  to `protodunevd`'s top branch. Revisit if it is meant to model the bottom CRP.
- **The protodunevd field response is assumed valid for the full FD CRP.**
  `protodunevd_FR_imbalance3p_260501` is a ProtoDUNE-VD CRP Garfield response;
  the FD is built from the same CRP design, so it is taken as transferable.
  Flag for a DUNE-VD owner if the FD CRP geometry differs.
- **The `wcls-*` (LArSoft) entry configs are NOT synced.** They take field
  response, noise and `response_plane` from external variables (`files_fields`,
  `files_noise`, `response_plane`) supplied by the dunereco fcl wrapper at
  runtime. Syncing those requires updating the dunereco `.fcl`, which is outside
  this config repository.
- **`params-10kt.jsonnet` is not driven by any checked-in entry config** — every
  `wct-*` entry config sets `ncrm = 24`, which selects `params.jsonnet`. The
  dual-drift sync was compile-verified with a standalone test (§6).

## 5. Out of scope

Not touched: `chndb-perfect.jsonnet` (borrowed from `dune10kt-1x2x6`; `dune-vd`
has no per-region coherent-noise response-kernel machinery like `protodunevd`'s
`chndb-resp-{bot,top}.jsonnet` — porting that is a separate task), `nf.jsonnet`,
`sp-filters.jsonnet`, and the stale "BIG FAT FIXME" comment at `sp.jsonnet:6`.

## 6. Verification

`wcsonnet` (`WIRECELL_PATH = wire-cell-data:cfg`):

| Entry config | Status |
|---|---|
| `wct-sim-fans.jsonnet` | ✅ compiles |
| `wct-sim-npz.jsonnet` | ✅ compiles |
| `wct-depo-sim-fans.jsonnet` | ✅ compiles |
| `wct-depo-sim-img-fans.jsonnet` | ✅ compiles |
| `wcls-nf-sp.jsonnet` | ⚠ extVar-blocked (unchanged from before the sync) |
| `wcls-sim-drift-simchannel.jsonnet` | ⚠ extVar-blocked (unchanged) |
| `wcls-sim-drift-simchannel-nf-sp.jsonnet` | ⚠ extVar-blocked (unchanged) |
| `wcls-sim-drift-simchannel-nf-sp-img.jsonnet` | ⚠ extVar-blocked (unchanged) |

The `wct-*` configs all evaluate. The dual-drift `params-10kt.jsonnet` path
(`tools` + `sim` + `sp`, 24 anodes spanning both drifts) was compile-verified
with a standalone test. The `wcls-*` configs fail only on LArSoft-wrapper
external variables (or a pre-existing wrapper-supplied-object error) — these
failures are **identical to before the sync**; nothing here introduces a new
error. "Compiles" means the Jsonnet evaluates; it is not a run-verified test.
The referenced data files were confirmed present in `wire-cell-data`.
