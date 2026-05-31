# Q-L Matching — hard-coded numbers to fix

Tracking document for the magic numbers and unit-less literals in the
charge-light (Q-L) matching code. Each entry records the current value, where it
lives, what it physically means, and the direction of the fix. Per the repo
convention, any new toggle/parameter defaults to the *current* behavior so
existing SBND production configs stay bit-identical.

Scope: `match/src/QLMatching.cxx`, `match/inc/WireCellMatch/QLMatching.h`,
`match/src/TimingTPCBundle.cxx`, `match/src/Opflash.cxx`.

Status legend: ☐ not started · ◐ in progress · ☑ done.

## Pass 1 — units + config (DONE)

Proper Wire-Cell units applied to every physical literal, and the tuning-relevant
numbers pulled up to jsonnet config keys on the `QLMatching` component (defaults =
the historical literals, so production is bit-identical — verified by rerunning the
`sbnd_xin` Q/L events: matched flash ids, t0, predicted PE, KS, χ²/ndf and LASSO
strengths all unchanged). Cross-class constants are forwarded from `QLMatching`:
the `Opflash` PE-error model (`Match::PEErr`) and the `TimingTPCBundle` thresholds
(`Match::BundleQualityParams` via `set_quality_params`). The full knob list lives in
`cfg/pgrapher/experiment/sbnd/qlmatching.jsonnet`.

Intentionally **not** changed (out of scope for a units/config pass): the §D6
`std::abs(x)` logic bug (preserved as-is, to keep behavior identical — fix later),
sentinels/structural constants (§D4 `KS==1`, §C8 warm-start `1.0`, §G2 default
`PE_err=1`), and the SBND layout tables / even-odd mask (§H1/H2/H4).

---

## A. Physical lengths written as bare numbers (NO Wire-Cell units)

These are the highest-priority "numbers without Wire-Cell units" the request
called out. They are physical lengths (millimetres) hard-coded as plain
literals instead of `... * units::mm`, `... / units::cm`, or — better —
derived from the anode / `IDetectorVolumes` geometry. They silently assume the
SBND active-volume dimensions and break the moment the code is reused for
another detector or the geometry shifts.

| # | File:line | Literal | Means | Issue / fix direction | St |
|---|-----------|---------|-------|-----------------------|----|
| # | Literal | Means | Resolution | St |
|---|---------|-------|------------|----|
| A1 | `lo_x_bound = (tpc==0) ? -2000 : 0` | drift-volume low X edge (mm) | `-m_x_bound` / `0`; key `x_bound` (`2000*wc.mm`) | ☑ |
| A2 | `hi_x_bound = (tpc==0) ? 0 : 2000` | drift-volume high X edge (mm); `0` = cathode | `0` / `m_x_bound`; key `x_bound`. Cathode seam stays literal `0` | ☑ |
| A3 | `std::abs(y) > 2000` | active-volume half-height in Y (mm) | `m_y_bound`; key `y_bound` (`2000*wc.mm`) | ☑ |
| A4 | `z < 0 \|\| z > 5000` | active-volume Z span (mm) | `m_z_min`/`m_z_max`; keys `z_min` (`0`), `z_max` (`5000*wc.mm`) | ☑ |
| A5 | `std::abs(x) > 1950` | "close to PMT/anode" distance (mm) | `m_pmt_dist`; key `pmt_dist` (`1950*wc.mm`) | ☑ |
| A6 | `x / 10., y / 10., z / 10.` | mm → cm for `SemiAnalyticalModel` | now `/ units::cm` (pure units, no key) | ☑ |

Note A2 also encodes the **cathode plane at X = 0** — the same `0` is reused as
both the high edge of TPC 0 and the low edge of TPC 1. It is deliberately kept a
literal `0` (the true origin). When the two-TPC joint match is built, this shared
seam is the value that has to become real geometry.

---

## B. Drift / timing / model constants in the header (bare, but configurable)

These are exposed to jsonnet, so they are less urgent, but the *defaults* are
still bare physical numbers with the units only living in a comment.

| # | Literal | Means | Resolution | St |
|---|---------|-------|------------|----|
| B1 | `m_drift_speed{1.563e-3}` | drift speed (mm/ns) | now `1.563*units::mm/units::us`; SBND passes `params.lar.drift_speed` | ☑ |
| B2 | `m_flash_mintime{-1500e3}` / `{1500e3}` | flash time window (ns) | now `∓1.5*units::ms` | ☑ |
| B3 | `m_beam_mintime{-5e3}` / `{5e3}` | beam-gate window (ns) | now `∓5*units::us` | ☑ |
| B4 | `m_flash_minPE{500}` | min flash PE gate | PE (no WCT unit); already configurable (SBND override 50) | ☑ |
| B5 | `m_QtoL{0.5}` | charge→light scale | dimensionless; already configurable | ☑ |
| B6 | `m_strength_cutoff{0.05}` | LASSO keep threshold | already configurable | ☑ |
| B7 | `vuv_absorption_length` default `85.0` (cm) | VUV absorption length | left: it is a JSON-model fallback (cm in the model's own unit system), not a WCT-units quantity | ☐ |

---

## C. LASSO / cost-function weights (dimensionless, hard-coded inline)

The regression weights are declared as four inline locals and reused across both
rounds. They directly set the relative pull of the light vs. charge constraints
and the shape penalty — i.e. they tune the matching — yet none is configurable.

| # | Literal | Means | Resolution | St |
|---|---------|-------|------------|----|
| C1 | `lambda = 0.1` | LASSO L1 regularization λ | key `lasso_lambda` (feeds both rounds via the local alias) | ☑ |
| C2 | `delta_charge = 0.01` | charge-constraint stiffness (`1/delta` rows) | key `delta_charge` | ☑ |
| C3 | `delta_light = 0.025` | per-flash background-light stiffness (rd 1) | key `delta_light` | ☑ |
| C4 | `delta_shape = 0.01` | KS shape-penalty weight (rd 2) | key `delta_shape` | ☑ |
| C5 | `weights(nbundle+k) = 0.5` | background-column weight (rd 1) | key `bkg_weight` | ☑ |
| C6 | `0.3` (knee) and `0.3` (floor) | PE-mismatch weight (rd 1) | keys `pe_mismatch_knee` / `pe_mismatch_floor` | ☑ |
| C7 | `0.3` (knee) and `0.3` (floor) | same PE-mismatch weight (rd 2) | same two keys (shared across both rounds) | ☑ |
| C8 | `initial(n) = 1.0` | LASSO warm-start strength | left as a named structural constant (not tuning) | ☐ |

---

## D. Pre-selection / bad-match gates (QLMatching loop)

| # | Literal | Means | Resolution | St |
|---|---------|-------|------------|----|
| D1 | `get_total_PE() > 5000` | MC saturated-PMT mask trigger (PE) | key `mc_saturation_pe` | ☑ |
| D2 | `npt_outside_drift > 0.25 * npt` | drop bundle if >25% pts drift out | key `drift_out_frac` | ☑ |
| D3 | `get_total_pred_light() < 10` | min predicted PE to keep bundle | key `min_pred_pe` | ☑ |
| D4 | `get_ks_dis() == 1` | KS==1 ⇒ uncomputable ⇒ reject | left: sentinel (KS default value), not tuning | ☐ |
| D5 | `chi2/ndf > 1e4` | pre-select χ²/ndf ceiling | key `preselect_chi2ndf_max` | ☑ |
| D6 | `if (std::abs(x) && ...)` | **suspected porting bug**: `std::abs(x)` is truthy for ~all points, so `flag_close_to_PMT` is set almost always. | left **unchanged** — preserving behavior; logic fix is a separate later task | ☐ |

---

## E. Out-of-beam QA cuts (`organize_bundles`)

| # | Literal | Means | Resolution | St |
|---|---------|-------|------------|----|
| E1 | `get_ks_dis() > 0.2` | out-of-beam KS reject | key `outbeam_ks_max` | ☑ |
| E2 | `chi2/ndf > 20` | out-of-beam χ²/ndf reject | key `outbeam_chi2ndf_max` | ☑ |
| E3 | `> 0.5 * get_total_PE()` | out-of-beam total-PE mismatch (50%) | key `outbeam_pe_frac` | ☑ |

---

## F. Bundle-quality constants (`TimingTPCBundle.cxx`)

These govern bundle merging and the "high-consistent" flag. Forwarded from
`QLMatching` config via `BundleQualityParams` → `set_quality_params()`.

| # | Literal | Means | Resolution | St |
|---|---------|-------|------------|----|
| F1 | `candidate_ks_dis < 0.2`, `chi2/ndf < 20` | merge-acceptance gate | keys `bundle_ks_merge_max` / `bundle_chi2ndf_merge_max` | ☑ |
| F2 | exponent `0.8` | `ks·(chi2/ndf)^0.8` tie-break in `add_bundle` | key `bundle_addmerge_exponent` | ☑ |
| F3 | `ks_dis < 0.06`, `ndf >= 3`, `chi2 < ndf*nvalidopdets` | "high-consistent" flag | keys `highconsist_ks_max` / `highconsist_min_ndf` (the `ndf*nvalidopdets` term left structural) | ☑ |
| F4 | `pe[j] < 1 && pred_pe[j] < 1` | per-opdet ndf inclusion knee (1 PE) | key `bundle_pe_ndf_knee` | ☑ |

---

## G. Optical PE-error model (`Opflash.cxx`)

Forwarded from `QLMatching` config via the `Match::PEErr` struct.

| # | Literal | Means | Resolution | St |
|---|---------|-------|------------|----|
| G1 | `PE_err = (PE<1) ? 0.3 : 0.3*PE` | per-channel PE error: 0.3 floor + 30% scale | `PEErr{floor,frac,knee}`; keys `pe_err_floor` / `pe_err_frac` / `pe_err_knee` | ☑ |
| G2 | `PE_err.resize(nchan, 1)` | default PE_err = 1 before init | left: transient pre-init fill, overwritten in the loop | ☐ |
| G3 | `Opflash(ff, 0.0, m_nchan)` | `threshold = 0.0` ⇒ `fired` if PE > 0 | key `flash_pe_threshold` | ☑ |

---

## H. Hard-coded SBND layout (large literals; geometry, not tuning)

Lower priority — these are detector descriptions, not tunable cuts — but they
are still hard-coded SBND-specific data baked into C++ rather than read from
config/geometry.

| # | What | Issue / fix direction | St |
|---|------|-----------------------|----|
| H1 | 312-entry `opdet_mask` PMT on/off pattern | **DONE** — mask now derived per-channel from the injected `OpDets` table: on iff `type ∈ active_opdet_types` (config, default `[1]` = PMTs only). The hard-coded 312-array is gone; the derived mask is byte-identical to the old one on the SBND geometry. | ☑ |
| H2 | even/odd `idet % 2` TPC split | **DONE** — TPC membership now from OpDet position (`center.x` vs `cathode_x`), the same same-TPC test the optical model uses. Reproduces the old split for SBND PMTs bit-identically (the 6 channels where parity≠position-sign are all Arapucas, masked off anyway) and is correct for any layout. | ☑ |
| H3 | `m_nchan{312}` OpDet count default | configurable; OK | ☑ |
| H4 | 312-entry `m_VUVEfficiency` / `m_VISEfficiency` | SBND efficiency tables; overridable via config | ☐ |
| H5 | `kFlashGidStride = 1000000` global-flash-id stride | deliberate, well-documented; leave as-is | ☑ |

---

## I. Cosmetic / logging literals (no physics; lowest priority)

Display-only truncation in `log->debug`: `int(x*100)/100.` (2-dp rounding) at
`:256-258, 336-340, 800-803`, and `int(flash_time)/100.` time scaling at `:256`.
These do not affect results; listed only so a future cleanup pass doesn't mistake
them for tunable parameters. **Do not touch.**

---

## Remaining work (after Pass 1)

- **§D6 porting bug** — confirm/repair the `std::abs(x)` boundary test (it makes
  `flag_close_to_PMT` true for almost all points). This is a *behavior* change, so
  it needs its own validation, separate from the units/config pass.
- **§A geometry sourcing** — Pass 1 made the X/Y/Z bounds configurable literals; a
  follow-up could source them from `m_dv`/`m_anode` so they track the real geometry
  (and the cathode seam `0`) instead of SBND-specific defaults. Prerequisite for
  cross-detector / two-TPC reuse.
- **§H** — H1 (opdet mask) and H2 (TPC split) are now derived from the injected
  `OpDets` metadata (see table above); the remaining SBND-specific item is H4
  (VUV/VIS efficiency tables), to revisit when generalizing beyond SBND.
- Minor leftovers kept as named constants by design: §B7, §C8, §D4, §G2.

Any behavior change must land behind a config knob whose default reproduces today's
number exactly (bit-identical production), validated by rerunning the Q/L events.
