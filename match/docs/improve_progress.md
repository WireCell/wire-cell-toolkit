# Q-L Matching — hard-coded numbers to fix

Tracking document for the magic numbers and unit-less literals in the
charge-light (Q-L) matching code. **This is a catalog only** — no code changes
yet. Each entry records the current value, where it lives, what it physically
means, and (for later) the direction of the fix. Per the repo convention, any
future toggle/parameter must default to the *current* behavior so existing SBND
production configs stay bit-identical.

Scope: `match/src/QLMatching.cxx`, `match/inc/WireCellMatch/QLMatching.h`,
`match/src/TimingTPCBundle.cxx`, `match/src/Opflash.cxx`.

Status legend: ☐ not started · ◐ in progress · ☑ done.

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
| A1 | `QLMatching.cxx:234` | `lo_x_bound = (tpc==0) ? -2000 : 0` | drift-volume low X edge (mm) | bare mm; should come from `m_dv`/`m_anode` face X, not a literal | ☐ |
| A2 | `QLMatching.cxx:235` | `hi_x_bound = (tpc==0) ? 0 : 2000` | drift-volume high X edge (mm); `0` = cathode | bare mm; derive from geometry | ☐ |
| A3 | `QLMatching.cxx:287` | `std::abs(y) > 2000` | active-volume half-height in Y (mm) | bare mm; derive from `m_dv` Y extent | ☐ |
| A4 | `QLMatching.cxx:287` | `z < 0 \|\| z > 5000` | active-volume Z span (mm) | bare mm; derive from `m_dv` Z extent | ☐ |
| A5 | `QLMatching.cxx:291` | `std::abs(x) > 1950` | "close to PMT/anode" distance (mm) | bare mm; should be `(active_halfwidth − gap)`, configurable | ☐ |
| A6 | `QLMatching.cxx:297` | `x / 10., y / 10., z / 10.` | mm → cm for `SemiAnalyticalModel` | hard-coded `/10.`; comment even says "units::cm == 10". Use `/ units::cm` | ☐ |

Note A2 also encodes the **cathode plane at X = 0** — the same `0` is reused as
both the high edge of TPC 0 and the low edge of TPC 1. When the two-TPC joint
match is built, this single shared literal is exactly the seam that has to
become a real geometry value.

---

## B. Drift / timing / model constants in the header (bare, but configurable)

These are exposed to jsonnet, so they are less urgent, but the *defaults* are
still bare physical numbers with the units only living in a comment.

| # | File:line | Literal | Means | Issue / fix direction | St |
|---|-----------|---------|-------|-----------------------|----|
| B1 | `QLMatching.h:69` | `m_drift_speed{1.563e-3}` | drift speed (mm/ns) | bare; configs should pass `params.lar.drift_speed`; default has no `units::` | ☐ |
| B2 | `QLMatching.h:60-61` | `m_flash_mintime{-1500e3}` / `m_flash_maxtime{1500e3}` | flash time window (ns) | bare ns; express with `units::us`/`units::ns` | ☐ |
| B3 | `QLMatching.h:62-63` | `m_beam_mintime{-5e3}` / `m_beam_maxtime{5e3}` | beam-gate window (ns) | bare ns; express with `units::us` | ☐ |
| B4 | `QLMatching.h:59` | `m_flash_minPE{500}` | min flash PE gate | PE (no WCT unit); fine, but document SBND override (50) | ☐ |
| B5 | `QLMatching.h:64` | `m_QtoL{0.5}` | charge→light scale | calibration constant; document, keep configurable | ☐ |
| B6 | `QLMatching.h:73` | `m_strength_cutoff{0.05}` | LASSO keep threshold | already configurable; OK | ☑ |
| B7 | `QLMatching.cxx:97` | `vuv_absorption_length` default `85.0` | VUV absorption length (cm) | bare cm fallback when JSON omits it | ☐ |

---

## C. LASSO / cost-function weights (dimensionless, hard-coded inline)

The regression weights are declared as four inline locals and reused across both
rounds. They directly set the relative pull of the light vs. charge constraints
and the shape penalty — i.e. they tune the matching — yet none is configurable.

| # | File:line | Literal | Means | Issue / fix direction | St |
|---|-----------|---------|-------|-----------------------|----|
| C1 | `QLMatching.cxx:419` | `lambda = 0.1` | LASSO L1 regularization λ | pull to config; same value feeds round 1 & 2 | ☐ |
| C2 | `QLMatching.cxx:420` | `delta_charge = 0.01` | charge-constraint stiffness (`1/delta` rows) | config | ☐ |
| C3 | `QLMatching.cxx:421` | `delta_light = 0.025` | per-flash background-light stiffness (rd 1) | config | ☐ |
| C4 | `QLMatching.cxx:422` | `delta_shape = 0.01` | KS shape-penalty weight (rd 2) | config | ☐ |
| C5 | `QLMatching.cxx:485` | `weights(nbundle+k) = 0.5` | background-column weight (rd 1) | bare; config | ☐ |
| C6 | `QLMatching.cxx:477-479` | `0.3` (knee) and `0.3` (floor) | PE-mismatch weight: `>0.3·meas` ? ratio : `0.3` (rd 1) | two distinct 0.3's; name & config | ☐ |
| C7 | `QLMatching.cxx:562-565` | `0.3` (knee) and `0.3` (floor) | same PE-mismatch weight (rd 2) | duplicate of C6; share one constant | ☐ |
| C8 | `QLMatching.cxx:496, 579` | `initial(n) = 1.0` | LASSO warm-start strength | OK as 1.0, but note it | ☐ |

---

## D. Pre-selection / bad-match gates (QLMatching loop)

| # | File:line | Literal | Means | Issue / fix direction | St |
|---|-----------|---------|-------|-----------------------|----|
| D1 | `QLMatching.cxx:251` | `get_total_PE() > 5000` | MC saturated-PMT mask trigger (PE) | bare PE; config (MC-only) | ☐ |
| D2 | `QLMatching.cxx:293` | `npt_outside_drift > 0.25 * npt` | drop bundle if >25% pts drift out | bare fraction; config | ☐ |
| D3 | `QLMatching.cxx:321` | `get_total_pred_light() < 10` | min predicted PE to keep bundle | bare PE; config | ☐ |
| D4 | `QLMatching.cxx:323` | `get_ks_dis() == 1` | KS==1 ⇒ uncomputable ⇒ reject | sentinel, document | ☐ |
| D5 | `QLMatching.cxx:327` | `chi2/ndf > 1e4` | pre-select χ²/ndf ceiling | bare; config | ☐ |
| D6 | `QLMatching.cxx:289` | `if (std::abs(x) && ...)` | **suspected porting bug**: `std::abs(x)` is truthy for ~all points, so `flag_close_to_PMT` is set almost always; the intended boundary test looks lost. Flag for review before any tuning. | ☐ |

---

## E. Out-of-beam QA cuts (`organize_bundles`)

| # | File:line | Literal | Means | Issue / fix direction | St |
|---|-----------|---------|-------|-----------------------|----|
| E1 | `QLMatching.cxx:812` | `get_ks_dis() > 0.2` | out-of-beam KS reject | bare; config | ☐ |
| E2 | `QLMatching.cxx:812` | `chi2/ndf > 20` | out-of-beam χ²/ndf reject | bare; config | ☐ |
| E3 | `QLMatching.cxx:817` | `> 0.5 * get_total_PE()` | out-of-beam total-PE mismatch (50%) | bare fraction; config | ☐ |

---

## F. Bundle-quality constants (`TimingTPCBundle.cxx`)

These govern bundle merging and the "high-consistent" flag — the same family of
KS / χ² thresholds as §E, duplicated here as literals.

| # | File:line | Literal | Means | Issue / fix direction | St |
|---|-----------|---------|-------|-----------------------|----|
| F1 | `TimingTPCBundle.cxx:109` | `candidate_ks_dis < 0.2`, `chi2/ndf < 20` | merge-acceptance gate | duplicate of E1/E2; share constants | ☐ |
| F2 | `TimingTPCBundle.cxx:117-119` | exponent `0.8` | `ks·(chi2/ndf)^0.8` tie-break in `add_bundle` | bare; name & config | ☐ |
| F3 | `TimingTPCBundle.cxx:193` | `ks_dis < 0.06`, `ndf >= 3`, `chi2 < ndf*nvalidopdets` | "high-consistent" flag | three bare thresholds; config | ☐ |
| F4 | `TimingTPCBundle.cxx:103, 187` | `pe[j] < 1 && pred_pe[j] < 1` | per-opdet ndf inclusion knee (1 PE) | bare PE; share with F-family | ☐ |

---

## G. Optical PE-error model (`Opflash.cxx`)

| # | File:line | Literal | Means | Issue / fix direction | St |
|---|-----------|---------|-------|-----------------------|----|
| G1 | `Opflash.cxx:24` | `PE_err = (PE<1) ? 0.3 : 0.3*PE` | per-channel PE error: 0.3 floor + 30% scale | the `1` knee and both `0.3`s are bare; this error feeds every χ² in §C–F. Centralize & config | ☐ |
| G2 | `Opflash.cxx:20` | `PE_err.resize(nchan, 1)` | default PE_err = 1 before init | bare; document | ☐ |
| G3 | `QLMatching.cxx:198` | `Opflash(ff, 0.0, m_nchan)` | `threshold = 0.0` ⇒ `fired` if PE > 0 | bare threshold passed positionally; name it | ☐ |

---

## H. Hard-coded SBND layout (large literals; geometry, not tuning)

Lower priority — these are detector descriptions, not tunable cuts — but they
are still hard-coded SBND-specific data baked into C++ rather than read from
config/geometry.

| # | File:line | What | Issue / fix direction | St |
|---|-----------|------|-----------------------|----|
| H1 | `QLMatching.cxx:157-168` | 312-entry `opdet_mask` PMT on/off pattern | SBND-specific; could move to JSON/geometry | ☐ |
| H2 | `QLMatching.cxx:239-240` | even/odd `idet % 2` TPC split | bakes TPC identity into PMT index parity; blocker for two-TPC joint match | ☐ |
| H3 | `QLMatching.h:51` | `m_nchan{312}` | OpDet count default | configurable; OK, document | ☑ |
| H4 | `QLMatching.h:77-82` | 312-entry `m_VUVEfficiency` / `m_VISEfficiency` | SBND efficiency tables | overridable via config; document the source | ☐ |
| H5 | `QLMatching.cxx:30` | `kFlashGidStride = 1000000` | global-flash-id stride | deliberate, well-documented; leave as-is | ☑ |

---

## I. Cosmetic / logging literals (no physics; lowest priority)

Display-only truncation in `log->debug`: `int(x*100)/100.` (2-dp rounding) at
`:256-258, 336-340, 800-803`, and `int(flash_time)/100.` time scaling at `:256`.
These do not affect results; listed only so a future cleanup pass doesn't mistake
them for tunable parameters. **Do not touch.**

---

## Suggested order of work (for when we do touch the code)

1. **§A unit literals first** — they are the stated request and the prerequisite
   for any cross-detector / two-TPC reuse. Replace `/10.` with `/units::cm`, and
   source the X/Y/Z bounds from `m_dv`/`m_anode` instead of `±2000 / 0 / 5000`.
2. **§D6 porting bug** — confirm/repair the `std::abs(x)` boundary test before
   tuning anything that depends on `flag_close_to_PMT`.
3. **§C + §G** — promote the LASSO weights and the PE-error model to config so
   the matching can be tuned without a rebuild; collapse the duplicated `0.3`
   and KS/χ² thresholds (§E/§F) into shared named constants.
4. **§B time/speed defaults** — wrap in `units::` and wire to common params.
5. **§H** — only if/when generalizing beyond SBND.

Every item must land behind a config knob whose default reproduces today's
number exactly (bit-identical production).
