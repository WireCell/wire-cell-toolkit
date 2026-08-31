# L1SP Filter Examination

Files examined:
- `inc/WireCellSigProc/L1SPFilter.h`
- `src/L1SPFilter.cxx` (692 lines)

---

## Algorithm

### Mathematical Formulation

The L1SP filter implements **L1-norm regularized sparse signal deconvolution**
(a LASSO / compressed sensing variant) for LArTPC wire signal recovery in
**shorted wire regions** where induction and collection plane wires are
electrically connected.

The optimization problem solved for each ROI segment:

    minimize  ||G * beta - W||_2^2  +  lambda * ||beta||_1

where:
- **W** (size N): observed ADC waveform in the ROI segment
- **beta** (size 2N): unknown signal coefficients split into:
  - `beta[0..N-1]`: collection-plane (W-plane) signal amplitudes
  - `beta[N..2N-1]`: induction-plane (V-plane) signal amplitudes
- **G** (size N x 2N): forward model encoding both detector response functions.
  Column j of the left half contains the W-plane response at each measurement time
  for a unit signal at time j. The right half contains V-plane response.
- **lambda**: L1 regularization promoting sparsity (most time bins have zero signal)

### Why L1 Sparsity?

In shorted wire regions, the measured ADC is a mix of both induction (bipolar)
and collection (unipolar) responses. Standard deconvolution fails because it
assumes a single response function. L1SP jointly fits both response types,
decomposing the mixed signal. The L1 penalty ensures the recovered signal is
sparse, which is physically expected for wire signals from ionization tracks.

### Processing Pipeline (`operator()`)

1. **ROI identification** (lines 232-329):
   - From deconvolved signal traces: identify ticks with positive charge
   - From raw ADC traces: identify ticks exceeding `raw_ROI_th_nsigma * noise_sigma`,
     padded by `raw_pad` ticks
   - Merge overlapping ROIs with gap tolerance `roi_pad`

2. **First L1 pass** (lines 342-354):
   For each ROI, call `L1_fit`. Returns flag:
   - 0 = no action
   - 1 = shorted (L1 applied successfully)
   - 2 = remove signal (artifact)
   Negative outputs are clamped to zero.

3. **Forward propagation** (lines 357-375):
   If a flag-1 (shorted) ROI is found, adjacent flag-2 ROIs within 20 ticks
   are re-fitted as shorted.

4. **Reverse propagation** (lines 377-398):
   Same logic scanning in reverse direction. **Note: has index bug (BUG-L1-1).**

5. **Cross-channel cleaning** (lines 404-448):
   Flag-2 ROIs are zeroed if neighboring channels (ch +/- 1) have a flag-1 ROI
   overlapping within 3 ticks.

6. **Output** (lines 451-458): packaged into output frame with configured tag.

### L1_fit Function (lines 468-692)

**Classification** (lines 525-536):
- `flag = 1` (do L1): `sum/rescaled_abs_sum > ratio_threshold` AND `abs_sum > sum_threshold`
  (significant net-positive signal -- likely real charge)
- `flag = 2` (zero): rescaled absolute sum below limit, or large decon signal with small raw ADC

**L1 Solve** (lines 538-591):
- ROI divided into sections of ~`l1_seg_length` ticks
- For each section: build response matrix G from `lin_W` and `lin_V` interpolators
- Response windowed to [-15us, +10us] relative to measurement time
- LASSO solver finds sparse beta coefficients

**Post-processing** (lines 593-680):
- Collection/induction components weighted by `l1_col_scale` / `l1_ind_scale`
- Signal convolved with configurable smearing filter
- Values below `l1_decon_limit / l1_scaling_factor` threshold zeroed
- Small isolated peaks (below `peak_threshold` max or `mean_threshold` average) removed

### Key Configuration Parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `l1_seg_length` | 120 | Max ticks per LASSO solve segment |
| `l1_scaling_factor` | 500 | Scales response matrix for numerical conditioning |
| `l1_lambda` | 5 | LASSO regularization strength |
| `l1_epsilon` | 0.05 | LASSO convergence tolerance |
| `l1_niteration` | 100000 | Max LASSO iterations |
| `l1_col_scale` | 1.15 | Weight for collection plane component |
| `l1_ind_scale` | 0.5 | Weight for induction plane component |
| `l1_decon_limit` | 100 | Min signal threshold (electrons) |
| `raw_ROI_th_nsigma` | 4 | N-sigma threshold for raw ROI detection |
| `roi_pad` / `raw_pad` | 3 / 15 | ROI extension padding (ticks) |
| `filter` | [] | Smearing filter kernel |
| `collect_time_offset` | 3.0 us | Collection plane time offset relative to induction |

---

## Potential Bugs

### BUG-L1-1: Reverse-scan loop index mismatch (HIGH)
**File:** L1SPFilter.cxx:378-397

The reverse propagation pass scans ROIs in reverse order (distance check uses
`rois_save.size() - 1 - i1`) but the flag lookup and `L1_fit` call use forward
index `i1`:
```cpp
if (map_ch_flag_rois[trace->channel()].at(i1) == 1) {   // uses i1, not reverse index
```
This means the reverse pass reads the flag for the **wrong ROI**. The entire
reverse-pass propagation logic is scrambled -- it checks the time distance to
one ROI but operates on a different one.

### BUG-L1-2: Division by zero in flag classification (MEDIUM)
**File:** L1SPFilter.cxx:527

```cpp
if (temp_sum / (temp1_sum * adc_sum_rescaling * 1.0 / nbin_fit) > adc_ratio_threshold &&
    temp1_sum > adc_sum_threshold) {
```
The division `temp_sum / (temp1_sum * ...)` is the **left** operand of `&&`. C++
short-circuit evaluation only skips the **right** operand. When `temp1_sum == 0`
(all ADC values below threshold), the division by zero **executes**, producing
`inf` or `NaN`. The behavior then depends on the comparison with
`adc_ratio_threshold` which is implementation-defined for NaN.

### BUG-L1-3: `ntot_ticks` used before fully computed (LOW)
**File:** L1SPFilter.cxx:275

`ntot_ticks` is updated inside the loop iterating over traces. Early traces may
see an incomplete `ntot_ticks` value for their bounds check. Only clips some
edge ticks for early-processed channels.

### BUG-L1-4: Potential out-of-bounds for neighboring channel flags (LOW)
**File:** L1SPFilter.cxx:420,430

`map_ch_flag_rois[ch + 1].at(i3)` accessed without checking that `ch+1` exists
in `map_ch_flag_rois`. The `operator[]` would default-construct an empty vector,
and `.at(i3)` would throw.

### BUG-L1-5: Raw pointer ownership (LOW)
**File:** L1SPFilter.h:54-55

`linterp<double>* lin_V` and `lin_W` are raw `new`/`delete` managed. Should be
`std::unique_ptr`.

### BUG-L1-6: Unused member `m_period` (INFO)
**File:** L1SPFilter.h:52

Set in `operator()` but never read anywhere. Dead state.

---

## Efficiency Issues

### EFF-L1-1: Config re-read on every `L1_fit` call (HIGH)
**File:** L1SPFilter.cxx:473-496

Every invocation of `L1_fit` calls `get()` on the JSON config ~20 times. These
values are constant after `configure()`. Should be cached as member variables.
For a frame with thousands of ROIs, this is significant overhead from repeated
JSON parsing.

### EFF-L1-2: Dense matrix allocation per section (MEDIUM)
**File:** L1SPFilter.cxx:564

```cpp
MatrixXd G = MatrixXd::Zero(temp_nbin_fit, temp_nbin_fit * 2);
```
For each section of each ROI of each channel, a potentially large matrix
(up to 120x240, ~230 KB) is allocated and zero-initialized. Pre-allocating
and reusing a scratch buffer would be more efficient.

### EFF-L1-3: Unnecessary waveform copy (LOW)
**File:** L1SPFilter.cxx:265

```cpp
Waveform::realseq_t tmp_charge = charges;
```
Entire charge vector copied just to compute percentiles. Could avoid if
`percentile` doesn't modify input.

### EFF-L1-4: Repeated map lookups (LOW)
**File:** L1SPFilter.cxx:347-396

`map_ch_flag_rois[trace->channel()]` looked up repeatedly in tight loops.
Should cache a reference.

### EFF-L1-5: Smearing vector parsed from JSON per call (MEDIUM)
**File:** L1SPFilter.cxx:496

```cpp
std::vector<double> smearing_vec = get<std::vector<double>>(m_cfg, "filter");
```
Parses JSON array and constructs new vector on every `L1_fit` call. Should be
cached once in `configure()`.

### EFF-L1-6: LASSO solver iteration limit (INFO)
Default `l1_niteration = 100000` could cause very long solve times for
ill-conditioned sections. No early termination reporting -- caller has no
indication if solver hit the limit without convergence.
