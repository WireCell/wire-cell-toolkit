# Algorithm Documentation for sig/

## Overview

The `sig` package implements **2D deconvolution** for Liquid Argon TPC (LArTPC) signal processing. The core idea: raw digitized waveforms from TPC wires are the convolution of the true ionization signal with detector response functions (field response, electronics response). Deconvolution in the frequency domain recovers the original signal.

The "2D" refers to performing the deconvolution in both **time** and **wire** (spatial) dimensions simultaneously, as opposed to the more common 1D (time-only) approach. This exploits the fact that the field response couples adjacent wires -- an ionization charge drifting past one wire induces signals on neighboring wires. The 2D approach deconvolves this cross-wire coupling.

### Data Representation

Data flows as `ITensorSet` objects containing 2D arrays:
- **Rows** = wire channels
- **Columns** = time ticks
- Additional tensors carry channel lists, summary info, and channel masking maps (CMM)

### Pipeline Architecture

The signal processing is split into two sequential `ITensorSetFilter` stages:

```
Input (raw waveforms, real, time domain)
    |
    v
[Decon2DResponse] -- 2D FFT deconvolution, response removal, wire/time shifts
    |
    v
(intermediate: complex, frequency domain along time, real along wire)
    |
    v
[Decon2DFilter] -- time-domain software filtering, inverse FFT, baseline restoration, CMM masking
    |
    v
Output (deconvolved waveforms, real, time domain)
```

---

## Decon2DResponse: Detailed Algorithm

File: `sig/src/Decon2DResponse.cxx`

### Step 1: Build the Overall Response Function (`init_overall_response`)

**Purpose:** Construct the combined detector response (field + electronics) at the coarse digitization sampling, for one wire plane.

1. **Load field response:** Retrieve the fine-grained field response from `IFieldResponse`. This describes the current induced on each wire as a function of time for charge drifting at various impact positions (distances from the wire center).

2. **Average over impact positions:** `Response::wire_region_average(fr)` averages the field response over all impact positions within each wire region. This produces one response waveform per wire (typically ~10-20 wires capturing the cross-wire induction pattern).

3. **Generate electronics response:** Create a `ColdElec` (cold electronics) response with configured gain and shaping time. This models the ASIC front-end amplifier/shaper. The response is scaled by `inter_gain * ADC_mV * (-1)` to convert from charge (fC) to ADC counts, with the sign flip accounting for signal polarity.

4. **Convolve field and electronics in frequency domain:**
   - FFT the averaged field response along time (axis 1)
   - FFT the electronics response
   - Multiply element-wise: `c_data(i,j) = field_fft(i,j) * elec_fft(j) * fine_period`
   - The `fine_period` factor accounts for the discrete-to-continuous normalization
   - Inverse FFT back to time domain

5. **Apply fine time offset:** Circular shift the response in time to account for sub-tick timing offsets (e.g., the T0 of the response relative to the digitization clock).

6. **Redigitize:** The field response is defined on a fine time grid (typically ~0.1 us bins). The detector digitizes at a coarser rate (typically ~0.5 us). Linear interpolation resamples the convolved response onto the coarse time grid.

### Step 2: Per-Channel Electronics Response Correction (optional)

**Purpose:** Correct for channel-to-channel variations in the electronics response.

If a `PerChannelResponse` is configured:
- For each channel, retrieve its measured electronics response
- FFT both the measured per-channel response and the nominal `ColdElec` response
- In frequency domain, divide out the measured response and multiply in the nominal one: `c_data *= nominal_elec / measured_elec`
- This normalizes all channels to the same effective electronics transfer function

### Step 3: 2D FFT Deconvolution

This is the core of the algorithm:

1. **Forward FFT on time (axis 1):** Transform input waveforms to frequency domain along time. After this, each column represents a frequency bin.

2. **Forward FFT on wire (axis 0):** Transform along the wire dimension. Now the data is in full 2D frequency space (wire-frequency x time-frequency).

3. **Build 2D response in same frequency space:**
   - Place the overall response waveforms (from Step 1) into an array with the same dimensions as the data
   - Forward FFT on time, then forward FFT on wire
   - This gives the 2D transfer function H(k_wire, k_time)

4. **Deconvolution by division:**
   ```
   S(k_w, k_t) = D(k_w, k_t) / H(k_w, k_t)
   ```
   where D is the measured data and S is the recovered signal in 2D frequency space.

5. **NaN/Inf cleanup:** Replace any NaN or Inf values from division by zero with 0.

6. **Apply wire-direction software filter:** Multiply by a high-frequency filter (`HfFilter`) along the wire dimension. This suppresses high-spatial-frequency noise amplified by the deconvolution. The filter type depends on the plane (induction vs collection).

### Step 4: Inverse FFT and Shifts

1. **Inverse FFT on wire (axis 0):** Return to wire-space (still in time-frequency domain).

2. **Inverse FFT on time (axis 1):** Return to full time-wire space.

3. **Circular shift in wire:** Shift by `(n_response_wires - 1) / 2` to center the deconvolved signal. The response function is defined with the primary wire at one end; this re-centers it.

4. **Circular shift in time:** Apply coarse time offset + intrinsic time offset (origin distance / drift speed) to align the deconvolved signal with the correct absolute time.

5. **Forward FFT on time again:** The output is left in the time-frequency domain (complex). This is because the downstream `Decon2DFilter` will apply additional frequency-domain filters before the final inverse FFT.

---

## Decon2DFilter: Detailed Algorithm

File: `sig/src/Decon2DFilter.cxx`

### Input

Receives the output of `Decon2DResponse`: a 2D complex array where the wire dimension is in real space and the time dimension is in frequency space.

### Step 1: Apply Time-Domain Software Filters

```
for each configured filter:
    wave = filter.filter_waveform(nticks)
    combined_filter *= wave  (element-wise)
```

Multiple `IFilterWaveform` components are chained multiplicatively. These typically include:
- Low-pass filters to suppress high-frequency noise
- Deconvolution-specific filters (e.g., Wiener filters) that depend on expected signal-to-noise

The combined filter is applied to each row (channel) identically:
```
c_data(row, col) *= combined_filter[col]
```

### Step 2: Inverse FFT on Wire, then Time

1. **Inverse FFT on wire (axis 1):** This was a `c2r` transform -- the input is the complex frequency-domain data along time, and this step inverts the wire-direction FFT from `Decon2DResponse`.
   
   *Note:* Looking more carefully, the code does `inv_c2r(m_dft, c_data_afterfilter, 1)` which is inverse complex-to-real along axis 1 (time). The wire-direction inverse was already done in `Decon2DResponse`. So the filter receives data that is real-along-wire, complex-along-time, and this step converts complex-along-time back to real.

2. **Extract the valid wire region:** `r_data = tm_r_data.block(m_pad_nwires, 0, m_nwires, m_nticks)` strips any wire padding. (Currently `m_pad_nwires=0`, so this is a no-op.)

### Step 3: Baseline Restoration (`restore_baseline`)

File: `sig/src/Util.cxx`

After deconvolution, the signal may have a floating baseline (DC offset) due to filtering artifacts. The baseline restoration uses a robust iterative median approach:

1. **Collect non-zero samples:** For each wire, gather all non-zero time-tick values. Zero values are assumed to be masked/dead regions and excluded.

2. **First median:** Compute the median of non-zero values. This gives a rough baseline estimate robust to signal outliers.

3. **Trim outliers:** Keep only values within +/- 500 (ADC units) of the first median. This removes large signal pulses that would bias the baseline.

4. **Second median:** Compute the median of the trimmed set. This is the final baseline estimate.

5. **Subtract:** Subtract the baseline from all non-zero samples. Zero samples remain zero (preserving the masking).

The hardcoded threshold of 500 assumes specific ADC units/scale. The two-iteration median is a form of **iterative sigma clipping** that converges quickly for typical LArTPC signals where true signals are sparse (few ticks have large signals) and baseline fluctuations are small.

### Step 4: Apply Channel Masking Map (CMM)

If CMM data is present (from noise filtering or known bad channels):

```
for each (channel, time_range) in CMM:
    r_data[channel, tmin:tmax] = 0
```

This zeros out known-bad regions, preventing noise or artifacts from propagating into downstream reconstruction (e.g., hit finding, clustering).

### Step 5: Package Output

The deconvolved real-valued waveforms are packaged back into an `ITensorSet` with appropriate metadata (tag, padding info, time bin) and the channel list tensor.

---

## Key Physical Parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `gain` | 14.0 mV/fC | Cold electronics amplifier gain |
| `shaping_time` | 2.2 us | Electronics shaping time constant |
| `inter_gain` | 1.2 | Intermediate gain stage factor |
| `ADC_mV` | 4096/(2000 mV) | ADC conversion factor |
| `ftoffset` | 0.0 us | Fine time offset for response alignment |
| `ctoffset` | -8.0 us | Coarse time offset |
| `tick` | 0.5 us | Digitization period |
| `iplane` | 0 | Wire plane index (0,1=induction, 2=collection) |

---

## Mathematical Summary

Let:
- `d(w, t)` = measured data (wire w, time t)
- `h(w, t)` = detector response (field convolved with electronics)
- `s(w, t)` = true signal (ionization)
- `n(w, t)` = noise

The measurement model: `d = h ** s + n` (where `**` is 2D convolution)

In 2D Fourier space: `D(k,f) = H(k,f) * S(k,f) + N(k,f)`

The deconvolution estimate: `S_hat(k,f) = D(k,f) / H(k,f) * F(k,f)`

where `F(k,f)` is the combined software filter (wire filter from `HfFilter` and time filters from `IFilterWaveform`). The filters serve as regularization to prevent noise amplification at frequencies where `H` is small.

The per-channel electronics correction normalizes the effective `H` to be uniform across channels before the 2D deconvolution, which assumes a spatially uniform response.
