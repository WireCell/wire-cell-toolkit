# Noise Modeling Group - Code Examination

## Summary of Findings

| File | Bugs | Efficiency | Notes |
|------|------|------------|-------|
| EmpiricalNoiseModel | 2 medium, 2 minor | 1 medium | Complex two-layer caching, potential division-by-zero |
| GroupNoiseModel | 1 minor | Clean | Static dummy returns, clean resampling pipeline |
| AddNoise / AddGroupNoise | 1 medium (legacy) | 1 medium (legacy) | AddNoise refactored well; AddGroupNoise is legacy with static vars |
| NoiseSource | 1 medium, 1 minor | 1 medium | Legacy code with static state, config key typo |
| Noise (base) | Clean | Clean | Well-structured template base |
| SilentNoise | 1 minor | Clean | Raw `new` allocation |
| PerChannelVariation | 1 minor | 1 minor | Ignores m_truncate config, uses cerr |

**Critical/high findings: 0**
**Medium findings: 6**
**Minor findings: 6**

---

## EmpiricalNoiseModel

**Files:** `gen/inc/WireCellGen/EmpiricalNoiseModel.h`, `gen/src/EmpiricalNoiseModel.cxx`

### Algorithm Overview

Models per-channel noise spectra based on empirically measured data. Noise amplitude spectra are loaded from a JSON file indexed by plane and wire length. At query time, the model:

1. Maps channel ID to a wire length (summing segment lengths from the anode).
2. Interpolates the noise spectrum between two bracketing wire-length entries.
3. Adjusts for per-channel gain and shaping time deviations (via `IChannelStatus`).
4. Adds a frequency-independent constant noise term in quadrature.
5. Caches results using a packed key encoding (plane, wire-length, gain, shaping).

### Potential Bugs

**B1 (Medium): Possible division by zero in shaping-time response ratio (line 409)**

```cpp
for (size_t ind = 0; ind < namps; ++ind) {
    amp[ind] *= resp1->second[ind] / resp2->second[ind];
}
```

If `resp2->second[ind]` is zero (which happens at high frequencies where the electronics response decays to zero), this produces infinity or NaN. The amplitude spectrum `amp` would then contain corrupted values. A guard `if (resp2->second[ind] > threshold)` or clamping is needed.

**B2 (Medium): Bit-packing overflow in PWGS cache key (lines 419-432)**

```cpp
((p & 0b111) << 29)           // 3 bits for plane
| ((w & 0b111111111111) << 16) // 12 bits for wire length
| ((g & 0xFF) << 8)            // 8 bits for gain
| (s & 0xFF)                   // 8 bits for shaping
```

The plane field uses 3 bits (values 0-7) shifted to bits 29-31 of a 32-bit `unsigned int`. The total is 3+12+8+8 = 31 bits. This fits in 32 bits. However, the comment says "0 - 65535 length steps" for the wire length field but only 12 bits are allocated (max 4095). Long wires (e.g. 800 cm / 1 cm resolution = 800 steps) fit fine for current detectors, but the comment is misleading and could mask a future issue. More importantly, if `ch_gain/m_gres` or `ch_shaping/m_sres` exceed 255, the `& 0xFF` mask will silently alias different gain/shaping values to the same cache key, producing incorrect noise spectra for affected channels.

**B3 (Minor): Memory leak of NoiseSpectrum objects (line 194)**

```cpp
NoiseSpectrum* nsptr = new NoiseSpectrum();
```

Raw `new` is used without corresponding `delete`. The `m_spectral_data` map stores raw pointers. If `configure()` is called multiple times, `m_spectral_data.clear()` at line 192 leaks all previously allocated `NoiseSpectrum` objects.

**B4 (Minor): Recursive tail call for cache population (line 438)**

```cpp
m_channel_amp_cache[packkey] = amp;
return channel_spectrum(chid);
```

After populating the cache, the function calls itself recursively to perform the cached lookup. This works correctly but is fragile -- if the cache insertion fails for any reason, it would infinite-loop. A direct return of `m_channel_amp_cache[packkey]` would be safer and avoids the extra map lookups.

### Efficiency Concerns

**E1 (Medium): Redundant map lookups after insertion (lines 318-322, 391-406)**

```cpp
amp_cache[ilen] = interpolate_wire_length(iplane, ilen * m_wlres);
// ...
lenamp = amp_cache.find(ilen);
```

After inserting into `amp_cache`, the code does a second `find()`. Using the return value of `insert()` or `emplace()` would avoid the redundant lookup. Same pattern at lines 391-406 with `m_elec_resp_cache`.

### Key Algorithmic Details

- **Resampling** (line 94-151): Linear interpolation between irregularly-spaced frequency bins of the input spectrum to the regular FFT frequency grid. Amplitude scaling by `sqrt(N_new / N_old)` preserves noise power spectral density under FFT length change.
- **Wire-length interpolation** (line 234-280): Linear interpolation between the two nearest wire-length entries. Flat extrapolation (clamp) for out-of-bounds lengths.
- **Electronics response correction** (lines 380-411): When per-channel shaping differs from nominal, multiplies by ratio of new/old electronics response magnitudes in frequency domain.
- **Constant noise** (lines 414-417): Added in quadrature: `sqrt(amp^2 + const^2)`, representing amplifier-independent (white) noise.
- **Frequency grid** (lines 57-71): `gen_elec_resp_default()` builds the frequency array with positive frequencies for indices `[0, N/2]` and mirrored (Nyquist-folded) frequencies for `[N/2+1, N-1]`. The `<=` comparison at line 63 with `N/2.` means the Nyquist bin is included in the "positive" half.

---

## GroupNoiseModel

**Files:** `gen/inc/WireCellGen/GroupNoiseModel.h`, `gen/src/GroupNoiseModel.cxx`

### Algorithm Overview

Provides noise spectra indexed by "group ID" rather than by individual channel. Each group of channels shares the same noise spectrum. This enables coherent noise modeling where channels in the same group (e.g., sharing a regulator or ASIC) have correlated noise.

Configuration loads two data structures:
1. A channel-to-group mapping (JSON with `group`/`groupID` and `channels` arrays).
2. A set of spectra (JSON with `freqs`, `amps`, `nsamples`, `period` per group).

Spectra are resampled from their original sampling to the desired `nsamples`/`tick` using `Spectrum::resample()` with Hermitian mirror symmetry.

### Potential Bugs

**B1 (Minor): Static local `dummy` in const methods (lines 165, 193)**

```cpp
static amplitude_t dummy;
```

Two different `static` locals named `dummy` in `channel_spectrum()` and `group_spectrum()`. If these methods are called from multiple threads, the returned reference to a static empty vector is safe (it's never written after initialization), but it prevents detecting the "channel not found" condition at the call site since any empty result is ambiguous.

### Efficiency Concerns

None significant. The `irrterp` interpolation and `Spectrum::resample` are done once at configuration time.

### Key Algorithmic Details

- **Spectrum resampling pipeline** (lines 127-152): Reads irregular frequency/amplitude pairs up to Nyquist, interpolates to a regular grid of `norig` bins using `irrterp`, applies Hermitian mirror symmetry, then resamples to the target `nsamples` size with a tick-ratio scaling factor.
- **Group fallback** (lines 89, 101-109): Supports both `"group"` and `"groupID"` keys in JSON, and falls back to sequential counting if neither is present.

---

## AddNoise (IncoherentAddNoise / CoherentAddNoise)

**Files:** `gen/inc/WireCellGen/AddNoise.h`, `gen/src/AddNoise.cxx`

### Algorithm Overview

Frame filters that add noise to existing traces. Two variants:

- **IncoherentAddNoise**: Each channel gets an independent noise waveform generated from its per-channel spectrum. Uses recycled random normals for performance.
- **CoherentAddNoise**: Channels in the same group share the same noise waveform (generated once per group from the group spectrum).

Both use `Aux::RandTools::GeneratorN` which generates noise by: multiplying the amplitude spectrum element-wise by random normal samples, then performing an inverse FFT.

### Potential Bugs

**B1 (Minor): `static bool warned` not thread-safe (lines 62-63, 152)**

```cpp
static bool warned = false;
static bool warned2 = false;
```

These are function-local statics shared across all instances and calls. In a multi-threaded context, concurrent writes are a data race (though the consequence is merely a duplicated warning, not data corruption).

### Efficiency Concerns

**E1: Coherent noise copies per channel (line 195)**

```cpp
Waveform::increase(charge, gwlu[grpid]);
```

The `gwlu[grpid]` lookup is done via `operator[]` on every channel in the group. Using `find()` once and caching the iterator or reference would avoid repeated hash lookups. However, the hash lookup is O(1) amortized, so this is minor.

### Key Algorithmic Details

- **Recycled randoms** (lines 57-58): `Normals::make_recycling` creates a buffer of `2*nsamples` normal random values. On each call, only `2*rep_percent` fraction is refreshed and the starting index is randomized. This trades strict independence for speed.
- **sqrt(2/pi) factor** (line 70): The `sqrt2opi` scaling converts from the amplitude spectrum (RMS per frequency bin) to the sigma parameter needed by `GeneratorN`. This is the ratio `E[|X|] / sqrt(E[X^2])` for a standard normal.
- **Bug 202 compatibility** (line 57): A configuration knob `bug202` can restore legacy behavior where the recycling percentage was incorrectly set, introducing spurious coherence between channels.
- **Multiple models** (line 79): Supports layering multiple independent noise models, each contributing additively to the same channel's waveform.

---

## AddGroupNoise (Legacy)

**Files:** `gen/inc/WireCellGen/AddGroupNoise.h`, `gen/src/AddGroupNoise.cxx`

### Algorithm Overview

An older, standalone implementation of coherent (group) noise addition. Unlike `CoherentAddNoise`, it does not use the `NoiseBase` framework. It directly loads spectra and channel-group mappings from JSON files, generates one complex-valued noise spectrum per group, and applies it via inverse FFT.

### Potential Bugs

**B1 (Medium): Static local vectors cause cross-frame and cross-instance state leakage (lines 92-95)**

```cpp
static std::vector<double> random_real_part;
static std::vector<double> random_imag_part;
```

These `static` vectors persist across calls and across instances. They are resized but never cleared between frames (only `m_grp2noise` is cleared at line 88). Combined with the `resize()` call (which preserves existing elements if the size doesn't change), the random values from the previous frame partially carry over. This is intentional for the "recycling" optimization but the implementation is unsafe in a multi-threaded context and couples separate instances.

**B2 (Minor): Hermitian symmetry construction accesses out-of-bounds at boundary (lines 108-111)**

```cpp
for (int i = int(value.size()) / 2; i < int(value.size()); i++) {
    noise_freq.at(i).real(noise_freq.at(int(value.size()) - i).real());
    noise_freq.at(i).imag((-1) * noise_freq.at(int(value.size()) - i).imag());
}
```

When `i == int(value.size()) / 2`, the mirror index is `value.size() - value.size()/2 = value.size()/2` (for even sizes), so it reads from itself after being zeroed at line 106-107. This is correct but yields zero for that bin. When `i == value.size() - 1`, the mirror index is 1, which is valid.

However, the loop range starts at `value.size()/2` which was just zeroed at lines 106-107, and the first iteration mirrors index `value.size()/2` from itself (zero). This means the Nyquist bin always contributes zero, which is standard practice.

### Efficiency Concerns

**E1 (Medium): Full copy of noise_freq per channel (line 122)**

```cpp
WireCell::Waveform::compseq_t noise_freq = m_grp2noise[groupID];
```

This copies the entire complex vector for every channel, then computes an inverse FFT on each copy. Since channels in the same group share the same spectrum, the inverse FFT result is identical for all of them. The code should compute `inv_c2r` once per group and reuse the resulting waveform. This is a significant inefficiency since FFTs are expensive.

**E2 (Minor): `spec_freq` loaded but unused (line 72)**

```cpp
auto spec_freq = jdata["freqs"];
```

The frequency data is loaded from JSON but never used.

### Key Algorithmic Details

- **Direct frequency-domain noise generation**: Multiplies random normals by amplitude spectrum, constructs Hermitian-symmetric complex spectrum, then does `inv_c2r`.
- **No resampling**: Unlike `GroupNoiseModel`, this component does not resample the input spectra to match the desired tick/nsamples. The spectra must already match.

---

## NoiseSource (Legacy)

**Files:** `gen/inc/WireCellGen/NoiseSource.h`, `gen/src/NoiseSource.cxx`

### Algorithm Overview

A frame source (not a filter) that generates frames filled entirely with noise. It iterates over all channels in an anode plane, generates a noise waveform per channel using `generate_spectrum()`, and assembles them into a frame. Produces frames from `start_time` to `stop_time` at `readout_time` intervals, then sends EOS.

The file header warns this is "crufty old code" that should be rewritten.

### Potential Bugs

**B1 (Medium): Configuration key typo (line 154)**

```cpp
m_nsamples = get<int>(cfg, "m_nsamples", m_nsamples);
```

The config key is `"m_nsamples"` but `default_configuration()` at line 122 sets `cfg["nsamples"]`. The prefix `m_` is the C++ member naming convention leaking into the config key. As a result, `m_nsamples` can never be set from configuration -- it always retains its default value of 9600.

**B2 (Minor): `generate_spectrum` has static state (lines 42-44)**

```cpp
static std::vector<double> random_real_part;
static std::vector<double> random_imag_part;
```

Same issue as `AddGroupNoise` -- static vectors shared across all call sites, not thread-safe.

### Efficiency Concerns

**E1 (Medium): `generate_spectrum` applies shift incorrectly, doing two loops (lines 72-82)**

```cpp
for (int i = shift; i < int(spec.size()); i++) {
    // ...
    noise_freq.at(i - shift).real(...);
}
for (int i = 0; i < shift; i++) {
    noise_freq.at(i + int(spec.size()) - shift).real(...);
}
```

This implements a circular shift of the random buffer to decorrelate channels. The two-loop approach works but the shift is applied to the _output_ index while the random values are read at the _loop_ index, which means the amplitude spectrum bins are being permuted. This mixes frequency bins with the wrong amplitude, subtly distorting the noise power spectrum. Each bin `noise_freq[k]` gets amplitude `spec[k]` but random values from index `k+shift`. Since the randoms are all i.i.d. normal, the statistical properties are preserved, but the implementation is confusing.

### Key Algorithmic Details

- **Frame generation loop**: Produces one frame per call until `m_time >= m_stop`, then sends EOS (nullptr). After EOS, returns `false` to signal exhaustion.
- Uses `cerr` for debug output instead of the logging framework.

---

## Noise (NoiseBase / NoiseBaseT)

**Files:** `gen/inc/WireCellGen/Noise.h`, `gen/src/Noise.cxx`

### Algorithm Overview

Base class template providing shared configuration for all noise-adding components. Handles:
- IRandom and IDFT component lookup
- Noise model name(s) parsing (supports single string or array, plus plural "models" key)
- `nsamples`, `replacement_percentage`, and `bug202` parameters

`NoiseBaseT<IModel>` is a CRTP-style template that additionally resolves model type-names to `IModel::pointer` instances during `configure()`.

### Potential Bugs

None found. Clean, well-structured base class.

### Efficiency Concerns

None. Configuration-time only.

### Key Algorithmic Details

- **Multiple model support** (lines 54-68): Accepts `"model"` as string or array, and additionally checks `"models"` (plural) for an array. All are merged into `m_model_tns` set.
- **Bug 202 emulation** (lines 75-78): When `bug202` config is non-zero, passes it as the recycling percentage to `Normals::make_recycling`, reproducing a historical bug where the recycling parameter was incorrectly passed as 0.04 instead of 0.0.

---

## SilentNoise

**Files:** `gen/inc/WireCellGen/SilentNoise.h`, `gen/src/SilentNoise.cxx`

### Algorithm Overview

A trivial frame source that produces frames with zero-valued traces. Used as a placeholder or for testing. Produces `noutputs` frames (or infinite if 0), each with `nchannels` flat-zero traces. After exhausting outputs, sends one EOS then returns false.

### Potential Bugs

**B1 (Minor): Raw `new` for SimpleFrame (line 58)**

```cpp
auto sfout = new SimpleFrame(m_count, m_count * 5.0 * units::ms, traces);
```

Uses raw `new` then wraps in `IFrame::pointer` at line 64. This works correctly but is inconsistent with the rest of the codebase which uses `make_shared`. If an exception were thrown between allocation and the `IFrame::pointer` assignment, the memory would leak.

### Efficiency Concerns

None significant -- this is a trivial component.

### Key Algorithmic Details

- **EOS protocol** (lines 44-52): When `m_count == m_noutputs`, outputs nullptr (EOS) and returns true. On the next call (`m_count > m_noutputs`), returns false to signal permanent exhaustion. This two-phase termination is the standard WCT protocol.
- **Channel IDs** (line 56): Channels are numbered 0 through `nchannels-1`. No anode-plane aware channel mapping.

---

## PerChannelVariation

**Files:** `gen/inc/WireCellGen/PerChannelVariation.h`, `gen/src/PerChannelVariation.cxx`

### Algorithm Overview

A frame filter that applies per-channel electronics response variation. For each trace, it "replaces" the nominal electronics response with a per-channel calibrated response using deconvolution/reconvolution in the frequency domain (`DftTools::replace`). This simulates misconfigured or varying amplifiers.

### Potential Bugs

**B1 (Minor): `m_truncate` config is loaded but never used (line 75, 107)**

```cpp
m_truncate = cfg["truncate"].asBool();  // line 75
// ...
wave.resize(charge.size());             // line 107
```

The `m_truncate` member is set during configuration but line 107 unconditionally truncates the output to the input size. The commented-out code at lines 103-104 shows this was once used with `Waveform::replace_convolve` which accepted a truncation flag. After the refactor to `DftTools::replace`, the truncation flag was lost and the output is always truncated. If `m_truncate` is false, the user expects the extended waveform but gets truncated output.

### Efficiency Concerns

**E1 (Minor): `channel_response()` called per trace per frame (line 102)**

```cpp
Waveform::realseq_t tch_resp = m_cr->channel_response(chid);
```

The per-channel response is fetched (and copied) for every trace on every frame. If the response is static, caching the result per channel would avoid repeated lookups and copies.

### Key Algorithmic Details

- **Replace convolution**: `DftTools::replace(dft, signal, new_response, old_response)` deconvolves `old_response` from the signal and convolves with `new_response` in the frequency domain. This is equivalent to `IFFT(FFT(signal) * FFT(new_resp) / FFT(old_resp))`.
- **Response matching** (lines 64-71): The per-channel response binning must exactly match the configured `tick`. No interpolation or resampling is performed.
- Uses `std::cerr` for warnings (lines 62, 86, 90) instead of the logging framework, inconsistent with other components.

---

## Cross-Cutting Observations

1. **Legacy vs. modern split**: `AddGroupNoise` and `NoiseSource` are legacy implementations with static variables, raw pointer management, and `cerr` output. `IncoherentAddNoise`, `CoherentAddNoise`, and `GroupNoiseModel` represent the modern approach using `NoiseBase`, `RandTools`, and proper logging.

2. **Static variable thread-safety**: Three files (`AddGroupNoise`, `NoiseSource`, `AddNoise`) use function-local `static` variables for random number recycling or warning suppression. These are data races in multi-threaded execution.

3. **sqrt(2/pi) constant**: Used in `AddNoise` (line 70), `AddGroupNoise` (line 102), and `NoiseSource` (line 74). This converts between the amplitude spectrum (which represents the expected magnitude of complex Gaussian noise per bin) and the standard deviation parameter. The factor arises from `E[|Z|] = sigma * sqrt(2/pi)` for a complex Gaussian with independent real/imaginary parts each having variance `sigma^2`.

4. **Header guard typo**: `EmpiricalNoiseModel.h` line 11 uses `EMPERICALNOISEMODEL` (misspelling of "empirical"). Harmless but notable.
