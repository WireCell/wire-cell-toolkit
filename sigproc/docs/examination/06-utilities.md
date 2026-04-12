# Utility Components Examination

Files examined:
- `src/PeakFinding.cxx` (474 lines), `inc/WireCellSigProc/PeakFinding.h`
- `src/ChannelSelector.cxx` (194 lines)
- `src/ChannelSplitter.cxx`
- `src/DBChannelSelector.cxx`
- `src/FrameMerger.cxx` (147 lines)
- `src/FrameSplitter.cxx`
- `src/Diagnostics.cxx` (144 lines)
- `src/Derivations.cxx` (~70 lines)
- `src/Undigitizer.cxx`
- `src/HfFilter.cxx`, `src/LfFilter.cxx`, `src/FilterResponse.cxx`

---

## Algorithm Details

### PeakFinding

A port of ROOT's `TSpectrum::SearchHighRes` algorithm:

1. **SNIP background estimation** (lines 143-212): Statistics-sensitive Non-linear
   Iterative Peak-clipping. Iteratively clips the spectrum to estimate the smooth
   background component.
2. **Markov chain smoothing** (lines 218-296): Optional pass applying statistical
   smoothing based on neighbor ratios (exponential weighting), followed by a second
   background subtraction.
3. **Gold deconvolution** (lines 298-396): Iterative Richardson-Lucy deconvolution
   against a Gaussian response kernel of configurable width `sigma`.
4. **Peak search** (lines 415-458): Local maxima in the deconvolved spectrum are
   identified, filtered by threshold on both deconvolved and original spectra,
   sorted by descending original amplitude.

Used by `ROI_refinement::BreakROI` to identify peaks within wide ROIs for splitting.

### ChannelSelector / ChannelSplitter / DBChannelSelector

- **ChannelSelector**: Filters a frame to keep only traces whose channel IDs are
  in a configured `std::set<int>`. Supports tag-based selection with tag rules.
- **ChannelSplitter**: Splits a frame into N output frames based on anode plane
  channel assignments. Each anode's channels go to a separate output port.
- **DBChannelSelector**: Extends ChannelSelector by querying a channel noise
  database for the channel set (selects "bad" or "misconfigured" channels).

### FrameMerger / FrameSplitter

- **FrameMerger**: Joins two input frames using configurable merge map.
  Two rules: `"replace"` (frame 1 wins) and `"include"` (concatenate all traces).
- **FrameSplitter**: Trivially duplicates input frame pointer to both output ports
  (no deep copy).

### Diagnostics

- **Partial** (lines 9-32): Checks if frequency spectrum has a "partial" pattern --
  DC+1 bin dominates among first `nfreqs` bins, average power exceeds `maxpower`.
- **Chirp** (lines 34-139): Detects chirping noise (Mike Mooney's algorithm).
  Divides signal into windows, computes RMS per window, identifies contiguous
  low-RMS regions, determines chirp based on "normal neighbor" fraction.

### Derivations

- **CalcRMS**: Two-pass robust RMS (sigma-clipping at 4.5 sigma).
- **CalcMedian**: Per-tick median across channels, excluding outliers > 5*avg_RMS.

### Undigitizer

Converts ADC values back to voltage:
```
Vout = ADC/ADCmax * (fullscale_high - fullscale_low) + fullscale_low
Vin  = (Vout - baseline_for_plane) / gain
```
Clamps ADC to `[0, 2^bits]`.

### Filters

- **HfFilter**: Generalized Gaussian low-pass: `exp(-0.5 * (freq/sigma)^power)`.
  Optional DC suppression.
- **LfFilter**: High-pass: `1 - exp(-(freq/tau)^2)`. Suppresses low frequencies.
- **FilterResponse**: Loads per-channel filter response from JSON file.

---

## Potential Bugs

### BUG-UTIL-1: Uninitialized raw pointers in PeakFinding (HIGH)
**File:** PeakFinding.h:35-39

`source`, `destVector`, `fPositionX`, `fPositionY` are raw `double*` with **no
initialization to nullptr**. The destructor calls `delete[]` on all four. If
`find_peak()` was never called, these are uninitialized garbage, causing undefined
behavior (likely crash) on destruction.

Additionally, if `find_peak()` is called twice, old arrays are never freed (memory
leak). Fix: initialize to `nullptr` in constructor, call `Clear()` at start of
`find_peak()`.

### BUG-UTIL-2: Potential NaN in Chirp RMS (MEDIUM)
**File:** Diagnostics.cxx:76

```cpp
runningAmpRMS = std::sqrt(runningAmpRMS - runningAmpMean * runningAmpMean);
```
Due to floating-point arithmetic, the argument can be slightly negative for
constant-value windows, producing NaN. Should clamp: `std::sqrt(std::max(0.0, ...))`.

### BUG-UTIL-3: Partial spectrum out-of-bounds risk (MEDIUM)
**File:** Diagnostics.cxx:17,19

Loop accesses `spec[ind + 1]` where `ind` goes up to `nfreqs`. If `nfreqs + 1 >= spec.size()`,
this is out of bounds. The loop condition checks `ind < spec.size()` but the access
is `spec[ind + 1]`, which is one past what the condition guards.

### BUG-UTIL-4: Unknown merge rule silently produces empty output (MEDIUM)
**File:** FrameMerger.cxx:93,114

If `rule` is neither `"replace"` nor `"include"`, no traces are added and the
output frame is silently empty. No error message or fallback. Should have an
`else` with error logging.

### BUG-UTIL-5: ADC max off-by-one in Undigitizer (MEDIUM)
**File:** Undigitizer.cxx:89

```cpp
const float adcmax = 1<<bits;  // 4096 for 12-bit
```
Maximum valid ADC is 4095 (0 to 2^12-1). Line 102 clamps with `std::min(adc, adcmax)`,
allowing ADC=4096 which is out-of-range. The conversion gives `Vout = fullscale_high`
for ADC=4096, exceeding the digitizer's actual full scale.

### BUG-UTIL-6: Undigitizer baselines not bounds-checked against planes (MEDIUM)
**File:** Undigitizer.cxx:104

`baselines[plane]` where `plane` comes from `anode->resolve(chid).index()`. No check
that `plane < baselines.size()`.

### BUG-UTIL-7: ChannelSelector summary index not bounds-checked (MEDIUM)
**File:** ChannelSelector.cxx:120

```cpp
auto threshold = summary.size() ? summary[trind] : -999;
```
`trind` may exceed `summary.size()` if the summary vector has a different size
than the tagged traces vector.

### BUG-UTIL-8: ChannelSplitter passes input-frame trace indices to output frame (MEDIUM)
**File:** ChannelSplitter.cxx:105

`traces` (from `in->tagged_traces(inttag)`) contains indices into the **input**
frame's trace vector. These are passed to the output frame where they are invalid
(the output has a different, smaller trace vector).

### BUG-UTIL-9: Diagnostics division by zero if numLowRMS==0 (LOW)
**File:** Diagnostics.cxx:118

```cpp
normalNeighborFrac = ((double) numNormalNeighbors) / ((double) numLowRMS);
```
Produces NaN/inf when `numLowRMS == 0`. The guard at line 119 prevents it from
being used, but the division still executes.

### BUG-UTIL-10: Derivations dangling reference (LOW)
**File:** Derivations.cxx:34

```cpp
for (auto it : chansig)  // copies map entry
    WireCell::IChannelFilter::signal_t& signal = it.second;  // ref to copy's .second
```
Iterating by value, then taking reference to the copy's `.second`. Should be
`auto& it : chansig`.

### BUG-UTIL-11: FilterResponse wire-to-channel identity assumption (LOW)
**File:** FilterResponse.cxx:62-68

`channel_response(int channel_ident)` looks up `m_cr[channel_ident]`, but `m_cr`
was populated with wire indices from JSON. If wire numbers differ from channel IDs,
lookup fails.

---

## Efficiency Issues

### EFF-UTIL-1: PeakFinding raw `new[]`/`delete[]` (LOW)
Large heap allocations using raw pointers. `std::vector<double>` would be
exception-safe and equivalent in performance.

### EFF-UTIL-2: ChannelSelector uses `std::set` for channel lookup (LOW)
O(log n) per lookup. `std::unordered_set` would be O(1).

### EFF-UTIL-3: ChannelSplitter uses `std::map` for channel-to-port (LOW)
Same: `std::unordered_map` would be faster.

### EFF-UTIL-4: Derivations CalcMedian per-bin temp vector allocation (LOW)
Creates a new `temp` vector per time bin in the inner loop. Should pre-allocate
and reuse.
