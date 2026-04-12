# Noise Filtering Subsystem Examination

Files examined:
- `src/OmnibusNoiseFilter.cxx` (294 lines)
- `src/OmnibusPMTNoiseFilter.cxx` (405 lines)
- `src/OmniChannelNoiseDB.cxx` (702 lines)
- `src/SimpleChannelNoiseDB.cxx` (522 lines)
- `src/NoiseModeler.cxx` (257 lines)
- `src/NoiseRanker.cxx` (47 lines)
- Corresponding headers

---

## Algorithm Details

### OmnibusNoiseFilter

Three-pass pipeline on each frame:

1. **Per-channel filters** (`m_perchan`): Applied independently to each channel.
   Typical: chirp detection, sticky-code repair, RC correction, spectral filtering,
   baseline subtraction.

2. **Grouped filters** (`m_multigroup_chanfilters`): Applied to groups of channels
   sharing electronics. Waveforms copied into a `channel_signals_t` map, the filter
   modifies them jointly (coherent noise subtraction), results copied back.

3. **Per-channel status** (`m_perchan_status`): Final pass for channel quality
   determination (e.g., noisy channel tagging).

Bad channels from the noise database are zero-filled before processing.

### OmnibusPMTNoiseFilter

Removes PMT cross-talk into TPC wires:

1. Categorize channels by wire plane, compute per-channel RMS.
2. Identify PMT ROIs on collection wires: scan for large negative excursions
   below `-threshold * rms`. Contiguous regions exceeding `min_window_length`
   become PMTNoiseROI objects.
3. Identify PMT signals on induction wires at the known peak times.
4. Sort associated wires into spatially contiguous groups.
5. Remove: replace signal in ROI region with linear interpolation between
   best-baseline endpoints.

### OmniChannelNoiseDB

Full JSON-configurable per-channel noise database. Stores per channel:
- Scalar parameters: baseline, gain correction, response offset, RMS cuts
- Frequency-domain filter spectra: RCRC, electronics reconfig, frequency masks, response

Features filter caching: RCRC by rounded RC constant, reconfig by packed int key,
response by wire plane ID. Supports multiple channel specification formats in JSON.

### SimpleChannelNoiseDB

Simpler programmatic alternative (registered as `testChannelNoiseDB`). Uses lazy
index allocation -- first-seen channel gets next available index. Parameters set
via C++ setter methods rather than JSON.

### NoiseModeler + NoiseRanker

NoiseModeler accumulates noise spectra across frames using a ranker to classify
traces as "noise" (rank >= threshold). The NoiseRanker computes
`rank = 1 - outliers/N` where outliers = samples with |deviation from median| > maxdev.
At finalization, mean spectra are interpolated and written to a JSON noise file.

---

## Potential Bugs

### BUG-NF-1: Global mutable `PMT_ROIs` vector (CRITICAL)
**File:** OmnibusPMTNoiseFilter.cxx:23

```cpp
std::vector<PMTNoiseROI*> PMT_ROIs;   // file-scope global
```
This is a **file-scope global variable** holding raw pointers.
- **Not thread-safe**: concurrent instances share and corrupt this vector.
- **Stateful across calls**: cleared at line 206, but if `operator()` throws, stale data persists.
- Should be local to `operator()` or a member variable.

### BUG-NF-2: `min` variable never updated in peak-finding loop (HIGH)
**File:** OmnibusPMTNoiseFilter.cxx:265-268

```cpp
float min = signal.at(start_bin);
peak_bin = start_bin;
for (int j = start_bin + 1; j != end_bin; j++) {
    if (signal.at(j) < min) peak_bin = j;   // min never updated!
```
`min` always stays at `signal.at(start_bin)`. So `peak_bin` is set to the LAST
sample below `start_bin`'s value, not the actual minimum. This produces incorrect
peak positions whenever the true minimum occurs before a later below-start sample.
Fix: add `min = signal.at(j);` inside the if-body.

### BUG-NF-3: Division by zero in RemovePMTSignal (MEDIUM)
**File:** OmnibusPMTNoiseFilter.cxx:393

When `end_bin == start_bin` (pad_window adjustments converge):
```cpp
float content = start_content + (end_content - start_content) * (j - start_bin) / (end_bin - start_bin * 1.0);
```
Note: the denominator is `end_bin - start_bin * 1.0` which equals `end_bin - start_bin`.

### BUG-NF-4: Reconfig cache key overflow (MEDIUM)
**File:** OmniChannelNoiseDB.cxx:238-240

```cpp
int key = int(round(10.0 * from_gain / (units::mV / units::fC))) << 24 |
          int(round(10.0 * from_shaping / units::us)) << 16 | ...
```
Each field gets 8 bits (0-255). Values exceeding 25.5 in natural units overflow
into adjacent fields, causing **silent cache key collisions** and returning wrong
filters. E.g., `from_gain = 30 mV/fC` gives `round(10*30) = 300 > 255`.

### BUG-NF-5: Empty dummy filter returned for null filter pointers (MEDIUM)
**File:** OmniChannelNoiseDB.cxx:661-697

Four functions return a reference to a static empty `filter_t dummy` when the
filter pointer is null. Callers expecting a spectrum of size `m_nsamples` get an
empty vector, which will cause out-of-bounds access. Contrast with `default_filter()`
which creates a properly-sized filter of all ones.

### BUG-NF-6: Entire channel group skipped if one channel missing (MEDIUM)
**File:** OmnibusNoiseFilter.cxx:226-232

If even ONE channel in a group is missing from the input frame, the ENTIRE group
is skipped for coherent noise filtering. This silently disables coherent noise
subtraction for all channels in that group. May be intentional (coherent
subtraction needs all channels for proper median) but is aggressive.

### BUG-NF-7: Linear search for bad channels (LOW)
**File:** OmnibusNoiseFilter.cxx:185

```cpp
if (find(bad_channels.begin(), bad_channels.end(), ch) == bad_channels.end())
```
O(N) per channel. `std::binary_search` or `std::unordered_set` would be better.

### BUG-NF-8: `chind()` mutates state from const methods (LOW)
**File:** SimpleChannelNoiseDB.cxx:508-517

`chind()` is `const` but mutates `m_ch2ind` (via `mutable`). Any query for an
unknown channel permanently pollutes the index map, creating entries that don't
correspond to any configured data.

### BUG-NF-9: NoiseModeler group 0 default-constructed without DFT (LOW)
**File:** NoiseModeler.cxx:137,142

Channels not in any group default to group 0, but group 0's Collector is
default-constructed (without DFT pointer or nfft) since the commented-out
`ugroups.insert(0)` line is disabled.

### BUG-NF-10: NoiseRanker division by zero for empty traces (LOW)
**File:** NoiseRanker.cxx:45

`1.0 - (double)outliers/(double)siz` -- if trace has 0 samples, divides by zero.

### BUG-NF-11: Static default_filter cache (LOW)
**File:** OmniChannelNoiseDB.cxx:148-149

```cpp
static shared_filter_t def = make_filter();
```
Initialized once on first call. If `m_nsamples` changes via reconfiguration,
the static filter retains the old size.

---

## Efficiency Issues

### EFF-NF-1: Double waveform copy per trace in OmnibusPMTNoiseFilter (MEDIUM)
Each trace's waveform is copied twice: once for RMS calculation (line 98),
once into the by-plane map (lines 105/108/111).

### EFF-NF-2: Range-for loops copy shared_ptrs/vectors/map entries (MEDIUM)
Throughout OmnibusNoiseFilter.cxx and OmnibusPMTNoiseFilter.cxx:
- `for (auto filter : m_perchan)` -- copies shared_ptr
- `for (auto group : mgcf.channelgroups)` -- copies entire channel group vector
- `for (auto cs : bychan_coll)` -- copies channel + waveform

All should use `const auto&`.

### EFF-NF-3: Waveform copied in+out of grouped filter map (MEDIUM)
**File:** OmnibusNoiseFilter.cxx:229,239-241

For every channel group, waveforms are copied into a `channel_signals_t` map
(~9600 floats per channel), the filter runs, then results are copied back.
Two full copies per channel per group.

### EFF-NF-4: Repeated `parse_channels` with wpid (LOW)
**File:** OmniChannelNoiseDB.cxx:131

Iterates ALL anode channels for each JSON `channel_info` entry that uses
wire plane ID specification. Repeated for every such entry.
