# Core Signal Processing Examination

Files examined:
- `src/OmnibusSigProc.cxx` (1,925 lines)
- `src/ROI_refinement.cxx` (3,155 lines)
- `src/ROI_formation.cxx` (937 lines)
- Corresponding headers in `inc/WireCellSigProc/` and `src/`

---

## Algorithm Details

### OmnibusSigProc

The master orchestrator. Key member data:
- `m_r_data[3]`: Eigen 2D arrays [wire x tick] holding per-plane waveform data
- `m_c_data[3]`: Complex FFT counterparts
- Filter names (`m_wiener_*`, `m_gauss_*`) resolved via WCT Factory

The `operator()` pipeline per plane:
1. `init_overall_response(in)` -- convolve field response with electronics, redigitize
2. `load_data(traces, plane)` -- populate m_r_data, zero bad channels
3. `decon_2D_init(plane)` -- 2D FFT deconvolution (time+wire)
4. `decon_2D_tightROI(plane)` / `decon_2D_looseROI(plane)` -- filtered versions for ROI finding
5. `ROI_formation` -- identify signal regions
6. `ROI_refinement` -- multi-stage iterative cleanup
7. Final Wiener/Gauss filtering within ROIs

### ROI Formation

Two-tier ROI identification:
- **Tight ROIs**: `threshold = th_factor * RMS + 1`. Contiguous above-threshold regions,
  padded by `pad` ticks, merged when overlapping.
- **Loose ROIs**: Rebin waveform by factor `rebin`, apply lower threshold
  `l_factor * RMS * rebin`. Extended to cover tight ROI boundaries.

RMS is computed via robust two-pass sigma-clipping: IQR-based initial estimate,
then standard deviation of samples within the IQR range.

`create_ROI_connect_info()` interpolates ROIs on gap channels between neighbors
(channels with no signal but flanked by channels that have ROIs).

### ROI Refinement

Maintains three relationship maps between SignalROI pointers:
- `front_rois[roi]`: ROIs on the next wire that overlap in time
- `back_rois[roi]`: ROIs on the previous wire that overlap in time
- `contained_rois[loose_roi]`: tight ROIs contained within a loose ROI

**Stages:**
1. `CleanUpROIs`: BFS flood-fill from tight ROIs. Any loose ROI not connected
   (directly or transitively) to a tight ROI is removed.
2. `generate_merge_ROIs`: If a tight ROI is not contained in any loose ROI,
   promote it (create a new loose ROI covering it).
3. `MP3ROI`/`MP2ROI`: Multi-plane geometric protection. Require that signal
   regions are consistent across 3 (or 2) wire planes.
4. `BreakROIs`: TSpectrum-like peak finding to split wide ROIs at valleys.
   Runs for `m_break_roi_loop` iterations (default 2).
5. `ShrinkROIs`: Tighten ROI boundaries. For each ROI, check if neighbor ROIs
   support the current extent; if not, shrink to the signal region.
6. `CleanUpCollectionROIs`/`CleanUpInductionROIs`: Remove likely-fake signals
   based on charge ratio and shape criteria.
7. `ExtendROIs`: Extend each ROI to cover the time range of overlapping
   neighbor ROIs (ensures consistent coverage for baseline subtraction).

---

## Potential Bugs

### BUG-CORE-1: Division by zero in `apply_roi` (MEDIUM)
**Files:** ROI_refinement.cxx:139,166; ROI_formation.cxx:76,106

When `start_bin == end_bin` (single-bin ROI), the baseline interpolation:
```cpp
(end_content - start_content) * (i - start_bin) / (end_bin - start_bin)
```
divides by zero. Single-bin ROIs can arise after ShrinkROI or BreakROI processing.

### BUG-CORE-2: Division by zero in `BreakROI` (MEDIUM)
**File:** ROI_refinement.cxx:2156,2312

When two consecutive valley positions coincide (`end_pos == start_pos`):
```cpp
start_content + (end_content - start_content) * (k - start_pos) / (end_pos - start_pos)
```

### BUG-CORE-3: Dead merge loop in `find_ROI_loose` (MEDIUM)
**File:** ROI_formation.cxx:865

```cpp
int flag_repeat = 0;
while (flag_repeat) {    // <-- never enters
    flag_repeat = 1;
    ...
```
`flag_repeat` initialized to 0, so the iterative ROI merge pass is never executed.
The variable name and loop structure strongly suggest this should be initialized to 1.
Adjacent loose ROIs that should be merged are left separate.

### BUG-CORE-4: Float-to-int truncation in baseline subtraction (LOW-MEDIUM)
**Files:** ROI_refinement.cxx:137,167,206; ROI_formation.cxx:74,106,136

```cpp
int content = r_data(irow, i) - (...);
```
The subtraction produces a float but is stored in `int`, silently truncating.
This introduces systematic rounding toward zero of up to 1 ADC count.

### BUG-CORE-5: Config key typos (LOW)
**File:** OmnibusSigProc.cxx

- Line 320: `"r_th_precent"` vs `"r_th_percent"` (default_configuration vs configure)
- Line 357: `"isWarped"` vs `"isWrapped"` (same mismatch)

Default config JSON has the wrong key names. Users relying on dumped defaults
will get stale keys.

### BUG-CORE-6: Unsigned underflow in `section_protected` (MEDIUM)
**File:** ROI_refinement.cxx:2265

```cpp
for (unsigned int j = 0; j < saved_b.size() - 1; ++j) {
```
If `saved_b` is empty, `saved_b.size() - 1` wraps to UINT_MAX, causing a massive
out-of-bounds loop. Triggered when a signal is entirely above threshold through
a BreakROI pass.

### BUG-CORE-7: Fixed-size stack arrays in `BreakROI` (LOW)
**File:** ROI_refinement.cxx:1927,2010

```cpp
float valley_pos[205];
```
Sized for `max_npeaks <= 200`. If `max_npeaks` is configured larger, this is
an out-of-bounds stack write.

### BUG-CORE-8: Float comparison for integer wrap boundary (LOW)
**File:** ROI_refinement.cxx:431

```cpp
if (chid != nwire_u + nwire_v + nwire_w / 2. || !isWrapped)
```
`nwire_w / 2.` produces a float. When `nwire_w` is odd, the comparison is always
true, making wrap protection ineffective.

### BUG-CORE-9: Possibly inverted interpolation weights (UNCLEAR)
**File:** OmnibusSigProc.cxx:924-925

In `init_overall_response`, the linear interpolation weights appear swapped:
the weight for sample `fcount-1` grows as `ctime` moves *away* from it.
This has been in the code for a long time and may be compensated elsewhere,
but is mathematically unusual.

### BUG-CORE-10: TestROIs boundary access (LOW -- dead code)
**File:** ROI_refinement.cxx:2901-2962

Accesses `rois_u_loose.at(chid+1)` when `chid == nwire_u - 1` and `.at(chid-1)`
when `chid == 0`. Would crash if invoked, but the function is not called in production.

---

## Efficiency Issues

### EFF-CORE-1: `cal_RMS` takes vector by value (HIGH)
**File:** ROI_formation.h:73, ROI_formation.cxx:364

```cpp
double cal_RMS(Waveform::realseq_t signal);  // copies entire waveform
```
Called once per wire per plane (thousands of times). Should be `const&` reference.

### EFF-CORE-2: Massive 3x code duplication per plane (HIGH -- maintenance)
Nearly every function in ROI_formation.cxx and ROI_refinement.cxx has identical
logic copy-pasted three times for U, V, W planes. Examples:
- `extend_ROI_self()`: 3 copies of the same loop
- `create_ROI_connect_info()`: 3 copies
- `apply_roi()`: 3 copies in both files
- `CleanUpROIs()`, `generate_merge_ROIs()`, `ExtendROIs()`: 2-3 copies each

This is a maintenance hazard (bugs must be fixed in 3 places) and increases
instruction cache pressure. A plane-parameterized helper would be better.

### EFF-CORE-3: Repeated map lookups in inner loops (MEDIUM)
**File:** ROI_refinement.cxx:291-307, ROI_formation.cxx:411-417,708-717

`bad_ch_map[irow + offset]` looked up 3 times per column in nested loops.
Each lookup is O(log n). Should cache the reference.

### EFF-CORE-4: Redundant c_data_afterfilter allocation (MEDIUM)
In `decon_2D_tightROI`, `decon_2D_looseROI`, etc., a full complex array is
allocated, filled, and IFFTed ~8 times per plane. The allocation could be reused.

### EFF-CORE-5: `get_mp2_rois()` non-const returns by value (LOW)
**File:** ROI_refinement.h:47

Returns a copy of the entire `MapMPROI` map. The const overload correctly
returns by reference. Likely should also return by reference.

### EFF-CORE-6: Redundant filter Factory lookups (LOW)
Filter waveform objects looked up via `Factory::find<IFilterWaveform>(...)` every
time `decon_2D_*` functions are invoked. These are the same named objects each time.

### EFF-CORE-7: Raw pointer ownership of SignalROI objects
All `SignalROI` objects are `new`ed and manually `delete`d. If any code path throws
or returns early, memory leaks. The `BreakROI1` and `ShrinkROI` functions delete
old ROIs and create new ones -- if construction throws, the old ROI is already
deleted and the maps are inconsistent. `std::unique_ptr` or a pool allocator
would be safer.
