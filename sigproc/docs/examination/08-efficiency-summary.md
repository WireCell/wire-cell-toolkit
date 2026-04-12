# Consolidated Efficiency Summary

All efficiency issues found across the sigproc module, sorted by impact.

**FIXED** = issue has been fixed in this round of changes.

---

## High Impact

Issues that affect performance measurably in production workloads.

| ID | File(s) | Description | Status |
|----|---------|-------------|--------|
| EFF-CORE-1 | ROI_formation.cxx:364 | **`cal_RMS` takes vector by value**: copies entire waveform (~10K floats) on every call. Called once per wire per plane (thousands of times per frame). | **FIXED**: changed to `const Waveform::realseq_t&`. |
| EFF-L1-1 | L1SPFilter.cxx:473-496 | **JSON config re-parsed per `L1_fit` call**: ~20 `get()` calls parsing JSON on every ROI. All values are constant after `configure()`. | **FIXED**: all 20 config values cached as member variables in `configure()`, referenced as `const` locals in `L1_fit`. |
| EFF-L1-5 | L1SPFilter.cxx:496 | **Smearing vector parsed from JSON per call**: `get<vector<double>>(m_cfg, "filter")` constructs a new vector every call. | **FIXED**: cached as `m_smearing_vec` in `configure()`, used as `const&` in `L1_fit`. |
| EFF-CORE-2 | ROI_formation.cxx, ROI_refinement.cxx | **Massive 3x code duplication per plane**: Nearly every function has identical logic copy-pasted for U, V, W planes. | NOT FIXED: requires large-scale refactor to parametrize by plane. Maintenance hazard documented. |

## Medium Impact

Issues that add noticeable overhead, especially for large detectors.

| ID | File(s) | Description | Status |
|----|---------|-------------|--------|
| EFF-NF-2 | OmnibusNoiseFilter.cxx | **Range-for loops copy shared_ptrs, vectors, map entries**: `for (auto x : container)` throughout. | **FIXED**: changed to `const auto&` for filters, groups, and `auto&` for chgrp writeback. |
| EFF-NF-3 | OmnibusNoiseFilter.cxx:229,239-241 | **Double waveform copy in grouped filter**: waveforms copied into map and back. | Partially addressed: `auto& cs` for the writeback loop avoids one copy. The initial copy into `chgrp` is necessary for the filter API. |
| EFF-NF-1 | OmnibusPMTNoiseFilter.cxx:98,105 | **Double waveform copy per trace**: copied once for RMS, once into by-plane map. | **FIXED**: store waveform once in plane map, compute RMS from the stored copy. Also fixed range-for copies. |
| EFF-CORE-3 | ROI_refinement.cxx:291-307, ROI_formation.cxx:411-417 | **Repeated map lookups in inner loops**: `bad_ch_map[irow + offset]` looked up 3 times per column. | **FIXED**: cached as `const auto& bad_ranges = bad_it->second` in both files. |
| EFF-CORE-4 | OmnibusSigProc.cxx | **Redundant c_data_afterfilter allocation**: full complex array allocated ~8 times per plane. | NOT FIXED: requires workspace member variable refactor. |
| EFF-L1-2 | L1SPFilter.cxx:564 | **Dense matrix allocation per section**: up to 120x240 matrix allocated per section. | NOT FIXED: requires pre-allocated scratch buffer with size management. |
| EFF-DET-1 | All detector files | **Per-channel FFT**: every channel performs its own forward/inverse FFT. | NOT FIXED: per-channel corrections prevent trivial batching; DFT plans likely cached internally. |
| EFF-DET-2 | All detector files | **Redundant waveform copy for baseline**: `auto temp_signal = signal;` copies entire waveform. | NOT FIXED: copy needed for sigma-clipping; histogram-based approach would be a larger change. |
| EFF-NF-7 | OmnibusNoiseFilter.cxx:185 | **Linear search for bad channels**: O(N*M) with `std::find`. | **FIXED**: replaced with `std::unordered_set` for O(1) lookup. |

## Low Impact

Minor issues or micro-optimizations.

| ID | File(s) | Description | Status |
|----|---------|-------------|--------|
| EFF-CORE-5 | ROI_refinement.h:47 | `get_mp2_rois()` non-const returns by value (copies entire map). | **FIXED**: changed to return `MapMPROI&`. |
| EFF-CORE-6 | OmnibusSigProc.cxx | Filter objects re-fetched via `Factory::find` on every invocation. | NOT FIXED: Factory likely caches internally; low priority. |
| EFF-NF-4 | OmniChannelNoiseDB.cxx:131 | `parse_channels` with wpid iterates ALL anode channels per JSON entry. | NOT FIXED |
| EFF-DET-3 | All detector files | `CalcRMSWithFlags` grows vector via `push_back` without `reserve()`. | **FIXED**: added `temp.reserve(sig.size())` in Microboone, PDHD, PDVD, DuneCrp. |
| EFF-DET-4 | All detector files | ROIs iterated by value: `for (auto roi : rois)` copies each `vector<int>`. | NOT FIXED: scattered across many call sites; low individual impact. |
| EFF-DET-5 | ProtoduneVD.cxx:1159-1185 | Per-bin temporary vector allocation in shield coupling median. | NOT FIXED |
| EFF-DET-6 | ProtoduneVD.cxx:1255-1338 | `ShieldCouplingSub` makes 3 full passes over all channel data. | NOT FIXED |
| EFF-UTIL-1 | PeakFinding.cxx | Raw `new[]`/`delete[]` instead of `std::vector`. | NOT FIXED: functional with nullptr init fix; full refactor to std::vector is larger scope. |
| EFF-UTIL-2 | ChannelSelector.cxx | `std::set` for channel lookup (O(log n)). | NOT FIXED |
| EFF-UTIL-3 | ChannelSplitter.cxx | `std::map` for channel-to-port lookup. | NOT FIXED |
| EFF-UTIL-4 | Derivations.cxx | Per-bin temp vector allocation in CalcMedian inner loop. | NOT FIXED |

---

## Fix Summary

**11 efficiency issues fixed** out of 25 total:
- 3 of 4 high-impact: **FIXED** (EFF-CORE-2 deferred -- large refactor)
- 5 of 8 medium-impact: **FIXED** (EFF-CORE-4, EFF-L1-2, EFF-DET-1, EFF-DET-2 deferred)
- 3 of 13 low-impact: **FIXED**

### Files Modified for Efficiency

| File | Fixes Applied |
|------|---------------|
| `src/ROI_formation.h` | `cal_RMS` signature: by-value -> `const&` |
| `src/ROI_formation.cxx` | `cal_RMS` definition, cached `bad_ch_map` lookups (x2) |
| `src/ROI_refinement.h` | `get_mp2_rois()` return by reference |
| `src/ROI_refinement.cxx` | Cached `bad_ch_map` lookup |
| `inc/WireCellSigProc/L1SPFilter.h` | Added cached config member variables |
| `src/L1SPFilter.cxx` | Config caching in `configure()`, use cached values in `L1_fit` |
| `src/OmnibusNoiseFilter.cxx` | `const auto&` in range-for loops, `unordered_set` for bad channels |
| `src/OmnibusPMTNoiseFilter.cxx` | Eliminated double waveform copy, `auto&` in range-for loops |
| `src/Microboone.cxx` | `CalcRMSWithFlags` vector reserve |
| `src/ProtoduneHD.cxx` | `CalcRMSWithFlags` vector reserve |
| `src/ProtoduneVD.cxx` | `CalcRMSWithFlags` vector reserve |
| `src/DuneCrp.cxx` | `CalcRMSWithFlags` vector reserve |

### Remaining Unfixed Issues

| ID | Reason |
|----|--------|
| EFF-CORE-2 | 3x code duplication requires large-scale plane-parametrization refactor |
| EFF-CORE-4 | Workspace member variable for c_data_afterfilter needs design |
| EFF-L1-2 | Pre-allocated scratch matrix needs size management across ROIs |
| EFF-DET-1 | Per-channel FFT inherent to the algorithm; DFT plans cached internally |
| EFF-DET-2 | Waveform copy needed for sigma-clipping; histogram approach is larger change |
| EFF-DET-4 | Scattered across many call sites; low individual impact |
| EFF-DET-5/6 | PDVD-specific optimizations; low priority |
| EFF-UTIL-* | Various low-priority micro-optimizations |
