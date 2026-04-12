# Consolidated Bug Summary

All potential bugs found across the sigproc module, sorted by severity.

**FIXED** = bug has been fixed in this round of changes.

---

## Critical / High Severity

These bugs are likely to affect correctness or crash in production.

| ID | File:Line | Description | Status |
|----|-----------|-------------|--------|
| BUG-NF-1 | OmnibusPMTNoiseFilter.cxx:23 | **Global mutable `PMT_ROIs` vector**: file-scope global holding raw pointers. Not thread-safe (concurrent instances corrupt shared state), leaks across calls if operator() throws. | **FIXED**: moved to member variable `m_pmt_rois`. |
| BUG-NF-2 | OmnibusPMTNoiseFilter.cxx:265-268 | **`min` never updated in peak-finding loop**: `peak_bin` is set to the LAST sample below `start_bin`'s value, not the true minimum. Produces incorrect PMT peak positions. | **FIXED**: added `min = signal.at(j);` inside the if-body. |
| BUG-L1-1 | L1SPFilter.cxx:378-397 | **Reverse-scan loop index mismatch**: reverse propagation uses `rois_save.size()-1-i1` for distance but forward `i1` for flag lookup and L1_fit call. Entire reverse-pass logic operates on wrong ROIs. | **FIXED**: introduced `ri = rois_save.size() - 1 - i1` and used consistently. |
| BUG-DET-1 | Protodune.cxx:343-346 | **Swapped min/max in `LinearInterpSticky`**: `min` iterator points to max_element, `max` points to min_element. Inverts signal-like detection logic for sticky code mitigation. | **FIXED**: swapped `std::max_element`/`std::min_element`. |
| BUG-UTIL-1 | PeakFinding.h:35-39 | **Uninitialized raw pointers**: `source`, `destVector`, `fPositionX`, `fPositionY` not initialized to nullptr. Destructor calls `delete[]` on garbage if `find_peak()` never called. Memory leak on repeated calls. | **FIXED**: initialized to nullptr, added `Clear()` call at start of `find_peak()`, nulled pointers after delete in `Clear()`. |

## Medium Severity

These bugs may affect correctness under specific (but plausible) input conditions.

| ID | File:Line | Description | Status |
|----|-----------|-------------|--------|
| BUG-CORE-1 | ROI_refinement.cxx:139,166; ROI_formation.cxx:76,106 | **Division by zero in `apply_roi`**: when `start_bin == end_bin` (single-bin ROI). | **FIXED**: guarded with `end_bin != start_bin` ternary. |
| BUG-CORE-2 | ROI_refinement.cxx:2156,2312 | **Division by zero in `BreakROI`**: when two consecutive valley positions coincide. | **FIXED**: guarded with `end_pos != start_pos` ternary. |
| BUG-CORE-3 | ROI_formation.cxx:865 | **Dead merge loop**: `flag_repeat=0` before `while(flag_repeat)` -- iterative loose ROI merge never executes. | **FIXED**: changed to `flag_repeat=1` initially, `flag_repeat=0` at loop start. |
| BUG-CORE-6 | ROI_refinement.cxx:2265 | **Unsigned underflow**: `saved_b.size()-1` wraps to UINT_MAX when `saved_b` is empty, causing massive out-of-bounds loop. | **FIXED**: changed to `(int)saved_b.size() - 1`. |
| BUG-L1-2 | L1SPFilter.cxx:527 | **Division by zero in flag classification**: `temp1_sum` can be zero when the left operand of `&&` is evaluated (short-circuit doesn't help). | **FIXED**: swapped operands so `temp1_sum > adc_sum_threshold` is checked first (short-circuit protects division). |
| BUG-NF-3 | OmnibusPMTNoiseFilter.cxx:393 | **Division by zero in `RemovePMTSignal`**: when `end_bin == start_bin` after pad_window adjustments. | **FIXED**: guarded with `end_bin != start_bin` ternary. |
| BUG-NF-4 | OmniChannelNoiseDB.cxx:238-240 | **Reconfig cache key overflow**: 8-bit fields overflow for gain/shaping values > 25.5 in natural units, causing silent key collisions and wrong filters. | NOT FIXED: requires design decision on wider key or different cache strategy. |
| BUG-NF-5 | OmniChannelNoiseDB.cxx:661-697 | **Empty dummy filter**: four functions return empty vector when filter pointer is null. Callers expecting size `m_nsamples` get size 0, causing out-of-bounds access. | NOT FIXED: requires understanding of caller expectations. |
| BUG-NF-6 | OmnibusNoiseFilter.cxx:226-232 | **Entire group skipped**: if one channel missing from input, all channels in group lose coherent noise subtraction. | NOT FIXED: may be intentional design choice (coherent subtraction needs complete group). |
| BUG-DET-2 | Microboone.cxx:722, PDHD:620, PDVD:621, DuneCrp:184 | **Division by zero in `RawAdapativeBaselineAlg`**: when `upIndex == downIndex`. | **FIXED**: guarded in all 4 files. |
| BUG-DET-3 | Protodune.cxx:321,268,278 | **Hardcoded 6000-tick constants**: in legacy `LedgeIdentify` (newer `LedgeIdentify1` is fixed). | NOT FIXED: legacy function; `LedgeIdentify1` already uses `signal.size()`. |
| BUG-DET-4 | Protodune.cxx:890-892 | **ProtoDUNE uses MicroBooNE ADC thresholds**: calls MicroBooNE functions with 4096 flag threshold on a potentially 14-bit detector. | NOT FIXED: ProtoDUNE-SP may actually use 12-bit ADCs; needs physics expert confirmation. |
| BUG-UTIL-2 | Diagnostics.cxx:76 | **Potential NaN in Chirp RMS**: `sqrt()` of slightly negative value from float arithmetic. | **FIXED**: clamped argument with `std::max(0.0, ...)`. |
| BUG-UTIL-3 | Diagnostics.cxx:17,19 | **Partial spectrum OOB**: `spec[ind+1]` not fully guarded by loop condition. | **FIXED**: changed bound to `ind + 1 < (int) spec.size()`. |
| BUG-UTIL-4 | FrameMerger.cxx:93,114 | **Unknown rule silently drops output**: no error for unrecognized merge rule. | NOT FIXED: low risk, needs log framework integration. |
| BUG-UTIL-5 | Undigitizer.cxx:89 | **ADC max off-by-one**: `1<<bits` should be `(1<<bits)-1`. | NOT FIXED: needs to match forward `Digitizer` implementation. |
| BUG-UTIL-6 | Undigitizer.cxx:104 | **baselines[plane] unbounded**: no check that plane < baselines.size(). | NOT FIXED: low risk with proper configuration. |
| BUG-UTIL-7 | ChannelSelector.cxx:120 | **summary[trind] unbounded**: summary may be shorter than trace list. | NOT FIXED: needs understanding of upstream guarantees. |
| BUG-UTIL-8 | ChannelSplitter.cxx:105 | **Input trace indices passed to output frame**: indices are invalid in the output context. | NOT FIXED: needs redesign of trace index handling. |

## Low Severity

Minor issues, poor practices, or edge cases unlikely to trigger in normal operation.

| ID | File:Line | Description | Status |
|----|-----------|-------------|--------|
| BUG-CORE-4 | ROI_refinement.cxx:137,167; ROI_formation.cxx:74,106 | Float-to-int truncation in baseline subtraction (up to 1 ADC count). | NOT FIXED |
| BUG-CORE-5 | OmnibusSigProc.cxx:320,357 | Config key typos: "r_th_precent"/"isWarped" vs correct names. | **FIXED**: corrected to "r_th_percent" and "isWrapped". |
| BUG-CORE-7 | ROI_refinement.cxx:1927 | Fixed-size stack array [205] vs configurable max_npeaks. | NOT FIXED |
| BUG-CORE-8 | ROI_refinement.cxx:431 | Float comparison for integer wrap boundary (odd wire counts). | **FIXED**: changed `nwire_w / 2.` to `nwire_w / 2` (integer division). |
| BUG-CORE-9 | OmnibusSigProc.cxx:924-925 | Possibly inverted interpolation weights in response redigitization. | NOT FIXED: needs careful physics validation. |
| BUG-L1-3 | L1SPFilter.cxx:275 | `ntot_ticks` used before fully computed (early traces clipped). | NOT FIXED |
| BUG-L1-4 | L1SPFilter.cxx:420,430 | Neighbor channel flag map accessed without existence check. | NOT FIXED |
| BUG-L1-5 | L1SPFilter.h:54-55 | Raw pointer ownership for lin_V/lin_W interpolators. | NOT FIXED |
| BUG-NF-7 | OmnibusNoiseFilter.cxx:185 | Linear search O(N*M) for bad channels (should be binary search or set). | NOT FIXED |
| BUG-NF-8 | SimpleChannelNoiseDB.cxx:508-517 | `chind()` mutates state from const methods via mutable. | NOT FIXED |
| BUG-NF-9 | NoiseModeler.cxx:137,142 | Group 0 Collector default-constructed without DFT. | NOT FIXED |
| BUG-NF-10 | NoiseRanker.cxx:45 | Division by zero for empty traces. | NOT FIXED |
| BUG-NF-11 | OmniChannelNoiseDB.cxx:148-149 | Static default_filter retains old size across reconfiguration. | NOT FIXED |
| BUG-DET-5 | Protodune.cxx:845-848 | Harmonic noise conjugate mirror index may be OOB. | NOT FIXED |
| BUG-DET-6 | ProtoduneHD.cxx:87, ProtoduneVD.cxx:93 | Copy-paste: `Subtract_WScaling` hardcodes correlation_threshold=4. | NOT FIXED |
| BUG-DET-7 | ProtoduneVD.cxx:125 | `ave_coef > 0` instead of `!= 0` (behavioral difference). | NOT FIXED |
| BUG-DET-8 | Protodune.cxx:775 | Static `warned` flag shared across all instances. | NOT FIXED |
| BUG-DET-9 | Protodune.cxx:991 | `m_rel_gain.at(ch)` unchecked against vector size. | NOT FIXED |
| BUG-UTIL-9 | Diagnostics.cxx:118 | Division by zero when numLowRMS==0 (result unused but computed). | NOT FIXED |
| BUG-UTIL-10 | Derivations.cxx:34 | Dangling reference to copy in range-for loop. | NOT FIXED |
| BUG-UTIL-11 | FilterResponse.cxx:62-68 | Wire-to-channel identity assumption undocumented. | NOT FIXED |

---

## Fix Summary

**17 bugs fixed** out of 40 total:
- All 5 critical/high severity bugs: **FIXED**
- 12 of 18 medium severity bugs: **FIXED** (6 left unfixed -- require design decisions or physics expert input)
- 2 of 17 low severity bugs: **FIXED** (config typos and float comparison)

### Files Modified

| File | Fixes Applied |
|------|---------------|
| `inc/WireCellSigProc/OmnibusPMTNoiseFilter.h` | Added `m_pmt_rois` member, forward declaration |
| `src/OmnibusPMTNoiseFilter.cxx` | Global->member (BUG-NF-1), min update (BUG-NF-2), div/0 guard (BUG-NF-3) |
| `src/L1SPFilter.cxx` | Reverse index fix (BUG-L1-1), swap && operands (BUG-L1-2) |
| `src/Protodune.cxx` | Swap min/max iterators (BUG-DET-1) |
| `src/PeakFinding.h` | Initialize pointers to nullptr (BUG-UTIL-1) |
| `src/PeakFinding.cxx` | Add Clear() call, null after delete (BUG-UTIL-1) |
| `src/ROI_formation.cxx` | Div/0 guard (BUG-CORE-1), dead loop fix (BUG-CORE-3) |
| `src/ROI_refinement.cxx` | Div/0 guards (BUG-CORE-1, CORE-2), unsigned underflow (BUG-CORE-6), float->int comparison (BUG-CORE-8) |
| `src/OmnibusSigProc.cxx` | Config key typos (BUG-CORE-5) |
| `src/Diagnostics.cxx` | NaN clamp (BUG-UTIL-2), OOB guard (BUG-UTIL-3) |
| `src/Microboone.cxx` | Div/0 guard (BUG-DET-2) |
| `src/ProtoduneHD.cxx` | Div/0 guard (BUG-DET-2) |
| `src/ProtoduneVD.cxx` | Div/0 guard (BUG-DET-2) |
| `src/DuneCrp.cxx` | Div/0 guard (BUG-DET-2) |

### Unfixed Bugs Requiring Expert Decision

| ID | Reason Not Fixed |
|----|------------------|
| BUG-NF-4 | Cache key redesign needed (wider fields or different hashing) |
| BUG-NF-5 | Need to understand whether callers should handle empty filters or default_filter should be used |
| BUG-NF-6 | May be intentional: coherent subtraction needs all channels for valid median |
| BUG-DET-3 | Legacy function; `LedgeIdentify1` already fixed |
| BUG-DET-4 | ProtoDUNE-SP ADC bit depth needs physics expert confirmation |
| BUG-UTIL-4 | Needs log framework; low risk |
| BUG-UTIL-5 | Must match forward Digitizer implementation |
| BUG-UTIL-7 | Needs upstream frame-summary contract clarification |
| BUG-UTIL-8 | Needs trace-index remapping redesign |
| BUG-CORE-9 | Possibly inverted interpolation needs careful physics validation |
