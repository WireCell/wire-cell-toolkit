# Digitization & Pipeline Components -- Code Examination

## Summary of Findings

This document examines 30 files (header + implementation pairs) in the `gen/` directory covering digitization, frame/depo fan-in/out, merging, chunking, filtering, and related pipeline utilities.

### Critical Bugs

| ID | File | Line(s) | Description |
|----|------|---------|-------------|
| B1 | Digitizer.cxx | 136-151 | Skipped channels leave nullptr in `adctraces` vector |
| B2 | FrameSummer.cxx | 54-57 | `newtwo` frame constructed but never used; summing uses original `two` |
| B3 | FrameFanin.cxx | 155,163 | Unsigned underflow when `m_ft.nrules()` returns 0 |
| B4 | Scaler.cxx | 111-158 | No EOS/nullptr check on `depo` causes null-pointer dereference; off-by-one in bounds clamping |
| B5 | DepoSetFilterYZ.cxx | 111 | Double unit application on position coordinates |
| B6 | DepoSetFilterYZ.cxx / Scaler.cxx | 126-148 / 142-148 | Off-by-one: clamped index `== n_bins` accesses out-of-bounds |
| B7 | DepoChunker.cxx | 47 | EOS nullptr pushed into `m_depos` after emit, leaked to next cycle |
| B8 | DepoSetRotate.cxx | 56-66 | When `m_rotate` is false, all depos are silently dropped (empty output) |
| B9 | DepoSetFilter.cxx | 66 | Logging message computes `nout-nin` (negative) instead of `nin-nout` |
| B10 | TruthTraceID.cxx | 153-218 | Wire loop nested inside depo loop: O(depos*wires) with cumulative bindiff state corruption |
| B11 | Retagger.cxx | 119-121 | Mutates const reference `summary` from input frame |
| B12 | ZSEndedTrace.cxx | 14 | UB when `m_nbins==0`: computes `min(bin, -1)` on unsigned-promoted value |
| B13 | Scaler.cxx | 152 | Unconditional `depo->prior()->id()` dereferences potentially null prior |

### Efficiency Concerns

| ID | File | Description |
|----|------|-------------|
| E1 | Digitizer.cxx | Dense array allocation for all channels, even sparse input |
| E2 | Scaler.cxx:156 | Sorts entire output queue on every single depo insertion |
| E3 | MultiDuctor.cxx:93-105 | Json::Value access inside inner depo-processing loop |
| E4 | TruthTraceID.cxx:127-133 | Wire/time filters regenerated in `process()` every readout |
| E5 | FrameFanout.cxx:117 | Copies all traces via `*in->traces()` for each output port |

---

## Per-File Analysis

---

### Digitizer (Digitizer.h / Digitizer.cxx)

**Algorithm overview:** Converts a voltage-domain frame to ADC counts. For each trace, looks up the wire plane via the anode to find the baseline, adds `gain * voltage + baseline`, then maps through a linear ADC transfer function with configurable resolution (default 12-bit), clamping to [0, adcmaxval]. Optionally rounds or floors.

**Bug B1 -- nullptr traces in output vector (line 131-151):**
```cpp
ITrace::vector adctraces(nrows);        // line 131: pre-sized with nrows elements
for (size_t irow = 0; irow < nrows; ++irow) {
    ...
    if (!wpid.valid()) {
        log->warn("...");
        continue;                        // line 140: skips, leaving adctraces[irow] as nullptr
    }
    ...
    adctraces[irow] = make_shared<SimpleTrace>(...);  // line 151
}
```
When `wpid` is invalid for a channel, the trace at that index remains a default-constructed `shared_ptr` (nullptr). The output frame then contains null trace pointers. Downstream consumers that dereference traces without null-checking will crash.

**Efficiency E1 (line 128):** A dense `nrows x ncols` Eigen array is allocated even if the input traces are sparse. For frames with many channels but short traces, this wastes memory. The `fill()` utility scatters data into this array and then each row is re-extracted sequentially -- a direct per-trace loop would avoid the intermediate dense buffer.

**Key detail:** The `digitize()` method (line 77-93) computes `adcmaxval = (1 << m_resolution) - 1` on every call. This is cheap but could be cached.

---

### Reframer (Reframer.h / Reframer.cxx)

**Algorithm overview:** Creates a "rectangular" frame from tagged traces. Initializes a map of channel -> waveform vectors (one per anode channel, sized to `m_nticks`, filled with `m_fill`), then overlays input traces by accumulating samples. Carries forward trace tags and summaries.

**Potential issue (line 111):** `*out_it += *in_it` accumulates rather than overwrites. If multiple input traces share a channel, their samples are summed. This is intentional but may surprise users expecting the last trace to win.

**Edge case (line 99-104):** The delta_tbin logic correctly handles truncation and padding, but if `m_tbin` is negative (requesting samples before the frame's reference time), `delta_tbin` will be large positive, potentially advancing `in_it` past `in_end`. The while-loop guard (line 110) protects against this, so no actual bug, but negative `m_tbin` silently produces empty traces.

**Efficiency:** The `std::map<int, vector<float>>` for `waves` involves per-channel heap allocation and tree traversal. An `unordered_map` or a flat vector indexed by channel offset would be faster for large channel counts.

---

### FrameFanin (FrameFanin.h / FrameFanin.cxx)

**Algorithm overview:** Merges N input frames into one output frame by concatenating all traces. Supports per-port tagging and tag rule transformation. Has both fixed and dynamic multiplicity modes.

**Bug B3 -- Unsigned underflow in tag rule lookup (lines 155, 163):**
```cpp
const size_t n = m_ft.nrules("frame");
auto fo = m_ft.transform(iport >= n ? n-1 : iport, "frame", fintags);
```
If no rules are configured (`n == 0`), then `n-1` underflows to `SIZE_MAX` for `size_t`. The `transform()` call would then access an invalid rule index. The same pattern appears on line 163 for trace rules. This is only reached when tag rules are configured with zero entries of a given type, which is unlikely but not impossible.

**Memory leak (line 196):** `auto sf = new Aux::SimpleFrame(...)` uses raw `new`. While it is later wrapped in `IFrame::pointer(sf)` on line 211, if an exception is thrown between lines 196-211, the frame leaks. Should use `make_shared`.

---

### FrameFanout (FrameFanout.h / FrameFanout.cxx)

**Algorithm overview:** Fans out one frame to N output ports. In trivial mode, passes the same shared pointer to all outputs. In rule-based mode, creates new SimpleFrame per port with transformed tags.

**Efficiency E5 (line 117):** In rule-based mode, `*in->traces()` is dereferenced and copied into a new SimpleFrame for every output port. Since traces are shared_ptr-based, the copy is of the vector of pointers, not the data. Still, N copies of a potentially large vector could be avoided if the output frames could share the trace vector.

**No bugs found.** Clean implementation with proper EOS handling.

---

### FrameSummer (FrameSummer.h / FrameSummer.cxx)

**Algorithm overview:** Joins two frames using `Aux::sum()`. Supports time alignment and offset between the two frames.

**Bug B2 -- Dead code / logic error (lines 53-57):**
```cpp
auto vtraces2 = two->traces();
ITrace::vector out_traces(vtraces2->begin(), vtraces2->end());
auto newtwo = std::make_shared<SimpleFrame>(two->ident(), t2, out_traces, two->tick());

out = Aux::sum(IFrame::vector{one, two}, one->ident());  // uses original 'two', not 'newtwo'!
```
A new frame `newtwo` is constructed with the adjusted time `t2`, but the sum on line 57 uses the original `two` (with original time), discarding the alignment/offset. The `newtwo` variable is unused. This means the `align` and `offset` configuration options have no effect.

**Missing member initialization:** `m_toffset` and `m_align` are initialized in the constructor body (lines 61-64) rather than in-class, but this is fine since default_configuration also returns these defaults.

---

### FrameSummerYZ (FrameSummerYZ.h / FrameSummerYZ.cxx)

**Algorithm overview:** Sums N frames (fan-in interface, not a 2-input joiner like FrameSummer). Uses `Aux::sum()` directly on the input vector.

**Minor issue (line 39):** Error message says "FrameFanin" instead of "FrameSummerYZ":
```cpp
THROW(ValueError() << errmsg{"FrameFanin multiplicity must be positive"});
```

**No alignment/offset support.** Unlike FrameSummer, there is no time alignment logic, so all frames must already be time-aligned. This is a design limitation, not a bug.

---

### MultiDuctor (MultiDuctor.h / MultiDuctor.cxx)

**Algorithm overview:** Routes depos to different sub-ductors based on configurable rules (wire bounds or boolean). Buffers output frames and extracts them in readout-time-aligned windows. Supports continuous and triggered modes.

**Efficiency E3 (lines 93-105):** The `Wirebounds` functor accesses `Json::Value` objects (`jargs`) in the hot path of every depo dispatch. The fixme comment on line 88 acknowledges this. Converting to a pre-parsed struct of vectors in `configure()` would significantly improve throughput.

**Known limitation (line 351-358):** Two self-documented bugs: (1) a sub-ductor frame larger than `m_readout_time` after splitting can produce an oversized output frame, and (2) overlapping traces on the same channel are not merged.

**Style issue:** Extensive use of `std::cerr` instead of the logging framework. This component does not inherit from `Aux::Logger`.

**Potential issue (line 219-226):** `start_processing()` logic: when `m_continuous` is false and `m_eos` is true and the depo is non-null (line 219), it sets `m_start_time` and returns `false`. But `m_eos` is only reset by the assignment at line 215 (`m_eos = true` when depo is null). There is no code that resets `m_eos` to false after the first EOS cycle, meaning after an EOS, the next non-null depo will re-capture the start time but subsequent processing may behave unexpectedly if there is a second stream of depos.

---

### Fourdee (Fourdee.h / Fourdee.cxx)

**Algorithm overview:** A monolithic application that chains depo source, drifter, optional depo filter, ductor, optional noise source, optional digitizer, optional filter, and sink into a manual pipeline. Has two execution strategies: `execute_new()` (drain-end-first pipeline) and `execute_old()` (simple loop with goto).

**Style issues:** Heavy use of `cerr`, `goto` (line 467), raw `new` for pipeline procs (lines 263-309). The `execute_old()` method uses `goto bail` for error handling.

**Memory leak risk:** Pipeline procs created with `new` (lines 263-309) are stored in a `Pipeline` vector but never explicitly deleted. If `Pipeline` does not own and delete its elements, these leak.

**No critical bugs** beyond the style and memory management concerns. This appears to be legacy/test code.

---

### DuctorFramer (DuctorFramer.h / DuctorFramer.cxx)

**Algorithm overview:** Wraps an IDuctor and IFrameFanin to convert a DepoSet into a single Frame. Feeds each depo (plus an EOS sentinel) to the ductor, collects all output frames, then merges them via the fanin.

**Potential issue (line 62-63):** `all_frames.back()` is checked for non-null to detect missing EOS. If the ductor produces no frames at all (not even an EOS), `all_frames` is empty and `.back()` is undefined behavior on an empty deque/vector.

**Design note:** The fanin receives all EOS markers from intermediate processing as the ductor may produce multiple frames. Stripping only the final EOS (line 65 `pop_back()`) is correct assuming the ductor passes through exactly one EOS.

---

### DepoFanout (DepoFanout.h / DepoFanout.cxx)

**Algorithm overview:** Trivially copies one depo pointer to N output ports. EOS passes through.

**No bugs found.** Simple and correct.

---

### DepoSetFanout (DepoSetFanout.h / DepoSetFanout.cxx)

**Algorithm overview:** Same as DepoFanout but for DepoSet objects. Trivially copies the pointer.

**No bugs found.** Clean implementation.

---

### DepoMerger (DepoMerger.h / DepoMerger.cxx)

**Algorithm overview:** Merge-sorts two streams of time-ordered depos into one output stream. Processes one depo at a time from whichever input has the earlier timestamp. When both inputs have reached EOS, outputs a single EOS.

**Potential issue (lines 40-57):** When `t0 == t1`, both depos are output (both `if` conditions are true). This is intentional for preserving all depos at coincident times.

**Design limitation:** After `m_eos` is set to true (line 79), subsequent calls return false (line 21). The component cannot be reused without external reset, but there is no reset method. This is typical for one-shot pipeline operation.

**Style issue:** Uses `std::cerr` instead of logging framework. No `Aux::Logger` base.

---

### DepoChunker (DepoChunker.h / DepoChunker.cxx)

**Algorithm overview:** Collects depos within a time window and emits a DepoSet when a depo falls outside the window. On EOS, flushes accumulated depos and resets the gate.

**Bug B7 -- Nullptr leaked into depo vector (lines 45-49):**
```cpp
if (!depo) {  // EOS
    emit(deposetqueue);      // emits m_depos, then clears m_depos
    m_depos.push_back(depo); // pushes nullptr into the now-empty m_depos
    m_gate = m_starting_gate;
    return true;
}
```
After `emit()` clears `m_depos`, the nullptr is pushed back in. On the next non-EOS cycle, this nullptr will be included in the accumulation and eventually emitted as part of a DepoSet. Downstream code iterating over the DepoSet's depos may dereference null.

**Edge case (line 60-66):** When a depo arrives past `m_gate.second`, the gate advances by exactly one window width. If a depo arrives more than two windows past the current gate, it will be placed in the immediately-next window but potentially much earlier than expected. There is no loop to advance the gate to the correct window.

---

### DepoBagger (DepoBagger.h / DepoBagger.cxx)

**Algorithm overview:** Collects all depos (optionally within a time gate) until EOS, then emits one DepoSet followed by an EOS.

**Design detail (line 55):** The gate check `m_gate.first == 0.0 and m_gate.second == 0.0` uses floating-point equality to detect the "no gate" sentinel. This works because the default is exactly `{0,0}`, but a user-configured gate of `[0.0, 0.0]` would be indistinguishable from "no gate" and would bag all depos.

**No bugs found.** Straightforward implementation.

---

### DeposOrBust (DeposOrBust.h / DeposOrBust.cxx)

**Algorithm overview:** A hydra node that routes non-empty DepoSets to output port 0 and creates EmptyFrame objects for empty DepoSets on output port 1. Provides a branch pattern.

**Design concern:** The `EmptyFrame` class (lines 19-50) returns `nullptr` from `traces()` and `0.0` from `tick()`. Downstream code that calls `frame->traces()->size()` without null-checking will crash. The `tick()==0` could also cause division-by-zero in time calculations.

**Asymmetric EOS handling:** On EOS, both output ports receive nullptr. But for non-empty depo sets, only port 0 gets data while port 1 gets nothing. For empty depo sets, only port 1 gets data. This asymmetry requires careful downstream synchronization (the header mentions FrameSync as the matching merge).

---

### Retagger (Retagger.h / Retagger.cxx)

**Algorithm overview:** Applies tag rule transformations (frame tags, trace tags, and merge rules) to rewrite tags on a frame without modifying traces. Also supports channel mask map remapping.

**Bug B11 -- Mutation of const reference from input frame (lines 117-121):**
```cpp
auto& summary = get<2>(ttls);  // reference to summary stored in stash
if (summary.empty()) {
    for (size_t ind = 0; ind < traces.size(); ++ind) {
        summary.push_back(0);  // MUTATES the stash entry
    }
}
```
The `summary` reference refers to data originally obtained from `inframe->trace_summary(inttag)` (line 85), stored in the stash tuple. The push_back mutates this stored copy. While this does not corrupt the input frame (it is a copy in the tuple), it does mean that if the same input tag maps to multiple output tags via merge rules, the second iteration sees a non-empty summary that was zero-padded by the first iteration, potentially corrupting the output. The zero-padding should only be accumulated into `osummary`, not into the stash.

---

### Scaler (Scaler.h / Scaler.cxx)

**Algorithm overview:** Scales depo charge based on a YZ position map loaded from a file. Looks up the depo's Y/Z bin in a 2D array and multiplies the charge by the scale factor.

**Bug B4 -- No null/EOS check (line 111-115):**
```cpp
bool Gen::Scaler::operator()(const input_pointer& depo, output_queue& outq)
{
    const double Qi = depo->charge();  // dereferences depo without null check
```
If `depo` is nullptr (EOS marker), this immediately dereferences a null pointer.

**Bug B6 -- Off-by-one in bounds clamping (lines 142-148):**
```cpp
if(depo_bin_y > n_ybin){
    depo_bin_y = n_ybin;       // yzmap[...] is sized n_ybin, valid indices 0..n_ybin-1
}
if(depo_bin_z > n_zbin){
    depo_bin_z = n_zbin;       // yzmap is sized n_zbin, valid indices 0..n_zbin-1
}
```
The clamp should use `>= n_ybin` / `>= n_zbin` and clamp to `n_ybin-1` / `n_zbin-1`. As written, when `depo_bin_y == n_ybin` or `depo_bin_z == n_zbin`, the condition is false, and line 150 accesses `yzmap[n_zbin][n_ybin]` which is out-of-bounds. Even when the condition is true, it clamps to the out-of-bounds index.

**Bug B13 -- Null prior dereference (line 152):**
```cpp
auto newdepo = make_shared<Aux::SimpleDepo>(depo->time(), depo->pos(), Qi*scale,
    depo->prior(), depo->extent_long(), depo->extent_tran(),
    depo->prior()->id(), depo->prior()->pdg(), depo->prior()->energy());
```
`depo->prior()` can be nullptr for root depos. The chained `->id()`, `->pdg()`, `->energy()` calls will crash.

**Efficiency E2 (line 156):** `std::sort(outq.begin(), outq.end(), scaler_by_time)` is called on every single depo insertion, making the overall complexity O(N^2 log N). A final sort after all depos are processed would reduce this to O(N log N).

---

### TruthSmearer (TruthSmearer.h / TruthSmearer.cxx)

**Algorithm overview:** A ductor variant that smears charge depositions in time and wire dimensions without full convolution. Applies Gaussian time smearing and discrete wire weighting to produce "truth" frames.

**Edge case (line 224-229):** `dist_to_wire` comparison uses strict inequalities. When `dist_to_wire == 0.5` exactly (depo at wire boundary), neither `< 0.5` nor `> 0.5` is true, so `wire_weight` remains at the initialized value of 1.0. This means boundary positions get unscaled weight rather than the intended smeared weight.

**Style issue:** Uses `std::cerr` for diagnostics instead of the logging framework.

---

### TruthTraceID (TruthTraceID.h / TruthTraceID.cxx)

**Algorithm overview:** Similar to TruthSmearer but applies frequency-domain wire and time filters. Creates wire filters (HfFilter) per plane and time filters per readout, then convolves charge spectra with these filters.

**Bug B10 -- Wire loop inside depo loop (lines 153-218):**
```cpp
for (auto depo : m_depos) {
    bindiff.add(depo, ...);         // adds depo to binned diffusion
    ...
    for (int iwire = 0; iwire < numwires; iwire++) {
        ...
        for (int imp = min_impact; imp <= max_impact; imp++) {
            auto id = bindiff.impact_data(imp);
            ...
        }
        bindiff.erase(0, min_impact);  // erases data that may be needed for later depos
    }
}
```
The wire loop is nested inside the depo loop. Each depo triggers a full scan of all wires, using `bindiff` that has been progressively erased by previous iterations. This means: (1) for the first depo, only that depo's contribution is processed across all wires; (2) `bindiff.erase()` removes data needed by subsequent depo iterations; (3) the overall result is incorrect -- depos should all be added first, then wires should be scanned once. Compare with TruthSmearer where `bindiff.add()` is in a separate loop from the wire iteration.

**Performance concern (line 175):** `std::cout << "Truth: charge spectrum extracted."` is an unconditional stdout print in the hot loop.

---

### Misconfigure (Misconfigure.h / Misconfigure.cxx)

**Algorithm overview:** Applies a "misconfiguration" to traces by deconvolving the assumed electronics response and reconvolving with a different one. Uses `DftTools::replace()` for the frequency-domain operation.

**Potential issue (line 92):** `wave.resize(charge.size())` always truncates regardless of the `m_truncate` config. The commented-out original code (line 89) would have respected the flag. The current code unconditionally truncates.

**No other bugs found.** Clean, concise implementation.

---

### StaticChannelStatus (StaticChannelStatus.h / StaticChannelStatus.cxx)

**Algorithm overview:** Provides per-channel preamp gain and shaping time via a lookup table. Channels not in the "deviants" map return nominal values.

**No bugs found.** Simple key-value lookup with proper defaults.

---

### ZSEndedTrace (ZSEndedTrace.h / ZSEndedTrace.cxx)

**Algorithm overview:** A trace implementation that accumulates charge by time bin and provides end-zero-suppressed output. Internally uses a `std::map<int, float>` for sparse storage.

**Bug B12 -- Undefined behavior when m_nbins==0 (line 14):**
```cpp
bin = min(bin, m_nbins - 1);  // m_nbins is int, m_nbins-1 = -1
```
When `m_nbins` is 0 (default), `m_nbins - 1` is -1. The `min()` with a potentially positive `bin` returns -1, which is then used as a map key. This is not UB per se (map keys can be negative), but it means the trace has a negative tbin, which may confuse downstream code. More critically, if the intention was to clamp to valid bins, this defeats the purpose.

**Thread safety:** The `charge()` method lazily populates `m_charge` (marked `mutable`). If called concurrently, this is a data race.

---

### TimeGatedDepos (TimeGatedDepos.h / TimeGatedDepos.cxx)

**Algorithm overview:** Filters depos by time gate. On EOS, advances the gate by a configurable period. Supports accept or reject mode.

**Bug in configure (line 58):**
```cpp
m_accept = cfg["accept"].asString() == "accept";
```
The default_configuration sets `cfg["accept"] = m_accept` (a boolean, default true). But configure() compares `asString()` against "accept". A boolean true converts to the string "true" via `asString()`, not "accept". So the default config would yield `m_accept = false` after round-tripping through configure, which is the opposite of the declared default.

**Edge case:** If `m_period` is 0 (default), the gate never advances on EOS. This is documented behavior but may surprise users who expect auto-advance.

---

### DepoSetFilter (DepoSetFilter.h / DepoSetFilter.cxx)

**Algorithm overview:** Filters depos in a DepoSet by checking if they fall within any anode face's sensitive bounding box.

**Bug B9 -- Wrong sign in log message (line 66):**
```cpp
log->debug("call={} depos: input={} output={} dropped={}",
           m_count, nin, nout, nout-nin);  // should be nin-nout
```
`nout - nin` is negative (or zero) since filtering can only reduce depo count. Should be `nin - nout`.

**No other bugs found.** Correct filtering logic with proper EOS handling.

---

### DepoSetFilterYZ (DepoSetFilterYZ.h / DepoSetFilterYZ.cxx)

**Algorithm overview:** Filters depos based on a YZ response map loaded from file. Keeps depos whose YZ bin maps to a specific response value.

**Bug B5 -- Double unit application (line 111-112):**
```cpp
double depo_y = idepo->pos().y()*units::mm;
double depo_z = idepo->pos().z()*units::mm;
```
Positions returned by `idepo->pos()` are already in the WCT internal unit system (which is mm-based). Multiplying by `units::mm` applies the unit conversion a second time, making the coordinates 1000x too large (since `units::mm` is approximately 1.0 in WCT, this may be a near-no-op depending on the unit system, but it is conceptually wrong and would break under non-default unit systems).

**Bug B6 -- Same off-by-one as Scaler (lines 120-130):**
```cpp
if(depo_bin_y > nbinsy)
    depo_bin_y = nbinsy;       // out-of-bounds index
if(depo_bin_z > nbinsz)
    depo_bin_z = nbinsz;       // out-of-bounds index
```
Identical to the Scaler bug. Should clamp to `nbinsy-1` and `nbinsz-1`, and the condition should use `>=`.

**Default config issue (lines 25-37):** Default configuration values are strings like `"BinWidth"`, `"TPCWidth"` etc. instead of numeric defaults. If a user does not provide these, `get<double>` will fail or return 0, leading to division by zero in bin calculations.

---

### DepoSetRotate (DepoSetRotate.h / DepoSetRotate.cxx)

**Algorithm overview:** Applies coordinate transpose and scaling to depo positions within a DepoSet. Supports axis swapping and sign flipping to simulate rotations.

**Bug B8 -- Depos dropped when rotate is false (lines 56-66):**
```cpp
for (auto depo: *(in->depos())) {
    if (m_rotate) {
        ...
        all_depos.push_back(newdepo);
    }
    // no else clause!
}
```
When `m_rotate` is false, no depos are added to `all_depos`. The output DepoSet will be empty. The component should pass through depos unchanged when rotation is disabled. This means configuring `rotate: false` silently drops all depos.

---

### DepoSetScaler (DepoSetScaler.h / DepoSetScaler.cxx)

**Algorithm overview:** Wraps a per-depo IDrifter (Scaler) to operate on DepoSets. Feeds each depo to the scaler and collects results.

**Potential issue (line 48-49):** The EOS flush (`in_depos.push_back(nullptr)`) is commented out. If the underlying scaler buffers state that needs an EOS to flush (as many IDrifter implementations do), some depos may be lost. However, the Scaler class processes depos immediately, so this is not currently a problem.

**Debug residue:** Lines 52, 73 show `charge_in` and `charge_out` initialized to 0 but never accumulated (the accumulation code is commented out), yet they are logged. The log always shows `in=0 out=0`.

---

### DumpDepos (DumpDepos.h / DumpDepos.cxx)

**Algorithm overview:** Debug sink that prints depo information to stderr.

**Style note:** In `WireCell` namespace directly, not `WireCell::Gen`. Uses `cerr` directly.

**No bugs found.** Trivial implementation.

---

### DumpFrames (DumpFrames.h / DumpFrames.cxx)

**Algorithm overview:** Debug sink that logs frame metadata (ident, time, trace count, tags).

**No bugs found.** Clean implementation using the logging framework.

---

## Cross-Cutting Concerns

### Inconsistent logging
Many older components (MultiDuctor, Fourdee, DepoMerger, DumpDepos, TruthSmearer, TruthTraceID) use `std::cerr` instead of the `Aux::Logger` framework. This makes it impossible to control log levels at runtime.

### Raw new usage
FrameFanin (line 196) and Retagger (line 45) use raw `new` for SimpleFrame. While wrapped in shared_ptr shortly after, exception safety would benefit from `make_shared`.

### EOS handling inconsistency
Some components (DepoBagger) pass EOS through, others (DepoChunker) reset state but accidentally leak nullptr into depo vectors, and others (DepoSetFilter) return `true` with `out=nullptr` on EOS without explicitly pushing an EOS marker.

### Repeated off-by-one pattern
The bounds-clamping pattern in Scaler and DepoSetFilterYZ (`if (idx > N) idx = N`) is an identical bug in two places. The correct pattern is `if (idx >= N) idx = N-1`.
