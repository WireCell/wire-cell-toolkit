# General I/O Components Examination

> **Bug Fix Status**: Bugs 1-3, 5-9, 11-14 fixed. Bug 4 (MagnifySource channel numbering) deferred — may be intentional for Magnify format. Bug 10 (MagnifySink do_shunt) deferred — complex ownership semantics need further review.
>
> **Efficiency Fix Status**: Efficiency 1 (CelltreeSource file reopening) fixed — file opened once, TTree* passed to read_traces/read_cmm. Efficiency 2 (MagnifySink shared_ptr copy) fixed — uses const auto&. Efficiency 3 (HistFrameSink Fill→SetBinContent) fixed — uses direct bin index computation. Efficiencies 4-7 deferred — low impact.

Files examined:
- `root/src/MagnifySource.cxx` / `root/inc/WireCellRoot/MagnifySource.h`
- `root/src/MagnifySink.cxx` / `root/inc/WireCellRoot/MagnifySink.h`
- `root/src/HistFrameSink.cxx` / `root/inc/WireCellRoot/HistFrameSink.h`
- `root/src/CelltreeSource.cxx` / `root/inc/WireCellRoot/CelltreeSource.h`
- `root/src/CelltreeFrameSink.cxx` / `root/inc/WireCellRoot/CelltreeFrameSink.h`
- `root/src/RootfileCreation.cxx` / `root/inc/WireCellRoot/RootfileCreation.h`

---

## 1. Potential Bugs

### Bug 1: `RootfileCreation_frames` registered with wrong interface type
- **File:** `RootfileCreation.cxx`, line 10
- **Severity:** HIGH
- **Description:** The `WIRECELL_FACTORY` macro for `RootfileCreation_frames` registers it as `WireCell::IDepoFilter` instead of `WireCell::IFrameFilter`. The class in the header (`RootfileCreation.h`, line 31) inherits from `IFrameFilter`, and its `operator()` (line 85) takes `IFrame::pointer` arguments. The factory registration is inconsistent with the actual type, which means the framework will not be able to look up this component as an `IFrameFilter`. It may silently fail to be found or cause a runtime type mismatch.

### Bug 2: `MagnifySource` does not check for null TFile
- **File:** `MagnifySource.cxx`, line 66
- **Severity:** HIGH
- **Description:** After `TFile::Open(url.c_str())`, there is no check for whether `tfile` is null or a zombie file. `TFile::Open` returns null if the file cannot be opened. The code immediately proceeds to `tfile->Get("Trun")` on line 72, which will crash with a null pointer dereference if the file does not exist or is inaccessible.

### Bug 3: `MagnifySource` does not check for null histogram
- **File:** `MagnifySource.cxx`, line 126
- **Severity:** HIGH
- **Description:** `tfile->Get(hist_name.c_str())` is cast to `TH2*` without a null check. If the histogram does not exist in the file, the code will crash on line 130 (`hist->GetNbinsX()`). The text comments at line 37-41 say "missing histograms raise an exception" but no exception is actually thrown; the code just dereferences null.

### Bug 4: `MagnifySource` uses wrong channel numbering
- **File:** `MagnifySource.cxx`, lines 121, 146-148
- **Severity:** MEDIUM
- **Description:** The variable `channel_number` is a simple sequential counter starting from 0 that increments across all planes and all frame tags. This means the channel number assigned to each trace is just a positional index (0, 1, 2, ...), not the actual detector channel number. This may work for some Magnify files where channels are stored in order, but is fragile. The channel number is not read from the file. If different frame tags are loaded, the counter continues across them, so the second tag's traces start with the channel count from where the first tag ended, which is almost certainly wrong.

### Bug 5: `CelltreeSource::read_traces` and `read_cmm` leak memory
- **File:** `CelltreeSource.cxx`, lines 79-81 (`read_traces`), lines 152-154 (`read_cmm`)
- **Severity:** MEDIUM
- **Description:** In `read_traces`, `channelid`, `channel_threshold`, and `esignal` are allocated with `new` but never deleted. Similarly in `read_cmm`, `badChannel`, `badBegin`, and `badEnd` are allocated with `new` but never deleted. These are leaked every time the function is called. Also, `tfile` in `read_traces` is closed but not deleted (line 140), leaking the `TFile` object.

### Bug 6: `CelltreeSource::operator()` leaks TFile
- **File:** `CelltreeSource.cxx`, line 200
- **Severity:** LOW
- **Description:** `tfile` is opened and closed (line 231) but never deleted. ROOT `TFile::Open` returns a heap-allocated object that should be deleted after use. Additionally, `read_traces` opens the same file again independently for each branch, so the file is opened N+1 times total (once in `operator()`, once per branch in `read_traces`).

### Bug 7: `CelltreeFrameSink` uses channel ID as vector index for thresholds
- **File:** `CelltreeFrameSink.cxx`, lines 226-235
- **Severity:** HIGH
- **Description:** The `channelThreshold` vector is resized to `ntot` (number of traces), but then indexed by `chid` (channel ID) at line 232: `channelThreshold->at(chid) = val`. For typical detectors, channel IDs can be in the thousands (e.g., 0-8255 for MicroBooNE), while `ntot` is just the trace count. This will cause an out-of-range exception (or undefined behavior with `operator[]`) whenever `chid >= ntot`. This is a data-corrupting bug in all realistic scenarios.

### Bug 8: `CelltreeFrameSink` skips multi-period masks but pushes channel
- **File:** `CelltreeFrameSink.cxx`, lines 271-281
- **Severity:** MEDIUM
- **Description:** When a channel has more than one dead period (`mask.size() != 1`), the code prints a warning and `continue`s (line 276). However, the channel ID has already been pushed into the `Channel` vector (line 271) before the size check. This means `Channel` will have an entry without corresponding `Begin`/`End` entries, causing the three vectors to become misaligned. Downstream consumers that read these vectors in parallel will associate wrong values.

### Bug 9: `CelltreeFrameSink` leaks heap-allocated vectors
- **File:** `CelltreeFrameSink.cxx`, lines 150, 158, 222, 260-262
- **Severity:** MEDIUM
- **Description:** Multiple objects are allocated with `new` and never deleted: `raw_channelId` (line 150), `sim_wf` (line 158), `channelThreshold` (line 222), `Channel`, `Begin`, `End` (lines 260-262). These leak on every call.

### Bug 10: `MagnifySink::do_shunt` deletes TFile that owns objects being moved
- **File:** `MagnifySink.cxx`, line 304
- **Severity:** MEDIUM
- **Description:** After shunting histograms from `input_tf` to `output_tf` via `hist->SetDirectory(output_tf)`, `input_tf` is deleted (line 304). For histograms, `SetDirectory` properly transfers ownership, so this is correct. However, for trees, `tree->CloneTree()` creates a clone whose directory is then set to `output_tf` (lines 289-290). The original `obj` from `Get()` is not the cloned tree, but the code does `delete input_tf` which will also try to delete the original tree. This is fine for the clone, but if the `dynamic_cast<TTree*>` fails and falls through to `dynamic_cast<TH1*>` which also fails, the fetched `obj` pointer is leaked (line 280 fetches it, but nothing handles the "unknown type" case except a warning on line 301).

### Bug 11: `MagnifySink` summary histogram has wrong bin range
- **File:** `MagnifySink.cxx`, line 424
- **Severity:** LOW
- **Description:** The TH1F for summaries is created with `Binning(chf - ch0 + 1, ch0, chf)`. The number of bins is `chf - ch0 + 1` but the range is `[ch0, chf]`, so the bin width is `(chf - ch0) / (chf - ch0 + 1)` which is less than 1. For integer channel numbers, bins should span `[ch0 - 0.5, chf + 0.5]` or `[ch0, chf + 1]` to get unit-width bins. This means `FindBin` on line 429 may map channels to wrong bins.

### Bug 12: `MagnifySink` ignores configured `root_file_mode` during `create_file`
- **File:** `MagnifySink.cxx`, line 166
- **Severity:** LOW
- **Description:** The `create_file()` method hardcodes `"RECREATE"` mode (line 166) rather than using `m_cfg["root_file_mode"]`. So regardless of the user's configuration, the file is always recreated during `configure()`. Then in `operator()`, the configured mode (line 322) is used. If user sets mode to `"UPDATE"`, the file is still recreated at configure time, which would destroy any pre-existing content.

### Bug 13: `MagnifySource` does not close or delete the TFile
- **File:** `MagnifySource.cxx`, line 66
- **Severity:** LOW
- **Description:** The `TFile` opened on line 66 is never closed or deleted. This means the file handle remains open until program exit. While ROOT will eventually clean up global file objects, this is a resource leak.

### Bug 14: Jsonnet config has misspelled tag
- **File:** `root/test/test-celltree-to-framefile.jsonnet`, line 14
- **Severity:** LOW
- **Description:** The `frame_tags` list contains `"weiner"` (line 14) which is likely a misspelling of `"wiener"`. However, `trace_tags` on line 8 correctly uses `"wiener"`. This mismatch may cause the wiener data to be tagged incorrectly.

---

## 2. Efficiency Improvements

### Efficiency 1: `CelltreeSource::read_traces` reopens file per branch
- **File:** `CelltreeSource.cxx`, lines 70-71 (called from line 245 in a loop)
- **Current approach:** Each call to `read_traces` opens the ROOT file, gets the tree, reads the entry, and closes the file. The main `operator()` also opens the same file separately.
- **Suggested improvement:** Open the file once in `operator()`, pass the `TTree*` pointer to `read_traces`, and close after all branches are read.
- **Expected impact:** MEDIUM -- eliminates redundant file open/close cycles, each of which involves filesystem I/O and ROOT metadata parsing.

### Efficiency 2: `MagnifySink::collate_byplane` copies traces by value
- **File:** `MagnifySink.cxx`, line 118
- **Current approach:** The function signature takes `const ITrace::vector& traces` and `ITrace::vector byplane[]`. Internally, it copies `ITrace::pointer` (shared_ptr) into the per-plane arrays. The outer loop on line 122 also iterates `for (auto trace : traces)` which copies each shared_ptr.
- **Suggested improvement:** Use `for (const auto& trace : traces)` to avoid incrementing/decrementing reference counts.
- **Expected impact:** LOW -- shared_ptr copy involves an atomic increment/decrement, which is cheap but unnecessary for many traces.

### Efficiency 3: `HistFrameSink` uses `Fill()` instead of `SetBinContent()`
- **File:** `HistFrameSink.cxx`, lines 119-122
- **Current approach:** Each sample is placed using `hist->Fill(t, fch, charge)` which involves ROOT's bin-finding algorithm (binary search) for every single sample.
- **Suggested improvement:** Compute bin indices directly and use `SetBinContent()`, similar to what `MagnifySink` does. The MagnifySink code comments (lines 369-374) explicitly note that `Fill()` produces files that are 2x larger than `SetBinContent()` due to Sumw2 error storage.
- **Expected impact:** MEDIUM -- avoids per-sample bin lookup overhead and prevents ROOT from automatically creating the Sumw2 error array, significantly reducing both CPU time and output file size.

### Efficiency 4: `CelltreeFrameSink` string comparison chain for tag mapping
- **File:** `CelltreeFrameSink.cxx`, lines 152-155 and 161-163
- **Current approach:** A chain of `if (!tag.compare("gauss"))` / `if (!tag.compare("wiener"))` / `if (!tag.compare("orig"))` is used to map tag names to branch names. If none match, `channelIdname` and `wfname` remain empty, which will create branches with empty names.
- **Suggested improvement:** Use a `std::unordered_map<std::string, std::string>` lookup or, better, make the branch naming configurable rather than hardcoded. Also add an `else` clause or default to handle unknown tags.
- **Expected impact:** LOW -- the current approach works for the three known tags but is not extensible.

### Efficiency 5: `MagnifySink::operator()` iterates traces multiple times
- **File:** `MagnifySink.cxx`, lines 326-384
- **Current approach:** For each frame tag, `tagged_traces()` retrieves the trace list, then `collate_byplane()` iterates all traces to find min/max channel/tbin and sort by plane, then the fill loop iterates traces again.
- **Suggested improvement:** Combine the collation and histogram filling into a single pass.
- **Expected impact:** LOW -- the trace lists are typically not huge, and the operations are simple.

### Efficiency 6: `MagnifySource` uses `std::cerr` instead of the logging framework
- **File:** `MagnifySource.cxx`, lines 54, 58, 61, 84, 92, 97, 125, 150-151, 159
- **Current approach:** All diagnostic output goes to `std::cerr` with no log level control. This produces output unconditionally.
- **Suggested improvement:** Use `WireCellUtil/Logging.h` logger (as `MagnifySink` and `CelltreeSource` do) to allow log-level filtering.
- **Expected impact:** LOW -- cosmetic but affects usability in production where stderr noise is undesirable.

### Efficiency 7: `CelltreeFrameSink` and `HistFrameSink` also use `std::cerr`
- **File:** `CelltreeFrameSink.cxx`, lines 113, 119, 148, 196, 219, 258, 274, 290, 292; `HistFrameSink.cxx`, lines 43, 46, 55, 81, 126
- **Current approach:** Same as above -- raw `std::cerr` instead of structured logging.
- **Suggested improvement:** Use the `Aux::Logger` mixin (as `CelltreeSource` does) or `Log::logger()`.
- **Expected impact:** LOW

---

## 3. Algorithm and Code Explanation

### 3.1 MagnifySource (`MagnifySource.cxx`)

**Purpose:** Reads a "Magnify" format ROOT file (used by the Magnify event display) and produces a single `IFrame` containing the waveform data. The source can be called at most twice: the first call returns the frame, the second returns EOS (end of stream).

**Key data structures:**
- `m_calls`: tracks call count for the two-shot source/EOS pattern
- `m_cfg`: JSON configuration holding filename, frame tags, and CMM tree mappings

**Algorithm walkthrough:**
1. On first call (`m_calls == 1`), opens the ROOT file and reads the `Trun` tree for event metadata (event number, number of ticks).
2. Reads channel mask maps (CMM) from configured TTrees. Each entry contains a channel ID, plane, start time, and end time defining masked regions.
3. For each configured frame tag (e.g., "raw"), iterates over three planes (u, v, w). For each plane, loads a 2D histogram named `h{u,v,w}_{tag}` where X bins are channels and Y bins are ticks.
4. Extracts the charge values from each histogram row (channel) and creates `SimpleTrace` objects. Traces are tagged with the frame tag and tracked by index.
5. Constructs a `SimpleFrame` with a hardcoded tick of 0.5 microseconds and the loaded CMM.

**Edge cases:**
- Missing `Trun` tree is tolerated (warning printed, defaults used)
- Second call returns EOS (null frame, `true` return)
- Third+ calls return `false` (past EOS)

### 3.2 MagnifySink (`MagnifySink.cxx`)

**Purpose:** Writes frame data to a Magnify-format ROOT file. Despite the name "Sink", it implements `IFrameFilter` so it passes the frame through, allowing it to sit in the middle of a pipeline chain.

**Key data structures:**
- `m_anode`: anode plane for channel-to-wire-plane resolution
- `m_nrebin`: rebinning factor to reduce time resolution in output
- `m_cfg`: full configuration including file names, tags, shunt list, and run info

**Algorithm walkthrough:**
1. On `configure()`, creates an empty output ROOT file (always with RECREATE mode).
2. On each `operator()` call:
   - Opens the output file with the configured mode (RECREATE or UPDATE).
   - For each configured frame tag, retrieves tagged traces, collates them by plane using the anode's channel resolution.
   - `collate_byplane()` determines channel and tick ranges for each plane, producing `Binning` objects for both channel (centered on integer IDs) and tick axes.
   - Creates a TH2F per plane with appropriate binning, fills it using `SetBinContent` (accumulating values for rebinned ticks).
   - Handles trace summaries: for each summary tag, creates a per-plane TH1F and fills it with summary values using either "set" or "sum" aggregation.
   - Writes channel mask maps to TTrees.
   - Calls `do_shunt()` to optionally copy objects (trees, histograms) from an input Magnify file to the output file. Can also construct a `Trun` tree from JSON runinfo or from a celltree input file.
   - Writes and closes the output file.

**Important details:**
- The "shunt" mechanism violates DFP encapsulation by reading an input file outside the data flow graph.
- The `create_file()` call during configure ensures the file exists before any operator calls, which is needed if multiple sinks share the file with UPDATE mode.

### 3.3 HistFrameSink (`HistFrameSink.cxx`)

**Purpose:** Saves each frame as a set of per-plane TH2F histograms in a ROOT file. The filename can contain a printf format specifier (%) which is formatted with the frame ID, allowing multiple frames to be written to separate files.

**Key data structures:**
- `m_anode`: anode plane for channel-to-plane resolution
- `m_filepat`: filename pattern with optional format specifier
- `m_units`: charge unit scaling factor

**Algorithm walkthrough:**
1. For each incoming frame, formats the filename using the frame's ident.
2. Collates traces by wire plane identity, tracking channel ranges and tick ranges per plane.
3. For each plane, creates a TH2F with time (in microseconds) on X and channel on Y.
4. Fills the histogram using `Fill()` with time computed from the frame's base time plus tick offsets.
5. Writes and closes the file.

**Important details:**
- Unlike MagnifySink, this sink uses physical time (microseconds) on the X axis rather than tick indices.
- Charge values are divided by `m_units` for output.
- Uses `wpid.ident()` for plane grouping rather than `wpid.index()`, so plane identifiers may not be 0,1,2.

### 3.4 CelltreeSource (`CelltreeSource.cxx`)

**Purpose:** Reads a "celltree" format ROOT file (produced by LArSoft Celltree module) and produces a single `IFrame`. Like MagnifySource, it uses the two-shot source/EOS pattern.

**Key data structures:**
- `m_calls`: call counter for source/EOS pattern
- `m_cfg`: configuration with filename, entry number, branch names, output tags, threshold branch names, and time scale

**Algorithm walkthrough:**
1. On first call, opens the file and reads event metadata (run/subrun/event numbers) from the `/Event/Sim` tree.
2. For each configured input branch, calls `read_traces()` which:
   - Opens the file, gets the `/Event/Sim` tree
   - Reads the channel ID vector and waveform TClonesArray for the specified entry
   - Optionally reads a threshold vector
   - Expands the waveform by `time_scale` factor (e.g., each celltree tick becomes 4 output ticks, with the value divided by time_scale to conserve total charge)
   - Creates `SimpleTrace` objects for each channel
3. Calls `read_cmm()` to load bad channel information (channel ID, begin/end tick ranges).
4. Constructs a `SimpleFrame` with tagged traces and threshold summaries.

**Important details:**
- The `time_scale` parameter handles the difference between celltree's compressed time binning and the full tick resolution needed downstream.
- Threshold values are also divided by `time_scale` to maintain consistency.

### 3.5 CelltreeFrameSink (`CelltreeFrameSink.cxx`)

**Purpose:** Writes frame data to a celltree-format ROOT file. Like MagnifySink, it implements `IFrameFilter` to pass frames through.

**Key data structures:**
- `m_anode`: anode plane for channel resolution
- `m_nsamples`: total number of time samples (must be configured)
- `m_nrebin`: rebinning factor
- `m_cfg`: configuration for output file, tags, summaries

**Algorithm walkthrough:**
1. Opens/creates the output ROOT file with an `Event/Sim` TTree structure.
2. If the Sim tree does not exist, creates it with dummy run/subrun/event branches.
3. For each configured frame tag:
   - Maps tags ("gauss", "wiener", "orig") to celltree branch naming conventions (e.g., "gauss" -> "calibGaussian_channelId", "calibGaussian_wf").
   - Creates a channelId branch and a TClonesArray of TH1F for waveforms.
   - For each trace, creates a TH1F with `nsamples/nrebin` bins and fills it from the trace's charge sequence, supporting rebinning.
   - Fills the branches.
4. Handles summaries (thresholds) -- but has the critical indexing bug described in Bug 7.
5. Writes channel mask (bad channel) information to vector branches.

**Important details:**
- The tag-to-branch-name mapping is hardcoded for exactly three tags. Any other tag produces empty branch names.
- The `cget_tagged_traces()` helper function first checks for tagged traces, then falls back to frame-tagged traces (returning all traces).

### 3.6 RootfileCreation (`RootfileCreation.cxx`)

**Purpose:** A minimal pass-through filter whose sole job is to create an empty ROOT file during configuration. This is used to ensure the output file exists before downstream components (like MagnifySink with UPDATE mode) try to write to it.

**Key data structures:**
- `m_cfg`: configuration with output filename and ROOT file mode

**Algorithm walkthrough:**
1. On `configure()`, opens the ROOT file with the configured mode, immediately closes it, and deletes the TFile object.
2. On each `operator()` call, simply passes the input through to the output unchanged.

**Important details:**
- Two variants exist: `RootfileCreation_depos` (for depo pipelines) and `RootfileCreation_frames` (for frame pipelines).
- The `_frames` variant has the factory registration bug (Bug 1) where it is registered as `IDepoFilter` instead of `IFrameFilter`.
- The `Close("R")` call uses the "R" option which means "also delete the file from ROOT's internal list", ensuring clean separation.
