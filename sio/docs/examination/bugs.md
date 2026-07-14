# SIO Package: Potential Bugs

## 1. DepoFileSource: Wrong variable in column check (HIGH)

**File:** `src/DepoFileSource.cxx:178`

The second column check re-tests `darr.cols()` instead of `iarr.cols()`:

```cpp
const size_t ndata = 7;
if (darr.cols() != ndata) {   // line 171 - checks data columns, correct
    ...
}

const size_t ninfo = 4;
if (darr.cols() != ndata) {   // line 178 - BUG: should be "iarr.cols() != ninfo"
    ...
}
```

The info array column count is never actually validated. A malformed info array
with the wrong number of columns would pass this check and lead to
out-of-bounds access at lines 199, 210-211, 218.

---

## 2. DepoFileSource: Wrong variable in error log message (LOW)

**File:** `src/DepoFileSource.cxx:116`

The ident mismatch error message logs `m_count` instead of `fmd1.ident`:

```cpp
log->error("call={}, ident mismatch {} != {} in file={}",
           m_count, m_count, fmd1.ident, fmd2.ident);
//                  ^^^^^^^ should be fmd1.ident
```

This means the error message would print the call count twice instead of
showing both mismatched idents.

---

## 3. TensorFileSink: Wrong class name in error message (LOW)

**File:** `src/TensorFileSink.cxx:38`

```cpp
const std::string msg = "ClusterFileSink: unsupported outname: " + m_outname;
//                       ^^^^^^^^^^^^^^^^ should be "TensorFileSink"
```

Copy-paste error. The error message refers to ClusterFileSink, but the code
is in TensorFileSink.

---

## 4. TensorFileSink: m_count not incremented in dump_mode (LOW-MEDIUM)

**File:** `src/TensorFileSink.cxx:87-91`

When `m_dump_mode` is true, the function returns early without incrementing
`m_count`. This is inconsistent - the debug log at line 83 increments with
`m_count++`, but that only happens on EOS. During dump mode:

```cpp
if(m_dump_mode) {
    log->debug("dumping tensor set ident={} at call {}",
               in->ident(), m_count);
    return true;   // m_count never incremented
}
```

The `m_count` will always show 0 for non-EOS calls in dump mode, making the
log confusing.

---

## 5. NumpyDepoTools: Raw pointer management risk (MEDIUM)

**File:** `src/NumpyDepoTools.cxx:59-84`

Raw `new` is used to create `SimpleDepo` objects, and they are stored as raw
pointers in a `std::vector<Aux::SimpleDepo*>`:

```cpp
auto sdepo = new Aux::SimpleDepo(...);
```

If an early `return false` is hit (e.g., at line 76), all previously
allocated `SimpleDepo` objects that have not been wrapped in a shared pointer
are leaked. The gen=0 depos remain as raw pointers in `sdepos` and are never
freed.

Compare with `DepoFileSource.cxx:190` which correctly uses
`std::make_shared<Aux::SimpleDepo>(...)`.

---

## 6. NumpyDepoTools: Prior depo linking overwrites without check (LOW)

**File:** `src/NumpyDepoTools.cxx:81`

```cpp
sdepos[other]->set_prior(idepo);
```

If `sdepos[other]` has already had its prior set (perhaps from a different
prior depo referencing the same child), this silently overwrites the previous
prior link. No check or warning is given. This is the same behavior in
`DepoFileSource.cxx:224`, so it may be intentional, but could mask data
issues.

---

## 7. JsonDepoSource: Uses raw `new`/`delete` for adapter (LOW-MEDIUM)

**File:** `src/JsonDepoSource.cxx:77, 120-146`

The `m_adapter` is managed with raw `new`/`delete`:

```cpp
m_adapter = new ElectronsAdapter(scale);
...
if (m_adapter) {
    delete m_adapter;
    m_adapter = nullptr;
}
```

If `configure()` throws after allocation but before completion, or if the
destructor doesn't call `delete m_adapter`, there is a memory leak. Should
use `std::unique_ptr<JsonRecombinationAdaptor>`.

Additionally, the destructor `~JsonDepoSource() {}` does NOT delete
`m_adapter`, so there IS a memory leak whenever the component is destroyed
without a second call to `configure()`.

---

## 8. JsonDepoSource: Uses `cerr` instead of logging framework (LOW)

**File:** `src/JsonDepoSource.cxx:133, 139, 166`

Multiple places use `cerr` for output instead of the toolkit's logging
framework. This component does not inherit from `Aux::Logger` (unlike most
other sio components), so it has no structured logging:

```cpp
cerr << "Sio::JsonDepoSource: using electrons with scale=" << scale << endl;
cerr << "Sio::JsonDepoSource::configure: unknown recombination model: ..."
cerr << "Sio::JsonDepoSource::configure: slurped in ..."
```

---

## 9. JsonDepoSource: Silent failure on missing recombination model (MEDIUM)

**File:** `src/JsonDepoSource.cxx:136-146`

If the model type is neither "electrons", "MipRecombination",
"BirksRecombination", nor "BoxRecombination", `m_adapter` remains `nullptr`.
Later in `jdepo2idepo()`, dereferencing `(*m_adapter)(jdepo)` will
segfault:

```cpp
const double q = (*m_adapter)(jdepo);  // nullptr dereference if model unknown
```

Even if the model is found but is an unrecognized type, the `if` chain at
lines 141-145 falls through without setting `m_adapter`.

---

## 10. JsonDepoSource: Duplicate include guard collision with BeeDepoSource (MEDIUM)

**File:** `inc/WireCellSio/JsonDepoSource.h:19` and `inc/WireCellSio/BeeDepoSource.h:29`

Both headers use the same include guard macro:

```cpp
#ifndef WIRECELLSIO_JSONDEPOSOURCE   // in BeeDepoSource.h
#define WIRECELLSIO_JSONDEPOSOURCE   // in BeeDepoSource.h
```

```cpp
#ifndef WIRECELLSIO_JSONDEPOSOURCE   // in JsonDepoSource.h
#define WIRECELLSIO_JSONDEPOSOURCE   // in JsonDepoSource.h
```

If both headers are included in the same translation unit, the second one
will be completely skipped. `BeeDepoSource.h` should use
`WIRECELLSIO_BEEDEPOSOURCE` instead.

---

## 11. NumpyFrameSaver: Saves beyond unique channel range (LOW)

**File:** `src/NumpyFrameSaver.cxx:148`

```cpp
cnpy::npz_save(fname, aname, channels.data(), {nrows}, mode);
```

`nrows` is computed from `std::distance(chbeg, chend)` where `chend` is the
result of `std::unique`. However, `channels.data()` still points to the full
original vector. The save uses `{nrows}` as the shape, so only `nrows`
elements are written, which is correct. But the underlying vector
`channels` still has its original full size, which is slightly misleading
(not a data bug, just a clarity issue).

---

## 12. FrameFileSource: No validation that tickinfo has >= 3 elements (MEDIUM)

**File:** `src/FrameFileSource.cxx:324, 337-338`

The tickinfo vector is accessed by index without bounds checking:

```cpp
const int tbin0 = (int)framelet.tickinfo[2];  // line 324
...
const double time = framelets[0].tickinfo[0];  // line 337
const double tick = framelets[0].tickinfo[1];  // line 338
```

If the tickinfo array in the file is malformed or has fewer than 3 elements,
this will access beyond vector bounds, causing undefined behavior.

---

## 13. ClusterFileSource::load_numpy: Dangling reference to pigenc data (HIGH)

**File:** `src/ClusterFileSource.cxx:180-197`

```cpp
pigenc::File pig;
pig.read(m_in);
...
const node_element_t* data = pig.as_type<node_element_t>();
nas.emplace(pf.code[0], boost::const_multi_array_ref<node_element_t, 2>(data, shape));
```

The `pig` object is local to the while loop body. When the loop iterates,
`pig` is destroyed, freeing its data buffer. But `nas` and `eas` store
`const_multi_array_ref` objects that hold raw pointers to `pig`'s data. After
the loop, these references point to freed memory. The call to
`to_cluster(nas, eas, m_anodes)` at line 235 then reads from dangling
pointers.

This is a use-after-free bug, though it may work in practice if `to_cluster`
copies the data before the ref becomes invalid within each iteration (needs
verification of `to_cluster` behavior).

---

## 14. FrameFileSink: Summary written in original trace order, not sorted-channel order (HIGH) — FIXED

**File:** `src/FrameFileSink.cxx`, `one_tag()` function

When serializing a tagged frame, `FrameFileSink` builds the 2D waveform array in
sorted-channel order (line 157: `std::sort(tmp.begin(), tmp.end())`), but writes
the trace `summary` (e.g. per-wire Wiener thresholds) in the original trace-iteration
order (osp-wire order in `OmnibusSigProc`):

```cpp
// channels written in sorted order:
std::sort(tmp.begin(), tmp.end());
channels.insert(channels.begin(), chbeg, chend);

// summary written in original (unsorted) order:
write(m_out, aname, summary);   // BUG: not reordered to match channels
```

When `FrameFileSource` reads the file back, it pairs `summary[i]` with the
i-th sorted channel, creating a systematic mismatch.  In `MagnifySink` this
caused the per-channel Wiener-threshold 1D histograms (`h{u,v,w}_threshold<N>`)
to show threshold = 0 on channels that are *not* dead, while the actually-dead
channels (in `T_bad`) received nonzero thresholds from adjacent good channels.
The mismatch was exact: per-plane counts of threshold-0 channels equalled the
dead-channel count, but the channel identities were wrong.

Affected any detector whose WAN (Wire-Attachment-Number) ordering is not
identical to WCT channel-ID ordering, which is the case for ProtoDUNE-VD and
ProtoDUNE-HD.

**Fix (committed):** Before writing the summary, build a `chid→summary_value`
map from the traces and re-emit summary values in the same order as the sorted
`channels` array.  This keeps the positional correspondence that
`FrameFileSource` relies on.

---

## 15. TensorFileSink header: Swapped metadata/array file extensions in comment (LOW)

**File:** `inc/WireCellSio/TensorFileSink.h:72-73`

The comment says:
```
<prefix>tensor_<ident>_<index>_metadata.npy   <- says .npy
<prefix>tensor_<ident>_<index>_array.json     <- says .json
```

But the actual code at `src/TensorFileSink.cxx:102-103` writes:
```cpp
jsonify(ten->metadata(), ppre + "_metadata.json");  // metadata is .json
numpyify(ten, ppre + "_array.npy");                 // array is .npy
```

The comment has the extensions swapped.
