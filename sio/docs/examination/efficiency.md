# SIO Package: Efficiency Observations

## 1. NumpyDepoTools: Unnecessary double-pass and raw pointer overhead

**File:** `src/NumpyDepoTools.cxx:49-93`

The function uses raw `new` for depo creation then wraps in shared pointers
in a second pass. The `DepoFileSource` counterpart uses
`std::make_shared<>` directly (one allocation instead of two per depo).
For large depo sets (100k+ depos), the extra allocation per depo is
measurable.

Additionally, the two-pass approach (create all, then filter active) can be
consolidated. The first pass creates all depos, and the second pass filters
and resolves priors. This doubles the iteration count over the depo array.

---

## 2. DepoFileSource: Extra copy in transpose correction

**File:** `src/DepoFileSource.cxx:154-161`

```cpp
if (darr.rows() != iarr.rows() and darr.cols() == iarr.cols()) {
    array_xxfrw dtmp = darr.transpose();
    darr = dtmp;
    array_xxirw itmp = iarr.transpose();
    iarr = itmp;
}
```

The transpose creates a temporary view, which is then materialized into
`dtmp`, then copied into `darr`. With Eigen, `.transpose()` returns a lazy
expression. Using `darr = darr.transpose().eval()` would avoid the explicit
temporary, though Eigen's aliasing rules make this tricky. The current
approach is safe but uses 2x memory temporarily.

For large depo arrays this extra allocation may matter.

---

## 3. DepoFileSource: Unnecessary shared_ptr copy in IDepo vector construction

**File:** `src/DepoFileSource.cxx:234-238`

```cpp
WireCell::IDepo::vector idepos;
for (auto sdepo: sdepos) {
    idepos.push_back(sdepo);
}
```

This copies shared pointers one by one. Could use:
```cpp
WireCell::IDepo::vector idepos(sdepos.begin(), sdepos.end());
```
or `std::move_iterator` to avoid incrementing reference counts.

---

## 4. FrameFileSink: String concatenation in tag-join loop

**File:** `src/FrameFileSink.cxx:94-99`

```cpp
std::string comma="";
std::string stags="";
for (const auto& tag : m_tags) {
    stags += comma + tag;
    comma = ",";
}
```

Minor: repeated string concatenation with temporary. For typical tag counts
(< 10) this is negligible, but a `std::ostringstream` or
`fmt::join` would be cleaner and avoid temporaries.

---

## 5. FrameFileSink::one_tag: Full array allocation even for sparse frames

**File:** `src/FrameFileSink.cxx:167`

```cpp
Array::array_xxf arr = Array::array_xxf::Zero(nrows, ncols) + m_baseline;
```

This allocates a full dense 2D array for all channels x all ticks, even when
the traces are sparse. For detectors with many channels (e.g., ~10k) and
many ticks (e.g., ~10k), this is a 400MB float array. The `+ m_baseline`
also creates a temporary array of the same size.

Using `Array::array_xxf::Constant(nrows, ncols, m_baseline)` would avoid
the temporary from the addition.

---

## 6. FrameFileSource::load: Per-row ChargeSequence copy

**File:** `src/FrameFileSource.cxx:329`

```cpp
ITrace::ChargeSequence charges(row.data(), row.data() + ncols);
```

Each row of the Eigen array is copied into a separate `std::vector<float>`.
For large frames (10k channels x 10k ticks), this means ~10k individual
heap allocations. This is inherent to the ITrace interface (each trace owns
its data), but is worth noting as the dominant cost of frame loading.

---

## 7. NumpyFrameSaver: Reads config values on every call

**File:** `src/NumpyFrameSaver.cxx:76-81`

```cpp
const float baseline = m_cfg["baseline"].asFloat();
const float scale = m_cfg["scale"].asFloat();
const float offset = m_cfg["offset"].asFloat();
const bool digitize = m_cfg["digitize"].asBool();
const std::string fname = m_cfg["filename"].asString();
```

On every `operator()` call, config values are re-parsed from JSON. Most
other components parse these once in `configure()` and store in member
variables. With large frame counts this adds unnecessary overhead from
repeated JSON map lookups and string conversions.

---

## 8. NumpyFrameSaver: Verbose debug logging on every frame

**File:** `src/NumpyFrameSaver.cxx:91-105`

A `std::stringstream` is constructed and populated with frame/trace tag
information on every single frame, even when debug logging is disabled. This
includes iterating over all frame tags and trace tags. Using the spdlog
lazy-evaluation pattern (like `log->debug(...)`) would avoid the string
construction when logging is disabled.

---

## 9. BeeDepoSource: Entire file loaded then reversed

**File:** `src/BeeDepoSource.cxx:40-57`

The entire JSON file is loaded, all depos are created in forward order, then
the vector is reversed so they can be served via `pop_back()`:

```cpp
m_depos.resize(ndepos, nullptr);
for (int idepo = 0; idepo < ndepos; ++idepo) {
    m_depos[idepo] = ...;
}
...
std::reverse(m_depos.begin(), m_depos.end());
```

Similarly, `m_filenames` is reversed at configure time. While this works, a
front-to-back index counter would avoid the `O(n)` reverse operations.

---

## 10. TensorFileSink::jsonify: Extra string copy via stringstream

**File:** `src/TensorFileSink.cxx:67-68`

```cpp
std::stringstream ss;
ss << md;
auto mdstr = ss.str();
```

The JSON is first written to a stringstream, then copied to a string via
`ss.str()`. The string is then written to the output stream. This creates
an extra copy of the JSON text. Writing directly to the output would avoid
this, though the custard protocol requires knowing the size beforehand, so
two passes (or a pre-serialization to string) is necessary.

---

## 11. ClusterFileSink: Blob validation loop on every cluster

**File:** `src/ClusterFileSink.cxx:147-157`

```cpp
for (auto vtx : mir(boost::vertices(gr))) {
    const auto& node = gr[vtx];
    if (node.code() == 'b') {
        auto iblob = std::get<IBlob::pointer>(gr[vtx].ptr);
        Aux::BlobCategory bcat(iblob);
        if (bcat.ok()) continue;
        log->warn("malformed blob: {}", bcat.str());
    }
}
```

This debug/validation loop iterates over ALL vertices in the cluster graph
on every call, even in production. The comment says "fixme: debugging" but
it remains enabled. For large clusters with thousands of vertices, this
adds overhead. The same loop exists in `ClusterFileSource.cxx:343-353`.

---

## 12. FrameFileSink: channels vector contains duplicates then unique'd

**File:** `src/FrameFileSink.cxx:156-160`

```cpp
auto tmp = Aux::channels(traces);
std::sort(tmp.begin(), tmp.end());
auto chbeg = tmp.begin();
auto chend = std::unique(chbeg, tmp.end());
channels.insert(channels.begin(), chbeg, chend);
```

The channels are extracted from all traces, sorted, unique'd, then
copied. For frames with many traces per channel this creates a potentially
large temporary vector. A `std::set` or `std::unordered_set` during
extraction could avoid the sort+unique, though the current approach is
reasonable for typical workloads.

---

## 13. NumpyDepoSaver: Entire depo stream buffered in memory

**File:** `src/NumpyDepoSaver.cxx:66-72`

All depos are buffered in `m_depos` until EOS:

```cpp
if (indepo) {
    outdepo = indepo;
    m_depos.push_back(indepo);
    return true;
}
```

For very large events this means the entire depo set is held in memory until
the end-of-stream marker arrives, preventing streaming/incremental output.
This is documented behavior ("This saver will buffer depos in memory until
EOS is received"), but worth noting for memory-constrained scenarios.
