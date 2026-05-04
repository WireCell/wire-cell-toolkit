# Efficiency Examination of iface

This document records efficiency observations for the `iface/` directory.
Since iface is predominantly an interface library (abstract base classes),
most performance-critical code lives in implementing packages.  However,
several default implementations and inline helpers here do affect runtime.

---

## EFF-1: In-place transpose uses O(N) auxiliary memory (IDFT.cxx:64-82)

```cpp
// In-place transpose
const int size = nrows*ncols;
std::vector<bool> visited(size);  // O(N) allocation every call
```

**Observation**: The in-place transpose allocates a `std::vector<bool>` of
size `nrows * ncols` on every call.  For large 2D arrays (e.g., a
frame with 2400 channels x 10000 ticks = 24M elements), this is a
non-trivial allocation.

**Note**: The comment says "Implementations, please override if you can
offer something faster."  This is a fallback default.  If concrete
implementations (e.g., FFTW, cuFFT) override the transpose, this
code is never called on the hot path.

**Recommendation**: If this default is used in production, consider a
cache-oblivious transpose algorithm or block-transposition to improve
cache behavior.  Alternatively, allocate the visited vector once and
reuse.

---

## EFF-2: `fwd1b`/`inv1b` with `axis=0` does two full transposes (IDFT.cxx:15-43)

```cpp
void IDFT::fwd1b(const complex_t* in, complex_t* out,
                 int nrows, int ncols, int axis) const
{
    if (axis) {
        // axis=1: simple row-wise FFT loop
    }
    else {
        this->transpose(in, out, nrows, ncols);       // O(N) copy
        this->fwd1b(out, out, ncols, nrows, 1);       // N FFTs
        this->transpose(out, out, ncols, nrows);       // O(N) in-place
    }
}
```

**Observation**: Performing column-wise batch FFT via
transpose-FFT-transpose involves two full matrix transposes.  Each
transpose has poor cache locality for large matrices.  The total
extra memory traffic is ~4x the matrix size (2 transposes, each
reading and writing every element).

**Note**: Again, this is a fallback.  The comment says to override
for batch optimization.  FFTW and GPU implementations should
provide native column-wise batch FFTs without transposing.

**Recommendation**: If the default is ever used for axis=0, consider
strided FFT (if the underlying `fwd1d` supports it) to avoid transposing.

---

## EFF-3: `cluster_node_t::ident()` creates vector of closures per call (ICluster.h:113-129)

```cpp
int ident() const
{
    using ident_f = std::function<int()>;
    std::vector<ident_f> ofs {
        [](){return 0;},
        [&](){return std::get<channel_t>(ptr)->ident();},
        [&](){return std::get<wire_t>(ptr)->ident();},
        [&](){return std::get<blob_t>(ptr)->ident();},
        [&](){return std::get<slice_t>(ptr)->ident();},
        [&](){return std::get<meas_t>(ptr)->ident();},
    };
    ...
    return ofs[ind]();
}
```

**Observation**: Every call to `ident()`:
1. Constructs 6 `std::function` objects (which may heap-allocate)
2. Constructs a `std::vector` (heap allocation)
3. Invokes one closure, then destroys all 6 + the vector

This is called in graph traversals where `ident()` may be invoked
thousands of times.

**Recommendation**: Replace with a `switch` statement or `std::visit`:

```cpp
int ident() const {
    return std::visit([](auto&& arg) -> int {
        using T = std::decay_t<decltype(arg)>;
        if constexpr (std::is_same_v<T, size_t>) return 0;
        else return arg->ident();
    }, ptr);
}
```

This is zero-allocation and likely inlineable.

---

## EFF-4: `cluster_node_t` hash uses pointer-to-integer cast (ICluster.h:136-162)

```cpp
std::size_t operator()(const WireCell::cluster_node_t& n) const
{
    size_t h = 0;
    switch (n.ptr.index()) {
    case 1: h = (std::size_t) std::get<1>(n.ptr).get(); break;
    case 2: h = (std::size_t) std::get<2>(n.ptr).get(); break;
    ...
    }
    return h;
}
```

**Observation**: Using a raw pointer value as a hash is fast but produces
poor hash distribution if pointers are aligned (e.g., 8-byte or
16-byte aligned objects have bottom 3-4 bits always zero).  For
hash tables with power-of-2 bucket counts, this can cause clustering.

**Recommendation**: Apply a finalizer mix (e.g., shift and XOR):

```cpp
h = h ^ (h >> 16);  // simple mix for better distribution
```

Or use `std::hash<void*>` which many implementations already apply
such mixing.

---

## EFF-5: `IBlobSet::shapes()` returns by value (IBlobSet.cxx:5-12)

```cpp
WireCell::RayGrid::blobs_t WireCell::IBlobSet::shapes() const
{
    RayGrid::blobs_t ret;
    for (const auto& ib : blobs()) {
        ret.push_back(ib->shape());
    }
    return ret;
}
```

**Observation**: Returns a vector by value, rebuilding it from scratch
every call.  No `reserve()` is called, so the vector may reallocate
multiple times during `push_back`.

**Recommendation**: Add `ret.reserve(blobs().size())` before the loop.
However, NRVO should optimize the return-by-value.

---

## EFF-6: `IHydraNodeVV::operator()` copies entire input queues (IHydraNode.h:209-231)

```cpp
virtual bool operator()(any_queue_vector& anyinqs, any_queue_vector& anyoutqs)
{
    const size_t isize = HydraInput::input_types().size();
    input_queues iqs(isize);
    for (size_t ind=0; ind<isize; ++ind) {
        iqs[ind] = boost::any_cast<input_queue>(anyinqs[ind]);  // COPY
    }
    ...
}
```

**Observation**: `boost::any_cast<input_queue>(anyinqs[ind])` extracts
the typed queue by **value**, copying all shared pointers in each
deque.  For large queues this is expensive.  Furthermore, any changes
the body function makes to `iqs` (popping consumed elements) are lost
since the typed queues are copies.

**Impact**: This means the hydra dispatch framework copies input queues
on every invocation.  For the VV variant, this also means consumed
elements are NOT removed from the `anyinqs`.

**Recommendation**: Use `boost::any_cast<input_queue&>` to get a reference
where possible, or document that the VV adapter requires the body
to manage the any queues directly.

---

## EFF-7: `IHydraNodeVT::operator()` iterates input elements individually (IHydraNode.h:293-314)

```cpp
for (size_t ind=0; ind<isize; ++ind) {
    for (auto any : anyinqs[ind]) {       // iterates each element
        iqs[ind].push_back(boost::any_cast<input_pointer>(any));
    }
}
```

**Observation**: This copies each element individually via `any_cast`,
rather than casting the entire deque at once.  This is the correct
approach when the any queue contains individual `boost::any` objects,
but it is slower than the VV approach for large queues because of
per-element `any_cast` overhead.

**Recommendation**: Minor.  The per-element approach is inherent to the
design where each element is wrapped in `boost::any`.

---

## EFF-8: `IRandom` uses `std::bind` for closures (IRandom.cxx)

```cpp
IRandom::int_func IRandom::make_binomial(int max, double prob)
{
    return std::bind(&IRandom::binomial, this, max, prob);
}
```

**Observation**: `std::bind` returns an implementation-defined type that is
stored via `std::function`, which may heap-allocate.  Modern C++
lambdas are preferred and often avoid the heap allocation.

**Impact**: Very minor; these factory functions are called once during setup.

---

## EFF-9: `boost::any` type erasure overhead in all node adapters

The entire node framework (`IFunctionNode`, `ISinkNode`, `ISourceNode`,
`IQueuedoutNode`, `IFaninNode`, `IFanoutNode`, `IHydraNode`) uses
`boost::any` for type erasure at the graph-execution boundary:

```cpp
virtual bool operator()(const boost::any& anyin, boost::any& anyout)
{
    const input_pointer& in = boost::any_cast<const input_pointer&>(anyin);
    output_pointer out;
    bool ok = (*this)(in, out);
    if (!ok) return false;
    anyout = out;
    return true;
}
```

**Observation**: Each data transfer through a graph edge involves:
1. `boost::any_cast` (RTTI + type comparison) on the input side
2. `boost::any` construction (potential heap allocation) on the output side

For high-throughput pipelines processing millions of frames/slices,
this per-message overhead accumulates.

**Note**: This is an architectural choice that enables the flexible
node-graph framework.  The overhead is typically dominated by the
actual computation in the node body.

**Recommendation**: Consider `std::variant` for the small set of
commonly-used types, or `std::any` (C++17) which may be slightly
more optimized in some standard library implementations.

---

## Summary

| ID    | File           | Impact | Description                              |
|-------|----------------|--------|------------------------------------------|
| EFF-1 | IDFT.cxx       | Medium | O(N) alloc per in-place transpose call   |
| EFF-2 | IDFT.cxx       | Medium | Double transpose for axis=0 batch FFT   |
| EFF-3 | ICluster.h     | Medium | Vector+closures alloc per ident() call   |
| EFF-4 | ICluster.h     | Low    | Poor pointer hash distribution           |
| EFF-5 | IBlobSet.cxx   | Low    | Missing reserve() in shapes()            |
| EFF-6 | IHydraNode.h   | Medium | Copies entire input queues per dispatch  |
| EFF-7 | IHydraNode.h   | Low    | Per-element any_cast in VT adapter       |
| EFF-8 | IRandom.cxx    | Low    | std::bind vs lambda                      |
| EFF-9 | Node framework | Low    | boost::any overhead per graph edge        |
