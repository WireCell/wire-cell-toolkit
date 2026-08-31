# Algorithm and Architecture: iface

This document explains the general design, architectural patterns, and
specific algorithmic details of the Wire-Cell Toolkit `iface/` package.

---

## 1. Overview

The `iface` package defines the **abstract interface layer** of the
Wire-Cell Toolkit (WCT).  It contains 134 header files and 12 source
files.  Its purpose is to declare contracts (pure virtual interfaces)
that concrete implementations in other packages (`sigproc`, `img`,
`gen`, `aux`, etc.) fulfill.

The package depends only on `WireCellUtil` and defines no heavy
computation itself.  It is the backbone of WCT's plugin architecture:
components are loaded at runtime by type name and communicate through
these interfaces.

---

## 2. Core Design Patterns

### 2.1 CRTP Data Interface (`IData<T>`)

**File**: `IData.h`

All data types (IDepo, IFrame, ITrace, IBlob, ISlice, ICluster, etc.)
inherit from `IData<T>` using the Curiously Recurring Template Pattern:

```cpp
class IDepo : public IData<IDepo> { ... };
```

`IData<T>` provides standardized type aliases:
- `pointer` = `std::shared_ptr<const T>` (immutable shared ownership)
- `vector` = `std::vector<pointer>`
- `shared_vector` = `std::shared_ptr<const vector>`

**Key design decision**: All data objects are immutable after creation.
They are always accessed through `shared_ptr<const T>`.  This enables
safe sharing across threads without locks.

### 2.2 Component Interface (`IComponent<T>`)

**Files**: Various (IDFT.h, IFieldResponse.h, IAnodePlane.h, etc.)

Service-like singletons (DFT engine, noise database, field response,
anode geometry) inherit from `IComponent<T>` (defined in
`WireCellUtil`).  This hooks into WCT's dependency injection system
where components are resolved by type and instance name at runtime.

### 2.3 Node-Graph Architecture

**Files**: `INode.h`, `ISourceNode.h`, `ISinkNode.h`, `IFunctionNode.h`,
`IQueuedoutNode.h`, `IFaninNode.h`, `IFanoutNode.h`, `ISplitNode.h`,
`IJoinNode.h`, `IHydraNode.h`

WCT processes data through a directed graph of nodes.  Each node has a
**category** that defines its I/O pattern:

```
sourceNode:     () -> T              One output, no input
sinkNode:       T -> ()              One input, no output
functionNode:   T1 -> T2             One-to-one transform
queuedoutNode:  T1 -> queue<T2>      One input, zero-or-more outputs
faninNode:      vector<T1> -> T2     N inputs (same type) to one output
fanoutNode:     T1 -> vector<T2>     One input to N outputs (same type)
joinNode:       tuple<T1,...> -> T    N inputs (different types) to one
splitNode:      T -> tuple<T1,...>    One input to N outputs (different types)
hydraNode:      queues -> queues     Most general: M input queues to N output queues
```

**Type erasure**: The graph executor sees only `INode` base pointers
and uses `boost::any` to pass data between nodes.  Each typed node
template (e.g., `IFunctionNode<IFrame, IFrame>`) provides an adapter
that casts between `boost::any` and the concrete pointer type:

```cpp
// IFunctionNode adapter (simplified)
bool operator()(const boost::any& anyin, boost::any& anyout) {
    auto in = boost::any_cast<input_pointer>(anyin);
    output_pointer out;
    bool ok = (*this)(in, out);  // call the typed implementation
    anyout = out;
    return ok;
}
```

This allows the graph framework to be type-agnostic while individual
nodes are fully typed.

---

## 3. Data Model

### 3.1 The Signal Processing Chain

The data model follows the physical signal processing chain in a
Liquid Argon Time Projection Chamber (LArTPC):

```
IDepo (charge deposition)
  |  drift, diffusion
  v
IFrame (digitized waveforms per channel)
  |  contains
  +---> ITrace (charge vs time on one channel)
  |  slice in time
  v
ISlice (channel activity in a time window)
  |  tiling / blob-finding
  v
IBlob (2D region in wire-plane space)
  |  group by time
  v
IBlobSet (collection of blobs in one time slice)
  |  clustering
  v
ICluster (graph of blobs, slices, channels, wires, measures)
```

### 3.2 IDepo (Charge Deposition)

**Files**: `IDepo.h`, `IDepo.cxx`

Represents a point-like energy deposition in the detector.  Key fields:
- `pos()`: 3D position
- `time()`: deposition time
- `charge()`: number of ionization electrons
- `energy()`: deposited energy (MeV)
- `prior()`: pointer to the original (pre-drift) deposition
- `extent_long()`, `extent_tran()`: Gaussian sigma of the charge cloud

The `prior()` chain allows tracing a drifted deposition back to its
origin.  `depo_chain()` walks this chain and returns the full history.

**Comparators**: `ascending_time`, `descending_time` sort depositions
by time with x-position tiebreaker.  `IDepoDriftCompare` sorts by
drift-adjusted time `t + x/v_drift`, useful for ordering depositions
as they would arrive at the anode.

### 3.3 IFrame and ITrace

**Files**: `IFrame.h`, `ITrace.h`

An `IFrame` represents a time window of digitized data across all channels.
It contains:
- A vector of `ITrace` objects (each trace is charge vs time on one channel)
- A reference time and tick period
- A tagging system: traces can be tagged (e.g., "raw", "gauss", "wiener")
  and retrieved by tag

An `ITrace` is a contiguous `std::vector<float>` of charge measurements
starting at time bin `tbin` on channel `channel()`.  Traces are sparse:
they only exist where there is data, and `tbin` locates them in the
frame's time axis.

### 3.4 ISlice

**Files**: `ISlice.h`

A time slice maps channels to (value, uncertainty) pairs over a time
span.  It uses an `unordered_map` keyed by `IChannel::pointer` with
custom hash/equality based on channel ident number (not pointer identity):

```cpp
typedef std::unordered_map<IChannel::pointer, value_t, IdentHash, IdentEq> map_t;
```

This is important: two different `IChannel` smart pointers pointing to
channels with the same `ident()` are considered the same key.

### 3.5 IBlob and IBlobSet

**Files**: `IBlob.h`, `IBlobSet.h`

A blob is a 2D region in the plane transverse to the drift direction,
defined by the intersection of wire ranges across multiple planes.
It carries:
- `value()`: associated charge
- `uncertainty()`: charge uncertainty
- `shape()`: the geometric description via `RayGrid::Blob`
- `face()`: the anode face it belongs to
- `slice()`: the time slice it was found in

An `IBlobSet` groups blobs from a single time slice.

### 3.6 ICluster (Graph-based Data Structure)

**Files**: `ICluster.h`, `ICluster.cxx`

The most complex data structure in iface.  An `ICluster` is a
**boost::adjacency_list** graph where nodes are typed using `std::variant`:

```cpp
using ptr_t = std::variant<size_t, channel_t, wire_t, blob_t, slice_t, meas_t>;
```

The variant index determines the node type, encoded as single characters:
`'c'` (channel), `'w'` (wire), `'b'` (blob), `'s'` (slice), `'m'` (measure).

**Graph type**: `boost::adjacency_list<setS, vecS, undirectedS, cluster_node_t>`
- `setS` for edge container: no parallel edges
- `vecS` for vertex container: O(1) vertex access by descriptor
- `undirectedS`: edges are bidirectional

**Helper functions**:
- `oftype<T>(graph)`: extract all vertices of a given type
- `neighbors_oftype<T>(graph, node)`: get typed neighbors of a node
- `neighbors(graph, vertex)`: get all neighbor vertex descriptors

These helpers iterate all vertices/edges and filter by type, which is
O(V) or O(degree) respectively.

---

## 4. Detector Geometry

### 4.1 Wire Plane Identification (`WirePlaneId`)

**Files**: `WirePlaneId.h`, `WirePlaneId.cxx`

Wire planes are identified by a bit-packed integer encoding three fields:

```
Bits:  [apa_shift+:apa] [face_shift:face] [0-2:layer]
       [4+:apa]         [3:face]          [0-2:layer_mask]
```

Layer values are **bitmasks**, not indices:
- `kUlayer = 1`, `kVlayer = 2`, `kWlayer = 4`, `kAllLayers = 7`

This allows OR-combination of layers (e.g., for selecting all wires).
The `index()` method converts layer enum to 0/1/2 index.

### 4.2 IAnodePlane and IAnodeFace

**Files**: `IAnodePlane.h`, `IAnodeFace.h`

`IAnodePlane` represents an entire anode assembly (e.g., one APA in DUNE).
It contains one or more `IAnodeFace` objects (DUNE has 2 faces per APA,
MicroBooNE has 1).

Each face provides:
- Wire planes and their geometry
- A `RayGrid::Coordinates` for ray-based tiling operations
- A bounding box for the sensitive volume
- Direction indicator (`dirx()`: +1 or -1) for drift direction

### 4.3 IWire and Wire Selectors

**Files**: `IWire.h`, `IWireSelectors.h`, `IWireSelectors.cxx`

`IWire` represents a physical wire segment with:
- Global ident, plane-local index, channel number
- Segment number (distance from channel input, for wrapped wires)
- Geometric ray (start and end points)

`WirePlaneSelector` is a functor that filters wires by layer/face/apa:
```cpp
bool operator()(IWire::pointer wire) {
    WirePlaneId ident = wire->planeid();
    if (layers && !(layers & ident.ilayer())) return false;
    if (apa >= 0 && ident.apa() != apa) return false;
    if (face >= 0 && ident.face() != face) return false;
    return true;
}
```

Pre-built selectors are provided: `select_u_wires`, `select_v_wires`,
`select_w_wires`, `select_all_wires`.

---

## 5. Signal Processing Interfaces

### 5.1 IDFT (Discrete Fourier Transform)

**Files**: `IDFT.h`, `IDFT.cxx`

The DFT interface provides 6 methods as the outer product of:
- Direction: `fwd` (forward, no normalization), `inv` (inverse, 1/N normalization)
- Dimensionality: `1d` (single transform), `1b` (batched along one axis),
  `2d` (full 2D transform)

All operate on `complex<float>` arrays in row-major (C) order.

**Default batch implementation** (`fwd1b`, `inv1b`):
- For `axis=1`: simple loop over rows, calling `fwd1d`/`inv1d` per row
- For `axis=0`: transpose -> batch-along-rows -> transpose back

**Default transpose** (`transpose`):
- Out-of-place: simple nested loop with index transformation
  `out[col*nrows + row] = in[row*ncols + col]`
- In-place: cycle-following algorithm from
  [Wikipedia](https://en.wikipedia.org/wiki/In-place_matrix_transposition)

Concrete implementations (FFTW, cuFFT) are expected to override these
defaults with optimized versions.

### 5.2 IChannelNoiseDatabase

**File**: `IChannelNoiseDatabase.h`

Provides per-channel noise characterization for signal processing:
- `rcrc(channel)`: RC coupling response filter (frequency domain)
- `config(channel)`: misconfiguration correction filter
- `noise(channel)`: noise attenuation filter
- `response(channel)`: nominal detector response spectrum
- `coherent_channels()`: channel groupings for coherent noise subtraction
- Various thresholds: `coherent_nf_decon_limit`, `min_rms_cut`, etc.

This interface is used by the noise filtering pipeline in `sigproc`.

### 5.3 IPlaneImpactResponse and IFieldResponse

**Files**: `IPlaneImpactResponse.h`, `IFieldResponse.h`

These interfaces provide access to the detector's field response:
- `IFieldResponse`: loads and shares the field response data
- `IPlaneImpactResponse`: provides frequency-domain response spectra
  at specific impact positions (positions along the pitch direction)

Key methods:
- `closest(relpitch)`: nearest impact response to a pitch position
- `bounded(relpitch)`: two impact responses bracketing a pitch position
  (used for interpolation in convolution/deconvolution)

---

## 6. Random Number Generation

**File**: `IRandom.h`, `IRandom.cxx`

Defines abstract random distribution methods:
- `binomial`, `poisson`, `normal`, `uniform`, `exponential`, `range`

Also provides factory methods (`make_binomial`, etc.) that return
`std::function` closures.  The factory pattern allows callers to
create a closure once and call it repeatedly without re-specifying
parameters.  Concrete implementations are encouraged to override
the factories to return closures that avoid repeated generator
construction.

---

## 7. Data Flow Graph Framework

### 7.1 Node Categories and Multiplicity

The `INode::NodeCategory` enum defines 10 node types.  The most commonly
used are:

- **functionNode**: 1-to-1 transforms (e.g., frame filtering, deconvolution)
- **queuedoutNode**: 1-to-many (e.g., frame slicing produces many slices)
- **faninNode**: N-to-1 with same type (e.g., merging blob sets from multiple faces)
- **fanoutNode**: 1-to-N with same type (e.g., sending a frame to multiple consumers)
- **hydraNode**: general M-to-N with queues (e.g., branch-merge patterns)

### 7.2 Hydra Node Design (IHydraNode.h)

The hydra node is the most flexible but also most complex node type.
It provides a 2x2 matrix of variants:

```
         Tuple (T)              Vector (V)
Input:   IHydraInputT           IHydraInputV
Output:  IHydraOutputT          IHydraOutputV
```

- **T (Tuple)**: ports have different types, cardinality is static
  (compile-time).  Access via `std::get<N>()`.
- **V (Vector)**: ports have the same type, cardinality is dynamic
  (runtime).  Access via `qs[N]`.

The four concrete combinations are:
- `IHydraNodeVV`: vector in, vector out (most common)
- `IHydraNodeTT`: tuple in, tuple out
- `IHydraNodeVT`: vector in, tuple out
- `IHydraNodeTV`: tuple in, vector out

Each provides an adapter from the untyped `any_queue_vector` interface
to the typed interface.

### 7.3 Concurrency Model

`INode::concurrency()` returns 1 by default.  The comments warn
strongly against setting higher concurrency due to:
1. Nodes typically have mutable state
2. Concurrent execution violates ordering assumptions
3. Historical issues (see issue #121)

---

## 8. Legacy / Deprecated Interfaces

### 8.1 ISequence (ISequence.h)

An iterator-based sequence abstraction using facade iterators.
It appears to be largely unused in favor of the simpler
`IData::shared_vector` pattern.  The `end()` methods contain bugs
(returning `cbegin()`) which further suggests this is dead code.

### 8.2 IWire-based interfaces

Several wire-centric interfaces (`IWireSource`, `IWireGenerator`,
`IWireSummarizer`) appear to be from an earlier design where wires
were first-class streaming data.  Modern WCT accesses wire geometry
through `IAnodeFace` and `IWirePlane`.

### 8.3 IDiffusion / IDiffuser

These appear to be from an older simulation approach.  Modern WCT
simulation uses `IDrifter` and `IDuctor` instead.

---

## 9. Configuration System

Interfaces that need configuration inherit from `IConfigurable`
(defined in `IConfigurable.h`), which requires:
- `default_configuration()`: return a JSON object with defaults
- `configure(cfg)`: accept a JSON configuration

This integrates with WCT's Jsonnet-based configuration system where
component parameters are specified in `.jsonnet` files and resolved
at runtime.

---

## 10. Key Architectural Properties

1. **Immutability**: All data objects are const-shared (`shared_ptr<const T>`).
   Once created, data never changes.  This is critical for thread safety.

2. **Plugin architecture**: Components are loaded by type name string.
   The interface package defines no concrete implementations, enabling
   clean separation of concerns.

3. **Minimal dependencies**: iface depends only on WireCellUtil and
   standard libraries (plus Boost for graph and any).  This keeps
   compile times down and avoids circular dependencies.

4. **Forward-compatible**: New node types and data types can be added
   by defining new interfaces without modifying existing ones.

5. **Graph-centric**: The data processing model is a directed graph
   of typed nodes, with type erasure at graph edges and full type
   safety within each node's implementation.
