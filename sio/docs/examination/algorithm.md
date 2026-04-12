# SIO Package: Algorithm Documentation

## Overview

The SIO ("Simple I/O") package provides ROOT-independent serialization and
deserialization for Wire Cell Toolkit (WCT) data objects. It handles four
primary data types: **depositions** (energy deposits in the detector),
**frames** (digitized waveforms), **tensors** (generic multi-dimensional
arrays), and **clusters** (graph-structured reconstruction objects).

All components follow the WCT dataflow programming model: each component is
a node in a directed dataflow graph. **Sources** produce data (read from
files), **Sinks** consume data (write to files), and **Filters** transform
data (pass-through with side-effect I/O).

## Data Formats

### Streaming archive format (tar + compression)

The file-based components (DepoFile, FrameFile, TensorFile, ClusterFile) use
a **streaming tar archive** format. The tar stream is managed by the
`custard` library and the `boost::iostreams::filtering_stream` framework.

The pipeline looks like:

```
Data objects -> Numpy/JSON serialization -> custard tar framing -> compression -> file
```

Compression is determined by file extension:
- `.tar` - no compression
- `.tar.bz2` - bzip2 compression
- `.tar.gz` - gzip compression
- `.tar.xz` / `.txz` - xz compression (if available)

Within the tar stream, individual entries are Numpy `.npy` files (via the
`pigenc` library) or `.json` files. The custard protocol frames each entry
with a simple text header: `name <filename>\nbody <bytecount>\n<data>`.

### Numpy NPZ format

The Numpy-prefix components (NumpyDepoLoader/Saver, NumpyFrameSaver) use
the standard `.npz` format (ZIP archive of `.npy` files). This uses the
`cnpy` library for direct read/write of Numpy-compatible arrays.

### JSON format

JsonDepoSource and BeeDepoSource read JSON files directly. The Bee format
is a simple dictionary with parallel arrays (x, y, z, q, t).

---

## Component Details

### DepoFileSink / DepoFileSource

**Purpose:** Stream deposition sets to/from tar archives.

**Data layout:** Each depo set is written as two Numpy arrays:
- `depo_data_<ident>.npy` - shape `(N, 7)`, float32: columns are
  `[time, charge, x, y, z, extent_long, extent_tran]`
- `depo_info_<ident>.npy` - shape `(N, 4)`, int32: columns are
  `[id, pdg, generation, child_index]`

**Algorithm (Sink):**
1. Receive an `IDepoSet`
2. Call `Aux::fill()` to flatten the depo tree (including prior/child
   relationships) into two 2D arrays
3. Write both arrays as `.npy` entries in the tar stream
4. Flush the stream

**Algorithm (Source):**
1. Read two consecutive `.npy` files from the tar stream (data + info)
2. Allow either ordering (data first or info first) by checking the "tag"
   field in the filename
3. Validate array shapes match (rows equal, 7 data columns, 4 info columns)
4. Handle transposed arrays (a legacy compatibility feature): if rows don't
   match but cols do, transpose both arrays
5. **Two-pass depo reconstruction:**
   - Pass 1: Create all depos from the arrays
   - Pass 2: Link prior depos (generation > 0) to their children using the
     child_index column. Collect generation-0 (active) depos as output.
6. Wrap in `SimpleDepoSet` and return

**Prior/child depo model:** Depos form a tree where "active" depos
(generation=0) are the final results, and each may have a chain of "prior"
depos representing their ancestry (e.g., before diffusion). The info array
column 2 holds the generation number, and column 3 holds the index of the
child depo in the flat array. This allows full tree reconstruction.

---

### FrameFileSink / FrameFileSource

**Purpose:** Stream detector frames (waveforms) to/from tar archives.

**Data layout per tag per frame:**
- `frame_<tag>_<ident>.npy` - shape `(nchannels, nticks)`, float32 or int16
- `channels_<tag>_<ident>.npy` - 1D array of channel IDs
- `tickinfo_<tag>_<ident>.npy` - 3-element array: `[time, tick_period, tbin0]`
- `summary_<tag>_<ident>.npy` - optional per-trace summary values
- `chanmask_<maskname>_<ident>.npy` - shape `(N, 3)`, int32: columns
  `[channel_id, begin_bin, end_bin]`

**Algorithm (Sink):**
1. Receive an `IFrame`
2. For each configured tag (or "*" for all traces):
   a. Extract matching traces (by trace tag, then frame tag)
   b. Determine channel range and time bin range
   c. Allocate a dense 2D array `(nchannels, nticks)`
   d. Fill array from traces using `Aux::fill()`: each trace's samples are
      placed at the correct row (channel) and column offset (tbin)
   e. Apply linear transform: `arr = (arr + baseline) * scale + offset`
   f. Optionally cast to int16 ("digitize") for compact storage
   g. Write the frame array, channel list, tick info, and optional summary
3. Write channel mask maps if configured

**Dense mode:** When `dense` config is given, the array spans a fixed
channel and tick range (chbeg:chend, tbbeg:tbend) regardless of which
channels have data. This is useful for producing uniformly-shaped arrays.

**Algorithm (Source):**
1. Read `.npy` entries from the tar stream one at a time
2. Parse each filename to determine its type (frame, channels, tickinfo,
   summary, chanmask) and tag
3. Group entries by `<ident>` (frame identifier). A read-ahead mechanism
   detects when the ident changes, signaling the end of the current frame.
4. Build "framelets" - one per tag - containing the trace array, channel
   list, tick info, and summary
5. Tag filtering: if tags are configured, only matching tags produce tagged
   traces; others are either untagged or excluded
6. Construct `SimpleFrame`:
   a. Convert each row of each framelet's trace array into a `SimpleTrace`
   b. Tag traces within each framelet with the framelet's tag
   c. Apply frame-level tags
   d. Include channel mask maps

---

### TensorFileSink / TensorFileSource

**Purpose:** Stream generic tensor sets to/from tar or NPZ archives.

**Data layout per tensor set:**
- `<prefix>tensorset_<ident>_metadata.json` - tensor set metadata
- `<prefix>tensor_<ident>_<index>_metadata.json` - per-tensor metadata
- `<prefix>tensor_<ident>_<index>_array.npy` - per-tensor array data

**Algorithm (Sink):**
1. Receive an `ITensorSet`
2. Write the tensor set metadata as JSON
3. For each tensor in the set:
   a. Write tensor metadata as JSON
   b. Write tensor array data as `.npy` (via `Stream::write`)
4. Support "dump mode" that skips actual file writes (for debugging)

**Algorithm (Source):**
1. Parse filenames to determine type (set metadata, tensor metadata, tensor
   array), ident, and index
2. Group by ident. When a new ident is encountered, emit the previous set.
3. Construct `TFSTensor` objects from pigenc data (deep-copies the array
   data into a `std::vector<std::byte>` buffer)
4. Merge metadata JSON with array data for each tensor
5. Handle "array-less tensors" (tensors with metadata but no array data)

**Prefix mechanism:** The `prefix` config allows multiple sink/source pairs
to share the same archive file by using different prefixes to namespace
their entries.

---

### ClusterFileSink / ClusterFileSource

**Purpose:** Serialize/deserialize cluster graphs to/from tar archives.

**Supported formats:**
- `json` - Full cluster graph in JSON "cluster graph schema"
- `numpy` - Cluster arrays: separate `.npy` files for each node and edge type
- `dot` - GraphViz DOT format (sink only)

**Numpy array schema:**
Clusters are graphs with typed nodes and edges. In the numpy format:
- Node arrays: `<prefix>_<clusterID>_<code>nodes.npy` where `<code>` is a
  single character identifying node type:
  - `c` = channel, `b` = blob, `m` = measure, `s` = slice, `w` = wire
- Edge arrays: `<prefix>_<clusterID>_<code>edges.npy` where `<code>` is a
  two-character pair identifying the connected node types

**Algorithm (Sink):**
1. Receive an `ICluster`
2. Run blob validation (debug: checks all blob nodes for well-formedness)
3. Dispatch to format-specific serializer:
   - JSON: call `Aux::jsonify()` to convert graph to JSON, write to stream
   - Numpy: call `Aux::ClusterArrays::to_arrays()` to convert graph to
     typed arrays, write each array to stream
   - DOT: call `Aux::dotify()` for visualization format

**Algorithm (Source):**
1. Read filenames from tar stream
2. Parse filename to determine format (JSON vs Numpy) and cluster ID
3. Dispatch to format-specific loader:
   - JSON: parse JSON, use `ClusterLoader` to reconstruct graph
   - Numpy: accumulate all node/edge arrays for one cluster ID, then call
     `to_cluster()` to reconstruct the graph from arrays
4. Requires `IAnodePlane` instances (configured via "anodes") to resolve
   wire/channel references during reconstruction

---

### NumpyDepoLoader / NumpyDepoSaver

**Purpose:** Load/save depositions from/to `.npz` files (ZIP of Numpy).

These are the older NPZ-based counterparts to the streaming
DepoFileSource/Sink.

**Algorithm (Saver):**
1. Buffer all incoming depos in memory until EOS
2. On EOS, flatten the depo tree via `Aux::fill()` to get `(N,7)` data and
   `(N,4)` info arrays
3. Append both arrays to the `.npz` file using `cnpy::npz_save` in append
   mode
4. Array naming: `depo_data_<count>`, `depo_info_<count>`

**Algorithm (Loader):**
1. Attempt to load `depo_data_<count>` and `depo_info_<count>` from the
   `.npz` file
2. Reconstruct depos using `NumpyDepoTools::load()`:
   - Same two-pass prior-linking algorithm as DepoFileSource
3. Buffer depos in a `std::deque` and serve them one at a time
4. Send EOS (nullptr) after each depo set is drained
5. Try to load the next set (increment count) until no more arrays exist

**NumpyDepoSetLoader** is a variant that wraps each loaded set into an
`IDepoSet` instead of emitting individual depos.

---

### NumpyFrameSaver

**Purpose:** Save frames to `.npz` files.

Same core algorithm as FrameFileSink but writes to NPZ instead of tar:
1. Extract traces by tag
2. Build dense 2D array from traces
3. Apply baseline/scale/offset transform, optionally digitize
4. Save frame array, channel array, and tickinfo array using `cnpy::npz_save`

Key difference from FrameFileSink: NumpyFrameSaver is a **filter**
(IFrameFilter) - it passes the frame through unchanged while saving as a
side effect. FrameFileSink is a **sink** that consumes the frame.

---

### BeeBlobSink

**Purpose:** Write blob data to Bee visualization format (ZIP archive).

**Algorithm:**
1. Receive `IBlobSet` objects
2. Track frame identity changes via ident:
   - Same ident: accumulate points
   - New ident: flush accumulated points to Bee store, start new
3. For each blob set, run configured "samplers" (`IBlobSampler`) to convert
   blobs into point clouds (`Bee::dump`)
4. Accumulate sampled points in `Bee::Points`
5. On flush, write the accumulated points to the Bee ZIP archive via
   `Bee::Sink`

The Bee format organizes data by "events" (identified by run/subrun/event
numbers).

---

### BeeDepoSource

**Purpose:** Load depositions from Bee-format JSON files.

**Algorithm:**
1. Configured with a list of JSON filenames
2. Process files in reverse order (LIFO via pop_back)
3. For each file:
   - Load JSON (expects parallel arrays: x, y, z, q, t)
   - Create SimpleDepo for each entry
   - Optionally insert EOS (nullptr) between files (unless "stream" policy)
4. Serve depos one at a time via pop_back

---

### JsonDepoSource

**Purpose:** Load depositions from generic JSON files with recombination
model support.

**Algorithm:**
1. Load a JSON file, navigate to the "depos" key (configurable via
   "jsonpath")
2. For each depo entry, compute charge using a recombination model:
   - `electrons` model: charge = n * e_charge (just unit conversion)
   - `MipRecombination`: charge = model(dE) (point energy deposit)
   - `BirksRecombination` / `BoxRecombination`: charge = model(dE, dX)
     (step energy deposit)
3. Sort depos by descending time
4. Serve depos one at a time (LIFO via pop_back, but already reversed by
   the descending sort, so effectively time-ordered output)

---

## Architecture Notes

### Stream Protocol

The tar-stream components use a two-layer protocol:
1. **custard** - Provides tar-level framing (filename + size headers) and
   compression filter setup
2. **pigenc** - Provides Numpy `.npy` format read/write for array data

The read-ahead pattern (used by FrameFileSource, TensorFileSource,
ClusterFileSource) is necessary because these components must detect when
the stream transitions from one logical unit (frame, tensor set, cluster) to
the next. Since the tar stream is sequential, the only way to know a unit
is complete is to read the first entry of the next unit and cache it.

### EOS (End of Stream) Protocol

All source components follow the WCT EOS protocol:
1. Return `true` with a valid pointer for each data object
2. Return `true` with `nullptr` to signal end-of-stream
3. Return `false` on subsequent calls after EOS

Sink components always return `true` (both for data and EOS). The EOS signal
triggers any buffered data to be flushed.

### Component Registration

All components are registered with the WCT factory system via
`WIRECELL_FACTORY` macros, specifying the interfaces they implement
(e.g., `IDepoSetSink`, `IConfigurable`, `ITerminal`). Components that
manage output streams implement `ITerminal` with a `finalize()` method
to properly close the stream.
