# BlobSetReframer

Projects each blob in a BlobSet into frame traces, distributing blob charge uniformly over its wire and time extent, with optional per-anode trace tagging and selectable charge source.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `BlobSetReframer` |
| Concrete class | `WireCell::Img::BlobSetReframer` |
| Node category | function |
| Primary interface | `IBlobSetFramer` |
| Input type(s) | `IBlobSet` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `frame_tag` | tag applied to the output frame (default: "") |
| `tick` | output frame sampling period in time units (default: 0.5 us) |
| `source` | selects charge source; "blob" uses blob value/error, "activity" uses slice channel activity (default: "blob") |
| `measure` | selects which scalar to use as charge; "value" or "error" (default: "value") |

## Time origin

`ISlice::start()` is relative to the time of the frame the slice was made
from (see `ISlice.h`).  `MaskSlicer` and `SumSlicer` both follow this.
`BlobSetReframer` therefore:

- sets the output frame time to the time of the blobs' slice frame
  (`slice->frame()->time()`), and raises an error if the blobs disagree;
- places each trace at `tbin = round(slice start / tick)`, with `tick` the
  configured output tick.

When the blobs come straight from a slicer in the same job, the output frame
lines up tick-for-tick with the frame the blobs were imaged from.

### History

A 2023 trial (`SLICE_START_TIME_IS_RELATIVE`, commit `7330008`) tried to make
slice start times absolute so that blobs could be reframed without the
original `IFrame`.  A parallel change to `MaskSlicer` made its starts relative
again, so the trial only took effect in `SumSlicer`, leaving the two slicers
inconsistent.  The trial was removed, and the relative convention, which most
consumers of `ISlice::start()` assume, is now used throughout.  Before that,
`BlobSetReframer` used the earliest slice start as the output frame time,
which dropped the real frame time.

## Known limitation: clusters read from file

`BlobSetReframer` is most useful in a job that reads clusters from file (for
example `img/test/tenio-img-brf.jsonnet`).  Cluster files do not store the
frame time or tick.  `Aux::ClusterArrays::to_cluster()` rebuilds each slice's
frame from only its ident, with time 0 and the default tick.  (It used to
take the earliest slice start as the frame time, which double counts that
start under the relative convention.)  So for file-loaded clusters the output
frame time is 0, and trace times are relative to the original, unknown,
frame time.  `IBlobSet` and `ICluster` carry no time of their own; the only
route to it is `ISlice::frame()`.

Suggested remedies, in rough order of preference:

### B. Store the frame header in cluster files

Save each frame ident's `time` and `tick` along with the cluster, mirroring
the `tickinfo` of frame files (`aux/docs/frame-files.org`, `tickinfo.org`).
For example, add a small per-frame array or a metadata entry to the cluster
npz, and an entry to the tensor data model's cluster metadata
(`aux/src/TensorDMcluster.cxx`).  `to_cluster()` would then rebuild frames
with their true time and tick, which fixes every consumer of loaded clusters
and not just reframers.  Old files would keep the current behavior (time 0),
with a warning.  This needs no change to `ISlice` or to the data-flow graph.

### C. Rejoin the frame header with the blobs in the graph

Add a join node, shaped like `IFrameJoiner`
(`IJoinNode<tuple<IFrame,IFrame>,IFrame>`), taking a blob set (or cluster)
plus a frame and producing a frame.  Feed it from the frame file the blobs
were imaged from, read header-only (no traces).  In a single job, a fanout
can carry the signal frame's header past the slicing and imaging.  This needs
no file format change and is a natural fit when the reframed output is to be
overlaid on or subtracted from the original frame.  The new cost is keeping
the two streams in step: matching idents, and EOS handling through fanouts
such as `BlobDeclustering` or per-anode splits.

### D. Make the frame header a first-class object

Each slice currently carries its own frame header copy (see the `fixme` about
`TracelessFrame` in `img/src/MaskSlice.cxx`).  Remedy B makes that header
persistent.  If C is built, the frame header could become its own small data
type passed through the graph instead of being reached only through
`ISlice::frame()`.
