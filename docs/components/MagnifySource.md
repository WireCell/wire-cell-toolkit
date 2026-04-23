# MagnifySource

Reads waveform data from a ROOT file in the Magnify format, loading per-plane 2D histograms (hu_*, hv_*, hw_*) and optionally bad-channel TTrees to reconstruct a tagged IFrame.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `MagnifySource` |
| Concrete class | `WireCell::Root::MagnifySource` |
| Node category | source |
| Primary interface | `IFrameSource` |
| Input type(s) | `(none)` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `filename` | path to the input Magnify ROOT file (required) |
| `frames` | array of frame/trace tags to load; each tag maps to histograms named h[uvw]_<tag> (default ["raw"]) |
| `cmmtree` | array of [cmmkey, treename] pairs mapping channel mask keys to input TTree names |
