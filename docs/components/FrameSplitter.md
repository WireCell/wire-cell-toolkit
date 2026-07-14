# FrameSplitter

Passes the same input frame pointer to both output ports, providing a simple no-copy split for use when the same frame data must flow to two downstream branches simultaneously.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `FrameSplitter` |
| Concrete class | `WireCell::SigProc::FrameSplitter` |
| Node category | split |
| Primary interface | `IFrameSplitter` |
| Input type(s) | `IFrame` |
| Output type(s) | `(IFrame, IFrame)` |
| Configurable | no |
