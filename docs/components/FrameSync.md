# FrameSync

Merges multiple input frame streams into a single output stream ordered by frame ident number, forwarding the frame with the smallest ident at each step.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `FrameSync` |
| Concrete class | `WireCell::Aux::FrameSync` |
| Node category | hydra |
| Primary interface | `IFrameMerge` |
| Input type(s) | `vector<IFrame>` |
| Output type(s) | `IFrame` |
| Configurable | no |
