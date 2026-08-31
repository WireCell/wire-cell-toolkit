# FrameMerger

Merges two input frames into a single output frame by combining traces according to a configurable merge map and rule; the "replace" rule keeps the last trace per channel, while "include" keeps all traces from both frames.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `FrameMerger` |
| Concrete class | `WireCell::SigProc::FrameMerger` |
| Node category | join |
| Primary interface | `IFrameJoiner` |
| Input type(s) | `(IFrame, IFrame)` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `mergemap` | array of [tag0, tag1, output_tag] triples specifying which tagged trace sets to merge and what output tag to assign |
| `rule` | merge algorithm, "replace" (frame 0 takes precedence) or "include" (all traces kept) (default "replace") |
