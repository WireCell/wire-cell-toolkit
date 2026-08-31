# BlipSource

Generates a sequence of point-like energy depositions (blips) with configurable spatial, time, and charge distributions, suitable for simulating diffuse backgrounds such as Ar-39 decays.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `BlipSource` |
| Concrete class | `WireCell::Gen::BlipSource` |
| Node category | source |
| Primary interface | `IDepoSource` |
| Input type(s) | `(none)` |
| Output type(s) | `IDepo` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `rng` | type:name of the IRandom service (default: "Random") |
| `charge` | object describing the charge per blip; "type" is "mono" (fixed value via "value") or "pdf" (sampled from "edges"/"pdf" arrays) |
| `time` | object with "type" ("decay"), "start", "stop", and "activity" (in Bq) controlling time distribution |
| `position` | object describing spatial distribution; "type" "box" with sub-key "extent" specifying a bounding volume |
