# DepoFileSource

Reads IDepoSet objects from a streaming archive file containing paired NumPy arrays (depo_data and depo_info), reconstructing depositions with prior-depo relationships and applying an optional charge scale factor.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `DepoFileSource` |
| Concrete class | `WireCell::Sio::DepoFileSource` |
| Node category | source |
| Primary interface | `IDepoSetSource` |
| Input type(s) | `(none)` |
| Output type(s) | `IDepoSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `inname` | input archive file path (required) |
| `scale` | multiplicative factor applied to each deposition charge (default 1.0) |

## Notes

Each set is the pair `depo_data_<ident>.npy` / `depo_info_<ident>.npy`, optionally followed by `depo_truth_<ident>.npy` (energy, orig_id, reserved) written by `DepoFileSink` with `truth_sidecar: true`; when present the depo energy is restored, otherwise it is the `SimpleDepo` default (1.0). Prior chains are rebuilt from the `child` column; legacy files whose prior rows reference themselves (written before 2026-10-21) are read through a "child = previous row" fallback with a warning.
