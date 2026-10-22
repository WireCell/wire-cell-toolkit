# DepoFileSink

Writes IDepoSet contents to a streaming archive file as pairs of NumPy arrays (depo_data and depo_info) encoding deposition kinematics and metadata, one pair per IDepoSet received.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `DepoFileSink` |
| Concrete class | `WireCell::Sio::DepoFileSink` |
| Node category | sink |
| Primary interface | `IDepoSetSink` |
| Input type(s) | `IDepoSet` |
| Output type(s) | `(none)` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `outname` | output archive file path; extension determines compression format (required) |
| `truth_sidecar` | if true, also write `depo_truth_<ident>.npy`, double (ndepos, 3) = (energy, orig_id, reserved), after `depo_info`; default false, output byte-identical to before the knob existed |

## Array layout

`depo_data_<ident>.npy` is float32 (ndepos, 7) = (t, q, x, y, z, long, tran); `depo_info_<ident>.npy` is int32 (ndepos, 4) = (id, pdg, gen, child). Rows are pushed depth-first: an active depo (gen 0) is followed by its prior (gen 1), that one's prior (gen 2) and so on; `child` of a gen>0 row is the row index of the depo whose `prior()` it is. Files written before 2026-10-21 carried `child == own row` on every prior row (defect in `Aux::fill`); readers detect that pattern and fall back to "child = previous row".
