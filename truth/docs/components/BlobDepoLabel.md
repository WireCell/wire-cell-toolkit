# BlobDepoLabel

`BlobDepoFill` duplicated into `WireCellTruth` with per-track bookkeeping (`Truth::depo_blob_assoc`). Same interface and knobs; with the extra-smearing knobs off the output cluster is bit-identical to `BlobDepoFill`'s, so it can replace it in an existing truth tier under a byte gate. The per-track table is emitted by `BlobDepoLabelTensor`.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `BlobDepoLabel` |
| Concrete class | `WireCell::Truth::BlobDepoLabel` |
| Node category | join |
| Primary interface | `IBlobDepoFill` |
| Input type(s) | `ICluster`, `IDepoSet` |
| Output type(s) | `ICluster` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `time_offset` | as BlobDepoFill: added to depo times before slicing (default 0) |
| `nsigma` | as BlobDepoFill: Gaussian extent considered, in sigma (default 3) |
| `speed` | as BlobDepoFill: drift speed converting the longitudinal extent to time (required) |
| `pindex` | as BlobDepoFill: primary plane the depos are diced onto (default 2) |
| `id_from_prior` | track key and energy taken from the depo's `prior()` (the undrifted G4 depo) when it exists (default true) |
| `extra_smear_long` | width in time added in quadrature to the longitudinal extent; 0 = off, integration exactly BlobDepoFill's (default 0) |
| `extra_smear_tran` | width in distance added in quadrature to the transverse extent; 0 = off (default 0) |

Design and parity gate: `wcp-porting-img/pr_ml/docs/19_native_truth_plumbing.md` sec 3.5.
