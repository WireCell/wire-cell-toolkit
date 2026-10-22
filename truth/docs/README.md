# WireCellTruth — native MC-truth plumbing

Carries Geant4 identity (track id, pdg, energy, undrifted position) from the
simulation depos through imaging and the point tree, so that per-blob and
per-point truth labels, particle-flow truth tables and label QA are produced
inside the toolkit without LArSoft.

The specification is `wcp-porting-img/pr_ml/docs/19_native_truth_plumbing.md`
(round 4 of the PR-chain brainstorm, 2026-10-21); this directory holds the
per-component documentation as the components land.

## Status

| Step (doc 19 sec 9) | Content | State |
|---|---|---|
| T1 | package skeleton: `wscript_build`, `CMakeLists.txt`, `Schema.h`, `wcdoctest-truth` | done 2026-10-21 |
| T2 | depo identity round trip: `Aux::fill` child-index defect fix, `DepoFileSource` legacy fallback, `DepoFileSink.truth_sidecar` | done 2026-10-21 (in `aux/`, `sio/`) |
| T3 | `DepoBlobAssoc`, `BlobDepoLabel`, `BlobDepoLabelTensor` (+ `iface/IClusterDepoTensor`) | done 2026-10-21 (`docs/components/`) |
| T4 | `TruthTableSource`, `python/dump_truth_table.py` | pending |
| T5 | `PFMerger`, `DepoPolyline`, `Matching` | pending |
| T6 | `TruthPointLabeler` | pending |
| T7–T9 | `validate_2d.py`, `validate_3d.py`, `qa_report.py` | pending |

## Rules

- Depends on `WireCellUtil`, `WireCellIface`, `WireCellAux` only; never on another plugin.
  Pure algorithms are header-only so `clus` visitors can include them without a link edge.
- Every behaviour change in production packages is behind a default-OFF knob; the only
  exception is the `Aux::fill` child-index defect fix (doc 19 sec 3.3, owner decision O2),
  which changes the `child` column of depo files written with priors.
- `Schema.h` is the contract: names, orders, sentinels (`tid_ghost = -1`, `tid_shared = -2`),
  class codes. Any change bumps `schema_version`.

## Depo file format (after T2)

`DepoFileSink` writes per depo set `depo_data_<ident>.npy` (N,7) float32
`t, q, x, y, z, long, tran` and `depo_info_<ident>.npy` (N,4) int32
`id, pdg, gen, child`; with `truth_sidecar: true` also `depo_truth_<ident>.npy`
(N,3) float64 `energy, orig_id, reserved`. Rows are pushed depth-first: an active
depo (gen 0) is followed by its prior (gen 1), that one's prior (gen 2), etc.; the
`child` of a gen>0 row is the row of the depo whose `prior()` it is. Files written
before 2026-10-21 have `child == own row` on every prior row; `DepoFileSource` and
`NumpyDepoSetLoader` detect that and use "child = previous row".
