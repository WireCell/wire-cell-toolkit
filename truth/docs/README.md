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
| T4 | `TruthTable`, `TruthTableSource`, `python/{truth_table,dump_truth_table}.py` | done 2026-10-21 (`docs/components/TruthTableSource.md`) |
| T5 | `PFMerger`, `DepoPolyline`, `Matching` (header-only, `inc/`) | done 2026-10-21 |
| T6 | `TruthPointLabel` (algorithm), `TruthPointLabeler` (+ `iface/ITruthPointLabeler`), `cfg/pgrapher/common/truth.jsonnet`, PDVD `clus.jsonnet` args | done 2026-10-21 (`docs/components/TruthPointLabeler.md`) |
| T7 | `python/validate_2d.py` (E0.5: native depo projection vs SP frames, the sec 6.3 timing self-check) | done 2026-10-21 |
| T8 | `python/validate_3d.py` (E0.6: per-blob from `blobtruth`, per-point from the `truth` cloud, count-based matching) | done 2026-10-21 |
| T9 | `python/qa_report.py` (sec 8 items; 3-5 and 8 placeholders) | done 2026-10-21 |

## Python

`python/truth_table.py` (file layout, column/process lists), `dump_truth_table.py` (Stage A dumper, gallery), `validate_2d.py`, `validate_3d.py`, `qa_report.py`; every script prints its usage with `-h`. The R1 campaign drivers that call them live in `wcp-porting-img/pr_ml/r1/`.

## Rules

- Depends on `WireCellUtil`, `WireCellIface`, `WireCellAux` only; never on another plugin.
  Pure algorithms are header-only so `clus` visitors can include them without a link edge.
- Every behaviour change in production packages is behind a default-OFF knob; the only
  exception is the `Aux::fill` child-index defect fix (doc 19 sec 3.3, owner decision O2),
  which changes the `child` column of depo files written with priors.
- `Schema.h` is the contract: names, orders, sentinels (`tid_ghost = -1`, `tid_shared = -2`),
  class codes, process codes. Any change bumps `schema_version` (2 since T4: the process code
  table). `python/truth_table.py` mirrors the column and process lists.

## Depo file format (after T2)

`DepoFileSink` writes per depo set `depo_data_<ident>.npy` (N,7) float32
`t, q, x, y, z, long, tran` and `depo_info_<ident>.npy` (N,4) int32
`id, pdg, gen, child`; with `truth_sidecar: true` also `depo_truth_<ident>.npy`
(N,3) float64 `energy, orig_id, reserved`. Rows are pushed depth-first: an active
depo (gen 0) is followed by its prior (gen 1), that one's prior (gen 2), etc.; the
`child` of a gen>0 row is the row of the depo whose `prior()` it is. Files written
before 2026-10-21 have `child == own row` on every prior row; `DepoFileSource` and
`NumpyDepoSetLoader` detect that and use "child = previous row".
