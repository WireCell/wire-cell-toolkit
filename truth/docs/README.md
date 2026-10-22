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
| WP0.2 G | `clus/src/SteinerTruthLabeler.cxx` (clus visitor: `steiner_truth` / `steiner_edges` cluster clouds from the `truth` cloud), `python/steiner_truth.py` (edge classes, flow, vertex flags, skeleton metrics; bands 15 / 40 mm) | started 2026-10-22 (pr_ml doc 23); exact particle search in the reference 2026-10-22 (doc 23 sec 12) |
| WP0.2 R | `clus/src/PRTruthExport.cxx` (clus visitor: `pr_segments` / `pr_seg_points`, the PR chain's segment -> cluster point assignment), `python/pr_truth.py` (segment / particle / interaction matching to PF nodes, SPINE / Pandora / PQ, direction, vertex, PID; truth beside every fit point) | FD-HD, SBND 2026-10-22 (pr_ml doc 23 sec 11, 13) |
| T6 fixes | `TruthPointLabeler` knobs `polyline_drift_volume_only`, `polyline_per_particle` (C++ default false; jsonnet `truth_polyline_drift_volume_only` / `truth_polyline_per_particle` on every clustering job, default true there -- truth-label jobs only), `build_particle_polylines` | done 2026-10-22 (pr_ml doc 23 sec 10) |

## Python

`python/truth_table.py` (file layout, column/process lists), `dump_truth_table.py` (Stage A dumper, gallery), `validate_2d.py`, `validate_3d.py`, `qa_report.py`, `xres_check.py` (doc 21), `parity_larsoft.py` (doc 22), `steiner_truth.py`, `pr_truth.py` (doc 23); every script prints its usage with `-h`. The R1 campaign drivers that call them live in `wcp-porting-img/pr_ml/r1/`.

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
