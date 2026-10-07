# `match_isFC` in `T_tagger`: implementation note for the next SBND production

For Haiwang, to prepare a new production round. Written 2026-10-07 at toolkit `469ae3be`
(master merged into `apply-pointcloud`; the change itself is `0bb05b3c`). Full study with the
gate records: `wcp-porting-validation` `sbnd/sbnd_xin/docs/133_match-isfc-per-candidate-in-root.md`.

## What changed

`tracking-pr.root` now has one more branch:

| tree | branch | type | meaning |
|---|---|---|---|
| `T_tagger` | `match_isFC` | `Float_t`, 0 or 1 | the neutrino candidate's main cluster is fully contained in the PR fiducial volume |

- One value per `T_tagger` row, so **one per neutrino candidate**. Join to the candidate with
  `nu_index` / `cluster_id` on the same row.
- It is the value the two top BDTs already read (numu input 70, nue input 0). Nothing about how it
  is computed changed; it was computed and used but never written.
- Nothing else in the file changes relative to the build just before this change: every other
  branch of every tree is identical before and after.
- Against **older production files** more will differ, none of it from this change. Compared with
  the 2026-09-25 production outputs (`work-*-m0925pr`), the new samples also carry the
  particle-link additions (`T_segment`, `T_kine.kine_particle_*`, `T_kine.kine_main_vertex_id`,
  `T_tagger.act_n_seg*`, `act_in_enu`), and one NCpi0 event (359980) has a different candidate
  cluster (`vertex_moved_cluster` 0 to 1, `numu_score` -1.1088 to -1.1177). All other common
  physics branches are identical on the 48 + 19 events. The build just before the merge
  reproduces the new NCpi0 sample exactly, so these come from production changes made between
  09-25 and this change.

## What you need to do in the configuration

**Nothing.** The branch is booked by `UbooneTaggerOutputVisitor` whenever its `nu_provenance` key
is set. SBND production already sets it:

- `cfg/pgrapher/experiment/sbnd/clus.jsonnet`, `pr()` default `root_nu_record=true`, passed to
  `tagger_output(nu_provenance=root_nu_record)`.
- The 2-step job (`sbnd-pr-stage.jsonnet`, `wct-pr.jsonnet`) and the 1-step `clus_pr` both go
  through `pr()`, so both get the branch with the defaults.
- Checked on the compiled 2-step job: `wcsonnet --tla-str input=qlpctree.tar.gz --tla-str
  reality=data pgrapher/experiment/sbnd/wct-pr.jsonnet` (and `reality=sim`) gives
  `UbooneTaggerOutputVisitor` with `nu_provenance: true`. The event runs below used the
  per-event job `wct-pr-perevt.jsonnet`; the 2-step job was compiled, not run.

To switch it off (pre-change schema): `root_nu_record=false`. That also drops the other
selection-provenance branches, as before.

## What to rebuild

Only `libWireCellRoot.so` changed for this branch (`root/src/UbooneTaggerOutputVisitor.cxx`, one
`Branch` call). A normal build of the merged toolkit is enough; no data file, weight file or
fcl change.

## How to check a new production file

```python
import uproot
t = uproot.open("tracking-pr.root")["T_tagger"]
a = t.arrays(["nu_index", "cluster_id", "match_isFC", "numu_score", "nue_score"], library="np")
assert set(a["match_isFC"]) <= {0.0, 1.0}
```

Stronger check, the one used here: re-infer the BDT scores from the file with
`bdt_roundtrip.py` (`/nfs/data/1/xning/wirecell-working/SBND/bdt_training/parameters_check/`).
The value its 0/1 scan infers must equal the stored branch on every candidate.

## Validation done

| check | sample | result |
|---|---|---|
| key on, before vs after the change, geometric vertex | 48 nueCC + 17 events (13 with two candidates) | only `match_isFC` added; all other branches identical |
| key off, before vs after | 12 nueCC | identical, no branch added |
| before vs after the master merge (PR 536), geometric vertex | 48 nueCC | identical except `Trun.toolkit_git` |
| before vs after the master merge, production defaults (DL vertex on) | 19 NCpi0 | identical |
| stored value vs BDT round-trip scan | 74 + 68 candidates | equal on all |
| clean cmake build of `469ae3be` with ROOT, aggregate `wcdoctest` | | 1006 / 1006 cases pass |
| clean cmake build of `469ae3be` in release mode (`-Werror`), all libraries incl. root and pytorch | | builds, 0 errors |

Reference samples made with the merged build at production defaults (DL vertex on), on wcgpu1:

- `/nfs/data/1/xqian/toolkit-dev/wcp-porting-img/sbnd/sbnd_xin/work-nuecc48-d133pr/pr_evt<ID>/tracking-pr.root`
  (48 nueCC events, 48 candidates, `match_isFC` = 1 on 40)
- `/nfs/data/1/xqian/toolkit-dev/wcp-porting-img/sbnd/sbnd_xin/work-ncpi0-d133pr/pr_evt<ID>/tracking-pr.root`
  (19 NCpi0 events, 20 candidates, `match_isFC` = 1 on 17)

Both are the PR stage only, run on the existing `m0925` charge-light matching products.

## Things to know

- **Old files do not have the branch.** Anything produced before `0bb05b3c` needs the PR stage
  rerun to get it. For old files the 0/1 scan of `bdt_roundtrip.py` recovers the same value.
- **`T_cluster.fc` and `T_tagger.act_fc` are not `match_isFC`.** They carry the earlier
  fully-contained tagger flag. They agree with `match_isFC` on 71 of 74 candidates, not all.
  For retraining use `match_isFC`.
- **To look up a candidate in `T_cluster`, join on `T_tagger.cluster_id`.** `T_cluster.is_main`
  marks every bundle's main cluster (9 to 27 rows per event), not the candidate's.
- **Other detectors.** The same writer key is on in the FD-HD job (`dune10kt-1x2x6/pr.jsonnet`) and
  the ICARUS job, so their `T_tagger` gains the branch too (not run). PDHD and PDVD do not set the
  key.
- **cmake notes** (not from this change): the `mcs` unit tests need
  `MCS_TEST_DATA=<tree>/mcs/test/data`, a directory that is git-ignored; and on wcgpu1 the
  programs that link libtorch fail to link against the system BLAS, so the aggregate test binary
  was built with `-DWITH_LIBTORCH=no`.
