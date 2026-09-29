# The SBND 2-step workflow: `wcls-img-clus-matching.jsonnet` + `wct-pr.jsonnet`

The LArSoft 1-step chain (`wcls-img-clus-matching-xin-lib.jsonnet`, see [`wcls-img-clus-matching-xin-chain.md`](wcls-img-clus-matching-xin-chain.md)) is split into two jobs (ai-helper issue 33):

1. **Step 1, LArSoft (`lar`):** imaging, clustering, charge-light (Q/L) matching and all-APA clustering. It writes the matched bundles, with the light, CTPC, metadata and MC truth, to an ITensorSet tar on disk.
2. **Step 2, standalone (`wire-cell`, no art):** pattern recognition on that tar. This covers taggers, track/shower separation, trajectory fitting, vertexing, PID, particle flow, energy reconstruction and the numu/nue scores. It can be re-run as often as needed without redoing step 1.

Validated on MC-9, NCpi0-19 and nueCC-48: the output is identical to the 1-step chain, in `tracking-pr.root` and in every Bee layer (issue 33, log b and c).

Both steps build their PR stage from one definition, `sbnd-pr-stage.jsonnet`, and the 1-step chain uses the same file. The steps cannot drift apart.

## 1. Workflow

```mermaid
flowchart LR
  reco1[("reco1 artROOT<br/>recob::Wire, recob::OpHit<br/>(+ MCTruth, MCParticle,<br/>SimEnergyDeposit on MC)")]
  subgraph S1 ["Step 1: lar -c wcls-img-clus-matching[-data].fcl<br/>(wcls-img-clus-matching.jsonnet = lib stage='ql')"]
    s1["imaging, per-APA clustering,<br/>hit flashes, Q/L matching,<br/>all-APA clustering,<br/>labeler_truth, truth attacher"]
  end
  tar[("qlpctree.tar.gz<br/>ITensorSet per event<br/>(all events of the lar job)")]
  bee1[("mabc.zip<br/>img / clustering / op /<br/>dead-area / truth Bee sets")]
  h5[("nugraph.h5")]
  subgraph S2 ["Step 2: wire-cell -c pgrapher/experiment/sbnd/wct-pr.jsonnet<br/>--tla-str input=qlpctree.tar.gz --tla-str reality=data or sim"]
    s2["PR MABC clus_pr<br/>16 visitors<br/>(sbnd-pr-stage.jsonnet)"]
  end
  trk[("pr_evt{E}/tracking-pr.root<br/>Trun, T_kine, T_tagger, T_cluster,<br/>T_rec_charge, T_bundle, T_flash, ...<br/>+ T_truth_nu, T_truth_pf (MC)")]
  bee2[("mabc-pr.zip<br/>clustering-pr / track_fit /<br/>shower_track / vertices /<br/>mc / tagger_* Bee sets")]

  reco1 --> s1
  s1 --> tar
  s1 --> bee1
  s1 --> h5
  tar --> s2
  s2 --> trk
  s2 --> bee2
  tar -. "re-run PR as often as needed<br/>(run-2step.pbs SKIP1=1)" .-> s2

  classDef art fill:#d6eaff,stroke:#1f5fa8
  classDef wct fill:#fde2b3,stroke:#b36b00
  classDef file fill:#eeeeee,stroke:#666
  class s1 art
  class s2 wct
  class reco1,tar,bee1,h5,trk,bee2 file
```

The two Bee zips together hold every layer the 1-step chain's single `mabc.zip` holds. `bee-split.py` in issue 33 merges them per event, for comparison or upload.

| | step 1 | step 2 |
|---|---|---|
| **job** | `sbnd/wcls-img-clus-matching.fcl` (MC), `-data.fcl` (data), in wcp-porting-validation | `cfg/pgrapher/experiment/sbnd/wct-pr.jsonnet` |
| **config** | `wcls-img-clus-matching.jsonnet` = `wcls-img-clus-matching-xin-lib.jsonnet(flash_source='hits', xtpc_sc1_light_gate=true, xtpc_sc1_overpred_max=2.9, stage='ql')` | TLAs `input`, `reality`, `output_dir` (`.`), `evt_subdir` (`pr_evt%1%`), `enable_tracking_root` (true), `bee_outname` (`mabc-pr.zip`), `pr_tensor_outname` (`''` = no post-PR tar) |
| **events** | any number per `lar` job, all in one tar | every tensor set in the tar, one `wire-cell` process |
| **run** | `lar -n <N> --nskip <k> -c wcls-img-clus-matching-data.fcl -s <reco1.root> --no-output` | `mkdir pr_evt<E>` for each event in the tar, then `wire-cell -c pgrapher/experiment/sbnd/wct-pr.jsonnet --tla-str input=qlpctree.tar.gz --tla-str reality=data` |
| **cost (Aurora, 2 cores)** | about 20–40 s per event | about 5–25 s per event |

## 2. Step 1: `wcls-img-clus-matching.jsonnet`

This is the 1-step graph up to and including the truth attacher: 106 nodes and 118 edges. The per-APA imaging is expanded in `wcls-img-clus-matching-xin-chain.md` section 2. The PR tail (`clus_pr`, `labeler_tagger`, `TensorFileSink:clus_all_apa`, 34 components) is replaced by one `TensorFileSink`.

```mermaid
flowchart TB
  sigs("wclsCookedFrameSource:sigs<br/>recob::Wire")
  fan["FrameFanout:sig_fanout"]
  subgraph apa0 ["APA0"]
    img0["imaging<br/>(img.jsonnet multi-3view, full_deghost)"]
    ptb0["PointTreeBuilding:apa0-0<br/>(3d, CTPC, dead winds)"]
    rse0("wclsTruthInformationAttacher:rse_apa0<br/>truth: false")
    mabc0[["MABC apa0-0<br/>16 visitors"]]
    fa0["FlashTensorToOpticalPCs:flash_attach_apa0"]
  end
  subgraph apa1 ["APA1"]
    img1["imaging (same as APA0)"]
    ptb1["PointTreeBuilding:apa1-0"]
    rse1("wclsTruthInformationAttacher:rse_apa1<br/>truth: false")
    mabc1[["MABC apa1-0<br/>16 visitors"]]
    fa1["FlashTensorToOpticalPCs:flash_attach_apa1"]
  end
  subgraph light ["light (hit-rebuilt flashes)"]
    oh0("wclsOpHitSource:tpc0") --> ff0["SBNDOpFlashFinder:tpc0"]
    oh1("wclsOpHitSource:tpc1") --> ff1["SBNDOpFlashFinder:tpc1"]
  end
  qlm["QLMatching:matching_joint<br/>(scenario-1 light gate)"]
  rseall("wclsTruthInformationAttacher:rse_all_apa<br/>truth: false")
  mabcall[["MABC clus_all_apa<br/>10 visitors"]]
  lt("wclsTensorSetLabeler:labeler_truth<br/>blob trackid, nu_* metadata,<br/>bee_pf_truth, truth/sed Bee, nugraph.h5")
  truth("wclsTruthInformationAttacher:truth<br/>truth: true<br/>+ truth/{E}/nu, truth/{E}/pf")
  sink["TensorFileSink:ql_pctree<br/>qlpctree.tar.gz<br/>prefix clustering_, real tensors"]
  bee[("BeeSink:mabc_shared<br/>mabc.zip")]

  sigs --> fan
  fan -->|0| img0
  fan -->|1| img1
  img0 -->|live| ptb0
  img0 -->|dead| ptb0
  img1 -->|live| ptb1
  img1 -->|dead| ptb1
  ptb0 --> rse0 --> mabc0 -->|clusters| fa0
  ptb1 --> rse1 --> mabc1 -->|clusters| fa1
  ff0 -->|opflash| fa0
  ff1 -->|opflash| fa1
  fa0 -->|apa0| qlm
  fa1 -->|apa1| qlm
  qlm --> rseall --> mabcall --> lt --> truth --> sink
  mabc0 -.-> bee
  mabc1 -.-> bee
  mabcall -.-> bee
  lt -.-> bee

  classDef mabc fill:#fde2b3,stroke:#b36b00,stroke-width:2px
  classDef art fill:#d6eaff,stroke:#1f5fa8
  classDef out fill:#eeeeee,stroke:#666
  class mabc0,mabc1,mabcall mabc
  class sigs,oh0,oh1,rse0,rse1,rseall,lt,truth art
  class sink out
```

- **Colours:** blue marks larwirecell components. These are art-event visitors, so each must be listed in the fcl `inputers`. Orange marks MABC nodes.
- **Visitor lists:** the three MABCs' EnsembleVisitor lists are the 1-step chain's, in `wcls-img-clus-matching-xin-chain.md` section 3.
- **`wclsTruthInformationAttacher`** replaces `wclsTensorSetMetadataAttacher`. With `truth: false` it only stamps `runNo` / `subRunNo` / `eventNo` into the set metadata, where each downstream MABC reads them (`rse_from_metadata`).
- **The `truth` instance** also appends the two MC truth tables described in section 3. On data it stamps the RSE only.

## 3. What travels in `qlpctree.tar.gz`

This is one ITensorSet per event, keyed by the set ident (the art event number). As measured on an MC tar:

| datapath | datatype | content | needed by |
|---|---|---|---|
| `pointtrees/<E>/live` | `pctree` (+159 `pcarray`, 23 `pcdataset`) | the matched live clusters, listed below | the whole PR stage |
| `pointtrees/<E>/dead` | `pctree` (+26 `pcarray`) | the dead (masked) grouping | `MakeFiducialUtils` |
| `truth/<E>/nu` | `truth_nu` | one row per neutrino (MC); columns in the tensor metadata | `T_truth_nu` |
| `truth/<E>/pf` | `truth_pf` | one row per particle of the truth particle flow (MC) | `T_truth_pf` |
| `truthtracks/<E>` | `truth_per_track` | the labeler's beam-nu primaries table (MC) | not read |

The live clusters carry:
- blob `3d` points, with `x_t0cor`, `y_cor`, `z_cor`, the per-plane charge, and `trackid` on MC;
- `ctpc_a{0,1}f0p{U,V,W}`, the 2D charge used for the Steiner graph and the trajectory fit;
- `dead_winds_*` and `dead_gap_*`;
- the light: `flash`, `light`, `flashlight`, `opflash`;
- the provenance in `perblob`: `real_cluster_*` and `assoc_cluster_*`;
- the flags and `matched_flash_gid` in `cluster_scalar`.

**Set metadata:** `runNo`, `subRunNo`, `eventNo`; the labeler's `n_nu` and `nu_*` arrays; and on MC `bee_pf_truth`, the truth particle tree the PR MABC merges into Bee `mc`.

**`truth_nu` columns:** `nu_idx, pdg, ccnc, mode, int_type, flavor, E, vtx_x, vtx_y, vtx_z, t, edep`.

**`truth_pf` columns:** `nu_row, trackid, parent_trackid, mother_trackid, pdg, process, E, KE, start_x/y/z/t, end_x/y/z/t, start_px/py/pz`. The rows are beam-neutrino-derived particles with KE > 10 MeV, each attached to its nearest kept ancestor, the same selection as the Bee `mc` tree. Units are cm, ns and GeV.

## 4. Step 2: `wct-pr.jsonnet`

Standalone: 3 nodes, 2 edges and 44 configured components.

```mermaid
flowchart LR
  src["TensorFileSource:ql_pctree<br/>inname = input (qlpctree.tar.gz)<br/>prefix clustering_"]
  mabcpr[["MultiAlgBlobClustering:clus_pr<br/>16 EnsembleVisitors (table below)<br/>rse_from_metadata, event_from_ident,<br/>reset_shower_ids_per_event,<br/>aux_datatypes = truth_nu, truth_pf"]]
  sink["TensorFileSink:pr<br/>dump_mode (no file)<br/>or pr_tensor_outname"]
  bee[("BeeSink:mabc_pr<br/>mabc-pr.zip")]
  trk[("pr_evt{E}/tracking-pr.root")]
  src --> mabcpr --> sink
  mabcpr -. "PR Bee layers,<br/>TaggerBeeVisitor" .-> bee
  mabcpr -. "SbndPrMagnifyTrackingVisitor,<br/>UbooneTaggerOutputVisitor" .-> trk
  classDef mabc fill:#fde2b3,stroke:#b36b00,stroke-width:2px
  classDef out fill:#eeeeee,stroke:#666
  class mabcpr mabc
  class bee,trk out
```

- **Multi-event:** one process runs every event of the tar.
  - The per-event outputs go to `output_dir/pr_evt<E>/`, where `<E>` is the set ident. The runner must create those directories.
  - `reset_shower_ids_per_event` restarts the shower-id counter at each event. Each event's ids then equal a one-event process's, which is required for identity with the 1-step chain.
- **RSE:** taken from the step-1 metadata, which MABC ranks above the ident. `event_from_ident` is only the required partner of `evt_subdir`.
- **Truth:** the MABC publishes the input tensors of datatype `truth_nu` / `truth_pf` on the Ensemble (`Ensemble::aux_tensor`) and forwards them to its output. `SbndPrMagnifyTrackingVisitor` writes them as `T_truth_nu` / `T_truth_pf`: one entry per row, plus `runNo`, `subRunNo` and `eventNo`, with identifier columns as `Int_t`. On data they are absent, so no truth trees are written.
- **Bee sink:** `BeeSink:mabc_pr` is listed explicitly in the job's config, because `pr()` names its sink without listing it in its `uses`.

### EnsembleVisitors of `MultiAlgBlobClustering:clus_pr`

Visitors 1–15 are exactly the 1-step chain's PR MABC (`sbnd-pr-stage.jsonnet`, `pr()` defaults = the SBND production operating point). Visitor 16 is appended in step 2 only.

| # | EnsembleVisitor (type:name) | what it does |
|---|---|---|
| 1 | `ClusteringSwitchScope:pr` | switch to T0Correction coords |
| 2 | `ClusteringUnmergeBundle:pr` | split flash-merged bundles (real_cluster_id provenance, from step 1) |
| 3 | `ClusteringUnmergeBundle:prassoc` | split isolated groupings (assoc_cluster_id / assoc_cluster_main, from step 1) |
| 4 | `CreateSteinerGraph:pr` | Steiner graph (uses the CTPC) |
| 5 | `MakeFiducialUtils:pr` | fiducial utils from the live + dead groupings |
| 6 | `TaggerCheckTGM:pr` | through-going muon tagger |
| 7 | `TaggerCheckSTM:pr` | stopping muon tagger |
| 8 | `TaggerCheckFC:pr` | fully-contained check |
| 9 | `ClusteringProtectBundle:pr` | Protect_Over_Clustering, after the cosmic taggers |
| 10 | `CreateSteinerGraph:prrefresh` | rebuild the Steiner products the split purged |
| 11 | `TaggerCheckNeutrino:pr` | neutrino tagger: vertex, trajectory fit, track/shower, PID, particle flow, energy |
| 12 | `UbooneNumuBDTScorer:pr` | numu score (XGBoost `numu_scalars_scores_0923.xml`) |
| 13 | `UbooneNueBDTScorer:pr` | nue score (XGBoost `XGB_nue_seed2_0923.xml`); the Bee `mc` particle flow is dumped after it |
| 14 | `SbndPrMagnifyTrackingVisitor:pr` | `pr_evt<E>/tracking-pr.root` (RECREATE), incl. `T_truth_nu` / `T_truth_pf` on MC |
| 15 | `UbooneTaggerOutputVisitor:pr` | adds `T_tagger` / `T_kine` (UPDATE) |
| 16 | `TaggerBeeVisitor:pr` | **step 2 only**: the `tagger_stm/tgm/fc/lm` Bee sets (in the 1-step, the art-side `labeler_tagger` writes them) |

## 5. The MC truth, from art to `tracking-pr.root`

```mermaid
flowchart LR
  art[("art::Event<br/>generator MCTruth<br/>largeant MCParticle + Assns<br/>ionandscint:priorSCE SimEnergyDeposit")]
  lab("labeler_truth<br/>(wclsTensorSetLabeler)")
  tia("wclsTruthInformationAttacher:truth")
  md["set metadata<br/>nu_* arrays, bee_pf_truth"]
  tnu["tensor truth/{E}/nu<br/>(truth_nu)"]
  tpf["tensor truth/{E}/pf<br/>(truth_pf)"]
  tar[("qlpctree.tar.gz")]
  ens["PR MABC:<br/>Ensemble::aux_tensor"]
  trk[("tracking-pr.root<br/>T_truth_nu, T_truth_pf")]
  beemc[("Bee mc:<br/>reco nu + truth tree")]
  art -->|visit| lab
  art -->|visit| tia
  lab --> md
  tia --> tnu
  tia --> tpf
  md --> tar
  tnu --> tar
  tpf --> tar
  tar -->|step 2| ens
  ens --> trk
  md -. "merge_metadata_key<br/>bee_pf_truth" .-> beemc
```

The labeler and the attacher read the same art products with the same particle selection. The G4 process codes come from one shared header, `larwirecell/aiml/G4ProcessCode.h`. issue 33's `check-truth.py` confirms that the two agree event by event.

## 6. Files

| file | repo | role |
|---|---|---|
| `cfg/pgrapher/experiment/sbnd/wcls-img-clus-matching-xin-lib.jsonnet` | toolkit | the chain as a function; `stage='1step'` (full 1-step) or `'ql'` (step 1) |
| `cfg/pgrapher/experiment/sbnd/wcls-img-clus-matching.jsonnet` | toolkit | step 1 top-level job |
| `cfg/pgrapher/experiment/sbnd/wct-pr.jsonnet` | toolkit | step 2 top-level job |
| `cfg/pgrapher/experiment/sbnd/sbnd-pr-stage.jsonnet` | toolkit | the PR stage, shared by the 1-step and step 2 |
| `cfg/pgrapher/experiment/sbnd/clus.jsonnet` | toolkit | `pr()`, including its `tagger_bee` entry |
| `clus/src/TaggerBeeVisitor.cxx` | toolkit | tagger Bee sets from inside the PR MABC |
| `clus/inc/WireCellClus/Facade_Ensemble.h`, `MultiAlgBlobClustering.*` | toolkit | auxiliary (truth) tensors |
| `root/src/SbndPrMagnifyTrackingVisitor.cxx` | toolkit | `T_truth_nu` / `T_truth_pf` |
| `larwirecell/Components/TruthInformationAttacher.*` | larwirecell | `wclsTruthInformationAttacher` |
| `larwirecell/aiml/G4ProcessCode.h`, `TensorSetLabeler.*` | larwirecell | shared process codes; per-event RNG re-seed, which makes multi-event step-1 jobs reproducible |
| `sbnd/wcls-img-clus-matching{,-data}.fcl`, `sbnd/wcls-img-clus-matching.jsonnet` | wcp-porting-validation | step 1 fcls and the re-export |

The harness and validation scripts are in `issues/33-sbnd-1step-to-2step/scripts/` of wire-cell-toolkit-ai-helper:
- `run-2step.pbs`: both steps and the comparison; `SKIP1=1,RUN_DIR=` re-runs only PR;
- `gate-2step-cfg.sh`: the configuration gates;
- `check-truth.py`, `bee-split.py`, `bee-diff.py`, `evt-table-2step.py`.
