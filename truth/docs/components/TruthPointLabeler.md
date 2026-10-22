# TruthPointLabeler

Adds the per-point `truth` local point cloud to a live point tree. Joins the tensor set of `PointTreeBuilding`, the `ICluster` it was sampled from and the drifted depos; runs `Truth::depo_blob_assoc` (blob labels identical to `BlobDepoLabel`), the `PFMerger` on the truth table with the measured per-track charge, and `DepoPolyline` per track; writes the cloud on every live blob node in the order `PointTreeBuilding::sample_live` created them, cross-checking each node against its blob through the `scalar` keys (`wpid`, `slice_index_min`, the wire index bounds). The `3d` and `scalar` clouds are untouched; the output set carries every other input tensor unchanged plus `blobtruth`, `blobtruth_sum`, `truth_nu` and `truth_pf`.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `TruthPointLabeler` |
| Concrete class | `WireCell::Truth::TruthPointLabeler` |
| Node category | join |
| Primary interface | `ITruthPointLabeler` |
| Input type(s) | `ITensorSet`, `ICluster`, `IDepoSet` |
| Output type(s) | `ITensorSet` |
| Configurable | yes |

## Configuration Parameters

The association knobs of `BlobDepoLabel` (`time_offset`, `nsigma`, `speed` (required), `pindex`, `id_from_prior`, `extra_smear_long`, `extra_smear_tran`) plus:

| Parameter | Description |
| --- | --- |
| `inpath` | datapath pattern of the live tree, `%d` = set ident (default `pointtrees/%d/live`) |
| `anode`, `face` | the AnodePlane and face of the tree, as `PointTreeBuilding` (anode required) |
| `detector_volumes` | DetectorVolumes giving `time_offset`, `drift_speed`, `tick` for the reconstructed frame, as `PointTreeBuilding`; empty = use `reco_time_offset`, `drift_speed`, `tick` |
| `truth_table` | truth file pattern, `%d` = cluster ident; absent or '' = particle table from the depos (`TruthTable::from_depos`) |
| `truth_prefix`, `blobtruth_path` | datapath prefixes of the truth tensors (default `truth/%d`, `clusters/%d`) |
| `id_mode` | `auto` (index mode when the table has `sed_ids`), `track`, `index` |
| `pfmerger` | `PFMergerParams`: `pf_ke_min`, `q_min`, `delta_len`, `compton_q`, `track_len_min`, `michel_own`, `photon_prop`, `void_radius` |
| `refine_by_nearest_depo` | per-point track from the nearest contributing depo (default true); false = the blob's largest fraction, `tid_shared` below `shared_frac` |
| `d_cap`, `end_window`, `kink_window`, `kink_angle`, `heat_sigma`, `shared_frac` | `PointLabelParams` (defaults 25 mm, 10 mm, 3 cm, 15 deg, 1 cm, 0.5) |
| `labels_meta` | file to write the labels_meta JSON to; the same JSON is in the output set metadata (default none) |
| `label_dead` | reserved: dead nodes are not labelled yet (default false; true warns) |
| `polyline_drift_volume_only` | build the trajectories from the depos of the face's drift volume (x between its anode and cathode planes, any y / z, so a track keeps its part over a neighbouring anode of the same drift side; not behind the response plane); default false = every depo of the event, so a track crossing the cathode gets a mirrored copy in this face's frame (pr_ml doc 23 sec 7.4) |
| `polyline_per_particle` | one polyline per Geant4 particle of a key (index mode: the `sed_ids` orig_trackid; track mode: the depo id before folding); `d_traj`/`dir`/`s_t`/`dedx_true`/`end` from the nearest particle, end/kink heat from all of them; default false = one time-ordered polyline per key, a zigzag across a key's folded EM daughters (pr_ml doc 23 sec 7.4) |

## The `truth` point cloud

Keys in `Truth::truth_pc_keys()` order, same major size as `3d`:

| key | type | meaning |
| --- | --- | --- |
| `tid`, `orig_tid` | int | G4 track key of the point (nearest contributing depo when refined), its pre-merge id through `sed_ids`; `-1` ghost, `-2` shared (unrefined) |
| `pid_frac`, `true_q_frac` | double | blob charge fraction of `tid`; blob true charge over solved charge |
| `instance_valid` | int | 0 for ghost, shared and LE points |
| `pdg`, `cls`, `iid`, `pf_node` | int | from the table and `PFMerger` (`cls`: 0 track, 1 shower, 2 michel, 3 delta, 4 le, -1 unknown; `iid` = interaction row) |
| `d_traj`, `dir_x/y/z`, `s_t`, `end` | double / int | distance to the track's depo polyline (capped at `d_cap`), its tangent, arc fraction, 1/2 within `end_window` of the start/end (with `polyline_per_particle`: of the nearest particle's polyline) |
| `v_heat_vtx/end/kink/conv` | double | Gaussian heat (`heat_sigma`) around the surviving primaries' starts, this track's end and kinks, and the `conv` electrons' starts |
| `dedx_true` | double | the key depo's `energy()` per cm along the polyline; meaningful only with the depo truth sidecar (`DepoFileSink.truth_sidecar`), else 1 per depo step |

Gate: PDVD held-out evt246 anode 0, clustering with and without the node (`wcp-porting-img/pr_ml/docs/19_native_truth_plumbing.md` sec 14.3): every Bee zip identical, the pctree dump identical on every pre-existing datapath, the `truth` arrays added.
