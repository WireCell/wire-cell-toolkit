# NeutrinoID → Toolkit Function Mapping

Mapping from the prototype `WCPPID::NeutrinoID` (in `prototype_pid/src/`) to the toolkit
`WireCell::Clus::PR::PatternAlgorithms` class (split across `clus/src/Neutrino*.cxx`).

---

## Structural Overview

| Concept | Prototype | Toolkit |
|---|---|---|
| Main class | `WCPPID::NeutrinoID` (monolithic ~5000-line god-class) | `WireCell::Clus::PR::PatternAlgorithms` (split across 8 source files) |
| Entry point | `NeutrinoID` constructor (runs everything) | `TaggerCheckNeutrino::visit()` (calls `PatternAlgorithms` methods) |
| Cluster | `WCPPID::PR3DCluster*` | `WireCell::Clus::Facade::Cluster&` + `WireCell::Clus::PR::Graph&` |
| Vertex | `WCPPID::ProtoVertex*` (raw ptr) | `WireCell::Clus::PR::VertexPtr` (`shared_ptr<PRVertex>`) |
| Segment | `WCPPID::ProtoSegment*` (raw ptr) | `WireCell::Clus::PR::SegmentPtr` (`shared_ptr<PRSegment>`) |
| Shower | `WCPPID::WCShower*` (raw ptr) | `WireCell::Clus::PR::ShowerPtr` (`shared_ptr<PRShower>`) |
| Tracking/fit engine | `PR3DCluster` member functions + global state | `WireCell::Clus::TrackFitting` (separate class, passed by ref) |
| Detector geometry | `ToyFiducial`, `WCPSst::GeomDataSource` | `IDetectorVolumes::pointer dv` |
| Particle data (PID tables) | hardcoded / `TPCParams` | `ParticleDataSet::pointer particle_data` |
| Recombination model | hardcoded | `IRecombinationModel::pointer recomb_model` |
| Accumulator IDs | `acc_vertex_id`, `acc_segment_id` member vars | passed explicitly; vertex IDs managed by graph |

Prototype member variable `main_vertex` (set as side-effect) → toolkit functions take `VertexPtr& main_vertex` as an **in-out reference** parameter.

Prototype implicit cluster context (accessed via `map_vertex_segments`, etc. member maps) → toolkit passes `Graph& graph` + `Facade::Cluster& cluster` explicitly.

---

## Prototype Source Files and Their Toolkit Equivalents

### `NeutrinoID_proto_vertex.h` → `NeutrinoPatternBase.cxx`

Core pattern recognition: initial segment finding, track breaking, structure helpers.

| Prototype (`WCPPID::NeutrinoID::`) | Toolkit (`PatternAlgorithms::`) | Notes |
|---|---|---|
| `find_proto_vertex(cluster, flag_break_track, nrounds, [flag_back_search])` | `find_proto_vertex(graph, cluster, track_fitter, dv, flag_break_track, nrounds, flag_back_search)` | Signature extended with explicit `graph`, `track_fitter`, `dv` |
| `init_point_segment(cluster)` | `init_point_segment(graph, cluster, track_fitter, dv)` | Same logic |
| `calc_PCA_main_axis(points)` | `calc_PCA_main_axis(points)` | Identical logic |
| `break_segments(segments, cluster, dis_cut)` | `break_segments(graph, track_fitter, dv, remaining_segments, dis_cut)` | |
| `find_other_segments(cluster, flag_break, range, scale)` | `find_other_segments(graph, cluster, track_fitter, dv, flag_break, range, scale)` | → `NeutrinoOtherSegments.cxx` |
| `find_vertex_other_segment(cluster, seg, flag_fwd, wcp, ...)` | `find_vertex_other_segment(graph, cluster, seg, flag_fwd, wcp, track_fitter, dv)` | → `NeutrinoOtherSegments.cxx` |
| `modify_vertex_isochronous(vtx, v1, sg, v2, cluster)` | `modify_vertex_isochronous(graph, cluster, vtx, v1, sg, v2, track_fitter, dv)` | → `NeutrinoOtherSegments.cxx` |
| `modify_segment_isochronous(sg1, v1, sg, v2, cluster, ...)` | `modify_segment_isochronous(graph, cluster, sg1, v1, sg, v2, track_fitter, dv, ...)` | → `NeutrinoOtherSegments.cxx` |
| `check_end_point(cluster, path, flag_front, ...)` | `check_end_point(graph, cluster, path, flag_front, ...)` | |
| `del_proto_vertex(pv)` / `del_proto_segment(ps)` | Graph edge/vertex removal via boost graph API | No direct wrapper; inlined |
| `add_proto_connection(pv, ps, cluster)` | `Graph::add_edge(v1, v2, segment)` equivalents | Inlined |
| `clean_up_maps_vertices_segments(cluster)` | `clean_up_graph(graph, cluster)` | |
| `find_vertices(sg)` → pair of `ProtoVertex*` | Graph incident vertices lookup | No direct wrapper |
| `find_other_vertex(sg, vtx)` | `find_other_vertex` in graph context | Likely inlined or graph helper |
| `do_rough_path(cluster, p1, p2)` | `do_rough_path(cluster, p1, p2)` | |
| `do_rough_path_reg_pc(...)` | `do_rough_path_reg_pc(cluster, p1, p2, graph_name)` | |
| `proto_extend_point(cluster, p, dir, dir_other, flag)` | `proto_extend_point(cluster, p, dir, dir_other, flag)` | |
| `proto_break_tracks(cluster, first_wcp, curr_wcp, last_wcp, list1, list2, flag)` | `proto_break_tracks(cluster, first_wcp, curr_wcp, last_wcp, list1, list2, flag)` | |
| `merge_nearby_vertices(cluster, ...)` | `merge_nearby_vertices(graph, cluster, track_fitter, dv)` | |
| `merge_two_segments_into_one(sg1, vtx, sg2)` | `merge_two_segments_into_one(graph, sg1, vtx, sg2, dv)` | |
| `merge_vertex_into_another(vtx_from, vtx_to)` | `merge_vertex_into_another(graph, vtx_from, vtx_to, dv)` | |
| `create_segment_for_cluster(cluster, dv, path, dir)` | `create_segment_for_cluster(cluster, dv, path, dir)` | |
| `create_segment_from_vertices(graph, cluster, v1, v2, dv)` | `create_segment_from_vertices(graph, cluster, v1, v2, dv)` | |
| `replace_segment_and_vertex(...)` | `replace_segment_and_vertex(...)` (two overloads) | |
| `break_segment_into_two(...)` | `break_segment_into_two(...)` | |
| `vertex_get_dir(vertex, dis_cut)` | `vertex_get_dir(vertex, graph, dis_cut)` | |
| `vertex_segment_get_dir(vertex, segment, dis_cut)` | `vertex_segment_get_dir(vertex, segment, graph, dis_cut)` | |
| `print_segs_info(cluster_id, main_vertex)` | `print_segs_info(graph, cluster, vertex)` | |
| `init_first_segment(cluster, main_cluster, ...)` | `init_first_segment(graph, cluster, main_cluster, track_fitter, dv, flag_back_search)` | |
| `transfer_info_from_segment_to_cluster(cluster)` | `transfer_info_from_segment_to_cluster(graph, cluster)` | |
| `calc_dir_cluster(graph, cluster, point, dis_cut)` | `calc_dir_cluster(graph, cluster, point, dis_cut)` | |

---

### `NeutrinoID_examine_structure.h` → `NeutrinoStructureExaminer.cxx`

Structural clean-up of the PR graph (merging linear vertices, etc.).

| Prototype | Toolkit | Notes |
|---|---|---|
| `examine_structure(cluster)` | `examine_structure(graph, cluster, track_fitter, dv)` | calls `_1` and `_2` |
| `examine_structure_1(cluster)` | `examine_structure_1(graph, cluster, track_fitter, dv)` | |
| `examine_structure_2(cluster)` | `examine_structure_2(graph, cluster, track_fitter, dv)` | |
| `examine_structure_3(cluster)` | `examine_structure_3(graph, cluster, track_fitter, dv)` | |
| `examine_structure_4(vtx, cluster, flag_final)` | `examine_structure_4(vtx, flag_final, graph, cluster, track_fitter, dv)` | param order changed |
| `examine_segment(cluster)` | `examine_segment(graph, cluster, track_fitter, dv)` | |
| `crawl_segment(sg, v1, cluster)` | `crawl_segment(graph, cluster, seg, vertex, track_fitter, dv)` | |
| `examine_partial_identical_segments(cluster)` | `examine_partial_identical_segments(graph, cluster, track_fitter, dv)` | |
| `examine_vertices(cluster)` | `examine_vertices(graph, cluster, track_fitter, dv, main_vertex)` | |
| `examine_vertices_1(cluster)` | `examine_vertices_1(graph, cluster, track_fitter, dv, main_vertex)` | |
| `examine_vertices_1(v1, v2, ...)` (many params) | `examine_vertices_1p(graph, v1, v2, track_fitter, dv)` | renamed with `p` suffix |
| `examine_vertices_2(cluster)` | `examine_vertices_2(graph, cluster, track_fitter, dv, main_vertex)` | |
| `examine_vertices_4(cluster)` | `examine_vertices_4(graph, cluster, track_fitter, dv, main_vertex)` | |
| `examine_vertices_4(v1, v2)` | `examine_vertices_4p(graph, v1, v2, track_fitter, dv)` | renamed with `p` suffix |
| `examine_vertices_3()` (no args, uses member vars) | `examine_vertices_3(graph, cluster, pair_vertices, track_fitter, dv)` | explicit params |
| `get_local_extension(cluster, wcp)` | `get_local_extension(cluster, wcp)` | |

**Final structure examination** (called from `improve_vertex` path):

| Prototype | Toolkit | Notes |
|---|---|---|
| `examine_structure_final(cluster)` | `examine_structure_final(graph, main_vertex, cluster, track_fitter, dv)` | |
| `examine_structure_final_1(cluster)` | `examine_structure_final_1(graph, main_vertex, cluster, track_fitter, dv)` | |
| `examine_structure_final_1p(cluster)` | `examine_structure_final_1p(graph, main_vertex, cluster, track_fitter, dv)` | |
| `examine_structure_final_2(cluster)` | `examine_structure_final_2(graph, main_vertex, cluster, track_fitter, dv)` | |
| `examine_structure_final_3(cluster)` | `examine_structure_final_3(graph, main_vertex, cluster, track_fitter, dv)` | |

---

### `NeutrinoID_track_shower.h` → `NeutrinoTrackShowerSep.cxx`

Track/shower classification for all segments.

| Prototype | Toolkit | Notes |
|---|---|---|
| `separate_track_shower(cluster)` | `separate_track_shower(graph, cluster)` | |
| `separate_track_shower()` (no-arg, runs all clusters) | called per-cluster in toolkit; no global version | |
| `determine_direction(cluster)` | `determine_direction(graph, cluster, particle_data, recomb_model)` | PID tables now explicit params |
| `shower_determing_in_main_cluster(cluster)` | `shower_determining_in_main_cluster(graph, cluster, particle_data, recomb_model, dv)` | typo fixed in toolkit name |
| `improve_maps_one_in(cluster, flag_strong)` | `improve_maps_one_in(graph, cluster, particle_data, recomb_model, flag_strong)` | |
| `improve_maps_shower_in_track_out(cluster, flag_strong)` | `improve_maps_shower_in_track_out(graph, cluster, particle_data, recomb_model, flag_strong)` | |
| `improve_maps_no_dir_tracks(cluster_id)` | `improve_maps_no_dir_tracks(graph, cluster, particle_data, recomb_model)` | |
| `improve_maps_multiple_tracks_in(cluster_id)` | `improve_maps_multiple_tracks_in(graph, cluster, particle_data, recomb_model)` | |
| `fix_maps_multiple_tracks_in(cluster_id)` | `fix_maps_multiple_tracks_in(graph, cluster)` | |
| `fix_maps_shower_in_track_out(cluster_id)` | `fix_maps_shower_in_track_out(graph, cluster)` | |
| `examine_maps(cluster)` | `examine_maps(graph, cluster)` | |
| `examine_good_tracks(cluster_id)` | `examine_good_tracks(graph, cluster, particle_data)` | |
| `judge_no_dir_tracks_close_to_showers(cluster_id)` | `judge_no_dir_tracks_close_to_showers(graph, cluster, particle_data, dv)` | |
| `examine_all_showers(cluster_id)` | `examine_all_showers(graph, cluster, particle_data)` | |
| `calculate_num_daughter_showers(vtx, sg, flag)` | `calculate_num_daughter_showers(graph, vtx, sg, flag)` | |
| `calculate_num_daughter_tracks(vtx, sg, flag_count_shower, length_cut)` | `calculate_num_daughter_tracks(graph, vtx, sg, flag_count_shower, length_cut)` | BFS count of non-shower (or all) segments beyond `sg` from `vtx` |
| `find_cont_muon_segment_nue(sg, vtx, flag_ignore_dQ_dx)` | `find_cont_muon_segment_nue(graph, sg, vtx, flag_ignore_dQ_dx)` | Like `find_cont_muon_segment` but 30 cm angle threshold instead of 50 cm; used in `bad_reconstruction` |
| `set_default_shower_particle_info(cluster)` | `set_default_shower_particle_info(graph, cluster, particle_data, recomb_model)` | → `NeutrinoPatternBase.cxx` |
| `change_daughter_type(vtx, sg, type, mass)` | `change_daughter_type(graph, vtx, sg, type, mass, particle_data, recomb_model)` | → `NeutrinoVertexFinder.cxx` |

---

### `NeutrinoID_improve_vertex.h` → `NeutrinoVertexFinder.cxx`

Vertex fitting and refinement.

| Prototype | Toolkit | Notes |
|---|---|---|
| `fit_vertex(vtx, sg_set, cluster)` | `fit_vertex(cluster, vtx, main_vertex, sg_set, track_fitter, dv)` | |
| `improve_vertex(cluster, flag_search_vtx, flag_final)` | `improve_vertex(graph, cluster, main_vertex, vertices_in_long_muon, segments_in_long_muon, track_fitter, dv, particle_data, recomb_model, flag_search_vtx, flag_final)` | many new explicit params |
| `search_for_vertex_activities(vtx, sg_set, cluster, range)` | `search_for_vertex_activities(graph, vtx, segments_set, cluster, track_fitter, dv, range)` | |
| `eliminate_short_vertex_activities(cluster, existing_segs)` | `eliminate_short_vertex_activities(graph, cluster, main_vertex, existing_segs, track_fitter, dv)` | |
| `get_dir(vtx, sg, dis)` | `vertex_segment_get_dir(vtx, sg, graph, dis)` | renamed |

---

### `NeutrinoID_final_structure.h` → `NeutrinoVertexFinder.cxx`

Main vertex determination across all clusters.

| Prototype | Toolkit | Notes |
|---|---|---|
| `determine_main_vertex(cluster)` / `determine_main_vertex(cluster, false)` | `determine_main_vertex(graph, cluster, main_vertex, vertices_in_long_muon, segments_in_long_muon, track_fitter, dv, particle_data, recomb_model)` | `false` flag (other cluster) removed; toolkit treats all clusters identically |
| `determine_overall_main_vertex()` | `determine_overall_main_vertex(graph, map_cluster_main_vertices, main_cluster, other_clusters, vertices_in_long_muon, segments_in_long_muon, track_fitter, dv, particle_data, recomb_model, flag_dev_chain)` | returns `VertexPtr` instead of setting member var |
| `examine_direction(main_vertex)` | `examine_direction(graph, vertex, main_vertex, vertices_in_long_muon, segments_in_long_muon, particle_data, recomb_model, flag_final)` | |
| `examine_main_vertex_candidate(vtx)` | `examine_main_vertex_candidate(graph, vtx)` | |
| `compare_main_vertices(candidates)` | `compare_main_vertices(graph, cluster, candidates)` | |
| `compare_main_vertices_all_showers(cluster, candidates, ...)` | `compare_main_vertices_all_showers(graph, cluster, candidates, track_fitter, dv, particle_data, recomb_model)` | |
| `compare_main_vertices_global(candidates, cluster, ...)` | `compare_main_vertices_global(graph, candidates, cluster, track_fitter, dv)` | |
| `calc_conflict_maps(vtx)` | `calc_conflict_maps(graph, vtx)` | |
| `find_cont_muon_segment(sg, vtx, flag)` | `find_cont_muon_segment(graph, sg, vtx, flag)` | |
| `swap_main_cluster(cluster)` | `swap_main_cluster(new_main, old_main, other_clusters)` | explicit args |
| `examine_main_vertices(...)` | `examine_main_vertices(graph, map_cluster_main_vertices, main_cluster, other_clusters)` | |
| `check_switch_main_cluster(...)` | `check_switch_main_cluster(graph, map_cluster_main_vertices, main_cluster, other_clusters, track_fitter, dv)` | |
| `check_switch_main_cluster_2(...)` | `check_switch_main_cluster_2(graph, temp_vertex, max_length_cluster, main_cluster, other_clusters)` | |
| `examine_main_vertices_local(vertices, ...)` | `examine_main_vertices_local(graph, vertices, particle_data, recomb_model)` | |

---

### `NeutrinoID_DL.h` → `NeutrinoVertexFinder.cxx`

Deep-learning vertex selection.

| Prototype | Toolkit | Notes |
|---|---|---|
| `determine_overall_main_vertex_DL()` (no args, uses member `flag_dl_vtx`, `dl_vtx_cut`) | `determine_overall_main_vertex_DL(graph, map_cluster_main_vertices, main_cluster, other_clusters, vertices_in_long_muon, segments_in_long_muon, track_fitter, dv, particle_data, recomb_model, dl_weights, dl_vtx_cut, dQdx_scale, dQdx_offset)` | returns `bool` (did DL change the vertex?); same semantics |

---

### `NeutrinoID_deghost.h` → `NeutrinoDeghoster.cxx`

Ghost hit removal across all clusters.

| Prototype | Toolkit | Notes |
|---|---|---|
| `deghosting()` | `deghosting(graph, map_cluster_main_vertices, all_clusters, track_fitter, dv)` | |
| `deghost_clusters()` | `deghost_clusters(graph, all_clusters, track_fitter, dv)` | |
| `deghost_segments()` | `deghost_segments(graph, map_cluster_main_vertices, all_clusters, track_fitter, dv)` | |
| `order_clusters(ordered, map_id_segs, map_total_len)` | `order_clusters(graph, ordered, map_cluster_segs, map_total_len)` | |
| `order_segments(ordered, segs)` | `order_segments(ordered, segs)` | |

---

### `NeutrinoID_shower_clustering.h` + `NeutrinoID_em_shower.h` → `NeutrinoShowerClustering.cxx`

Shower building and pi0 identification.

| Prototype | Toolkit | Notes |
|---|---|---|
| `shower_clustering_with_nv()` (no args, uses member vars) | `shower_clustering_with_nv(acc_segment_id, pi0_showers, map_shower_pio_id, map_pio_id_showers, map_pio_id_mass, map_pio_id_saved_pair, pio_kine, vertices_in_long_muon, segments_in_long_muon, graph, main_vertex, showers, main_cluster, other_clusters, map_cluster_main_vertices, map_vertex_in_shower, map_segment_in_shower, map_vertex_to_shower, used_shower_clusters, track_fitter, dv, particle_data, recomb_model)` | all state now explicit |
| `shower_clustering_with_nv_in_main_cluster()` | `shower_clustering_with_nv_in_main_cluster(graph, main_vertex, showers, ...)` | |
| `shower_clustering_connecting_to_main_vertex()` | `shower_clustering_connecting_to_main_vertex(graph, main_vertex, showers, ...)` | |
| `shower_clustering_with_nv_from_main_cluster()` | `shower_clustering_with_nv_from_main_cluster(graph, main_vertex, main_cluster, ...)` | |
| `shower_clustering_with_nv_from_vertices()` | `shower_clustering_with_nv_from_vertices(graph, main_vertex, main_cluster, other_clusters, ..., track_fitter, dv, ...)` | |
| `shower_clustering_in_other_clusters()` | `shower_clustering_in_other_clusters(graph, main_vertex, showers, ...)` | |
| `examine_merge_showers()` | `examine_merge_showers(showers, main_vertex, ..., track_fitter, dv, ...)` | |
| `examine_shower_1()` | `examine_shower_1(graph, main_vertex, showers, ..., track_fitter, dv, ...)` | |
| `examine_showers()` | `examine_showers(graph, main_vertex, showers, ..., track_fitter, dv, ...)` | |
| `id_pi0_with_vertex()` | `id_pi0_with_vertex(acc_segment_id, pi0_showers, ..., graph, main_vertex, ...)` | |
| `id_pi0_without_vertex()` | `id_pi0_without_vertex(acc_segment_id, pi0_showers, ..., graph, main_vertex, ...)` | |
| `update_shower_maps(...)` | `update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, map_vertex_to_shower, used_shower_clusters)` | |

---

### `NeutrinoID_energy_reco.h` → `NeutrinoEnergyReco.cxx`

Charge-based energy reconstruction.

| Prototype | Toolkit | Notes |
|---|---|---|
| `collect_2D_charges()` (member var side-effect) | `collect_charge_maps(track_fitter)` → populates `m_charge_2d_u/v/w`, `m_map_apa_ch_plane_wires` | called once at start of `shower_clustering_with_nv` |
| `cal_kine_charge(shower)` | `cal_kine_charge(shower, graph, track_fitter, dv)` (convenience) or fast overload with pre-collected maps | |
| `cal_kine_charge(sg)` | `cal_kine_charge(segment, graph, track_fitter, dv)` | |
| `cal_corr_factor(point, offset_u, slope_yu, ...)` | `cal_corr_factor(point, track_fitter, dv)` | geometry params now come from `dv` |
| `calculate_shower_kinematics(...)` | `calculate_shower_kinematics(showers, vertices_in_long_muon, segments_in_long_muon, graph, track_fitter, dv, particle_data, recomb_model)` | |

---

## Not Yet Ported (Empty Stub Files)

The following prototype headers have corresponding **empty** toolkit `.cxx` files.
These are the next porting targets.

### `NeutrinoID_cosmic_tagger.h` + parts of `NeutrinoID_nue_tagger.h` → `NeutrinoTaggerCosmic.cxx` ✓ PORTED

| Prototype (`WCPPID::NeutrinoID::`) | Toolkit (`PatternAlgorithms::`) | Notes |
|---|---|---|
| `cosmic_tagger()` | `cosmic_tagger(graph, main_vertex, showers, map_segment_in_shower, map_vertex_to_shower, segments_in_long_muon, main_cluster, all_clusters, dv, ti)` | 10 cosmic rejection flags; fills `TaggerInfo` BDT features |
| `bad_reconstruction(shower)` | `bad_reconstruction(graph, main_vertex, shower, flag_fill, ti)` | 3 sub-checks: long stem, muon continuation via `find_cont_muon_segment_nue`, track-like continuation near shower start |

**Key translation notes for this file:**
- `TVector3::Theta()` / `TVector3::Phi()` → `vec_theta(dir)` / `vec_phi(dir)` (static helpers, since `D3Vector` has no spherical angle methods)
- `Shower::get_last_segment_vertex_long_muon(IndexedSegmentSet)` takes a plain `std::set<SegmentPtr>` (no comparator); convert from `IndexedSegmentSet` at call site
- `main_cluster->get_cluster_id()` used for flag 9 cluster-PCA block and flag 10 front-face vertex check



### `NeutrinoID_numu_tagger.h` → `NeutrinoTaggerNuMu.cxx` (EMPTY)

Entry points:
- `std::pair<bool, double> WCPPID::NeutrinoID::numu_tagger()`
- `std::pair<int, int> WCPPID::NeutrinoID::count_daughters(ProtoSegment*)` / `count_daughters(WCShower*)`

### `NeutrinoID_nue_tagger.h` → `NeutrinoTaggerNuE.cxx` (EMPTY)

Entry point: `bool WCPPID::NeutrinoID::nue_tagger(double muon_length)`

Helper functions also in this file (may become member functions or free functions):
- `low_energy_michel`, `single_shower`, `angular_cut`, `track_overclustering`

`NeutrinoID_nue_functions.h` (additional nue helpers) also maps here:
- `stem_direction`, `multiple_showers`, `other_showers`, `stem_length`, `vertex_inside_shower`, `compare_muon_energy`

### `NeutrinoID_pio_tagger.h` → `NeutrinoTaggerPi0.cxx` (EMPTY)

*Note: file exists in prototype but is empty in the grep output above — may be subsumed by shower clustering.*

### `NeutrinoID_singlephoton_tagger.h` → `NeutrinoTaggerSinglePhoton.cxx` (EMPTY)

Entry point: `bool WCPPID::NeutrinoID::singlephoton_tagger(double muon_length)`

Helper: `low_energy_michel_sp`

### `NeutrinoID_ssm_tagger.h` → `NeutrinoTaggerSSM.cxx` (EMPTY)

Entry point: `bool WCPPID::NeutrinoID::ssm_tagger()`

### `NeutrinoID_numu_bdts.h` (no toolkit file yet)

BDT score calculation for numu selection:
- `cal_numu_bdts_xgboost()`, `cal_numu_bdts()`
- Individual BDT sub-scores: `cal_cosmict_2_4_bdt`, `cal_cosmict_3_5_bdt`, `cal_cosmict_6_bdt`, `cal_cosmict_7_bdt`, `cal_cosmict_8_bdt`, `cal_cosmict_10_bdt`, `cal_numu_1_bdt`, `cal_numu_2_bdt`, `cal_numu_3_bdt`

### `NeutrinoID_nue_bdts.h` (no toolkit file yet)

BDT score calculation for nue selection:
- `cal_bdts_xgboost()`, `cal_bdts()`
- Many individual BDT sub-scores: `cal_mipid_bdt`, `cal_gap_bdt`, `cal_hol_lol_bdt`, `cal_cme_anc_bdt`, `cal_mgo_mgt_bdt`, `cal_br1_bdt`, `cal_br3_bdt`, `cal_br3_3_bdt`, `cal_br3_5_bdt`, `cal_br3_6_bdt`, `cal_stemdir_br2_bdt`, `cal_trimuon_bdt`, `cal_br4_tro_bdt`, `cal_mipquality_bdt`, `cal_pio_1_bdt`, `cal_pio_2_bdt`, `cal_stw_spt_bdt`, `cal_vis_1_bdt`, and more.

### `NeutrinoID_kine.h` → `NeutrinoKinematics.cxx`

| Prototype (`WCPPID::NeutrinoID::`) | Toolkit (`PatternAlgorithms::`) | Notes |
|---|---|---|
| `init_tagger_info()` (NeutrinoID.cxx:2217) | `init_tagger_info(TaggerInfo& ti)` | Body: `ti = TaggerInfo{}`; C++ default-member-initializers on `TaggerInfo` replace 1200-line assignment list. Data structs live in `NeutrinoTaggerInfo.h` |
| `fill_kine_tree(KineInfo& ktree)` | `fill_kine_tree(main_vertex, showers, pio_kine, graph, track_fitter, dv, geom_helper, particle_data, recomb_model) → KineInfo` | Returns by value instead of out-param. `geom_helper` (IClusGeomHelper::pointer) added for SCE vertex correction (pass nullptr to skip). `pio_kine` (Pi0KineFeatures) replaces `kine_pio_*` member vars. |

**Data format** (`clus/inc/WireCellClus/NeutrinoTaggerInfo.h`, namespace `WireCell::Clus::PR`):
- `struct KineInfo` — reconstructed neutrino kinematics output (~20 fields)
- `struct TaggerInfo` — ~500+ BDT input features for all sub-taggers (cosmic/numu/nue/ssm/singlephoton), all with correct C++ in-class defaults

**Key prototype→toolkit translation in `fill_kine_tree`**:
- `map_vertex_segments[vtx]` → `boost::out_edges(vtx->get_descriptor(), graph)` + `graph[*ei].segment`
- `find_other_vertex(seg, vtx)` → `find_other_vertex(graph, seg, vtx)`
- `shower->get_start_segment()` → `shower->start_segment()`
- `shower->get_start_segment()->get_particle_type()` → `shower->get_particle_type()`
- `seg->get_particle_type()` → `seg->particle_info()->pdg()`
- `seg->get_kine_best()` → `seg->particle_info()->kinetic_energy()`
- `seg->get_particle_mass()` → `seg->particle_info()->mass()`
- `cal_kine_charge(seg)` → `cal_kine_charge(seg, graph, track_fitter, dv)`
- `seg->cal_kine_dQdx()` → `segment_cal_kine_dQdx(seg, recomb_model)`
- `seg->cal_kine_range()` → `cal_kine_range(segment_track_length(seg), pdg, particle_data)`
- `shower->get_start_vertex()` → `shower->get_start_vertex_and_type()`
- SCE: `mp.func_pos_SCE_correction(nu_vtx)` → `geom_helper->get_corrected_point(nu_vtx, IClusGeomHelper::SCE, apa, face)`

---

## Call Order in `TaggerCheckNeutrino::visit()` vs. Prototype Constructor

The prototype `NeutrinoID` constructor runs the full chain in sequence. The toolkit
`TaggerCheckNeutrino::visit()` covers through shower clustering. Everything marked
**[NOT IN visit()]** is prototype-only and needs a downstream consumer.

```
Prototype NeutrinoID constructor          Toolkit TaggerCheckNeutrino::visit()
─────────────────────────────────────────────────────────────────────────────────
For main cluster:
  find_proto_vertex(main, true, 2)      → pattern_algos.find_proto_vertex(..., true, 2, true)
  clustering_points(main)               → pattern_algos.clustering_points(...)
  separate_track_shower(main)           → pattern_algos.separate_track_shower(...)
  determine_direction(main)             → pattern_algos.determine_direction(...)
  shower_determing_in_main_cluster(main)→ pattern_algos.shower_determining_in_main_cluster(...)
  determine_main_vertex(main)           → pattern_algos.determine_main_vertex(...)

For each other cluster (length > 6 cm):
  same 5 calls                          → same 5 calls (identical structure)

For each other cluster (length ≤ 6 cm):
  find_proto_vertex OR init_point_seg   → same
  same 4 calls                          → same

deghosting()                            → pattern_algos.deghosting(...)

if flag_dl_vtx:
  determine_overall_main_vertex_DL()   → pattern_algos.determine_overall_main_vertex_DL(...)
  (fallback) determine_overall_main_vertex() → pattern_algos.determine_overall_main_vertex(...)

improve_vertex(main, true, true)        → pattern_algos.improve_vertex(...)
clustering_points(main)                 → pattern_algos.clustering_points(...)  [again]
examine_direction(main_vertex)          → pattern_algos.examine_direction(...)
separate_track_shower()    [global]     → (not re-called in toolkit)
collect_2D_charges()                    → (called internally in shower_clustering_with_nv)
shower_clustering_with_nv()             → pattern_algos.shower_clustering_with_nv(...)

─────── NOT YET IN TOOLKIT visit() ─────────────────────────────────────────────
init_tagger_info(tagger_info)           → NeutrinoKinematics.cxx  [PORTED]
fill_kine_tree(kine_info)               → NeutrinoKinematics.cxx  [PORTED]
cosmic_tagger()                         → NeutrinoTaggerCosmic.cxx [PORTED]
numu_tagger()                           → NeutrinoTaggerNuMu.cxx  [EMPTY]
ssm_tagger()                            → NeutrinoTaggerSSM.cxx   [EMPTY]
nue_tagger(muon_length)                 → NeutrinoTaggerNuE.cxx   [EMPTY]
singlephoton_tagger(muon_length)        → NeutrinoTaggerSinglePhoton.cxx [EMPTY]
cal_numu_bdts_xgboost()                 → (no toolkit file yet)
cal_bdts_xgboost() / cal_bdts()         → (no toolkit file yet)
```

---

## Key Differences to Remember When Porting

1. **No member-variable side effects**: prototype methods mutate `main_vertex`, `showers`, etc. as class members.  Toolkit methods take these as explicit by-ref parameters.

2. **Graph is explicit**: any operation that reads/writes the PR graph must receive `Graph& graph` as a parameter.

3. **`acc_segment_id` / `acc_vertex_id`**: in the prototype these are NeutrinoID member counters. In the toolkit, vertex IDs are managed by the `PRVertex` constructor / graph; `acc_segment_id` is passed explicitly to shower functions.

4. **`separate_track_shower()` (no-arg global)**: called once after the main vertex is found in the prototype (line ~239 of NeutrinoID.cxx). This re-runs track/shower on all clusters. The toolkit does **not** repeat this call in `TaggerCheckNeutrino::visit()`. Keep this in mind when porting tagging logic that depends on a "final" track/shower classification.

5. **Prototype `determine_main_vertex(cluster, false)`**: the `false` flag means "this is an other cluster" (not the main cluster). In the toolkit this distinction is removed — all calls use the same signature.

6. **Typo fixed**: prototype `shower_determing_in_main_cluster` → toolkit `shower_determining_in_main_cluster` (extra `i`).

7. **`tagger_info` struct**: prototype fills a flat `tagger_info` struct with many float members that become BDT inputs. Toolkit equivalent is still to be designed.
