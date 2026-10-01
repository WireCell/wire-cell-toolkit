// dune10kt-1x2x6/pr.jsonnet -- neutrino pattern recognition (PR) for DUNE FD-HD 1x2x6.
//
// Counterpart of the SBND neutrino tail (sbnd/clus.jsonnet pr(), production chain of
// sbnd-pr-stage.jsonnet) on the PDHD-style APA geometry of dune10kt-1x2x6/clus.jsonnet.
// Forked BY DUPLICATION (CLAUDE.md M10): sbnd/ and pdhd/ files are untouched.  wcfm/docs/24 (scope)
// and wcfm/docs/25 (this job) in wcp-porting-validation.
//
// What is taken from where:
//   * geometry and plumbing (12 APAs x 2 live faces, wrapped U/V, DetectorVolumes, the per-(anode,
//     face) retile samplers, the PdvdPrMagnifyTrackingVisitor writer) -- PDHD (pdhd/pr.jsonnet),
//     the same APA type;
//   * the TaggerCheckNeutrino operating point, the steiner / protect_bundle / FC settings and the
//     Bee layers -- SBND production (doc sbnd_xin/120 sec 4, 118, 119), the only neutrino tail
//     validated end to end (doc sbnd_xin/115).  sbnd_tcn_knobs below is SBND's compiled
//     TaggerCheckNeutrino configuration of toolkit 0319ea67 (wct-pr.jsonnet), minus the keys the
//     helper sets itself and minus the FD overrides (fd_tcn_overrides).
// What is different for FD, and why:
//   * NO LIGHT.  Q/L matching normally makes the main + associated bundle the neutrino PR takes.
//     ClusteringBundleNoLight (stage 'bundle_no_light', first) makes the whole simulated event one
//     bundle: longest cluster = main, the rest associated, gid 0, cluster_t0 0 (the true t0 of
//     the simulation).  Not for data or pile-up.  The beam window is the whole readout and
//     contains t0 = 0, so every beam-window gate keeps the bundle.  unmerge_bundle /
//     unmerge_assoc (Q/L provenance), the cosmic taggers TGM / STM, the uBooNE-trained BDT
//     scorers and the DL (SCN, uBooNE-trained) vertex are not used.
//   * x = 0 IS AN ANODE PLANE (the 12 APAs; cathodes at +-3629 mm).  Every PR knob that holds a
//     single cathode plane at cathode_x (default 0) would act on the APA plane: the kink and MCS
//     cathode cuts, the long-muon cathode bridge and protect_bundle's cathode re-join are off.
//   * recombination: Modified Box A 0.93, B 0.212 at 0.5 kV/cm = the LArG4 model the FD depos were
//     made with, paired with particle_dataset.jsonnet (same model and field).
//   * track fitting: dune10kt_track_fitting.json (FD diffusion, PDHD transverse seeds, electron
//     lifetime correction OFF for these flash-less clusters until its FD behaviour is understood).
//   * cosmic_y_* re-anchored to the FD top (600.12 cm) as pdhd/wct-pr-perevt.jsonnet does;
//     vertex_z_prior_scale by the uBooNE rule 200 cm x 1392.6 / 1037 = 268.6 cm.
// Not calibrated on FD: mip_dqdx / mip_dqdx_median, muon_dqdx_curve, the kine_* recombination
// factors, the transverse fit widths (all SBND/PDHD values at the same 0.5 kV/cm field).
local wc = import 'wirecell.jsonnet';
local g = import 'pgraph.jsonnet';
local clus = import 'pgrapher/common/clus.jsonnet';
local P = import 'pgrapher/experiment/dune10kt-1x2x6/clus_params.jsonnet';

// SBND production TaggerCheckNeutrino knobs (compiled from toolkit 0319ea67
// cfg/pgrapher/experiment/sbnd/wct-pr.jsonnet), minus helper-set and FD-overridden keys.
// long_muon_cathode_bridge_* sub-knobs are inert here (the bridge is off, fd_tcn_overrides).
local sbnd_tcn_knobs = {
    assoc_clear_on_merge: true,
    assoc_full_recluster: true,
    assoc_reassign_orphans: true,
    break_seg_orient: true,
    broken_muon_cluster_id_count: true,
    conn3_stitch_max: 1,
    cont_muon_dir3_30cm: true,
    cosmic_companion_min_length: 15,
    cosmic_consistent_fv: true,
    daughter_count_proto_examine_showers: true,
    daughter_count_proto_main_vertex: true,
    daughter_shower_angle_reclass_straight_guard: true,
    dir_track_median_local: true,
    dir_weak_use_score: true,
    dl_vtx_dual_chain: true,
    dual_chain_mode: "snap",
    dual_chain_transfer: true,
    dual_chain_transfer_max: 2,
    endpoint_trim_retry: true,
    es3_stub_guard: true,
    esva_ignore_empty_2d: true,
    examine_direction_dirsign_shower_in_guard: true,
    excl_t0_frame: true,
    fit_blob_coverage: 0,
    fit_exclusion: true,
    fit_vertex_min_seg_length: 1,
    iso_endpoint: true,
    iso_snap_min_dir_mag: 4,
    kine_charge_dedup: true,
    kine_charge_rebuild: true,
    kine_charge_track_ctx: true,
    kine_count_conn4_near: true,
    kine_count_guard_freed: true,
    kine_count_near_cross_cluster: true,
    kine_count_orphan_tracks: true,
    kine_dqdx_skip_zero_dx: true,
    kine_drop_stray_satellites: true,
    kine_guard_freed_impact: 20,
    kine_guard_freed_miss_deg: 30,
    kine_hadronic_dqdx: true,
    kine_long_muon_mode: 2,
    kine_mainvtx_used_guard: true,
    kine_mass_rules: true,
    kine_near_pointing_impact: 200,
    kine_near_pointing_miss_deg: 30,
    kine_proton_recom_factor: 0.51,
    kine_recom_factor: 0.87,
    kine_shower_fudge_factor: 0.86,
    kine_shower_pdg_live: true,
    kine_shower_recom_factor: 0.58,
    kink_break_protect: true,
    kink_walk_dqdx_stop: true,
    long_muon_angle_relax_long: true,
    long_muon_cathode_bridge_lever: 15,
    long_muon_cathode_bridge_short_gap: 8,
    long_muon_cathode_bridge_tail_min_len: 20,
    long_muon_cathode_bridge_track_partner: true,
    long_muon_cathode_bridge_track_types: true,
    long_muon_members_geometry: true,
    long_muon_range_empty_chain_fallback: true,
    long_muon_stub_bridge: true,
    long_muon_stub_bridge_len: 7.5,
    main_vertex_candidate_flag: true,
    main_vertex_graph_audit: true,
    main_vertex_require_descriptor: true,
    main_vertex_swap_apply: true,
    mcs_bridged_members: true,
    mcs_enable: true,
    mcs_muon_source: "long_muon_else_pf",
    mcs_point_source: "muon_segments",
    mcs_range_comparator_chain: true,
    michel_stem_michel_check: true,
    michel_stem_muon_rescue: true,
    mip_dqdx: 56000,
    mip_dqdx_median: 48000,
    muon_dqdx_curve: [0.8826, 1.0587, 18, 0.4745],
    muon_multi_proton_pion: true,
    mvfit_robust: true,
    mvga_ac_chord_max: 30,
    mvga_ac_no_cascade: true,
    mvga_approach_collapse: 15,
    mvga_dup_frac: 0.8,
    mvga_dup_starved_asym: 0.55,
    mvga_dup_starved_mip: 0.8,
    mvga_dup_starved_span: 0.5,
    mvga_interposed: true,
    mvga_interposed_angle: 130,
    mvga_interposed_deg1: true,
    mvga_interposed_fallback: true,
    mvga_interposed_fallback_min_angle: 45,
    mvga_interposed_len: 10,
    mvga_op1_dup_frac: 0.7,
    mvga_op1_post: true,
    mvga_op1_radius: -1,
    mvga_passthru: 4,
    mvga_proj_angle: 25,
    mvga_proj_dqdx_ratio: 0.55,
    mvga_proj_dup_frac: 0.7,
    mvga_reseat_angle: 0,
    mvga_sat_dup_frac: 0.7,
    mvga_satellite: 3,
    mvga_splice_straighten: 5,
    mvga_straighten_radius: 1,
    mvga_stub: 2.5,
    mvga_stub_pts: 3,
    neutrino_type_bitmask: true,
    nu_bundle_flash_group: true,
    nu_fallback_demoted_mains: true,
    nu_per_bundle: true,
    nu_per_bundle_demoted_acts: true,
    nu_per_bundle_min_length: 15,
    nu_provenance: true,
    nu_selected_as_main: true,
    nu_selected_as_main_snapshot_all: true,
    nu_skip_cosmic: true,
    nu_skip_cosmic_bundle: true,
    nu_skip_cosmic_bundle_min_length: 15,
    nue_sp_consistent_fv: true,
    oov_prototype_parity: true,
    other_seg_empty_2d_guard: true,
    other_seg_keep_isolated: true,
    other_seg_keep_isolated_len_admit: 30,
    pi0_admit_muon_showers: true,
    pi0_admit_type3: true,
    pi0_attached_partner_min_mev: 29,
    pi0_bp_vertex_miss_cm: 8,
    pi0_id_shared_allocator: true,
    pi0_nc_floor_mev: 5,
    pi0_nc_frag_merge: true,
    pi0_nc_pf_assoc_deg: 20,
    pi0_nc_sig_angle_deg: 15,
    pi0_prefer_main_vertex: true,
    pi0_readmit_retyped: true,
    pid_flag_reconcile: true,
    proton_dir_vote: true,
    reclass_never_computed_ke_floor: true,
    reclass_preserve_4mom: true,
    sccc_bridge_body: true,
    sccc_kink_max: 18,
    sccc_max_gap: 10,
    sgp_max_sep: 3,
    sgp_weak_qref: 6000,
    sgp_weak_scale: 5,
    shower_absorb_track_guard: true,
    shower_absorb_unreachable_main: true,
    shower_accept_pid_guard: true,
    shower_bragg_protect_start_segment: true,
    shower_cone_absorb_guard: true,
    shower_conn3_unreachable: true,
    shower_connect_from_vertices_straight_guard: true,
    shower_connect_main_vertex_straight_guard: true,
    shower_connect_start_seg_straight_guard: true,
    shower_dedup_start_seg: true,
    shower_detach_track_stem: true,
    shower_em_collinear_deg: 10,
    shower_em_collinear_dis_cm: 120,
    shower_endpoint_exclude_start_vertex: true,
    shower_endpoint_skip_orphan_vtx: true,
    shower_ex1_dedup_rehome: true,
    shower_flag_pdg_electron: true,
    shower_ghost_member_drop: true,
    shower_hadronic_growth_max: 0.7,
    shower_hadronic_stem_ratio: 2.8,
    shower_hadronic_tag: true,
    shower_in_cascade_guard: true,
    shower_less_id_tiebreak: true,
    shower_long_muon_keep_type: true,
    shower_merge_relax: true,
    shower_merge_relax_continuity: true,
    shower_nv_bridge_track: true,
    shower_pass3_backfill_guard_len: 15,
    shower_pass3_cone_guard_len: 15,
    shower_pass4_best_owner: true,
    shower_pass4_prefilter_v1_escape: true,
    shower_pass4_prefilter_v1_max_v2: 90,
    shower_pass4_prox_guard_len: 50,
    shower_pass4_prune_detached: true,
    shower_pass4_prune_gap2: 25,
    shower_pass4_track_guard_len: 50,
    shower_pdg_exact_muon_test: true,
    shower_pdg_from_shower_type: true,
    shower_pdg_from_start_segment: true,
    shower_proton_daughter_pion: true,
    shower_proton_daughter_pion_dissolve: true,
    shower_reclass_case_b_dqdx_guard: true,
    shower_reclass_dqdx_guard: true,
    shower_samevtx_track_absorb: true,
    shower_satellite_absorb: true,
    shower_split: true,
    shower_split_em_start: true,
    shower_stem_backfill: true,
    shower_topo_demote_len: 50,
    shower_topo_dqdx_guard: true,
    shower_topo_reexam_straight_guard: true,
    shower_topo_reset: true,
    shower_traj_michel_stem: true,
    shower_traj_recheck_parity: true,
    shower_traj_straight_guard: true,
    shower_vote_track_pid_counts: true,
    shower_walk_visited_parity: true,
    single_muon_long_muon_claim: true,
    single_muon_proton_chain_veto: true,
    skip_cosmic_companions: true,
    sp_dedx_use_recomb_model: true,
    sp_mean_dedx_cut: 2.23,
    sp_photon_flag: true,
    steiner_gap_penalty: 2,
    stem_backfill_back_dvtx: 45,
    stem_backfill_back_guard: true,
    stem_endpoint_wcpt_parity: true,
    straight_cont_cross_cluster: true,
    swap_orphan_dup_audit: true,
    tagger_ordered_segment_sets: true,
    teb_bragg_veto_turn: 30,
    teb_turn_min_arm_frac: 0.4,
    track_comp_empty_abstain: true,
    track_pid_persist_4mom: true,
    track_pid_persist_dqdx: true,
    track_pid_persist_dqdx_electron_guard: true,
    two_end_break: true,
    v3_extension_guard: true,
    vertex_dir_use_fit_point: true,
    vertex_junction_snap: true,
    vertex_kink_snap: true,
    vjs_override_kink_snap: true,
};

// FD values of the SBND knobs that carry SBND geometry (wcfm/docs/25 audit table).
local fd_tcn_overrides(pr_y_top, vertex_z_prior_scale) = {
    // cathode at x = 0 in SBND; the FD APA plane here => all cuts at that plane off.
    cathode_kink_xcut: 0,            // SBND 5 cm
    cathode_wide_kink_angle: 0,      // SBND 25 deg
    mcs_cathode_xcut: 0,             // SBND 5 cm
    long_muon_cathode_bridge: false, // SBND true
    // the cosmic tagger's y cuts, same offsets below the top face as SBND / PDHD
    cosmic_y_top_main: pr_y_top - 17,
    cosmic_y_top_strict: pr_y_top - 15,
    cosmic_y_top_loose: pr_y_top - 37,
    cosmic_y_small_piece: pr_y_top - 67,
    // (z - min_z) / scale main-vertex penalty, cm
    vertex_z_prior_scale: vertex_z_prior_scale,
    // FD U/V channels wrap over BOTH live faces: map a 2-D charge cell to the wire nearest the
    // object, not the channel's first wire (C++ default false; wcfm/docs/25).
    kine_charge_all_wires: true,
};

function(output_dir='', runNo=1, subRunNo=1, eventNo=1)
{
    local clus_maker = (import 'pgrapher/experiment/dune10kt-1x2x6/clus.jsonnet')(
        output_dir=output_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo),
    local out_prefix = if output_dir == '' then '' else output_dir + '/',
    local bee_dir = if output_dir == '' then 'data' else output_dir,

    // anodes: the run's tools.anodes objects.  pipeline_names: the stage list (wct-pr-perevt.jsonnet
    // default_pipeline).  pds: particle_dataset.jsonnet().
    pr(anodes, pds, pipeline_names,
       dump=true, tensor_outname='', pr_bee=true,
       trackfitting_config_file='pgrapher/experiment/dune10kt-1x2x6/dune10kt_track_fitting.json',
       // beam window on cluster_t0 (WCT units); must contain the bundle's t0 (0).
       beam_window=[-1e4 * wc.us, 1e4 * wc.us],
       pr_y_top=600.12,                 // cm, top of the active volume (6001.2 mm)
       vertex_z_prior_scale=268.6,      // cm
       nticks=P.nticks,
       bee_detector=P.bee_detector,
       bundle_t0=0.0) :: {
        local dv = clus_maker.detector_volumes(anodes),
        local pcts = clus_maker.pc_transforms(dv),
        local cm_old = clus.clustering_methods(
            prefix='pr', detector_volumes=dv, pc_transforms=pcts, fiducial=dv, coords=clus_maker.scope_coords),
        local cm = clus.clustering_methods(
            prefix='pr', detector_volumes=dv, pc_transforms=pcts, fiducial=dv, coords=clus_maker.t0cor_coords),

        // Practical units (kV/cm, (kV/cm)(g/cm^2)/MeV, g/cm^3), as sbnd_box_recomb (doc 88).
        local fd_recomb = {
            type: 'PracticalBoxRecombination',
            name: 'fd_box_recomb',
            data: { A: 0.93, B: 0.212, Efield: 0.5, rho: 1.38, Wi: 23.6e-6 },
        },
        // The active box (no inset); the taggers' margins come in through fv_tolerance.
        local fd_pr_fv = {
            type: 'BoxFiducial',
            name: 'fd_pr_fv',
            data: {
                bounds: {
                    tail: { x: P.active_box.xmin, y: P.active_box.ymin, z: P.active_box.zmin },
                    head: { x: P.active_box.xmax, y: P.active_box.ymax, z: P.active_box.zmax },
                },
            },
        },
        // SBND production margins (x 2.5, y 3, downstream z 5, upstream z 3 cm).
        local fd_pr_fv_margins = [-2.5 * wc.cm, -2.5 * wc.cm, -3 * wc.cm, -3 * wc.cm, -5 * wc.cm, -3 * wc.cm],
        local beam_gate = beam_window[1] > beam_window[0],

        // Retiler for the steiner stage: the 'charge_stepped' retile sampler (SBND/PDHD production),
        // one per (anode, face) -- both faces of every FD anode are live.
        local improve2 = cm.improve_cluster_2(
            anodes=anodes,
            samplers=[clus.sampler(clus_maker.live_sampler(a, f, strategy_name='charge_stepped'),
                                   apa=a.data.ident, face=f)
                      for a in anodes for f in [0, 1]]),
        // SBND production steiner settings (doc sbnd_xin/118 flip: prefer3 / 0.5 / tree+path).
        local steiner_args = {
            retiler: improve2, perf: true, require_beam_flash: false,
            beam_window_only: beam_gate, beam_window_low: beam_window[0], beam_window_high: beam_window[1],
            terminal_wire_tol: 1, terminal_adjacent_slice: true, edge_charge_forward_dead_mix: true,
            terminal_blank_plane_mode: 'prefer3', base_weight_blank_alpha: 0.5, base_weight_scope: 'tree+path',
        },
        local tracking_pr_root = out_prefix + 'tracking-pr.root',

        local cm_by_name = {
            // MUST be first: the bundle and its t0 feed switch_scope's T0 correction.
            bundle_no_light: {
                type: 'ClusteringBundleNoLight',
                name: 'pr',
                data: { grouping: 'live', bundle_gid: 0, cluster_t0: bundle_t0 },
            },
            switch_scope: cm_old.switch_scope(),
            steiner: cm.steiner(retiler=steiner_args.retiler, perf=steiner_args.perf,
                                require_beam_flash=steiner_args.require_beam_flash,
                                beam_window_only=steiner_args.beam_window_only,
                                beam_window_low=steiner_args.beam_window_low,
                                beam_window_high=steiner_args.beam_window_high,
                                terminal_wire_tol=steiner_args.terminal_wire_tol,
                                terminal_adjacent_slice=steiner_args.terminal_adjacent_slice,
                                edge_charge_forward_dead_mix=steiner_args.edge_charge_forward_dead_mix,
                                terminal_blank_plane_mode=steiner_args.terminal_blank_plane_mode,
                                base_weight_blank_alpha=steiner_args.base_weight_blank_alpha,
                                base_weight_scope=steiner_args.base_weight_scope),
            fiducialutils: cm.fiducialutils(),
            tagger_check_fc: cm.tagger_check_fc(
                evaluate_demoted_mains=true,
                fiducial=wc.tn(fd_pr_fv),
                fv_tolerance=fd_pr_fv_margins,
                require_in_scope=true,
                beam_window_only=beam_gate,
                beam_window_low=beam_window[0],
                beam_window_high=beam_window[1]),
            // SBND production settings; the cathode re-join is OFF (x = 0 is the APA plane).
            protect_bundle: cm.protect_bundle(
                graph_name='relaxed_strict_img_2d_rescue_long_wtrack',
                beam_window_only=beam_gate,
                beam_window_low=beam_window[0],
                beam_window_high=beam_window[1],
                open_convicted_bundles=true),
            // rebuilds only the clusters protect_bundle purged (replace=false); always after it.
            steiner_refresh: cm.steiner(name='refresh', replace=false,
                                retiler=steiner_args.retiler, perf=steiner_args.perf,
                                require_beam_flash=steiner_args.require_beam_flash,
                                beam_window_only=steiner_args.beam_window_only,
                                beam_window_low=steiner_args.beam_window_low,
                                beam_window_high=steiner_args.beam_window_high,
                                terminal_wire_tol=steiner_args.terminal_wire_tol,
                                terminal_adjacent_slice=steiner_args.terminal_adjacent_slice,
                                edge_charge_forward_dead_mix=steiner_args.edge_charge_forward_dead_mix,
                                terminal_blank_plane_mode=steiner_args.terminal_blank_plane_mode,
                                base_weight_blank_alpha=steiner_args.base_weight_blank_alpha,
                                base_weight_scope=steiner_args.base_weight_scope),
            tagger_check_neutrino: cm.tagger_check_neutrino(
                trackfitting_config_file=trackfitting_config_file,
                particle_dataset=wc.tn(pds.particle_dataset),
                recombination_model=wc.tn(fd_recomb),
                perf=true,
                dl_weights='',               // geometric vertex (SCN is uBooNE-trained; M4)
                dl_vtx_min_accept_score=10,  // SBND value; inert without DL weights
                beam_window_low=beam_window[0],
                beam_window_high=beam_window[1],
                fiducial=wc.tn(fd_pr_fv),
                fv_tolerance=fd_pr_fv_margins,
                knobs=sbnd_tcn_knobs + fd_tcn_overrides(pr_y_top, vertex_z_prior_scale)),
            // Magnify-tracking ROOT dump; channel scheme from every anode's wires (wrapped U/V
            // channels counted once), as PDHD/PDVD production.
            tracking_visitor: {
                type: 'PdvdPrMagnifyTrackingVisitor',
                name: 'pr',
                data: {
                    grouping: 'live',
                    output_filename: tracking_pr_root,
                    runNo: runNo,
                    subRunNo: subRunNo,
                    eventNo: eventNo,
                    anodes: [wc.tn(a) for a in anodes],
                    detector_volumes: wc.tn(dv),
                    dQdx_scale: 0.1,
                    dQdx_offset: -1000.0,
                    flag_skip_vertex: false,
                    nticks: nticks,
                    save_in_scope: true,
                },
            },
            // T_tagger / T_kine (UbooneTaggerOutputVisitor, geometry-free); UPDATE mode, so after
            // tracking_visitor.  SBND production booking.
            tagger_output: cm.tagger_output(output_filename=tracking_pr_root,
                                            neutrino_type_bitmask=true,
                                            nu_per_bundle=true,
                                            mcs_output=true,
                                            nu_provenance=true,
                                            nu_particle_links=true),
            // PR event-display calib dump (docs/pr/26); after tagger_check_neutrino.
            pr_display: {
                type: 'PrDisplayDump',
                name: 'pr',
                data: {
                    grouping: 'live',
                    output_filename: out_prefix + 'calib-pr-evt' + std.toString(eventNo) + '.json',
                    runNo: runNo,
                    subRunNo: subRunNo,
                    eventNo: eventNo,
                    anodes: [wc.tn(a) for a in anodes],
                    detector_volumes: wc.tn(dv),
                    dQdx_scale: 0.1,
                    dQdx_offset: -1000.0,
                    nticks: nticks,
                    pseudo_shower_track_paint: true,
                },
            },
        },
        local cm_pipeline = [cm_by_name[n] for n in pipeline_names],
        local tcn_on = std.member(pipeline_names, 'tagger_check_neutrino'),
        local pr_visitor = 'TaggerCheckNeutrino:pr',
        local uses_extra = (if tcn_on then [fd_recomb] + pds.all else [])
                           + (if tcn_on || std.member(pipeline_names, 'tagger_check_fc') then [fd_pr_fv] else []),

        local mabc = g.pnode({
            type: 'MultiAlgBlobClustering',
            name: 'clus_pr',
            data: {
                inpath: 'pointtrees/%d',
                outpath: 'pointtrees/%d',
                perf: true,
                bee_dir: bee_dir,
                bee_zip: if pr_bee then out_prefix + 'mabc-pr.zip' else '',
                bee_detector: bee_detector,
                initial_index: 0,
                use_config_rse: true,
                runNo: runNo,
                subRunNo: subRunNo,
                eventNo: eventNo,
                reset_shower_ids_per_event: true,
                save_deadarea: true,
                dead_area_version: 2,
                save_opflash: false,
                anodes: [wc.tn(a) for a in anodes],
                detector_volumes: wc.tn(dv),
                cluster_id_order: 'tree',
                bee_points_sets: [
                    {
                        name: 'clustering-pr',
                        detector: bee_detector,
                        algorithm: 'clustering',
                        pcname: '3d',
                        coords: clus_maker.t0cor_coords,
                        individual: false,
                    },
                ] + (if tcn_on then [
                    {
                        name: 'track_fit', visitor: pr_visitor, grouping: 'live', detector: bee_detector,
                        algorithm: 'track_fit', pcname: '3d', coords: ['x', 'y', 'z'], individual: false,
                        dQdx_scale: 0.1, dQdx_offset: -1000.0, include_vertex_points: true, require_pr_graph: true,
                    },
                    {
                        name: 'shower_track', visitor: pr_visitor, grouping: 'live', detector: bee_detector,
                        algorithm: 'shower_track', pcname: '3d', coords: ['x', 'y', 'z'], individual: false,
                        particle_ids: true, pseudo_shower_track_paint: true, require_pr_graph: true,
                        use_associate_points: true,
                    },
                    {
                        name: 'vertices', visitor: pr_visitor, grouping: 'live', detector: bee_detector,
                        algorithm: 'vertices', pcname: '3d', coords: ['x', 'y', 'z'], individual: false,
                        require_pr_graph: true, use_graph_vertices: true,
                    },
                ] else []),
                // Particle-flow "mc" tree; SBND production flags, emitted after TaggerCheckNeutrino.
                bee_pf: if tcn_on then [
                    {
                        name: 'mc', visitor: pr_visitor, grouping: 'live',
                        em_ke_min: 5, np_ke_min: 3, emit_empty: true, prototype_names: true,
                        merge_metadata_key: 'bee_pf_truth', merge_node_text: 'reco nu',
                        pf_conn4_near_candidate: true, pf_direct_when_touching: true,
                        pf_drop_stray_satellites: true, pf_orphan_audit_only: true,
                        pf_orphan_confident_track: true, pf_orphan_guard_freed: true,
                        pf_orphan_near_cross_cluster: true, pf_orphan_track_parentage: true,
                        pf_pdg_name_prototype_fallback: true, pf_pi0_node_per_id: true,
                        pf_pseudo_gap_from_main: true, pf_shower_parent_precedence: true,
                        pf_shower_vertex_barrier: true, pf_track_bridged_clusters: true,
                        pf_track_main_cluster_only: true, pf_track_owns_loose_vertex: true,
                        pf_unique_node_ids: true,
                    },
                ] else [],
                pipeline: wc.tns(cm_pipeline),
            },
        }, nin=1, nout=1, uses=anodes + [dv, pcts] + cm_pipeline + uses_extra),

        local sink = g.pnode({
            type: 'TensorFileSink',
            name: 'clus_pr',
            data: {
                outname: if tensor_outname == '' then out_prefix + 'trash-pr.tar.gz' else tensor_outname,
                prefix: 'clustering_',
                dump_mode: tensor_outname == '',
            },
        }, nin=1, nout=0),

        ret: if dump then g.pipeline([mabc, sink]) else mabc,
    }.ret,
}
