// This file provides some helper functions to configure components from WCT
// "clus/" sub-package.  In particular, to configure MultiAlgBlobClustering
// (MABC) and its pipeline of "clustering method" components.

local wc = import "wirecell.jsonnet";

{
    /// Create a "factory" object for creating Clustering* "method" components
    /// (eg ClusteringLiveDead).
    ///
    /// The clustering_methods() function takes a number of "general" arguments
    /// with default values.  Some are common to all Clustering* method
    /// components (like "prefix" to which individual object names are appended)
    /// while others may be ignored by some Clustering* method components.
    ///
    /// This function returns an object with a number of elements, each
    /// providing a function to construct a specific Clustering* component.
    /// Each of these constructor functions accept the set of "specific"
    /// arguments with default values that are relevant to the particular
    /// Clustering* component.  The "specific" arguments are named to match the
    /// names of the configuration parameters that they pass.  Eg,
    /// "dead_live_overlap_offset".
    ///
    /// Users may override either the general or specific default values as
    /// needed for their particular needs.
    ///
    /// Users note: The factory object keywords are generally matching the name
    /// of their implementation (.cxx) source file name.  This generally (with
    /// some exceptions) the class name with "Clustering" part removed and the
    /// remaining name converted from CamelCase to snake_case.
    ///
    /// Developers note: As new Clustering* components are developed, developers
    /// should extend the factory object.
    ///
    /// Example use:
    ///
    /// local cm = clus.cluster_methods("all", dv, pcts);
    /// local cm_objs = [
    ///   cm.live_dead(),     // "ClusteringLiveDead:all", defaults okay
    ///   cm.regular("one"), // "ClusteringRegular:allone", must make names unique
    ///                      // "ClusteringRegular:alltwo", because we have a second one: 
    ///   cm.regular("two", length_cut=30*wc.cm, flag_enable_extend=true),
    ///                      // Use generic() if config support not yet added.
    ///                      // This makes a tn of "ClusterNewType:allnew".
    /// ];
    /// local mabc = g.pnode({
    ///    type: "MultiAlgBlobClustering",
    ///    data: {
    ///        clustering_methods = wc.tns(cm_objs);
    ///        ...
    ///    },
    /// }, nin=1, nout=1, uses=cm_objs + [...]); // include objects that MABC "uses" directly

    clustering_methods(prefix="", detector_volumes=null, pc_transforms=null, fiducial=null,
                       pc_name="3d", coords=["x", "y", "z"] ) :: {
        // abbreviations covering commonalities across different Clustering* method components.
        local dv_tn = wc.tn(detector_volumes),
        local dv_cfg = {detector_volumes: dv_tn},
        local fiducial_tn = wc.tn(fiducial),
        local fiducial_cfg = { fiducial: fiducial_tn },
        local pcts_tn = wc.tn(pc_transforms),
        local pcts_cfg = {pc_transforms: pcts_tn},
        local scope_cfg = {pc_name: pc_name, coords: coords},

        // Use "parent" inside of a function to call sibling functions.
        local parent = self,

        tagger_flag_transfer(name="", enable_debug=false) :: {
            type: "ClusteringTaggerFlagTransfer",
            name: prefix+name,
            data: {
                enable_debug: enable_debug,
            },
        },

clustering_recovering_bundle(name="", graph_name="relaxed") :: {
            type: "ClusteringRecoveringBundle",
            name: prefix + name,
            data: dv_cfg + pcts_cfg + scope_cfg + {
                grouping: "live",           // Which grouping to process
                array_name: "isolated",     // Array name for pcarray lookup
                pcarray_name: "perblob",    // PCArray name for blob separation
                graph_name: graph_name,     // Graph flavor for connected_blobs examine step
            },
            uses: [detector_volumes, pc_transforms],
        },

        // require_in_scope (default false): also require each candidate main to
        // pass the default-scope filter set by switch_scope, i.e. to have blobs
        // whose T0-corrected points land in the active volume.  switch_scope
        // separates the out-of-volume blobs into their own cluster that keeps an
        // inherited flag_main_cluster; without this the taggers evaluate those
        // non-physical shards (which are outside the FV by construction, so they
        // satisfy the TGM CASE-A test almost automatically).  Key emitted only
        // when true so existing compiled configs stay byte-identical.
        tagger_check_stm(name="", trackfitting_config_file="", particle_dataset="", recombination_model="",
                         require_in_scope=false) :: {
            type: "TaggerCheckSTM",
            name: prefix + name,
            data: {
                grouping: "live",           // Which grouping to process
                trackfitting_config_file: trackfitting_config_file, 
                particle_dataset: particle_dataset,
                recombination_model: recombination_model,
            } + dv_cfg + pcts_cfg
              + (if require_in_scope then { require_in_scope: true } else {})
        },

        // Through-going-muon tagger (port of prototype check_tgm).  fiducial
        // names the IFiducial for the inside/outside-FV tests (e.g. a
        // BoxFiducial spanning ALL TPCs so cathode crossers are not exiters);
        // fv_tolerance = [x_lo,x_hi,y_lo,y_hi,z_lo,z_hi] margins (negative =
        // inset).  check_neutrino_candidate (C++ default false): enable the
        // ported prototype Dijkstra path-topology neutrino veto so
        // in-beam-window bundles may be tagged; when false (default, key
        // omitted => byte-identical pre-port config) in-beam bundles are
        // never tagged through the protected branches.
        // require_chord_charge (C++ default false; keys omitted when off =>
        // byte-identical pre-fix config): before a pair of extreme points may
        // tag, require the cluster to carry charge ALONG the chord between
        // them -- no contiguous stretch longer than chord_max_gap without a
        // cluster point within chord_support_radius.  Guards against the
        // flash-time merge in examine_bundles(use_flash_t0) putting a detached
        // fragment into the tagged Cluster: its chord to the real track passes
        // the FV-only test trivially.  C++ defaults 6 cm / 30 cm.
        // chord_charge_mode (C++ default "chord"; key omitted then =>
        // byte-identical): "chord" samples the STRAIGHT segment between the
        // extremes (rejects genuinely curved tracks -- SBND evt285185
        // cluster 16 bows 10 cm off its 480 cm chord and lost a real TGM);
        // "path" instead requires a piecewise charge path through the
        // cluster's own points with no jump longer than chord_max_gap
        // (chord_support_radius unused in that mode).
        // component_extremes (C++ default false; key omitted when off): find
        // the 8 extreme points PER connected component and union them.  The
        // global scan gives each slot (max z, min x, ...) to whichever merge
        // component reaches furthest, hiding the other component's own
        // wall-exit, so a genuine through-goer inside a merged bundle can never
        // form its pair.  USE TOGETHER WITH require_chord_charge: the union
        // also creates cross-component pairs, and the chord test is what
        // rejects those.  C++ component_min_length default 10 cm.
        // component_rescue (C++ default false; key omitted when off): a
        // component SHORTER than component_min_length still donates its
        // extremes when it is path-connected (30 cm-step charge path, the
        // path-mode chord rule) to a component that passed the length cut --
        // a genuine track end that fragments into a sub-10 cm piece behind
        // small gaps keeps its wall exit (SBND evt286681 cluster 7), while a
        // detached merge-grafted speck stays dropped (path-disconnected).
        // Only consulted when component_extremes is on.
        // rescue_chord_check (C++ default false; key omitted when off): a
        // pair whose end was donated by a RESCUED component must also pass
        // the STRAIGHT-chord support test even in path mode -- a genuine
        // fragmented track end lies on its own pair's chord, but path mode
        // alone lets a rescued speck pair across TWO merged cosmics through
        // an L-shaped charge detour (SBND evt288727 cluster 6).  Only
        // consulted when component_rescue is on.
        tagger_check_tgm(name="", fiducial="", fv_tolerance=[], beam_window_low=0, beam_window_high=0, length_limit_frac=0.45, enable_case_b=true, require_in_scope=false, check_neutrino_candidate=false, require_chord_charge=false, chord_support_radius=null, chord_max_gap=null, chord_charge_mode="chord", component_extremes=false, component_min_length=null, component_rescue=false, rescue_chord_check=false) :: {
            type: "TaggerCheckTGM",
            name: prefix + name,
            data: {
                grouping: "live",
                fv_tolerance: fv_tolerance,
                beam_window_low: beam_window_low,
                beam_window_high: beam_window_high,
                length_limit_frac: length_limit_frac,
                enable_case_b: enable_case_b,
            } + dv_cfg + pcts_cfg + (if fiducial == "" then {} else { fiducial: fiducial })
              + (if require_in_scope then { require_in_scope: true } else {})
              + (if check_neutrino_candidate then { check_neutrino_candidate: true } else {})
              + (if require_chord_charge then { require_chord_charge: true } else {})
              + (if chord_support_radius == null then {} else { chord_support_radius: chord_support_radius })
              + (if chord_max_gap == null then {} else { chord_max_gap: chord_max_gap })
              // Key only when the guard is ON and the mode is non-default, so
              // a runner may pass the mode unconditionally without disturbing
              // the knob-off compiled config.
              + (if require_chord_charge && chord_charge_mode != "chord" then { chord_charge_mode: chord_charge_mode } else {})
              + (if component_extremes then { component_extremes: true } else {})
              + (if component_min_length == null then {} else { component_min_length: component_min_length })
              + (if component_rescue then { component_rescue: true } else {})
              + (if rescue_chord_check then { rescue_chord_check: true } else {}),
        },

        // Fully-contained (FC) tagger.  Records Facade::cluster_fc_check's
        // verdict as the cluster flag "FC" -- the tagger-computed sibling of
        // "TGM"/"STM", i.e. the prototype's event_type bit 2 / match_isFC.
        // Needs the steiner and fiducialutils stages ahead of it (without
        // them cluster_fc_check returns is_fc=false for every cluster).
        // require_in_scope (C++ default false; key omitted when off =>
        // byte-identical): evaluate only clusters passing switch_scope's
        // active-volume filter.
        // fiducial / fv_tolerance (both C++-default absent => keys omitted =>
        // byte-identical): redirect the DIRECT containment tests from
        // FiducialUtils to this IFiducial with these margins, so FC and TGM
        // judge containment identically.  Pass the same values as
        // tagger_check_tgm.  The dead-region / signal-processing checks keep
        // using FiducialUtils either way, exactly as TaggerCheckTGM does.
        tagger_check_fc(name="", fiducial="", fv_tolerance=[], require_in_scope=false) :: {
            type: "TaggerCheckFC",
            name: prefix + name,
            data: {
                grouping: "live",
            } + dv_cfg + pcts_cfg
              + (if fiducial == "" then {} else { fiducial: fiducial, fv_tolerance: fv_tolerance })
              + (if require_in_scope then { require_in_scope: true } else {}),
        },

        tagger_check_neutrino(name="", trackfitting_config_file="", particle_dataset="", recombination_model="", perf=false, dl_weights="", dQdx_scale=0.1, dQdx_offset=-1000.0, clus_geom_helper="", dl_vtx_rerank=true, dl_vtx_top_k=5, dl_vtx_min_accept_score=4.0, dl_vtx_score_scale=1000.0, beam_window_low=0, beam_window_high=0) :: {
            type: "TaggerCheckNeutrino",
            name: prefix + name,
            data: {
                grouping: "live",           // Which grouping to process
                trackfitting_config_file: trackfitting_config_file,
                particle_dataset: particle_dataset,
                recombination_model: recombination_model,
                perf: perf,
                dl_weights: dl_weights,     // path to SCN vertex .pth file (empty = DL disabled)
                dQdx_scale: dQdx_scale,     // scale factor for dQ passed to SCN network
                dQdx_offset: dQdx_offset,   // offset for dQ passed to SCN network
                clus_geom_helper: clus_geom_helper, // type/name of SimpleClusGeomHelper; empty = no SCE
                dl_vtx_rerank: dl_vtx_rerank,           // true → top-K + soft re-rank; false → legacy single argmax
                dl_vtx_top_k: dl_vtx_top_k,             // number of top DL voxels to re-rank (only when dl_vtx_rerank==true)
                dl_vtx_min_accept_score: dl_vtx_min_accept_score,  // min composite score to accept re-ranked DL vertex
                dl_vtx_score_scale: dl_vtx_score_scale, // scale factor on raw DL score (1.0=unscaled; ~1000 for typical ~0.005 scores)
                beam_window_low: beam_window_low,   // beam window on cluster_t0 (matched flash time); low >= high
                beam_window_high: beam_window_high, // (default) disables the gate = uBooNE single-main selection
            } + dv_cfg + pcts_cfg
        },

        // Run pattern recognition (find_proto_vertex) on the main cluster.
        // mode is passed for future use (e.g. "multiple" for multi-track mode).
        do_tracking(name="", mode="", perf=false, clus_geom_helper="") :: $.tagger_check_neutrino(name=name, perf=perf, clus_geom_helper=clus_geom_helper),

        // Run numu CC BDT scoring (TMVA/xgboost).
        // Must run AFTER tagger_check_neutrino in the visitor list.
        // XML weight files should be resolved from wire-cell-data uboone/weights/.
        // Pass empty strings to disable (scorer will skip booking and EvaluateMVA).
        numu_bdt_scorer(name="",
                        numu1_weights_xml="",
                        numu2_weights_xml="",
                        numu3_weights_xml="",
                        cosmict10_weights_xml="",
                        numu_xgboost_xml="") :: {
            type: "UbooneNumuBDTScorer",
            name: prefix + name,
            data: {
                grouping: "live",
                numu1_weights_xml:    numu1_weights_xml,
                numu2_weights_xml:    numu2_weights_xml,
                numu3_weights_xml:    numu3_weights_xml,
                cosmict10_weights_xml: cosmict10_weights_xml,
                numu_xgboost_xml:     numu_xgboost_xml,
            }
        },

        // Run nueCC BDT scoring (TMVA/xgboost).
        // Must run AFTER tagger_check_neutrino (and after numu_bdt_scorer) in the visitor list.
        // All 30 sub-BDT XML files plus the top-level XGB combiner must be resolved from
        // wire-cell-data uboone/weights/.  Pass empty strings to disable individual sub-BDTs.
        nue_bdt_scorer(name="",
                       mipid_weights_xml="",
                       gap_weights_xml="",
                       hol_lol_weights_xml="",
                       cme_anc_weights_xml="",
                       mgo_mgt_weights_xml="",
                       br1_weights_xml="",
                       br3_weights_xml="",
                       br3_3_weights_xml="",
                       br3_5_weights_xml="",
                       br3_6_weights_xml="",
                       stemdir_br2_weights_xml="",
                       trimuon_weights_xml="",
                       br4_tro_weights_xml="",
                       mipquality_weights_xml="",
                       pio_1_weights_xml="",
                       pio_2_weights_xml="",
                       stw_spt_weights_xml="",
                       vis_1_weights_xml="",
                       vis_2_weights_xml="",
                       stw_2_weights_xml="",
                       stw_3_weights_xml="",
                       stw_4_weights_xml="",
                       sig_1_weights_xml="",
                       sig_2_weights_xml="",
                       lol_1_weights_xml="",
                       lol_2_weights_xml="",
                       tro_1_weights_xml="",
                       tro_2_weights_xml="",
                       tro_4_weights_xml="",
                       tro_5_weights_xml="",
                       nue_xgboost_xml="") :: {
            type: "UbooneNueBDTScorer",
            name: prefix + name,
            data: {
                grouping: "live",
                mipid_weights_xml:       mipid_weights_xml,
                gap_weights_xml:         gap_weights_xml,
                hol_lol_weights_xml:     hol_lol_weights_xml,
                cme_anc_weights_xml:     cme_anc_weights_xml,
                mgo_mgt_weights_xml:     mgo_mgt_weights_xml,
                br1_weights_xml:         br1_weights_xml,
                br3_weights_xml:         br3_weights_xml,
                br3_3_weights_xml:       br3_3_weights_xml,
                br3_5_weights_xml:       br3_5_weights_xml,
                br3_6_weights_xml:       br3_6_weights_xml,
                stemdir_br2_weights_xml: stemdir_br2_weights_xml,
                trimuon_weights_xml:     trimuon_weights_xml,
                br4_tro_weights_xml:     br4_tro_weights_xml,
                mipquality_weights_xml:  mipquality_weights_xml,
                pio_1_weights_xml:       pio_1_weights_xml,
                pio_2_weights_xml:       pio_2_weights_xml,
                stw_spt_weights_xml:     stw_spt_weights_xml,
                vis_1_weights_xml:       vis_1_weights_xml,
                vis_2_weights_xml:       vis_2_weights_xml,
                stw_2_weights_xml:       stw_2_weights_xml,
                stw_3_weights_xml:       stw_3_weights_xml,
                stw_4_weights_xml:       stw_4_weights_xml,
                sig_1_weights_xml:       sig_1_weights_xml,
                sig_2_weights_xml:       sig_2_weights_xml,
                lol_1_weights_xml:       lol_1_weights_xml,
                lol_2_weights_xml:       lol_2_weights_xml,
                tro_1_weights_xml:       tro_1_weights_xml,
                tro_2_weights_xml:       tro_2_weights_xml,
                tro_4_weights_xml:       tro_4_weights_xml,
                tro_5_weights_xml:       tro_5_weights_xml,
                nue_xgboost_xml:         nue_xgboost_xml,
            }
        },

        // Write T_tagger and T_kine trees into the existing tracking output ROOT file.
        // Must run AFTER numu_bdt_scorer and nue_bdt_scorer (BDT scores must be filled).
        // Must run AFTER UbooneMagnifyTrackingVisitor (file must already exist to UPDATE).
        tagger_output(name="", output_filename="tracking_proj.root") :: {
            type: "UbooneTaggerOutputVisitor",
            name: prefix + name,
            data: {
                grouping: "live",
                output_filename: output_filename,
            }
        },

        pointed(name="", groupings=["live"]) :: {
            type: "ClusteringPointed",
            name: prefix+name,
            data: {
                groupings: groupings,
            },
        },
        
        test(name="") :: {
            type: "ClusteringTest",
            name: prefix+name,
            data: dv_cfg + pcts_cfg,
            uses: [detector_volumes, pc_transforms],
        },

        ctpointcloud(name="") :: {
            type: "ClusteringCTPointcloud",
            name: prefix+name,
            data: dv_cfg + pcts_cfg,
            uses: [detector_volumes, pc_transforms],
        },

        live_dead(name="", dead_live_overlap_offset=2, use_flash_t0=false, flash_t0_window=80*wc.ns) :: {
            type: "ClusteringLiveDead",
            name: prefix+name,
            data: {
                dead_live_overlap_offset: dead_live_overlap_offset,
                use_flash_t0: use_flash_t0,
                flash_t0_window: flash_t0_window,
            } + dv_cfg + scope_cfg,
            uses: [detector_volumes],
        },

        extend(name="", flag=0, length_cut=150*wc.cm, num_try=0, length_2_cut=3*wc.cm, num_dead_try=3, use_flash_t0=false, flash_t0_window=80*wc.ns) :: {
            type: "ClusteringExtend",
            name: prefix+name,
            data: {
                flag: flag,
                length_cut: length_cut,
                num_try: num_try,
                length_2_cut: length_2_cut,
                num_dead_try: num_dead_try,
                use_flash_t0: use_flash_t0,
                flash_t0_window: flash_t0_window,
            } + dv_cfg + scope_cfg,
            uses: [detector_volumes],
        },


        regular(name="",  length_cut=45*wc.cm, flag_enable_extend=true, use_flash_t0=false, flash_t0_window=80*wc.ns) :: {
            type: "ClusteringRegular",
            name: prefix+name,
            data: {
                length_cut: length_cut,
                flag_enable_extend: flag_enable_extend,
                use_flash_t0: use_flash_t0,
                flash_t0_window: flash_t0_window,
            } + dv_cfg + scope_cfg,
            uses: [detector_volumes],
        },

        parallel_prolong(name="", length_cut=35*wc.cm, use_flash_t0=false, flash_t0_window=80*wc.ns) :: {
            type: "ClusteringParallelProlong",
            name: prefix+name,
            data: {
                length_cut: length_cut,
                use_flash_t0: use_flash_t0,
                flash_t0_window: flash_t0_window,
            } + dv_cfg + scope_cfg,
            uses: [detector_volumes],
        },

        close(name="", length_cut=1*wc.cm, use_flash_t0=false, flash_t0_window=80*wc.ns) :: {
            type: "ClusteringClose",
            name: prefix+name,
            data: {
                length_cut: length_cut,
                use_flash_t0: use_flash_t0,
                flash_t0_window: flash_t0_window,
            } + scope_cfg,
        },

        // SBND cathode-crossing connector (default-OFF, retireable; see
        // clus/docs/cathode-crossing-clustering.md).  Connects the two halves of a
        // cathode-crossing track left unmerged by the generic passes, using a narrow
        // cathode-specific cut set (collinear + opposite TPCs + both ends at the
        // cathode + same drift depth).  The 3D closest-point distance is handled in two
        // regimes: below dis_cut, accept on the local (Hough) track collinearity alone;
        // from dis_cut to max_dis (large transverse / in-cathode-plane travel), the local
        // Hough direction can be unreliable for a blobby half, so the cluster PCA axis is
        // added as an ALTERNATIVE direction (Hough OR PCA), and the p1->p2 connection
        // vector must align with the track within conn_far_cut (rejects parallel-offset
        // cosmics).  Cannot fire within a single TPC, so it is safe to add to the all-APA
        // pipeline only.  cathode_x is the cathode position in the T0-corrected frame.
        // use_flash_t0 (default true) gates pairs on flash-time coincidence; set false
        // on detectors without flash matching (e.g. PDHD) where the gate would veto
        // every pair.
        cathode_connect(name="", drift_cut=5*wc.cm, dis_cut=5*wc.cm, max_dis=25*wc.cm,
                        angle_cut=10.0, conn_far_cut=30.0, cathode_x=0.0,
                        cathode_x_cut=3.5*wc.cm, hough_radius=20*wc.cm,
                        min_length=10*wc.cm, min_length_short=null,
                        short_dir_len=null, conn_short_cut=30.0,
                        tip_touch_cut=null, tip_touch_angle_cut=null,
                        use_flash_t0=true, flash_t0_window=80*wc.ns,
                        crosser_conn_relax=null, crosser_pca_angle=null,
                        cathode_band_dis=null) :: {
            type: "ClusteringCathodeConnect",
            name: prefix+name,
            data: {
                drift_cut: drift_cut,
                dis_cut: dis_cut,
                max_dis: max_dis,
                angle_cut: angle_cut,
                conn_far_cut: conn_far_cut,
                cathode_x: cathode_x,
                cathode_x_cut: cathode_x_cut,
                hough_radius: hough_radius,
                min_length: min_length,
                // null => C++ defaults min_length_short to min_length (symmetric gate)
                [if min_length_short != null then "min_length_short"]: min_length_short,
                // null => C++ defaults short_dir_len to 0 (short-stub prolongation OFF)
                [if short_dir_len != null then "short_dir_len"]: short_dir_len,
                conn_short_cut: conn_short_cut,
                // null => C++ defaults tip_touch_cut to 0 (cc_pca tip-touch relaxation OFF)
                [if tip_touch_cut != null then "tip_touch_cut"]: tip_touch_cut,
                // null => C++ defaults tip_touch_angle_cut to angle_cut (local-Hough fallback OFF)
                [if tip_touch_angle_cut != null then "tip_touch_angle_cut"]: tip_touch_angle_cut,
                // null => C++ defaults crosser_conn_relax to 0 (6cm-cathode cc_pca relaxation OFF)
                [if crosser_conn_relax != null then "crosser_conn_relax"]: crosser_conn_relax,
                // null => C++ defaults crosser_pca_angle to 0 (6cm-cathode tt_pca bound raise OFF)
                [if crosser_pca_angle != null then "crosser_pca_angle"]: crosser_pca_angle,
                // null => C++ defaults cathode_band_dis to 0 (near-cathode closest-approach retry OFF)
                [if cathode_band_dis != null then "cathode_band_dis"]: cathode_band_dis,
                use_flash_t0: use_flash_t0,
                flash_t0_window: flash_t0_window,
            } + scope_cfg,
        },

        extend_loop(name="", num_try=0, use_flash_t0=false, flash_t0_window=80*wc.ns) :: {
            type: "ClusteringExtendLoop",
            name: prefix+name,
            data: {
                num_try: num_try,
                use_flash_t0: use_flash_t0,
                flash_t0_window: flash_t0_window,
            } + dv_cfg + scope_cfg,
            uses: [detector_volumes],
        },

        // max_hull_points: cap on points for the separation-decision convex hull
        // (Cluster::get_hull). -1 (default) uses Constants::MaxHullPoints (10000),
        // i.e. bit-identical to prior behavior; raise it to let large full-detector
        // overclusters be considered for separation.
        separate(name="", use_ctpc=true, max_hull_points=-1, sbnd_boundary_tag=false,
                 collinear_recover=false, collinear_interior=false,
                 collinear_member_merge=false,
                 track_repartition=false, band_merge_back=false, band_recarve=false,
                 drift_side_fv_x=false,
                 far_point_x_cut=null, far_point_mid_dis=null, track_recarve=false,
                 dec1_guard_main_angle=null, iso_slab_split=false, tag_family=false,
                 collinear_global_merge=false) :: {
            type: "ClusteringSeparate",
            name: prefix+name,
            data: {
                use_ctpc: use_ctpc,
                max_hull_points: max_hull_points,
                // SBND-only two-track upstream-boundary tag; key omitted when false
                // so existing (non-SBND) configs stay bit-identical.
                [if sbnd_boundary_tag then 'sbnd_boundary_tag']: sbnd_boundary_tag,
                // Post-separation refinements (PDVD/PDHD): recover stranded
                // collinear track tips / re-carve two crossing isochronous bands.
                // Keys omitted when false so existing configs stay bit-identical.
                [if collinear_recover then 'collinear_recover']: collinear_recover,
                // Interior-bite reclaim extension of collinear_recover: at a track
                // crossing, Separate_2's 5 cm relink can absorb an interior segment
                // of one track into the other's cluster; reclaim it.  Only effective
                // when collinear_recover is also on.  Key omitted when false.
                [if collinear_interior then 'collinear_interior']: collinear_interior,
                // Rejoin a single straight track the carve cut into long thin
                // touching collinear pieces.  Key omitted when false.
                [if collinear_member_merge then 'collinear_member_merge']: collinear_member_merge,
                // Pairwise k=2 3D repartition of two crossing thin-track family
                // members: fixes a mid-track segment of one track fused into the
                // other's cluster at the crossing.  Key omitted when false.
                [if track_repartition then 'track_repartition']: track_repartition,
                // Re-assemble a single isochronous band that the carve hatched
                // into interleaved pieces (keeps distinct parallel and genuinely
                // crossing bands apart).  Key omitted when false.
                [if band_merge_back then 'band_merge_back']: band_merge_back,
                [if band_recarve then 'band_recarve']: band_recarve,
                // Drift-side FV x-range for common-face multi-APA scopes (drift
                // groups): the out-of-time apparent-x test uses the group's drift
                // side instead of the cryostat overall x.  Key omitted when false
                // so existing configs stay bit-identical.
                [if drift_side_fv_x then 'drift_side_fv_x']: drift_side_fv_x,
                // Drift-x deviation promoting a boundary point to a "far" point in
                // JudgeSeparateDec_2's two-endpoint test.  null (default) keeps the
                // prototype-exact 140 cm (effectively dead); PDHD/PDVD set the
                // evidently intended 14 cm.  Key omitted when null so existing
                // configs stay bit-identical.
                [if far_point_x_cut != null then 'far_point_x_cut']: far_point_x_cut,
                // Midpoint-to-cluster cap in the same far-point test.  null
                // (default) keeps the prototype-exact 25 cm; raise it so two
                // diverging/forking tracks keep their far-point evidence.
                [if far_point_mid_dis != null then 'far_point_mid_dis']: far_point_mid_dis,
                // Post-separation k=2 3D-line self-split of a member holding two
                // long crossing track arms (an "X" that pure connectivity cannot
                // hold apart).  Key omitted when false: bit-identical.
                [if track_recarve then 'track_recarve']: track_recarve,
                // Dec_1 drift-aligned protection guard applies only when the
                // cluster MAIN axis is within this angle (deg) of drift.  null
                // (default) keeps the legacy unconditional guard, which wide
                // isochronous/multi-track complexes trip by accident.
                [if dec1_guard_main_angle != null then 'dec1_guard_main_angle']: dec1_guard_main_angle,
                // x-slab-aware split of a member mixing isochronous bands (one
                // dense narrow x-slab each) with drift-direction tracks that
                // chain them together under pure connectivity.  Key omitted
                // when false: bit-identical.
                [if iso_slab_split then 'iso_slab_split']: iso_slab_split,
                // Stamp final family members with a "sep_family" cluster scalar
                // so a later same-stage pass (connect1 respect_separate_family)
                // can decline to undo the split.  Key omitted when false.
                [if tag_family then 'tag_family']: tag_family,
                // Grouping-wide end-to-end stitch of two long thin collinear
                // clusters (same gates as collinear_member_merge, applied
                // beyond one separation family).  Key omitted when false.
                [if collinear_global_merge then 'collinear_global_merge']: collinear_global_merge,
            } + dv_cfg + pcts_cfg + scope_cfg,
            uses: [detector_volumes, pc_transforms],
        },

        // iso_max_dis (default null/-1 == OFF, byte-identical): upper bound on the
        // actual cluster-to-cluster closest-point distance for the isochronous-relaxed
        // connection branch, which otherwise merges two separated isochronous tracks on
        // the (misleadingly small) infinite-line distance.  SBND sets a finite value.
        // allow_mixed_faces (default null == same-face required): waive the same-face
        // requirement of the multi-wpid drift-group validation when running at the
        // per-drift-group scope (PDVD: both faces of a CRP share one drift volume).
        // respect_separate_family (default false == byte-identical): refuse to
        // reconnect two clusters that a same-stage separate(tag_family=true)
        // deliberately split apart.
        connect1(name="", use_flash_t0=false, flash_t0_window=80*wc.ns, iso_max_dis=null, allow_mixed_faces=null, respect_separate_family=false) :: {
            type: "ClusteringConnect1",
            name: prefix+name,
            data: {
                use_flash_t0: use_flash_t0,
                flash_t0_window: flash_t0_window,
                [if iso_max_dis != null then 'iso_max_dis']: iso_max_dis,
                [if allow_mixed_faces != null then 'allow_mixed_faces']: allow_mixed_faces,
                [if respect_separate_family then 'respect_separate_family']: respect_separate_family,
            } + dv_cfg + scope_cfg,
            uses: [detector_volumes],
        },

        // allow_mixed_faces: as for connect1() above.
        // empty_view_unique (default null/OFF == byte-identical): group-stage
        // semantics for empty per-(face,apa) 2D indices — a point whose view holds
        // no other cluster counts as unique evidence instead of a bogus overlap
        // (else the first cluster of each not-yet-seeded volume is wrongly
        // destroyed).  Set true on per-drift-group instances.
        deghost(name="", use_ctpc=true, length_cut=0, allow_mixed_faces=null, empty_view_unique=null) :: {
            type: "ClusteringDeghost",
            name: prefix+name,
            data: {
                use_ctpc: use_ctpc,
                length_cut: length_cut,
                [if allow_mixed_faces != null then 'allow_mixed_faces']: allow_mixed_faces,
                [if empty_view_unique != null then 'empty_view_unique']: empty_view_unique,
            } + dv_cfg + pcts_cfg + scope_cfg,
            uses: [detector_volumes, pc_transforms],
        },

        // length_cut / range_cut (default null): small/big classification
        // thresholds. When null the key is omitted and ClusteringIsolated falls
        // back to its built-in defaults (20 cm / 150), so existing configs stay
        // byte-identical. Set to opt into a tighter/looser threshold.
        isolated(name="", use_flash_t0=false, flash_t0_window=80*wc.ns, length_cut=null, range_cut=null) :: {
            type: "ClusteringIsolated",
            name: prefix+name,
            data: {
                use_flash_t0: use_flash_t0,
                flash_t0_window: flash_t0_window,
                [if length_cut != null then 'length_cut']: length_cut,
                [if range_cut != null then 'range_cut']: range_cut,
            } + dv_cfg + scope_cfg,
        },

        // flags_from_longest (default false): on the flash-time merge, take the
        // merged cluster's flags from the same representative member that donates
        // its flash instead of from an arbitrary (last-visited) member, so a
        // matched main cannot lose flag_main_cluster to a co-merged fragment.
        // Key emitted only when true so existing compiled configs stay
        // byte-identical.  See merge_clusters() in clus/inc/.../ClusteringFuncs.h.
        examine_bundles(name="", graph_name="relaxed", use_flash_t0=false, flash_t0_window=80*wc.ns,
                        flags_from_longest=false) :: {
            type: "ClusteringExamineBundles",
            name: prefix+name,
            data: dv_cfg + pcts_cfg + scope_cfg + {
                graph_name: graph_name,
                use_flash_t0: use_flash_t0,
                flash_t0_window: flash_t0_window,
            } + (if flags_from_longest then { flags_from_longest: true } else {}),
            uses: [detector_volumes, pc_transforms],
        },

        // allow_mixed_faces (default false): waive the same-face requirement on
        // multi-wpid groupings (NOT the identical-FV_x-metadata one) for
        // detectors where both faces of an anode share one drift volume
        // (PDVD: faces are the y-halves of one CRP).  Key emitted only when
        // true so existing compiled configs stay byte-identical.
        examine_x_boundary(name="", allow_mixed_faces=false) :: {
            type: "ClusteringExamineXBoundary",
            name: prefix+name,
            data: dv_cfg + scope_cfg + {
                [if allow_mixed_faces then "allow_mixed_faces"]: true,
            },
            uses: [detector_volumes],
        },

        protect_overclustering(name="") :: {
            type: "ClusteringProtectOverclustering",
            name: prefix+name,
            data: dv_cfg + pcts_cfg + scope_cfg,
            uses: [detector_volumes, pc_transforms],
        },

        // protect_iso_band (default false == byte-identical): decline to merge
        // an isochronous band (narrow drift slab, large y-z footprint) with a
        // non-band cluster unless the two genuinely touch — the extended-cloud
        // prolongations otherwise bridge tens of cm.
        neutrino(name="", num_try=1, use_flash_t0=false, flash_t0_window=80*wc.ns, protect_iso_band=false) :: {
            type: "ClusteringNeutrino",
            name: prefix+name,
            data: {
                num_try: num_try,
                use_flash_t0: use_flash_t0,
                flash_t0_window: flash_t0_window,
                [if protect_iso_band then 'protect_iso_band']: protect_iso_band,
            } + dv_cfg + scope_cfg,
            uses: [detector_volumes],
        },

        switch_scope(name="", correction_name="T0Correction") :: {
            type: "ClusteringSwitchScope",
            name: prefix+name,
            data: {
                correction_name: correction_name,
            } + pcts_cfg + scope_cfg,
            uses: [pc_transforms],
        },

        // This configures RetileCluster, a per-cluster helper for
        // ClusteringRetile as well as others.  Use the sampler() function to
        // provide properly formed elements to the array-of-object argument
        // "samplers".
        retiler(name="", anodes=[], samplers=[], cut_time_low=-1e9, cut_time_high=1e9) :: {
            local sampler_objs = [s.sobj for s in samplers],
            local sampler_cfgs = [{name:wc.tn(s.sobj), apa:s.apa, face:s.face} for s in samplers],
            type: "RetileCluster",
            name: prefix+name,
            data: {
                cut_time_low: cut_time_low,
                cut_time_high: cut_time_high,
                anodes: wc.tns(anodes),
                samplers: sampler_cfgs,
            } + dv_cfg + pcts_cfg,
            uses: [detector_volumes, pc_transforms]+anodes+sampler_objs,
        },

        // Use the sampler() function to provide properly formed elements to the
        // array-of-object argument "samplers".
        retile(name="", retiler={}) :: {
            // local sampler_objs = [s.sobj for s in samplers],
            // local sampler_cfgs = [{name:wc.tn(s.sobj), apa:s.apa, face:s.face} for s in samplers],
            // local rc = parent.retiler(name, anodes, samplers, cut_time_low, cut_time_high),
            type: "ClusteringRetile",
            name: prefix+name,
            data: {
                retiler: wc.tn(retiler),
            } + scope_cfg,
            uses: [retiler],
        },

        improve_cluster_1(name="", anodes=[], samplers=[]) :: {
            local sampler_objs = [s.sobj for s in samplers],
            local sampler_cfgs = [{name:wc.tn(s.sobj), apa:s.apa, face:s.face} for s in samplers],
            type: "ImproveCluster_1",
            name: prefix+name,
            data: {
                anodes: wc.tns(anodes),
                samplers: sampler_cfgs,
            } + dv_cfg + pcts_cfg,
            uses: [detector_volumes, pc_transforms]+anodes+sampler_objs,
        },

        // This configures ImproveCluster_2, which inherits from ImproveCluster_1
        // and adds advanced Steiner tree improvements.
        improve_cluster_2(name="", anodes=[], samplers=[], verbose=true) :: {
            local sampler_objs = [s.sobj for s in samplers],
            local sampler_cfgs = [{name:wc.tn(s.sobj), apa:s.apa, face:s.face} for s in samplers],
            type: "ImproveCluster_2",
            name: prefix+name,
            data: {
                anodes: wc.tns(anodes),
                samplers: sampler_cfgs,
                verbose: verbose,
            } + dv_cfg + pcts_cfg,
            uses: [detector_volumes, pc_transforms]+anodes+sampler_objs,
        },
        

        // Use an ImproveCluster_1 retiler for clustering retile operations
        improve_retile_1(name="", improver={}) :: {
            type: "ClusterImprove_1",
            name: prefix+name,
            data: {
                retiler: wc.tn(improver),
            } + scope_cfg,
            uses: [improver],
        },

        // Use an ImproveCluster_2 retiler for clustering retile operations
        improve_retile_2(name="", improver={}) :: {
            type: "ClusterImprove_2",
            name: prefix+name,
            data: {
                retiler: wc.tn(improver),
            } + scope_cfg,
            uses: [improver],
        },

        // Run steiner-related on clusters in grouping, saving graph to them of the given name.
        // require_beam_flash=true (uBooNE): only beam_flash-flagged clusters; false
        // (post-QL-matching detectors without that flag): every scope-passing cluster.
        steiner(name="", retiler={}, grouping="live", graph="steiner", perf=true, require_beam_flash=true) :: {
            type: "CreateSteinerGraph",
            name: prefix+name,
            data: {
                grouping: grouping,
                graph: graph,
                retiler: wc.tn(retiler),
                perf: perf,
                require_beam_flash: require_beam_flash,
            } + dv_cfg + pcts_cfg,
            uses: [detector_volumes, pc_transforms, retiler]
        },

        // Add a "FiducialUtils" to a grouping.  
        fiducialutils(name="", live_grouping="live", dead_grouping="dead", target_grouping="live") :: {
            type: "MakeFiducialUtils",
            name: prefix+name,
            data: {
                live: live_grouping,
                dead: dead_grouping,
                target: target_grouping,
            } + dv_cfg + fiducial_cfg + pcts_cfg
        },

    }, // clustering_methods(),

    /// Use this function to provide the elements of retile's "samplers"
    /// array-of-objects parameter.  It requires a configuration object
    /// for an IBlobSampler component as first argument.
    sampler(sampler_object, apa=0, face=0) :: { sobj:sampler_object, apa: apa, face: face},

    test: {
        cm : $.clustering_methods(detector_volumes={type:"DetectorVolumes", name:""}),
        ld : self.cm.live_dead(),
    }

}
