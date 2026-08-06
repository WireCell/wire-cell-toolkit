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

        // Restore the prototype "main cluster + associated clusters" data
        // product on a post-Q/L tree: split every flag_main_cluster cluster
        // that the flash-time merge (examine_bundles use_flash_t0) built back
        // into its pre-merge main (retained, keeps the flags and the cluster
        // scalars) and one flag_associated_cluster cluster per other member.
        // Without it the STM tagger fits a bundle of detached cosmics as one
        // track and check_other_clusters() has no companions left to count.
        //
        // mode (C++ default "real"): "real" = per-blob real_cluster_main /
        // real_cluster_id flash-merge provenance (exact; needs a pctree saved
        // with MultiAlgBlobClustering save_real_cluster_id -- a cluster
        // WITHOUT it is left alone with a warning, NOT proxied, because a
        // connectivity split is a clustering decision and breaks cathode
        // crossers).  "component" = longest connected component, opt-in only.
        // require_in_scope (C++ default true): skip the out-of-volume shards
        // switch_scope splits off.
        //
        // MUST be placed BEFORE the steiner stage: separate() does not carry
        // node-local point clouds, and the retained main would otherwise keep
        // a steiner_pc built from the pre-split blob set.
        // id_aname / main_aname (C++ defaults "real_cluster_id" /
        // "real_cluster_main"; keys omitted when null => byte-identical pre-knob
        // config): which per-blob provenance pair to split on.  Pass
        // "assoc_cluster_id" / "assoc_cluster_main" for a SECOND instance that
        // undoes the isolated GROUPING (doc 52 Stage 3) rather than the flash
        // merge.  Run the flash one first (outer) and the isolated one second
        // (inner): that order reproduces the prototype's main_cluster +
        // additional_clusters layout, and both must precede the steiner stage.
        // restore_demoted_mains (C++ default false; key omitted when null =>
        // byte-identical pre-knob config): additionally tag a split-off part
        // that was ITSELF a matched Q/L bundle main before the flash-time merge
        // with flag_demoted_main, read from the per-blob "real_cluster_was_main"
        // array (examine_bundles save_bundle_main_provenance must be on, or the
        // visitor warns and flags nothing).  The part keeps
        // flag_associated_cluster and deliberately does NOT get
        // flag_main_cluster back -- nu_skip_cosmic_bundle builds its veto set
        // from main_cluster.  Doc pr/20 Part I P2.
        // require_provenance (C++ default false; key omitted when null =>
        // byte-identical pre-knob config): with restore_demoted_mains on, a
        // pctree with NO wasmain array is a LEGACY tree (pre pr/20 Part I) and
        // the C++ aborts instead of warn-and-skip -- the warn path runs to
        // completion silently reproducing pre-pr/20 behaviour (doc pr/23
        // sec 4.2).  Intentional legacy runs pass false explicitly.
        unmerge_bundle(name="", mode="real", graph_name="relaxed", require_in_scope=true,
                       id_aname=null, main_aname=null, restore_demoted_mains=null,
                       require_provenance=null) :: {
            type: "ClusteringUnmergeBundle",
            name: prefix + name,
            data: dv_cfg + pcts_cfg + {
                grouping: "live",
                mode: mode,
                pcarray_name: "perblob",
                graph_name: graph_name,
                require_in_scope: require_in_scope,
                [if id_aname != null then 'id_aname']: id_aname,
                [if main_aname != null then 'main_aname']: main_aname,
                [if restore_demoted_mains != null then 'restore_demoted_mains']: restore_demoted_mains,
                [if require_provenance != null then 'require_provenance']: require_provenance,
            },
            uses: [detector_volumes, pc_transforms],
        },

        // PR-stage overclustering protection: the toolkit counterpart of the
        // prototype's second graph-examination round after Q/L matching
        // (WCPPID::Protect_Over_Clustering -> PR3DCluster::Examine_graph,
        // ProtectOverClustering.cxx:6-160, wire-cell-prod-nue.cxx:1322).
        // Splits each bundle cluster at graph-component boundaries; the
        // longest component keeps the retained cluster (and, for a bundle
        // main, its main_cluster role), fragments become associated clusters
        // fit separately downstream.  Place AFTER unmerge_assoc and BEFORE
        // the steiner stage.  The stage only acts when named in
        // pipeline_names, so existing pipelines are byte-identical.
        // beam_window_only + beam_window_low/high (C++ defaults false/0/0):
        // restrict to the beam-coincident bundle(s), the prototype's scope --
        // same gate and key names as CreateSteinerGraph.
        // cathode_x / cathode_rejoin_xcut / cathode_rejoin_dyz /
        // cathode_rejoin_dis (C++ defaults 0 / 0 = pass disabled / 4 cm /
        // 8 cm; keys omitted when null => prototype-faithful behavior):
        // re-unite component pairs that meet across the cathode band before
        // splitting -- the relaxed graph does not join the two halves of a
        // cathode crosser (SBND geometry uboone did not have; doc pr/20,
        // pr/23).  Internal units.
        // cathode_rejoin_perp / cathode_rejoin_angle / cathode_rejoin_dir_radius
        // / cathode_rejoin_dir_npts (C++ defaults 0 = fallback disabled / 20.0
        // deg / 15 cm / 20 pts; keys omitted when null => prototype-faithful,
        // i.e. dyz-only, behavior): direction-agreement fallback for a pair
        // that fails ONLY cathode_rejoin_dyz.  dyz is a frame-aligned bound on
        // transverse offset and is right for a near-drift-parallel crosser,
        // but for an OBLIQUE crosser the transverse offset across the cathode
        // x-gap is real track travel and dyz rejects it by construction (doc
        // pr/25, SBND evt 489327).  cathode_rejoin_angle is degrees; the rest
        // are internal units like the block above (NOT cm, unlike
        // cathode_kink_xcut elsewhere -- the doc pr/20 trap).
        protect_bundle(name="", graph_name="relaxed",
                       beam_window_only=null, beam_window_low=null, beam_window_high=null,
                       skip_convicted=null,
                       cathode_x=null, cathode_rejoin_xcut=null,
                       cathode_rejoin_dyz=null, cathode_rejoin_dis=null,
                       cathode_rejoin_perp=null, cathode_rejoin_angle=null,
                       cathode_rejoin_dir_radius=null, cathode_rejoin_dir_npts=null) :: {
            type: "ClusteringProtectBundle",
            name: prefix + name,
            data: dv_cfg + pcts_cfg + {
                grouping: "live",
                graph_name: graph_name,
                pcarray_name: "perblob",
                [if beam_window_only != null then 'beam_window_only']: beam_window_only,
                [if beam_window_low != null then 'beam_window_low']: beam_window_low,
                [if beam_window_high != null then 'beam_window_high']: beam_window_high,
                // C++ default true (doc pr/23 ordering): a TGM/STM/lm-convicted
                // in-window main does not open its bundle.  Key omitted when
                // null => C++ default.
                [if skip_convicted != null then 'skip_convicted']: skip_convicted,
                [if cathode_x != null then 'cathode_x']: cathode_x,
                [if cathode_rejoin_xcut != null then 'cathode_rejoin_xcut']: cathode_rejoin_xcut,
                [if cathode_rejoin_dyz != null then 'cathode_rejoin_dyz']: cathode_rejoin_dyz,
                [if cathode_rejoin_dis != null then 'cathode_rejoin_dis']: cathode_rejoin_dis,
                [if cathode_rejoin_perp != null then 'cathode_rejoin_perp']: cathode_rejoin_perp,
                [if cathode_rejoin_angle != null then 'cathode_rejoin_angle']: cathode_rejoin_angle,
                [if cathode_rejoin_dir_radius != null then 'cathode_rejoin_dir_radius']: cathode_rejoin_dir_radius,
                [if cathode_rejoin_dir_npts != null then 'cathode_rejoin_dir_npts']: cathode_rejoin_dir_npts,
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
        // save_stm_fit (C++ default false; key omitted when off => byte-identical):
        // persist the per-pass STM track fits as cluster PCs
        // (stm_fit/stm_pass/stm_eval) and a grouping-level "stm" TrackFitting
        // for downstream writers (Bee stm_fit layer, SbndMagnifyTrackingVisitor).
        // mip_dqdx (C++ default 50e3 = the MicroBooNE value; null here => key
        // omitted => byte-identical pre-knob config): the MIP dQ/dx scale in
        // ELECTRONS PER CM (not per internal length unit -- unlike the
        // MIP_dQdx default arguments in PRSegmentFunctions.h, which are
        // written 50000/units::cm).  One number, two roles inside
        // TaggerCheckSTM: it sets the flat "MIP-like, not stopping" reference
        // curves the fitted dQ/dx is KS-compared against, AND it is the
        // divisor that normalizes measured dQ/dx for ~33 hard-coded cut
        // thresholds.  Raising it rescales the reference and tightens all of
        // those cuts on unchanged measured charge, so it is a behavior change,
        // not a recalibration of the reference alone.
        // fiducial / fv_tolerance (both C++ default unset; keys omitted then =>
        // byte-identical pre-fix config): the fiducial volume for the
        // cluster_fc_check containment gate that decides whether an STM fit is
        // attempted at all.  Unset means the historical fallback, FiducialUtils
        // -> DetectorVolumes::contained() = the union of per-(apa,face)
        // sensitive volumes with NO margin -- which is NOT the volume
        // tagger_check_tgm / tagger_check_fc use, and is more permissive at
        // every wall (see clus/src/TaggerCheckSTM.cxx configure()).  Pass the
        // SAME IFiducial name and the SAME margin vector given to
        // tagger_check_tgm so "contained" means one thing across all three
        // verdicts; that also matches the prototype, where check_stm and
        // check_tgm are members of one ToyFiducial and share its boundaries.
        // fv_tolerance = [x_lo,x_hi,y_lo,y_hi,z_lo,z_hi], negative = inset.
        // beam_window_only (C++ default false; keys omitted when off =>
        // byte-identical pre-knob config): evaluate ONLY the beam-coincident
        // bundle, i.e. mains whose matched flash time (cluster_t0) is in
        // [beam_window_low, beam_window_high).  Same gate as the steiner stage
        // and tagger_check_{tgm,fc}; a surviving main's verdict is unchanged.
        // accept_guards (C++ default false; key omitted when off =>
        // byte-identical pre-knob config): the doc-63 round-1 pass-level
        // acceptance guards (charge-desert one-objectness, spike-not-ramp
        // vertex veto, eval ratio2 normalization cap).  Guard thresholds keep
        // their C++ defaults (measured on the doc-62 owner baseline); only
        // the master switch is threaded here.
        // proton_muon_guard (C++ default false; key omitted when off =>
        // byte-identical): the doc-63 round-2 muon-consistency guard on
        // detect_proton's end-proton branches -- an end region matching the
        // muon hypothesis in shape and normalization is not called a proton.
        tagger_check_stm(name="", trackfitting_config_file="", particle_dataset="", recombination_model="",
                         require_in_scope=false, evaluate_demoted_mains=false,
                         save_stm_fit=false, mip_dqdx=null,
                         fiducial=null, fv_tolerance=[],
                         beam_window_only=false, beam_window_low=0, beam_window_high=0,
                         accept_guards=false, proton_muon_guard=false,
                         cathode_guard=false, anode_dist_fix=false,
                         second_track_guard=false, deficit_guard=false,
                         vertex_kink_guard=false,
                         michel_res_length_cut=null, proton_tm_max=null,
                         proton_b_ks2_max=null, proton_c_peak_max=null) :: {
            type: "TaggerCheckSTM",
            name: prefix + name,
            data: {
                grouping: "live",           // Which grouping to process
                trackfitting_config_file: trackfitting_config_file,
                particle_dataset: particle_dataset,
                recombination_model: recombination_model,
            } + dv_cfg + pcts_cfg
              + (if require_in_scope then { require_in_scope: true } else {})
              + (if evaluate_demoted_mains then { evaluate_demoted_mains: true } else {})
              + (if save_stm_fit then { save_stm_fit: true } else {})
              + (if mip_dqdx != null then { mip_dqdx: mip_dqdx } else {})
              + (if fiducial != null then { fiducial: fiducial } else {})
              + (if std.length(fv_tolerance) > 0 then { fv_tolerance: fv_tolerance } else {})
              + (if beam_window_only then {
                     beam_window_only: true,
                     beam_window_low: beam_window_low,
                     beam_window_high: beam_window_high,
                 } else {})
              + (if accept_guards then { accept_guards: true } else {})
              + (if proton_muon_guard then { proton_muon_guard: true } else {})
              + (if cathode_guard then { cathode_guard: true } else {})
              // C++ default false.  Key omitted when off => byte-identical
              // pre-fix config (doc-63 round 4a dist_to_anode face fix).
              + (if anode_dist_fix then { anode_dist_fix: true } else {})
              // C++ default false.  Key omitted when off => byte-identical
              // pre-knob config (doc-63 round 4b/4c second-track vetoes).
              + (if second_track_guard then { second_track_guard: true } else {})
              // C++ defaults false.  Keys omitted when off => byte-identical
              // pre-knob config (doc-63 round 5 stop-region vetoes).
              + (if deficit_guard then { deficit_guard: true } else {})
              + (if vertex_kink_guard then { vertex_kink_guard: true } else {})
              // doc-66 sec 12 diffusion-margin cut package.  C++ defaults are
              // the prototype constants (6 cm / 1.0 / 0.05 / 4.3); null here
              // omits the key => byte-identical legacy config.
              + (if michel_res_length_cut != null then { michel_res_length_cut: michel_res_length_cut } else {})
              + (if proton_tm_max != null then { proton_tm_max: proton_tm_max } else {})
              + (if proton_b_ks2_max != null then { proton_b_ks2_max: proton_b_ks2_max } else {})
              + (if proton_c_peak_max != null then { proton_c_peak_max: proton_c_peak_max } else {})
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
        // main_component_pairs (C++ default false; key omitted when off): a
        // pair may tag only when at least one end lies in the cluster's MAIN
        // charge component (the 30 cm-step path component with the most
        // points; a cathode crosser is ONE such component).  On a
        // flash-merged Cluster a merged-in fragment that is itself
        // through-going otherwise tags the whole bundle on its own
        // within-component pair -- which the chord guard deliberately allows
        // (SBND evt289343 cluster 9: in-beam bundle tagged TGM by a 26 cm
        // corner-clipping cosmic fragment 450 cm from the main track).  The
        // TGM verdict follows the main cluster instead.
        // main_component_mode (C++ default "path"; key omitted then =>
        // byte-identical): with main_component_pairs on, "path" identifies the
        // main as the largest 30 cm path component (a proxy); "real" reads the
        // per-blob "real_cluster_main" flash-merge provenance (exact, needs a
        // pctree saved with save_real_cluster_id; falls back to the proxy when
        // the array is absent).
        // beam_window_only (C++ default false; key omitted when off =>
        // byte-identical pre-knob config): evaluate ONLY the beam-coincident
        // bundle, i.e. mains whose matched flash time (cluster_t0) is in the
        // SAME [beam_window_low, beam_window_high) window the in-beam
        // protection already uses.  Same gate as the steiner stage and
        // tagger_check_{stm,fc}; verdicts on surviving mains are unchanged.
        tagger_check_tgm(name="", fiducial="", fv_tolerance=[], beam_window_low=0, beam_window_high=0, beam_window_only=false, length_limit_frac=0.45, enable_case_b=true, require_in_scope=false, check_neutrino_candidate=false, require_chord_charge=false, chord_support_radius=null, chord_max_gap=null, chord_charge_mode="chord", component_extremes=false, component_min_length=null, component_rescue=false, rescue_chord_check=false, main_component_pairs=false, main_component_mode="path", interior_fv_tolerance=[], evaluate_demoted_mains=false, exempt_demoted_main_pairs=false) :: {
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
              + (if beam_window_only then { beam_window_only: true } else {})
              + (if require_in_scope then { require_in_scope: true } else {})
              + (if evaluate_demoted_mains then { evaluate_demoted_mains: true } else {})
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
              + (if rescue_chord_check then { rescue_chord_check: true } else {})
              + (if main_component_pairs then { main_component_pairs: true } else {})
              // Key only when the guard is ON and the mode is non-default, so
              // a runner may pass the mode unconditionally without disturbing
              // the knob-off compiled config.
              + (if main_component_pairs && main_component_mode != "path" then { main_component_mode: main_component_mode } else {})
              // C++ default false. Key omitted when off => byte-identical.
              // Exempts a cluster carrying flag_demoted_main from
              // main_pair_rejects (doc pr/25, SBND evt 320029): its own
              // real_cluster_main/path-component provenance is all-zero by
              // construction, so with main_component_pairs on this guard
              // rejected every demoted-main pair unconditionally.
              + (if exempt_demoted_main_pairs then { exempt_demoted_main_pairs: true } else {})
              // C++ default empty (= use fv_tolerance). Key omitted when
              // empty => byte-identical pre-knob config.  Separate tolerance
              // for the CASE-A interior-support tests only (doc 32 caveat:
              // endpoint-only widening of the downstream-z inset).
              + (if std.length(interior_fv_tolerance) > 0 then { interior_fv_tolerance: interior_fv_tolerance } else {}),
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
        // beam_window_only (C++ default false; keys omitted when off =>
        // byte-identical pre-knob config): evaluate ONLY the beam-coincident
        // bundle, i.e. mains whose matched flash time (cluster_t0) is in
        // [beam_window_low, beam_window_high).  Same gate as the steiner stage
        // and tagger_check_{tgm,stm}.
        // evaluate_demoted_mains (C++ default false; key omitted when off =>
        // byte-identical pre-knob config): also evaluate a cluster carrying
        // flag_demoted_main -- a split part that was ITSELF a matched Q/L bundle
        // main before the flash-time merge demoted it (unmerge_bundle
        // restore_demoted_mains, doc pr/20 Part I P2/P3).  Today no cosmic
        // tagger ever looks at one.  Needs that knob upstream, or nothing
        // carries the flag and this is inert.
        tagger_check_fc(name="", fiducial="", fv_tolerance=[], require_in_scope=false,
                        beam_window_only=false, beam_window_low=0, beam_window_high=0,
                        evaluate_demoted_mains=false) :: {
            type: "TaggerCheckFC",
            name: prefix + name,
            data: {
                grouping: "live",
            } + dv_cfg + pcts_cfg
              + (if fiducial == "" then {} else { fiducial: fiducial, fv_tolerance: fv_tolerance })
              + (if require_in_scope then { require_in_scope: true } else {})
              + (if evaluate_demoted_mains then { evaluate_demoted_mains: true } else {})
              + (if beam_window_only then {
                     beam_window_only: true,
                     beam_window_low: beam_window_low,
                     beam_window_high: beam_window_high,
                 } else {}),
        },

        tagger_check_neutrino(name="", trackfitting_config_file="", particle_dataset="", recombination_model="", perf=false, dl_weights="", dQdx_scale=0.1, dQdx_offset=-1000.0, clus_geom_helper="", dl_vtx_rerank=true, dl_vtx_top_k=5, dl_vtx_min_accept_score=4.0, dl_vtx_score_scale=1000.0, beam_window_low=0, beam_window_high=0, nu_skip_cosmic=false, nu_skip_cosmic_bundle=false, nu_skip_cosmic_bundle_min_length=0, dir_weak_use_score=false, mip_dqdx=null, mip_dqdx_median=null, proton_dir_vote=false, proton_dir_score_max=null, proton_dir_asym_min=null, endpoint_trim_retry=false, fit_vertex_min_seg_length=null, cathode_x=null, cathode_kink_xcut=null, shower_topo_demote_len=null, iso_endpoint=false, iso_endpoint_min_length=null, iso_endpoint_max_xext=null, iso_endpoint_xext_frac=null, iso_endpoint_xext_quantile=null, iso_endpoint_tube_radius=null, iso_endpoint_min_aspect=null, v3_extension_guard=false, v3_extension_min_gain=null, cosmic_y_top_main=null, cosmic_y_top_strict=null, cosmic_y_top_loose=null, cosmic_y_small_piece=null, vertex_z_prior_scale=null, ssm_target_dir=null, ssm_absorber_dir=null, kine_fudge_factor=null, kine_recom_factor=null, kine_shower_fudge_factor=null, kine_shower_recom_factor=null, kine_proton_recom_factor=null, kine_plane_weights=null, kine_plane_asym_switch=null, kine_w_value=null, kine_shower_pdg_live=false, muon_dqdx_curve=null, sp_dedx_use_recomb_model=false, sp_mean_dedx_cut=null, dl_vtx_cut=null, skip_cosmic_companions=false, cosmic_companion_min_length=null, sp_photon_flag=false, fit_exclusion=false, graph_endpoint_strict=false, graph_endpoint_tol=null, oov_prototype_parity=false, first_seg_local_pca=null, other_seg_relaxed_accept=null, shower_topo_proto_dir=false, vertex_dir_use_fit_point=false, shower_traj_recheck_parity=false, main_vertex_require_descriptor=false, main_vertex_candidate_flag=false, cont_muon_dir3_30cm=false, track_comp_empty_abstain=false, shower_topo_reset=false, reclass_preserve_4mom=false, dir_track_median_local=false, examine_showers_vertex_by_index=false, fiducial=null, fv_tolerance=[], sp_sce_correction=false, tagger_ordered_segment_sets=false, stem_endpoint_wcpt_parity=false, broken_muon_cluster_id_count=false, neutrino_type_bitmask=false, daughter_count_proto_main_vertex=false, daughter_count_proto_examine_showers=false, shower_pdg_from_start_segment=false, shower_pdg_from_shower_type=false, shower_pdg_exact_muon_test=false, pi0_id_shared_allocator=false, shower_flag_pdg_electron=false, shower_less_id_tiebreak=false) :: {
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
              // C++ default false.  Key omitted when off => byte-identical pre-knob config (uBooNE).
              // When on (beam gate only): skip in-window mains with flag_TGM/flag_STM/lm_flag>0.
              + (if nu_skip_cosmic then { nu_skip_cosmic: true } else {})
              // C++ default false.  Key omitted when off => byte-identical pre-knob config.
              // When on (needs nu_skip_cosmic): the cosmic verdict vetoes the whole flash
              // bundle -- every in-window main sharing matched_flash_gid with a tagged main
              // is skipped, so no PR runs on a bundle that holds a cosmic (docs/pr/3 sec. 8).
              + (if nu_skip_cosmic_bundle then { nu_skip_cosmic_bundle: true } else {})
              // C++ default 0 (cm).  Key omitted when 0 => byte-identical pre-knob config.
              // When > 0 (needs nu_skip_cosmic_bundle): an UNTAGGED in-window main at least
              // this long survives the bundle veto -- the taggers examined it and declined to
              // tag it; its cosmic-tagged bundle-mate stays excluded regardless (docs/pr/16
              // design A: keeps vetoing unexamined shards like SBND evt 52195's 1.3 cm main).
              + (if nu_skip_cosmic_bundle_min_length > 0 then { nu_skip_cosmic_bundle_min_length: nu_skip_cosmic_bundle_min_length } else {})
              // skip_cosmic_companions / cosmic_companion_min_length (cm):
              // doc pr/20 Part I P4.  Drop a TGM- or STM-tagged COMPANION from
              // the neutrino's other_clusters, so its charge cannot become a
              // particle in the flow tree or enter kine_reco_Enu.  The length
              // floor keeps a SHORT tagged companion in regardless of verdict.
              // Nothing tags a companion unless the taggers run with
              // evaluate_demoted_mains, so this is inert without P3.  C++
              // defaults false / 0; keys omitted when off => byte-identical.
              + (if skip_cosmic_companions then { skip_cosmic_companions: true } else {})
              + (if cosmic_companion_min_length != null then { cosmic_companion_min_length: cosmic_companion_min_length } else {})
              // sp_photon_flag: store singlephoton_tagger()'s verdict in
              // TaggerInfo::photon_flag, the way prototype NeutrinoID.cxx:271
              // does.  The port already runs that tagger and fills its ~90
              // shw_sp_* BDT features; only the verdict was dropped, leaving
              // the uBooNE tagger ntuple's photon_flag branch a constant 0
              // (sbnd_xin/docs/pr/26 sec. 8.2).  Nothing in the chain reads
              // the field, so this changes one output branch and no
              // reconstruction.  C++ default false; key omitted when off =>
              // byte-identical pre-knob config.
              + (if sp_photon_flag then { sp_photon_flag: true } else {})
              // C++ default false.  Key omitted when off => byte-identical pre-knob config (uBooNE).
              // When on: direction-weakness reads use segment_is_dir_weak() (score thresholds),
              // the faithful port of prototype ProtoSegment::is_dir_weak() -- see sbnd_xin/docs/pr/6.
              + (if dir_weak_use_score then { dir_weak_use_score: true } else {})
              // MIP dQ/dx scales in e/cm.  C++ defaults 50000 (flat-template role) /
              // 43000 (median-threshold role) = the uBooNE legacy values; keys omitted
              // when null => byte-identical pre-knob config.  SBND: 56000 / 48000
              // (owner 2026-07-30; docs pr/7 sec 5 + pr/8).
              + (if mip_dqdx != null then { mip_dqdx: mip_dqdx } else {})
              + (if mip_dqdx_median != null then { mip_dqdx_median: mip_dqdx_median } else {})
              // Proton-template direction vote (doc pr/8).  C++ default false = legacy
              // muon-vs-flat-only direction; keys omitted when off => byte-identical.
              // Threshold keys emitted only when explicitly set (C++ defaults 0.25 / 1.3
              // are initial values pending the pr/8 sec 6 calibration).
              + (if proton_dir_vote then { proton_dir_vote: true } else {})
              + (if proton_dir_score_max != null then { proton_dir_score_max: proton_dir_score_max } else {})
              + (if proton_dir_asym_min != null then { proton_dir_asym_min: proton_dir_asym_min } else {})
              // Endpoint-trim retry (doc pr/9 sec 6 F1).  C++ default false = legacy
              // abstention; key omitted when off => byte-identical.  When on: on
              // double abstention, retry the dQ/dx PID once with exactly 1 sample
              // excluded at each orientation's hypothesized stopping end (runs
              // before the proton_dir_vote fallback).
              + (if endpoint_trim_retry then { endpoint_trim_retry: true } else {})
              // fit_vertex short-segment exclusion (doc pr/9 sec 11 F3c), value in cm.
              // C++ default 0 = all segments enter the vertex position fit; key omitted
              // when null => byte-identical.  When set: segments with wcpt-path length
              // below the cut do not count as vertex-fit legs -- >=3 surviving legs
              // fit on the survivors, <=2 skip the fit (the two-leg position was
              // already fit; refitting on stub-contaminated re-tracked clouds
              // reproduces the drag).  SBND: 1.0 cm (owner 2026-07-30, deliberate
              // prototype divergence).
              + (if fit_vertex_min_seg_length != null then { fit_vertex_min_seg_length: fit_vertex_min_seg_length } else {})
              // Cathode kink veto (doc sbnd_xin/docs/pr/20 Part II B0), both in cm.
              // segment_search_kink skips candidate fit points within
              // cathode_kink_xcut of the cathode plane, where the ~2 cm transverse
              // cathode mismatch fakes a turn while the para_angle guard is wide
              // open (the crossing is drift-x dominated).  C++ default 0 => no
              // point is ever skipped => keys omitted when null are byte-identical.
              + (if cathode_x != null then { cathode_x: cathode_x } else {})
              + (if cathode_kink_xcut != null then { cathode_kink_xcut: cathode_kink_xcut } else {})
              // shower_topo_demote_len (cm, doc pr/25 sec 3): a segment the
              // topology test would flag kShowerTopology whose geometric
              // length exceeds this is demoted to a track, so it gets real
              // track PID instead of the hard-coded pdg=11/score=100.  The
              // test's only measurement axis is drift-aligned for an
              // isochronous segment, where it reads a 0.313 cm slice lattice
              // against a 0.4 cm cut.  C++ default 0 => the guard never fires
              // => key omitted when null is byte-identical.
              + (if shower_topo_demote_len != null then { shower_topo_demote_len: shower_topo_demote_len } else {})
              // ---- doc sbnd_xin/docs/pr/30 §11 port-fidelity knobs -------------
              // fit_exclusion (P1): pass flag_exclusion=true to the 27 knobbed
              // do_multi_tracking call sites, as 28 of the 30 live prototype
              // sites do.  With it on, form_map_graph calls update_association
              // and strips from each segment's 2-D associations the (wire,tick)
              // cells belonging to OTHER segments.  break_segments' two sites
              // and the single-segment local fitter are never knobbed -- they
              // already match the prototype.  C++ default false.
              + (if fit_exclusion then { fit_exclusion: true } else {})
              // graph_endpoint_strict (P8): REFUSE a PR::add_segment whose
              // vertices do not sit within graph_endpoint_tol of the segment's
              // two wcpt ends, as the prototype's add_proto_connection does.
              // The WARN and the counter are unconditional in C++; only the
              // refusal is gated.  C++ defaults false / 0.3 cm.
              + (if graph_endpoint_strict then { graph_endpoint_strict: true } else {})
              + (if graph_endpoint_tol != null then { graph_endpoint_tol: graph_endpoint_tol } else {})
              // oov_prototype_parity (F2, was P9): make a point outside every
              // TPC vote the way the prototype's own helper answers for it, at
              // all three sites -- bad (not connected) in
              // modify_segment_isochronous, not-dead in examine_vertices_1p,
              // unique (segment kept) in examine_vertices_3.  Today all three
              // vote the opposite way.  C++ default false.
              + (if oov_prototype_parity then { oov_prototype_parity: true } else {})
              // shower_topo_proto_dir (doc pr/31 sec 11, F2 was P2): skip the
              // stage-3 segment_determine_shower_direction call so a topology
              // shower keeps the direction segment_is_shower_topology set --
              // the prototype's state, since its determine_dir_shower_topology
              // does not touch flag_dir and its determine_shower_direction()
              // runs only in stage 4.  C++ default false.  Key omitted when
              // off => byte-identical pre-pr/31 config.
              + (if shower_topo_proto_dir then { shower_topo_proto_dir: true } else {})
              // doc sbnd_xin/docs/pr/32 sec 11 -- the four kept findings of the
              // stage-4 (neutrino vertex identification) port audit.  Every
              // C++ default is false = today's path, and every key is omitted
              // when off => byte-identical pre-pr/32 config.
              //
              // vertex_dir_use_fit_point (F1, was P1): measure the
              // calc_conflict_maps direction vectors and the all-showers PCA
              // projection / z tie-breaks / Steiner path endpoints from the
              // CONTINUOUS fit, as the prototype's get_fit_pt() does, instead
              // of from the fit snapped to the nearest Steiner node.  Eleven
              // expressions; NOT byte-identical when on.
              + (if vertex_dir_use_fit_point then { vertex_dir_use_fit_point: true } else {})
              // shower_traj_recheck_parity (F2, was P3): restore the
              // prototype's improve_vertex shower-trajectory recheck -- outer
              // gates read the STORED kShowerTrajectory flag, inner test
              // recomputes at 10 cm with the mip_dqdx scale, and
              // segment_is_shower_trajectory re-caches the flag (clearing it
              // when the test says no) the way the prototype does.  The three
              // move together: fixing only the inner parameters makes the
              // block dead code.
              + (if shower_traj_recheck_parity then { shower_traj_recheck_parity: true } else {})
              // main_vertex_require_descriptor (F3, was P7): drop
              // invalid-descriptor candidates before compare_main_vertices
              // scores, so the min_z scan, the fiducial term, the conflict
              // penalty and the argmax see the same candidate set as the two
              // blocks that already guard.  Expected byte-identical (the path
              // looks unreachable) -- the drop is counted, not assumed.
              + (if main_vertex_require_descriptor then { main_vertex_require_descriptor: true } else {})
              // main_vertex_candidate_flag (F4, was P12): set
              // VertexFlags::kMainCandidate on each per-cluster main-vertex
              // candidate, the prototype's map_cluster_main_candidate_vertices.
              // DIAGNOSTIC ONLY -- only PrDisplayDump reads it.
              + (if main_vertex_candidate_flag then { main_vertex_candidate_flag: true } else {})
              // doc sbnd_xin/docs/pr/31 sec 12 -- the sec 10.12 port-fidelity
              // round (topology/PID/direction audit survivors).  Every C++
              // default is false = today's path, and every key is omitted when
              // off => byte-identical pre-pr/31-sec-12 config.
              //
              // cont_muon_dir3_30cm (F5, was P6): find_cont_muon_segment_nue's
              // hoisted dir3 always at 30 cm, as the prototype computes it,
              // instead of falling back to the 15 cm dir1 for a short
              // reference segment.
              + (if cont_muon_dir3_30cm then { cont_muon_dir3_30cm: true } else {})
              // track_comp_empty_abstain (F6, was P7): an empty dQ/dx
              // comparison window ABSTAINS from the direction gate (the
              // prototype's degenerate answer, verified by execution) instead
              // of declaring the orientation confirmed.
              + (if track_comp_empty_abstain then { track_comp_empty_abstain: true } else {})
              // shower_topo_reset (F3, was P13): segment_is_shower_topology
              // clears kShowerTopology and dirsign at entry, before its early
              // returns, as the prototype does -- no stale flag survives a
              // re-test, no stale direction survives an empty-cloud return.
              + (if shower_topo_reset then { shower_topo_reset: true } else {})
              // reclass_preserve_4mom (F1, was P1+P3a+P4): the 15
              // reclassification sites preserve the existing 4-momentum and
              // recompute only where the prototype's get_particle_4mom(3)>0
              // guard passes.  Moves kine_reco_Enu directly.
              + (if reclass_preserve_4mom then { reclass_preserve_4mom: true } else {})
              // dir_track_median_local (F4, was P8): determine_dir_track's
              // median dQ/dx over the SAME local vector the PID receives,
              // as the prototype's nth_element does, instead of the filtered
              // helper rebuild.
              + (if dir_track_median_local then { dir_track_median_local: true } else {})
              // examine_showers_vertex_by_index (F7, was P5): order
              // examine_all_showers' vertex pair by graph index before the
              // asymmetric 165/150-degree branches.  DELIBERATELY DORMANT --
              // stays off pending pr/30 F4's find_vertices adjudication.
              + (if examine_showers_vertex_by_index then { examine_showers_vertex_by_index: true } else {})
              // first_seg_local_pca (P2) and other_seg_relaxed_accept (P4) are
              // the two knobs whose C++ default is TRUE, because the behaviour
              // they gate is already production.  null => key omitted => the
              // C++ default => byte-identical; pass false to restore the
              // prototype's narrower behaviour for measurement.
              + (if first_seg_local_pca != null then { first_seg_local_pca: first_seg_local_pca } else {})
              + (if other_seg_relaxed_accept != null then { other_seg_relaxed_accept: other_seg_relaxed_accept } else {})
              // Isochronous first-segment endpoint finding (doc pr/24 round 2,
              // SBND evt 271851): for a long cluster whose quantile-trimmed
              // drift-x extent is small (a filled 2-D sheet), pick the first
              // PR segment's endpoints from the sheet's principal axis instead
              // of the wire-footprint boundary metric (which degenerates to
              // sheet corners), and skip the local-PCA endpoint refinement on
              // that branch.  C++ defaults false / 40 cm / 25 cm / 0.35 /
              // 0.02; keys omitted when off/null => byte-identical.
              + (if iso_endpoint then { iso_endpoint: true } else {})
              + (if iso_endpoint_min_length != null then { iso_endpoint_min_length: iso_endpoint_min_length } else {})
              + (if iso_endpoint_max_xext != null then { iso_endpoint_max_xext: iso_endpoint_max_xext } else {})
              + (if iso_endpoint_xext_frac != null then { iso_endpoint_xext_frac: iso_endpoint_xext_frac } else {})
              + (if iso_endpoint_xext_quantile != null then { iso_endpoint_xext_quantile: iso_endpoint_xext_quantile } else {})
              // C++ defaults 4 cm / 0.12 (doc pr/24 round 3).  Keys omitted when
              // null => byte-identical compiled config.
              + (if iso_endpoint_tube_radius != null then { iso_endpoint_tube_radius: iso_endpoint_tube_radius } else {})
              + (if iso_endpoint_min_aspect != null then { iso_endpoint_min_aspect: iso_endpoint_min_aspect } else {})
              // doc pr/24 round 5.  examine_vertices_3's get_local_extension
              // recovery step can retract a track endpoint instead of
              // extending it; C++ default false / -1.0 cm => unconditional
              // accept, byte-identical when off.
              + (if v3_extension_guard then { v3_extension_guard: true } else {})
              + (if v3_extension_min_gain != null then { v3_extension_min_gain: v3_extension_min_gain } else {})
              // Detector-extent literals (docs/pr/2 sec 2e(iv)), all in cm.  C++
              // defaults = the uBooNE prototype values, so keys omitted when null =>
              // byte-identical pre-knob config.
              //
              // cosmic_y_*: cosmic_tagger()'s four "reaches the top of the detector"
              // tests.  uBooNE's active volume tops out at y = +117 cm
              // (prototype_base/pid/apps/wire-cell-prod-nue.cxx:417), so the legacy
              // 100 / 102 / 80 / 50 cm are 17 / 15 / 37 / 67 cm BELOW THE TOP FACE --
              // an offset a downward cosmic entering the top satisfies whatever the
              // detector height is, hence a taller detector moves the cuts up rather
              // than scaling them.  Roles: top_main = the MAIN cluster's own highest
              // point (relaxes the vertical-angle cut 20 -> 30 deg); top_strict =
              // event highest point in the single-cosmic branch; top_loose = event
              // highest point, global gate on the whole decision; small_piece = PCA
              // centre of a <3 cm cluster counted as cosmic debris.  SBND (top
              // y = +200 cm): 183 / 185 / 163 / 133.
              + (if cosmic_y_top_main != null then { cosmic_y_top_main: cosmic_y_top_main } else {})
              + (if cosmic_y_top_strict != null then { cosmic_y_top_strict: cosmic_y_top_strict } else {})
              + (if cosmic_y_top_loose != null then { cosmic_y_top_loose: cosmic_y_top_loose } else {})
              + (if cosmic_y_small_piece != null then { cosmic_y_small_piece: cosmic_y_small_piece } else {})
              // vertex_z_prior_scale (cm): denominator of the upstream-z penalty
              // (z - min_z)/scale in compare_main_vertices and
              // compare_main_vertices_global, competing with the +0.25-per-track
              // bonuses.  uBooNE's 200 cm spans ~5.2 penalty units over its 1037 cm
              // detector; a shorter detector keeps the same per-cm trade-off but
              // loses dynamic range, which bites hardest in the _global comparison
              // (candidates from DIFFERENT clusters of the beam bundle).  SBND:
              // 100 cm ~ 200 x 501/1037.  Must be > 0.
              + (if vertex_z_prior_scale != null then { vertex_z_prior_scale: vertex_z_prior_scale } else {})
              // SSM beam-line reference directions [x,y,z] in the detector frame,
              // feeding the 8 ssm_*_angle_{target,absorber} BDT features.  C++
              // defaults = the prototype's uBooNE BNB-target (0.46,0.05,0.885) and
              // NuMI-absorber (0.33,0.75,-0.59) vectors; keys omitted when null =>
              // byte-identical pre-knob config.  SBND has no value for either yet
              // (docs/pr/2 sec 2e(i)) -- the numbers are merely reachable now, not
              // fixed.  Note they are deliberately NOT unit vectors.
              + (if ssm_target_dir != null then { ssm_target_dir: ssm_target_dir } else {})
              + (if ssm_absorber_dir != null then { ssm_absorber_dir: ssm_absorber_dir } else {})
              // Max distance (mm; C++ default 25.0 = 2.5 cm) from the DL SCN
              // prediction to accept a candidate vertex.  Coupled to the
              // uBooNE-trained net (docs/pr/2 sec 7.4) -- threaded here so it
              // has a proper configuration path.  Key omitted when null =>
              // byte-identical pre-knob config.
              + (if dl_vtx_cut != null then { dl_vtx_cut: dl_vtx_cut } else {})
              // Charge -> kinetic-energy calibration of NeutrinoEnergyReco
              // (docs/pr/2 sec 2e(iii)).  C++ defaults are the uBooNE-tuned
              // literals these keys replaced -- 0.95/0.7 (track), 0.8/0.5
              // (shower-flagged), 0.35 (proton recombination only), plane
              // weights [0.25,0.25,1.0] for [U,V,W], asymmetry switch 0.04,
              // W-value 23.6 eV -- so keys omitted when null => byte-identical
              // pre-knob config.  None has been re-derived for SBND: the
              // recombination/fudge pair is field- and calibration-dependent
              // and the plane weights encode uBooNE induction-plane quality.
              // Exposed so a calibration can move them without a code change;
              // no detector sets one yet.
              //
              // E = sum_p(w_p Q_p)/sum(w) / recom / fudge * kine_w_value * 1e-6 MeV,
              // with the 3-plane average replaced by the (median,min) pair when
              // the two largest plane charges differ by more than the switch.
              // A single zero plane weight is accepted; an all-zero triple is
              // rejected by C++ with a warning.
              + (if kine_fudge_factor != null then { kine_fudge_factor: kine_fudge_factor } else {})
              + (if kine_recom_factor != null then { kine_recom_factor: kine_recom_factor } else {})
              + (if kine_shower_fudge_factor != null then { kine_shower_fudge_factor: kine_shower_fudge_factor } else {})
              + (if kine_shower_recom_factor != null then { kine_shower_recom_factor: kine_shower_recom_factor } else {})
              + (if kine_proton_recom_factor != null then { kine_proton_recom_factor: kine_proton_recom_factor } else {})
              + (if kine_plane_weights != null then { kine_plane_weights: kine_plane_weights } else {})
              + (if kine_plane_asym_switch != null then { kine_plane_asym_switch: kine_plane_asym_switch } else {})
              + (if kine_w_value != null then { kine_w_value: kine_w_value } else {})
              // doc pr/35 sec 10.2 (F1 = P1+P8): read the shower PDG live from
              // the start segment at the four fill_kine_tree sites (prototype
              // kine.h:53 :67 :175 :187) instead of Shower's cached field,
              // whose refresh path is incomplete.  C++ default false.  Key
              // omitted when off => byte-identical pre-knob config.
              + (if kine_shower_pdg_live then { kine_shower_pdg_live: true } else {})
              // Muon median-dQ/dx-vs-length envelope [c0, c1, pivot_cm, power]:
              // cut = c0 + c1*(pivot/L)^power, a multiple of mip_dqdx_median,
              // used at NINE tagger sites (numu x2, vertex-finder, nue x4, ssm,
              // cosmic -- docs/pr/2 sec 2e(iv)).  C++ default
              // [0.8866, 0.9533, 18, 0.4234] = the prototype's empirical uBooNE
              // stopping-muon refit; key omitted when null => byte-identical
              // pre-knob config.  pivot and power must be > 0 (C++ rejects with
              // a warning).
              + (if muon_dqdx_curve != null then { muon_dqdx_curve: muon_dqdx_curve } else {})
              // Single-photon stem dE/dx conversion (docs/pr/2 sec 2e(i) third
              // correctness item).  C++ default false = the inline float inverse
              // Modified-Box frozen at uBooNE's 0.273 kV/cm; key omitted when off
              // => byte-identical.  When on: shw_sp_vec_{median,mean}_dedx go
              // through the component's configured recombination_model instead.
              + (if sp_dedx_use_recomb_model then { sp_dedx_use_recomb_model: true } else {})
              // The hard shw_sp_vec_mean_dedx cut in MeV/cm.  C++ default 2.3 =
              // the legacy literal (compared as float => bit-identical); key
              // omitted when null.  Tuned against the INLINE (uBooNE-field)
              // dE/dx scale -- retune together with the knob above.
              + (if sp_mean_dedx_cut != null then { sp_mean_dedx_cut: sp_mean_dedx_cut } else {})
              // ---- doc sbnd_xin/docs/pr/36 sec 10 tagger-stage knobs -------
              // F1 (= P1): fiducial + margins for the match_isFC recompute --
              // the SAME objects tagger_check_{stm,tgm,fc} receive, restoring
              // one containment definition across the stage.  C++: an absent
              // "fiducial" key keeps the historical FiducialUtils fallback
              // (cluster_fc_check's nullptr path, documented bit-for-bit) =>
              // byte-identical pre-knob config.  Same key pair as
              // tagger_check_stm above.
              + (if fiducial != null then { fiducial: fiducial } else {})
              + (if std.length(fv_tolerance) > 0 then { fv_tolerance: fv_tolerance } else {})
              // F3 (= P2): single-photon SCE correction gate.  A SEPARATE
              // bool from clus_geom_helper so the two consumers of that key
              // (kine + single-photon) stay independently gateable.  C++
              // default false.  Key omitted when off => byte-identical.
              + (if sp_sce_correction then { sp_sce_correction: true } else {})
              // F4 (= P3+P5): iterate the three tagger accumulation sets in
              // graph-index order (M4 house-rule determinism fix; prototype
              // n/a -- it is address-ordered too).  C++ default false.
              + (if tagger_ordered_segment_sets then { tagger_ordered_segment_sets: true } else {})
              // F5 (= P6): prototype wcpt-identity stem-endpoint rule at the
              // 18 seg_endpoint_near sites.  C++ default false.
              + (if stem_endpoint_wcpt_parity then { stem_endpoint_wcpt_parity: true } else {})
              // F6 (= P8): broken_muon_id counts distinct cluster IDS
              // (prototype) instead of pointers.  C++ default false.
              + (if broken_muon_cluster_id_count then { broken_muon_cluster_id_count: true } else {})
              // F7 (= P4): compute the prototype's neutrino_type verdict
              // bitmask.  The matching T_tagger branch is booked by
              // tagger_output under the same key.  C++ default false.
              + (if neutrino_type_bitmask then { neutrino_type_bitmask: true } else {})
              // ---- doc sbnd_xin/docs/pr/33 §10 EM-shower-clustering knobs.
              // All C++ default false.  Key omitted when off => byte-identical
              // pre-knob config.
              // F1 (= P1): restore the prototype's calculate_num_daughter_tracks
              // callee at the two sites that call _showers.  Two knobs: the
              // sites err in opposite directions.
              + (if daughter_count_proto_main_vertex then { daughter_count_proto_main_vertex: true } else {})
              + (if daughter_count_proto_examine_showers then { daughter_count_proto_examine_showers: true } else {})
              // F2 (= P2): read the PDG off the object the prototype reads.
              // NOTE: prototype parity at the :170 site needs BOTH
              // shower_pdg_from_start_segment AND shower_pdg_exact_muon_test
              // on; from_start_segment alone is neither tree's behavior there.
              + (if shower_pdg_from_start_segment then { shower_pdg_from_start_segment: true } else {})
              + (if shower_pdg_from_shower_type then { shower_pdg_from_shower_type: true } else {})
              + (if shower_pdg_exact_muon_test then { shower_pdg_exact_muon_test: true } else {})
              // F3 (= P3): the two pi0 finders share one id allocation
              // stream (prototype member semantics), preventing pio_id
              // collisions in the tagger pi0 block and the Bee mc.json
              // grouping.  Scoped to the finders only.
              + (if pi0_id_shared_allocator then { pi0_id_shared_allocator: true } else {})
              // F4 (= P6): is_shower at the center-point site gains the
              // prototype's missing abs(pdg)==11 disjunct.
              + (if shower_flag_pdg_electron then { shower_flag_pdg_electron: true } else {})
              // F5 (= P12): shower_less same-index tie-break by stable
              // shower id instead of pointer address (house-rule determinism
              // fix; prototype n/a -- it is address-ordered everywhere).
              + (if shower_less_id_tiebreak then { shower_less_id_tiebreak: true } else {}),
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
        tagger_output(name="", output_filename="tracking_proj.root", neutrino_type_bitmask=false) :: {
            type: "UbooneTaggerOutputVisitor",
            name: prefix + name,
            data: {
                grouping: "live",
                output_filename: output_filename,
            }
              // doc sbnd_xin/docs/pr/36 sec 10.8 (F7 = P4): book the
              // neutrino_type/I branch (prototype
              // wire-cell-prod-nue-port.cxx:1486).  Same key as
              // tagger_check_neutrino's computing knob.  C++ default false =
              // branch not booked; key omitted when off => byte-identical
              // pre-knob config AND schema.
              + (if neutrino_type_bitmask then { neutrino_type_bitmask: true } else {}),
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

        // SBND cathode BUNDLE rescue (isolated patch for the flash-reco
        // absorbing-window defect; sbnd_xin/docs/pr/14).  Joins a cathode
        // crosser whose two halves sit in DIFFERENT flash bundles (which
        // cathode_connect's flash gate correctly refuses): finds a beam-window
        // half with a tip at the cathode plus an out-of-beam partner within
        // [-rescue_t0_early, +rescue_t0_late], decides which flash the joined
        // crosser belongs to from the bundle structure, merges the pair and
        // re-stamps the merged cluster's flash identity.  Place AFTER
        // cathode_connect and BEFORE examine_bundles.  Retire when the light
        // reconstruction fixes the absorbing window.
        //
        // C++ defaults leave the beam window EMPTY (low == high == 0), so the
        // component is a no-op unless beam_window_low/high are set.  The
        // geometric knobs mirror cathode_connect (same names, same C++
        // defaults; the PDVD-specific crosser_* relaxations are not offered).
        cathode_bundle_rescue(name="", beam_window_low=0.0, beam_window_high=0.0,
                              rescue_t0_early=8*wc.us, rescue_t0_late=13*wc.us,
                              require_far_out_of_beam=true,
                              drift_cut=5*wc.cm, dis_cut=5*wc.cm, max_dis=25*wc.cm,
                              angle_cut=10.0, conn_far_cut=30.0, cathode_x=0.0,
                              cathode_x_cut=3.5*wc.cm, hough_radius=20*wc.cm,
                              min_length=10*wc.cm, min_length_short=null,
                              short_dir_len=null, conn_short_cut=30.0,
                              tip_touch_cut=null, tip_touch_angle_cut=null,
                              cathode_band_dis=null,
                              rescue_unmatched=false,
                              unmatched_min_length=null, unmatched_min_npts=null,
                              adopt_nu_fragments=false, adopt_dis=null,
                              adopt_xcut=null, adopt_frag_max_length=null,
                              adopt_min_npts=null, adopt_beam_min_length=null) :: {
            type: "ClusteringCathodeBundleRescue",
            name: prefix+name,
            data: dv_cfg + pcts_cfg + scope_cfg + {
                beam_window_low: beam_window_low,
                beam_window_high: beam_window_high,
                rescue_t0_early: rescue_t0_early,
                rescue_t0_late: rescue_t0_late,
                require_far_out_of_beam: require_far_out_of_beam,
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
                // null => C++ defaults short_dir_len to 0 (both-long/short-stub branches OFF)
                [if short_dir_len != null then "short_dir_len"]: short_dir_len,
                conn_short_cut: conn_short_cut,
                // null => C++ defaults tip_touch_cut to 0 (tip-touch relaxation OFF)
                [if tip_touch_cut != null then "tip_touch_cut"]: tip_touch_cut,
                // null => C++ defaults tip_touch_angle_cut to angle_cut (local-Hough fallback OFF)
                [if tip_touch_angle_cut != null then "tip_touch_angle_cut"]: tip_touch_angle_cut,
                // null => C++ defaults cathode_band_dis to 0 (near-cathode retry OFF)
                [if cathode_band_dis != null then "cathode_band_dis"]: cathode_band_dis,
                // Unmatched-cluster adoption pass (sbnd_xin/docs/pr/17): a
                // flashless cluster that geometrically continues a beam-window
                // cluster across the cathode is merged into the beam bundle.
                // C++ default false.  Key omitted when off => byte-identical
                // pre-knob config.
                [if rescue_unmatched then "rescue_unmatched"]: true,
                // null => C++ default floors (30 cm / 200 pts) on the adopted cluster.
                [if unmatched_min_length != null then "unmatched_min_length"]: unmatched_min_length,
                [if unmatched_min_npts != null then "unmatched_min_npts"]: unmatched_min_npts,
                // Near-cathode fragment adoption pass (sbnd_xin/docs/pr/19):
                // a small flashless fragment reaching within adopt_xcut of the
                // cathode is merged into a beam-window cluster when its raw
                // closest approach under the beam-T0 hypothesis is within
                // adopt_dis.  Companion of clustering_isolated's
                // cathode_guard_xcut.  C++ default false.  Key omitted when
                // off => byte-identical pre-knob config.
                [if adopt_nu_fragments then "adopt_nu_fragments"]: true,
                // nulls => C++ defaults (10 cm / 30 cm / 60 cm / 5 pts / 10 cm).
                [if adopt_dis != null then "adopt_dis"]: adopt_dis,
                [if adopt_xcut != null then "adopt_xcut"]: adopt_xcut,
                [if adopt_frag_max_length != null then "adopt_frag_max_length"]: adopt_frag_max_length,
                [if adopt_min_npts != null then "adopt_min_npts"]: adopt_min_npts,
                [if adopt_beam_min_length != null then "adopt_beam_min_length"]: adopt_beam_min_length,
            },
            uses: [detector_volumes, pc_transforms],
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
                 collinear_global_merge=false, vertex_veto=false) :: {
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
                // Un-split a neutrino-vertex "V": both dominant separation
                // pieces END at their mutual closest approach, so the cluster
                // is one interaction, not two crossing cosmics (doc pr/15,
                // SBND run 18255 evt 56463).  C++ default false.  Key omitted
                // when false => byte-identical pre-fix config.
                [if vertex_veto then 'vertex_veto']: vertex_veto,
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
        // save_assoc_id (C++ default false; key omitted when off => byte-identical
        // pre-knob config): also record this pass's main + associated partition
        // into the dedicated per-blob pair "assoc_cluster_id" /
        // "assoc_cluster_main" ("perblob").  merge_clusters carries that pair
        // across every later merge and clustering_switch_scope carves it through
        // the scope rebuild, so ClusteringUnmergeBundle(id_aname=..., main_aname=...)
        // can undo the grouping before the taggers run -- the prototype's
        // main_cluster + additional_clusters layout.  Needs
        // MultiAlgBlobClustering save_assoc_cluster_id=true to survive the pctree
        // tarball.  Doc 52.
        // cathode_guard_xcut (default null => key omitted => guard OFF,
        // byte-identical pre-knob config; sbnd_xin/docs/pr/19): decline the
        // angle-less small->big absorb for a small cluster that reaches within
        // this distance of the cathode plane AND is farther from every big
        // cluster than from the cathode -- it plausibly belongs to activity in
        // the other drift volume, invisible to this per-APA pass.  cathode_x
        // (default null => C++ 0) is the cathode plane position in this
        // pass's raw frame.
        isolated(name="", use_flash_t0=false, flash_t0_window=80*wc.ns, length_cut=null, range_cut=null,
                 save_assoc_id=false, cathode_guard_xcut=null, cathode_x=null,
                 cathode_guard_dis_floor=null) :: {
            type: "ClusteringIsolated",
            name: prefix+name,
            data: {
                use_flash_t0: use_flash_t0,
                flash_t0_window: flash_t0_window,
                [if length_cut != null then 'length_cut']: length_cut,
                [if range_cut != null then 'range_cut']: range_cut,
                [if save_assoc_id then 'save_assoc_id']: true,
                [if cathode_guard_xcut != null then 'cathode_guard_xcut']: cathode_guard_xcut,
                [if cathode_x != null then 'cathode_x']: cathode_x,
                // null => C++ 0 (no floor): the guard declines regardless of how
                // close the big cluster is.  Set to keep nearby (< floor) absorbs.
                [if cathode_guard_dis_floor != null then 'cathode_guard_dis_floor']: cathode_guard_dis_floor,
            } + dv_cfg + scope_cfg,
        },

        // flags_from_longest (default false): on the flash-time merge, take the
        // merged cluster's flags from the same representative member that donates
        // its flash instead of from an arbitrary (last-visited) member, so a
        // matched main cannot lose flag_main_cluster to a co-merged fragment.
        // Key emitted only when true so existing compiled configs stay
        // byte-identical.  See merge_clusters() in clus/inc/.../ClusteringFuncs.h.
        // save_bundle_main_provenance (C++ default false; key omitted when off
        // => byte-identical pre-knob config): on the flash-time merge, also
        // write the per-blob "real_cluster_was_main" array -- 1 on the rows of
        // EVERY member that carried flag_main_cluster on the way in, i.e. the
        // bundle mains this merge demotes.  Distinct from "real_cluster_main",
        // which marks only the single representative member.  Consumed by
        // unmerge_bundle's restore_demoted_mains.  Doc pr/20 Part I P1.
        examine_bundles(name="", graph_name="relaxed", use_flash_t0=false, flash_t0_window=80*wc.ns,
                        flags_from_longest=false, save_bundle_main_provenance=false) :: {
            type: "ClusteringExamineBundles",
            name: prefix+name,
            data: dv_cfg + pcts_cfg + scope_cfg + {
                graph_name: graph_name,
                use_flash_t0: use_flash_t0,
                flash_t0_window: flash_t0_window,
            } + (if flags_from_longest then { flags_from_longest: true } else {})
              + (if save_bundle_main_provenance then { save_bundle_main_provenance: true } else {}),
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
        neutrino(name="", num_try=1, use_flash_t0=false, flash_t0_window=80*wc.ns, protect_iso_band=false, protect_iso_band_xext=null) :: {
            type: "ClusteringNeutrino",
            name: prefix+name,
            data: {
                num_try: num_try,
                use_flash_t0: use_flash_t0,
                flash_t0_window: flash_t0_window,
                [if protect_iso_band then 'protect_iso_band']: protect_iso_band,
                // C++ default 0 (off).  Key omitted when null => byte-identical
                // pre-knob config.  When set (needs protect_iso_band on), a
                // band/non-band merge is refused even on genuine touch if the
                // non-band partner spans more than this in drift x (doc pr/18,
                // SBND evt 10550).
                [if protect_iso_band_xext != null then 'protect_iso_band_xext']: protect_iso_band_xext,
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
        // beam_window_only (C++ default false; keys omitted when off =>
        // byte-identical pre-knob config): build steiner graphs ONLY for the
        // beam-coincident bundle -- the main clusters whose matched flash time
        // (cluster_t0) is in [beam_window_low, beam_window_high) plus the
        // companions sharing their matched_flash_gid.  Pass the same window the
        // taggers get: with tagger_check_{tgm,stm,fc} gated the same way, the
        // clusters that lose their graph are exactly the ones no tagger reads.
        // terminal_wire_tol / terminal_adjacent_slice (doc pr/29 D1, D12) fix
        // the Steiner TERMINAL filter only -- get_extreme_wcps shares the same
        // C++ helper and keeps the exact, no-slack, no-fallback behaviour the
        // prototype gives it.  Both default to the historical toolkit values
        // and their keys are omitted when off => byte-identical pre-knob config.
        steiner(name="", retiler={}, grouping="live", graph="steiner", perf=true, require_beam_flash=true,
                beam_window_only=false, beam_window_low=0, beam_window_high=0, replace=null,
                terminal_wire_tol=0, terminal_adjacent_slice=false,
                edge_charge_forward_dead_mix=false) :: {
            type: "CreateSteinerGraph",
            name: prefix+name,
            data: {
                grouping: grouping,
                graph: graph,
                retiler: wc.tn(retiler),
                perf: perf,
                require_beam_flash: require_beam_flash,
                // C++ default true.  false = keep existing steiner products
                // (the doc pr/23 second pass rebuilds ONLY the clusters
                // protect_bundle purged; replacing an existing graph
                // erase+emplaces the store node and dangles any
                // GraphAlgorithms made in the tagger stage).  Key omitted
                // when null => byte-identical pre-knob config.
                [if replace != null then 'replace']: replace,
                // C++ default 0.  1 = the prototype's one wire of slack on
                // both sides of all three planes in the terminal filter
                // (PR3DCluster_steiner.h:285-290); get_extreme_wcps is
                // unaffected.  Key omitted when 0 => byte-identical pre-knob
                // config.
                [if terminal_wire_tol != 0 then 'terminal_wire_tol']: terminal_wire_tol,
                // C++ default false.  true = the adjacent-slice fallback steps
                // by the face's ticks-per-slice, which is what makes it resolve
                // at all (the time_blob_map key is in ticks).  Key omitted when
                // false => byte-identical pre-knob config.
                [if terminal_adjacent_slice then 'terminal_adjacent_slice']: true,
                // C++ default false.  true = the edge-weight charges honour
                // the disable_dead_mix_cell create_steiner_tree was called
                // with (false in this chain), as the prototype does
                // (PR3DCluster_steiner.h:514,:521); false keeps the dropped
                // argument and its always-true downstream default.  Key
                // omitted when false => byte-identical pre-knob config.
                [if edge_charge_forward_dead_mix then 'edge_charge_forward_dead_mix']: true,
            } + dv_cfg + pcts_cfg
              + (if beam_window_only then {
                     beam_window_only: true,
                     beam_window_low: beam_window_low,
                     beam_window_high: beam_window_high,
                 } else {}),
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
