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
                       skip_convicted=null, open_convicted_bundles=null,
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
                // doc sbnd_xin/docs/pr/94 round 3.  When true a convicted main
                // still OPENS its bundle, so the bundle's other, unconvicted
                // members -- the secondary activity a cosmic-tagged main would
                // otherwise silence -- get the same graph examination every
                // member of an unconvicted bundle gets.  The convicted cluster
                // itself is still never split (the per-member skip_convicted
                // guard is untouched).  C++ default false; key omitted when
                // null => byte-identical pre-round-3 config.
                [if open_convicted_bundles != null then 'open_convicted_bundles']: open_convicted_bundles,
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
                         descent_guard=false, guard_descent_cos_y=null,
                         guard_descent_min_cm=null,
                         vertex_hadron_guard=false, guard_hadron_len_cm=null,
                         guard_hadron_mip=null,
                         entry_rise_guard=false, guard_entry_frac=null,
                         guard_entry_min_cm=null, guard_entry_max_cm=null,
                         guard_entry_min_len_cm=null, guard_entry_kink_deg=null,
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
              // doc-94 round 1 descent veto: a cosmic stopping muon must have
              // travelled DOWNWARD to reach its stop.  C++ default false and
              // cos_y default +1.01 (above the feature range = pure probe);
              // keys omitted when off => byte-identical pre-knob config.
              + (if descent_guard then { descent_guard: true } else {})
              + (if guard_descent_cos_y != null then { guard_descent_cos_y: guard_descent_cos_y } else {})
              + (if guard_descent_min_cm != null then { guard_descent_min_cm: guard_descent_min_cm } else {})
              // doc-94 round 1 vertex-hadron veto: a long, heavily-ionizing prong
              // off the fitted main is a neutrino vertex hadron, not the
              // second muon check_other_tracks looks for.  C++ defaults
              // false/12/1.5; keys omitted when off => byte-identical.
              + (if vertex_hadron_guard then { vertex_hadron_guard: true } else {})
              + (if guard_hadron_len_cm != null then { guard_hadron_len_cm: guard_hadron_len_cm } else {})
              + (if guard_hadron_mip != null then { guard_hadron_mip: guard_hadron_mip } else {})
              // doc-94 round 2 entry-rise veto: a contiguous elevated dQ/dx run
              // ANCHORED at the boundary end of the fit and decaying to the
              // body level is two particles, one of which left the detector.
              // The only predicate in TaggerCheckSTM that reads the end where
              // the fit STARTS.  Round 3 adds the owner's KINK: the path must
              // also TURN somewhere, which is what tells a two-particle vertex
              // from a delta-ray fluctuation on a straight muon.
              // C++ defaults false/1.3/5/60/70/22; keys omitted when off =>
              // byte-identical pre-knob config.
              + (if entry_rise_guard then { entry_rise_guard: true } else {})
              + (if guard_entry_frac != null then { guard_entry_frac: guard_entry_frac } else {})
              + (if guard_entry_min_cm != null then { guard_entry_min_cm: guard_entry_min_cm } else {})
              + (if guard_entry_max_cm != null then { guard_entry_max_cm: guard_entry_max_cm } else {})
              + (if guard_entry_min_len_cm != null then { guard_entry_min_len_cm: guard_entry_min_len_cm } else {})
              + (if guard_entry_kink_deg != null then { guard_entry_kink_deg: guard_entry_kink_deg } else {})
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

        tagger_check_neutrino(name="", trackfitting_config_file="", particle_dataset="", recombination_model="", perf=false, dl_weights="", dQdx_scale=0.1, dQdx_offset=-1000.0, clus_geom_helper="", dl_vtx_rerank=true, dl_vtx_top_k=5, dl_vtx_min_accept_score=4.0, dl_vtx_score_scale=1000.0, beam_window_low=0, beam_window_high=0, nu_skip_cosmic_bundle_min_length=0, dir_weak_use_score=false, mip_dqdx_median=null, proton_dir_vote=false, proton_dir_score_max=null, proton_dir_asym_min=null, endpoint_trim_retry=false, fit_vertex_min_seg_length=null, nu_per_bundle=false, nu_per_bundle_demoted_acts=false, nu_per_bundle_min_length=null, fiducial=null, fv_tolerance=[], teb_min_len=null, teb_min_arm=null, teb_min_arm_pts=null, teb_stub_max=null, teb_accept_range=null, teb_rise_r1=null, teb_rise_r2=null, teb_abs_end_min=null, teb_dip_floor=null, teb_score_cap_r1=null, teb_score_cap_r2=null, teb_turn_angle=null, teb_turn_baseline=null, teb_turn_skirt=null, kink_dqdx_hot_ratio=null, knobs={}) :: {
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
              // C++ default 0 (cm).  Key omitted when 0 => byte-identical pre-knob config.
              // When > 0 (needs nu_skip_cosmic_bundle): an UNTAGGED in-window main at least
              // this long survives the bundle veto -- the taggers examined it and declined to
              // tag it; its cosmic-tagged bundle-mate stays excluded regardless (docs/pr/16
              // design A: keeps vetoing unexamined shards like SBND evt 52195's 1.3 cm main).
              + (if nu_skip_cosmic_bundle_min_length > 0 then { nu_skip_cosmic_bundle_min_length: nu_skip_cosmic_bundle_min_length } else {})
              // nu_per_bundle / nu_per_bundle_demoted_acts (doc
              // sbnd_xin/docs/pr/94 Phase 2): run the PR chain once per
              // in-beam-window flash bundle instead of once per event, so a
              // cosmic-convicted activity can no longer take a co-bundled
              // neutrino candidate down with it (SBND 18255/395148).  The
              // per-main cosmic veto is kept; only the event-level BUNDLE veto
              // is dropped in this mode.  nu_per_bundle_demoted_acts mirrors
              // the cosmic taggers' evaluate_demoted_mains admission gate so
              // the per-activity act_evaluated column is exact -- wire it from
              // the SAME variable that feeds TaggerCheck{TGM,STM,FC}.  C++
              // defaults false; keys omitted when off => byte-identical
              // pre-pr/94 config.
              + (if nu_per_bundle then { nu_per_bundle: true } else {})
              + (if nu_per_bundle && nu_per_bundle_demoted_acts then { nu_per_bundle_demoted_acts: true } else {})
              // doc pr/94 Phase 5b round 2 -- length floor (cm) for a
              // per-bundle candidate, exempting the one activity the legacy
              // event-wide selector would itself have picked.  C++ default 0
              // = no floor.  Key omitted when null or when nu_per_bundle is
              // off => byte-identical compiled config.
              + (if nu_per_bundle && nu_per_bundle_min_length != null then { nu_per_bundle_min_length: nu_per_bundle_min_length } else {})
              // C++ default false.  Key omitted when off => byte-identical pre-knob config (uBooNE).
              // When on: direction-weakness reads use segment_is_dir_weak() (score thresholds),
              // the faithful port of prototype ProtoSegment::is_dir_weak() -- see sbnd_xin/docs/pr/6.
              + (if dir_weak_use_score then { dir_weak_use_score: true } else {})
              // MIP dQ/dx scales in e/cm.  C++ defaults 50000 (flat-template role) /
              // 43000 (median-threshold role) = the uBooNE legacy values; keys omitted
              // when null => byte-identical pre-knob config.  SBND: 56000 / 48000
              // (owner 2026-07-30; docs pr/7 sec 5 + pr/8).
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
              // doc pr/48 -- back-to-back track fixes.  Bools: C++ default
              // false, key omitted when off => byte-identical.  teb_*
              // numerics: null = C++ defaults (10/1.8/4/4/15 cm-ish floors,
              // rise 1.3/1.15, abs 1.7x median, caps 0.6/0.9, turn 25 deg /
              // 35 cm / 3 cm), inert while two_end_break is off.
              + (if teb_min_len != null then { teb_min_len: teb_min_len } else {})
              + (if teb_min_arm != null then { teb_min_arm: teb_min_arm } else {})
              + (if teb_min_arm_pts != null then { teb_min_arm_pts: teb_min_arm_pts } else {})
              + (if teb_stub_max != null then { teb_stub_max: teb_stub_max } else {})
              + (if teb_accept_range != null then { teb_accept_range: teb_accept_range } else {})
              + (if teb_rise_r1 != null then { teb_rise_r1: teb_rise_r1 } else {})
              + (if teb_rise_r2 != null then { teb_rise_r2: teb_rise_r2 } else {})
              + (if teb_abs_end_min != null then { teb_abs_end_min: teb_abs_end_min } else {})
              + (if teb_dip_floor != null then { teb_dip_floor: teb_dip_floor } else {})
              + (if teb_score_cap_r1 != null then { teb_score_cap_r1: teb_score_cap_r1 } else {})
              + (if teb_score_cap_r2 != null then { teb_score_cap_r2: teb_score_cap_r2 } else {})
              + (if teb_turn_angle != null then { teb_turn_angle: teb_turn_angle } else {})
              + (if teb_turn_baseline != null then { teb_turn_baseline: teb_turn_baseline } else {})
              + (if teb_turn_skirt != null then { teb_turn_skirt: teb_turn_skirt } else {})
              // doc pr/90 round 4 -- chain-topology gate admission (D1),
              // route R3 turn/activity thresholds (deg / x mip median, D3),
              // R2 bragg veto turn (deg, D4).  C++ defaults false/0 =
              // legacy; false/null = key omitted => byte-identical pre-fix
              // config.
              + (if kink_dqdx_hot_ratio != null then { kink_dqdx_hot_ratio: kink_dqdx_hot_ratio } else {})
              // The knob bag: every default-OFF pattern-recognition knob whose only
              // job is to reach this component as a key.  The caller emits the key or
              // omits it -- an absent key is the C++ default, exactly as the per-knob
              // `+ (if x then {x: ...})` clauses this replaces did.  SBND builds it in
              // wct-pr-perevt.jsonnet (doc 77 round 2); uBooNE builds it inline.
              + knobs,
        },

        // Run pattern recognition (find_proto_vertex) on the main cluster.
        // mode is passed for future use (e.g. "multiple" for multi-track mode).
        do_tracking(name="", mode="", perf=false, clus_geom_helper="") :: $.tagger_check_neutrino(name=name, perf=perf, clus_geom_helper=clus_geom_helper),

        // Run numu CC BDT scoring (TMVA/xgboost).
        // Must run AFTER tagger_check_neutrino in the visitor list.
        // XML weight files should be resolved from wire-cell-data uboone/weights/.
        // Pass empty strings to disable (scorer will skip booking and EvaluateMVA).
        // fast_xgb_forest: C++ default false = TMVA::Reader books the XGB
        // combiner (legacy).  true = TmvaGradForest, the same score from a
        // compact scan of the same XML (sbnd_xin/docs/76 round 2).  Key
        // omitted when off => byte-identical pre-knob config.
        numu_bdt_scorer(name="",
                        numu1_weights_xml="",
                        numu2_weights_xml="",
                        numu3_weights_xml="",
                        cosmict10_weights_xml="",
                        numu_xgboost_xml="",
                        fast_xgb_forest=false) :: {
            type: "UbooneNumuBDTScorer",
            name: prefix + name,
            data: {
                grouping: "live",
                numu1_weights_xml:    numu1_weights_xml,
                numu2_weights_xml:    numu2_weights_xml,
                numu3_weights_xml:    numu3_weights_xml,
                cosmict10_weights_xml: cosmict10_weights_xml,
                numu_xgboost_xml:     numu_xgboost_xml,
                [if fast_xgb_forest then 'fast_xgb_forest']: true,
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
                       nue_xgboost_xml="",
                       // see numu_bdt_scorer: C++ default false; key omitted when off
                       fast_xgb_forest=false) :: {
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
                [if fast_xgb_forest then 'fast_xgb_forest']: true,
            }
        },

        // Write T_tagger and T_kine trees into the existing tracking output ROOT file.
        // Must run AFTER numu_bdt_scorer and nue_bdt_scorer (BDT scores must be filled).
        // Must run AFTER UbooneMagnifyTrackingVisitor (file must already exist to UPDATE).
        tagger_output(name="", output_filename="tracking_proj.root", neutrino_type_bitmask=false, nu_per_bundle=false, mcs_output=false) :: {
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
              + (if neutrino_type_bitmask then { neutrino_type_bitmask: true } else {})
              // doc sbnd_xin/docs/pr/94 Phase 1: book the per-bundle identity
              // + per-activity cosmic-flag branches (cluster_id,
              // matched_flash_gid, nu_index, act_*).  Plumbing only -- nothing
              // populates them yet.  C++ default false = branches not booked;
              // key omitted when off => byte-identical pre-knob config AND
              // schema.
              + (if nu_per_bundle then { nu_per_bundle: true } else {})
              // doc 80 round 3: book the five kine_mcs_* T_kine branches (MCS
              // muon momentum).  C++ default false = branches not booked; key
              // omitted when off => byte-identical pre-knob config AND schema.
              + (if mcs_output then { mcs_output: true } else {}),
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
                              adopt_min_npts=null, adopt_beam_min_length=null,
                              rescue_allow_in_beam_far=false,
                              rescue_geom_first=false, geom_first_dis=null,
                              rescue_pierce_test=false, pierce_cut=null,
                              conn_drift_frac=null, conn_min_dis=null,
                              rescue_dest_beam_for_new=false,
                              far_contain_tol=null,
                              rescue_beam_main_only=false) :: {
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
                // ---- round 2 (sbnd_xin/docs/73) -------------------------------
                // Four independent openings of a measured blocker, each C++
                // default false / inert.  Keys omitted when off => compiled
                // config byte-identical to pre-round-2.
                //
                // rescue_allow_in_beam_far: K_far may itself be in the beam
                // window provided it is in a different flash bundle (both halves
                // matched, each to its own side's in-beam flash; the T0
                // hypothesis then moves by <= 0.31 cm, sub-blob).
                [if rescue_allow_in_beam_far then "rescue_allow_in_beam_far"]: true,
                // rescue_geom_first: test a pair the [-rescue_t0_early,
                // +rescue_t0_late] window rejected, but demand a TIGHTENED
                // geometry (geom_first_dis + pierce_cut + collinearity) -- the
                // wrong flash is measured up to 855 us away, beyond any time
                // prior.  Purely additive to the legacy accept path.
                [if rescue_geom_first then "rescue_geom_first"]: true,
                // null => C++ default 8 cm (a genuine one-sided crosser's
                // tip-to-tip separation runs to 6.7 cm, evt65053).
                [if geom_first_dis != null then "geom_first_dis"]: geom_first_dis,
                // rescue_pierce_test: where the tip-to-tip vector is dominated
                // by drift (it is then the cathode dead gap, and the conn angle
                // silently becomes a cut on |dir_x| > cos(conn_far_cut)) or too
                // short to define a direction, substitute the cathode-PIERCING
                // agreement -- a fixed transverse bound, not one that scales
                // with the tip separation.
                [if rescue_pierce_test then "rescue_pierce_test"]: true,
                // nulls => C++ defaults (8 cm / 0.8 / 8 cm).
                [if pierce_cut != null then "pierce_cut"]: pierce_cut,
                [if conn_drift_frac != null then "conn_drift_frac"]: conn_drift_frac,
                [if conn_min_dis != null then "conn_min_dis"]: conn_min_dis,
                // rescue_dest_beam_for_new: a pair admitted ONLY by one of the
                // three knobs above adopts the beam bundle, instead of the
                // length-based a/b/c/d rule which can send it to the cosmic one
                // when the beam-side donor is still a pre-collapse stub.  Legacy
                // pairs keep a/b/c/d unchanged.
                [if rescue_dest_beam_for_new then "rescue_dest_beam_for_new"]: true,
                // far_contain_tol: how far a round-2 pair's FAR half may end up
                // on the wrong side of the cathode once re-materialized under the
                // destination T0.  The tip tests say nothing about the rest of
                // the cluster, and a cosmic matched hundreds of us away lands
                // tens of cm inside the other TPC.  null => C++ default 1 cm
                // (a direct bound on the T0-hypothesis error: 1 cm = 6.4 us;
                // good merges overshoot <= 0.55 cm, false ones by 33.5).
                // Inert unless a round-2 knob is on.
                [if far_contain_tol != null then "far_contain_tol"]: far_contain_tol,
                // Round 3 (sbnd_xin/docs/73 sec 12): the beam-side donor must
                // BE its bundle's matched main (Flags::main_cluster).  On SBND
                // data evt 51128 a 3.8 cm associated fragment donated the beam
                // T0 to a 283.9 cm cosmic and the merge displaced the bundle's
                // real 57.7 cm main.  C++ default false.  Key omitted when
                // off => byte-identical pre-round-3 config.
                [if rescue_beam_main_only then "rescue_beam_main_only"]: true,
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
        // graph_name (C++ default "ctpc"; key omitted when null => byte-identical
        // pre-knob config): doc 79 round 2 -- "ctpc_fast" arms the busy-cluster
        // lazy walk of the >30 cm skeleton shortest-path graph build.
        deghost(name="", use_ctpc=true, length_cut=0, allow_mixed_faces=null, empty_view_unique=null, graph_name=null) :: {
            type: "ClusteringDeghost",
            name: prefix+name,
            data: {
                use_ctpc: use_ctpc,
                length_cut: length_cut,
                [if allow_mixed_faces != null then 'allow_mixed_faces']: allow_mixed_faces,
                [if empty_view_unique != null then 'empty_view_unique']: empty_view_unique,
                [if graph_name != null then 'graph_name']: graph_name,
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

        // busy_num_threshold (C++ default 0 = off; key omitted when null =>
        // byte-identical pre-knob config): doc 78 round 3 -- busy-cluster
        // lazy walk in the overclustering-protection pair loop, the same
        // mode as the "relaxed_fast" graph flavor.
        protect_overclustering(name="", busy_num_threshold=null) :: {
            type: "ClusteringProtectOverclustering",
            name: prefix+name,
            data: dv_cfg + pcts_cfg + scope_cfg
              + (if busy_num_threshold != null then { busy_num_threshold: busy_num_threshold } else {}),
            uses: [detector_volumes, pc_transforms],
        },

        // protect_iso_band (default false == byte-identical): decline to merge
        // an isochronous band (narrow drift slab, large y-z footprint) with a
        // non-band cluster unless the two genuinely touch — the extended-cloud
        // prolongations otherwise bridge tens of cm.
        neutrino(name="", num_try=1, use_flash_t0=false, flash_t0_window=80*wc.ns, protect_iso_band=false, protect_iso_band_xext=null, record_band_veto=false) :: {
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
                // Record each iso-band refusal as a per-blob "nu_band_veto_role"
                // provenance array so the SEPARATE all-APA clustering chain --
                // which has no iso-band guard of its own -- declines to
                // re-merge the pair (doc pr/66).  C++ default false; key
                // omitted when off => byte-identical compiled config.  Only
                // meaningful with protect_iso_band on.
                [if record_band_veto then 'record_band_veto']: record_band_veto,
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
