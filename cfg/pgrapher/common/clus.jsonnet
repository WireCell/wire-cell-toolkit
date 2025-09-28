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

        clustering_recovering_bundle(name="") :: {
            type: "ClusteringRecoveringBundle",
            name: prefix + name,
            data: {
                grouping: "live",           // Which grouping to process
                array_name: "isolated",     // Array name for pcarray lookup  
                pcarray_name: "perblob",    // PCArray name for blob separation
            },
        },

        tagger_check_stm(name="", trackfitting_config_file="", linterp_function="", recombination_model="") :: {
            type: "TaggerCheckSTM",
            name: prefix + name,
            data: {
                grouping: "live",           // Which grouping to process
                trackfitting_config_file: trackfitting_config_file, 
                linterp_function: linterp_function,
                recombination_model: recombination_model,
            } + dv_cfg + pcts_cfg
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

        live_dead(name="", dead_live_overlap_offset=2) :: {
            type: "ClusteringLiveDead",
            name: prefix+name,
            data: {
                dead_live_overlap_offset: dead_live_overlap_offset,
            } + dv_cfg + scope_cfg,
            uses: [detector_volumes],
        },

        extend(name="", flag=0, length_cut=150*wc.cm, num_try=0, length_2_cut=3*wc.cm, num_dead_try=3) :: {
            type: "ClusteringExtend",
            name: prefix+name,
            data: {
                flag: flag,
                length_cut: length_cut,
                num_try: num_try,
                length_2_cut: length_2_cut,
                num_dead_try: num_dead_try,
            } + dv_cfg + scope_cfg,
            uses: [detector_volumes],
        },
        

        regular(name="",  length_cut=45*wc.cm, flag_enable_extend=true) :: {
            type: "ClusteringRegular",
            name: prefix+name,
            data: {
                length_cut: length_cut,
                flag_enable_extend: flag_enable_extend,
            } + dv_cfg + scope_cfg,
            uses: [detector_volumes],
        },

        parallel_prolong(name="", length_cut=35*wc.cm) :: {
            type: "ClusteringParallelProlong",
            name: prefix+name,
            data: {
                length_cut: length_cut,
            } + dv_cfg + scope_cfg,
            uses: [detector_volumes],
        },

        close(name="", length_cut=1*wc.cm) :: {
            type: "ClusteringClose",
            name: prefix+name,
            data: {
                length_cut: length_cut,
            } + scope_cfg,
        },

        extend_loop(name="", num_try=0) :: {
            type: "ClusteringExtendLoop",
            name: prefix+name,
            data: {
                num_try: num_try,
            } + dv_cfg + scope_cfg,
            uses: [detector_volumes],
        },

        separate(name="", use_ctpc=true) :: {
            type: "ClusteringSeparate",
            name: prefix+name,
            data: {
                use_ctpc: use_ctpc,
            } + dv_cfg + pcts_cfg + scope_cfg,
            uses: [detector_volumes, pc_transforms],
        },

        connect1(name="") :: {
            type: "ClusteringConnect1",
            name: prefix+name,
            data: dv_cfg + scope_cfg,
            uses: [detector_volumes],
        },

        deghost(name="", use_ctpc=true, length_cut=0) :: {
            type: "ClusteringDeghost",
            name: prefix+name,
            data: {
                use_ctpc: use_ctpc,
                length_cut: length_cut,
            } + dv_cfg + pcts_cfg + scope_cfg,
            uses: [detector_volumes, pc_transforms],
        },

        isolated(name="") :: {
            type: "ClusteringIsolated",
            name: prefix+name,
            data: dv_cfg + scope_cfg,
        },

        examine_bundles(name="") :: {
            type: "ClusteringExamineBundles",
            name: prefix+name,
            data: dv_cfg + pcts_cfg + scope_cfg,
            uses: [detector_volumes, pc_transforms],
        },

        examine_x_boundary(name="") :: {
            type: "ClusteringExamineXBoundary",
            name: prefix+name,
            data: dv_cfg + scope_cfg,
            uses: [detector_volumes],
        },

        protect_overclustering(name="") :: {
            type: "ClusteringProtectOverclustering",
            name: prefix+name,
            data: dv_cfg + pcts_cfg + scope_cfg,
            uses: [detector_volumes, pc_transforms],
        },

        neutrino(name="", num_try=1) :: {
            type: "ClusteringNeutrino",
            name: prefix+name,
            data: {
                num_try: num_try,
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
        improve_cluster_2(name="", anodes=[], samplers=[], verbose=false) :: {
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
        steiner(name="", retiler={}, grouping="live", graph="steiner") :: {
            type: "CreateSteinerGraph",
            name: prefix+name,
            data: {
                grouping: grouping,
                graph: graph,
                retiler: wc.tn(retiler),
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
