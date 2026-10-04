// dune10kt-1x2x6/clus.jsonnet -- clustering module for DUNE FD-HD 1x2x6 (12 APAs, both faces live).
//
// Promoted VERBATIM from wcp-porting-validation wcfm/clus.jsonnet (wcfm/docs/02 sec 3, docs/13,
// docs/25), itself forked BY DUPLICATION from pdhd/clus.jsonnet (untouched; CLAUDE.md M10).  The
// wcfm file is now a re-export of this one.  What differs from PDHD, and why:
//   * geometry from clus_params.jsonnet: drift speed 1.6 mm/us, the `dvm` DetectorVolumes
//     metadata generated for ALL 12 anodes x 2 faces (PDHD lists anodes 0-3 by hand), the
//     overall FV box of the 1x2x6 module, bee_detector;
//   * topology: PDHD merges the two faces of an APA (stage 2) and then groups x-aligned APAs
//     by ident parity.  Here BOTH faces of every APA are live and are OPPOSITE drift volumes
//     (validate_drift_group refuses a mixed-face scope, and ident parity is the y row), so
//     there is no per-APA stage: per-(anode, face) stage 1 (with the PDHD stage-2
//     deghost/protect_overclustering appended: one (anode, face) is one imaging view, as one
//     PDHD APA is) -> per-FACE group over all run anodes (stage 3) -> all-TPC
//     (stage 4) with switch_scope only.  cathode_connect is dropped: the cathodes are the
//     outer walls of this module, nothing crosses them;
//   * no Q/L, no flash keys (there is no light here).
// The visitor lists and every cut value are the PDHD production ones.
// The live sampler also takes the PR job's retile-sampler options (pr.jsonnet), with defaults
// that compile to the clustering job's sampler unchanged.
local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";
local f = import 'pgrapher/common/funcs.jsonnet';
local clus = import "pgrapher/common/clus.jsonnet";
local P = import "pgrapher/experiment/dune10kt-1x2x6/clus_params.jsonnet";

function (output_dir='', runNo=1, subRunNo=1, eventNo=1,
          time_offset=0 * wc.us, trigger_offset=0 * wc.us)

local drift_speed = P.drift_speed;
local initial_index = "0";
local index = std.parseInt(initial_index);

local common_coords = ["x", "y", "z"];
local common_corr_coords = ["x_t0cor", "y", "z"];

local aname(anode) = "anode%d" % anode.data.ident;

// Per-face drift-volume block: face 0 = +x volume [W plane, cathode], face 1 = mirror.
local face_block(face) = {
    drift_speed: drift_speed,
    tick: P.tick,
    tick_drift: self.drift_speed * self.tick,
    time_offset: time_offset,
    trigger_offset: trigger_offset,
    nticks_live_slice: P.tick_span,
    FV_xmin: if face == 0 then P.fv_x.wire else -P.fv_x.cathode,
    FV_xmax: if face == 0 then P.fv_x.cathode else -P.fv_x.wire,
    FV_xmin_margin: 2 * wc.cm,
    FV_xmax_margin: 2 * wc.cm,
};

local dvm = {
    overall: P.fv_overall {
        FV_xmin_margin: 2 * wc.cm,
        FV_xmax_margin: 2 * wc.cm,
        FV_ymin_margin: 2.5 * wc.cm,
        FV_ymax_margin: 2.5 * wc.cm,
        FV_zmin_margin: 3 * wc.cm,
        FV_zmax_margin: 3 * wc.cm,
        vertical_dir: [0, 1, 0],
        beam_dir: [0, 0, 1],
    },
} + {
    ["a%df%dpA" % [a, face]]: face_block(face)
    for a in std.range(0, 11)
    for face in [0, 1]
};

// Component name suffix of a DetectorVolumes / PCTransformSet: the anode idents joined by "-", then the face as
// "-f<face>" when the volume is face-specific.  The face used to be appended as a bare "-<face>", so the face-1 volume
// of anode 0 ("dv-apa0-1") and the anode-pair volume of anodes 0 and 1 ("dv-apa0-1") were ONE component configured
// twice; the pair volume then carried anode 0's metadata only and T0Correction::forward threw std::out_of_range
// (map::at on apa 1) for every event whose active anodes are exactly {0, 1}.  The "f" makes the two names distinct;
// component names do not reach any output (wcfm doc 33: pctree archives byte-identical on 8 events).
local anodes_name(anodes, face="") =
    std.join("-", [std.toString(a.data.ident) for a in anodes]) + if face == "" then "" else "-f" + std.toString(face);

// Both faces of every anode are always registered (Grouping::fill_dv_cache reads every face).
local detector_volumes(anodes, face="") = {
    "type": "DetectorVolumes",
    "name": "dv-apa" + anodes_name(anodes, face),
    "data": {
        "anodes": [wc.tn(anode) for anode in anodes],
        metadata:
            {overall: dvm["overall"]} +
            { ["a" + std.toString(a.data.ident) + "f0pA"]: dvm["a" + std.toString(a.data.ident) + "f0pA"] for a in anodes } +
            { ["a" + std.toString(a.data.ident) + "f1pA"]: dvm["a" + std.toString(a.data.ident) + "f1pA"] for a in anodes }
    },
    uses: anodes
};

local pctransforms(dv) = {
    type: "PCTransformSet",
    name: dv.name,
    data: { detector_volumes: wc.tn(dv) },
    uses: [dv]
};

// strategy_name 'charge_stepped' (+ wire_product, charge_threshold) is for the PR job's retile
// samplers only (pr.jsonnet improve2), as pdhd/clus.jsonnet and sbnd/clus.jsonnet bs_live_face: the
// port of the prototype's retile sampler (CalcPoints.cxx:75-160).  It carries its own name ('live-cs-')
// so it can never replace the clustering sampler of the same (anode, face) in one config.
// Default 'stepped' => the clustering sampler, compiled config unchanged.
local bs_live_face(apa, face, wrapped_channel_charge=true,
                   strategy_name='stepped', wire_product=null, charge_threshold=null) = {
    assert strategy_name == 'stepped' || strategy_name == 'charge_stepped' :
        "bs_live_face: strategy_name must be 'stepped' or 'charge_stepped'",
    type: "BlobSampler",
    name: if strategy_name == 'charge_stepped' then "live-cs-%s-%d"%[apa, face] else "live-%s-%d"%[apa, face],
    data: {
        drift_speed: drift_speed,
        time_offset: time_offset,
        strategy: if strategy_name == 'charge_stepped' then
            [{name: "charge_stepped",
              disable_mix_dead_cell: false,
              [if wire_product != null then 'max_wire_product_threshold']: wire_product,
              [if charge_threshold != null then 'charge_threshold_max']: charge_threshold,
              [if charge_threshold != null then 'charge_threshold_min']: charge_threshold,
              [if charge_threshold != null then 'charge_threshold_other']: charge_threshold}]
        else ["stepped"],
        extra: [".*wire_index", ".*charge_val", ".*charge_unc", "wpid"],
        // FD APAs wrap their induction planes exactly like PDHD's (1149/1148 wires on 800
        // channels); PDHD production ON since 2026-09-06 (pdhd/clus.jsonnet).
        wrapped_channel_charge: wrapped_channel_charge,
    }
};
local bs_dead_face(apa, face) = {
    type: "BlobSampler",
    name: "dead-%s-%d"%[apa, face],
    data: {
        strategy: ["center"],
        extra: [".*"] // want all the extra
    }
};

// Stage 1: one (anode, face), ending with the PDHD stage-2 steps.  wcfm/docs/13 sec 1: the
// stage-2 deghost (empty_view_unique false, the C++ default, as PDHD runs it) seeds its 2-D
// view with the longest cluster only; in a multi-anode scope every other anode's view is empty,
// an empty view counts as an overlap and every cluster of every other anode is removed.  Until
// doc 13 it ran at the head of the multi-anode stage 3 and kept one anode only.
local clus_per_face (
    anode,
    face,
    dump = true,
    bee_dir = "data",
    runNo = 1,
    subRunNo = 1,
    eventNo = 1,
    wrapped_channel_charge = true,
    ) =
{
    local an = aname(anode),
    local dv = detector_volumes([anode], face),
    local pcts = pctransforms(dv),

    local cluster_scope_filter_live = g.pnode({
        type: "ClusterScopeFilter",
        name: "csf-live-%s-%d"%[an, face],
        data: { face_index: face },
    }, nin=1, nout=1, uses=[]),

    local cluster_scope_filter_dead = g.pnode({
        type: "ClusterScopeFilter",
        name: "csf-dead-%s-%d"%[an, face],
        data: { face_index: face },
    }, nin=1, nout=1, uses=[]),

    local bsl = bs_live_face(an, face, wrapped_channel_charge=wrapped_channel_charge),
    local bsd = bs_dead_face(an, face),

    local ptb = g.pnode({
        type: "PointTreeBuilding",
        name: "%s-%d"%[an, face],
        data:  {
            samplers: { "3d": wc.tn(bsl), "dead": wc.tn(bsd) },
            multiplicity: 2,
            tags: ["live", "dead"],
            anode: wc.tn(anode),
            face: face,
            detector_volumes: wc.tn(dv),
        }
    }, nin=2, nout=1, uses=[bsl, bsd, dv]),

    local cluster2pct = g.intern(
        innodes = [cluster_scope_filter_live, cluster_scope_filter_dead],
        centernodes = [],
        outnodes = [ptb],
        edges = [
            g.edge(cluster_scope_filter_live, ptb, 0, 0),
            g.edge(cluster_scope_filter_dead, ptb, 0, 1)
        ]
    ),

    local face_name = "%s-%d"%[an, face],

    local cm = clus.clustering_methods(prefix=face_name,
                                       detector_volumes=dv,
                                       pc_transforms=pcts,
                                       coords=common_coords),
    local cm_pipeline = [
        cm.pointed(),
        cm.live_dead(dead_live_overlap_offset=2),
        cm.extend(flag=4, length_cut=60*wc.cm, num_try=0, length_2_cut=15*wc.cm, num_dead_try=1),
        cm.regular(name="-one", length_cut=60*wc.cm, flag_enable_extend=false),
        cm.regular(name="_two", length_cut=30*wc.cm, flag_enable_extend=true),
        cm.parallel_prolong(length_cut=35*wc.cm),
        cm.close(length_cut=1.2*wc.cm),
        cm.extend_loop(num_try=3),
        cm.connect1(),
        // PDHD stage 2 (per APA, one live face there == one (anode, face) here)
        cm.deghost(),
        cm.protect_overclustering(),
    ],

    local mabc = g.pnode({
        local name = "%s-%d"%[an, face],
        type: "MultiAlgBlobClustering",
        name: name,
        data:  {
            inpath: "pointtrees/%d",
            outpath: "pointtrees/%d",
            perf: true,
            bee_dir: bee_dir,
            bee_zip: "%s/mabc-%s-face%d.zip"%[bee_dir, an, face],
            bee_detector: P.bee_detector,
            initial_index: index,
            use_config_rse: true,
            runNo: runNo,
            subRunNo: subRunNo,
            eventNo: eventNo,
            save_deadarea: true,
            dead_area_version: 2,
            anodes: [wc.tn(anode)],
            face: face,
            detector_volumes: wc.tn(dv),
            bee_points_sets: [
                {
                    name: "clustering",
                    detector: P.bee_detector,
                    algorithm: "clustering",
                    pcname: "3d",
                    coords: ["x", "y", "z"],
                    individual: true
                }
            ],
            pipeline: wc.tns(cm_pipeline),
        }
    }, nin=1, nout=1, uses=[dv, anode, pcts]+cm_pipeline),

    local sink = g.pnode({
        type: "TensorFileSink",
        name: "clus_per_face-%s-%d"%[an, face],
        data: {
            outname: "%s/trash-%s-face%d.tar.gz"%[bee_dir, an, face],
            prefix: "clustering_",
            dump_mode: true,
        }
    }, nin=1, nout=0),

    local end = if dump then g.pipeline([mabc, sink]) else g.pipeline([mabc]),

    ret :: g.pipeline([cluster2pct, end], "clus_per_face-%s-%d"%[an, face])
}.ret;

// Stage 3: one drift volume = all run anodes seen through ONE face: the PDHD stage-3 list (the
// PDHD stage-2 deghost + protect_overclustering now end stage 1, per (anode, face); doc 13 sec 1).
local clus_per_group (
    anodes,
    group_name,
    face,
    dump = true,
    bee_dir = "data",
    runNo = 1,
    subRunNo = 1,
    eventNo = 1,
    save_assoc_id = false,
    ) = {
    local nanodes = std.length(anodes),
    local pcmerging = g.pnode({
        type: "PointTreeMerging",
        name: "clus_per_group-%s"%group_name,
        data:  {
            multiplicity: nanodes,
            inpath: "pointtrees/%d",
            outpath: "pointtrees/%d",
            tolerate_missing: true,
        }
    }, nin=nanodes, nout=1),

    local dv = detector_volumes(anodes, face),
    local pcts = pctransforms(dv),

    local cm = clus.clustering_methods(prefix=group_name,
                                       detector_volumes=dv,
                                       pc_transforms=pcts,
                                       coords=common_coords),
    local cm_pipeline = [
        // PDHD stage 3
        cm.extend(flag=4, length_cut=60*wc.cm, num_try=0, length_2_cut=15*wc.cm, num_dead_try=1),
        cm.regular(name="1", length_cut=60*wc.cm, flag_enable_extend=false),
        cm.regular(name="2", length_cut=30*wc.cm, flag_enable_extend=true),
        cm.parallel_prolong(length_cut=35*wc.cm),
        cm.close(length_cut=1.2*wc.cm),
        cm.extend_loop(num_try=3),
        cm.separate(use_ctpc=true, max_hull_points=1000000, collinear_recover=true, collinear_interior=true,
                    collinear_member_merge=true,
                    track_repartition=true, band_merge_back=true, band_recarve=true, drift_side_fv_x=true,
                    far_point_x_cut=14*wc.cm, far_point_mid_dis=60*wc.cm, track_recarve=true, dec1_guard_main_angle=45,
                    iso_slab_split=true, tag_family=true, collinear_global_merge=true),
        cm.connect1(respect_separate_family=true),
        cm.deghost(empty_view_unique=true),
        cm.examine_x_boundary(),
        cm.neutrino(protect_iso_band=true),
        cm.isolated(save_assoc_id=save_assoc_id),
    ],

    local mabc = g.pnode({
        type: "MultiAlgBlobClustering",
        name: "clus_per_group-%s"%group_name,
        data:  {
            inpath: "pointtrees/%d",
            outpath: "pointtrees/%d",
            perf: true,
            bee_dir: bee_dir,
            bee_zip: "%s/mabc-%s.zip"%[bee_dir, group_name],
            bee_detector: P.bee_detector,
            initial_index: index,
            use_config_rse: true,
            runNo: runNo,
            subRunNo: subRunNo,
            eventNo: eventNo,
            save_deadarea: true,
            dead_area_version: 2,
            [if save_assoc_id then 'save_assoc_cluster_id']: true,
            anodes: [wc.tn(a) for a in anodes],
            detector_volumes: wc.tn(dv),
            bee_points_sets: [
                {
                    name: "clustering",
                    detector: P.bee_detector,
                    algorithm: "clustering",
                    pcname: "3d",
                    coords: ["x", "y", "z"],
                    individual: false,
                }
            ],
            pipeline: wc.tns(cm_pipeline),
        }
    }, nin=1, nout=1, uses=anodes+[dv, pcts]+cm_pipeline),

    local sink = g.pnode({
        type: "TensorFileSink",
        name: "clus_per_group-%s"%group_name,
        data: {
            outname: "%s/trash-%s.tar.gz"%[bee_dir, group_name],
            prefix: "clustering_",
            dump_mode: true,
        }
    }, nin=1, nout=0),

    local end = if dump then g.pipeline([mabc, sink]) else g.pipeline([mabc]),

    ret :: g.intern(
        innodes = [pcmerging],
        centernodes = [],
        outnodes = [end],
        edges = [ g.edge(pcmerging, end, 0, 0) ]
    ),
}.ret;

// Stage 4: all-TPC, merging the two face groups; switch_scope only (no cathode crossers, no light).
local clus_all_tpc (
    anodes,
    ngroups = 2,
    dump = true,
    bee_dir = "data",
    runNo = 1,
    subRunNo = 1,
    eventNo = 1,
    tensor_outname = '',
    save_assoc_id = false,
    ) = {
    local pcmerging = g.pnode({
        type: "PointTreeMerging",
        name: "clus_all_tpc",
        data:  {
            multiplicity: ngroups,
            inpath: "pointtrees/%d",
            outpath: "pointtrees/%d",
            tolerate_missing: true,
        }
    }, nin=ngroups, nout=1),

    local dv = detector_volumes(anodes),
    local pcts = pctransforms(dv),

    local cm_old = clus.clustering_methods(prefix="all",
                                           detector_volumes=dv,
                                           pc_transforms=pcts,
                                           coords=common_coords),
    local cm_pipeline = [
        cm_old.switch_scope(),
    ],

    // Dead-area Bee groups by drift side (the PDHD apa_drift_groups role).
    local dead_groups = [
        { name: "groupf0", apas: [a.data.ident for a in anodes], face: 0 },
        { name: "groupf1", apas: [a.data.ident for a in anodes], face: 1 },
    ],

    local mabc = g.pnode({
        type: "MultiAlgBlobClustering",
        name: "clus_all_tpc",
        data:  {
            inpath: "pointtrees/%d",
            outpath: "pointtrees/%d",
            perf: true,
            bee_dir: bee_dir,
            bee_zip: "%s/mabc-all-apa.zip"%[bee_dir],
            bee_detector: P.bee_detector,
            initial_index: index,
            use_config_rse: true,
            runNo: runNo,
            subRunNo: subRunNo,
            eventNo: eventNo,
            save_deadarea: true,
            [if save_assoc_id then 'save_assoc_cluster_id']: true,
            dead_area_version: 2,
            dead_apa_groups: dead_groups,
            anodes: [wc.tn(a) for a in anodes],
            detector_volumes: wc.tn(dv),
            bee_points_sets: [
                {
                    name: "clustering",
                    detector: P.bee_detector,
                    algorithm: "clustering",
                    pcname: "3d",
                    coords: common_corr_coords,
                    individual: false
                },
                {
                    name: "img",
                    detector: P.bee_detector,
                    algorithm: "img",
                    pcname: "3d",
                    coords: ["x", "y", "z"],
                    individual: false,
                }
            ],
            pipeline: wc.tns(cm_pipeline),
        },
    }, nin=1, nout=1, uses=anodes+[dv, pcts]+cm_pipeline),

    local sink = g.pnode({
        type: "TensorFileSink",
        name: "clus_all_tpc",
        data: {
            outname: if tensor_outname == '' then "%s/trash-all-apa.tar.gz"%[bee_dir] else tensor_outname,
            prefix: "clustering_",
            dump_mode: tensor_outname == '',
        }
    }, nin=1, nout=0),
    local end = if dump then g.pipeline([mabc, sink]) else g.pipeline([mabc]),
    ret :: g.intern(
        innodes = [pcmerging],
        centernodes = [],
        outnodes = [end],
        edges = [ g.edge(pcmerging, end, 0, 0) ]
    ),
}.ret;


{
    local bee_dir = if output_dir == '' then 'data' else output_dir,
    per_face(anode, face=0, dump=true, wrapped_channel_charge=true) :: clus_per_face(anode, face=face, dump=dump, bee_dir=bee_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo, wrapped_channel_charge=wrapped_channel_charge),
    per_group(anodes, group_name, face, dump=true, save_assoc_id=false) :: clus_per_group(anodes, group_name, face, dump=dump, bee_dir=bee_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo, save_assoc_id=save_assoc_id),
    all_tpc(anodes, ngroups=2, dump=true, tensor_outname='', save_assoc_id=false) :: clus_all_tpc(anodes, ngroups=ngroups, dump=dump, bee_dir=bee_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo, tensor_outname=tensor_outname, save_assoc_id=save_assoc_id),
    detector_volumes(anodes, face="") :: detector_volumes(anodes, face),
    pc_transforms(dv) :: pctransforms(dv),
    live_sampler(anode, face, wrapped_channel_charge=true, strategy_name='stepped', wire_product=null, charge_threshold=null) ::
        bs_live_face(aname(anode), face, wrapped_channel_charge=wrapped_channel_charge,
                     strategy_name=strategy_name, wire_product=wire_product, charge_threshold=charge_threshold),
    drift_speed :: drift_speed,
    time_offset :: time_offset,
    scope_coords :: common_coords,
    t0cor_coords :: common_corr_coords,
}
