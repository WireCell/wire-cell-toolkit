local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";
local f = import 'pgrapher/common/funcs.jsonnet';
// local params = import "pgrapher/experiment/dune-vd/params.jsonnet";
// local tools_maker = import 'pgrapher/common/tools.jsonnet';
// local tools = tools_maker(params({}));
// local anodes = tools.anodes;

local clus (
    anode,
    face = 0,
    drift_speed = 1.6 * wc.mm / wc.us,
    time_offset = 0 * wc.us,
    bee_dir = "data",
    bee_zip = "mabc.zip",
    initial_index = "0",
    initial_runNo = "1",
    initial_subRunNo = "1",
    initial_eventNo = "1") =
{
    local index = std.parseInt(initial_index),

    local LrunNo = std.parseInt(initial_runNo),
    local LsubRunNo = std.parseInt(initial_subRunNo),
    local LeventNo  = std.parseInt(initial_eventNo),

    // Note, the "sampler" must be unique to the "sampling".
    local bs_live = {
        type: "BlobSampler",
        name: "%s-%d"%[anode.name, face],
        data: {
            drift_speed: drift_speed,
            time_offset: time_offset,
            strategy: [
                // "center",
                // "corner",
                // "edge",
                // "bounds",
                "stepped",
                // {name:"grid", step:1, planes:[0,1]},
                // {name:"grid", step:1, planes:[1,2]},
                // {name:"grid", step:1, planes:[2,0]},
                // {name:"grid", step:2, planes:[0,1]},
                // {name:"grid", step:2, planes:[1,2]},
                // {name:"grid", step:2, planes:[2,0]},
            ],
            // extra: [".*"] // want all the extra
            extra: [".*wire_index"] //
            // extra: [] //
        }},
    local bs_dead = {
        type: "BlobSampler",
        name: "%s-%d"%[anode.name, face],
        data: {
            strategy: [
                "center",
            ],
            extra: [".*"] // want all the extra
        }},


    local detector_volumes = 
    {
        "type": "DetectorVolumes",
        "name": "",
        "data": {
            "anodes": [wc.tn(anode)],
        }
    },
    
    local geom_helper = {
        type: "SimpleClusGeomHelper",
        name: "uboone",
        data: {
            a0f0: {
                face: 0,
                pitch_u: 3 * wc.mm,
                pitch_v: 3 * wc.mm,
                pitch_w: 3 * wc.mm,
                angle_u: 1.0472,    // 60 degrees
                angle_v: -1.0472,   // -60 degrees
                angle_w: 0,         // 0 degrees
                drift_speed: drift_speed*1,
                tick: 0.5 * wc.us,  // 0.5 mm per tick
                tick_drift: self.drift_speed * self.tick,
                time_offset: time_offset,
                nticks_live_slice: 4,
                FV_xmin: 1 * wc.cm,
                FV_xmax: 255 * wc.cm,
                FV_ymin: -99.5 * wc.cm,
                FV_ymax: 101.5 * wc.cm,
                FV_zmin: 15 * wc.cm,
                FV_zmax: 1022 * wc.cm,
                FV_xmin_margin: 2 * wc.cm,
                FV_xmax_margin: 2 * wc.cm,
                FV_ymin_margin: 2.5 * wc.cm,
                FV_ymax_margin: 2.5 * wc.cm,
                FV_zmin_margin: 3 * wc.cm,
                FV_zmax_margin: 3 * wc.cm
            },
            a1f0: self.a0f0 {},
            a2f0: self.a0f0 {},
            a3f0: self.a0f0 {},
            a4f0: self.a0f0 {},
            a5f0: self.a0f0 {},
            a6f0: self.a0f0 {},
            a7f0: self.a0f0 {},
            a8f0: self.a0f0 {},
            a9f0: self.a0f0 {},
            a10f0: self.a0f0 {},
            a11f0: self.a0f0 {},
            a12f0: self.a0f0 {},
            a13f0: self.a0f0 {},
            a14f0: self.a0f0 {},
            a15f0: self.a0f0 {},
            a16f0: self.a0f0 {},
            a17f0: self.a0f0 {},
            a18f0: self.a0f0 {},
            a19f0: self.a0f0 {},
            a20f0: self.a0f0 {},
            a21f0: self.a0f0 {},
            a22f0: self.a0f0 {},
            a23f0: self.a0f0 {},
        }
    },

    local ptb = g.pnode({
        type: "PointTreeBuilding",
        name: "%s-%d"%[anode.name, face],
        data:  {
            samplers: {
                "3d": wc.tn(bs_live),
                "dead": wc.tn(bs_dead),
            },
            multiplicity: 2,
            tags: ["live", "dead"],
            anode: wc.tn(anode),
            face: face,
            geom_helper: wc.tn(geom_helper),
            detector_volumes: "DetectorVolumes",
        }
    }, nin=2, nout=1, uses=[bs_live, bs_dead, detector_volumes]),

    local mabc = g.pnode({
        type: "MultiAlgBlobClustering",
        name: "%s-%d"%[anode.name, face],
        data:  {
            inpath: "pointtrees/%d",
            outpath: "pointtrees/%d",
            // grouping2file_prefix: "grouping%s-%d"%[anode.name, face],
            perf: true,
            bee_dir: bee_dir, // "data/0/0", // not used
            bee_zip: "mabc-%s-face%d.zip"%[anode.name, face],
            bee_detector: "sbnd",
            initial_index: index,   // New RSE configuration
            use_config_rse: true,  // Enable use of configured RSE
            runNo: LrunNo,
            subRunNo: LsubRunNo,
            eventNo: LeventNo,
            save_deadarea: true, 
            anode: wc.tn(anode),
            detector_volumes: "DetectorVolumes",
            face: face,
            geom_helper: wc.tn(geom_helper),
            func_cfgs: [
                {name: "clustering_test", detector_volumes: "DetectorVolumes"},
                {name: "clustering_live_dead", dead_live_overlap_offset: 2},
                {name: "clustering_extend", flag: 4, length_cut: 60 * wc.cm, num_try: 0, length_2_cut: 15 * wc.cm, num_dead_try: 1},
                {name: "clustering_regular", length_cut: 60*wc.cm, flag_enable_extend: false},
                {name: "clustering_regular", length_cut: 30*wc.cm, flag_enable_extend: true},
                {name: "clustering_parallel_prolong", length_cut: 35*wc.cm},
                {name: "clustering_close", length_cut: 1.2*wc.cm},
                {name: "clustering_extend_loop", num_try: 3},
                {name: "clustering_separate", use_ctpc: true},
                {name: "clustering_connect1"},
                {name: "clustering_deghost"},
                {name: "clustering_examine_x_boundary"},
                {name: "clustering_protect_overclustering"},
                {name: "clustering_neutrino"},
                {name: "clustering_isolated"},
            ],
        }
    }, nin=1, nout=1, uses=[geom_helper]),

    local sink = g.pnode({
        type: "TensorFileSink",
        name: "%s-%d"%[anode.name, face],
        data: {
            outname: "clus-%s-face%d.tar.gz"%[anode.name, face],
            prefix: "clustering_", // json, numpy, dummy
            dump_mode: true,
        }
    }, nin=1, nout=0),

    clus_pipe(dump=true) ::
    if dump then
        g.pipeline([ptb, mabc, sink], "clus_pipe-%s-%d"%[anode.name, face])
    else
        g.pipeline([ptb, mabc], "clus_pipe-%s-%d"%[anode.name, face]),
};

function () {
    per_volume(anode, face=0, dump=true) :: clus(anode, face=face).clus_pipe(dump),
}