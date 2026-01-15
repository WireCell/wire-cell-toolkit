/// This provides a (1 deposet) -> 3 x (per-view ROI tensor) for one TPC
/// (drifted depos)->splat->reframer->totdm->resample->threshold->(rois x 3 views)

local wc = import "wirecell.jsonnet";
local real_pg = import "pgraph.jsonnet";

local sim_js = import "sim.jsonnet";
local frame_js = import "frame.jsonnet";
local tpc_js = import "tpc.jsonnet";
local decon_js = import "decon.jsonnet";

function(tpc, control, downsample_factor=4, splat_threshold=100.0, pg=real_pg) 

    local sim = sim_js(tpc, control);
    local frame = frame_js(control);
    local tpcmod = tpc_js(tpc, control, pg);
    local decon = decon_js(tpc, control);

    // [1]frame -> tensor[ngroups]
    local to_tdm = frame.to_tdm(tpc, extra_name="_splat");
    local to_tens = frame.tensorset_unpacker(tpc, "_splat");

    local frame_source = pg.pipeline([
        sim.splat,
        sim.reframer("_splat"),
        
        to_tdm,                 // frame->tensor SET
        to_tens                 // tensor SET to per group tensors
    ]);

    // [ngroups]tensor->tensor[nviews]
    local frame_port_nodes = [pg.oport_node(frame_source, oport.port) for oport in frame_source.oports];
    local view_splats = decon.collect_groups_by_view(frame_port_nodes);

    local upstream = pg.intern(innodes=[frame_source], outnodes=view_splats);


    local one_view(view_index) = 
        local this_name = tpc.name + "_splat_v" + std.toString(view_index);
        local downsampler = pg.pnode({
            type:'SPNGResampler',
            name: this_name,
            data: {
                ratio: 1.0/downsample_factor,
            } + control
        }, nin=1, nout=1);
        local thresholder = pg.pnode({
            type: 'SPNGThreshold',
            name: this_name,
            data: {
                nominal: splat_threshold,
                tag: "roi",
                datapath_format: "/traces/Threshold/" + this_name,
            } + control
        }, nin=1, nout=1);
        pg.pipeline([downsampler, thresholder]);
    local downstream = pg.crossline([one_view(view_index) for view_index in [0,1,2]]);


    // After this we have crossline and not atomic nodes
    local subgraph = pg.shuntlines([
        upstream,
        downstream,
    ]);
    subgraph
