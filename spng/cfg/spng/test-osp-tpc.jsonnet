local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local fio = import "fileio.jsonnet";

local io = import "spng/io.jsonnet";
local tio = import "spng/torchio.jsonnet";

local detconf = import "spng/detconf.jsonnet";
local detector = import "spng/detector.jsonnet";
local control_mod = import "spng/control.jsonnet";

function(input, output='test-osp-tpc-%(tier)s.npz',
         tpcid=0, 
         detname='pdhd',
         engine='Pgrapher', device='cpu')
    // FIXME: device is ignored

    local controls = control_mod(device=device);

    local det = detector.subset(detconf[detname], [wc.intify(tpcid)]);
    local tpc = det.tpcs[tpcid];

    local source = io.frame_array_source(input);

    local sp = tpc.osp_subgraph.sp;
    local sp_tags = [
        "loose_lf" + std.toString(tpcid),
        "mp2_roi" + std.toString(tpcid),
        "mp3_roi" + std.toString(tpcid),
        "decon_charge" + std.toString(tpcid)];
        
    local sptap = fio.tap('FrameFanout', io.frame_array_sink(output % "osp", tags=sp_tags));
    local roi = tpc.osp_subgraph.dnnroi;

    local sink = io.frame_array_sink(output % "sig");

    local graph = pg.pipeline([source, sp, sptap, roi, sink]);
    pg.main(graph, engine,
            plugins=["WireCellSpng", "WireCellSigProc", "WireCellPytorch"],
            uses = controls.uses)


