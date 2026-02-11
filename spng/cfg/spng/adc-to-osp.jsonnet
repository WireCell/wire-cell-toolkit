// This job converts ADC to OSP signals
// FIXME: needs NF for real data.

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

local io = import "spng/io.jsonnet";
local detconf = import "spng/detconf.jsonnet";
local detector = import "spng/detector.jsonnet";
local control_js = import "spng/control.jsonnet";

local sg_js = import "spng/subgraphs.jsonnet";

// The TLAs:
//
// @param input The name of a file in WCT "depo file" format, usually .npz.
// @param output The output file name.
// @param detname The name of a supported detector, default "pdhd".
// @param tpcid The TPC ID number.
// @param engine The name of the graph execution engine, default Pgrapher or TbbFlow.
// @param device The name of the device for SPNG nodes, default "cpu" or "gpu", "gpu1", etc. 
// @param verbosity The verbosity level for additional logging.
//
// The only required TLA is "input".
//
// Notes:
//
// - Input depos should be arranged to populate the given tpcid.

function(input,
         output="osp.npz",
         detname='pdhd',
         tpcid=3,
         engine='Pgrapher',
         device='cpu')
    
    local controls = control_js(device=device);
    local det = detconf.get(detname, [wc.numberify(tpcid)]);
    local tpc = det.tpcs[0];

    local source = io.frame_array_source(input);

    local sp = tpc.osp_subgraphs.sp;
    local sp_tags = [
        "loose_lf" + std.toString(tpcid),
        "mp2_roi" + std.toString(tpcid),
        "mp3_roi" + std.toString(tpcid),
        "decon_charge" + std.toString(tpcid)];
        
    local roi = tpc.osp_subgraphs.dnnroi;

    local sink = io.frame_array_sink(output);

    local graph = pg.pipeline([source, sp, roi, sink]);
    pg.main(graph, engine,
            plugins=["WireCellSpng", "WireCellSigProc", "WireCellPytorch"],
            uses = controls.uses)


