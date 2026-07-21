local wc = import "wirecell.jsonnet";
local pg = import 'pgraph.jsonnet';
local det_mod = import "spng/det.jsonnet";
local io = import "spng/io.jsonnet";

local control_mod = import "spng/control.jsonnet";

local detconf = import "spng/detconf.jsonnet";

// Test det.jsonnet by having a TLA name the function to use to provide the
// non-IO part of the graph.
//
// The "input_type" names the I/O source to use.  It should be one of:
// - depo for WCT "depo file"
// - frame_tensor for WCT "tensor file"
// - frame_array for WCT frame file"
//
// The "job" names an entry in det.jsonnet.
//
function(input, output="test-det-%(tier)s-tpc%(tpcid)d.npz",
         input_type="depo",
         job="depos_to_adc",
         detname='pdhd', tpcs="", engine='Pgrapher', device='cpu', verbosity=0)

    local tpcids = wc.intlistify(tpcs);
    local controls = control_mod(device=device, verbosity=wc.intify(verbosity));
    local det = detconf.get(detname, tpcids);

    // make source node
    local source = io[input_type + "_source"](input);

    // Sink a node that has one oport per TPC across the det.
    local sink_det(src, tier) = pg.shuntline(src,  pg.crossline(
        [io.frame_array_sink(output % {tpcid: tpc.ident, tier: tier})
         for tpc in det.tpcs]));


    local guts = det_mod(det, controls.config)[job];

    // fixme, while most in det.jsonnet require a single sink, some fanout to
    // many outputs.  For example, kitchen_sink returns .{depos_sink,
    // {splat,osp,spng}_source}.  Need to provide custom output subgraph for
    // them.  For now, any multi-out will break us.
    // local head = pg.pipeline([source, guts.depo_sink]);
    // local splat = sink_det(guts.splat_source, "splat");
    // local osp = sink_det(guts.osp_source, "osp");
    // local spng = sink_det(guts.spng_source, "spng");

    local sink = sink_det(guts, job);

    local graph = pg.shuntlines([source, sink]);

    pg.main(graph, app=engine,
            plugins=["WireCellSpng", "WireCellSigProc", "WireCellGen", "WireCellPytorch"],
            uses=controls.uses)
