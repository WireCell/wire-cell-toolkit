// Build everything.  depos to splat, osp and spng 

local pg = import 'pgraph.jsonnet';
local det_mod = import "spng/det.jsonnet";
local io = import "spng/io.jsonnet";

local control_mod = import "spng/control.jsonnet";

local detconf = import "spng/detconf.jsonnet";

function(input, output, detname='pdhd', tpcids=[], engine='Pgrapher', device='cpu', verbosity=0)

    local controls = control_mod(device=device, verbosity=verbosity);


    local det = detconf.get(detname, tpcids);

    // make source node
    local source = io.depo_source(input);

    // Sink a node that has one oport per TPC across the det.
    local sink_det(src, tier) = pg.shuntline(src,  pg.crossline(
        [io.frame_array_sink(output % {tpcid: tpc.ident, tier: tier})
         for tpc in det.tpcs]));


    // returns .{depos_sink, {splat,osp,spng}_source}.
    local guts = det_mod(det, controls.config).kitchen_sink;

    local head = pg.pipeline([source, guts.depo_sink]);
    local splat = sink_det(guts.splat_source, "splat");
    local osp = sink_det(guts.osp_source, "osp");
    local spng = sink_det(guts.spng_source, "spng");


    local graph = pg.components([head, splat, osp, spng]);

    pg.main(graph, app=engine,
            plugins=["WireCellSpng", "WireCellSigProc", "WireCellGen", "WireCellPytorch"],
            uses=controls.uses)

