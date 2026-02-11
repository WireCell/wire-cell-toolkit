/// Process depo files to dnnroi-truth (ground truth ROI masks).
///
/// This is a single-TPC job producing only the splat aka "truth" tier.
///
/// Minimal example usage:
///
/// wire-cell spng/cfg/spng/depos-to-splat.jsonnet -A input=depos.npz 

local wc = import "wirecell.jsonnet";
local pg = import 'pgraph.jsonnet';
local io = import "spng/io.jsonnet";
local tio = import "spng/torchio.jsonnet";
local drift_js = import "spng/drift.jsonnet";
local subgraphs_js = import "spng/subgraphs.jsonnet";
local control_js = import "spng/control.jsonnet";
local frame_js = import "spng/frame.jsonnet";
local hacks = import "spng/hacks.jsonnet";
local detconf = import "spng/detconf.jsonnet";

// The TLAs:
//
// @param input The name of a file in WCT "depo file" format, usually .npz.
// @param output The output file name.
// @param detname The name of a supported detector, default "pdhd".
// @param engine The name of the graph execution engine, default Pgrapher or TbbFlow.
// @param tpcid The TPC ID number.
//
// The only required TLA is "input".
//
// Notes:
//
// - Input depos should be arranged to populate the given tpcid.
function(input,
         output="splat.npz",
         detname='pdhd',
         engine='Pgrapher',
         tpcid=3)


    local controls = control_js();
    local control = controls.config;

    // We focus here on just one TPC
    local det = detconf.get(detname, [wc.numberify(tpcid)]);
    local tpc = det.tpcs[0];

    local frame = frame_js(control);

    // Source of depos
    local source = io.depo_source(input);

    // Common drifter
    local drift = drift_js(det, control).drifter;

    local upstream = pg.pipeline([source, drift]);


    local sg = subgraphs_js(tpc, control);

    local truth = sg.splat(extra_name="_TRUTH");

    local sink = io.frame_sink(output);
    
    local graph = pg.pipeline([upstream, truth, sink]);

    local result = pg.main(graph, app=engine,
                           plugins=["WireCellSpng", "WireCellGen", "WireCellSio"],
                           uses=controls.uses);
    result


