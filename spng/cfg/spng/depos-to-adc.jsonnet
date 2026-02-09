/// Process depo files to ADC.
///
/// This is a single-TPC job.
///
/// Minimal example usage:
///
/// wire-cell spng/cfg/spng/depos-to-adc.jsonnet -A input=depos.npz

local wc = import "wirecell.jsonnet";
local pg = import 'pgraph.jsonnet';
local io = import "spng/io.jsonnet";
local tio = import "spng/torchio.jsonnet";
local det_js = import "spng/det.jsonnet";
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
         output="adc.npz",
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

    local detmod = det_js(det, control);

    local sim = detmod.inducer;

    local sink = io.frame_sink(output);


    local graph = pg.pipeline([upstream, sim, sink]);

    // fixme: strictly, only need HIO if saving to HDF.
    local result = pg.main(graph, app=engine,
                           plugins=["WireCellSpng", "WireCellGen", "WireCellHio"],
                           uses=controls.uses);
    result


