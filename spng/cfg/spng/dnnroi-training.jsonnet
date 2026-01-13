/// Process depo files to dnnroi-input and dnnroi-truth.
///
/// This is a single-TPC jobs

local wc = import "wirecell.jsonnet";
local pg = import 'pgraph.jsonnet';
local io = import "spng/io.jsonnet";
local torchio = import "spng/torchio.jsonnet";
local det_js = import "spng/det.jsonnet";
local drift_js = import "spng/drift.jsonnet";
local splatroi_js = import "spng/splatroi.jsonnet";
local roifodder_js = import "spng/roifodder.jsonnet";
local control_js = import "spng/control.jsonnet";

local detconf = import "spng/detconf.jsonnet";

// The TLAs:
//
// @param input The name of a file in WCT "depo file" format, usually .npz.
// @param outpat The output file name pattern with format variables.
// @param detname The name of a supported detector, default "pdhd".
// @param engine The name of the graph execution engine, default Pgrapher or TbbFlow.
// @param device The name of the device for SPNG nodes, default "cpu" or "gpu", "gpu1", etc. 
//
// The only required TLA is "input".  
//
// The outpat must include these format variables:
// - %(tier)s will be filled with the label "input" or "truth".
//
// Note, this hard-wires use of TPC ID 0.  Input depos should be arranged to
// populate this TPC.
//
function(input,
         outpat="dnnroi-training-%(tier)s-%(view)s.pkl",
         detname='pdhd',
         engine='Pgrapher',
         device='cpu',
         downsample_factor=4,
         tpcid=0,
         verbosity=0)


    local controls = control_js(device=device, verbosity=wc.intify(verbosity));

    // We focus here on just one TPC
    local det = detconf.get(detname, [tpcid]);
    local tpc = det.tpcs[0];

    // Source of depos
    local source = io.depo_source(input);

    // Common drifter
    local drift = drift_js(det, controls.config).drifter;

    local upstream = pg.pipeline([source, drift]);

    // Fan out depos to the splat and sim+sigproc subgraphs
    local depo_fan = pg.pnode({
        type:'DepoSetFanout',
        name: det.name + "_depo",
        data: {
            multiplicity:2
        }
    },nin=1, nout=2);


    local truth = splatroi_js(tpc, controls.config, downsample_factor);
    local detmod = det_js(det, controls.config);

    local sim = detmod.inducer;

    local fodder = roifodder_js(tpc, controls.config);
    local simfodder = pg.pipeline([sim, fodder]);

    local body = pg.intern(innodes=[upstream], outnodes=[truth, simfodder],
                           centernodes=[depo_fan],
                           edges=[
                               pg.edge(upstream, depo_fan),
                               pg.edge(depo_fan, truth, 0, 0),
                               pg.edge(depo_fan, simfodder, 1, 0)]);


    local truth_sinks = [
        torchio.pickle_tensor_set(outpat % {tier:"truth", view:view})
        for view in [0,1,2]];
    local fodder_sinks = [
        torchio.pickle_tensor_set(outpat % {tier:"fodder", view:view})
        for view in [0,1]];
    local graph = pg.intern(centernodes=[body]+truth_sinks+fodder_sinks,
                            edges=[
                                pg.edge(truth, truth_sinks[view], view, 0),
                                for view in [0,1,2]
                            ] + [
                                pg.edge(fodder, fodder_sinks[view], view, 0),
                                for view in [0,1]]);

    // local graph = pg.components([body, sinks]); 
    //local graph = fodder;
    //local graph = pg.components([sim, fodder]);

    local result = pg.main(graph, app=engine,
                           plugins=["WireCellSpng", "WireCellGen"],
                           uses=controls.uses);
    result

