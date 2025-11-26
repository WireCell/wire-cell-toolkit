// This produces a final graph to input depos and output crossview tensors

local jobs = import "spng/jobs.jsonnet";
local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local io = import "spng/io.jsonnet";
local tio = import "spng/torchio.jsonnet";
local tpc_nodes = import "spng/tpc.jsonnet";



local detconf = import "spng/detconf.jsonnet";
local detector = import "spng/detector.jsonnet";
local control_module = import "spng/control.jsonnet";


/// Top-level arguments:
///
/// input :: file name pattern where %d receives anode ID.
/// output :: if given, a file name with '%d' to receive anode ID number
/// view_crossed :: which views should be subject to crossviews
function(input, output, view_crossed=[1,1,0],
         detname='pdhd', tpcids=[], engine='Pgrapher', device='cpu', verbosity=0)
    
    local control = control_module.bundle(device=device, verbosity=wc.intify(verbosity));
    local det = detector.subset(detconf[detname], tpcids);

    // Iframe -> cv tensor set per TPC
    local tpc_pipes = [
        pg.pipeline([
            // fixme: refactor to include detsim
            io.frame_array_source(input%tpc.ident),
            tdm(tpc, control).graph_crossview,
            // fixme: refactor to make a single output file
            tio.pickle_tensor_set(output%tpc.ident),
        ])
        for tpc in det.tpcs];

    local graph = pg.components(tpc_pipes);
    pg.main(graph, engine, plugins=["WireCellSpng"])



//// original content
// local cv = import "crossviews.jsonnet";
// local decon = import "decon.jsonnet";
// local detconf = import "detconf.jsonnet";
// local det = detconf.pdhd;

// local tests = {
//     local tpc = det.tpcs[0],
//     tpc: tpc,
//     tv: cv.threshold_views(tpc),

//     cv0: cv.crossview(tpc, 0),

//     cf: cv.crossfan(tpc, $.tv),
// };
// function(which="")
//     if which ==""
//     then tests
//     else tests[which]

