
local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';
local util = import 'pgrapher/experiment/pdhd/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';

local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';

local params = import 'pgrapher/experiment/pdhd/simparams.jsonnet';

local tools = tools_maker(params);


local sim_maker = import 'pgrapher/experiment/pdhd/sim.jsonnet';
local sim = sim_maker(params, tools);

local wcls_maker = import "pgrapher/ui/wcls/nodes.jsonnet";
local wcls = wcls_maker(params, tools);
local wcls_input = {
    depos: wcls.input.depos(name="", art_tag="IonAndScint"),
};

local bagger = sim.make_bagger("bagger%d");

// local output_file=std.extVar('outname');

local sinks = g.pnode({
    name: "deposink",
    type: "DepoFileSink",
    data: {
        // outname: output_file,
    }
}, nin=1, nout=0);

local save = g.pnode({
    type: "NumpyDepoSaver",
    name: 'numpysave',
    data: {
        // filename: fname
    }
}, nin=1, nout=1);
local dump = g.pnode({type: 'DumpDepos', name: 'dump_depos', data: {}}, nin=1, nout=0);

local graph = g.pipeline([wcls_input.depos, save, dump], 'mygraph');
// g4 sim as input
// local graph = g.intern(
//     innodes=[wcls_input.depos], centernodes=[bagger], outnodes=[sinks],
//     edges = 
//         [
//             g.edge(wcls_input.depos, bagger),
//             g.edge(bagger, sinks, 0, 0),

//         ],
// );

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
};
g.uses(graph) + [app]
