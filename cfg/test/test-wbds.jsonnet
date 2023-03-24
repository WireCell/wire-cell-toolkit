local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local params = import "pgrapher/experiment/uboone/simparams.jsonnet";
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local tools = tools_maker(params);
local sim_maker = import "pgrapher/experiment/uboone/sim.jsonnet";
local sim = sim_maker(params, tools);
local drifter = sim.deposet_filter(sim.drifter, "drifter");
local graph = pg.pipeline([drifter, sim.signal_sets]);
local app = {
    type: 'Pgrapher',
    name: "",
    data: {
        edges: pg.edges(graph),
    },
};

local tests = [
    std.assertEqual(wc.tn(drifter), "DepoSetDrifter:drifter"),
    std.assertEqual(wc.tn(sim.drifter), "Drifter"),
];

pg.uses(graph) + [app]
    

