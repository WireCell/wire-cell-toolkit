local g = import "pgraph.jsonnet";
local f = import "pgrapher/common/funcs.jsonnet";
local wc = import "wirecell.jsonnet";

local tools_maker = import 'pgrapher/common/tools.jsonnet';

// added Ewerton 2023-09-06 
local reality = 'sim'; // <- fixme put this in fcl!!
local params_maker =
if reality == 'data' then import 'params.jsonnet'
else import 'simparams.jsonnet';

local base = import 'pgrapher/experiment/sbnd/simparams.jsonnet';
local params = base {
  lar: super.lar { // <- super.lar overrides default values
    // Longitudinal diffusion constant
    DL: std.extVar('DL') * wc.cm2 / wc.s,
    // Transverse diffusion constant
    DT: std.extVar('DT') * wc.cm2 / wc.s,
    // Electron lifetime
    lifetime: std.extVar('lifetime') * wc.ms,
    // Electron drift speed, assumes a certain applied E-field
    drift_speed: std.extVar('driftSpeed') * wc.mm / wc.us,
  },
};


local tools_all = tools_maker(params);
local tools = tools_all {anodes: [tools_all.anodes[n] for n in [0,1]]}; //added Ewerton 2023-09-08


// added Ewerton 2023-09-08
local wiener_label = std.extVar("input_wiener_label");
local gauss_label = std.extVar("input_gauss_label");
local badmasks_label = std.extVar("badmasks_inputTag");
local threshold_label = std.extVar("threshold_tag");

// must match name used in fcl
local wcls_input = g.pnode({
    type: 'wclsCookedFrameSource', //added wcls Ewerton 2023-07-27
    name: 'sigs',
    data: {
        nticks: params.daq.nticks,
        wiener_inputTag: wiener_label,        // input recob::Wire (wiener)
        gauss_inputTag: gauss_label,          // input recob::Wire (gauss)
        badmasks_inputTag: [badmasks_label],    // input bad masks
        threshold_inputTag: threshold_label,  // input threshold
        frame_tags: ["orig"],                 // frame tags (only one frame in this module)
        cmm_tag: "bad",                        // single label for now
        scale: 50,                             // scale up input recob::Wire by this factor
        debug_channel: 6800                   // debug purposes. Deleteme later. 
    },
}, nin=0, nout=1);

local img = import 'pgrapher/experiment/sbnd/img.jsonnet';
local img_maker = img();
local img_pipes = [img_maker.per_anode(a, "multi", add_dump = true) for a in tools.anodes];

local fanout_apa_rules =
[
    {
        frame: {
            //'.*': 'number%d' % n,
            //'.*': 'gauss%d' % n,
            //'.*': 'framefanout%d ' % n,
            '.*': 'orig%d' % n,
        },
        trace: {
            // fake doing Nmult SP pipelines
            //orig: ['wiener', 'gauss'],
            gauss: 'gauss%d' % n, //uncommented Ewerton 2023-09-27
            wiener: 'wiener%d' % n, //created Ewerton 2023-09-27
            //'.*': 'orig',
        },
    }
    for n in std.range(0, std.length(tools.anodes) - 1)
];
local parallel_graph = f.fanout("FrameFanout", img_pipes, "parallel_graph", fanout_apa_rules);

local graph = g.pipeline([wcls_input, parallel_graph], "main"); // added Ewerton 2023-09-08

local app = {
  type: 'Pgrapher', //Pgrapher, TbbFlow
  data: {
    edges: g.edges(graph),
  },
};

local cmdline = {
    type: "wire-cell",
    data: {
        plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio", "WireCellSigProc", "WireCellRoot", "WireCellTbb", "WireCellImg"],
        apps: ["Pgrapher"] //TbbFlow
    }
};

[cmdline] + g.uses(graph) + [app]
