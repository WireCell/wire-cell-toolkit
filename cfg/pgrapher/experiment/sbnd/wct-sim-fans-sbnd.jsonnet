local g = import "pgraph.jsonnet";
local f = import "pgrapher/experiment/sbnd/funcs.jsonnet";
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
local art_label = std.extVar("art_label");
local wiener_label = std.extVar("input_wiener_label");
local gauss_label = std.extVar("input_gauss_label");
local badmasks_label = std.extVar("badmasks_inputTag");
local threshold_label = std.extVar("threshold_tag");

// must match name used in fcl
local wcls_input = g.pnode({
    type: 'wclsCookedFrameSource', //added wcls Ewerton 2023-07-27
    name: 'sigs',
    data: {
        art_tag: art_label,                   // art inputTag for recob::Wire (no need to setup it up if imaging=true below)
        nticks: params.daq.nticks,
        wiener_inputTag: wiener_label,        // input recob::Wire (wiener)
        gauss_inputTag: gauss_label,          // input recob::Wire (gauss)
        //badmasks_inputTag: badmasks_label,    // input bad masks
        badmasks_inputTag: ["wctsp:badmasks"], //["wctsp:bad0masks","wctsp:bad1masks"],
        threshold_inputTag: threshold_label,  // input threshold
        frame_tags: ["orig"],                 // frame tags (only one frame in this module)
        cmm_tag: "bad",                        // single label for now
        scale: 50,                             // scale up input recob::Wire by this factor
        debug_channel: 6800                   // debug purposes. Deleteme later. 
    },
}, nin=0, nout=1);

/* uncomment this only if used in this jsonnet
local anode = tools.anodes[0]; //added Ewerton 2024-01-22 
local mag_test = g.pnode({
          type: 'MagnifySink',
          //name: 'magimgtest0', //commented Ewerton 2024-01-22
          name: 'magimgtest', //added Ewerton 2024-01-22
          data: {
                output_filename: "magoutput-cooked.root",
                root_file_mode: 'UPDATE',
                frames: ['gauss', 'wiener'],
                summaries: ['wiener'],
                summary_operator: { ['wiener']: 'set' },
                trace_has_tag: true,
                //anode: wc.tn(tools.anodes), //commented Ewerton 2024-01-22
                anode: wc.tn(anode), //added Ewerton 2024-01-22
          },
}, nin=1, nout=1);
*/

// Parallel part //////////////////////////////////////////////////////////////////////////////

local img = import 'pgrapher/experiment/sbnd/img.jsonnet'; //added Ewerton 2023-06-09
local img_maker = img();
local img_pipes = [img_maker.per_anode(a) for a in tools.anodes];

// added Ewerton 2023-09-08 
  local parallel_pipes = [
    g.pipeline([
                  img_pipes[n],
            ],
            'parallel_pipe_%d' % n)
  for n in std.range(0, std.length(tools.anodes) - 1)];


// added Ewerton 2023-09-10
local parallel_graph = f.fanpipe2('FrameFanout', parallel_pipes, 'FrameFanin', 'fanout'); //added Ewerton 2023-09-14


// Final pipeline //////////////////////////////////////////////////////////////////////////////

local graph = g.pipeline([wcls_input, parallel_graph], "main"); // added Ewerton 2023-09-08
//local graph = g.pipeline([wcls_input, mag_test, parallel_graph], "main"); // added Ewerton 2023-09-08

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
