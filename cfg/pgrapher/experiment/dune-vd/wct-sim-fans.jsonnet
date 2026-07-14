local g = import "pgraph.jsonnet";
local f = import "pgrapher/common/funcs.jsonnet";
local wc = import "wirecell.jsonnet";

local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local response_plane = std.extVar('response_plane')*wc.cm;

// local fcl_params = {
//     response_plane: 18.92*wc.cm,
//     nticks: 8500,
//     ncrm: 320,
//     wires: 'dunevd10kt_3view_30deg_v6_refactored.json.bz2',
//     use_dnnroi: false,
//     process_crm: 'all',
// };
local fcl_params = {
    response_plane: 18.1*wc.cm, // synced to protodunevd
    nticks: 8500,
    wires: 'dunevd10kt_3view_30deg_v5_refactored_1x8x6ref.json.bz2',
    ncrm: 24,
    use_dnnroi: false,
    process_crm: 'test1',
};
local params_maker =
if fcl_params.ncrm ==320 then import 'pgrapher/experiment/dune-vd/params-10kt.jsonnet'
else import 'pgrapher/experiment/dune-vd/params.jsonnet';

local params = params_maker(fcl_params) {
  lar: super.lar {
    // Longitudinal diffusion constant
    DL: 6.2e-9 * wc.cm2 / wc.ns,
    // Transverse diffusion constant
    DT: 16.3e-9 * wc.cm2 / wc.ns,
    // Electron lifetime
    lifetime: 10.4e3 * wc.us,
    // Electron drift speed, assumes a certain applied E-field
    drift_speed: 1.60563 * wc.mm / wc.us,
  },
  files: super.files {
      wires: fcl_params.wires,
      fields: [ 'protodunevd_FR_imbalance3p_260501.json.bz2', ], // synced to protodunevd
      noise: 'pdvd-top-noise-spectra-v3.json.bz2',               // synced to protodunevd
  },
};

local tools_all = tools_maker(params);
local tools =
if fcl_params.process_crm == "partial"
then tools_all {anodes: [tools_all.anodes[n] for n in std.range(32, 79)]}
else if fcl_params.process_crm == "test1"
then tools_all {anodes: [tools_all.anodes[n] for n in [0,1,4,5]]}
// then tools_all {anodes: [tools_all.anodes[n] for n in [0]]}
else if fcl_params.process_crm == "test2"
then tools_all {anodes: [tools_all.anodes[n] for n in std.range(0, 7)]}
else tools_all;

local sim_maker = import 'pgrapher/experiment/dune-vd/sim.jsonnet';
local sim = sim_maker(params, tools);

// Deposit and drifter ///////////////////////////////////////////////////////////////////////////////

local tracklist = [
    {
        time: 0*wc.ms,
        charge: -5000,
        ray: {
            tail: wc.point(200, -500, 50, wc.cm),
            head: wc.point(200, -500, 200, wc.cm),
        }
    },
    {
        time: 0*wc.ms,
        charge: -5000,
        ray: {
            tail: wc.point(200, -550, 50, wc.cm),
            head: wc.point(200, -550, 200, wc.cm),
        }
    },
];
// local tracklist = [
//     {
//         time: 0*wc.ms,
//         charge: -5000,
//         ray: {
//             tail: wc.point(100, -550, 50, wc.cm),
//             head: wc.point(300, -550, 200, wc.cm),
//         }
//     },
// ];

local depos = g.pnode({
        type: 'TrackDepos',
        data: {
            step_size: 1.0*wc.mm,
            tracks: tracklist
        },
}, nin=0, nout=1);

local drifter = sim.drifter;
local bagger = sim.make_bagger();

// Parallel part //////////////////////////////////////////////////////////////////////////////


// local sn_pipes = sim.signal_pipelines;
local sn_pipes = sim.splusn_pipelines;

local sp_maker = import 'pgrapher/experiment/dune-vd/sp.jsonnet';
local sp_override = {
    sparse: true,
    use_roi_debug_mode: false,
    use_multi_plane_protection: false,
    mp_tick_resolution: 10,
};
local sp = sp_maker(params, tools, sp_override);
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];

local img = import 'pgrapher/experiment/dune-vd/img.jsonnet';
local img_maker = img();
local img_pipes = [img_maker.per_anode(a) for a in tools.anodes];

local magoutput = 'mag-sim-sp.root';
local magnify = import 'pgrapher/experiment/dune-vd/magnify-sinks.jsonnet';
local sinks = magnify(tools, magoutput);
local frame_tap = function(name, outname, tags, digitize) {
    ret: g.fan.tap('FrameFanout',  g.pnode({
        type: "FrameFileSink",
        name: name,
        data: {
            outname: outname,
            tags: tags,
            digitize: digitize,
        },  
    }, nin=1, nout=0), name),
}.ret;
local frame_sink = function(name, outname, tags, digitize) {
    ret: g.pnode({
        type: "FrameFileSink",
        name: name,
        data: {
            outname: outname,
            tags: tags,
            digitize: digitize,
        },
    }, nin=1, nout=0),
}.ret;

local dnnroi = import 'pgrapher/experiment/dune-vd/dnnroi.jsonnet';
local ts = {
    type: "TorchService",
    name: "dnnroi",
    data: {
        model: "ts-model/unet-l23-cosmic500-e50.ts",
        device: "gpucpu",
        concurrency: 1,
    },
};

local parallel_pipes = [
  g.pipeline([ 
                sn_pipes[n],
                // frame_tap(
                //     name="orig%d"%tools.anodes[n].data.ident,
                //     outname="frame-orig%d.tar.bz2"%tools.anodes[n].data.ident,
                //     tags=["orig%d"%tools.anodes[n].data.ident],
                //     digitize=true
                // ),
                // sinks.orig_pipe[n],
                sp_pipes[n],
                // frame_tap(
                //     name="gauss%d"%tools.anodes[n].data.ident,
                //     outname="frame-gauss%d.tar.bz2"%tools.anodes[n].data.ident,
                //     tags=["gauss%d"%tools.anodes[n].data.ident],
                //     digitize=false
                // ),
                sinks.decon_pipe[n],
                // sinks.debug_pipe[n], // use_roi_debug_mode=true in sp.jsonnet
                // dnnroi(tools.anodes[n], ts, output_scale=1.2),
                // sinks.dnnroi_pipe[n],
                // g.pnode({type: "DumpFrames", name: "dumpframes-%d"%tools.anodes[n].data.ident}, nin = 1, nout=0),
                img_pipes[n],
          ], 
          'parallel_pipe_%d' % n) 
  for n in std.range(0, std.length(tools.anodes) - 1)];

local outtags = [];
local tag_rules = {
    frame: {
        '.*': 'framefanin',
    },
    trace: {['gauss%d' % anode.data.ident]: ['gauss%d' % anode.data.ident] for anode in tools.anodes}
        + {['wiener%d' % anode.data.ident]: ['wiener%d' % anode.data.ident] for anode in tools.anodes}
        + {['threshold%d' % anode.data.ident]: ['threshold%d' % anode.data.ident] for anode in tools.anodes}
        + {['dnnsp%d' % anode.data.ident]: ['dnnsp%d' % anode.data.ident] for anode in tools.anodes},
};

// local parallel_graph = f.multifanpipe('DepoSetFanout', parallel_pipes, 'FrameFanin', [1,4], [4,1], [1,4], [4,1], 'sn_mag', outtags, tag_rules);
// The per-anode pipelines terminate in imaging (img per_anode ends in
// dump sinks), so they cannot fan back into a FrameFanin.  Use
// multifanout: fan out to the terminal imaging pipelines, no fanin.
local parallel_graph =
if fcl_params.process_crm == "test1"
then f.multifanout('DepoSetFanout', parallel_pipes, [1,4], [4,1], 'sn_mag', tag_rules)
else if fcl_params.process_crm == "test2"
then f.multifanout('DepoSetFanout', parallel_pipes, [1,8], [8,1], 'sn_mag', tag_rules)
else f.multifanout('DepoSetFanout', parallel_pipes, [1,2,8,32], [2,4,4,10], 'sn_mag', tag_rules);


// Only one sink ////////////////////////////////////////////////////////////////////////////


local sink = sim.frame_sink;


// Final pipeline //////////////////////////////////////////////////////////////////////////////

local graph = g.pipeline([depos, drifter, bagger, parallel_graph], "main");
// local graph = g.pipeline([depos, drifter, bagger, parallel_graph], "main");
// local graph = g.pipeline([depos, drifter, bagger, parallel_pipes[0]], "main");

local app = {
  type: 'TbbFlow', //Pgrapher, TbbFlow
  data: {
    edges: g.edges(graph),
  },
};

local cmdline = {
    type: "wire-cell",
    data: {
        plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio", "WireCellSigProc", "WireCellRoot", "WireCellTbb", "WireCellImg", "WireCellPytorch"],
        apps: ["TbbFlow"]
    }
};

[cmdline] + g.uses(graph) + [app]
