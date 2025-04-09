local g = import "pgraph.jsonnet";
local f = import "pgrapher/common/funcs.jsonnet";
local wc = import "wirecell.jsonnet";

local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local params = import 'pgrapher/experiment/pdhd/simparams.jsonnet';
local fcl_params = {
    use_dnnroi: false,
};

local tools = tools_maker(params);
// local tools_all = tools_maker(params);
// local tools = tools_all {anodes: [tools_all.anodes[n] for n in [0,1,2,3]]};

local sim_maker = import 'pgrapher/experiment/pdhd/sim.jsonnet';
local sim = sim_maker(params, tools);

// Deposit and drifter ///////////////////////////////////////////////////////////////////////////////

local thetaXZ = 45*wc.deg;

local stubby = {
  tail: wc.point(-100, 100, 100, wc.cm),
  head: wc.point(-100*(1 + std.tan(thetaXZ)), 100, 100*(1+1), wc.cm),
};

local tracklist = [

  {
    time: 0 * wc.us,
    charge: -500, // 5000 e/mm
    ray: stubby, // params.det.bounds,
  },

];

// local depos = sim.tracks(tracklist, step=1.0 * wc.mm);
local depos = sim.tracks(tracklist, step=0.1 * wc.mm);


local drifter = sim.drifter;
local bagger = sim.make_bagger();

// Parallel part //////////////////////////////////////////////////////////////////////////////


// local sn_pipes = sim.signal_pipelines;
local sn_pipes = sim.splusn_pipelines;

local sp_maker = import 'pgrapher/experiment/pdhd/sp.jsonnet';
local sp_override = {
    sparse: false,
    use_roi_debug_mode: false,
    use_multi_plane_protection: true,
    mp_tick_resolution: 4,
};
local sp = sp_maker(params, tools, sp_override);
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];

local img = import 'pgrapher/experiment/pdhd/img.jsonnet';
local img_maker = img();
local img_pipes = [img_maker.per_anode(a) for a in tools.anodes];

local magoutput = 'protodunehd-sim-fans.root';
local magnify = import 'pgrapher/experiment/pdhd/magnify-sinks.jsonnet';
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

// local dnnroi = import 'pgrapher/experiment/pdhd/dnnroi.jsonnet';
// local ts = {
//     type: "TorchService",
//     name: "dnnroi",
//     data: {
//         model: "ts-model/unet-l23-cosmic500-e50.ts",
//         device: "gpucpu",
//         concurrency: 1,
//     },
// };

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
                // sinks.decon_pipe[n],
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

// local parallel_graph = f.multifanout('DepoSetFanout', parallel_pipes, [1,4], [4,1], 'sn_mag', tag_rules);

local nanodes = std.length(tools.anodes);
local parallel_graph = f.multifanout('DepoSetFanout', parallel_pipes, [1,nanodes], [nanodes,1], 'sn_mag', tag_rules);



// Only one sink ////////////////////////////////////////////////////////////////////////////


local sink = sim.frame_sink;


// Final pipeline //////////////////////////////////////////////////////////////////////////////

// local graph = g.pipeline([depos, drifter, bagger, parallel_graph, sink], "main");
local graph = g.pipeline([depos, drifter, bagger, parallel_graph], "main");
// local graph = g.pipeline([depos, drifter, bagger, parallel_pipes[0]], "main");

local app = {
  type: 'Pgrapher', //Pgrapher, TbbFlow
  data: {
    edges: g.edges(graph),
  },
};

local cmdline = {
    type: "wire-cell",
    data: {
        // plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio", "WireCellSigProc", "WireCellRoot", "WireCellTbb", "WireCellImg", "WireCellPytorch"],
        plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio", "WireCellSigProc", "WireCellRoot", "WireCellTbb", "WireCellImg"],
        apps: ["Pgrapher"]
    }
};

[cmdline] + g.uses(graph) + [app]
