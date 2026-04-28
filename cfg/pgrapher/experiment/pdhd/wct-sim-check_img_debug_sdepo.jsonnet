# usage: wire-cell -l stdout wct-sim-check-with-img-clustering.jsonnet
# Complete simulation + signal processing + imaging + clustering pipeline

local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';

local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local params = import 'pgrapher/experiment/pdhd/simparams.jsonnet';

local tools = tools_maker(params);

local sim_maker = import 'pgrapher/experiment/pdhd/sim.jsonnet';
local sim = sim_maker(params, tools);

// local track= {
// head: wc.point(0.0359088, 467.694, 35.7493, wc.cm),
// tail: wc.point(17.8166, 328.196, 446.613, wc.cm),
// };
local track= {
head: wc.point(10, 467.694, 35.7493, wc.cm),
tail: wc.point(10, 328.196, 446.613, wc.cm),
};

local tracklist = [
  {
    time: 0 * wc.us,
    charge: -500, // negative means # electrons per step (see below configuration) 
    // ray: params.det.bounds,
    ray: track,

  },
];

local depos = sim.tracks(tracklist, step=0.1 * wc.mm); // MIP <=> 5000e/mm

local nanodes = std.length(tools.anodes);
local anode_iota = std.range(0, nanodes-1);
local anode_idents = [anode.data.ident for anode in tools.anodes];

local drifter = sim.drifter;
local bagger = sim.make_bagger();

local depo_saver = g.pnode({
  type: "NumpyDepoSaver", 
  name: "depo_saver",
  data: {
    filename: "depositions.npz"
  }
}, nin=1, nout=1);

// signal plus noise pipelines
local sn_pipes = sim.splusn_pipelines;

local perfect = import 'pgrapher/experiment/pdhd/chndb-base.jsonnet';
local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  data: perfect(params, tools.anodes[n], tools.field, n){dft:wc.tn(tools.dft)},
  uses: [tools.anodes[n], tools.field, tools.dft],
} for n in anode_iota];

local nf_maker = import 'pgrapher/experiment/pdhd/nf.jsonnet';
local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n) for n in std.range(0, std.length(tools.anodes) - 1)];

local sp_override = {
    sparse: false,
    use_roi_debug_mode: true,
    use_multi_plane_protection: true,
    process_planes: [0, 1, 2]
};

local sp_maker = import 'pgrapher/experiment/pdhd/sp.jsonnet';
local sp = sp_maker(params, tools, sp_override);
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];

// IMAGING PIPELINE
local img = import 'pgrapher/experiment/pdhd/img_simple.jsonnet';
local img_maker = img();
local img_pipes = [img_maker.per_anode(a,'single') for a in tools.anodes];

// CLUSTERING PIPELINE (Optional - can be enabled/disabled)
// local clus = import 'pgrapher/experiment/pdhd/clus.jsonnet';
// local clus_maker = clus();
// local clus_pipes = [clus_maker.per_apa(tools.anodes[n], dump=true) for n in std.range(0, std.length(tools.anodes) - 1)];

local magoutput = 'protodunehd-sim-check-full.root';
local magnify = import 'pgrapher/experiment/pdhd/magnify-sinks.jsonnet';
local magnifyio = magnify(tools, magoutput);



local parallel_pipes = [
  g.pipeline([
               sn_pipes[n],
               // magnifyio.orig_pipe[n],
               nf_pipes[n],
              //  magnifyio.raw_pipe[n],
               sp_pipes[n],
               // magnifyio.debug_pipe[n],
              //  magnifyio.decon_pipe[n],
               img_pipes[n],
             ],
             'parallel_pipe_%d' % n)
  for n in std.range(0, std.length(tools.anodes) - 1)
];

local outtags = ['raw%d' % n for n in std.range(0, std.length(tools.anodes) - 1)];
// local parallel_graph = f.fanpipe('DepoSetFanout', parallel_pipes, 'FrameFanin', 'sn_mag_nf_img', outtags);
local nanodes = std.length(tools.anodes);
local parallel_graph = f.multifanout('DepoSetFanout', parallel_pipes, [1,nanodes], [nanodes,1], 'sn_mag_nf_img', outtags);

// // Add clustering sink for final output
// local clustering_sink = if use_clustering then
//   g.pnode({
//     type: "TensorFileSink",
//     name: "clustering_sink",
//     data: {
//       outname: "protodunehd-sim-check-clusters.tar.gz",
//       prefix: "clustering_",
//       dump_mode: true,
//     }
//   }, nin=1, nout=0)
// else null;

local sink = sim.frame_sink;

// Build final graph
// local graph = g.pipeline([depos, drifter, bagger, parallel_graph, sink]);
// local graph = g.pipeline([depos, drifter, bagger, depo_saver,parallel_graph]);
local graph = g.pipeline([depos, drifter,  depo_saver,bagger, parallel_graph]);

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
};

// local cmdline = {
//     type: "wire-cell",
//     data: {
//         plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio", "WireCellSigProc", "WireCellRoot"],
//         apps: ["Pgrapher"]
//     }
// };

local cmdline = {
    type: "wire-cell",
    data: {
        // plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio", "WireCellSigProc", "WireCellRoot", "WireCellTbb", "WireCellImg", "WireCellPytorch"],
        plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio", "WireCellSigProc", "WireCellRoot", "WireCellTbb", "WireCellImg"],
        apps: ["Pgrapher"]
    }
};

// Finally, the configuration sequence which is emitted.

[cmdline] + g.uses(graph) + [app]