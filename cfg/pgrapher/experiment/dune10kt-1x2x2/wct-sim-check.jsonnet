// This is a main entry point for configuring a wire-cell CLI job to
// simulate protoDUNE-SP.  It is simplest signal-only simulation with
// one set of nominal field response function.  It excludes noise.
// The kinematics are a mixture of Ar39 "blips" and some ideal,
// straight-line MIP tracks.
//
// Output is a Python numpy .npz file.

local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';

local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local params_maker = import 'pgrapher/experiment/dune10kt-1x2x2/simparams.jsonnet';
local params = params_maker({});

local tools = tools_maker(params);

local sim_maker = import 'pgrapher/experiment/dune10kt-1x2x2/sim.jsonnet';
local sim = sim_maker(params, tools);

local stubby = {
  tail: wc.point(100.0, 300.0,   0.0, wc.cm),
  head: wc.point(100.0, 300.0, 230.0, wc.cm),
};

local tracklist = [
  {
    time: 0,
    charge: -5000,
    ray: stubby, // params.det.bounds,
  },
];
local output = 'wct-sim-ideal-sig.npz';


//local depos = g.join_sources(g.pnode({type:"DepoMerger", name:"BlipTrackJoiner"}, nin=2, nout=1),
//                             [sim.ar39(), sim.tracks(tracklist)]);
local depos = sim.tracks(tracklist);


//local deposio = io.numpy.depos(output);
local drifter = sim.drifter;
local bagger = sim.make_bagger();
// signal plus noise pipelines
//local sn_pipes = sim.signal_pipelines;
local sn_pipes = sim.splusn_pipelines;

local multimagnify = import 'pgrapher/experiment/dune10kt-1x2x2/multimagnify.jsonnet';
local magoutput = 'dune10kt-1x2x2-sim-check.root';
// please remove the root file before you generate a new one


local rootfile_creation_depos = g.pnode({
    type: "RootfileCreation_depos",
    name: "origmag",
    data: {
        output_filename: magoutput,
        root_file_mode: "RECREATE",
    },
}, nin=1, nout=1);


local multi_magnify = multimagnify('orig', tools, magoutput);
local magnify_pipes = multi_magnify.magnify_pipelines;
local multi_magnify2 = multimagnify('raw', tools, magoutput);
local magnify_pipes2 = multi_magnify2.magnify_pipelines;
local multi_magnify3 = multimagnify('gauss', tools, magoutput);
local magnify_pipes3 = multi_magnify3.magnify_pipelines;
local multi_magnify4 = multimagnify('wiener', tools, magoutput);
local magnify_pipes4 = multi_magnify4.magnify_pipelines;
local multi_magnify5 = multimagnify('threshold', tools, magoutput);
local magnify_pipes5 = multi_magnify5.magnifysummaries_pipelines;

local perfect = import 'chndb-perfect.jsonnet';
local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  data: perfect(params, tools.anodes[n], tools.field){dft:wc.tn(tools.dft)},
  uses: [tools.anodes[n], tools.field, tools.dft],
} for n in std.range(0, std.length(tools.anodes) - 1)];

//local chndb_maker = import 'pgrapher/experiment/pdsp/chndb.jsonnet';
//local noise_epoch = "perfect";
//local noise_epoch = "after";
//local chndb_pipes = [chndb_maker(params, tools.anodes[n], tools.fields[n]).wct(noise_epoch)
//                for n in std.range(0, std.length(tools.anodes)-1)];
local nf_maker = import 'pgrapher/experiment/dune10kt-1x2x2/nf.jsonnet';
// local nf_pipes = [nf_maker(params, tools.anodes[n], chndb_pipes[n]) for n in std.range(0, std.length(tools.anodes)-1)];
local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n) for n in std.range(0, std.length(tools.anodes) - 1)];

local sp_maker = import 'pgrapher/experiment/dune10kt-1x2x2/sp.jsonnet';
local sp = sp_maker(params, tools);
//local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];
local sp_pipes = [sp.make_sigproc(tools.anodes[n], n) for n in std.range(0, std.length(tools.anodes) - 1)];

//local parallel_pipes = [g.pipeline([sn_pipes[n], magnify_pipes[n],
//                                nf_pipes[n], magnify_pipes2[n] ], "parallel_pipe_%d"%n)
//                    for n in std.range(0, std.length(tools.anodes)-1)];
local parallel_pipes = [
  g.pipeline([
               sn_pipes[n],
               magnify_pipes[n],
               // nf_pipes[n],
               // magnify_pipes2[n],
               sp_pipes[n],
               // magnify_pipes3[n],
               // magnify_pipes4[n],
               // magnify_pipes5[n],
             ],
             'parallel_pipe_%d' % n)
  for n in std.range(0, std.length(tools.anodes) - 1)
];
local outtags = ['raw%d' % n for n in std.range(0, std.length(tools.anodes) - 1)];
local parallel_graph = f.fanpipe('DepoSetFanout', parallel_pipes, 'FrameFanin', 'sn_mag_nf', outtags);


//local frameio = io.numpy.frames(output);
local sink = sim.frame_sink;

local graph = g.pipeline([depos, drifter, bagger, parallel_graph, sink]);

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
};

local cmdline = {
    type: "wire-cell",
    data: {
        plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio", "WireCellSigProc", "WireCellRoot"],
        apps: ["Pgrapher"]
    }
};

// Finally, the configuration sequence which is emitted.

[cmdline] + g.uses(graph) + [app]
