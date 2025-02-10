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

local input = std.extVar('input');
local output = std.extVar('output');

local tools_maker = import 'pgrapher/common/tools.jsonnet';
local base = import 'pgrapher/experiment/pdsp/simparams.jsonnet';
local params = base {
  lar: super.lar {
    // Longitudinal diffusion constant
    // DL: std.extVar('DL') * wc.cm2 / wc.s,
    DL: 4.0 * wc.cm2 / wc.s,
    // Transverse diffusion constant
    // DT: std.extVar('DT') * wc.cm2 / wc.s,
    DT: 8.8 * wc.cm2 / wc.s,
    // Electron lifetime
    // lifetime: std.extVar('lifetime') * wc.ms,
    lifetime: 10.4 * wc.ms,
    // Electron drift speed, assumes a certain applied E-field
    // drift_speed: std.extVar('driftSpeed') * wc.mm / wc.us,
    drift_speed: 1.565 * wc.mm / wc.us,
  },
};

local tools_all = tools_maker(params);
local tools = tools_all {
    anodes: [tools_all.anodes[0]]
};

local sim_maker = import 'pgrapher/experiment/pdsp/sim.jsonnet';
local sim = sim_maker(params, tools);

local nanodes = std.length(tools.anodes);
local anode_iota = std.range(0, nanodes - 1);

//local deposio = io.numpy.depos(output);
local bagger = sim.make_bagger();

local drifter = sim.drifter;
local setdrifter = g.pnode({
            type: 'DepoSetDrifter',
            data: {
                drifter: "Drifter"
            }
        }, nin=1, nout=1,
        uses=[drifter]);

// signal plus noise pipelines
//local sn_pipes = sim.signal_pipelines;
local sn_pipes = sim.splusn_pipelines;

local perfect = import 'pgrapher/experiment/pdsp/chndb-perfect.jsonnet';
local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  data: perfect(params, tools.anodes[n], tools.field, n){dft:wc.tn(tools.dft)},
  uses: [tools.anodes[n], tools.field, tools.dft],
} for n in anode_iota];

//local chndb_maker = import 'pgrapher/experiment/pdsp/chndb.jsonnet';
//local noise_epoch = "perfect";
//local noise_epoch = "after";
//local chndb_pipes = [chndb_maker(params, tools.anodes[n], tools.fields[n]).wct(noise_epoch)
//                for n in std.range(0, std.length(tools.anodes)-1)];
local nf_maker = import 'pgrapher/experiment/pdsp/nf.jsonnet';
// local nf_pipes = [nf_maker(params, tools.anodes[n], chndb_pipes[n]) for n in std.range(0, std.length(tools.anodes)-1)];
local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n) for n in anode_iota];

local sp_maker = import 'pgrapher/experiment/pdsp/sp.jsonnet';
local sp_override = {
    sparse: false,
    use_roi_debug_mode: false,
    m_save_negative_charge: false,
    use_multi_plane_protection: false,
    mp_tick_resolution: 10,
};
local sp = sp_maker(params, tools, sp_override);
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];


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

local magoutput = 'mag.root';
local magnify = import 'pgrapher/experiment/dune-vd/magnify-sinks.jsonnet';
local sinks = magnify(tools, magoutput);

local sio_sinks = g.pnode({
        type: "FrameFileSink",
        data: {
            outname: output, // "frames.tar.bz2",
            tags: ["orig", "gauss"],
            digitize: false,
        },
    }, nin=1, nout=0);

local hio_sp = [g.pnode({
      type: 'HDF5FrameTap',
      name: 'hio_sp%d' % n,
      data: {
        anode: wc.tn(tools.anodes[n]),
        trace_tags: ['loose_lf%d' % n 
        , 'tight_lf%d' % n 
        , 'cleanup_roi%d' % n 
        , 'break_roi_1st%d' % n 
        , 'break_roi_2nd%d' % n 
        , 'shrink_roi%d' % n 
        , 'extend_roi%d' % n 
        , 'mp3_roi%d' % n 
        , 'mp2_roi%d' % n 
        , 'decon_charge%d' % n 
        , 'gauss%d' % n], 
        filename: "g4-rec-%d.h5" % n,
        // chunk: [1, 1], // ncol, nrow
        gzip: 2,
        high_throughput: true,
      },
    }, nin=1, nout=1),
    for n in std.range(0, std.length(tools.anodes) - 1)
];

local multipass = [
    g.pipeline([
        sn_pipes[n],
        // sinks.orig_pipe[n],
        // sio_sinks[n],
        // nf_pipes[n],
        sp_pipes[n],
        // hio_sp[n],
        sinks.decon_pipe[n],
        sinks.debug_pipe[n],
        // sinks.threshold_pipe[n],
        // dnnroi(tools.anodes[n], ts, output_scale=1.0),
    ], 'multipass%d' % n)
  for n in anode_iota
];
local outtags = ['orig%d' % n for n in anode_iota];
local bi_manifold = f.fanpipe('DepoSetFanout', multipass, 'FrameFanin', 'sn_mag_nf', outtags);
// local bi_manifold = f.fanpipe('DepoFanout', multipass, 'FrameFanin', 'sn_mag_nf', outtags);

local retagger = g.pnode({
  type: 'Retagger',
  data: {
    // Note: retagger keeps tag_rules an array to be like frame fanin/fanout.
    tag_rules: [{
      // Retagger also handles "frame" and "trace" like fanin/fanout
      // merge separately all traces like gaussN to gauss.
      frame: {
        '.*': 'orig',
      },
      merge: {
        'orig\\d': 'daq',
      },
    }],
  },
}, nin=1, nout=1);

//local frameio = io.numpy.frames(output);
local sink = sim.frame_sink;


local depo_source  = g.pnode({
    type: 'DepoFileSource',
    data: { inname: input } // "depos.tar.bz2"
}, nin=0, nout=1);

local graph = g.pipeline([depo_source, setdrifter, bi_manifold, retagger, sio_sinks]);
local plugins = [ "WireCellSio", "WireCellGen", "WireCellSigProc","WireCellApps", "WireCellPgraph", "WireCellTbb", "WireCellRoot", "WireCellHio", "WireCellPytorch"];

// Pgrapher or TbbFlow
local engine = "Pgrapher";
local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
};

local cmdline = {
    type: "wire-cell",
    data: {
        plugins: plugins,
        apps: ["Pgrapher"],
    }
};


// Finally, the configuration sequence which is emitted.
[cmdline] + g.uses(graph) + [app]
