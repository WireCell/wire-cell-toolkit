// This is used by test-pdsp-simsn.bats.

local pg = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';
local io = import 'pgrapher/common/fileio.jsonnet';

// for channel noise "database"
local nf_maker = import 'pgrapher/experiment/pdsp/nf.jsonnet';
local sp_maker = import 'pgrapher/experiment/pdsp/sp.jsonnet';

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
  // cfg/pgrapher/common/params.jsonnet leaves elec.gain and
  // adc.resolution commented out as "read from fcl" — larsoft injects
  // them at runtime.  Standalone wcsonnet has no fcl injection, so
  // restore the historical PDSP defaults here.
  elec: super.elec {
    gain: 14.0 * wc.mV / wc.fC,
  },
  adc: super.adc {
    resolution: 12,
  },
};

local tools = tools_maker(params);

// We only use first anode.
local anode = tools.anodes[0];
local pirs = tools.pirs[0];

//
// Simulation
//
local sim_maker = import "pgrapher/common/sim/nodes.jsonnet";
local sim = sim_maker(params, tools);

local drifter = sim.drifter;
local setdrifter = pg.pnode({
    type: 'DepoSetDrifter',
    data: {
        drifter: "Drifter"
    }
}, nin=1, nout=1, uses=[drifter]);

local transform =  pg.pnode({
    type:'DepoTransform',
    name: "",
    data: {
        rng: wc.tn(tools.random),
        dft: wc.tn(tools.dft),
        anode: wc.tn(anode),
        pirs: std.map(function(pir) wc.tn(pir), pirs),
        fluctuate: false,   // minimize randomness
        drift_speed: params.lar.drift_speed,
        first_frame_number: 100, // just to be different....
        readout_time: params.sim.ductor.readout_time,
        start_time: params.sim.ductor.start_time,
        tick: params.daq.tick,
        nsigma: 3,
    },
}, nin=1, nout=1, uses=[anode, tools.random, tools.dft] + pirs);
local reframer = pg.pnode({
    type: 'Reframer',
    name: '',
    data: {
        anode: wc.tn(anode),
        tbin: params.sim.reframer.tbin, // 120
        nticks: params.sim.reframer.nticks,
    },
}, nin=1, nout=1);
local noise_model = {
    type: "EmpiricalNoiseModel",
    name: "",
    data: {
        anode: wc.tn(anode),
        dft: wc.tn(tools.dft),
        chanstat: "",
        spectra_file: params.files.noise,
        nsamples: params.daq.nticks,
        period: params.daq.tick,
        wire_length_scale: 1.0*wc.cm, // optimization binning
    },
    uses: [anode, tools.dft],
};
local addnoise = pg.pnode({
    type: "AddNoise",
    name: "",
    data: {
        rng: wc.tn(tools.random),
        dft: wc.tn(tools.dft),
        model: wc.tn(noise_model),
	nsamples: params.daq.nticks,
        replacement_percentage: 0.02, // random optimization
    }}, nin=1, nout=1, uses=[tools.random, tools.dft, noise_model]);

local digitizer = pg.pnode({
    type: "Digitizer",
    name: "",
    data : params.adc {
        anode: wc.tn(anode),
        frame_tag: "orig0",
    }
}, nin=1, nout=1, uses=[anode]);

// Break up pipeline between signal voltage and adding noise so we can
// save out individual.
local vlt = pg.pipeline([setdrifter, transform, reframer]);
local noi = pg.pipeline([addnoise, digitizer]);

//
// Noise Filter
//
local chndb_data = import 'pgrapher/experiment/pdsp/chndb-base.jsonnet';

local chndb = {
  type: 'OmniChannelNoiseDB',
  name: '',
    data: chndb_data(params, anode, tools.field, 0){
        dft:wc.tn(tools.dft),
    },
  uses: [anode, tools.field, tools.dft]};


local nf = nf_maker(params, anode, chndb, 0, "");

//
// Signal Processing
//
local sp = sp_maker(params, tools).make_sigproc(anode);

// I/O

local depo_source(input)  = pg.pnode({
    type: 'DepoFileSource',
    name: wc.basename(input),
    data: { inname: input } // "depos.tar.bz2"
}, nin=0, nout=1);

local sio_sink(output, tags=[]) = pg.pnode({
    type: "FrameFileSink",
    name: wc.basename(output),
    data: {
        outname: output, // "frames.tar.bz2",
        tags: tags,
        digitize: false,
    },
}, nin=1, nout=0);

local sio_tap(output, tags=[]) = pg.fan.tap('FrameFanout', sio_sink(output, tags), wc.basename(output));

local plugins = [ "WireCellSio", "WireCellGen", "WireCellSigProc","WireCellApps", "WireCellPgraph", "WireCellTbb"];

function(input="input.npz",
         vltoutput="frames-vlt.npz", adcoutput="frames-adc.npz",
         nfoutput="frames-nf.npz", spoutput="frames-sp.npz")

    local graph = pg.pipeline([depo_source(input),
                               vlt, sio_tap(vltoutput),
                               noi, sio_tap(adcoutput, ["orig0"]),
                               nf, sio_tap(nfoutput, ["raw0"]),
                               sp, sio_sink(spoutput, ["gauss0","wiener0"])
                              ]);

    // Pgrapher or TbbFlow
    local engine = "Pgrapher";
    local app = {
        type: 'Pgrapher',
        data: {
            edges: pg.edges(graph),
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
    [cmdline] + pg.uses(graph) + [app]
