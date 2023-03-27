// Test the "roundtrip":
//
// noise: file -> sim -> modeler -> file

// CAUTION: this Jsonnet file does not follow good configuration
// style.  In order to be be self-contained it is overly verbose,
// underly generic and not at all extensible.  DO NOT USE IT AS AN
// EXAMPLE.

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

// These define the noise spectra and groups.  We may either pass
// JSON/Jsonnet file NAMES or we may provide the data structures
// directly.  Since this tests algorithmically generates everything,
// it's easier to bring these inline.  Also, we happen to use the same
// spectra for both grouped coherent and grouped incoherent.  In real
// world, such a coincidence is very unlikely.
local tns = import "test-noise-spectra.jsonnet";
local tng = {
    inco: import "test-noise-groups-incoherent.jsonnet",
    cohe: import "test-noise-groups-coherent.jsonnet",
    empno: import "test-noise-groups-incoherent.jsonnet",
};

// Size of original (fictional) waveforms
local nsamples_modelin = 4096;
// and how many subsamples to store.
local nsamples_modelinsave = 64;

// Size of generated waveforms.
//local nsamples_generate = 6000;
local nsamples_generate = 4096;

// As above but for saving out the model.
local nsamples_modelout = 4096;
local nsamples_modeloutsave = 64;

// services
local svcs = {rng: { type: 'Random' }, dft: { type: 'FftwDFT' }};
local wires = {
    type: "WireSchemaFile",
    data: {
        filename: "protodune-wires-larsoft-v4.json.bz2",
    }};
local isnoise = {
    type: "NoiseRanker",
    data: {
        maxdev: 10*wc.mV,
    },
};

// pdsp apa 0, dumped from "real" jsonnet
local faces = [{
    "anode": -3578.36,
    "cathode": -1.5875,
    "response": -3487.8875
}, {
    "anode": -3683.14,
    "cathode": -7259.9125,
    "response": -3773.6125
}];

local anode = {
    type: "AnodePlane",
    name: "",
    data: {
        ident: 0,
        wire_schema: wc.tn(wires),
        faces: faces,
    },
    uses: [wires]};

// From nothing comes all.
local absurd = pg.pnode({
    type: 'SilentNoise',
    data: {
        noutputs: 1,
        nchannels: 2560,
    }}, nin=0, nout=1);

local reframer = pg.pnode({
    type: 'Reframer',
    data: {
        nticks: nsamples_generate,
        anode: wc.tn(anode),
    }}, nin=1, nout=1, uses=[anode]);


// same data for digitizing and undigitizing
local digidata = {
    anode: wc.tn(anode),
    gain: 1.0,
    resolution: 12,
    baselines: [wc.volt,wc.volt,wc.volt],
    fullscale: [0, 2*wc.volt],
};


local tick = 0.5*wc.us;

// The graph will have two major pipelines split based on coherent and
// incoherent noise.  Here are the short nick names to ID each:
local group_nicks = ["inco", "cohe"];
local adder_types = {
    inco: "IncoherentAddNoise",
    cohe: "CoherentAddNoise",
    empno: "IncoherentAddNoise",
};
local models = {
    [one]: {
        type: "GroupNoiseModel",
        name: one,
        data: {
            // This can also be given as a JSON/Jsonnet file
            spectra: tns(nsamples=nsamples_modelin, nsave=nsamples_modelinsave,
                         rms=1*wc.mV), // PDSP is eg 1.3 mV
            groups: tng[one],
            nsamples: nsamples_generate,
            tick: tick,
        }
    }
    for one in group_nicks} + {
        empno: {
            type: "EmpiricalNoiseModel",
            name: "empno",
            data: {
                anode: wc.tn(anode),
                chanstat: "",
                spectra_file: "protodune-noise-spectra-v1.json.bz2",
                nsamples: nsamples_generate,
                period: tick,
                wire_length_scale: 1*wc.cm,
            }, uses: [anode]
        },
    };


local pipes(round=true, bug202 = 0.0) = {
    [one]: [

        // add noise
        pg.pnode({
            type: adder_types[one],
            name: one,
            data: {
                dft: wc.tn(svcs.dft),
                rng: wc.tn(svcs.rng),
                model: wc.tn(models[one]),
                nsamples: nsamples_generate,
                bug202: bug202,
            }}, nin=1, nout=1, uses=[models[one], svcs.dft, svcs.rng]),

        // tap input to modeler
        pg.fan.tap('FrameFanout',  pg.pnode({
            type: "FrameFileSink",
            name: one+'vin',
            data: {
                outname: "test-noise-roundtrip-%s-vin.npz"%one, 
                digitize: false,
            },
        }, nin=1, nout=0), one+'vin'),

        // digitize
        pg.pnode({
            type: "Digitizer",
            name: one,
            data: digidata + {round: round}
        }, nin=1, nout=1, uses=[anode]),

        // tap generated noise
        pg.fan.tap('FrameFanout',  pg.pnode({
            type: "FrameFileSink",
            name: one+'adc',
            data: {
                outname: "test-noise-roundtrip-%s-adc.npz"%one, 
                digitize: true,
            },
        }, nin=1, nout=0), one+'adc'),

        // undigitize
        pg.pnode({
            type: "Undigitizer",
            name: one,
            data: digidata,
        }, nin=1, nout=1, uses=[anode]),

        // tap input to modeler
        pg.fan.tap('FrameFanout',  pg.pnode({
            type: "FrameFileSink",
            name: one+'dac',
            data: {
                outname: "test-noise-roundtrip-%s-dac.npz"%one, 
                digitize: false,
            },
        }, nin=1, nout=0), one+'dac'),
        
        // and finally the modeler itself
        pg.pnode({
            type: "NoiseModeler",
            name: one,
            data: {
                dft: wc.tn(svcs.dft),
                isnoise: wc.tn(isnoise),
                threshold: 0.9,
                groups: tng[one],
                outname: "test-noise-roundtrip-%s-spectra.json.bz2"%one
            },
        }, nin=1, nout=0, uses=[svcs.dft, isnoise])] // one pipe
    for one in std.objectFields(models)
};                              // pipes

function(round=true, bug202=0.0) 
    local bool_round = round == true || round == "true" || round == "round";
    local number_bug202 = if std.type(bug202) == "number" then bug202 else std.parseJson(bug202);

    local nicks_to_use = ["empno", "cohe", "inco"];
    local fanout = pg.fan.fanout('FrameFanout', [
        pg.pipeline(pipes(bool_round, number_bug202)[nick]) for nick in nicks_to_use
    ]);
    local graph = pg.pipeline([absurd, reframer, fanout]);

    local app = {
        type: 'Pgrapher',
        data: {
            edges: pg.edges(graph),
        },
    };
    local cmdline = {
        type: "wire-cell",
        data: {
            plugins: ["WireCellAux", "WireCellSigProc", "WireCellGen",
                      "WireCellApps", "WireCellPgraph", "WireCellSio"],
            apps: [app.type]
        }
    };
    [cmdline] + pg.uses(graph) + [app]

                                 
