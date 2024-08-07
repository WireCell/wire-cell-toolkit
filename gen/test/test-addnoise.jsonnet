// A test of using AddNoise, backwards compatible to before the noise refactoring.

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

local tick = 0.5*wc.us;

// services
local svcs = {rng: { type: 'Random' }, dft: { type: 'FftwDFT' }};
local wires = {
    type: "WireSchemaFile",
    data: {
        filename: "protodune-wires-larsoft-v4.json.bz2",
    }};


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

local absurd = pg.pnode({
    type: 'SilentNoise',
    data: {
        noutputs: 1,
        nchannels: 2560,
    }}, nin=0, nout=1);

local reframer(nsamples) = pg.pnode({
    type: 'Reframer',
    data: {
        nticks: nsamples,
        anode: wc.tn(anode),
    }}, nin=1, nout=1, uses=[anode]);

local pdsp_noise_file = "protodune-noise-spectra-v1.json.bz2";

local noise_models(nsamples) = {

    empno: {
        type: "EmpiricalNoiseModel",
        name: "empno",
        data: {
            anode: wc.tn(anode),
            chanstat: "",
            spectra_file: pdsp_noise_file,
            nsamples: nsamples,
            period: tick,
            wire_length_scale: 1*wc.cm,
        }, uses: [anode]
    },
    inco: {
        type: "GroupNoiseModel",
        name: "inco",
        data: {
            // This can also be given as a JSON/Jsonnet file
            spectra: pdsp_noise_file,
            groups: import "test-noise-groups-pdsp.jsonnet",
            nsamples: nsamples,
            tick: tick,
        }
    }
};


local addnoise(model, nsamples) =
    pg.pnode({
        type: 'AddNoise',
        name: "",
        data: {
            dft: wc.tn(svcs.dft),
            rng: wc.tn(svcs.rng),
            model: wc.tn(model),
            nsamples: nsamples,
        }}, nin=1, nout=1, uses=[model, svcs.dft, svcs.rng]);

local digi = pg.pnode({
    type: "Digitizer",
    name: "",
    data: {
        anode: wc.tn(anode),
        gain: 1.0,
        resolution: 12,
        baselines: [wc.volt,wc.volt,wc.volt],
        fullscale: [0, 2*wc.volt],
        round: false,           // post 0.23
    }
}, nin=1, nout=1, uses=[anode]);

function(outfile="test-addnoise.tar.bz2", noise="empno", nsamples=4096)
    local int_nsamples = if std.type(nsamples) == "number" then nsamples else std.parseInt(nsamples);

    local nmodel = noise_models(int_nsamples)[noise];

    local sink = pg.pnode({
        type: "FrameFileSink",
        name: outfile,
        data: {
            outname: outfile,
            digitize: true,
        },
    }, nin=1, nout=0);
    local graph = pg.pipeline([absurd, reframer(int_nsamples), addnoise(nmodel, int_nsamples), digi, sink]);
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
