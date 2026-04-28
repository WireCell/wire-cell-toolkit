local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";
local f = import "pgrapher/common/funcs.jsonnet";
local sim_maker = import "pgrapher/common/sim/nodes.jsonnet";


// return some nodes, includes base sim nodes.
function(params, tools) {
    local sim = sim_maker(params, tools),

    local nanodes = std.length(tools.anodes),


    local ductors = [sim.make_depotransform("ductor%d"%anode.data.ident, anode, tools.pirs[0]) for anode in tools.anodes],

    local reframers = [
        g.pnode({
            type: 'Reframer',
            name: 'reframer%d'%anode.data.ident,
            data: {
                anode: wc.tn(anode),
                tags: [],           // ?? what do?
                fill: 0.0,
                tbin: params.sim.reframer.tbin,
                toffset: 0,
                nticks: params.sim.reframer.nticks,
            },
        }, nin=1, nout=1) for anode in tools.anodes],

    // fixme: see https://github.com/WireCell/wire-cell-gen/issues/29
    local make_noise_model = function(anode, csdb=null) {
        type: "EmpiricalNoiseModel",
        name: "empericalnoise%s"% anode.name,
        data: {
            anode: wc.tn(anode),
            dft: wc.tn(tools.dft),
            chanstat: if std.type(csdb) == "null" then "" else wc.tn(csdb),
            spectra_file: params.files.noise,
            nsamples: params.daq.nticks,
            period: params.daq.tick,
            wire_length_scale: 1.0*wc.cm, // optimization binning
        },
        uses: [anode, tools.dft] + if std.type(csdb) == "null" then [] else [csdb],
    },
    local noise_models = [make_noise_model(anode) for anode in tools.anodes],


    local add_noise = function(model) g.pnode({
        type: "AddNoise",
        name: "addnoise%s"%[model.name],
        data: {
            rng: wc.tn(tools.random),
            dft: wc.tn(tools.dft),
            model: wc.tn(model),
	    nsamples: params.daq.nticks,
            replacement_percentage: 0.02, // random optimization
        }}, nin=1, nout=1, uses=[tools.random, tools.dft, model]),

    local noises = [add_noise(model) for model in noise_models],

    local digitizers = [
        sim.digitizer(anode, name="digitizer%d"%anode.data.ident, tag="orig%d"%anode.data.ident)
        for anode in tools.anodes],

    ret : {

        signal_pipelines: [g.pipeline([ductors[n], reframers[n],  digitizers[n]],
                                      name="simsigpipe%d"%tools.anodes[n].data.ident,) for n in std.range(0, nanodes-1)],

        splusn_pipelines: [g.pipeline([ductors[n], reframers[n], noises[n],  digitizers[n]],
                                      name="simsigpipe%d"%tools.anodes[n].data.ident) for n in std.range(0, nanodes-1)],


    } + sim,      
}.ret
