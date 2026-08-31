local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";
local f = import "pgrapher/common/funcs.jsonnet";
local sim_maker = import "pgrapher/common/sim/nodes.jsonnet";


// return some nodes, includes base sim nodes.
function(params, tools) {
    local sim = sim_maker(params, tools),

    local nanodes = std.length(tools.anodes),

    // Drift-volume index for an anode: bottom CRMs (ident < half) -> 0, top
    // -> 1.  Clamped to the number of field responses so single-drift params
    // (one field/elec/noise) map every anode to index 0.
    local noise_files = if std.objectHas(params.files, "noises")
                        then params.files.noises else [params.files.noise],
    local ndrift = std.length(noise_files),
    local drift_idx = function(anode)
        local raw = if anode.data.ident < std.length(params.det.volumes) / 2
                    then 0 else 1;
        if raw < ndrift then raw else ndrift - 1,

    local ductors = [sim.make_depotransform("ductor%d"%anode.data.ident, anode, tools.pirs[drift_idx(anode)]) for anode in tools.anodes],

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
    local make_noise_model = function(anode, spectra, csdb=null) {
        type: "EmpiricalNoiseModel",
        name: "empericalnoise%s"% anode.name,
        data: {
            anode: wc.tn(anode),
            dft: wc.tn(tools.dft),
            chanstat: if std.type(csdb) == "null" then "" else wc.tn(csdb),
            spectra_file: spectra,
            nsamples: params.daq.nticks,
            period: params.daq.tick,
            wire_length_scale: 1.0*wc.cm, // optimization binning
        },
        uses: [anode, tools.dft] + if std.type(csdb) == "null" then [] else [csdb],
    },

    // Electronics-noise spectra synced to protodunevd: top and bottom drift
    // volumes use different spectra, so there is one shared noise model per
    // drift (noise_files = [bottom, top]; single-drift params has just one).
    local mega_anode = function(d) {
        type: 'MegaAnodePlane',
        name: 'meganodes%d' % d,
        data: {
            anodes_tn: [wc.tn(a) for a in tools.anodes if drift_idx(a) == d],
        },
        uses: [a for a in tools.anodes if drift_idx(a) == d],
    },

    local add_noise = function(model, name="") g.pnode({
        type: "AddNoise",
        name: "addnoise%s"%[name],
        data: {
            rng: wc.tn(tools.random),
            dft: wc.tn(tools.dft),
            model: wc.tn(model),
            nsamples: params.daq.nticks,
            replacement_percentage: 0.02, // random optimization
        }}, nin=1, nout=1, uses=[tools.random, tools.dft, model]),

    local noises = [
        add_noise(make_noise_model(mega_anode(drift_idx(anode)),
                                   noise_files[drift_idx(anode)]), anode.name)
        for anode in tools.anodes],

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
