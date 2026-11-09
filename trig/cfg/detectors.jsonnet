// Per-detector setup for the SNB TA study: the parameters, the one anode
// simulated and its simulation chain (DepoTransform, Reframer, optional
// AddNoise, Digitizer), with defaults for the trigger.
//
//   local det = (import "detectors.jsonnet")[name](anode, nticks, lifetime, noise)
//
// gives {params, tools, anode, drifter, simpipe, thresholds, plane_offsets_us,
// geometry}.  Names:
//
// - hd :: DUNE FD horizontal drift, one APA of dune10kt-1x2x6 (both faces,
//   2560 channels, 14-bit ADC).  Anode index 0..11, default 0.
// - vd :: DUNE FD vertical drift, one top-drift CRU of the v6 geometry
//   (two CRMs with jumpered induction strips, 1536 channels, 12-bit ADC,
//   TDE electronics at 500 ns).  Anode ident 160..319, default 160.

local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";
local tools_maker = import "pgrapher/common/tools.jsonnet";
local sim_nodes = import "pgrapher/common/sim/nodes.jsonnet";

// The drifter as a component (the common one is only a pnode).
local make_drifter(params, tools) = {
    type: "Drifter",
    data: params.lar {
        rng: wc.tn(tools.random),
        xregions: std.filter(function(x) x != null,
                             wc.unique_list(std.flattenArrays([v.faces for v in params.det.volumes]))),
        time_offset: params.sim.depo_toffset,
        drift_speed: params.lar.drift_speed,
        fluctuate: params.sim.fluctuate,
        DL: params.lar.DL,
        DT: params.lar.DT,
        lifetime: params.lar.lifetime,
    },
    uses: [tools.random],
};

// DepoTransform -> Reframer -> [AddNoise] -> Digitizer for one anode.
local make_simpipe(params, tools, anode, pirs, noise_file, noise) =
    local sim = sim_nodes(params, tools);
    local name = anode.name;
    local reframer = g.pnode({
        type: "Reframer",
        name: "reframer" + name,
        data: {
            anode: wc.tn(anode),
            tags: [],
            fill: 0.0,
            tbin: params.sim.reframer.tbin,
            toffset: 0,
            nticks: params.sim.reframer.nticks,
        },
    }, nin=1, nout=1, uses=[anode]);
    local model = {
        type: "EmpiricalNoiseModel",
        name: "noise" + name,
        data: {
            anode: wc.tn(anode),
            dft: wc.tn(tools.dft),
            chanstat: "",
            spectra_file: noise_file,
            nsamples: params.daq.nticks,
            period: params.daq.tick,
            wire_length_scale: 1.0*wc.cm,
        },
        uses: [anode, tools.dft],
    };
    local addnoise = g.pnode({
        type: "AddNoise",
        name: "addnoise" + name,
        data: {
            rng: wc.tn(tools.random),
            dft: wc.tn(tools.dft),
            model: wc.tn(model),
            nsamples: params.daq.nticks,
            replacement_percentage: 0.02,
        },
    }, nin=1, nout=1, uses=[tools.random, tools.dft, model]);
    g.pipeline([sim.make_depotransform("ductor" + name, anode, pirs), reframer]
               + (if noise then [addnoise] else [])
               + [sim.digitizer(anode, name="digitizer" + name, tag="orig")]);

local with_lifetime(params, lifetime) =
    if lifetime == null then params
    else params { lar: super.lar { lifetime: lifetime } };

{
    hd(anode=0, nticks=6000, lifetime=null, noise=true)::
        local params0 = (import "pgrapher/experiment/dune10kt-1x2x6/simparams.jsonnet")({});
        local params = with_lifetime(params0 { daq: super.daq { nticks: nticks } }, lifetime);
        local tools = tools_maker(params);
        local the_anode = tools.anodes[anode];
        {
            name: "hd",
            geometry: "DUNE FD HD 1x2x6, APA %d" % anode,
            params: params,
            tools: tools,
            anode: the_anode,
            drifter: make_drifter(params, tools),
            simpipe: make_simpipe(params, tools, the_anode, tools.pirs[0], params.files.noise, noise),
            // dunetrigger HD SimpleThreshold
            thresholds: [70, 70, 60],
            // measured with the offset_est column (v1)
            plane_offsets_us: [59, 61.5, 64],
            drift_length: 3.63*wc.m,
        },

    vd(anode=160, nticks=6000, lifetime=null, noise=true)::
        local params0 = (import "pgrapher/experiment/dune-vd/params-10kt.jsonnet")({
            response_plane: 18.1*wc.cm, nticks: nticks, ncrm: 320});
        local params = with_lifetime(params0 {
            files: super.files { wires: "dunevd10kt_3view_30deg_v6_refactored.json.bz2" },
        }, lifetime);
        local tools = tools_maker(params);
        local the_anode = tools.anodes[anode];
        local top = anode >= 160;
        {
            name: "vd",
            geometry: "DUNE FD VD v6, %s CRU %d" % [if top then "top" else "bottom", anode],
            params: params,
            tools: tools,
            anode: the_anode,
            drifter: make_drifter(params, tools),
            simpipe: make_simpipe(params, tools, the_anode, tools.pirs[if top then 1 else 0],
                                  params.files.noises[if top then 1 else 0], noise),
            // dunetrigger VD SimpleThreshold (4 sigma)
            thresholds: [36, 36, 36],
            // measured with the offset_est column
            plane_offsets_us: [112, 114.6, 118.3],
            drift_length: 6.5*wc.m,
        },
}
