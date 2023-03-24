// Microboone-specific functions for simulation related things.  The
// structure of function is parameterized on params and tools.

local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";
local simnodes = import "pgrapher/common/sim/nodes.jsonnet";
local f = import "pgrapher/common/funcs.jsonnet";

function(params, tools)
{
    local sim = simnodes(params, tools),

    local sr = params.shorted_regions,

    // A depo set oriented subgraph for WireBoundedDepos.
    local wbdepos_sets = [
        sim.deposet_filter(
            sim.make_wbdepo(tools.anode, sr.uv + sr.vy, "reject", "nominal"),
            "nominal"),
        sim.deposet_filter(
            sim.make_wbdepo(tools.anode, sr.uv, "accept", "shorteduv"),
            "shoteduv"),
        sim.deposet_filter(
            sim.make_wbdepo(tools.anode, sr.vy, "accept", "shortedvy"),
            "shortedvy"),
    ],

    // The older per-depo WireBoundedDepos subgraph, best to use the
    // depo set version above
    local wbdepos = [
        sim.make_wbdepo(tools.anode, sr.uv + sr.vy, "reject", "nominal"),
        sim.make_wbdepo(tools.anode, sr.uv, "accept", "shorteduv"),
        sim.make_wbdepo(tools.anode, sr.vy, "accept", "shortedvy"),
    ],

    local baggers = [sim.make_bagger("bagger%d"%n) for n in [0,1,2]],
        
    local zippers = [sim.make_depozipper("depozip%d"%n, tools.anode, tools.pirs[n]) for n in [0,1,2]],
    local transforms = [sim.make_depotransform("depotran%d"%n, tools.anode, tools.pirs[n]) for n in [0,1,2]],

    //local depos2frames = zippers,
    local depos2frames = transforms,

    local pipelines = [g.pipeline([wbdepos[n], baggers[n], depos2frames[n]], "ubsigpipe%d"%n)
                       for n in [0,1,2]],
    local pipelines_sets = [g.pipeline([wbdepos_sets[n], depos2frames[n]], "ubsigpipe%d"%n)
                            for n in [0,1,2]],

    
    local ubsigtags = ['ubsig%d'%n for n in [0,1,2]],

    local signal = g.pipeline([f.fanpipe('DepoFanout', pipelines, 'FrameFanin', 'ubsigraph', ubsigtags),
                               sim.make_reframer('ubsigrf', tools.anode, ubsigtags)], 'ubsignal'),

    local signal_sets = g.pipeline([f.fanpipe('DepoSetFanout', pipelines_sets, 'FrameFanin', 'ubsigraph', ubsigtags),
                                    sim.make_reframer('ubsigrf', tools.anode, ubsigtags)], 'ubsignal'),

    //
    // Noise:
    //
    
    // A channel status with no misconfigured channels.
    local empty_csdb = {
        type: "StaticChannelStatus",
        name: "uniform",
        data: {
            nominal_gain: params.elec.gain,
            nominal_shaping: params.elec.shaping,
            deviants: [
                //// This is what elements of this array look like:
                //// One entry per "deviant" channel.
                // {
                //     chid: 0,               // channel number
                //     gain: 4.7*wc.mV/wc.fC, // deviant gain
                //     shaping: 1*wc.us,      // deviant shaping time
                // }
            ],
        }
    },


    // A channel status configured with nominal misconfigured channels.
    local miscfg_csdb = {
        type: "StaticChannelStatus",
        name: "misconfigured",
        data: {
            nominal_gain: params.elec.gain,
            nominal_shaping: params.elec.shaping,
            deviants: [ { chid: ch,
                                gain: params.nf.misconfigured.gain,
                                shaping: params.nf.misconfigured.shaping,
                              } for ch in params.nf.misconfigured.channels ],
        },
    },

        
    // Make a noise model bound to an anode and a channel status
    local make_noise_model = function(anode, csdb) {
        type: "EmpiricalNoiseModel",
        name: "empericalnoise%s"% csdb.name,
        data: {
            anode: wc.tn(anode),
            dft: wc.tn(tools.dft),
            chanstat: wc.tn(csdb),
            spectra_file: params.files.noise,
            nsamples: params.daq.nticks,
            period: params.daq.tick,
            wire_length_scale: 1.0*wc.cm, // optimization binning
        },
        uses: [anode, csdb, tools.dft],
    },


    // this is a frame filter that adds noise
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

    ret: {
        signal : signal,
        signal_sets : signal_sets,
        empty_csdb : empty_csdb,
        miscfg_csdb : miscfg_csdb,
        make_noise_model :: make_noise_model,
        add_noise:: add_noise,
    } + sim,
    


}.ret
