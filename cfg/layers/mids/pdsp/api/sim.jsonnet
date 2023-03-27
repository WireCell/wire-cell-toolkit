// The .sim sub-API produces a few sub graphs related to simulation so
// end user can achieve different goals.

// Nothing in here should encode any variants.  Put any variant
// information in variants/*.jsonnet.  Information from one of those
// files is passed as "params"

local low = import "../../../low.jsonnet";
// some short hands
local wc = low.wc;
local pg = low.pg;
local idents = low.util.idents;

local frs = import "frs.jsonnet";

function(services, params, options={}) {

    // Signal binning may be extended from nominal.
    local sig_binning = params.ductor.binning,

    local fr = frs(params).sim,

    local cer = {
        type: "ColdElecResponse",
        data: params.elec + sig_binning,
    },

    local rc = {
        type: 'RCResponse',
        data: sig_binning { width: params.rc.width }
    },

    // some have more than one
    local short_responses = [ cer, ],
    local long_responses = [ rc, ],
        
    local pirs = [{
        type: 'PlaneImpactResponse',
        name: std.toString(plane),
        uses: [fr, services.dft] + short_responses + long_responses,
        data: sig_binning {
            plane: plane,
            dft: wc.tn(services.dft),
            field_response: wc.tn(fr),
            short_responses: [wc.tn(sr) for sr in short_responses],
            long_responses: [wc.tn(lr) for lr in long_responses],
            overall_short_padding: 100*wc.us,
            long_padding: 1.5*wc.ms,
        },
    } for plane in [0,1,2]],

    // API method sim.track_depos: subgraph making some depos from
    // ideal tracks in the detector.
    track_depos :: function(tracklist = [{
        time: 0,
        charge: -5000,         
        ray:  {
            tail: wc.point(-4.0, 0.0, 0.0, wc.m),
            head: wc.point(+4.0, 6.1, 7.0, wc.m),
        }}])
        pg.pipeline([
            low.gen.track_depos(tracklist),
            low.gen.bagger(params.ductor.readout_time,
                           params.ductor.start_time)
        ]),
        

    // API method sim.reframer
    reframer :: function(anode)
        pg.pnode({
            type: 'Reframer',
            name: idents(anode),
            data: {
                anode: wc.tn(anode),
                tags: [],
                fill: 0.0,
                toffset: 0,
                tbin: params.ductor.tbin,
                nticks: sig_binning.nticks - self.tbin,
            },
        }, nin=1, nout=1, uses=[anode]),
    
    // API method sim.signal: subgraph making pure signal voltage from
    // depos.
    signal :: function(anode)
        pg.pipeline([
            pg.pnode({
                type: 'DepoTransform',
                name: idents(anode),
                data: {
                    rng: wc.tn(services.random),
                    dft: wc.tn(services.dft),
                    anode: wc.tn(anode),
                    pirs: [wc.tn(p) for p in pirs],
                    fluctuate: true,
                    drift_speed: params.lar.drift_speed,
                    first_frame_number: 0,
                    readout_time: params.ductor.readout_time,
                    start_time: params.ductor.start_time,
                    tick: sig_binning.tick,
                    nsigma: 3,
                },
            }, nin=1, nout=1, uses = pirs + [anode, services.random, services.dft]),
            $.reframer(anode)]),

    // API method sim.noise: subgraph adding noise to voltage
    noise :: function(anode)
        local model = {
            type: 'EmpiricalNoiseModel',
            name: idents(anode),
            data: params.noise.model {
                anode: wc.tn(anode),
                dft: wc.tn(services.dft),
                chanstat:"",    // must explicitly empty
            },
            uses: [anode, services.dft]
        };
        pg.pnode({
            type: 'AddNoise',
            name: idents(anode),
            data: {
                rng: wc.tn(services.random),
                dft: wc.tn(services.dft),
                model: wc.tn(model),
                nsamples: params.noise.model.nsamples,
                replacement_percentage: params.noise.replacement_percentage,
            }}, nin=1, nout=1, uses=[services.random, services.dft, model]),

    // API method sim.digitizer: return subgraph adding digitization
    // of voltage to produce ADC
    digitizer :: function(anode)
        pg.pnode({
            type: 'Digitizer',
            name: idents(anode),
            data: params.digi {
                anode: wc.tn(anode),
                frame_tag: "orig" + idents(anode),
            }
        }, nin=1, nout=1, uses=[anode]),

}
