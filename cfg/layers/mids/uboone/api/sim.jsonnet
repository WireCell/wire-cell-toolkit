// The trifurcation requires the depo stream to split based on depo
// location relative to the wires.  Thus, we MUST assume individial
// depos

local low = import "../../../low.jsonnet";
local pg = low.pg;
local wc = low.wc;

local frs = import "frs.jsonnet";
local sr = import "shorted-regions.jsonnet";

local idents = low.util.idents;

function(services, params, options={}) {

    // Signal binning may be extended from nominal.
    local sig_binning = params.ductor.binning,

    // Must be *3* field responses
    local fr = frs(params).sims,

    local broken_detector = [{
        regions: sr.uv+sr.vy,
        mode:"reject",
        name:"nominal",
    },{
        regions: sr.uv,
        mode:"accept",
        name:"shorteduv",
    },{
        regions: sr.vy,
        mode:"accept",
        name:"shortedvy"
    }],

    local cer = {
        type: "ColdElecResponse",
        data: params.elec + sig_binning,
    },

    local rc = {
        type: 'RCResponse',
        data: sig_binning { width: params.rc.width }
    },


    // "sys status" is false for nominal.  No sys_resp.  
    local short_responses = [ cer, ],
    local long_responses = [ rc, rc ],

    local make_pirs = function(index) [{
        type: 'PlaneImpactResponse',
        name: "p"+std.toString(plane)+"i"+std.toString(index),
        uses: [fr[index], services.dft] + short_responses + long_responses,
        data: sig_binning {
            plane: plane,
            dft: wc.tn(services.dft),
            field_response: wc.tn(fr[index]),
            short_responses: [wc.tn(sr) for sr in short_responses],
            long_responses: [wc.tn(lr) for lr in long_responses],
            overall_short_padding: 100*wc.us,
            long_padding: 1.5*wc.ms,
        },
    } for plane in [0,1,2]],

    local transforms = function(anode, index, pirs) {
        local common = {
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
        local uses = pirs + [anode, services.random, services.dft],
        standard: 
            pg.pnode({
                type: 'DepoTransform',
                name: "a"+idents(anode)+"i"+std.toString(index),
                data: common,
            }, nin=1, nout=1, uses = uses),
        // BIG FAT WARNING: this will not work in standard WCT and
        // requires std.extVar's to be set.....
        overlay: 
            pg.pnode({
                type:'wclsReweightedDepoTransform',
                name: "a"+idents(anode)+"i"+std.toString(index),
                data: common + params.ductor.overlay,
            }, nin=1, nout=1, uses = uses),
    },
                
    local make_branch = function(anode, index)
        local pirs = make_pirs(index);
        pg.pipeline([
            pg.pnode({
                type: 'WireBoundedDepos',
                name: broken_detector[index].name,
                data: {
                    anode: wc.tn(anode),
                    regions: broken_detector[index].regions,
                    mode: broken_detector[index].mode,
                },
            }, nin=1, nout=1, uses = [anode]),
            pg.pnode({
                type: 'Bagger',
                name: broken_detector[index].name,
                data: {
                    gate: [params.ductor.start_time,
                           params.ductor.start_time+params.ductor.readout_time],
                },
            }, nin=1, nout=1),
            transforms(anode, index, pirs)[params.ductor.transform],
        ]),


    local sys = {
        type: 'ResponseSys',
        data: sig_binning {
            // overall_short_padding should take into account this offset "start".
            // currently all "start" should be the same cause we only have an overall time offset
            // compensated in reframer
            // These values correspond to files.fields[0, 1, 2]
            // e.g. normal, shorted U, and shorted Y
            start: [-10*wc.us, -10*wc.us, -10*wc.us], 
            magnitude: [1.0, 1.0, 1.0],
            time_smear: [0.0*wc.us, 0.0*wc.us, 0.0*wc.us],
        }
    },

    // Make some depos.  Note, this does NOT have a bagger.
    track_depos :: function(tracklist = [{
        time: 0,
        charge: -5000,         
        ray:  {
            tail: wc.point(-4.0, 0.0, 0.0, wc.m),
            head: wc.point(+4.0, 6.1, 7.0, wc.m),
        }}])
        low.gen.track_depos(tracklist),

    local ubsigtags = ['ubsig%d'%n for n in [0,1,2]],

    // API method sim.signal: subgraph making pure signal voltage from
    // depos.  Note, you MUST provide a stream of individial depos,
    // not a stream of deposets!
    signal :: function(anode) pg.pipeline([
        pg.fan.pipe('DepoFanout', [make_branch(anode, index) for index in [0,1,2]],
                    'FrameFanin', 'ubsigraph', ubsigtags),
        pg.pnode({
            type: 'Reframer',
            name: idents(anode),
            data: params.reframer {
                anode: wc.tn(anode),
                tags: ubsigtags,
                fill: 0.0,
            },
        }, nin=1, nout=1, uses=[anode])]),

    
    local csdb = {
        type: "StaticChannelStatus",
        name: "",
        data: {
            nominal_gain: params.elec.gain,
            nominal_shaping: params.elec.shaping,
            deviants:  [ {
                chid: ch,
                gain: params.misconfigured.gain,
                shaping: params.misconfigured.shaping,
            } for ch in params.misconfigured.channels ],
        }
    },

    noise :: function(anode) 
        local model = {
            type: 'EmpiricalNoiseModel',
            name: idents(anode),
            data: params.noise.model {
                anode: wc.tn(anode),
                dft: wc.tn(services.dft),
                chanstat:wc.tn(csdb),
            },
            uses: [anode, csdb, services.dft]
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
