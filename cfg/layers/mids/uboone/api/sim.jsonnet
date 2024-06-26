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

    // local res = low.resps(params).sim,

    // Signal binning may be extended from nominal.
    local sig_binning = params.ductor.binning,

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
        // A set of multiple FRs.  0-element is "nominal"
        local uboone_frs = frs(params).sims,
        type: 'PlaneImpactResponse',
        name: "p"+std.toString(plane)+"i"+std.toString(index),
        uses: [uboone_frs[index], services.dft] + short_responses + long_responses,
        data: sig_binning {
            plane: plane,
            dft: wc.tn(services.dft),
            field_response: wc.tn(uboone_frs[index]),
            short_responses: [wc.tn(sr) for sr in short_responses],
            long_responses: [wc.tn(lr) for lr in long_responses],
            overall_short_padding: 100*wc.us,
            long_padding: 1.5*wc.ms,
        },
    } for plane in [0,1,2]],

    // Make a DepoTransform parameterized by its FR.
    local transforms = function(anode, pirs, index=0) {
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
        standard: pg.pnode({
            type: 'DepoTransform',
            name: "a"+idents(anode)+"i"+std.toString(index),
            data: common,
        }, nin=1, nout=1, uses = uses),
        // BIG FAT WARNING: this will not work in standard WCT and
        // requires std.extVar's to be set.....
        overlay: pg.pnode({
            type:'wclsReweightedDepoTransform',
            name: "a"+idents(anode)+"i"+std.toString(index),
            data: common + params.ductor.overlay,
        }, nin=1, nout=1, uses = uses),
    },
                

    local sys = {
        type: 'ResponseSys',
        data: sig_binning {
            // overall_short_padding should take into account this offset
            // "start".  currently all "start" should be the same cause we only
            // have an overall time offset compensated in reframer. These values
            // correspond to files.fields[0, 1, 2] e.g. normal, shorted U, and
            // shorted Y
            start: [-10*wc.us, -10*wc.us, -10*wc.us], 
            magnitude: [1.0, 1.0, 1.0],
            time_smear: [0.0*wc.us, 0.0*wc.us, 0.0*wc.us],
        }
    },

    // Make some depos.  Note, this does NOT have a bagger.
    track_depos : function(tracklist = [{
        time: 0,
        charge: -5000,         
        ray:  {
            tail: wc.point(-4.0, 0.0, 0.0, wc.m),
            head: wc.point(+4.0, 6.1, 7.0, wc.m),
        }}])
        low.gen.track_depos(tracklist),

    // Simple signal simulation with a single FR.
    signal_single(anode, name) : pg.pipeline([
        local pirs = make_pirs(index=0);
        pg.pipeline([
            low.drifter(services.random,
                        low.util.driftsToXregions(params.geometry.volumes),
                        params.lar, name=name),
            transforms(anode, pirs)[params.ductor.transform],
            low.reframer(params, anode, name=name),
        ])]),
        
    // For signal_multiple().
    local make_branch = function(anode, index)
        local pirs = make_pirs(index);
        pg.pipeline([
            pg.pnode({
                // despite the name, just calls whatever its "drifter" is
                type: 'DepoSetDrifter',
                name: broken_detector[index].name,
                data: {
                    drifter: 'WireBoundedDepos:' + broken_detector[index].name,
                }
            }, nin=1, nout=1, uses=[pg.pnode({
                type: 'WireBoundedDepos',
                name: broken_detector[index].name,
                data: {
                    anode: wc.tn(anode),
                    regions: broken_detector[index].regions,
                    mode: broken_detector[index].mode,
                },
            }, nin=1, nout=1, uses = [anode])]),
            transforms(anode, pirs, index)[params.ductor.transform],
        ]),

    // More complex signal subgraph needed for "actual".
    signal_multiple(anode, name) :
        local ubsigtags = ['ubsig%d'%n for n in [0,1,2]];
        pg.pipeline([
            // Super position of multiple FRs
            pg.fan.pipe('DepoSetFanout',
                    [make_branch(anode, index) for index in [0,1,2]],
                    'FrameFanin', 'ubsigraph', ubsigtags),
        // A special reframer with a magic number time offset.
        pg.pnode({              
            type: 'Reframer',
            name: name,
            data: {
                anode: wc.tn(anode),
                local wtf_is_this = 81*wc.us,
                tbin: wtf_is_this/(params.ductor.binning.tick),
                nticks: params.ductor.binning.nticks,
                toffset: params.ductor.drift_dt - wtf_is_this,
                tags: ubsigtags,
                fill: 0.0,
            },
        }, nin=1, nout=1, uses=[anode])]),


    signal : self["signal_" + params.ductor.universe],


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

    
    noise(anode, name) :
        local model = {
            type: 'EmpiricalNoiseModel',
            name: name,
            data: params.noise.model {
                anode: wc.tn(anode),
                dft: wc.tn(services.dft),
                chanstat:wc.tn(csdb),
            },
            uses: [anode, csdb, services.dft]
        };
        pg.pnode({
            type: 'AddNoise',
            name: name,
            data: {
                rng: wc.tn(services.random),
                dft: wc.tn(services.dft),
                model: wc.tn(model),
                nsamples: params.noise.model.nsamples,
                replacement_percentage: params.noise.replacement_percentage,
            }}, nin=1, nout=1, uses=[services.random, services.dft, model]),



}
