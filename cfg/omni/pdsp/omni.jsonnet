// The omni object(s) for ProtoDUNE Single-Phase (PDSP) detector.

local wc = import 'wirecell.jsonnet';

local params = import 'simparams.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local tools = tools_maker(params);

local omniapi = import "../../omni.jsonnet";

local detectors = import "detectors.jsonnet";
local det = detectors.pdsp;

local nf = import "nf.jsonnet";
local sp = import "sp.jsonnet";


// Assume we may have variants, so make an initial default omni object as a
// local to which we may override.
local omni = omniapi + nf + sp + {
    name: "pdsp",
    detname: "pdsp",
    
    // The detector parameters are from wcls-sim-drift-simchannel.fcl and
    // wcls-sim-drift-simchannel.jsonnet under pgrapher's pdsp.

    lar: super.lar + {
        lifetime : 10.4*wc.ms,
        drift_speed : 1.565*wc.mm/wc.us,
        DL: 4.0 * wc.cm2 / wc.s,
        DT: 8.8 * wc.cm2 / wc.s,
    },

    adc: super.adc + {
        // per tdr, chapter 2
        // induction plane: 2350 ADC, collection plane: 900 ADC
        baselines: [1003.4*wc.millivolt,1003.4*wc.millivolt,507.7*wc.millivolt],

        // check this.  The tdr says, "The ADC ASIC has an input
        // buffer with offset compensation to match the output of the
        // FE ASIC.  The input buffer first samples the input signal
        // (with a range of 0.2 V to 1.6 V)..."
        fullscale: [0.2*wc.volt, 1.6*wc.volt],        
    },

    // from original simparams.jsonnet.
    local apa_cpa = 3.63075*wc.m,
    local cpa_thick = 3.175*wc.mm,
    local apa_w2w = 85.725*wc.mm,
    local plane_gap = 4.76*wc.mm,
    local apa_g2g = 114.3*wc.mm,

    local apa_plane = 0.5*apa_g2g - plane_gap,
    local res_plane = 0.5*apa_w2w + $.detdata.response_plane,
    local cpa_plane = apa_cpa - 0.5*cpa_thick,

    // This follows original params.jsonnet
    volumes:  [ {
        local sign = 2*(n%2)-1,
        local centerline = sign*apa_cpa,
        wires: n,       // anode number
        name: "apa%d"%n,
        faces:
            // top, front face is against cryo wall
            if sign > 0
            then [
                null,           // insensitive
                {
                    anode: centerline - apa_plane,
                    response: centerline - res_plane,
                    cathode: centerline - cpa_plane, 
                }
            ]
            // bottom, back face is against cryo wall
            else [
                {
                    anode: centerline + apa_plane,
                    response: centerline + res_plane,
                    cathode: centerline + cpa_plane, 
                },
                null            // insensitive
            ],
    } for n in std.range(0,5)],

    binning: super.binning {
        local tick = 500*wc.ns,
        local nticks = 6000,
        default_time_binning: {
            spacing: tick,
            number: nticks,
            start: 0.0*wc.ns
        },
        // The "absolute" time (ie, in G4 time) that the lower edge of
        // of final readout tick #0 should correspond to.  This is a
        // "fixed" notion.
        local tick0_time = -250*wc.us,
        
        // Open the ductor's gate a bit early.
        local response_time_offset = $.detdata.response_plane / $.lar.drift_speed,
        local response_nticks = wc.roundToInt(response_time_offset / tick),

        sim: {
            spacing: tick,
            number: nticks + response_nticks,
            start: tick0_time - response_time_offset,
        },
        // fixme: check sp, splat
    },


    // pdsp needs unique er and rc.  fr is copy-pasted from parent.
    responses(anode, kind, binning=$.binning.sim) :: {

        // Note, for now we do not distinguish based on "kind".

        fr: [{
            type: 'FieldResponse',
            name: $.detname,
            data: { filename: $.detdata.field } }],

        er: [{
            type: 'ColdElecResponse',
            name: $.detname,
            data: {
                start: binning.start,
                tick: binning.spacing,
                nticks: binning.number,
                shaping : 2.2*wc.us,
                gain: 14.0*wc.mV/wc.fC, 
                postgain: 1.1365,
            }}],

        local rcr = {
            type: 'RCResponse',
            name: $.detname,
            data: {
                width: 1.1*wc.ms,
            }
        },
        // "RCRC" response
        rc: [rcr, rcr]
    },



    // local img = import "img.jsonnet",
    // img :: img,
        
    

};

omni

