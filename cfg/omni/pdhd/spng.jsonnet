// This file provides support for an SPNG variant of PDHD.
//
// It is cribbed from code in spng/test/wct-framesource.jsonnet
//
// For now, SPNG configuration needs a major conceptual cleanup and so we must
// write this file specific to PDHD.

local pg = import "pgraph.jsonnet";
local wc = import "wirecell.jsonnet";

local filters = import "spng-filters.jsonnet";

// Used by TorchFRERSpectrum but we really should remove it.
local should_not_need(plane) = {
    default_nchans : [800, 800, 480, 480][plane],
    default_nticks: 6000,
    default_period: 500*wc.ns,
};
    

// Little helper
local apname(anode, plane) = anode.name + std.toString(plane);

{
    // Return the FR+ER response component bundling the individuals.  The fr and
    // er are as provided by omni.responses(), the adc is omni.adc.    
    frer(anode, plane, fr, er, adc) :: {
        type: 'TorchFRERSpectrum',
        name: apname(anode, plane),
        data: {
            field_response: wc.tn(fr),
            elec_respnse: wc.tn(er),
            fr_plane_id: plane,
            // Should rename to like "ADC_count_per_voltage".  Units are NOT ADC/mV.
            ADC_mV: 1.0/adc.lsb_voltage,
            // should rename to "gain"
            inter_gain: adc.gain,
            // should remove
            anode_num: anode.data.ident,
        } + should_not_need(plane),
        uses: [fr, er]
    }

    // Return a per-plane SPNG plan deconvolution node
    decon(anode, plane, frer) :: pg.pnode({
        local wf = if plane > 1
                   then spng_filters.torch_wire_filters[1]
                   else spng_filters.torch_wire_filters[0],
        type: 'SPNGDecon',
        name: apname(anode, plane),
        data: {
            frer_spectrum: wc.tn(frer),
            wire_filter: wc.tn(wf) ,
            coarse_time_offset: 1000.0,
            pad_wire_domain: (plane > 1), #Non-periodic planes get padded
        },
    } nin=1, nout=1, uses=[wf, frer]),

    // A 4-way fanout treating wrapped U and V as one group and W as two.
    fanout_uvww(anode) :: pg.pnode({
        type: 'FrameToTorchSetFanout',
        name: anode.name,
        data: {
            anode: wc.tn(anode),

            /// fixme: why do we give this number?
            expected_nticks: 6000,

            output_groups: [
                [wc.WirePlaneId(wc.Ulayer, 0, anode.data.ident),
                 wc.WirePlaneId(wc.Ulayer, 1, anode.data.ident)],
                
                [wc.WirePlaneId(wc.Vlayer, 0, anode.data.ident),
                 wc.WirePlaneId(wc.Vlayer, 1, anode.data.ident)],
                
                [wc.WirePlaneId(wc.Wlayer, 0, anode.data.ident)],
                
                [wc.WirePlaneId(wc.Wlayer, 1, anode.data.ident)],
            ],
        }
    }, nin=1, nout=4, uses=[anode]),

    // Return a "sigproc subgraph" for the context of one anode.
    // For now, this gives a 1-in/4-out subgraph
    sigproc(anode, resps, adc) ::
        local fr = resps.fr[0];
        local er = resps.er[0];
        local fout = fanout_uvww(anode);
        // local fin = fanin_uvww(anode);
        local decons = [$.decon(anode, plane, $.frer(anode, plane, fr, er, adc)) for plane in [0,1,2,3]];
        // Not finished, 
        pg.intern(
            innodes=[fout],
            outnodes=decons,
            centernodes=[],
            edges = [
                pg.edge(fout, spng_decons[0], 0),
                pg.edge(fout, spng_decons[1], 1),
                pg.edge(fout, spng_decons[2], 2),
                pg.edge(fout, spng_decons[3], 3),
            ])
        
