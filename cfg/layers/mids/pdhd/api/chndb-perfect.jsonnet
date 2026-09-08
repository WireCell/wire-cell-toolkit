// This produces a ChannelNoiseDb for a "perfect" PDSP detector.  See
// also chndb-actual.jsonnet.

local low = import "../../../low.jsonnet";

local groups = import "channel-groups.jsonnet";

function(anode, binning, fr, dft) {
    type: 'OmniChannelNoiseDb',
    name: low.util.idents(anode),
    uses: [anode, fr, dft],
    data: {
        anode: low.wc.tn(anode),
        field_response: low.wc.tn(fr),
        tick: binning.tick,
        nsamples: binning.nticks,
        groups: groups(anode.data.ident, binning.nticks),

    }
}

    
    
