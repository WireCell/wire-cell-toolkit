// This produces a ChannelNoiseDb with actual bad channels in all
// their gory detail.  See also chndb-perfect.jsonnet.

local low = import "../../../low.jsonnet";

local ci = import "channel-info.jsonnet";

// strictly, these span more than any given anode.
local bc = import "actual-bad-channels.jsonnet";

local groups = import "channel-groups.jsonnet";

function(anode, binning, fr, dft) {
    type: 'OmniChannelNoiseDB',
    name: low.util.idents(anode),
    uses: [anode, fr, dft],
    data: {
        anode: low.wc.tn(anode),
        field_response: low.wc.tn(fr),
        tick: binning.tick,
        nsamples: binning.nticks,
        groups: groups(anode.data.ident),
        bad: bc,
        channel_info: ci(anode.data.ident, binning.nticks)
    }
}

    
    
