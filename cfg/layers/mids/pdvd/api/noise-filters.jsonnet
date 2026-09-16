// All possible noise filters for PDSP.  If additional variants are
// needed for any given type, make a fully new entry in this object.
// Do not play games with "if".  Then, in variants/*.jsonnet add
// whichever key names found here are needed to the corresponding
// nf.filters.<category> arrays.
//
// NOTE: the key names here are used by different variants.  The key
// *means* its value.  If you tweak a value in-place you are changing
// many variants.  

local low = import "../../../low.jsonnet";
local tn = low.wc.tn;

function(services, params, anode, chndb)
{
    local ident = low.util.idents(anode),

    single: {
        type: 'pdOneChannelNoise',
        name: ident,
        uses: [services.dft, chndb, anode],
        data: {
            noisedb: tn(chndb),
            anode: tn(anode),
            dft: tn(services.dft),
            resmp: [
                {channels: std.range(2128, 2175), sample_from: 5996},
                {channels: std.range(1520, 1559), sample_from: 5996},
                {channels: std.range( 440,  479), sample_from: 5996},
            ],
        },
    },
    grouped: {
        type: 'mbCoherentNoiseSub',
        name: ident,
        uses: [services.dft, chndb, anode],
        data: {
            noisedb: tn(chndb),
            anode: tn(anode),
            dft: tn(services.dft),
            rms_threshold: 0.0,
        },
    },
    sticky: {
        type: 'pdStickyCodeMitig',
        name: ident,
        uses: [services.dft, chndb, anode],
        data: {
            extra_stky: [
                {channels: std.range(anode.data.ident * 2560, (anode.data.ident + 1) * 2560 - 1), bits: [0,1,63]},
                {channels: [4], bits: [6]  },
                {channels: [159], bits: [6]  },
                {channels: [164], bits: [36] },
                {channels: [168], bits: [7]  },
                {channels: [323], bits: [24] },
                {channels: [451], bits: [25] },
            ],
            noisedb: tn(chndb),
            anode: tn(anode),
            dft: tn(services.dft),
            stky_sig_like_val: 15.0,
            stky_sig_like_rms: 2.0,
            stky_max_len: 10,
        },
    },
    gaincalib: {
        type: 'pdRelGainCalib',
        name: ident,
        uses: [chndb, anode],
        data: {
            noisedb: tn(chndb),
            anode: tn(anode),
            rel_gain: import "rel-gain.jsonnet",
        },
    }
}

