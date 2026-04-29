// This provides some noise filtering related pnodes,

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';
local sp_filters = import 'pgrapher/experiment/pdhd/sp-filters.jsonnet';

local default_dft = { type: 'FftwDFT' };

function(params, anode, chndbobj, n, name='', dft=default_dft) {
    local single = {
        type: 'PDHDOneChannelNoise',
        name: name,
        uses: [dft, chndbobj, anode],
        data: {
            noisedb: wc.tn(chndbobj),
            anode: wc.tn(anode),
            dft: wc.tn(dft),
            // adaptive_baseline left at C++ default (false): PDHD cold
            // electronics is DC-coupled, so the IS_RC partial-RC gate that
            // fronts the adaptive baseline (see Microboone.cxx:963-1047) has
            // no physical meaning here. Side effect: the lf_noisy mask
            // emitted under is_partial in ProtoduneHD.cxx is no longer
            // produced; any bad-channel info will be supplied separately.
        },
    },
    local grouped = {
        type: 'PDHDCoherentNoiseSub',
        name: name,
        uses: [dft, chndbobj, anode] + sp_filters,
        data: {
            noisedb: wc.tn(chndbobj),
            anode: wc.tn(anode),
            dft: wc.tn(dft),
            rms_threshold: 0.0,
            time_filters:
              if anode.data.ident == 0
              then ['Wiener_tight_U_APA1', 'Wiener_tight_V_APA1', 'Wiener_tight_W_APA1']
              else ['Wiener_tight_U', 'Wiener_tight_V', 'Wiener_tight_W'],
            lf_tighter_filter: 'ROI_tighter_lf',
            lf_loose_filter: 'ROI_loose_lf',
        },
    },

    local obnf = g.pnode({
        type: 'OmnibusNoiseFilter',
        name: name,
        data: {

            // Nonzero forces the number of ticks in the waveform
            nticks: 0,

            // channel bin ranges are ignored
            // only when the channelmask is merged to `bad`
            // maskmap: {sticky: "bad", ledge: "bad", noisy: "bad"},
            // maskmap: {noisy:"bad", lf_noisy: "bad"},
            channel_filters: [
                wc.tn(single),
            ],
            grouped_filters: [
                wc.tn(grouped),
            ],
            channel_status_filters: [
            ],
            noisedb: wc.tn(chndbobj),
            // intraces: 'orig%d' % n,  // frame tag get all traces
            // intraces: 'orig',        // use when orig frames have tag 'orig'
            intraces: '',  // '' means use all traces (wildcard); HD orig frames use '*' tag
            outtraces: 'raw%d' % n,
        },
    }, uses=[chndbobj, anode, single, grouped], nin=1, nout=1),


    pipe: g.pipeline([obnf], name=name),
}.pipe
