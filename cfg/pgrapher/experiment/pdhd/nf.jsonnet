// This provides some noise filtering related pnodes,

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';
local sp_filters = import 'pgrapher/experiment/pdhd/sp-filters.jsonnet';

local default_dft = { type: 'FftwDFT' };

function(params, anode, chndbobj, n, name='', dft=default_dft,
         debug_dump_path='', debug_dump_groups=[],
         maskmap=null) {
    // pr_ml doc 32: params.apa0_nominal (absent = false = production) treats APA 0 like APAs 1-3 (the training-sample
    // simulation with a nominal APA 0 field response, pdhd/d31 apa0_nominal); absent => compiled config unchanged.
    local apa0_special = anode.data.ident == 0 && !(std.objectHas(params, 'apa0_nominal') && params.apa0_nominal),
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
              if apa0_special
              then ['Wiener_tight_U_APA1', 'Wiener_tight_V_APA1', 'Wiener_tight_W_APA1']
              else ['Wiener_tight_U', 'Wiener_tight_V', 'Wiener_tight_W'],
            lf_tighter_filter: 'ROI_tighter_lf',
            lf_loose_filter: 'ROI_loose_lf',
            debug_dump_path: debug_dump_path,
            debug_dump_groups: debug_dump_groups,
        },
    },
    local fembfilt = {
         type: 'PDHDFEMBNoiseSub',
         name: name,
         uses: [anode],
         data: {
             anode: wc.tn(anode),
             width: 50.0,
             pad_nticks: 20,
             nsigma: 3.5,
         },
     },

    local obnf = g.pnode({
        type: 'OmnibusNoiseFilter',
        name: name,
        data: {

            // Nonzero forces the number of ticks in the waveform
            nticks: 0,

            // maskmap routes NF tag names to output tag names.
            // null → use OmnibusNoiseFilter C++ defaults (chirp+noisy → bad).
            // Pass an explicit object to rename each tag independently, e.g.:
            //   {chirp:"chirp", noisy:"noisy", femb_noise:"femb_noise"}
            // or merge all to "bad":
            maskmap: { noisy:"bad", femb_noise:"bad"},
            // [if maskmap != null then 'maskmap']: maskmap,
            channel_filters: [
                wc.tn(single),
            ],
            grouped_filters: [
            //     wc.tn(grouped),
            ],
            multigroup_chanfilters: [
             {
               channelgroups: chndbobj.data.femb_negpulse_groups,
               filters: [wc.tn(fembfilt)],
             },
            {
               channelgroups: chndbobj.data.groups,
               filters: [wc.tn(grouped)],
             },
            ],

            channel_status_filters: [
            ],
            noisedb: wc.tn(chndbobj),
            // intraces: 'orig%d' % n,  // frame tag get all traces
            // intraces: 'orig',        // use when orig frames have tag 'orig'
            intraces: '',  // '' means use all traces (wildcard); HD orig frames use '*' tag
            outtraces: 'raw%d' % n,
        },
    }, uses=[chndbobj, anode, single, grouped, fembfilt], nin=1, nout=1),


    pipe: g.pipeline([obnf], name=name),
}.pipe
