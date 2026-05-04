// This provides some noise filtering related pnodes,

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';
local sp_filters = import 'pgrapher/experiment/protodunevd/sp-filters.jsonnet';

function(params, anode, chndbobj, n, name='',
         debug_dump_path='', debug_dump_groups=[])
  {
    local single = {
      type: 'PDVDOneChannelNoise',
      name: name,
      data: {
        noisedb: wc.tn(chndbobj),
        anode: wc.tn(anode),
        // adaptive_baseline left at C++ default (false): PDVD hardware has
        // no AC coupling, so the IS_RC partial-RC gate that fronts the
        // adaptive baseline (see Microboone.cxx:963-1047) has no physical
        // meaning here.
      },
    },
    // Top (_t) vs bottom (_b) anode filter suffix.  Bottom = ident 0..3,
    // top = ident 4..7.  See sp-filters.jsonnet for the registered names.
    local sfx = if anode.data.ident < 4 then '_b' else '_t',
    local grouped = {
      type: 'PDVDCoherentNoiseSub',
      name: name,
      uses: [chndbobj, anode] + sp_filters,
      data: {
        noisedb: wc.tn(chndbobj),
        anode: wc.tn(anode),
        rms_threshold: 0.0,
        time_filters: ['Wiener_tight_U' + sfx,
                       'Wiener_tight_V' + sfx,
                       'Wiener_tight_W' + sfx],
        lf_tighter_filter: 'ROI_tighter_lf' + sfx,
        lf_loose_filter:   'ROI_loose_lf'   + sfx,
        debug_dump_path: debug_dump_path,
        debug_dump_groups: debug_dump_groups,
      },
    },
    local shieldcoupling_grouped = {
         type: 'PDVDShieldCouplingSub',
         name: name,
         uses: [anode],
         data: {
             anode: wc.tn(anode),
             noisedb: wc.tn(chndbobj),
             strip_length: params.files.strip_length,
             rms_threshold: 0.0,
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
        maskmap: {sticky: "bad", ledge: "bad", noisy: "bad"},
        channel_filters: [
          wc.tn(single),
        ],
             multigroup_chanfilters: [

            {
               channelgroups: chndbobj.data.groups,
               filters: [wc.tn(grouped)],
             },

            ]
            // only apply to top
            +if anode.data.ident > 3 then[
              {
                channelgroups: chndbobj.data.top_u_groups,
                filters: [wc.tn(shieldcoupling_grouped)],
              },
            ]else []
            ,       
        grouped_filters: [
          // wc.tn(grouped),
        ],
        channel_status_filters: [
        ],
        noisedb: wc.tn(chndbobj),
        intraces: 'orig',  // frame tag get all traces
        outtraces: 'raw%d' % anode.data.ident,
      },
    }, uses=[chndbobj, anode, single, grouped,shieldcoupling_grouped], nin=1, nout=1),


    pipe: g.pipeline([obnf], name=name),
  }.pipe
