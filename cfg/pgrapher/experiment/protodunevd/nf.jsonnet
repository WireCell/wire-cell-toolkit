// This provides some noise filtering related pnodes,

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

function(params, anode, chndbobj, n, name='')
  {
    local single = {
      type: 'PDVDOneChannelNoise',
      name: name,
      data: {
        noisedb: wc.tn(chndbobj),
        anode: wc.tn(anode),
      },
    },
    local grouped = {
      type: 'PDVDCoherentNoiseSub',
      name: name,
      data: {
        noisedb: wc.tn(chndbobj),
        anode: wc.tn(anode),
        rms_threshold: 0.0,
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
        ]else [],       
        grouped_filters: [
          // wc.tn(grouped),
        ],
        channel_status_filters: [
        ],
        noisedb: wc.tn(chndbobj),
        // intraces: 'orig%d' % n,  // frame tag get all traces
        intraces: 'orig',  // frame tag get all traces
        outtraces: 'raw%d' % anode.data.ident,
      },
    }, uses=[chndbobj, anode, single, grouped,shieldcoupling_grouped], nin=1, nout=1),


    pipe: g.pipeline([obnf], name=name),
  }.pipe
