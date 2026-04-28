// This provides some noise filtering related pnodes,

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

function(params, anode, chndbobj, n, name='')
  {
    local status = {
      type: std.trace("Warning MB in DUNE?", 'mbOneChannelStatus'),
      name: name,
      data: {
        Threshold: 3.5,
        Window: 5,
        Nbins: 250,
        Cut: 14,
        anode: wc.tn(anode),
      },
    },
    local single = {
      type: std.trace("Warning MB in DUNE?", 'mbOneChannelNoise'),
      name: name,
      data: {
        noisedb: wc.tn(chndbobj),
        anode: wc.tn(anode),
      },
    },
    local grouped = {
      type: std.trace("Warning MB in DUNE?", 'mbCoherentNoiseSub'),
      name: name,
      data: {
        noisedb: wc.tn(chndbobj),
        anode: wc.tn(anode),
      },
    },

    local obnf = g.pnode({
      type: 'OmnibusNoiseFilter',
      name: name,
      data: {

        // This is the number of bins in various filters
        nsamples: params.nf.nsamples,

        //maskmap: { chirp: "bad", noisy: "bad" },
        channel_filters: [
          wc.tn(single)
        ],
        grouped_filters: [
          // wc.tn(grouped),
        ],
        channel_status_filters: [
          // wc.tn(status),
        ],
        noisedb: wc.tn(chndbobj),
        intraces: 'orig%d' % n,  // frame tag get all traces
        outtraces: 'raw%d' % n,
      },
    // }, uses=[chndbobj, anode, single, grouped, status], nin=1, nout=1),
    }, uses=[chndbobj, anode, single], nin=1, nout=1),


    pipe: g.pipeline([obnf], name=name),
  }.pipe
