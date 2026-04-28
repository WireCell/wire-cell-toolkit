local reality = std.extVar('reality');
local charge_input_label = std.extVar('charge_input_label');  // "gauss", "dnnsp"
local wiener_input_label = std.extVar('wiener_input_label');  // "wiener"

local wc = import 'wirecell.jsonnet';
local f = import "pgrapher/common/funcs.jsonnet";
// local util = import 'pgrapher/experiment/pdsp/funcs.jsonnet';
local g = import 'pgraph.jsonnet';

local data_params = import 'pgrapher/experiment/pdsp/params.jsonnet';
local simu_params = import 'pgrapher/experiment/pdsp/simparams.jsonnet';
local params = if reality == 'data' then data_params else simu_params;

local tools_maker = import 'pgrapher/common/tools.jsonnet';
local tools = tools_maker(params);

local wcls_maker = import 'pgrapher/ui/wcls/nodes.jsonnet';
local wcls = wcls_maker(params, tools);
// Collect the WC/LS input converters for use below.  Make sure the
// "name" argument matches what is used in the FHiCL that loads this
// file.  In particular if there is no ":" in the inputer then name
// must be the emtpy string.
local wcls_input = {
  charge_input: g.pnode({
    type: 'wclsCookedFrameSource',
    name: 'charge',
    data: {
      recobwire_tags: [charge_input_label], summary_tags: [''], trace_tags: ['gauss'],
      frame_tags: ['gauss'],  // this is a WCT designator
      // nticks: params.daq.nticks,
    },
  }, nin=0, nout=1),

  wiener_input: g.pnode({
    type: 'wclsCookedFrameSource',
    name: 'wiener',
    data: {
      recobwire_tags: [wiener_input_label], summary_tags: ['wclsdatasp:thresholdsummary'], trace_tags: ['wiener'],
      frame_tags: ['wiener'],  // this is a WCT designator
      // nticks: params.daq.nticks,
    },
  }, nin=0, nout=1),
};

local mega_anode = {
  type: 'MegaAnodePlane',
  name: 'meganodes',
  data: {
    anodes_tn: [wc.tn(anode) for anode in tools.anodes],
  },
};

local img = import 'pgrapher/experiment/pdsp/img.jsonnet';
local img_maker = img();
local img_pipes = [img_maker.per_anode(a) for a in tools.anodes];

local chsel_pipes = [
  g.pnode({
    type: 'ChannelSelector',
    name: 'chsel%d' % n,
    data: {
      channels: std.range(2560 * n, 2560 * (n + 1) - 1),
      tags: ["gauss","wiener"],
      // tag_rules: [{
      //   frame: {'.*': 'sigproc',},
      //   trace: {
      //     'gauss': 'gauss%d' %n,
      //     'wiener': 'wiener%d' %n,
      //   },
      // }],
    },
  }, nin=1, nout=1)
  for n in std.range(0, std.length(tools.anodes) - 1)
];

local imgfanpipes = [
  g.pipeline([
               chsel_pipes[n],
               img_pipes[n],
             ],
             'imgfan_pipe_%d' % n)
  for n in std.range(0, std.length(tools.anodes) - 1)
];

// local fanpipe = util.fanpipe('FrameFanout', imgfanpipes, 'FrameFanin', 'sn_mag_nf');
local fanout_tag_rules = [ 
          {
            frame: {
              '.*': 'orig%d' % tools.anodes[n].data.ident,
            },
            trace: {
              // fake doing Nmult SP pipelines
              //orig: ['wiener', 'gauss'],
              //'.*': 'orig',
            },
          }
          for n in std.range(0, std.length(tools.anodes) - 1)
        ];

local anode_ident = [tools.anodes[n].data.ident for n in std.range(0, std.length(tools.anodes) - 1)];
local fanin_tag_rules = [
          {
            frame: {
              //['number%d' % n]: ['output%d' % n, 'output'],
              '.*': 'framefanin',
            },
            trace: {
              ['gauss%d'%ind]:'gauss%d'%ind,
              ['wiener%d'%ind]:'wiener%d'%ind,
              ['thresholdsummary%d'%ind]:'thresholdsummary%d'%ind,
            },

          }
          for ind in anode_ident
        ];

// local fanpipe = util.fanpipe('FrameFanout', imgfanpipes, 'FrameFanin', 'nfsp', [], fanout_tag_rules, fanin_tag_rules);
local nanodes = std.length(tools.anodes);
local fanpipe = f.multifanout('FrameFanout', imgfanpipes, [1,nanodes], [nanodes,1], 'sn_mag', fanin_tag_rules);

local retagger = g.pnode({
  type: 'Retagger',
  data: {
    tag_rules: [{
      frame: {
        '.*': 'retagger',
      },
      merge: {
        'gauss\\d': 'gauss',
        'wiener\\d': 'wiener',
      },
    }],
  },
}, nin=1, nout=1);

// duo-head source input: charge input ("gauss" or "dnnsp") & wiener input
local sources = [wcls_input.charge_input, wcls_input.wiener_input];

local joiner = g.pnode({ 
      type: 'FrameFanin',
      name:"joins",
      data:{
          multiplicity: std.length(sources),
          tags: ["gauss", "wiener"],
          tag_rules: [
            {
              frame: {'.*': 'joins',},
              trace: {
                "gauss" : "gauss",
                "wiener" : "wiener",
              },
            }
          ],
      } 
      }, nin=std.length(sources), nout=1);

local joined_source = g.join_sources(joiner, sources, n=2);
local graph = g.pipeline([joined_source, fanpipe]);
// local sink = g.pnode({ type: 'DumpFrames' }, nin=1, nout=0);

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
};

// Finally, the configuration sequence
g.uses(graph) + [app]
