
local reality = std.extVar('reality');
local sigoutform = std.extVar('signal_output_form');  // eg "sparse" or "dense"
local save_tradsp = true;

local wc = import 'wirecell.jsonnet';
local f = import "pgrapher/common/funcs.jsonnet";
local g = import 'pgraph.jsonnet';

local raw_input_label = std.extVar('raw_input_label');  // eg "daq"


local data_params = import 'params.jsonnet';
local simu_params = import 'simparams.jsonnet';
local params = if reality == 'data' then data_params else simu_params;


local tools_maker = import 'pgrapher/common/tools.jsonnet';
local tools = tools_maker(params);

local wcls_maker = import 'pgrapher/ui/wcls/nodes.jsonnet';
local wcls = wcls_maker(params, tools);

//local chndb_maker = import "pgrapher/experiment/pdsp/chndb.jsonnet";

local sp_maker = import 'pgrapher/experiment/pdhd/sp.jsonnet';

local use_resampler = (reality == 'data');

// Collect the WC/LS input converters for use below.  Make sure the
// "name" argument matches what is used in the FHiCL that loads this
// file.  In particular if there is no ":" in the inputer then name
// must be the emtpy string.
local wcls_input = {
  adc_digits: g.pnode({
    type: 'wclsRawFrameSource',
    name: '',
    data: {
      art_tag: raw_input_label,
      frame_tags: ['orig'],  // this is a WCT designator
      // nticks: params.daq.nticks,
      tick: 512*wc.ns,
    },
  }, nin=0, nout=1),

};

// Collect all the wc/ls output converters for use below.  Note the
// "name" MUST match what is used in theh "outputers" parameter in the
// FHiCL that loads this file.
local mega_anode = {
  type: 'MegaAnodePlane',
  name: 'meganodes',
  data: {
    anodes_tn: [wc.tn(anode) for anode in tools.anodes],
  },
};

local wcls_output = {
  // The noise filtered "ADC" values.  These are truncated for
  // art::Event but left as floats for the WCT SP.  Note, the tag
  // "raw" is somewhat historical as the output is not equivalent to
  // "raw data".
  nf_digits: g.pnode({
    type: 'wclsFrameSaver',
    name: 'nfsaver',
    data: {
      // anode: wc.tn(tools.anode),
      anode: wc.tn(mega_anode),
      digitize: true,  // true means save as RawDigit, else recob::Wire
      frame_tags: ['raw'],
      // nticks: params.daq.nticks,
      chanmaskmaps: ['bad'],
    },
  }, nin=1, nout=1, uses=[mega_anode]),


  // The output of signal processing.  Note, there are two signal
  // sets each created with its own filter.  The "gauss" one is best
  // for charge reconstruction, the "wiener" is best for S/N
  // separation.  Both are used in downstream WC code.
  sp_signals: g.pnode({
    type: 'wclsFrameSaver',
    name: 'spsaver',
    data: {
      // anode: wc.tn(tools.anode),
      anode: wc.tn(mega_anode),
      digitize: false,  // true means save as RawDigit, else recob::Wire
      frame_tags: ['gauss', 'wiener'],
      frame_scale: [0.001, 0.001],
      // nticks: params.daq.nticks,
      chanmaskmaps: [],
      nticks: -1,
    },
  }, nin=1, nout=1, uses=[mega_anode]),

  dnn_signals: g.pnode({
    type: 'wclsFrameSaver',
    name: 'dnnsaver',
    data: {
      anode: wc.tn(mega_anode),
      digitize: false,  // true means save as RawDigit, else recob::Wire
      frame_tags: ['dnnsp'],
      frame_scale: [0.001],
      nticks: -1,

    },
  }, nin=1, nout=1, uses=[mega_anode]),
};

// local perfect = import 'chndb-perfect.jsonnet';
local base = import 'chndb-base.jsonnet';
local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  // data: perfect(params, tools.anodes[n], tools.field, n) { dft:wc.tn(tools.dft) },
  data: base(params, tools.anodes[n], tools.field, n) { dft:wc.tn(tools.dft) },
  uses: [tools.anodes[n], tools.field, tools.dft],
} for n in std.range(0, std.length(tools.anodes) - 1)];

local nf_maker = import 'pgrapher/experiment/pdhd/nf.jsonnet';
local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n) for n in std.range(0, std.length(tools.anodes) - 1)];

local sp_override = { // assume all tages sets in base sp.jsonnet
    sparse: sigoutform == 'sparse',
    // wiener_tag: "",
    // gauss_tag: "",
    use_roi_refinement: true,
    use_roi_debug_mode: true,
    save_negtive_charge: false, // no negative charge in gauss
    tight_lf_tag: "",
    // loose_lf_tag: "",
    cleanup_roi_tag: "",
    break_roi_loop1_tag: "",
    break_roi_loop2_tag: "",
    shrink_roi_tag: "",
    extend_roi_tag: "",
    // decon_charge_tag: "",
    use_multi_plane_protection: true,
    do_not_mp_protect_traditional: true, // do_not_mp_protect_traditional to 
                                         // make a clear ref, defualt is false
    mp_tick_resolution: 10,
};
//local sp = sp_maker(params, tools, { sparse: sigoutform == 'sparse' });
local sp = sp_maker(params, tools, sp_override);
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];

local img = import 'pgrapher/experiment/pdhd/img.jsonnet';
local img_maker = img({use_dnn_img: true});
local img_pipes = [img_maker.per_anode(a) for a in tools.anodes];

//local util = import 'pgrapher/experiment/pdhd/funcs.jsonnet';
local chsel_pipes = [
  g.pnode({
    type: 'ChannelSelector',
    name: 'chsel%d' % n,
    data: {
      channels: std.range(2560 * n, 2560 * (n + 1) - 1),
      // tags: ['orig%d' % n], // traces tag
    },
  }, nin=1, nout=1)
  for n in std.range(0, std.length(tools.anodes) - 1)
];

local dnnroi = import 'pgrapher/experiment/pdhd/dnnroi.jsonnet';
local ts = {
    type: "TorchService",
    name: "dnnroi",
    data: {
        model: "ts-model/unet-cosmic390-newwc-depofluxsplat-pdhd.ts",
        device: "cpu", // "gpucpu",
        concurrency: 1,
    },
};

local resamplers_config = import 'pgrapher/common/resamplers.jsonnet';
local load_resamplers = resamplers_config(g, wc, tools);
local resamplers = load_resamplers.resamplers;

local magoutput = 'protodunehd-data-check.root';
local magnify = import 'pgrapher/experiment/pdhd/magnify-sinks.jsonnet';
local magio = magnify(tools, magoutput);

local use_magnify = std.extVar("use_magnify");
local nfsp_pipes = [
  g.pipeline(
    [chsel_pipes[n]] +
    (if use_resampler then [resamplers[n]] else []) +
    (if use_magnify =='true' then [
      magio.orig_pipe[n],
      nf_pipes[n],
      magio.raw_pipe[n],
      sp_pipes[n],
      magio.decon_pipe[n],
      dnnroi(tools.anodes[n], ts, output_scale=1.0, nticks=params.daq.nticks, nchunks=1),
      magio.dnnsp_pipe[n],
      // magio.threshold_pipe[n],
      // magio.debug_pipe[n], // use_roi_debug_mode=true in sp.jsonnet
      img_pipes[n],
    ]
    else [
      nf_pipes[n],
      sp_pipes[n],
      dnnroi(tools.anodes[n], ts, output_scale=1.0, nticks=params.daq.nticks, nchunks=1),
      img_pipes[n],
    ]),
    'nfsp_pipe_%d' % n)
    for n in std.range(0, std.length(tools.anodes) - 1)
];

// local fanpipe = util.fanpipe('FrameFanout', nfsp_pipes, 'FrameFanin', 'sn_mag_nf');
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
              ['threshold%d'%ind]:'threshold%d'%ind,
              // ['tight_lf%d'%ind]:'tight_lf%d'%ind,
              ['loose_lf%d'%ind]:'loose_lf%d'%ind,
            },

          }
          for ind in anode_ident
        ];
// local fanpipe = util.fanpipe('FrameFanout', nfsp_pipes, 'FrameFanin', 'nfsp', [], fanout_tag_rules, fanin_tag_rules);

local nanodes = std.length(tools.anodes);
local fanpipe = f.multifanout('FrameFanout', nfsp_pipes, [1,nanodes], [nanodes,1], 'sn_mag', fanin_tag_rules);

local retagger = g.pnode({
  type: 'Retagger',
  name: 'dnnout',
  data: {
    tag_rules: [{
      frame: {'.*': 'dnnretagger',},
      merge: {'dnnsp\\d': 'dnnsp',},
    }],
  },
}, nin=1, nout=1);

local sink = g.pnode({ type: 'DumpFrames' }, nin=1, nout=0);
// local graph = g.pipeline([wcls_input.adc_digits, fanpipe, retagger, wcls_output.dnn_signals, sink]);
local graph = g.pipeline([wcls_input.adc_digits, fanpipe], retagger);

// Build an incomplete subgraph ending to be spliced for saving out frames 
local ofanin = g.pnode({ 
      type: 'FrameFanin',
      name:"outfanin",
      data:{
          multiplicity: std.length(tools.anodes),
          tag_rules: [
            {
              frame: {'.*': 'outfanin',},
              trace: {
                ['gauss%d' % n]: ['gauss%d' % n],
                ['wiener%d' % n]: ['wiener%d' % n],
                // ['threshold%d' % n]: ['threshold%d' % n],
              },
            }
            for n in std.range(0, std.length(tools.anodes) - 1)
          ],
      } 
      }, nin=std.length(tools.anodes), nout=1);
local osink = g.pnode({ type: 'DumpFrames', name:"outsink", data:{} }, nin=1, nout=0);
// local outsgr = g.intern(innodes=[ofanin,], centernodes = [osink],
//                       edges=[ g.edge(ofanin, osink) ], name="outsgr");
local outretagger = g.pnode({
  type: 'Retagger',
  name: 'spout',
  data: {
    tag_rules: [{
      frame: {'.*': 'spretagger',},
      merge: {
        'gauss\\d': 'gauss',
        'wiener\\d': 'wiener',
        // 'threshold\\d': 'threshold',
      },
    }],
  },
}, nin=1, nout=1);
local outgr = g.pipeline([ofanin, outretagger, wcls_output.sp_signals, osink]);

local edge_selector(e) = std.startsWith(e.tail.node, "OmnibusSigProc:");
local fanout_factory(n,e) = { type:'FrameFanout', name:"splice%d"%n, data:{multiplicity: 2} }; // "2-wire" splice

local spliced_graph =
  if save_tradsp then
    g.splice(graph, outgr, edge_selector, fanout_factory)
  else
    graph;

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(spliced_graph),
  },
};

// Finally, the configuration sequence
g.uses(spliced_graph) + [app]
