// Part 2 of a two-part split of wcls-nf-sp-img.jsonnet.
//
// This file runs only the imaging step, reading the SP output frames
// (gauss/wiener) from a tar file produced by FrameFileSink in the
// wcls-nf-sp.jsonnet step (via the frame_tap below).
//
// The tar format (FrameFileSink / FrameFileSource) preserves real channel
// IDs and trace summaries (RMS), which are required by pre_proc
// (ChargeErrorFrameEstimator) and MaskSlices.
//
// Run standalone with:
//   wire-cell -l stdout -L debug -c wct-img.jsonnet
//
// Or with an explicit input file:
//   wire-cell -l stdout -L debug \
//     --tla-str input="protodune-sp-frames.tar.bz2" \
//     -c wct-img.jsonnet
//
// To produce the input file, add a FrameFileSink tap to wcls-nf-sp.jsonnet
// after the sp_pipes step, e.g. via the frame_tap helper in
// wct-sim-sigproc-fans.jsonnet:
//
//   local frame_tap = function(name, outname, tags, digitize) { ... };
//   // insert frame_tap(...) after sp_pipes[n] in nfsp_pipes

local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';

// Use simulation params for geometry/tools
local params = import 'pgrapher/experiment/protodunevd/simparams.jsonnet';

local tools_maker = import 'pgrapher/common/tools.jsonnet';
local tools_all = tools_maker(params);
local tools = tools_all
              {
  anodes: [
    // tools_all.anodes[4]
    tools_all.anodes[5],
    // tools_all.anodes[6],
    // tools_all.anodes[7]
  ],
};

// Top-level function to allow overriding input file via TLA
function(
  input='protodune-sp-frames.tar.bz2'
)

  // ---- FrameFileSource: read SP frames from tar file -----------------------
  // FrameFileSource preserves real channel IDs and trace summaries (RMS),
  // unlike MagnifySource which reassigns sequential channel numbers from 0.
  // The tags list selects which tagged trace sets to load as tagged traces.

  local anode_ident = tools.anodes[0].data.ident;  // e.g. 3

  local frame_source = g.pnode({
    type: 'FrameFileSource',
    name: 'frame_source',
    data: {
      inname: input,
      tags: [
        'gauss%d' % anode_ident,
        'wiener%d' % anode_ident,
      ],
      // Optionally apply frame-level tags to the produced frame:
      // frame_tags: [],
    },
  }, nin=0, nout=1);

  // ---- Imaging pipeline per anode ------------------------------------------
  local img = import 'pgrapher/experiment/protodunevd/img.jsonnet';
  local img_maker = img();
  local img_pipes = [img_maker.per_anode(a) for a in tools.anodes];

  // ---- Full pipeline: source -> imaging ------------------------------------
  local nanodes = std.length(tools.anodes);

  local graph =
    if nanodes == 1
    then g.pipeline([frame_source, img_pipes[0]], 'img_graph')
    else error 'wct-img.jsonnet: only single-anode mode supported with FrameFileSource';

  // ---- App and cmdline configuration ---------------------------------------

  local app = {
    type: 'Pgrapher',
    data: {
      edges: g.edges(graph),
    },
  };

  local cmdline = {
    type: 'wire-cell',
    data: {
      plugins: [
        'WireCellGen',
        'WireCellPgraph',
        'WireCellSio',  // for FrameFileSource
        'WireCellSigProc',
        'WireCellImg',
        'WireCellClus',
      ],
      apps: ['Pgrapher'],
    },
  };

  [cmdline] + g.uses(graph) + [app]
