// Part 2 of a two-part split of wcls-nf-sp-img.jsonnet.
//
// Processes all anodes (or a subset) in one wire-cell run.
// Reads per-anode SP frames from files produced by wcls-nf-sp-out.jsonnet:
//   protodune-sp-frames-anode{N}.tar.bz2
//
// Each anode gets an independent FrameFileSource -> imaging pipeline.
// Pgrapher executes all independent source->sink chains.
//
// Run standalone (all anodes):
//   wire-cell -l stdout -L debug -c wct-img-all.jsonnet
//
// With a custom file prefix:
//   wire-cell -l stdout -L debug \
//     --tla-str input_prefix="protodune-sp-frames" \
//     -c wct-img-all.jsonnet
//
// To select a subset of anodes by index into tools_all.anodes:
//   wire-cell -l stdout -L debug \
//     --tla-code anode_indices='[4,5]' \
//     -c wct-img-all.jsonnet

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

local params = import 'pgrapher/experiment/protodunevd/simparams.jsonnet';

local tools_maker = import 'pgrapher/common/tools.jsonnet';
local tools_all = tools_maker(params);

// Top-level function: parameters overridable via --tla-str / --tla-code
function(
  // Prefix for per-anode input files: "{input_prefix}-anode{N}.tar.bz2"
  input_prefix = 'protodune-sp-frames',
  // Indices into tools_all.anodes to process; default = all
  anode_indices = std.range(0, std.length(tools_all.anodes) - 1)
)

  local anodes = [tools_all.anodes[i] for i in anode_indices];

  local img = import 'pgrapher/experiment/protodunevd/img.jsonnet';
  local img_maker = img();

  // Build one FrameFileSource + imaging pipeline per anode
  local per_anode_graph(anode) =
    local aid = anode.data.ident;
    local src = g.pnode({
      type: 'FrameFileSource',
      name: 'frame_source_anode%d' % aid,
      data: {
        inname: '%s-anode%d.tar.bz2' % [input_prefix, aid],
        tags: ['gauss%d' % aid, 'wiener%d' % aid],
      },
    }, nin=0, nout=1);
    g.pipeline([src, img_maker.per_anode(anode)],
               'img_graph_anode%d' % aid);

  local graphs = [per_anode_graph(a) for a in anodes];

  // Collect edges and component nodes from all per-anode subgraphs.
  // Pgrapher runs all connected components, so independent chains all execute.
  local all_edges = std.foldl(function(acc, gr) acc + g.edges(gr), graphs, []);
  local all_uses  = std.foldl(function(acc, gr) acc + g.uses(gr),  graphs, []);

  local app = {
    type: 'Pgrapher',
    data: { edges: all_edges },
  };

  local cmdline = {
    type: 'wire-cell',
    data: {
      plugins: [
        'WireCellGen',
        'WireCellPgraph',
        'WireCellSio',
        'WireCellSigProc',
        'WireCellImg',
        'WireCellClus',
      ],
      apps: ['Pgrapher'],
    },
  };

  [cmdline] + all_uses + [app]
