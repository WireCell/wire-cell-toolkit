// Standalone imaging for pdhd.
// Reads per-anode SP frames from archives produced upstream:
//   protodunehd-sp-frames-anode{N}.tar.bz2
//
// Each anode gets an independent FrameFileSource -> imaging pipeline.
// Pgrapher executes all independent source->sink chains.
//
// Run standalone (all anodes):
//   wire-cell -l stdout -L debug -c wct-img-all.jsonnet
//
// With a custom file prefix:
//   wire-cell -l stdout -L debug \
//     --tla-str input_prefix="protodunehd-sp-frames" \
//     -c wct-img-all.jsonnet
//
// To select a subset of anodes by index into tools_all.anodes:
//   wire-cell -l stdout -L debug \
//     --tla-code anode_indices='[0,1]' \
//     -c wct-img-all.jsonnet
//
// To write cluster output files to a specific directory:
//   wire-cell -l stdout -L debug \
//     --tla-str output_dir="work/027409_1" \
//     -c wct-img-all.jsonnet

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

// Imaging uses simparams for its 500 ns tick (data is resampled 512->500 ns
// before imaging).  But simparams redefines det.volumes with BOTH APA faces
// sensitive ("cryostat side included"), which makes the dead/masked fork tile
// phantom dead blobs on the driftless face.  Restore the data geometry's
// det.volumes (one face per APA nulled = non-sensitive) so GridTiling suppresses
// the driftless face; the sensitive face geometry is identical between the two.
local simparams = import 'pgrapher/experiment/pdhd/simparams.jsonnet';
local dataparams = import 'pgrapher/experiment/pdhd/params.jsonnet';
local params = simparams { det+: { volumes: dataparams.det.volumes } };

local tools_maker = import 'pgrapher/common/tools.jsonnet';
local tools_all = tools_maker(params);

// Top-level function: parameters overridable via --tla-str / --tla-code
function(
  // Prefix for per-anode input files: "{input_prefix}-anode{N}.tar.bz2"
  input_prefix = 'protodunehd-sp-frames',
  // Indices into tools_all.anodes to process; default = all
  anode_indices = std.range(0, std.length(tools_all.anodes) - 1),
  // Directory for output cluster files ('' means current directory)
  output_dir = '',
  // Reframer output length in ticks; 0 = fall back to params.daq.nticks.
  // Readout length can vary run to run -- run_img_evt.sh probes the input
  // frame and passes the actual value (see
  // pdvd/docs/sp-img-readout-window-truncation.md).
  nticks = 0
)

  local anodes = [tools_all.anodes[i] for i in anode_indices];

  // PDHD slices on any positive charge (1e-6 = charge>0 surrogate, see the
  // nthreshold note in the cfg img.jsonnet); the cfg default is 3.6 sigma.
  local img = import 'pgrapher/experiment/pdhd/img.jsonnet';
  local img_maker = img({nthreshold: [1e-6, 1e-6, 1e-6]});

  // Build one FrameFileSource + imaging pipeline per anode.
  // Reframer densifies sparse IFrames (where gauss/wiener may have
  // different row counts) to the full anode channel set before imaging.
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
    local reframer = g.pnode({
      type: 'Reframer',
      name: 'reframer_anode%d' % aid,
      data: {
        anode: wc.tn(anode),
        tags: ['gauss%d' % aid, 'wiener%d' % aid],
        tbin: 0,
        nticks: if nticks > 0 then nticks else params.daq.nticks,
        fill: 0.0,
        // Carry the SP frame's bad-channel mask onto the reframed frame so the
        // dead/masked imaging fork can tile dead regions; without this the CMM
        // is dropped and the dead Bee layer comes out empty.  Safe here because
        // tbin=0 keeps the mask tick ranges aligned.
        keep_masks: true,
      },
    }, nin=1, nout=1, uses=[anode]);
    g.pipeline([src, reframer, img_maker.per_anode(anode, "multi", output_dir)],
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
