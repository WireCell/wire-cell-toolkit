// CANONICAL in-tree SBND imaging job — single source of truth.
// Promoted from sbnd_xin 2026-07-27 (sbnd_xin/docs/64_cfg-sync.md);
// sbnd_xin/wct-img-all.jsonnet is now a one-line re-export of this file.
// Runnable directly:
//   wire-cell -c pgrapher/experiment/sbnd/wct-img-all.jsonnet --tla-...
// (WIRECELL_PATH must contain toolkit/cfg + wire-cell-data).
//
// NOTE full_deghost defaults TRUE here, which is the SBND production setting;
// the img.jsonnet module's own per_anode default stays false so other callers
// are unaffected.
//
// Run 3D imaging on SBND SP frames — standalone (no LArSoft).
//
// Input:  sp-frames.tar.bz2 — single-event archive (trace tag 'dnnsp', both anodes)
// Output: icluster-apa<N>-active.npz, icluster-apa<N>-masked.npz  (in output_dir)
//
// The shared imaging graph (pgrapher/experiment/sbnd/img.jsonnet) expects per-anode
// 'gauss<N>' and 'wiener<N>' trace tags.  A FrameFanout with per-anode trace rules
// duplicates the single 'dnnsp' tag into both 'gauss<N>' and 'wiener<N>' for each
// anode N.  The 'bad' channel mask flows through as-is.
//
// Usage (called from run_img_evt.sh):
//   wire-cell \
//     --tla-str  input=work/evt2/sp-frames.tar.bz2 \
//     --tla-code anode_indices='[0,1]' \
//     --tla-str  output_dir=work/evt2 \
//     -c wct-img-all.jsonnet

local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';

local params = import 'pgrapher/experiment/sbnd/simparams.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local tools_all = tools_maker(params);

function(
  input         = 'sp-frames.tar.bz2',
  anode_indices = [0, 1],
  output_dir    = '',
  full_deghost  = true,   // matches uBooNE chain (ProjectionDeghosting x2 + 3 ChargeSolving + 3 InSliceDeghosting); pass --tla-code full_deghost=false to revert to simple-solving
)
  local anodes  = [tools_all.anodes[i] for i in anode_indices];
  local nanodes = std.length(anodes);

  local img = import 'pgrapher/experiment/sbnd/img.jsonnet';
  local img_maker = img();

  local img_pipes = [img_maker.per_anode(anodes[n], 'multi-3view', add_dump=false, full_deghost=full_deghost)
                     for n in std.range(0, nanodes - 1)];

  // ClusterFileSink helpers — port 0 = active (live), port 1 = masked (dead)
  local prefix(aid) = if output_dir == '' then '' else output_dir + '/';
  local cfsink(fname) = g.pnode({
    type: 'ClusterFileSink',
    name: fname,
    data: { format: 'numpy', outname: fname },
  }, nin=1, nout=0);

  local cfsinks_active = [cfsink('%sicluster-apa%d-active.npz' % [prefix(anodes[n].data.ident), anodes[n].data.ident])
                          for n in std.range(0, nanodes - 1)];
  local cfsinks_masked = [cfsink('%sicluster-apa%d-masked.npz' % [prefix(anodes[n].data.ident), anodes[n].data.ident])
                          for n in std.range(0, nanodes - 1)];

  // Wire each per-anode imaging pipeline's ports 0/1 to the two ClusterFileSinks.
  // (img.jsonnet's internal chsel_pipes already restricts each branch to APA N's
  // 5638-channel range, so no extra pre-filter is needed here.)
  local img_dump_pipe = [g.intern(
    innodes=  [img_pipes[n]],
    centernodes= [],
    outnodes= [cfsinks_active[n], cfsinks_masked[n]],
    edges=    [
      g.edge(img_pipes[n], cfsinks_active[n], 0, 0),
      g.edge(img_pipes[n], cfsinks_masked[n], 1, 0),
    ]
  ) for n in std.range(0, nanodes - 1)];

  // FrameFanout rules: rename 'dnnsp' → 'gauss<N>' and 'wiener<N>' for anode N.
  // Producing two output tags from one input duplicates the trace for both imaging inputs.
  local fanout_rules = [
    {
      frame: { '.*': 'orig%d' % anodes[n].data.ident },
      trace: { dnnsp: ['gauss%d' % anodes[n].data.ident, 'wiener%d' % anodes[n].data.ident] },
    }
    for n in std.range(0, nanodes - 1)
  ];

  local src = g.pnode({
    type: 'FrameFileSource',
    name: 'frame_source',
    data: { inname: input, tags: ['dnnsp'] },
  }, nin=0, nout=1);

  local fanout_graph = f.fanout('FrameFanout', img_dump_pipe, 'img_fanout', fanout_rules);
  local graph = g.pipeline([src, fanout_graph], 'main');

  local app = {
    type: 'Pgrapher',
    data: { edges: g.edges(graph) },
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
        'WireCellRoot',
      ],
      apps: ['Pgrapher'],
    },
  };

  [cmdline] + g.uses(graph) + [app]
