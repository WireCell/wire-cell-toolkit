// Pure WireCell (no LArSoft) NF+SP pipeline.
//
// Part 2 of a two-part split of wcls-nf-sp-out.jsonnet:
//   Part 1: wcls-nf-sp-out.jsonnet  (runs in art/LArSoft)
//     - reads RawDigits, runs ChannelSelector, saves per-anode orig frames:
//         protodune-orig-frames-anode{N}.tar.bz2
//   Part 2: this file  (runs standalone with wire-cell CLI)
//     - reads per-anode orig frames, runs [Resampler ->] NF -> SP
//     - saves NF frames: protodune-nf-frames-anode{N}.tar.bz2
//     - saves SP frames: protodune-sp-frames-anode{N}.tar.bz2
//
// Run example (all anodes):
//   wire-cell -l stdout -L debug \
//     --tla-str orig_prefix="protodune-orig-frames" \
//     --tla-str sp_prefix="protodune-sp-frames" \
//     --tla-str reality="data" \
//     -c pgrapher/experiment/protodunevd/wct-nf-sp.jsonnet
//
// To process a subset of anodes:
//   wire-cell ... --tla-code anode_indices='[4,5]' -c wct-nf-sp.jsonnet
//
// Data vs simulation:
//   reality='data' (default) inserts a Resampler (512 ns -> 500 ns) on the
//   four bottom-drift anodes (n<4) before NF, and re-stamps the top-drift
//   anodes (n>=4) to 500 ns at read time -- the upstream extraction
//   mislabels the (physically 500 ns) top frames as 512 ns. Pass
//   reality='sim' to skip both, for simulated input already at 500 ns.

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

local params = import 'pgrapher/experiment/protodunevd/params.jsonnet';

local tools_maker = import 'pgrapher/common/tools.jsonnet';
local tools_all = tools_maker(params);

function(
  orig_prefix   = 'protodune-orig-frames',  // input prefix; reads {prefix}-anode{N}.tar.bz2
  raw_prefix    = 'protodune-sp-frames-raw',   // output prefix for NF (raw) frames
  sp_prefix     = 'protodune-sp-frames',    // output prefix for SP frames
  reality       = 'data',                   // 'data' enables the 512->500 ns Resampler on bottom anodes (n<4); 'sim' disables it
  sigoutform    = 'dense',                 // 'sparse' or 'dense'
  anode_indices = std.range(0, std.length(tools_all.anodes) - 1),
  use_freqmask  = true,                    // apply per-channel frequency mask in NF; override with --tla-code use_freqmask=false
  debug_dump_path = '',                    // when non-empty, PDVDCoherentNoiseSub dumps per-group .npz under this dir (default OFF)
  debug_dump_groups = [],                  // optional whitelist of group ids (= first-channel idents). [] = all groups
  shield_dump_path = '',                   // when non-empty, PDVDShieldCouplingSub dumps diagnostic npz per group here (default OFF)
  // L1SP defaults: tagger ON in dump mode, LASSO writeback OFF.  Users can
  // validate the ROI tagger via per-event NPZ dumps before the per-region
  // kernel files are generated.  Switch to 'process' (with kernels_file
  // populated in cfg/.../protodunevd/sp.jsonnet) for full L1SP replacement.
  l1sp_pd_mode = 'dump',                   // 'dump' (default; tagger only) / 'process' (full L1SP) / '' (OFF)
  l1sp_pd_dump_path = '',                  // dump directory (used when l1sp_pd_mode='dump'); pass via -c flag
  l1sp_pd_wf_dump_path = '',               // per-ROI waveform dump dir (process mode); pass via -w flag
  l1sp_pd_planes = [0, 1],                 // plane indices to process (0=U, 1=V; skip W)
  l1sp_pd_adj_enable = true,               // cross-channel adjacency expansion (default ON)
  l1sp_pd_adj_max_hops = 3,                // adjacency hop cap (default 3 = +/-3 channels from any donor)
  // Special debug mode: also dump the pre-Wire-filter, pre-ROI deconvolved
  // waveform (h{u,v,w}_rawdecon<ident> in the magnify ROOT) for offline
  // software-filter tuning.  OFF in production.  Pass via -R in run_nf_sp_evt.sh.
  dump_rawdecon = false,
  // ROI-debug mode: run SP with use_roi_debug_mode + multi-plane protection so
  // OmnibusSigProc emits the DNN-ROI input tags (loose_lf / tight_lf / mp2_roi /
  // mp3_roi / decon_charge, alongside gauss) into the SP frame archive.
  // OFF in production.
  roi_debug = false,
)

  local tools = tools_all;
  local use_resampler = (reality == 'data');

  local base = import 'pgrapher/experiment/protodunevd/chndb-base.jsonnet';
  local chndb = [{
    type: 'OmniChannelNoiseDB',
    name: 'ocndbperfect%d' % n,
    data: base(params, tools.anodes[n], tools.field, tools.anodes[n].data.ident, use_freqmask=use_freqmask) { dft: wc.tn(tools.dft) },
    uses: [tools.anodes[n], tools.field, tools.dft],
  } for n in std.range(0, std.length(tools.anodes) - 1)];

  local nf_maker = import 'pgrapher/experiment/protodunevd/nf.jsonnet';
  local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], tools.anodes[n].data.ident, name='nf%d' % tools.anodes[n].data.ident,
                             debug_dump_path=debug_dump_path, debug_dump_groups=debug_dump_groups,
                             shield_dump_path=shield_dump_path)
                    for n in std.range(0, std.length(tools.anodes) - 1)];

  local sp_maker = import 'pgrapher/experiment/protodunevd/sp.jsonnet';
  // Diagnostic: dump the 2D (wire x time-freq) input/response/decon spectra
  // from inside OmnibusSigProc::decon_2D_init() into NPZ files.  Off in
  // production; enabled here for the V-plane pole investigation.
  local sp = sp_maker(params, tools, { sparse: sigoutform == 'sparse',
                                       dump_2d_spectra: true,
                                       dump_2d_prefix: 'dumps_data/sp_dump' }
    + (if roi_debug
       then { use_roi_debug_mode: true, use_multi_plane_protection: true, mp_tick_resolution: 4 }
       else {}));
  local sp_pipes = [sp.make_sigproc(a,
                                    l1sp_pd_mode=l1sp_pd_mode,
                                    l1sp_pd_dump_path=l1sp_pd_dump_path,
                                    l1sp_pd_wf_dump_path=l1sp_pd_wf_dump_path,
                                    l1sp_pd_planes=l1sp_pd_planes,
                                    l1sp_pd_adj_enable=l1sp_pd_adj_enable,
                                    l1sp_pd_adj_max_hops=l1sp_pd_adj_max_hops,
                                    dump_rawdecon=dump_rawdecon)
                    for a in tools.anodes];

  local resamplers_config = import 'pgrapher/common/resamplers.jsonnet';
  local load_resamplers = resamplers_config(g, wc, tools);
  local resamplers = load_resamplers.resamplers;

  // Tap: save NF output (raw) frame per anode
  local raw_frame_tap = function(anode_ident)
    g.fan.tap('FrameFanout',
      g.pnode({
        type: 'FrameFileSink',
        name: 'rawframesink%d' % anode_ident,
        data: {
          outname: '%s-anode%d.tar.bz2' % [raw_prefix, anode_ident],
          tags: ['raw%d' % anode_ident],
          digitize: false,
          masks: false,
        },
      }, nin=1, nout=0),
      'rawframetap%d' % anode_ident);

  // Tap: save SP output (gauss+wiener) frame per anode
  local frame_tap = function(anode_ident)
    g.fan.tap('FrameFanout',
      g.pnode({
        type: 'FrameFileSink',
        name: 'spframesink%d' % anode_ident,
        data: {
          outname: '%s-anode%d.tar.bz2' % [sp_prefix, anode_ident],
          tags: ['gauss%d' % anode_ident, 'wiener%d' % anode_ident]
                + (if dump_rawdecon then ['rawdecon%d' % anode_ident] else [])
                + (if roi_debug then [t % anode_ident for t in
                     ['loose_lf%d', 'tight_lf%d', 'mp2_roi%d', 'mp3_roi%d', 'decon_charge%d']]
                   else []),
          digitize: false,
          masks: true,
        },
      }, nin=1, nout=0),
      'spframetap%d' % anode_ident);

  // Build one source -> [resampler ->] NF -> SP pipeline per anode
  local per_anode_graph(n) =
    local anode = tools.anodes[n];
    local aid = anode.data.ident;

    local src = g.pnode({
      type: 'FrameFileSource',
      name: 'origframesrc%d' % aid,
      data: {
        inname: '%s-anode%d.tar.bz2' % [orig_prefix, aid],
        tags: ['orig'],
      } + (
        // The upstream extraction stamps every anode's orig frames with a
        // single tick (512 ns in data mode).  That is correct for the
        // bottom CRP (n<4, handled by the Resampler below) but wrong for
        // the top CRP (n>=4), which is physically digitized at 500 ns.
        // Re-stamp the top tick at read time so NF/SP run on the correct
        // 500 ns grid.  This is a label correction, not a resample.
        if use_resampler && n >= 4 then { tick: 500 * wc.ns } else {}
      ),
    }, nin=0, nout=1);

    local sink = g.pnode({ type: 'DumpFrames', name: 'dump%d' % aid }, nin=1, nout=0);

    g.pipeline(
      [src]
      + (if use_resampler && n < 4 then [resamplers[n]] else [])
      + [nf_pipes[n]]
      + [raw_frame_tap(aid)]
      + [sp_pipes[n]]
      + [frame_tap(aid)]
      + [sink],
      'nfsp_pipe_%d' % n);

  local graphs = [per_anode_graph(n) for n in anode_indices];

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
        'WireCellAux',
      ],
      apps: ['Pgrapher'],
    },
  };

  [cmdline] + all_uses + [app]
