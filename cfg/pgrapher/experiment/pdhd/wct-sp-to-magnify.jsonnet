// Convert per-anode SP frame archives (protodunehd-sp-frames-anode{N}.tar.bz2)
// into a single Magnify ROOT file containing all anodes.
//
// Each anode produces: hu/hv/hw_gauss<N>, hu/hv/hw_wiener<N> (TH2F),
//   hu/hv/hw_threshold<N> (TH1F per-channel threshold), T_bad<N> (cmmtree).
// One Trun tree is written by the first anode's sink.
//
// When include_raw=true, the raw frame archive (protodunehd-sp-frames-raw-anode{N}.tar.bz2)
// is also read and hu/hv/hw_raw<N> TH2F histograms are appended to the same file.
//
// When include_orig=true, the pre-NF orig frame archive
// (protodunehd-orig-frames-anode{N}.tar.bz2) is also read and
// hu/hv/hw_orig<N> TH2F histograms are appended to the same file.
// pdhd orig archives store frames as frame_*_<evt>.npy (literal '*' in the tag
// slot). WCT FrameFileSource treats '*' the same as empty string — both are
// "untagged" — so no trace tag is attached.  Because there is nothing to rename,
// no Retagger is used; MagnifySink is configured with trace_has_tag:false so it
// writes h[uvw]_orig<N> from all traces in the frame keyed only by the frames label.
//
// All pipelines run in the same wire-cell process to avoid MagnifySink::create_file()
// wiping the output file between passes.
//
// Typical usage (called from run_sp_to_magnify_evt.sh):
//   wire-cell \
//     --tla-str  input_prefix=/path/to/protodunehd-sp-frames \
//     --tla-code anode_indices='[0,1,2,3]' \
//     --tla-str  output_file=magnify-run027409-evt1.root \
//     --tla-code run=27409 --tla-code subrun=0 --tla-code event=339870 \
//     -c wct-sp-to-magnify.jsonnet
//
// With raw and orig:
//   wire-cell ... --tla-code include_raw=true \
//     --tla-str raw_input_prefix=/path/to/protodunehd-sp-frames-raw \
//     --tla-code include_orig=true \
//     --tla-str orig_input_prefix=/path/to/protodunehd-orig-frames \
//     -c wct-sp-to-magnify.jsonnet

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

local params = import 'pgrapher/experiment/pdhd/simparams.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local tools_all = tools_maker(params);

function(
  input_prefix      = 'protodunehd-sp-frames',
  anode_indices     = std.range(0, std.length(tools_all.anodes) - 1),
  output_file       = 'magnify.root',
  run               = 0,
  subrun            = 0,
  event             = 0,
  // Used only for Trun's total_time_bin metadata; TH2F binning is driven by the
  // actual frame shape from FrameFileSource.  Callers should pass the real frame
  // tick count extracted from frame_gauss0_*.npy or frame_dnnsp0_*.npy.
  nticks            = 6000,
  include_raw       = true,
  raw_input_prefix  = 'protodunehd-sp-frames-raw',
  include_orig      = false,
  orig_input_prefix = 'protodunehd-orig-frames',
  // Special-mode pre-Wire-filter pre-ROI deconvolved waveform.  Set true
  // when the SP frame archive carries the rawdecon%d tag (i.e. wct-nf-sp.jsonnet
  // ran with dump_rawdecon=true).
  include_rawdecon  = false,
  // Set true when the dnnroi archive was produced with L1SP off (carries dnnsp%d
  // tags instead of gauss%d/wiener%d).  dnnsp%d is re-tagged to gauss%d so
  // MagnifySink writes the standard hu/hv/hw_gauss%d histograms.
  use_dnnsp         = false,
)
  local anodes  = [tools_all.anodes[i] for i in anode_indices];
  local nanodes = std.length(anodes);

  local runinfo = {
    runNo:          run,
    subRunNo:       subrun,
    eventNo:        event,
    total_time_bin: nticks,
  };

  // decon_trun_pipe = Trun-aware decon sinks (RECREATE + runinfo on the first
  // sink, threshold summaries, T_bad cmm trees); plain decon_pipe would
  // silently drop Trun/threshold/cmm output.
  local mag = (import 'pgrapher/experiment/pdhd/magnify-sinks.jsonnet')(
    { anodes: anodes }, output_file, runinfo=runinfo);

  local per_anode_graph(n, anode) =
    local aid = anode.data.ident;
    local src = g.pnode({
      type: 'FrameFileSource',
      name: 'frame_source_anode%d' % aid,
      data: {
        inname: '%s-anode%d.tar.bz2' % [input_prefix, aid],
        tags: if use_dnnsp
              then ['dnnsp%d' % aid]
              else ['gauss%d' % aid, 'wiener%d' % aid]
                   + (if include_rawdecon then ['rawdecon%d' % aid] else []),
      },
    }, nin=0, nout=1);
    // Retagger: for L1SP-on, duplicate wiener<N> → [wiener<N>, threshold<N>].
    // For L1SP-off (use_dnnsp), re-tag dnnsp<N> → gauss<N> so MagnifySink writes
    // the standard hu/hv/hw_gauss<N> histograms.
    local retag = g.pnode({
      type: 'Retagger',
      name: 'retagger_anode%d' % aid,
      data: {
        tag_rules: [{
          frame: { '.*': '' },
          trace: if use_dnnsp
                 then { ['dnnsp%d' % aid]: 'gauss%d' % aid }
                 else {
                   ['gauss%d' % aid]: 'gauss%d' % aid,
                   ['wiener%d' % aid]: ['wiener%d' % aid, 'threshold%d' % aid],
                 } + (if include_rawdecon then { ['rawdecon%d' % aid]: 'rawdecon%d' % aid } else {}),
        }],
      },
    }, nin=1, nout=1);
    local dump = g.pnode({
      type: 'DumpFrames',
      name: 'dump_anode%d' % aid,
    }, nin=1, nout=0);
    g.pipeline([src, retag, mag.decon_trun_pipe[n], dump],
               'magnify_graph_anode%d' % aid);

  // Per-anode raw pipeline: FrameFileSource(raw) → MagnifySink(UPDATE) → DumpFrames.
  // Always UPDATE (decon anode-0 sink already RECREATEd the file).
  local raw_anode_graph(anode) =
    local aid = anode.data.ident;
    local src = g.pnode({
      type: 'FrameFileSource',
      name: 'raw_source_anode%d' % aid,
      data: {
        inname: '%s-anode%d.tar.bz2' % [raw_input_prefix, aid],
        tags: ['raw%d' % aid],
      },
    }, nin=0, nout=1);
    local sink = g.pnode({
      type: 'MagnifySink',
      name: 'magraw%d' % aid,
      data: {
        output_filename: output_file,
        root_file_mode: 'UPDATE',
        frames: ['raw%d' % aid],
        trace_has_tag: true,
        anode: wc.tn(anode),
      },
    }, nin=1, nout=1, uses=[anode]);
    local dump = g.pnode({
      type: 'DumpFrames',
      name: 'rawdump_anode%d' % aid,
    }, nin=1, nout=0);
    g.pipeline([src, sink, dump], 'raw_graph_anode%d' % aid);

  // Per-anode orig pipeline: FrameFileSource(orig) → MagnifySink(UPDATE) → DumpFrames.
  // pdhd orig archives have literal '*' in the filename tag slot, which WCT treats as
  // "untagged" — no trace tags are attached.  Use trace_has_tag:false so MagnifySink
  // writes h[uvw]_orig<N> from all traces without needing a trace-tag lookup.
  local orig_anode_graph(anode) =
    local aid = anode.data.ident;
    local src = g.pnode({
      type: 'FrameFileSource',
      name: 'orig_source_anode%d' % aid,
      data: {
        inname: '%s-anode%d.tar.bz2' % [orig_input_prefix, aid],
        tags: [],  // load all traces regardless of tag (pdhd uses literal '*' as tag)
      },
    }, nin=0, nout=1);
    local sink = g.pnode({
      type: 'MagnifySink',
      name: 'magorig%d' % aid,
      data: {
        output_filename: output_file,
        root_file_mode: 'UPDATE',
        frames: ['orig%d' % aid],
        trace_has_tag: false,  // pdhd '*'-tagged traces carry no trace tag; take all traces
        anode: wc.tn(anode),
      },
    }, nin=1, nout=1, uses=[anode]);
    local dump = g.pnode({
      type: 'DumpFrames',
      name: 'origdump_anode%d' % aid,
    }, nin=1, nout=0);
    g.pipeline([src, sink, dump], 'orig_graph_anode%d' % aid);

  local decon_graphs = [per_anode_graph(n, anodes[n])
                        for n in std.range(0, nanodes - 1)];
  local raw_graphs   = if include_raw  then [raw_anode_graph(anodes[n])
                                             for n in std.range(0, nanodes - 1)]
                       else [];
  local orig_graphs  = if include_orig then [orig_anode_graph(anodes[n])
                                             for n in std.range(0, nanodes - 1)]
                       else [];
  // Node ordering controls execution order via reverse Kahn-sort:
  //   orig_graphs (lowest instance) → run LAST  (UPDATE)
  //   raw_graphs  (middle instance) → run SECOND (UPDATE)
  //   decon_graphs (highest instance) → run FIRST (RECREATE, creates the file)
  local graphs = orig_graphs + raw_graphs + decon_graphs;

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
        'WireCellRoot',
      ],
      apps: ['Pgrapher'],
    },
  };

  [cmdline] + all_uses + [app]
