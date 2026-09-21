// Convert per-anode SP frame archives (protodune-sp-frames-anode{N}.tar.bz2)
// into a single Magnify ROOT file containing all anodes.
//
// Each anode produces: hu/hv/hw_gauss<N>, hu/hv/hw_wiener<N> (TH2F),
//   hu/hv/hw_threshold<N> (TH1F per-channel threshold), T_bad<N> (cmmtree).
// One Trun tree is written by the first anode's sink.
//
// When include_raw=true, the raw frame archive (protodune-sp-frames-raw-anode{N}.tar.bz2)
// is also read and hu/hv/hw_raw<N> TH2F histograms are appended to the same file.
//
// When include_orig=true, the pre-NF orig frame archive
// (protodune-orig-frames-anode{N}.tar.bz2) is also read and
// hu/hv/hw_orig<N> TH2F histograms are appended to the same file.
// The orig archives carry a uniform trace tag 'orig'; a Retagger renames it
// to 'orig<N>' per anode before MagnifySink to avoid histogram-name collisions.
//
// All pipelines run in the same wire-cell process to avoid MagnifySink::create_file()
// wiping the output file between passes.
//
// Typical usage (called from run_sp_to_magnify_evt.sh):
//   wire-cell \
//     --tla-str  input_prefix=/path/to/protodune-sp-frames \
//     --tla-code anode_indices='[0,1,2,3,4,5,6,7]' \
//     --tla-str  output_file=magnify-run039324-evt1.root \
//     --tla-code run=39324 --tla-code subrun=0 --tla-code event=339870 \
//     -c wct-sp-to-magnify.jsonnet
//
// With raw and orig:
//   wire-cell ... --tla-code include_raw=true \
//     --tla-str raw_input_prefix=/path/to/protodune-sp-frames-raw \
//     --tla-code include_orig=true \
//     --tla-str orig_input_prefix=/path/to/protodune-orig-frames \
//     -c wct-sp-to-magnify.jsonnet

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

local params = import 'pgrapher/experiment/protodunevd/simparams.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local tools_all = tools_maker(params);

function(
  input_prefix      = 'protodune-sp-frames',
  anode_indices     = std.range(0, std.length(tools_all.anodes) - 1),
  output_file       = 'magnify.root',
  run               = 0,
  subrun            = 0,
  event             = 0,
  nticks            = 6000,  // fallback only; readout varies by run, run_sp_to_magnify_evt.sh probes and passes the actual value
  include_raw       = true,
  raw_input_prefix  = 'protodune-sp-frames-raw',
  include_orig      = false,
  orig_input_prefix = 'protodune-orig-frames',
  // Special-mode pre-Wire-filter pre-ROI deconvolved waveform.  Set true
  // when wct-nf-sp.jsonnet ran with dump_rawdecon=true so the SP frame
  // archive carries the rawdecon%d tag.
  include_rawdecon  = false,
)
  local anodes  = [tools_all.anodes[i] for i in anode_indices];
  local nanodes = std.length(anodes);

  local runinfo = {
    runNo:          run,
    subRunNo:       subrun,
    eventNo:        event,
    total_time_bin: nticks,
  };

  local mag = (import 'pgrapher/experiment/protodunevd/magnify-sinks.jsonnet')(
    { anodes: anodes }, output_file, runinfo=runinfo);

  // Per-anode decon pipeline: FrameFileSource(gauss,wiener) → Retagger →
  // MagnifySink.  The DNN-ROI archive (protodune-sp-dnnroi-frames) is
  // structurally a standard SP archive -- the L1SP-after-DNN envelope's
  // dnn_relabel writes gauss%d (= DNN-ROI output) and wiener%d (alias
  // carrying the per-channel threshold summary) -- so the same pipeline
  // serves both the plain-SP and DNN-ROI inputs.
  local sp_anode_graph(n, anode) =
    local aid = anode.data.ident;
    local src = g.pnode({
      type: 'FrameFileSource',
      name: 'frame_source_anode%d' % aid,
      data: {
        inname: '%s-anode%d.tar.bz2' % [input_prefix, aid],
        tags: ['gauss%d' % aid, 'wiener%d' % aid]
              + (if include_rawdecon then ['rawdecon%d' % aid] else []),
      },
    }, nin=0, nout=1);
    // Retagger duplicates wiener<N> → [wiener<N>, threshold<N>] so
    // MagnifySink can write h[uvw]_threshold<N> TH1F (Magnify expects that name).
    local retag = g.pnode({
      type: 'Retagger',
      name: 'retagger_anode%d' % aid,
      data: {
        tag_rules: [{
          frame: { '.*': '' },
          trace: {
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

  // Per-anode orig pipeline: FrameFileSource(orig) → Retagger(orig→orig<N>)
  // → MagnifySink(UPDATE) → DumpFrames.
  // The upstream tag in the archive is the uniform string 'orig' (not 'orig<N>'),
  // so a Retagger is needed to give each anode a distinct histogram name in Magnify.
  local orig_anode_graph(anode) =
    local aid = anode.data.ident;
    local src = g.pnode({
      type: 'FrameFileSource',
      name: 'orig_source_anode%d' % aid,
      data: {
        inname: '%s-anode%d.tar.bz2' % [orig_input_prefix, aid],
        tags: ['orig'],
      },
    }, nin=0, nout=1);
    local retag = g.pnode({
      type: 'Retagger',
      name: 'orig_retagger_anode%d' % aid,
      data: {
        tag_rules: [{
          frame: { '.*': '' },
          trace: { orig: 'orig%d' % aid },
        }],
      },
    }, nin=1, nout=1);
    local sink = g.pnode({
      type: 'MagnifySink',
      name: 'magorig%d' % aid,
      data: {
        output_filename: output_file,
        root_file_mode: 'UPDATE',
        frames: ['orig%d' % aid],
        trace_has_tag: true,
        anode: wc.tn(anode),
      },
    }, nin=1, nout=1, uses=[anode]);
    local dump = g.pnode({
      type: 'DumpFrames',
      name: 'origdump_anode%d' % aid,
    }, nin=1, nout=0);
    g.pipeline([src, retag, sink, dump], 'orig_graph_anode%d' % aid);

  local decon_graphs = [sp_anode_graph(n, anodes[n])
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
