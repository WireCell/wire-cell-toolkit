// Per-anode per-plane DNN-ROI subgraph for ProtoDUNE-VD.
//
// Uses two single-plane DNNROIFinding nodes (U, V) sharing the same
// TorchService.  Each call feeds the model an (1, 6, ~476, 1600) tensor;
// the two calls are serialized by the TorchService semaphore so peak
// activation memory is for one plane at a time.
//
// The PDVD .ts was trained on stacked U+V at (1, 6, ~952, 1600).  The
// MobileNetV3-large UNet is fully convolutional on the channel axis, so
// per-plane inference is structurally compatible — analogous to PDHD's
// per-plane mode which has been validated by
// DNN_ROI_SP/scripts/test_per_plane_ts.py.  An equivalent per-plane vs
// stacked side-by-side validation for PDVD is still pending; treat
// per-plane PDVD output as scientifically un-cross-checked until that
// comparison lands.
//
// Per-anode policy: U and V via model, W from standard SP gauss
// passthrough.  PDVD has no APA0 special case — all 8 anodes use the
// same model for both U and V.
//
// PDVD's traced UNet has 5 stride-2 levels post-rebin and needs the
// post-rebin tick width divisible by 32, i.e. input ticks divisible by
// 128.  Default tick_pad_multiple=128 so the C++ node zero-pads any
// input length to the next 128-multiple before inference (output is
// cropped back to input_ticks).
//
// Output frame trace tags per APA N:
//   - dnnspNu : DNN-ROI on U  (model output)
//   - dnnspNv : DNN-ROI on V  (model output)
//   - dnnspNw : standard SP gauss on W  (always passthrough)
// Unified frame tag "dnnspN" via the final Retagger.

local pg = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

function(anode, ts, prefix='dnnroi',
         output_scale=1.0,
         nticks=6000,
         tick_per_slice=4,
         nchunks=1,
         mask_thresh=0.2,
         nchan=6,
         debugfile='',
         tick_pad_multiple=128)

  local apaid = anode.data.ident;
  local prename = prefix + std.toString(apaid);

  // Input trace tags fed to the model, in the order the model expects.
  // PDVD ships 6-channel models (DAGMan 287); order must match the
  // training im_tags exactly.
  assert nchan == 6 : 'pdvd dnnroi nchan must be 6, got %d' % nchan;
  local apa_intags = [
    'loose_lf%d' % apaid,
    'mp2_roi%d' % apaid,
    'mp3_roi%d' % apaid,
    'tight_lf%d' % apaid,
    'decon_charge%d' % apaid,
    'gauss%d' % apaid,
  ];

  // Per-anode per-call debug filename suffix so U / V calls don't collide.
  local _dbg(suffix) =
    if debugfile == '' then '' else '%s_%s' % [debugfile, suffix];

  local dnn_node(plane, outtag, name_suffix, dbg_suffix) = pg.pnode({
    type: 'DNNROIFinding',
    name: prename + name_suffix,
    data: {
      anode: wc.tn(anode),
      plane: plane,
      intags: apa_intags,
      decon_charge_tag: 'decon_charge%d' % apaid,
      outtag: outtag,
      // Propagate the per-channel Wiener threshold (carried as the
      // trace_summary of wiener%d in OmnibusSigProc output) onto the
      // output dnnsp%d* trace tag so the downstream Magnify pipeline
      // can build hu/hv/hw_threshold{N} TH1F.
      summary_tag: 'wiener%d' % apaid,
      // PDVD 6-ch models bake per-channel z-scale normalization into
      // the .ts, so they run with input_scale 1.0.
      input_scale: 1.0,
      output_scale: output_scale,
      mask_thresh: mask_thresh,
      forward: wc.tn(ts),
      nticks: nticks,
      tick_per_slice: tick_per_slice,
      tick_pad_multiple: tick_pad_multiple,
      nchunks: nchunks,
      debugfile: _dbg(dbg_suffix),
    },
  }, nin=1, nout=1, uses=[ts, anode]);

  // Gauss passthrough for W: re-tag gauss%d → outtag.  summary_tags
  // pulls the wiener%d trace_summary (= per-channel threshold) alongside
  // the gauss waveforms so the threshold survives into the downstream
  // Magnify TH1F.
  local pass_through(plane, outtag, name_suffix) = pg.pnode({
    type: 'PlaneSelector',
    name: prename + name_suffix,
    data: {
      anode: wc.tn(anode),
      plane: plane,
      tags: ['gauss%d' % apaid],
      summary_tags: ['wiener%d' % apaid],
      tag_rules: [{
        frame: { '.*': 'DNNROIFinding' },
        trace: { ['gauss%d' % apaid]: outtag },
      }],
    },
  }, nin=1, nout=1, uses=[anode]);

  local u_node = dnn_node(0, 'dnnsp%du' % apaid, 'u', 'u');
  local v_node = dnn_node(1, 'dnnsp%dv' % apaid, 'v', 'v');
  local w_node = pass_through(2, 'dnnsp%dw' % apaid, '_wpass');

  local sub_nodes = [u_node, v_node, w_node];
  local mult = std.length(sub_nodes);  // always 3

  local dnnfanout = pg.pnode({
    type: 'FrameFanout',
    name: prename,
    data: { multiplicity: mult },
  }, nin=1, nout=mult);

  local plane_letters = ['u', 'v', 'w'];

  local dnnfanin = pg.pnode({
    type: 'FrameFanin',
    name: prename,
    data: {
      multiplicity: mult,
      tag_rules: [
        { frame: { '.*': 'dnnsp%d%s' % [apaid, plane] },
          trace: { '.*': 'dnnsp%d%s' % [apaid, plane] } }
        for plane in plane_letters
      ],
    },
  }, nin=mult, nout=1);

  local retagger = pg.pnode({
    type: 'Retagger',
    name: 'dnnroi%d' % apaid,
    data: {
      tag_rules: [{
        frame: { '.*': 'dnnsp%d' % apaid },
        merge: { '.*': 'dnnsp%d' % apaid },
      }],
    },
  }, nin=1, nout=1);

  pg.intern(
    innodes=[dnnfanout],
    outnodes=[retagger],
    centernodes=sub_nodes + [dnnfanin],
    edges=[pg.edge(dnnfanout, sub_nodes[i], i, 0) for i in std.range(0, mult - 1)]
          + [pg.edge(sub_nodes[i], dnnfanin, 0, i) for i in std.range(0, mult - 1)]
          + [pg.edge(dnnfanin, retagger, 0, 0)],
  )
