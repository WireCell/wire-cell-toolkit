// Per-anode per-plane DNN-ROI subgraph for PDHD.
//
// Uses two stock single-plane DNNROIFinding nodes (U, V) sharing the
// same CP43.ts (via TorchService).  Each call feeds the model an
// (1, 3, 800, 1500) tensor; the two calls are serialized by the
// TorchService semaphore so peak activation memory is for one plane
// at a time — roughly half the stacked (1, 3, 1600, 1500) peak.
//
// The training shape was (1, 3, 1600, 1500) (U+V stacked) but the
// MobileNetV3-Large UNet is fully convolutional on the channel axis,
// so the same .ts is shape-flexible.  Empirically verified by
// DNN_ROI_SP/scripts/test_per_plane_ts.py.
//
// Per-anode policy:
//   APA0:    U via model; V and W from standard SP gauss (passthrough).
//   APA1-3:  U and V via model; W from standard SP gauss.
//
// Output frame trace tags per APA N:
//   - dnnspNu : DNN-ROI on U  (model output)
//   - dnnspNv : DNN-ROI on V  (model on APA1-3; gauss passthrough on APA0)
//   - dnnspNw : standard SP gauss on W  (always)
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
         debugfile='')

  local apaid = anode.data.ident;
  local prename = prefix + std.toString(apaid);

  // Input trace tags fed to the model, in the order the model expects.
  // The 3-ch model (CP43.ts) uses the first three; the 6-ch KD/QAT models
  // add tight_lf/decon_charge/gauss.  Order must match the training im_tags.
  assert nchan == 3 || nchan == 6 : 'dnnroi nchan must be 3 or 6, got %d' % nchan;
  local apa_intags = [
    'loose_lf%d' % apaid,
    'mp2_roi%d' % apaid,
    'mp3_roi%d' % apaid,
    'tight_lf%d' % apaid,
    'decon_charge%d' % apaid,
    'gauss%d' % apaid,
  ][0:nchan];

  local emit_v_from_model = (apaid != 0);

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
      // can build hu/hv/hw_threshold0 TH1F.
      summary_tag: 'wiener%d' % apaid,
      // 6-ch models bake per-channel z-scale normalization into the .ts,
      // so they run with input_scale 1.0; the 3-ch model needs 1/4000.
      input_scale: if nchan == 6 then 1.0 else 1.0 / 4000,
      output_scale: output_scale,
      mask_thresh: mask_thresh,
      forward: wc.tn(ts),
      nticks: nticks,
      tick_per_slice: tick_per_slice,
      nchunks: nchunks,
      debugfile: _dbg(dbg_suffix),
    },
  }, nin=1, nout=1, uses=[ts, anode]);

  // Gauss passthrough for a non-DNN plane: re-tag gauss%d → outtag.
  // summary_tags pulls the wiener%d trace_summary (= per-channel
  // threshold) alongside the gauss waveforms so the threshold survives
  // into the downstream Magnify TH1F.
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
  local v_node =
    if emit_v_from_model
    then dnn_node(1, 'dnnsp%dv' % apaid, 'v', 'v')
    else pass_through(1, 'dnnsp%dv' % apaid, '_vpass');
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
