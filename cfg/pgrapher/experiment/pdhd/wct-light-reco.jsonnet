// WCT-native PDHD light reconstruction for one event (stage 2; see
// toolkit flash/docs/stage2-reconstruction.md):
//   PDHDOpWaveformSource (raw snippets) -> OpDecon (Wiener, SPE
//   template) -> OpHitFinder (SlidingWindow) -> OpFlashFinder
//   (OpFlashAlg) -> opflash_pdhd-wct.tar.gz
// plus a dump of the deconvolved snippet frames for validation
// (light-frames-wct.tar.bz2, tags "raw"/"decon").
//
// Compile/run:
//   wcsonnet -A input_file=... -A output_dir=... -S run=27305 -S event=150 \
//            -o cfg.json wct-light-reco.jsonnet
//   wire-cell -l stderr -L debug -c cfg.json

local g = import 'pgraph.jsonnet';
local flash = import 'pgrapher/experiment/pdhd/flash.jsonnet';

// offset_us: per-event readout-vs-trigger offset (tc - rd_timestamp, ~250) read
// from the ROOT trigoff tree by run_light_evt.sh and stamped into the opflash
// metadata so downstream charge clustering / Q-L matching can place charge on the
// trigger time base.  Default 0 (no offset).
function(input_file, output_dir='.', run=27305, event=150, offset_us=0)

  local run_n = if std.type(run) == 'string' then std.parseInt(run) else run;
  local evt_n = if std.type(event) == 'string' then std.parseInt(event) else event;
  local off_us = if std.type(offset_us) == 'string' then std.parseJson(offset_us) else offset_us;
  // Single node instance reused below so all references share identity.
  local opflash_finder = flash.opflash_finder(offset_us=off_us);

  local fanout = g.pnode({
    type: 'FrameFanout',
    name: 'light_reco',
    data: { multiplicity: 2 },
  }, nin=1, nout=2);

  local graph = g.intern(
    innodes=[flash.opwaveform_source(input_file, run_n, evt_n)],
    centernodes=[
      flash.opdecon(),
      fanout,
      flash.ophit(),
      opflash_finder,
    ],
    outnodes=[
      flash.waveform_sink('%s/light-frames-wct.tar.bz2' % output_dir, tags=['raw', 'decon'], name='wct'),
      flash.opflash_sink('%s/opflash_pdhd-wct.tar.gz' % output_dir, name='wct'),
    ],
    edges=[
      g.edge(flash.opwaveform_source(input_file, run_n, evt_n), flash.opdecon()),
      g.edge(flash.opdecon(), fanout),
      g.edge(fanout, flash.waveform_sink('%s/light-frames-wct.tar.bz2' % output_dir, tags=['raw', 'decon'], name='wct'), 0, 0),
      g.edge(fanout, flash.ophit(), 1, 0),
      g.edge(flash.ophit(), opflash_finder),
      g.edge(opflash_finder, flash.opflash_sink('%s/opflash_pdhd-wct.tar.gz' % output_dir, name='wct')),
    ],
  );

  local app = {
    type: 'Pgrapher',
    data: { edges: g.edges(graph) },
  };

  local cmdline = {
    type: 'wire-cell',
    data: {
      plugins: [
        'WireCellRoot',
        'WireCellFlash',
        'WireCellGen',
        'WireCellSio',
        'WireCellAux',
        'WireCellPgraph',
      ],
      apps: ['Pgrapher'],
    },
  };

  [cmdline] + g.uses(graph) + [app]
