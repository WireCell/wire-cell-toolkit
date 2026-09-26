// Convert one event of a temporary PDHD light-data ROOT file into WCT
// formats (see toolkit flash/docs/design.md):
//   - opflash_pdhd.tar.gz     opflash tensor set (LArSoft flashes, ophits,
//                             summaries; trigger-relative ns, 160 OpDets)
//   - light-frames.tar.bz2    "raw"/"deconv" waveform snippet frames
//
// Compile/run:
//   wcsonnet -A input_file=... -A output_dir=... -S run=27305 -S event=150 \
//            -o cfg.json wct-light-convert.jsonnet
//   wire-cell -l stderr -L debug -c cfg.json

local g = import 'pgraph.jsonnet';
local flash = import 'pgrapher/experiment/pdhd/flash.jsonnet';

function(input_file, output_dir='.', run=27305, event=150)

  local run_n = if std.type(run) == 'string' then std.parseInt(run) else run;
  local evt_n = if std.type(event) == 'string' then std.parseInt(event) else event;

  local opflash_pipe = g.pipeline([
    flash.opflash_source(input_file, run_n, evt_n),
    flash.opflash_sink('%s/opflash_pdhd.tar.gz' % output_dir),
  ]);

  local waveform_pipe = g.pipeline([
    flash.opwaveform_source(input_file, run_n, evt_n),
    flash.waveform_sink('%s/light-frames.tar.bz2' % output_dir),
  ]);

  local graphs = [opflash_pipe, waveform_pipe];
  local all_edges = std.foldl(function(acc, gr) acc + g.edges(gr), graphs, []);
  local all_uses = std.foldl(function(acc, gr) acc + g.uses(gr), graphs, []);

  local app = {
    type: 'Pgrapher',
    data: { edges: all_edges },
  };

  local cmdline = {
    type: 'wire-cell',
    data: {
      plugins: [
        'WireCellRoot',
        'WireCellSio',
        'WireCellPgraph',
      ],
      apps: ['Pgrapher'],
    },
  };

  [cmdline] + all_uses + [app]
