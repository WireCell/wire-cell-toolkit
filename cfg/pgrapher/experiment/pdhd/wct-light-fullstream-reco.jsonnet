// WCT-native PDHD light reconstruction for the -x FULL-STREAM PDs (opch
// 120-159), one event.  Same chain as wct-light-reco.jsonnet but on the
// continuous 343808-sample stream instead of 1024-tick self-trigger snippets:
//   PDHDOpWaveformSource (one long "snippet"/ch from fullstream_to_decoana.py)
//   -> OpDecon (samples=343808, FIXED Wiener filter fixed_snr -> same filter as
//      the snippet path) -> OpHitFinder (SlidingWindow over the full stream)
//   -> OpFlashFinder -> opflash_pdhd-fullstream-wct.tar.gz
// plus a dump of the deconvolved full-stream frame for tuning/validation
// (light-frames-fullstream-wct.tar.bz2, tags "raw"/"decon").
//
// The input is the converted ROOT file (decoana layout) produced by
// pdhd/pd_plot/fullstream_to_decoana.py.  The full-stream channels have no
// LArSoft OpHits, so PDHDOpWaveformSource falls back to t_first = rd_timestamp
// (the stream origin) -- frame_time = (rd_timestamp - tc_time)*16ns, the same
// trigger-relative clock as the snippet path (same tc_time), so the full-stream
// OpHits/flashes are directly time-comparable for coincidence.
//
// Compile/run (see run_light_fullstream_evt.sh):
//   wcsonnet -A input_file=fs.root -A output_dir=... -S run=27980 -S event=8 \
//            -S offset_us=249.808 -o cfg.json wct-light-fullstream-reco.jsonnet
//   wire-cell -l stderr -L debug -c cfg.json

local g = import 'pgraph.jsonnet';
local flash = import 'pgrapher/experiment/pdhd/flash.jsonnet';

// Full-stream record length (rawdump/raw_waveform nsamples for opch 120-159).
local FULLSTREAM_SAMPLES = 343808;

function(input_file, output_dir='.', run=27980, event=8, offset_us=0, fixed_snr=-1, spe_file='')

  local run_n = if std.type(run) == 'string' then std.parseInt(run) else run;
  local evt_n = if std.type(event) == 'string' then std.parseInt(event) else event;
  local off_us = if std.type(offset_us) == 'string' then std.parseJson(offset_us) else offset_us;
  local snr = if std.type(fixed_snr) == 'string' then std.parseJson(fixed_snr) else fixed_snr;
  // Single node instances reused below so all references share identity.
  // fixed_snr<=0 keeps flash.jsonnet's PDHD default (0.005); pass >0 to sweep.
  // spe_file='' keeps the default 2024 averages; pass pdhd-spe-templates-tuned.json
  // for the per-channel tuned FBK templates (pdhd/docs/pdhd-spe-template-tuning.md).
  local source = flash.opwaveform_source(input_file, run_n, evt_n);
  local decon = if snr > 0 then flash.opdecon(samples=FULLSTREAM_SAMPLES, fixed_snr=snr, spe_file=spe_file)
                else flash.opdecon(samples=FULLSTREAM_SAMPLES, spe_file=spe_file);
  // ROI cleaning of the decon ("decon" -> "decon_roi"): high-pass to find ROIs
  // above ~5 sigma (padded), then on the ORIGINAL decon zero everything outside
  // ROIs and linear-baseline each ROI to start/end at zero; ringing channels are
  // zeroed.  This SUPERSEDES the robust_baseline DC/ringing handling below (it
  // owns both); see pdhd-fullstream-light-reco.md.
  // opch 135 and 147 are bad data-quality PDs (hand-scan of pdhd/pics/pd/wf_ch*.png;
  // see pdhd-fullstream-light-reco.md): zero them outright so they raise no OpHits,
  // independent of the per-event MAD ringing veto.
  local roi = flash.oproi(veto_channels=[135, 147]);
  // Raised hit threshold (~5 sigma of the decon noise floor): the full-stream
  // scans 5.5 ms continuously, so a snippet-mode 3.0 threshold (~1.3 sigma here)
  // would integrate noise into ~thousands of spurious flashes.
  // The OpHit finder reads the ROI-cleaned 'decon_roi' traces.  The hysteresis
  // ROIs hug each pulse, so the in-ROI samples are signal-dominated and a robust
  // median/MAD over them would close the start gate -- use the known clean noise
  // floor (fixed_ped_sigma ~ the HPF rms 0.02 decon -> 2 scaled) with ped_mean=0
  // (the ROIs are endpoint-zeroed).  OpRoi already zeroes ringing channels, so no
  // robust_baseline veto is needed.  (Snippets keep the head method, wct-light-reco.jsonnet.)
  local hit = flash.ophit(hit_threshold=11.0, intag='decon_roi', fixed_ped_sigma=2.0);
  local opflash_finder = flash.opflash_finder(offset_us=off_us);
  local wf_sink = flash.waveform_sink('%s/light-frames-fullstream-wct.tar.bz2' % output_dir,
                                      tags=['raw', 'decon', 'decon_roi'], name='fswct');
  local fl_sink = flash.opflash_sink('%s/opflash_pdhd-fullstream-wct.tar.gz' % output_dir, name='fswct');

  local fanout = g.pnode({
    type: 'FrameFanout',
    name: 'light_fullstream_reco',
    data: { multiplicity: 2 },
  }, nin=1, nout=2);

  local graph = g.intern(
    innodes=[source],
    centernodes=[decon, roi, fanout, hit, opflash_finder],
    outnodes=[wf_sink, fl_sink],
    edges=[
      g.edge(source, decon),
      g.edge(decon, roi),
      g.edge(roi, fanout),
      g.edge(fanout, wf_sink, 0, 0),
      g.edge(fanout, hit, 1, 0),
      g.edge(hit, opflash_finder),
      g.edge(opflash_finder, fl_sink),
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
