// WCT-native PDHD light reconstruction over ALL 160 PDs in ONE processing:
// the self-trigger snippet PDs (opch 0-119) and the -x full-stream PDs (opch
// 120-159) are reconstructed in two branches that MERGE at the OpHit level, so
// a single OpFlashFinder builds flashes spanning all 160 PDs.  In particular a
// -x flash now covers the whole 80-159 wall (both readout halves) instead of
// only one half -- the product Q/L matching consumes.
//
//   snippet branch:    PDHDOpWaveformSource(snip decoana 0-119)
//                      -> OpDecon(samples=1024) -> OpHitFinder(head method)
//   full-stream branch: PDHDOpWaveformSource(fs decoana 120-159)
//                      -> OpDecon(samples=343808) -> OpRoi(veto 135/147)
//                      -> OpHitFinder(intag=decon_roi, fixed_ped_sigma=2)
//   both ophits -> OpHitMerge -> OpFlashFinder -> opflash_pdhd-allpd-wct.tar.gz
//
// Both branches use the SAME fixed Wiener filter (fixed_snr) so the 1024-tick
// snippets and the 343808-tick stream deconvolve identically and their OpHit
// times share the trigger-relative clock (same tc_time) -- a prerequisite for
// OpFlashFinder's 1us accumulator to group them.  This reproduces, branch for
// branch, wct-light-reco.jsonnet (snippets) and wct-light-fullstream-reco.jsonnet
// (full stream); only the flash step is shared.  The two inputs are the decoana
// files produced by pdhd/pd_plot/fullstream_to_decoana.py (run_light_allpd_evt.sh).
//
// Compile/run (see run_light_allpd_evt.sh):
//   wcsonnet -A snip_file=snip.root -A fs_file=fs.root -A output_dir=... \
//            -S run=27980 -S event=8 -S offset_us=249.808 \
//            -o cfg.json wct-light-allpd-reco.jsonnet
//   wire-cell -l stderr -L debug -c cfg.json

local g = import 'pgraph.jsonnet';
local flash = import 'pgrapher/experiment/pdhd/flash.jsonnet';

// Full-stream record length (rawdump/raw_waveform nsamples for opch 120-159).
local FULLSTREAM_SAMPLES = 343808;

function(snip_file, fs_file, output_dir='.', run=27980, event=8, offset_us=0, fixed_snr=-1, spe_file='', sat_pad=1024)

  local run_n = if std.type(run) == 'string' then std.parseInt(run) else run;
  local evt_n = if std.type(event) == 'string' then std.parseInt(event) else event;
  local off_us = if std.type(offset_us) == 'string' then std.parseJson(offset_us) else offset_us;
  local snr = if std.type(fixed_snr) == 'string' then std.parseJson(fixed_snr) else fixed_snr;

  // --- snippet branch (opch 0-119), identical to wct-light-reco.jsonnet ---
  // detect_saturation/veto_saturation drop snippets whose raw ADC railed at the
  // 14-bit ceiling (16383); a clipped flat-top would otherwise over-integrate
  // and fragment into spurious flashes.  See run29107-evt1015-light-anomaly.md.
  local snip_src = flash.opwaveform_source(snip_file, run_n, evt_n, name='snip');
  local sat_pad_n = if std.type(sat_pad) == 'string' then std.parseInt(sat_pad) else sat_pad;
  local snip_decon = flash.opdecon(name='snip', detect_saturation=true, saturation_pad=sat_pad_n);
  local snip_hit = flash.ophit(name='snip', veto_saturation=true);

  // --- full-stream branch (opch 120-159), identical to wct-light-fullstream-reco.jsonnet ---
  // Saturation veto is sub-range granular (OpDecon flags contiguous railed
  // runs; OpHitFinder drops only hits overlapping them), so a clipped region in
  // the 343808-sample stream is removed without losing the rest of the channel.
  local fs_src = flash.opwaveform_source(fs_file, run_n, evt_n, name='fs');
  local fs_decon = if snr > 0 then flash.opdecon(name='fs', samples=FULLSTREAM_SAMPLES, fixed_snr=snr, spe_file=spe_file, detect_saturation=true, saturation_pad=sat_pad_n)
                   else flash.opdecon(name='fs', samples=FULLSTREAM_SAMPLES, spe_file=spe_file, detect_saturation=true, saturation_pad=sat_pad_n);
  local fs_roi = flash.oproi(name='fs', veto_channels=[135, 147]);
  local fs_hit = flash.ophit(name='fs', hit_threshold=11.0, intag='decon_roi', fixed_ped_sigma=2.0, veto_saturation=true);

  // --- merge OpHits and build flashes once over all 160 PDs ---
  local merge = flash.ophit_merge(name='allpd', multiplicity=2, meta_port=0);
  local opflash_finder = flash.opflash_finder(offset_us=off_us);
  local fl_sink = flash.opflash_sink('%s/opflash_pdhd-allpd-wct.tar.gz' % output_dir, name='allpd');

  local graph = g.intern(
    innodes=[snip_src, fs_src],
    centernodes=[snip_decon, snip_hit, fs_decon, fs_roi, fs_hit, merge, opflash_finder],
    outnodes=[fl_sink],
    edges=[
      g.edge(snip_src, snip_decon),
      g.edge(snip_decon, snip_hit),
      g.edge(fs_src, fs_decon),
      g.edge(fs_decon, fs_roi),
      g.edge(fs_roi, fs_hit),
      g.edge(snip_hit, merge, 0, 0),
      g.edge(fs_hit, merge, 0, 1),
      g.edge(merge, opflash_finder),
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
