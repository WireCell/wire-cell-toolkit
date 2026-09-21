// WCT-native PDVD light reconstruction over ALL 40 OpDets (51 DAPHNE
// channels) in ONE processing: three population branches from the same
// raw_waveform file MERGE at the OpHit level, so a single OpFlashFinder
// builds all-PD flashes (bin_width 1 us; DAPHNE channels ganged to the
// 40 OpDet PE columns via pdvd-opch-map.json).
//
//   cathode branch (opch 10xx, full stream 468800/468864):
//     PDVDOpWaveformSource -> OpDecon(WI sigma=1.25, samples=468864)
//     -> OpRoi -> OpHitFinder(intag=decon_roi, fixed_ped_sigma)
//   membrane branch (opch 20xx, 1024-tick snippets):
//     PDVDOpWaveformSource -> OpDecon(WI sigma=1.0) -> OpHitFinder
//   PMT branch (opch 30xx, 1024-tick snippets):
//     PDVDOpWaveformSource -> OpDecon(WI sigma=3.5) -> OpHitFinder
//   all three ophits -> OpHitMerge -> OpFlashFinder
//     -> opflash_pdvd-wct.tar.gz
//
// All branches use the Wiener-INSPIRED filter (F(0)=1: PE normalization
// exact and filter-independent) and share the source's common tick
// origin (the event's earliest record over ALL channels), so their
// OpHit times are directly coincident -- pdvd/docs/pdvd-light-filter.md
// and pdvd-flash-dt.md.
//
// Compile/run (see run_light_evt.sh):
//   wcsonnet -A input_file=raw.root -A output_dir=... \
//            -S run=39252 -S event=298567 \
//            -o cfg.json wct-light-reco.jsonnet
//   wire-cell -l stderr -L debug -c cfg.json

local g = import 'pgraph.jsonnet';
local flash_v1 = import 'pgrapher/experiment/protodunevd/flash.jsonnet';

// Cathode full-stream record length: 468864 = max(nsamp); the shorter
// 468800-sample records are zero-padded by OpDecon.
local FULLSTREAM_SAMPLES = 468864;

function(input_file, output_dir='.', run=39252, event=298567, offset_us=0,
         offset_bot_us=null, offset_top_us=null,
         cath_thresh=3.7, mem_thresh=4.0, pmt_thresh=2.2,
         cath_ped_sigma=0.75, sat_pad=1024,
         // Drop OpHits overlapping a DAPHNE 14-bit rail (+- sat_pad).
         // Default true = legacy behavior, compiled config byte-identical.
         // false keeps hits from railed pulses (clipped => PE underestimated
         // instead of exactly 0).  See docs/qlmatch/
         // pdvd-pd-mapping-investigation.md sec.6-7 (42% of bright flashes
         // lose >=1 cathode channel to the veto).
         veto_saturation=true,
         // Keep-and-mark chain (docs/qlmatch/pdvd-saturation-recovery.md):
         // flag_saturation appends the ophit rail flag -> OpFlashFinder emits
         // the per-flash flash_sat tensor -> QLMatching use_saturation_flag
         // masks those channels per flash.  saturation_repair additionally
         // fills railed runs with the two-sided exp bridge before decon
         // (better flash totals; the mask still applies).  Both C++-default
         // false, keys suppressed => byte-identical when off.
         flag_saturation=false, saturation_repair=false,
         // Per-trace livetime rows -> OpFlashFinder flash_cov tensor ->
         // QLMatching use_coverage_flag masks self-trigger channels with no
         // snippet over a flash's window (membrane XA / PMT, duty ~5-30%;
         // docs/qlmatch/pdvd-lightpattern-sp-investigation.md).  Single knob
         // for all three populations (partial emission would zero the
         // others' coverage).  C++ default false, key suppressed =>
         // byte-identical when off.
         emit_coverage=false,
         // Validated SPE templates v2 (pd_plot/spe_v2.py): repairs the
         // broken v1 templates (ch2020 negative area, ch1051 multi-PE mode
         // latch, ch3010/3020 contaminated tails) that corrupt those
         // channels' measured-PE scale.  false => the v1 file, byte-identical.
         spe_v2=false,
         // Remap floor-pinned OVERFLOW runs to the rail before OpDecon's rail
         // scan, so the existing detect/flag/repair chain sees the membrane
         // self-trigger snippets whose over-range pulses pin at ADC 0 rather
         // than clamping at the ceiling (invisible to detect_saturation as-is).
         // Requires detect_saturation (already on for all three branches).
         // C++ default false, key suppressed => byte-identical when off.
         // NOT a production default: the encoding mechanism is unconfirmed --
         // docs/qlmatch/14_pdvd-lightpattern-sp-investigation.md ("The
         // zero-run").
         overflow_to_rail=false,
         // Membrane wide-hit handling (docs/qlmatch/25_pdvd-wall-xa-usability.md
         // §3/§7): a slow pulse spanning a 16.4-us membrane snippet becomes ONE
         // OpHit at its PEAK time and OpFlashFinder books the whole snippet PE
         // to that peak's 1-us bin (74% of wall-XA PE lands on the wrong or no
         // flash).  'start' books the full integral at the pulse onset (the
         // cathode/total-light convention, comparable to a photon-library
         // prediction); 'slice' books per 1-us slice (time-faithful, but a
         // flash then carries only its prompt share).  C++/toolkit defaults ''
         // = off (byte-identical); STUDY knob, not a production default -- the
         // wall XAs stay masked in Q/L either way.
         mem_wide_hit_mode='', mem_wide_hit_min_width_us=2.0, mem_slice_width_us=1.0,
         // Same handling for the PMT branch (also 16.4-us self-trigger
         // snippets: 1% of PMT hits are wide but they carry 22% of the PMT
         // PE, 46% of it unassigned).  Same C++ knob, same default off.
         pmt_wide_hit_mode='',
         // Slow-tail flash merge (docs/qlmatch/26, pdvd doc 23 §7d): absorb
         // the late flash member made of wide cathode-XA slow-tail hits on
         // the seed's own lit PDs; merged flash keeps the seed (fast-peak)
         // time.  C++ default false, keys suppressed => byte-identical off.
         tail_merge=false, tail_window_us=3.0, tail_min_width_us=1.0,
         tail_pe_frac=0.7)

  local run_n = if std.type(run) == 'string' then std.parseInt(run) else run;
  local evt_n = if std.type(event) == 'string' then std.parseInt(event) else event;
  local off_us = if std.type(offset_us) == 'string' then std.parseJson(offset_us) else offset_us;
  // Per-crate light<->charge offsets from the rawwf trigoff tree (see
  // run_light_evt.sh): chain_light_t0 - charge_{bde,tde}_window_start, us.
  // Stamped verbatim into the archive metadata for downstream Q/L matching.
  local off_bot_us = if std.type(offset_bot_us) == 'string' then std.parseJson(offset_bot_us) else offset_bot_us;
  local off_top_us = if std.type(offset_top_us) == 'string' then std.parseJson(offset_top_us) else offset_top_us;
  local sat_pad_n = if std.type(sat_pad) == 'string' then std.parseInt(sat_pad) else sat_pad;
  local cath_th = if std.type(cath_thresh) == 'string' then std.parseJson(cath_thresh) else cath_thresh;
  local mem_th = if std.type(mem_thresh) == 'string' then std.parseJson(mem_thresh) else mem_thresh;
  local pmt_th = if std.type(pmt_thresh) == 'string' then std.parseJson(pmt_thresh) else pmt_thresh;
  local cath_ps = if std.type(cath_ped_sigma) == 'string' then std.parseJson(cath_ped_sigma) else cath_ped_sigma;
  local veto_sat = if std.type(veto_saturation) == 'string' then std.parseJson(veto_saturation) else veto_saturation;
  local flag_sat = if std.type(flag_saturation) == 'string' then std.parseJson(flag_saturation) else flag_saturation;
  local emit_cov = if std.type(emit_coverage) == 'string' then std.parseJson(emit_coverage) else emit_coverage;
  local spe_v2b = if std.type(spe_v2) == 'string' then std.parseJson(spe_v2) else spe_v2;
  local mem_wh_minw = if std.type(mem_wide_hit_min_width_us) == 'string' then std.parseJson(mem_wide_hit_min_width_us) else mem_wide_hit_min_width_us;
  local mem_slice_w = if std.type(mem_slice_width_us) == 'string' then std.parseJson(mem_slice_width_us) else mem_slice_width_us;
  // spe_v2 swaps the OpDecon template file for the validated v2 set; the
  // toolkit flash.jsonnet default stays the v1 file (byte-identical off).
  local flash = flash_v1 + (if spe_v2b then {
      spe_file:: 'pgrapher/experiment/protodunevd/pdvd-spe-templates-v2.json',
  } else {});
  local sat_rep = if std.type(saturation_repair) == 'string' then std.parseJson(saturation_repair) else saturation_repair;
  local ovf_rail = if std.type(overflow_to_rail) == 'string' then std.parseJson(overflow_to_rail) else overflow_to_rail;
  local tmerge = if std.type(tail_merge) == 'string' then std.parseJson(tail_merge) else tail_merge;
  local twin = if std.type(tail_window_us) == 'string' then std.parseJson(tail_window_us) else tail_window_us;
  local tminw = if std.type(tail_min_width_us) == 'string' then std.parseJson(tail_min_width_us) else tail_min_width_us;
  local tfrac = if std.type(tail_pe_frac) == 'string' then std.parseJson(tail_pe_frac) else tail_pe_frac;

  // --- cathode branch (opch 10xx, continuous full streams) ---
  local cath_src = flash.opwaveform_source(input_file, run_n, evt_n, opch_lo=1000, opch_hi=1999, name='cath');
  local cath_decon = flash.opdecon(name='cath', samples=FULLSTREAM_SAMPLES, wi_sigma=1.25,
                                   detect_saturation=true, saturation_pad=sat_pad_n,
                                   saturation_repair=sat_rep, overflow_to_rail=ovf_rail);
  local cath_roi = flash.oproi(name='cath');
  local cath_hit = flash.ophit(name='cath', hit_threshold=cath_th, intag='decon_roi',
                               fixed_ped_sigma=cath_ps, veto_saturation=veto_sat,
                               flag_saturation=flag_sat, emit_coverage=emit_cov);

  // --- membrane XA branch (opch 20xx, snippets; top+bottom walls share
  //     sigma=1.0, the top-wall pickup sets the threshold) ---
  local mem_src = flash.opwaveform_source(input_file, run_n, evt_n, opch_lo=2000, opch_hi=2999, name='mem');
  local mem_decon = flash.opdecon(name='mem', wi_sigma=1.0,
                                  detect_saturation=true, saturation_pad=sat_pad_n,
                                  saturation_repair=sat_rep, overflow_to_rail=ovf_rail);
  local mem_hit = flash.ophit(name='mem', hit_threshold=mem_th, veto_saturation=veto_sat,
                              flag_saturation=flag_sat, emit_coverage=emit_cov,
                              wide_hit_mode=mem_wide_hit_mode,
                              wide_hit_min_width_us=mem_wh_minw,
                              slice_width_us=mem_slice_w);

  // --- PMT branch (opch 30xx, snippets) ---
  local pmt_src = flash.opwaveform_source(input_file, run_n, evt_n, opch_lo=3000, opch_hi=3999, name='pmt');
  local pmt_decon = flash.opdecon(name='pmt', wi_sigma=3.5,
                                  detect_saturation=true, saturation_pad=sat_pad_n,
                                  saturation_repair=sat_rep, overflow_to_rail=ovf_rail);
  local pmt_hit = flash.ophit(name='pmt', hit_threshold=pmt_th, veto_saturation=veto_sat,
                              flag_saturation=flag_sat, emit_coverage=emit_cov,
                              wide_hit_mode=pmt_wide_hit_mode,
                              wide_hit_min_width_us=mem_wh_minw,
                              slice_width_us=mem_slice_w);

  // --- merge OpHits and build flashes once over all 40 OpDets ---
  local merge = flash.ophit_merge(name='allpd', multiplicity=3, meta_port=0);
  local opflash_finder = flash.opflash_finder(offset_us=off_us,
                                              offset_bot_us=off_bot_us,
                                              offset_top_us=off_top_us,
                                              tail_merge=tmerge,
                                              tail_window_us=twin,
                                              tail_min_width_us=tminw,
                                              tail_pe_frac=tfrac);
  local fl_sink = flash.opflash_sink('%s/opflash_pdvd-wct.tar.gz' % output_dir, name='allpd');

  local graph = g.intern(
    innodes=[cath_src, mem_src, pmt_src],
    centernodes=[cath_decon, cath_roi, cath_hit,
                 mem_decon, mem_hit, pmt_decon, pmt_hit,
                 merge, opflash_finder],
    outnodes=[fl_sink],
    edges=[
      g.edge(cath_src, cath_decon),
      g.edge(cath_decon, cath_roi),
      g.edge(cath_roi, cath_hit),
      g.edge(mem_src, mem_decon),
      g.edge(mem_decon, mem_hit),
      g.edge(pmt_src, pmt_decon),
      g.edge(pmt_decon, pmt_hit),
      g.edge(cath_hit, merge, 0, 0),
      g.edge(mem_hit, merge, 0, 1),
      g.edge(pmt_hit, merge, 0, 2),
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
