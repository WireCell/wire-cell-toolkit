// PDVD light (photon detector) helper: builders for the light-data
// I/O and reconstruction graph nodes.  Mirrors
// pgrapher/experiment/pdhd/flash.jsonnet with the PDVD conventions:
//
//  - 40 OpDets, 51 live DAPHNE OpChannels (block codes 10xx cathode
//    X-ARAPUCA / 20xx membrane XA / 30xx PMT).  The X-ARAPUCAs have
//    TWO channels per OpDet, so traces/ophits carry the DAPHNE
//    OpChannel and OpFlashFinder gangs them to OpDet PE columns via
//    pdvd-opch-map.json.
//  - positive raw polarity (input_polarity +1), 16 ns ticks.
//  - cathode PDs are continuous full streams (468800/468864 samples,
//    deconvolved at samples=468864 with benign zero-padding);
//    membrane XA + PMT are 1024-tick self-trigger snippets.
//  - deconvolution uses the WIENER-INSPIRED filter (F(0)=1, PE
//    normalization exact, AutoScale skipped): sigma 1.25 MHz cathode /
//    1.0 membrane / 3.5 PMT, power 2 -- pdvd/docs/pdvd-light-filter.md.
//  - single all-PD flash, bin_width 1 us (pdvd/docs/pdvd-flash-dt.md);
//    group_by_side/flash_refine off (their x>=0 split and y/z grid are
//    PDHD horizontal-drift assumptions).

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

{
    nchan:: 40,

    spe_file:: 'pgrapher/experiment/protodunevd/pdvd-spe-templates.json',
    geom_file:: 'pgrapher/experiment/protodunevd/pdvd-opdet-geom.json',
    opch_map_file:: 'pgrapher/experiment/protodunevd/pdvd-opch-map.json',

    // raw_waveform TTree -> IFrame with "raw" tagged traces (channel =
    // DAPHNE OpChannel; tick origin = the event's earliest record over
    // ALL channels, shared by every per-population instance).
    // opch_lo/hi select the population: cathode 1000-1999, membrane
    // 2000-2999, PMT 3000-3999; -1 = unbounded.
    opwaveform_source(filename, run, event, opch_lo=-1, opch_hi=-1, name='')::  g.pnode({
        type: 'PDVDOpWaveformSource',
        name: name,
        data: {
            filename: filename,
            run: run,
            event: event,
            opch_lo: opch_lo,
            opch_hi: opch_hi,
        },
    }, nin=0, nout=1),

    // Opflash tensor-set archive writer (schema of flash/docs/design.md
    // §3.4 with nchan=40; consumed later by TensorFileSource +
    // FlashTensorToOpticalPCs{nchan: 40} + QLMatching{nchan: 40}).
    opflash_sink(outname, name='')::  g.pnode({
        type: 'TensorFileSink',
        name: name,
        data: {
            outname: outname,
            prefix: 'opflash_',
        },
    }, nin=1, nout=0),

    // Waveform frame dump (one dense array per tag).
    waveform_sink(outname, tags=['raw', 'decon'], name='')::  g.pnode({
        type: 'FrameFileSink',
        name: name,
        data: {
            outname: outname,
            tags: tags,
            digitize: false,
        },
    }, nin=1, nout=0),

    local dft = { type: 'FftwDFT' },

    // Wiener-inspired deconvolution of "raw" traces -> "decon" traces.
    // G = conj(H) F(f) / (|H|^2 + eps), F = exp(-0.5 (f/sigma)^2) with
    // F(0)=1 exactly: a 1-PE pulse deconvolves to unit area for any
    // sigma, so AutoScale is a no-op and the PE scale cannot move with
    // filter tuning.  The Gauss post-filter is disabled (F *is* the
    // band filter).  Per-population sigma from the real-noise scan
    // (pdvd-light-filter.md): cathode 1.25, membrane XA 1.0 (one branch
    // top+bottom), PMT 3.5 MHz.
    // samples: 1024 for snippets; 468864 for the cathode full streams
    // (the 468800-sample records are zero-padded by OpDecon).
    // DAPHNE 14-bit rail saturation flags as in PDHD.
    opdecon(name='', samples=1024, wi_sigma=1.0, detect_saturation=false, saturation_pad=0,
            saturation_repair=false)::  g.pnode({
        type: 'OpDecon',
        name: name,
        data: {
            dft: wc.tn(dft),
            spe_file: $.spe_file,
            noise_file: '',
            samples: samples,
            input_polarity: 1,
            wiener_inspired: true,
            wi_sigma_mhz: wi_sigma,
            wi_power: 2.0,
            apply_postfilter: false,
            // Perf knob (validated 2026-07-08, flash/docs/light-perf-round1.md):
            // true real FFTs -- same physics to float round-off, half the FFT work.
            use_real_dft: true,
            [if detect_saturation then 'detect_saturation']: detect_saturation,
            [if detect_saturation && saturation_pad != 0 then 'saturation_pad']: saturation_pad,
            // Two-sided exp bridge over railed runs before decon.  C++ default
            // false.  Key omitted when off => byte-identical pre-fix config.
            // See pdvd/docs/qlmatch/pdvd-saturation-recovery.md.
            [if saturation_repair then 'saturation_repair']: true,
        },
    }, nin=1, nout=1, uses=[dft]),

    // ROI identification + cleaning of the full-stream cathode decon
    // ("decon" -> "decon_roi"); same algorithm and prototype-tuned
    // defaults as PDHD (see pdhd/flash.jsonnet for the full rationale).
    // veto_sigma is exposed: it is a cap on the per-channel decon MAD
    // and must track the PDVD cathode noise floor, not PDHD's.
    oproi(name='', veto_sigma=0.1, veto_channels=[])::  g.pnode({
        type: 'OpRoi',
        name: name,
        data: {
            dft: wc.tn(dft),
            hpf_tau_mhz: 0.05,
            roi_seed_nsigma: 5.0,
            roi_ext_nsigma: 3.0,
            roi_pad_pre: 50,
            roi_post_peak: 300,
            veto_sigma: veto_sigma,
            veto_channels: veto_channels,
            use_real_dft: true,  // perf knob, round-off-level change only
        },
    }, nin=1, nout=1, uses=[dft]),

    // SlidingWindow hit finding on "decon" traces -> ophits tensor set.
    // Same split-pulse settings as PDHD.  hit_threshold in scaled decon
    // units (100 = 1 PE/tick); starting values from the WI noise floors
    // (pdvd-light-filter.md): cathode ~3.7, membrane ~4 (top-wall 5
    // sigma), PMT ~2.2.
    ophit(name='', hit_threshold=3.0, robust_baseline=false, intag='decon', fixed_ped_sigma=0, veto_saturation=false, flag_saturation=false, emit_coverage=false)::  g.pnode({
        type: 'OpHitFinder',
        name: name,
        data: {
            hit_threshold: hit_threshold,
            robust_baseline: robust_baseline,
            [if intag != 'decon' then 'intag']: intag,
            [if fixed_ped_sigma > 0 then 'fixed_ped_sigma']: fixed_ped_sigma,
            [if veto_saturation then 'veto_saturation']: veto_saturation,
            // Keep-and-mark alternative to the veto: 10th ophit column ->
            // OpFlashFinder flash_sat tensor -> QLMatching per-flash channel
            // mask.  C++ default false.  Key omitted when off => byte-identical.
            [if flag_saturation then 'flag_saturation']: true,
            // Per-trace livetime rows -> OpFlashFinder flash_cov tensor ->
            // QLMatching per-flash no-data mask (membrane XA / PMT channels
            // are 16.4-us self-trigger snippets; without this an uncovered
            // channel is scored measured = 0).  C++ default false.  Key
            // omitted when off => byte-identical.
            [if emit_coverage then 'emit_coverage']: true,
            algo: {
                split_enable: true,
                split_min_prominence: 0.4,
                split_min_prominence_abs: 100.0,
                split_min_peak: 3.0,
                split_min_separation: 2,
            },
        },
    }, nin=1, nout=1),

    // N-to-1 fan-in row-concatenating the "ophits" tensors of the three
    // population branches (cathode / membrane / PMT) so one
    // OpFlashFinder builds all-PD flashes.
    ophit_merge(name='', multiplicity=3, meta_port=0)::  g.pnode({
        type: 'OpHitMerge',
        name: name,
        data: {
            multiplicity: multiplicity,
            meta_port: meta_port,
        },
    }, nin=multiplicity, nout=1),

    // OpFlashAlg flash assembly -> opflash tensor set.  Single all-PD
    // flash (pdvd-flash-dt.md: all four PD groups coincide within
    // +-0.6 us, the double-sided cathode XAs see both drift volumes so
    // a per-volume split is not meaningful the PDHD way).  Hit channels
    // (DAPHNE) are ganged to the 40 OpDet PE columns via the channel
    // map.  bin_width 1000 ns (PDHD parity, validated on data).
    // Quality cuts scaled to 40 PDs (PDHD: 5/20 on 160).
    // offset_us: legacy scalar key, stays 0 for PDVD (inert).  The real
    // per-event light<->charge offsets are PER CRATE (the TDE/BDE charge
    // windows open up to ~32 us apart, each jittering vs the trigger):
    //   offset_bot_us = light_chain_t0 - charge_bde_start  (bottom volume)
    //   offset_top_us = light_chain_t0 - charge_tde_start  (top volume)
    // (negative, ~-2.1..-2.5 ms; ADD to a flash time to land on that crate's
    // charge time base).  Measured per event from the rawwf trigoff tree by
    // run_light_evt.sh and stamped verbatim into the archive metadata for
    // run_clus_evt.sh / QLMatching trigger_offsets.
    opflash_finder(name='', offset_us=0, min_fired_pds=2, min_total_pe=10.0,
                   offset_bot_us=null, offset_top_us=null)::  g.pnode({
        type: 'OpFlashFinder',
        name: name,
        data: {
            nchan: $.nchan,
            geom_file: $.geom_file,
            channel_map_file: $.opch_map_file,
            group_by_side: false,
            flash_refine: false,
            offset_us: offset_us,
            min_fired_pds: min_fired_pds,
            min_total_pe: min_total_pe,
            min_fired_pe: 1.0,
        } + if offset_bot_us == null && offset_top_us == null then {} else {
            metadata_extra: {
                offset_bot_us: offset_bot_us,
                offset_top_us: offset_top_us,
            },
        },
    }, nin=1, nout=1),
}
