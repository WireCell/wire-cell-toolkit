// PDHD light (photon detector) helper: builders for the light-data
// I/O and reconstruction graph nodes.  See flash/docs/design.md for
// the conventions (trigger-relative ns times, OpChannel == OpDet,
// nchan = 160) and the opflash tensor-set schema.
//
// The ROOT sources read the *temporary* light exchange format; when
// the format changes only the source builders here (and the two
// components in root/) need replacing.

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

{
    // PDHD optical channel count: 160 OpDets (4 APAs x 10 bars x 4
    // windows); offline OpChannel == OpDet.
    nchan:: 160,

    // LArSoft-reconstructed flashes (opflashana trees) -> opflash tensor set.
    opflash_source(filename, run, event, name='')::  g.pnode({
        type: 'PDHDOpFlashSource',
        name: name,
        data: {
            filename: filename,
            run: run,
            event: event,
            nchan: $.nchan,
        },
    }, nin=0, nout=1),

    // decoana waveform snippets -> IFrame with "raw"/"deconv" tagged traces.
    opwaveform_source(filename, run, event, name='')::  g.pnode({
        type: 'PDHDOpWaveformSource',
        name: name,
        data: {
            filename: filename,
            run: run,
            event: event,
        },
    }, nin=0, nout=1),

    // Opflash tensor-set archive writer (same convention as the SBND
    // opflash_apa<n>.tar.gz consumed by TensorFileSource +
    // FlashTensorToOpticalPCs{nchan: 160} + QLMatching{nchan: 160}).
    opflash_sink(outname, name='')::  g.pnode({
        type: 'TensorFileSink',
        name: name,
        data: {
            outname: outname,
            prefix: 'opflash_',
        },
    }, nin=1, nout=0),

    // Waveform frame dump (one dense array per tag).
    waveform_sink(outname, tags=['raw', 'deconv'], name='')::  g.pnode({
        type: 'FrameFileSink',
        name: name,
        data: {
            outname: outname,
            tags: tags,
            digitize: false,
        },
    }, nin=1, nout=0),

    // --- WCT-native reconstruction (DUNE-method port, see
    // flash/docs/stage2-reconstruction.md for fcl correspondence) ---

    // Wiener deconvolution of "raw" snippets -> "decon" traces.
    // Defaults follow protodunehd_deconvolution.  The SPE templates in
    // pdhd-spe-templates.json default to the 2024 NP04 FBK/HPK *average*
    // templates: they are DC-balanced and deconvolve to a flat tail,
    // whereas the run28368 v1 per-channel templates (kept in
    // pdhd-spe-templates-v1.json) over-subtract the slow tail below zero
    // (see pdhd/pics/pd/README.md).  Flat Wiener N^2 (noise_file empty);
    // the run27950 per-channel spectra in pdhd-noise-templates.json are a
    // second-order effect and slightly less flat with these templates.
    // fixed_snr: PDHD uses a FIXED Wiener filter (S2/N^2 = fixed_snr) instead
    // of the component-default adaptive S2 (taken from each snippet's own peak,
    // hence signal- and record-length dependent).  Fixing the S2/N^2 ratio makes
    // the filter independent of pulse amplitude AND record length, so the
    // 1024-tick self-trigger snippets (ch 0-119) and the 343808-tick full-stream
    // PDs (ch 120-159) are deconvolved with the *same* filter -- a prerequisite
    // for comparing their OpHits.  0.005 ~ a 20:1-amplitude reference pulse on the
    // FBK PDs at the 1024-tick reference (see pdhd-fullstream-light-reco.md).
    // samples: FFT/record length; 1024 for snippets, 343808 for the full stream.
    fixed_snr:: 0.005,
    local dft = { type: 'FftwDFT' },
    // spe_file: override the SPE-template file.  Default '' keeps the OpDecon.h
    // default (pdhd-spe-templates.json, 2024 FBK/HPK averages) -> existing
    // configs stay bit-identical.  Pass pdhd-spe-templates-tuned.json to use the
    // per-channel TUNED templates for the full-stream FBK channels whose average
    // template over-subtracts the post-pulse tail (pdhd/docs/pdhd-spe-template-tuning.md).
    opdecon(name='', samples=1024, fixed_snr=$.fixed_snr, spe_file='')::  g.pnode({
        type: 'OpDecon',
        name: name,
        data: {
            dft: wc.tn(dft),
            noise_file: '',
            samples: samples,
            fixed_snr: fixed_snr,
            [if spe_file != '' then 'spe_file']: spe_file,
        },
    }, nin=1, nout=1, uses=[dft]),

    // ROI identification + ROI-based cleaning of the full-stream decon
    // ("decon" -> "decon_roi"), mirroring the charge SP induction ROI chain.
    // A high-pass filter H(f) = 1 - exp(-(f/tau)^2) (tau ~ 1/20us, scintillation
    // is < 20us) gives a ROI-finding wave; ROIs are HYSTERESIS runs -- contiguous
    // above roi_ext_nsigma*MAD (low) that reach roi_seed_nsigma*MAD (high, ~ the
    // hit threshold so a real pulse is required).  Each ROI is padded roi_pad_pre
    // before the run and extended to at least roi_post_peak ticks past the decon
    // pulse peak (the late-light window ~3x the 1.6us LAr tau); a brighter tail
    // above ext runs further on its own, and overlapping windows merge.
    // The ROIs are applied to the ORIGINAL decon: everything outside is zeroed and
    // each ROI gets a linear endpoint-zeroing baseline.  Ringing channels (MAD >
    // veto_sigma) are zeroed entirely.  Hysteresis (vs a single threshold + fixed
    // blanket pad) keeps a bright pulse's ROI from ballooning to 60-75us and
    // ramping its baseline away the pulse PE.
    // Only the full-stream path uses this (between OpDecon and OpHitFinder, with
    // ophit(intag='decon_roi', fixed_ped_sigma=2.0)); the snippet path and every
    // other config never instantiate it.  Defaults are the prototype-tuned values
    // (pdhd/pd_plot/fullstream_roi_proto.py, pdhd/docs/pdhd-fullstream-light-reco.md).
    // OpRoi builds its high-pass spectrum from the actual trace length at
    // runtime, so (unlike opdecon) it needs no 'samples' parameter.
    // veto_channels: hard per-channel veto (OpChannel ids) zeroed unconditionally,
    // for known-bad data-quality channels a hand scan flags but whose MAD does not
    // reliably clear veto_sigma.  Default [] -> bit-identical (only the MAD veto).
    oproi(name='', veto_channels=[])::  g.pnode({
        type: 'OpRoi',
        name: name,
        data: {
            dft: wc.tn(dft),
            hpf_tau_mhz: 0.05,
            roi_seed_nsigma: 5.0,
            roi_ext_nsigma: 3.0,
            roi_pad_pre: 50,
            roi_post_peak: 300,
            veto_sigma: 0.1,
            veto_channels: veto_channels,
        },
    }, nin=1, nout=1, uses=[dft]),

    // SlidingWindow hit finding on "decon" traces -> ophits tensor set.
    // PDHD enables overlapping-pulse splitting: a second flash riding on
    // the first's scintillation tail is recovered as its own OpHit instead
    // of being absorbed (the component default is off / bit-identical).
    // hit_threshold: keep pulses whose deconvolved peak >= this (scaled decon
    // units, scale=100).  Default 3.0 for the triggered self-trigger snippets,
    // where every 1024-tick window is placed on a real pulse.  The full-stream
    // path scans the whole 5.5 ms continuously, so it raises this to ~5 sigma of
    // the decon noise floor (~11, set by the noisier FBK PDs) to reject sub-PE
    // noise excursions -- ~5 sigma also ~ a 1-PE peak (see
    // pdhd-fullstream-light-reco.md).
    // robust_baseline: per-channel median/MAD pedestal for the CONTINUOUS full
    // stream, where the OpHitFinder head-pedestal (a few leading samples) cannot
    // track a per-channel DC offset or a ringing channel, so two bad channels
    // over-produce flashes (opch 121 DC-offset, opch 147 ringing; see
    // pdhd-fullstream-light-reco.md s6).  Default false -> head method,
    // bit-identical to the self-trigger snippet path and every existing config.
    // The full-stream chain turns it on; it removes the DC offset (median ped)
    // and vetoes ringing channels (MAD above robust_veto_sigma), cutting
    // full-stream OpHits ~60% (evt 8) while retaining clean-channel real hits.
    // intag: input trace tag (default 'decon'; the full-stream ROI path passes
    // 'decon_roi' to read the OpRoi-cleaned traces).  fixed_ped_sigma: the known
    // clean noise floor (scaled units) used for ROI-cleaned input -- the OpRoi
    // hysteresis ROIs hug each pulse, so the in-ROI samples are signal-dominated
    // and a median/MAD over them would close the start gate; > 0 sets ped_mean=0
    // (ROIs are endpoint-zeroed) and ped_sigma to this (start gate robust_nsigma *
    // it).  Ringing channels are already zeroed by OpRoi, so robust_baseline is
    // not needed alongside it.  Conditional keys keep the snippet path and every
    // existing config byte-identical.
    ophit(name='', hit_threshold=3.0, robust_baseline=false, intag='decon', fixed_ped_sigma=0)::  g.pnode({
        type: 'OpHitFinder',
        name: name,
        data: {
            hit_threshold: hit_threshold,
            robust_baseline: robust_baseline,
            [if intag != 'decon' then 'intag']: intag,
            [if fixed_ped_sigma > 0 then 'fixed_ped_sigma']: fixed_ped_sigma,
            algo: {
                split_enable: true,
                split_min_prominence: 0.4,
                // absolute valley-depth floor (scaled decon units; 100 = 1
                // PE/tick) so genuine second pulses split out but ripples on
                // the slow scintillation tail do not over-fragment.
                split_min_prominence_abs: 100.0,
                split_min_peak: 3.0,
                split_min_separation: 2,
            },
        },
    }, nin=1, nout=1),

    // OpFlashAlg flash assembly -> opflash tensor set (design.md §3.4).
    // group_by_side builds flashes per drift volume (the opaque cathode
    // makes the two volumes optically independent); the component default
    // is off / all-OpDet.  On the 2024 single-side readout this is
    // bit-identical to all-TPC, but it is the physically correct grouping
    // once both volumes are instrumented.
    // offset_us: per-event readout-vs-trigger offset (tc - rd_timestamp from the
    // ROOT trigoff tree, ~250) stamped verbatim into the opflash metadata so the
    // downstream charge clustering / Q-L matching can place charge on the trigger
    // time base.  Default 0 (no offset).  run_light_evt.sh supplies the real value.
    // flash_refine: merge a later, dim, few-PD flash into an earlier one whose
    // lit OpDets are physically adjacent -- the same flash over-split by the
    // per-channel OpHit splitter + 1us accumulators into a bright primary plus
    // small satellites on its scintillation tail.  Component default off
    // (bit-identical to larana); on for PDHD with the run-27305 data-tuned cuts
    // (refine_pe_ratio / refine_max_fired; see pdhd-light-raw-data.md §4).
    opflash_finder(name='', offset_us=0)::  g.pnode({
        type: 'OpFlashFinder',
        name: name,
        data: {
            nchan: $.nchan,
            group_by_side: true,
            offset_us: offset_us,
            flash_refine: true,
            refine_window_us: 8.0,
            refine_pe_ratio: 0.5,
            refine_max_fired: 2,
            refine_fired_pe: 0.5,
            // subset escape: also merge a later dim flash that lights only PDs
            // the parent already lights (a tail of a bright extended flash whose
            // light spreads over >max_fired of the SAME OpDets; e.g. evt150
            // 818us 10-PD parent + its 6-/5-PD tail fragments). Additive to the
            // few-PD merges above; see pdhd-light-raw-data.md §4.4.
            refine_subset_merge: true,
        },
    }, nin=1, nout=1),
}
