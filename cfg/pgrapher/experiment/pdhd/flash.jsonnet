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
    local dft = { type: 'FftwDFT' },
    opdecon(name='')::  g.pnode({
        type: 'OpDecon',
        name: name,
        data: {
            dft: wc.tn(dft),
            noise_file: '',
        },
    }, nin=1, nout=1, uses=[dft]),

    // SlidingWindow hit finding on "decon" traces -> ophits tensor set.
    // PDHD enables overlapping-pulse splitting: a second flash riding on
    // the first's scintillation tail is recovered as its own OpHit instead
    // of being absorbed (the component default is off / bit-identical).
    ophit(name='')::  g.pnode({
        type: 'OpHitFinder',
        name: name,
        data: {
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
    opflash_finder(name='', offset_us=0)::  g.pnode({
        type: 'OpFlashFinder',
        name: name,
        data: {
            nchan: $.nchan,
            group_by_side: true,
            offset_us: offset_us,
        },
    }, nin=1, nout=1),
}
