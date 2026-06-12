// PDHD light (photon detector) helper: builders for the light-data
// I/O and reconstruction graph nodes.  See flash/docs/design.md for
// the conventions (trigger-relative ns times, OpChannel == OpDet,
// nchan = 160) and the opflash tensor-set schema.
//
// The ROOT sources read the *temporary* light exchange format; when
// the format changes only the source builders here (and the two
// components in root/) need replacing.

local g = import 'pgraph.jsonnet';

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
}
