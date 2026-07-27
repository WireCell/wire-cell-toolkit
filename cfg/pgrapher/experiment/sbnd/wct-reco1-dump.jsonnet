// CANONICAL in-tree SBND reco1 art-file dump job — single source of truth.
// Promoted from sbnd_xin 2026-07-27 (sbnd_xin/docs/64_cfg-sync.md);
// sbnd_xin/wct-reco1-dump.jsonnet is now a one-line re-export of this file.
// Runnable directly:
//   wire-cell -c pgrapher/experiment/sbnd/wct-reco1-dump.jsonnet --tla-...
// (WIRECELL_PATH must contain toolkit/cfg + wire-cell-data).
//
// Dump an SBND reco1 art/LArSoft ROOT file into the standalone-chain
// interchange formats -- directly with the toolkit, no LArSoft.
//
// Toolkit-only counterpart of yuhw's LArSoft dumps
// (input_files/wcls-frame-dump.jsonnet + wcls-flash-dump.jsonnet):
//   SBNDReco1FrameSource   -> FrameFileSink   frames-dnn.tar.bz2
//   SBNDReco1OpFlashSource -> TensorFileSink  opflash_apa{0,1}.tar.gz
//
// By default ALL events in the art file are streamed into the combined
// archives (yuhw sample-dir layout); pass entry >= 0 to dump one event.
//
// caf_offset_mode ('none'|'product'|'auto'|'override'): whether the
// opflash tensor-set metadata carries the per-event frame_apply_at_caf
// (ns) flash-time re-reference.  'product' reads the authoritative
// FrameShiftInfo::fFrameApplyAtCaf written by the sbndcode FrameShift
// producer (requires a *_frameshift.root input); 'auto' approximates it
// from the ported FrameShift derivation (pre-FrameShift files; low by
// ~0.26 us on the dev sample); 'none' omits the key =>
// FlashTensorToOpticalPCs no-op, byte-identical semantics to the older
// 10-event dumps.  See toolkit/root/docs/sbnd-reco1-source.md.
//
// Run (from sbnd_xin/): see run_reco1_dump.sh, or directly:
//   wire-cell -l stderr -L info \
//     --tla-str input=<reco1.root> --tla-str output_dir=<dir> \
//     -c wct-reco1-dump.jsonnet

function(input, output_dir='.', entry='-1', caf_offset_mode='none', caf_offset_override='0')
// caf_offset_mode: none | product | auto | override (validated in C++)

local g = import 'pgraph.jsonnet';

local entry_num = std.parseInt(entry);
local caf_override_ns = std.parseJson(caf_offset_override);

local frame_src = g.pnode({
    type: 'SBNDReco1FrameSource',
    name: 'reco1frames',
    data: {
        filename: input,
        entry: entry_num,
        // wire/badmask/summary products, tag, scales, nticks, tick:
        // C++ defaults match the SBND data reco1 files (sptpc2d, dnnsp).
    },
}, nin=0, nout=1);

local frame_sink = g.pnode({
    type: 'FrameFileSink',
    name: 'reco1frames',
    data: {
        outname: output_dir + '/frames-dnn.tar.bz2',
        tags: ['dnnsp'],
        digitize: false,  // float32 output
        masks: true,      // include chanmask_bad_<ident>.npy
    },
}, nin=1, nout=0);

local flash_srcs = [
    g.pnode({
        type: 'SBNDReco1OpFlashSource',
        name: 'tpc%d' % n,
        data: {
            filename: input,
            entry: entry_num,
            flash_product: 'recob::OpFlashs_opflashtpc%d__Reco1.' % n,
            caf_offset_mode: caf_offset_mode,
            // Only meaningful in override mode; C++ default 0.  Key omitted
            // otherwise => compiled config byte-identical to pre-override runs.
            [if caf_offset_mode == 'override' then 'caf_offset_override']: caf_override_ns,
        },
    }, nin=0, nout=1)
    for n in [0, 1]
];

local flash_sinks = [
    g.pnode({
        type: 'TensorFileSink',
        name: 'opflash_sink_apa%d' % n,
        data: {
            outname: output_dir + '/opflash_apa%d.tar.gz' % n,
            prefix: 'opflash_',
        },
    }, nin=1, nout=0)
    for n in [0, 1]
];

local pipes = [g.pipeline([frame_src, frame_sink], 'frames')] +
              [g.pipeline([flash_srcs[n], flash_sinks[n]], 'flash%d' % n) for n in [0, 1]];

local app = {
    type: 'Pgrapher',
    data: {
        edges: std.foldl(function(acc, p) acc + g.edges(p), pipes, []),
    },
};

local cmdline = {
    type: 'wire-cell',
    data: {
        plugins: ['WireCellRoot', 'WireCellSio', 'WireCellPgraph', 'WireCellAux'],
        apps: ['Pgrapher'],
    },
};

[cmdline] + std.foldl(function(acc, p) acc + g.uses(p), pipes, []) + [app]
