// PDVD DNN-ROI subgraph for one anode.  Adapted from spng/detconfs/pdhd/dnnroi.jsonnet.
//
// U and V planes go through DNNROIFinding; W (collection) is shunted through a
// PlaneSelector.  Output is a unified frame tagged "dnnspN" (N = anode ident).
//
// NOTE: the TorchService model path below is a placeholder copied from the PDHD
// config; PDVD has no dedicated trained model yet.  This subgraph exists mainly
// so detector.subset() can build osp_subgraphs; the check-dunereco-config SP
// comparison does not exercise the model.

local wc = import 'wirecell.jsonnet';
local pg = import 'pgraph.jsonnet';

function(tpc, prefix='dnnroi', output_scale=1.0, nticks=6000, tick_per_slice=10, nchunks=1, device='cpu')
    local ts = {
        type: 'TorchService',
        name: 'dnnroi',
        data: {
            // FIXME: placeholder — no PDVD-specific model yet.
            model: '/nfs/data/1/abashyal/spng/model_files/Pytorch-UNet/ts-model-2.3/unet-l23-cosmic500-e50.ts',
            device: device,
            concurrency: 1,
        },
    };
    local anode = tpc.anode;
    local apaid = anode.data.ident;
    local prename = prefix + std.toString(apaid);
    local intags = ['loose_lf%d' % apaid, 'mp2_roi%d' % apaid, 'mp3_roi%d' % apaid];

    local dnnroi_u = pg.pnode({
        type: 'DNNROIFinding',
        name: prename + 'u',
        data: {
            anode: wc.tn(anode),
            plane: 0,
            intags: intags,
            decon_charge_tag: 'decon_charge%d' % apaid,
            outtag: 'dnnsp%du' % apaid,
            output_scale: output_scale,
            forward: wc.tn(ts),
            tick_per_slice: tick_per_slice,
            nticks: nticks,
            nchunks: nchunks,
        },
    }, nin=1, nout=1, uses=[ts, anode]);
    local dnnroi_v = pg.pnode({
        type: 'DNNROIFinding',
        name: prename + 'v',
        data: {
            anode: wc.tn(anode),
            plane: 1,
            intags: intags,
            decon_charge_tag: 'decon_charge%d' % apaid,
            outtag: 'dnnsp%dv' % apaid,
            output_scale: output_scale,
            forward: wc.tn(ts),
            tick_per_slice: tick_per_slice,
            nticks: nticks,
            nchunks: nchunks,
        },
    }, nin=1, nout=1, uses=[ts, anode]);

    local dnnroi_w = pg.pnode({
        type: 'PlaneSelector',
        name: prename + 'w',
        data: {
            anode: wc.tn(anode),
            plane: 2,
            tags: ['gauss%d' % apaid],
            tag_rules: [{
                frame: { '.*': 'DNNROIFinding' },
                trace: { ['gauss%d' % apaid]: 'dnnsp%dw' % apaid },
            }],
        },
    }, nin=1, nout=1, uses=[anode]);

    local dnnpipes = [dnnroi_u, dnnroi_v, dnnroi_w];
    local dnnfanout = pg.pnode({
        type: 'FrameFanout',
        name: prename,
        data: { multiplicity: 3 },
    }, nin=1, nout=3);

    local dnnfanin = pg.pnode({
        type: 'FrameFanin',
        name: prename,
        data: {
            multiplicity: 3,
            tag_rules: [{
                frame: { '.*': 'dnnsp%d%s' % [apaid, plane] },
                trace: { '.*': 'dnnsp%d%s' % [apaid, plane] },
            } for plane in ['u', 'v', 'w']],
        },
    }, nin=3, nout=1);

    local retagger = pg.pnode({
        type: 'Retagger',
        name: 'dnnroi%d' % apaid,
        data: {
            tag_rules: [{
                frame: { '.*': 'dnnsp%d' % apaid },
                merge: { '.*': 'dnnsp%d' % apaid },
            }],
        },
    }, nin=1, nout=1);

    pg.intern(innodes=[dnnfanout],
              outnodes=[retagger],
              centernodes=dnnpipes + [dnnfanin],
              edges=[pg.edge(dnnfanout, dnnpipes[ind], ind, 0) for ind in [0, 1, 2]] +
                    [pg.edge(dnnpipes[ind], dnnfanin, 0, ind) for ind in [0, 1, 2]] +
                    [pg.edge(dnnfanin, retagger, 0, 0)])
