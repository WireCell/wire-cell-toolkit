// This produces a function for the pre-processing, inference and post-processing of DNN-ROI

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

function(ts, plane='u', prefix="dnnroi", config)
    local dnnpreprocess = (
        pg.pnode({
            type: 'SPNGDNNROIPreProcess',
            name: 'dnnroi_preprocess_%s'%[plane],
            data:{
                plane: plane,
                nchunks: config.nchunks,
                input_scale: config.input_scale,
                input_offset: config.input_offset,
                ntick: config.nticks,
                tick_per_slice: config.tick_per_slice,
                preprocess_all: config.all_preprocess,
            },
        }, nin=1, nout=1)
    );

    local dnnprocess = (
        pg.pnode({
            type: 'SPNGDNNROIProcess',
            name: 'dnnroi_%s'%[plane],
            data: {
                plane: plane,
                forward: wc.tn(ts),
            },
        }, nin=1, nout=1, uses=[ts])
    );

    local dnnpostprocess = (
        pg.pnode({
            type: 'SPNGDNNROIPostProcess',
            name: 'dnnroi_postprocess_%s'%[plane],
            data:{
                plane: plane,
                output_scale: config.output_scale,
                nchunks: config.nchunks,
                ntick: config.nticks,
                tick_per_slice: config.tick_per_slice,
            },
        }, nin=1, nout=1)
    );

    local torch_to_tensor = (
        pg.pnode({
            type: 'TorchToTensor',
            name: 'torch_to_tensor_%s'%[plane],
            data:{},
        }, nin=1, nout=1)
    );

    pg.intern(
        innodes = [dnnpreprocess],
        centernodes = [dnnprocess] + [dnnpostprocess],
        outnodes = [torch_to_tensor],
        edges = [
            pg.edge(dnnpreprocess, dnnprocess),
            pg.edge(dnnprocess, dnnpostprocess),
            pg.edge(dnnpostprocess, torch_to_tensor)
        ]
    )



