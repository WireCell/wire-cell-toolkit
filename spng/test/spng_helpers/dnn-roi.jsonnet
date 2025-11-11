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

    pg.intern(
        innodes = [dnnpreprocess],
        centernodes = [dnnprocess],
        outnodes = [dnnpostprocess],
        edges = [
            pg.edge(dnnpreprocess, dnnprocess),
            pg.edge(dnnprocess, dnnpostprocess),
        ]
    )



