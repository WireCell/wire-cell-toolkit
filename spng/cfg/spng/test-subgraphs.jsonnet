
// This produces a final graph to input depos and output crossview tensors

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

local detconf = import "spng/detconf.jsonnet";
local detector = import "spng/detector.jsonnet";
local control_js = import "spng/control.jsonnet";

local sg_js = import "spng/subgraphs.jsonnet";

function(detname='pdhd', tpcids=[], engine='Pgrapher', device='cpu', verbosity=0)
    
    local controls = control_js(device=device, verbosity=wc.intify(verbosity));
    local det = detector.subset(detconf[detname], tpcids);
    local sg = sg_js(det.tpcs[0], controls.config, pg);

    // True if care about cross view info for the view.
    local crossed_views = [1,1,0];
    local all_views = wc.iota(std.length(crossed_views));
    local rebin = 4;

    local dnnroi_model_file = "unet-l23-cosmic500-e50.ts";

    local sg0 = sg.frame_to_tdm(extra_name="_TOTDM");

    // For a training job we just need sg0+sg1_train
    // local sg1_train = sg.dnnroi_training_preface(rebin=rebin, extra_name="_TRAIN");

    local sg1_infer = sg.dnnroi_inference(modelfile=dnnroi_model_file,
                                          rebin=rebin, crossed_views=crossed_views);

    local sg1_sink = sg.attach_file_sink_views(sg1_infer, "test-subgraphs.npz",
                                               extra_name="_SINK", prefix="signals");

    local graph = pg.pipeline([sg0, sg1_sink]);
    pg.main(graph, engine, plugins=["WireCellSpng"])


