
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

    // // For inference, sg1_infer is just the beginning.  Must feed .fodder to
    // // dnnroi forward and then that plus .rois plus .gauss into apply_roi.
    // local sg1_infer = sg.dnnroi_inference_preface(rebin=rebin, extra_name="_INFER");

    // local sg1_feed = pg.shuntline(sg0, sg1_infer.decon_sink);

    // local dnnroifwd = sg.dnnroi_forward_views(dnnroi_model_file, 
    //                                           crossed_views=crossed_views, extra_name="_FORWARD");
    // local sg1_dnnroi = pg.shuntline(sg1_infer.fodder, dnnroifwd);

    // // Merge all sources of ROIs.
    // ///
    // // BIG FAT WARNING: this should be done with knowledge of crossed_views to
    // // get the port ordering correct.  As is, it assumes dnnroi rois all come
    // // before the views that have only initial rois.
    // local fat_rois = pg.crossline([sg1_dnnroi, sg1_infer.rois]);
    // local unrebin = sg.unrebin_views(rebin=rebin, views=all_views, extra_name="_ROIS");
    // local rois = pg.shuntline(fat_rois, unrebin);

    // /// Gives .roi_sink, .dense_sink and .signal_source
    // local applyrois = sg.applyroi_views(views=all_views, extra_name="_APPLYROIS");
    // local rois_cap = pg.shuntline(rois, applyrois.roi_sink);
    // local dense_cap = pg.shuntline(sg1_infer.gauss, applyrois.dense_sink);
    // local rollup = pg.components([sg1_feed, rois_cap, dense_cap, applyrois.signal_source]);
    

    //local graph = sg1_feed;
    //local graph = dnnroifwd;
    //local graph = rollup;
    local graph = pg.pipeline([sg0, sg1_infer]);
    pg.main(graph, engine, plugins=["WireCellSpng"])
    




    // // [1]frame -> tensor[3]
    // local sg1 = sg.frame_decon(extra_name="_FRAME_DECON");
    // local sg01 = pg.pipeline([sg0, sg1]);

    // // Fanout the decon for training (3*wiener+2*dense)
    // local fan_training = sg.fanout_for_dnnroi_training(crossed_views,extra_name="FOO");
    // local sg01_connection = pg.shuntline(sg01, fan_training.sink);

    // // FIXME: we carry some schizophrenia for traning vs inference.  After
    // // working this out, bundle each branch into a subgraphs function.
    // // Until then, the extra_name is required to avoid crossing the streams.
    // // Fanout the decon for inference (3*wiener+2*dense+3*gauss)
    // // local fan_inference = sg.fanout_for_dnnroi_inference(crossed_views, "_INFERENCE");
    // // local sg1_inference = pg.shuntline(sg1, fan_inference.sink);
    // // local sg1_inference = pg.shuntline(sg1, fan_inference.sink);

    // // [3]tensor -> tensor[3]
    // local sg2 = sg.tight_roi(extra_name="_TIGHT_ROI");
    // //local sg2_connection = pg.shuntline(fan_training.targets.wiener, sg2);

    // // [3]tensor -> tensor[2]
    // local sg3 = sg.cellviews_tensors(out_views=dnnroi_views, chunk_size=0, extra_name="_CELLVIEWS");

    // // [3]tensor -> tensor[2] combo of two above
    // local sg23 = pg.shuntline(sg2, sg3);
    // local sg23_connection = pg.shuntline(fan_training.targets.wiener, sg23);

    // // [2]tensor -> tensor[2]
    // local sg4 = sg.dnnroi_dense_views(views=dnnroi_views, extra_name="_DNNROI_DENSE");
    // local sg4_connection = pg.shuntline(fan_training.targets.dense, sg4);

    // // mp_sink:[3]tensor + dense_sink:[2]tensor(decon) -> source:tensor[2]
    // local sg5 = sg.connect_dnnroi_stack(sg3, sg4, views=dnnroi_views, extra_name="_DNNROI_STACK");

    // // [3]tensor->tensor[3]
    // local sg6 = sg.gauss_dense_views(extra_name="_GAUSS");

    // local graph = pg.components([sg01_connection, sg23_connection, sg4_connection, sg5]);

