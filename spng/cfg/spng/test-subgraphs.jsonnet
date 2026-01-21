
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

    local dnnroi_views = [0,1];

    local crossed_views = [1,1,0];

    // Normally we do not need extra_name, but add it here to check how it shows up.

    // [1]frame -> tensor[3]
    local sg1 = sg.frame_decon(extra_name="_FRAME_DECON");

    local fan_training = sg.fanout_for_training(crossed_views, "_TRAINING");
    local sg1_training = pg.shuntline(sg1, fan_training.sink);


    // [3]tensor -> tensor[3]
    local sg2 = sg.tight_roi(extra_name="_TIGHT_ROI");
    // [3]tensor -> tensor[2]
    local sg3 = sg.cellviews_tensors(out_views=dnnroi_views, chunk_size=0, extra_name="_CELLVIEWS");
    // [2]tensor -> tensor[2]
    local sg4 = sg.dnnroi_dense_views(views=dnnroi_views, extra_name="_DNNROI_DENSE");
    // [3]tensor(rois)+[2]tensor(decon) -> tensor(dnnroi fodder)[2]
    local sg5 = sg.connect_dnnroi_stack(sg3, sg4, views=dnnroi_views, extra_name="_DNNROI_STACK");
    // [3]tensor->tensor[3]
    local sg6 = sg.gauss_dense_views(extra_name="_GAUSS");

    //local graph = pg.shuntlines([sg1,sg2]);
    //local graph = pg.components([sg1, sg2, sg5, sg6]);
    //local graph = sg3;
    local graph = sg1_training;

    pg.main(graph, engine, plugins=["WireCellSpng"])
    // graph
