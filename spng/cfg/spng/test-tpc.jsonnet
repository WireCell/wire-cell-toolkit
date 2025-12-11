// build a graph in the scope of one tpc.

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

local fans_mod = import "fans.jsonnet";
local frame_mod = import "frame.jsonnet";
local decon_mod = import "decon.jsonnet";
local cv_mod = import "crossviews.jsonnet";
// local dnnroi_mod = import "dnnroi.jsonnet";
local controls_mod = import "control.jsonnet";

local detconf = import "detconf.jsonnet";

function(detname="pdhd", tpcid="tpc0", stage="decon_pack")
    local det = detconf[detname];
    local tpc = det.tpc[tpcid];

    local controls = controls_mod(device='gpu', verbosity=2, semaphore=1);
    local control = controls.config;

    local fans = fans_mod(control);
    local frame = frame_mod(control);
    local decon = decon_mod(control);
    local cv = cv_mod(control);
    // local dnnroi = dnnroi_mod(control);

    local stages = {
        // 1 frame -> 1 frame tensor set + ngroup tensors
        input_data: frame.set_plus_groups(tpc),

        // 1 frame tensor set -> ngroup tensors
        ts_unpack: frame.tensorset_unpacker(tpc),

        // ngroup tensors -> nview tensors
        decon_pack: decon.group_decon_view(tpc),

        // 1 frame tensor set -> nview tensors
        frame_decon: pg.shuntline($.ts_unpack, $.decon_pack),

        gauss_filter_one: decon.time_filter_one(tpc, "gauss", 0),

        local filter_names = ["gauss", "wiener", "dnnroi"],
        local filters = {
            [filter]: decon.time_filter_views(tpc, filter)
            for filter in filter_names
        },
        wiener: filters.wiener,

        frame_filters: fans.tensor_fanout_shuntline($.frame_decon, [
            filters.gauss, filters.wiener, filters.dnnroi
        ]),

        // 3->3
        threshold: cv.threshold_views(tpc),

        crossview0: cv.crossview(tpc, 0),

        // 3->3
        crossfan: cv.crossfan(tpc),

        // combine filter and threshold into single 3->3
        filter_threshold: pg.shuntline(filters.wiener, $.threshold),

        // Attach that to input of crossviews
        threshold_cv: pg.shuntline($.filter_threshold, $.crossfan),

        // dnnroi_prepare: dnnroi.prepare(tpc.name, filters.dnnroi, $.crossfan),

        // all: pg.components([$.frame_filters, $.threshold_cv, $.dnnroi_prepare]),

    };
    local graph = stages[stage];
    pg.main(graph, 'Pgrapher', plugins=["WireCellSpng"])
    // graph





        // // connect groups
        // // 1 frame -> 1frame tensor set + nview tenors
        // input_pack: pg.intern(innodes=[self.input_data],
        //                       outnodes=[self.input_data, self.decon_pack],
        //                       oports=[self.input_data.oports[0]]+self.decon_pack.oports,
        //                       edges=[
        //                           pg.edge(self.input_data, self.decon_pack, g+1, g)
        //                           for g in wc.iota(std.length(tpc.view_groups))]),

        // // short circuit convolution->gauss filter.
        // // More generally we must fanout to 3 filter types.
        // input_gauss: pg.intern(innodes=[self.input_pack],
        //                        outnodes=[self.input_pack, filters.gauss],
        //                        edges=[
        //                            pg.edge(self.input_pack, filters.gauss, v+1, v)
        //                            for v in [0,1,2]]),
        
