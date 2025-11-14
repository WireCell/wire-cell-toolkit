// Flexible ways to build an SPNG related graph or final config.

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

local fans = import "fans.jsonnet";
local frame = import "frame.jsonnet";
local decon = import "decon.jsonnet";
local cv = import "crossviews.jsonnet";
local dnnroi = import "dnnroi.jsonnet";

local detconf = import "detconf.jsonnet";

/// Return an object with all possible subgraphs for one tpc
local one_tpc(tpc) = {
    tpc: tpc,

    // Node: 1->1: IFrame -> frame tensor set.  This converts an ADC IFrame into
    // corresponding TDM tensor set.
    frame_to_tdm: frame.to_tdm(tpc),

    // Node: 1->2: frame tensor set fanout eg to send one to output of graph for
    // merging with SPNG results.
    frame_set_fanout: fans.fanout(tpc.name+"frame"),

    // Node: 1->ngroup: frame tensor set -> per-group tensors.  This breaks up the
    // ADC frame into individual ADC tensors grouped by channel groups. 
    frame_set_unpack: frame.tensorset_unpacker(tpc),

    // Node: ngroup->nview: apply per-group response decon and merge any split groups
    // to produce a subgraph with 3 per-view tensors of decon signal.  Note,
    // only FR*ER response channel filter is decon'ed.  Not time filter nor "RC"
    // response decon is included.
    response_decon: decon.group_decon_view(tpc),

    // Nodes: 3->3, each.  Apply the three types of time filters intended to follow
    // response decon.  These names are canonical.  Each detector must define
    // them in terms of filter function parameter set.
    //
    // They are interpreted as:
    //
    // - wiener :: provides input to thresholding prior to cross views
    //
    // - dnnroi :: provides the non-cross-views image for DNNROI style inference.
    //
    // - gauss :: provides final signals to which ROIs are applied.
    local filter_names = ["gauss", "wiener", "dnnroi"],
    local filters = {
        [filter]: decon.time_filter_views(tpc, filter)
        for filter in filter_names
    },

    // Node: 3->3*3. Fanout each view's decon result to each filter.  This provides a
    // single Pnode which we won't use except to add as a component. We will use
    // the individual filters.XXX defined above to make downstream connections.
    decon_filters: fans.fanout_shuntline($.response_decon, [
        filters.gauss, filters.wiener, filters.dnnroi
    ]),

    // Node: 3->3.  Connect decons directly to wiener filters for running a job
    // that ends does not include dnnroi nor gauss.
    crossview_filters: pg.shuntline($.response_decon, filters.wiener),

    // Node: 3->3.  Connect wiener filter and thresholding into single 3->3.
    filter_threshold: pg.shuntline(filters.wiener, cv.threshold_views(tpc)),
    
    // Node: 3->3.  3 inputs are fanned to three crossviews each of three inputs
    // making a fully-connected net.
    crossfan: cv.crossfan(tpc),

    // Node: 3->3. Connect threshold to the cross fan.
    threshold_crossviews: pg.shuntline($.filter_threshold, $.crossfan),
    
    // Node: 1->3. A subgraph consuming IFrame and producing 3 crossview
    // tensors.  It excludes gauss and dnnroi filtered processing.
    graph_crossview: pg.components([
        $.frame_to_tdm,
        $.frame_set_unpack,
        $.response_decon,
        $.crossview_filters,
        $.filter_threshold,
        $.threshold_crossviews
    ]),

    // Node: 2x3 -> 3.  Connect the dnnroi filters and output of crossviews
    dnnroi_prepare: dnnroi.prepare(tpc.name, filters.dnnroi, $.crossfan),
    
    // Node: 3 -> 3.  Connect dnnroi pre with forward inference.
    dnnroi_forward: dnnroi.forward(tpc.name, $.dnnroi_prepare),

    // Node: 2x3 -> 3.  Apply rois.
    dnnroi_apply: dnnroi.apply(tpc.name, filters.gauss, $.dnnroi_forward),

    // Node: 3->1.  Collect view tensors into a tensor set.
    frame_set_repack: frame.tensorset_view_repacker(tpc.name),

    // Node: 2->1.  Collect original tensor set and the repack with signals.
    frame_set_fanin: fans.fanin(tpc.name+"fanin", type='TensorSet'),

    // Node: 1->1
    frame_from_tdm: frame.from_tdm(tpc),


    // Node: 2->2.  The "bypass" is a bit weird in that it seems to produce a cycle:
    // But, this is a mirage as bypass is a subgraph. See the function comments
    // for details and a drawing that shows the inner workings.
    frame_bypass: frame.connect_bypass(frame.bypass(tpc.name+"bypass"),
                                       sources=[$.frame_to_tdm, $.frame_set_repack],
                                       targets=[$.frame_set_unpack, $.frame_from_tdm]),


    all: pg.components([$.frame_filters, $.threshold_cv, $.dnnroi_prepare]),



};

function(detname="pdhd", tpcid="tpc0", stage="decon_pack", main=true)
    local det = detconf[detname];
    local tpc = det.tpc[tpcid];
    // fixme: support "alltpcs" as loop over what follows

    local stages = one_tpc(tpc);
    local graph = stages[stage];
    if main == true || main == "true"
    then pg.main(graph, 'Pgrapher', plugins=["WireCellSpng"])
    else graph





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
        
