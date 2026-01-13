// A (adc)->tdm->decon+roi->(tensors) to produce input ("fodder") for dnnroi

// FIXME that probably will never happen: this is a copy+paste+prune of spng.jsonnet.
// It omits the bypass and anything downstream and including the dnnroi "forward" node.

local wc = import "wirecell.jsonnet";
local real_pg = import "pgraph.jsonnet";

local frame_js = import "spng/frame.jsonnet";
local fans_js = import "spng/fans.jsonnet";
local roi_js = import "spng/roi.jsonnet";
local tpc_js = import "spng/tpc.jsonnet";


/// Return a 1->1 node consuming ADC IFrame from one anode and producing a signal IFrame.
function(tpc, control, view_crossed=[1,1,0], pg=real_pg)

    // A list of view indices for views that have crossviews calculated for DNNROI type ROIs.
    local crossed_view_indices = [x.index for x in wc.enumerate(view_crossed) if x.value == 1];
    // A list of view indices for views where ROIs are calculated with simple thresholds.
    local threshold_view_indices = [x.index for x in wc.enumerate(view_crossed) if x.value == 0];

    local tpcmod = tpc_js(tpc, control, pg=pg);
    local roimod = roi_js(control, pg=pg);
    local fanmod = fans_js(control);

    local to_tdm = frame_js(control).to_tdm(tpc);

    // 1[TensorSet]->[Tensor]4 ngroups of tensors.  Its input is implicitly
    // connected in the bypass.  Later we meet repack which is the bookend to
    // unpack.
    local unpack = tpcmod.frame_set_unpack;


    // 4->3
    local decon = tpcmod.response_decon;

    // Only need to apply the "dnnroi" filter for those that we want crossed results
    local decon_fans = fanmod.fanout_select(tpc.name+"_decon_", N=3,
                                          targets_list=[
                                              if crossed == 1
                                              then ["wiener","dnnroi"]
                                              else ["wiener"]
                                              for crossed in view_crossed]);

    // 0->0 the deconvolution stage
    local decon_stage = pg.shuntlines([
        to_tdm,
        unpack,
        decon,
        decon_fans.sink,        // virtual sink
    ]);
    

    // Next, decon_fans.targets.* provide multiview nodes for input to each
    // filter chain.

    // rest are 3->3
    local wiener_filter = tpcmod.time_filter("wiener");
    local wiener_rebin = tpcmod.downsampler("wiener");
    local threshold = tpcmod.threshold;

    local crossview_wiener = pg.shuntline(decon_fans.targets.wiener, wiener_filter);
    local threshold_wfull = pg.pnode({
        type: 'SPNGThreshold',
        name: tpc.name + 'wfull',
        data: tpc.crossview_thresholds[2] + control
    }, nin=1, nout=1);
                            
    // 3->ncrossed cross views.
    local crossviews = tpcmod.crossfan(view_crossed);

    // 3->2 decon->thresholded initial ROIs.
    local crossviews_stage = pg.shuntlines([
        crossview_wiener,
        wiener_rebin,
        threshold,
        crossviews,
    ]);
    // crossviews.oports[2] is W ROI.

    // 
    local dense_filter = tpcmod.time_filter("dnnroi", views=crossed_view_indices);

    // 2->2
    local dense_stage = pg.shuntlines([
        decon_fans.targets.dnnroi,
        dense_filter,
    ]);


    // We will only use nodes prior to actual forward


    local downstream = [
        local sg = roimod.dnnroi_subgraph(tpc.name + "v" + std.toString(vi.value), forward=null);
        local cv = pg.oport_node(crossviews_stage, vi.value);
        local df = pg.oport_node(dense_stage, vi.index);
        pg.intern(centernodes=[ cv, df, sg.crossviews, sg.dense],
                  outnodes=[sg.nodes.stack],
                  edges=[
                      pg.edge(cv, sg.crossviews),
                      pg.edge(df, sg.dense),
                  ])
        for vi in wc.enumerate(crossed_view_indices)
    ];

    local results = pg.intern(innodes=[decon_stage], outnodes=downstream);

    //pg.components([decon_stage, crossviews_stage] + connects)
    results
