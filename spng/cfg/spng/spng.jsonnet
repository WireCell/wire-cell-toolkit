local wc = import "wirecell.jsonnet";
local real_pg = import "pgraph.jsonnet";

// Fixme: should we be passing these in like we do with "pg"?
local fans_mod = import "spng/fans.jsonnet";
local frame_mod = import "spng/frame.jsonnet";
local roi_mod = import "spng/roi.jsonnet";
local tpc_mod = import "spng/tpc.jsonnet";


/// Return a 1->1 node consuming ADC tensor from one anode and producing a signal tensor.
function(tpc, control, view_crossed=[1,1,0], pg=real_pg)

    // A list of view indices for views that have crossviews calculated for DNNROI type ROIs.
    local crossed_view_indices = [x.index for x in wc.enumerate(view_crossed) if x.value == 1];
    // A list of view indices for views where ROIs are calculated with simple thresholds.
    local threshold_view_indices = [x.index for x in wc.enumerate(view_crossed) if x.value == 0];


    local frame_nodes = frame_mod(control);

    local tpc_nodes = tpc_mod(tpc, control, pg=pg);
    local roi_nodes = roi_mod(control, pg=pg);
    local fan_nodes = fans_mod(control);

    // 1->4
    local unpack = tpc_nodes.frame_direct_unpack;
    // 4->3
    local decon = tpc_nodes.response_decon;

    // The targets_list implicitly reflects the view_crossed array, so FIXME:
    // parameterize its construction on view_crossed.
    local decon_fans = fan_nodes.fanout_select(tpc.name+"_decon_", N=3,
                                          targets_list=[
                                              ["wiener","gauss","dnnroi"],
                                              ["wiener","gauss","dnnroi"],
                                              ["wiener","gauss"]]);
    // 0->0 the deconvolution stage
    local decon_stage = pg.shuntlines([
        unpack,
        decon,
        decon_fans.sink,        // virtual sink
    ]);
    

    // Next, decon_fans.targets.* provide multiview nodes for input to each
    // filter chain.

    // rest are 3->3
    local wiener_filter = tpc_nodes.time_filter("wiener");
    local wiener_rebin = tpc_nodes.downsampler("wiener");
    local threshold = tpc_nodes.threshold;

    // W breaks the symmetry shared by U+V in that it feeds crossviews but does
    // not use crossviews to determine its own ROIs.  Instead, it gets ROI
    // directly from its wiener and with no downsample.  So, we must "forkout"
    // this edge so we can treat is special while letting it continue to be with
    // U/V wiener for crossviews.  This is dependent on view_crossed but here we
    // implicitly assume W is the special one.
    local ww = fan_nodes.forkout_oport(tpc.name+'ww', wiener_filter, 2);
    local crossview_wiener = pg.shuntline(decon_fans.targets.wiener, ww[0]);
    local threshold_wfull = pg.pnode({
        type: 'SPNGThreshold',
        name: tpc.name + 'wfull',
        data: tpc.crossview_thresholds[2] + control
    }, nin=1, nout=1);
    local roi_w = pg.intern(centernodes=[ww[1]],
                            outnodes=[threshold_wfull],
                            edges=[ pg.edge(ww[1], threshold_wfull) ]);
                            
    // 3->ncrossed cross views.
    local crossviews = tpc_nodes.crossfan(view_crossed);

    // 3->2 decon->thresholded initial ROIs.
    local crossviews_stage = pg.shuntlines([
        crossview_wiener,
        wiener_rebin,
        threshold,
        crossviews,
    ]);
    // crossviews.oports[2] is W ROI.

    local dnnroi_filter = tpc_nodes.time_filter("dnnroi", views=crossed_view_indices);

    // 2->2
    local dnnroi_filter_stage = pg.shuntlines([
        decon_fans.targets.dnnroi,
        dnnroi_filter,
    ]);

    local gauss_filter = tpc_nodes.time_filter("gauss", views=[0,1,2]);

    local gauss_filter_stage = pg.shuntlines([
        decon_fans.targets.gauss,
        gauss_filter,
    ]);

    // FIXME: implement this
    local forward = {
        type:'SPNGTensorForwardTS',
        name:'',
        data: {
            ts_filename:"unet-l23-cosmic500-e50.ts"
        } + control
    };

    // 3 source pnodes producing "signal" tensors
    local view_signals=[
        roi_nodes.connect_dnnroi(tpc.name + "v" + std.toString(vi.value), forward,
                                 pg.oport_node(crossviews_stage, vi.value),
                                 pg.oport_node(gauss_filter_stage, vi.value),
                                 pg.oport_node(dnnroi_filter_stage, vi.index))
        for vi in wc.enumerate(crossed_view_indices)
    ] + [
        // fixme: I don't want to hardwire "W" but rather go off threshold_view_indices.
        roi_nodes.connect_threshold(tpc.name + "v2", roi_w, 
                                    pg.oport_node(gauss_filter_stage, 2))
    ];
    local signals = pg.crossline(view_signals);

    // until here which is 3->1
    local repack = tpc_nodes.frame_set_repack;


    local end_stage = pg.shuntlines([
        signals,
        repack,
    ]);
    //pg.components([decon_stage, crossviews_stage, end_stage])
    pg.intern(innodes=[decon_stage], outnodes=[end_stage], centernodes=[crossviews_stage])

    


 
   
