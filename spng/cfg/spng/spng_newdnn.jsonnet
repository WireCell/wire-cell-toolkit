local wc = import "wirecell.jsonnet";
local real_pg = import "pgraph.jsonnet";

// Fixme: should we be passing these in like we do with "pg"?
local fans_mod = import "spng/fans.jsonnet";
local roi_mod = import "spng/roi.jsonnet";
local tpc_mod = import "spng/tpc.jsonnet";


/// Return a 1->1 node consuming ADC IFrame from one anode and producing a signal IFrame.
function(tpc, control, view_crossed=[1,1,0], pg=real_pg)

    // A list of view indices for views that have crossviews calculated for DNNROI type ROIs.
    local crossed_view_indices = [x.index for x in wc.enumerate(view_crossed) if x.value == 1];
    // A list of view indices for views where ROIs are calculated with simple thresholds.
    local threshold_view_indices = [x.index for x in wc.enumerate(view_crossed) if x.value == 0];

    local tpc_nodes = tpc_mod(tpc, control, pg=pg);
    local roi_nodes = roi_mod(control, pg=pg);
    local fan_nodes = fans_mod(control);

    // This node seems to introduce some magic but it is a mirage.  See
    // tpc.jsonnet comments for clarification.  For here, see "unpack" and
    // "repack".
    local bypass = tpc_nodes.frame_bypass;

    // 1[TensorSet]->[Tensor]4 ngroups of tensors.  Its input is implicitly
    // connected in the bypass.  Later we meet repack which is the bookend to
    // unpack.
    local unpack = tpc_nodes.frame_set_unpack;


    // // 4->3
    local decon = tpc_nodes.response_decon;

    // // The targets_list implicitly reflects the view_crossed array, so FIXME:
    // // parameterize its construction on view_crossed.
    local decon_fans = fan_nodes.fanout_select(
        tpc.name+"_decon_", N=3,
        targets_list=[
            ["dnnroi","gauss"],
            ["dnnroi","gauss"],
            ["dnnroi","gauss"]
        ]
    );
    // // 0->0 the deconvolution stage
    // local decon_stage = pg.shuntlines([
    //     unpack,
    //     decon,
    //     decon_fans.sink,        // virtual sink
    // ]);
    

    // // Next, decon_fans.targets.* provide multiview nodes for input to each
    // // filter chain.

    // // // rest are 3->3
    // // local wiener_filter = tpc_nodes.time_filter("wiener");
    // // local wiener_rebin = tpc_nodes.downsampler("wiener");
    // // local threshold = tpc_nodes.threshold;

    // // // W breaks the symmetry shared by U+V in that it feeds crossviews but does
    // // // not use crossviews to determine its own ROIs.  Instead, it gets ROI
    // // // directly from its wiener and with no downsample.  So, we must "forkout"
    // // // this edge so we can treat is special while letting it continue to be with
    // // // U/V wiener for crossviews.  This is dependent on view_crossed but here we
    // // // implicitly assume W is the special one.
    // // local ww = fan_nodes.forkout_oport(tpc.name+'ww', wiener_filter, 2);
    // // local crossview_wiener = pg.shuntline(decon_fans.targets.wiener, ww[0]);
    // // local threshold_wfull = pg.pnode({
    // //     type: 'SPNGThreshold',
    // //     name: tpc.name + 'wfull',
    // //     data: tpc.crossview_thresholds[2] + control
    // // }, nin=1, nout=1);
    // // local roi_w = pg.intern(centernodes=[ww[1]],
    // //                         outnodes=[threshold_wfull],
    // //                         edges=[ pg.edge(ww[1], threshold_wfull) ]);
                            
    // // // 3->ncrossed cross views.
    // // local crossviews = tpc_nodes.crossfan(view_crossed);

    // // // 3->2 decon->thresholded initial ROIs.
    // // local crossviews_stage = pg.shuntlines([
    // //     crossview_wiener,
    // //     wiener_rebin,
    // //     threshold,
    // //     crossviews,
    // // ]);
    // // // crossviews.oports[2] is W ROI.

    local dnnroi_filter = tpc_nodes.time_filter("dnnroi", views=[0,1,2]);

    // 2->2
    // local dnnroi_filter_stage = pg.shuntlines([
    //     decon_fans.targets.dnnroi,
    //     dnnroi_filter,
    // ]);

    local gauss_filter = tpc_nodes.time_filter("gauss", views=[0,1,2]);

    // local gauss_filter_stage = pg.intern([
    //     decon_fans.targets.gauss,
    //     gauss_filter,
    // ]);

    local group_gauss = pg.pnode({
        type: 'SPNGReduce',
        name: 'group_gauss',
        data: {
            multiplicity: 3,
            operation: "cat",
            dim: -2         // concatenate along channel dimension
        } + control
    }, nin=3, nout=1);

    // local gauss_branch = pg.intern(
    //     innodes=[gauss_filter_stage],
    //     outnodes=[group_gauss],
    //     edges = [
    //         pg.edge(gauss_filter_stage, group_gauss, 0, 0),
    //         pg.edge(gauss_filter_stage, group_gauss, 0, 1),
    //         pg.edge(gauss_filter_stage, group_gauss, 0, 2),
    //     ]
    // );

    local group_dnnroi = pg.pnode({
        type: 'SPNGReduce',
        name: 'group_dnnroi',
        data: {
            multiplicity: 3,
            operation: "cat",
            dim: -2         // concatenate along channel dimension
        } + control
    }, nin=3, nout=1);


    /// 1->1 rebin the "dense" image
    /// Forgive me Father, for I have rebinned (by 10)
    local rebin=10;
    local rebinner = pg.pnode({
        type: "SPNGRebinner",
        name: "all_rebin",
        data: {
            norm: "interpolation",
            factor: rebin,
        } + control
    }, nin=1, nout=1);

    local unbinner = pg.pnode({
        type: "SPNGRebinner",
        name: "all_unbin",
        data: {
            norm: "interpolation",
            factor: -rebin,
        } + control
    }, nin=1, nout=1);
    local this_name = 'spng';
    /// 1 -> 1.  DNNROI model requires (ntick,nchan) while WCT convention is
    /// (nchan,ntick) so we must bookend the forward with a transpose.
    local pre = pg.pnode({
        type: 'SPNGTransform',
        name: this_name+'_pre',
        data: {
            operations: [
                { operation: "transpose", dims: [-2,-1] }
            ]
        } + control
    }, nin=1, nout=1);

    /// Do inference (eg, DNNROI) or other type of "forward"
    local forward = {
        type:'SPNGTensorForwardTS',
        name:'',
        data: {
            ts_filename:"unet-l23-cosmic500-e50.ts" #FIXME
        } + control
    };
    local fwd = pg.pnode({
        type: 'SPNGTensorForward',
        name: this_name,
        data: {
            forward: wc.tn(forward),
            nbatch: 1,
        } + control
    }, nin=1, nout=1, uses=[forward]);
    local post = pg.pnode({
        type: 'SPNGTransform',
        name: this_name+'_post',
        data: {
            operations: [
                /// DNNROI traditionally works in (ntick,nchan) for last two
                /// dims while WCT works in (nchan,ntick) so we must bookend
                /// the forward with a transpose.
                { operation: "transpose", dims: [-2,-1] },
                // in principle could threshold here but see comments, next.

                // We have to force nbatch=1 on TensorForward because DNNROI
                // expects a batch dim.  But we don't actually want it on
                // output so remove it now.  Also, DNNROI keeps the "feature
                // dimension".
                { operation: "squeeze", dims: [1, 0] },
            ],
        } + control
    }, nin=1, nout=1);

    // Convert DNNROI output to Boolean mask.  Note, this can be done in a
    // Transform but we break it out specific as we may want to
    // resample prior to threshold for more precision.
    // https://github.com/WireCell/spng/issues/32
    local thresh = pg.pnode({
        type: 'SPNGThreshold',
        name: this_name+'_roi',
        data: {
            nominal: 0.5,   // FIXME: a study is needed to best set this
            tag: "roi",
            datapath_format: "/traces/Threshold/" + this_name+'_roi'
        } + control
    }, nin=1, nout=1);

    local TFTT = pg.pipeline([
        pre, fwd, post, thresh,
    ]);

    // // FIXME: implement this
    // local forward = {
    //     type:'SPNGTensorForwardTS',
    //     name:'',
    //     data: {
    //         ts_filename:"unet-l23-cosmic500-e50.ts"
    //     } + control
    // };

    // // 3 source pnodes producing "signal" tensors
    // local view_signals=[
    //     roi_nodes.connect_dnnroi(tpc.name + "v" + std.toString(vi.value), forward,
    //                             //  pg.oport_node(crossviews_stage, vi.value),
    //                              pg.oport_node(gauss_filter_stage, vi.value),
    //                              pg.oport_node(dnnroi_filter_stage, vi.index))
    //     for vi in wc.enumerate(crossed_view_indices)
    // ] + [
    //     // fixme: I don't want to hardwire "W" but rather go off threshold_view_indices.
    //     roi_nodes.connect_threshold(tpc.name + "v2", roi_w, 
    //                                 pg.oport_node(gauss_filter_stage, 2))
    // ];
    // local signals = pg.crossline(view_signals);

    // // [3]Tensor -> TensorSet[1].  This is the bookend to unpack.  Its output is
    // // implicitly connected inside the bypass.
    local repack = tpc_nodes.frame_set_repack_1;


    local mul = pg.pnode({
        type: 'SPNGReduce',
        name: this_name+"_applyroi",
        data: {
            operation: 'mul',
        } + control
    }, nin=2, nout=1);

    local rbl = pg.pnode({
        type: 'SPNGRebaseliner',
        name: this_name,
        data: {
            tag: "signal",
            datapath_format: "/traces/Rebaseliner/" + this_name
        } + control
    }, nin=1, nout=1);

    // local end_stage = pg.shuntlines([
    //     repack,
    // ]);
    //pg.components([decon_stage, crossviews_stage, end_stage])
    pg.intern(
        innodes=[tpc_nodes.frame_to_tdm],
        outnodes=[tpc_nodes.frame_from_tdm],
        centernodes=[
            decon, unpack, decon_fans.sink,
            decon_fans.targets.gauss, group_gauss,
            decon_fans.targets.dnnroi, group_dnnroi,
            rebinner,
            TFTT,
            unbinner,
            mul,
            rbl,
            repack,
        ],
        edges=[
            pg.edge(tpc_nodes.frame_to_tdm, unpack),
            pg.edge(unpack, decon, 0, 0),
            pg.edge(unpack, decon, 1, 1),
            pg.edge(unpack, decon, 2, 2),
            pg.edge(unpack, decon, 3, 3),

            pg.edge(decon, decon_fans.sink, 0, 0),
            pg.edge(decon, decon_fans.sink, 1, 1),
            pg.edge(decon, decon_fans.sink, 2, 2),

            pg.edge(decon_fans.targets.gauss, gauss_filter, 0, 0),
            pg.edge(decon_fans.targets.gauss, gauss_filter, 1, 1),
            pg.edge(decon_fans.targets.gauss, gauss_filter, 2, 2),

            pg.edge(gauss_filter, group_gauss, 0, 0),
            pg.edge(gauss_filter, group_gauss, 1, 1),
            pg.edge(gauss_filter, group_gauss, 2, 2),

            pg.edge(decon_fans.targets.dnnroi, dnnroi_filter, 0, 0),
            pg.edge(decon_fans.targets.dnnroi, dnnroi_filter, 1, 1),
            pg.edge(decon_fans.targets.dnnroi, dnnroi_filter, 2, 2),

            pg.edge(dnnroi_filter, group_dnnroi, 0, 0),
            pg.edge(dnnroi_filter, group_dnnroi, 1, 1),
            pg.edge(dnnroi_filter, group_dnnroi, 2, 2),

            pg.edge(group_dnnroi, rebinner),
            pg.edge(rebinner, TFTT),
            pg.edge(TFTT, unbinner),
            
            pg.edge(unbinner, mul, 0, 0),
            pg.edge(group_gauss, mul, 0, 1),
            
            pg.edge(mul, rbl),

            pg.edge(rbl, repack),

            pg.edge(repack, tpc_nodes.frame_from_tdm)
        ]
    )
