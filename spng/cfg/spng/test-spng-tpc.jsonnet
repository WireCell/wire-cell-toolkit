// This covers SPNG from ADC input to signal output.

local jobs = import "spng/jobs.jsonnet";
local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local io = import "spng/io.jsonnet";
local tio = import "spng/torchio.jsonnet";
local fans = import "spng/fans.jsonnet";
local frame = import "spng/frame.jsonnet";
local cv = import "spng/crossviews.jsonnet";
local roi = import "spng/roi.jsonnet";
local tpc_mod = import "spng/tpc.jsonnet";



local detconf = import "spng/detconf.jsonnet";
local detector = import "spng/detector.jsonnet";
local control_module = import "spng/control.jsonnet";


/// Top-level arguments:
///
/// input :: name of input file
///
/// output :: name for output file
///
/// view_crossed :: which views should be subject to crossviews
///
/// wrap :: if given, it provides a filename pattern with "%(tier)s" and "%(tpcid)d".
/// dump :: if given, a list of tiers to dump.  See pipeline definition below for labels
function(input, output, tpcid=0, view_crossed=[1,1,0],
         wrap="tensors-%(tier)s-%(tpcid)d.pkl",
         dump="",
         detname='pdhd', engine='Pgrapher', device='cpu', verbosity=0)

    // A list of view indices for views that have crossviews calculated for DNNROI type ROIs.
    local crossed_view_indices = [x.index for x in wc.enumerate(view_crossed) if x.value == 1];
    // A list of view indices for views where ROIs are calculated with simple thresholds.
    local threshold_view_indices = [x.index for x in wc.enumerate(view_crossed) if x.value == 0];

    local control = control_module.bundle(device=device, verbosity=wc.intify(verbosity));
    local det = detector.subset(detconf[detname], [wc.intify(tpcid)]);
    local tpc = det.tpcs[0];

    local dump_views_maybe(tier, node) =
        local nodump = std.length(std.findSubstr(tier, dump)) == 0;
        if nodump || wrap == ""
        then node
        else
            local name = tpc.name + '_' + tier;
            local filename = wrap % {tier: tier, tpcid: tpc.ident};
            local sink = tio.pickle_tensor_set(filename);
            local tap = frame.sink_taps(name, std.length(node.oports), sink, control=control);
            pg.shuntline(node, tap);
            
            
    local tpc_nodes = tpc_mod(tpc, control);

    local source = io.frame_array_source(input);
    // 1->4
    local unpack = dump_views_maybe("unpack", tpc_nodes.frame_direct_unpack);
    // 4->3
    local decon = dump_views_maybe("decon", tpc_nodes.response_decon);

    // The targets_list implicitly reflects the view_crossed array, so FIXME:
    // parameterize its construction on view_crossed.
    local decon_fans = fans.fanout_select(tpc.name+"_decon_", N=3,
                                          targets_list=[
                                              ["wiener","gauss","dnnroi"],
                                              ["wiener","gauss","dnnroi"],
                                              ["wiener","gauss"]]);
    // 0->0 the deconvolution stage
    local decon_stage = pg.shuntlines([
        source,
        unpack,
        decon,
        decon_fans.sink,        // virtual sink
    ]);
    

    // Next, decon_fans.targets.* provide multiview nodes for input to each
    // filter chain.

    // rest are 3->3
    local wiener_filter = dump_views_maybe("filterwiener", tpc_nodes.time_filter("wiener"));
    local wiener_rebin = dump_views_maybe("rebinwiener", tpc_nodes.downsampler("wiener"));
    local threshold = dump_views_maybe("threshold", tpc_nodes.threshold);

    // W breaks the symmetry shared by U+V in that it feeds crossviews but does
    // not use crossviews to determine its own ROIs.  Instead, it gets ROI
    // directly from its wiener and with no downsample.  So, we must "forkout"
    // this edge so we can treat is special while letting it continue to be with
    // U/V wiener for crossviews.  This is dependent on view_crossed but here we
    // implicitly assume W is the special one.
    local ww = fans.forkout_oport(tpc.name+'ww', wiener_filter, 2);
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

    local dnnroi_filter = dump_views_maybe("filterdnnroi",
                                           tpc_nodes.time_filter("dnnroi", views=crossed_view_indices));

    // 2->2
    local dnnroi_filter_stage = pg.shuntlines([
        decon_fans.targets.dnnroi,
        dnnroi_filter,
    ]);

    local gauss_filter = dump_views_maybe("filterdgauss", tpc_nodes.time_filter("gauss", views=[0,1,2]));

    local gauss_filter_stage = pg.shuntlines([
        decon_fans.targets.gauss,
        gauss_filter,
    ]);

    // FIXME: implement this
    local forward={type:'SPNGTorchScriptForward', name:'', data:{filename:"foo.ts"}};

    // 3 source pnodes producing "signal" tensors
    local view_signals=[
        roi.connect_dnnroi(tpc.name + "v" + std.toString(vi.value), forward,
                            pg.oport_node(crossviews_stage, vi.value),
                            pg.oport_node(gauss_filter_stage, vi.value),
                            pg.oport_node(dnnroi_filter_stage, vi.index))
        for vi in wc.enumerate(crossed_view_indices)
    ] + [
        // fixme: I don't want to hardwire "W" but rather go off threshold_view_indices.
        roi.connect_threshold(tpc.name + "v2", roi_w, 
                              pg.oport_node(gauss_filter_stage, 2))
    ];
    local signals=pg.crossline(view_signals);

    // until here which is 3->1
    local repack = tpc_nodes.frame_set_repack;
    local sink = tio.pickle_tensor_set(output);


    local end_stage = pg.shuntlines([
        signals,
        repack,
        sink
    ]);
    local graph = pg.components([decon_stage, crossviews_stage, end_stage]);

    // FIXME: still need to bring original frame tensor set to output to form a
    // full frame tensor set and then feed to TdmToFrame
    

    pg.main(graph, engine, plugins=["WireCellSpng"])

