// This covers SPNG from ADC input to signal output.

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local io = import "spng/io.jsonnet";
local tio = import "spng/torchio.jsonnet";

local fans_mod = import "spng/fans.jsonnet";
local frame_mod = import "spng/frame.jsonnet";
local roi_mod = import "spng/roi.jsonnet";
local tpc_mod = import "spng/tpc.jsonnet";



local detconf = import "spng/detconf.jsonnet";
local detector = import "spng/detector.jsonnet";
local control_module = import "spng/control.jsonnet";



// An experiment to inject a generic wrapper.
local wrapped_pnode(class_wrappers)=
    function(inode, nin=0, nout=0, uses=[], name=null) // pnode
        //local itype=std.trace('inode.type=%s'%inode.type, inode.type);
        local itype=inode.type;
        local w = std.get(class_wrappers, itype, function(inode, pnode) pnode);
        w(inode, pg.pnode(inode, nin=nin, nout=nout, uses=uses, name=name));


/// Top-level arguments:
///
/// input :: name of input file
///
/// output :: name for output file
///
/// view_crossed :: which views should be subject to crossviews
///
/// dump :: if given, a list of tiers to dump.  See pipeline definition below for labels
function(input, output, tpcid=0, view_crossed=[1,1,0],
         dump="",
         detname='pdhd',
         // fixme: add seeds=[1,2,3,4]
         engine='Pgrapher', device='cpu', verbosity=0, semaphore=1)

    // A list of view indices for views that have crossviews calculated for DNNROI type ROIs.
    local crossed_view_indices = [x.index for x in wc.enumerate(view_crossed) if x.value == 1];
    // A list of view indices for views where ROIs are calculated with simple thresholds.
    local threshold_view_indices = [x.index for x in wc.enumerate(view_crossed) if x.value == 0];

    local tpcid_num = wc.intify(tpcid);

    local det = detector.subset(detconf[detname], [tpcid_num]);
    local tpc = det.tpcs[tpcid_num];

    local controls = control_module(device=device,
                                    verbosity=wc.intify(verbosity),
                                    semaphore=wc.intify(semaphore));
    local control = controls.config;

    local frame_nodes = frame_mod(control);

    local dump_views_maybe(tier, pnode) =
        local wrap="tensors-%(tier)s-%(tpcid)s.pkl";
        local nodump = std.length(std.findSubstr(tier, dump)) == 0;
        if nodump
        then pnode
        else
            local name = 'wrap_tpc' + std.toString(tpcid) + '_' + tier;
            local filename = wrap % {tier: tier, tpcid: std.toString(tpcid)};
            local sink = tio.pickle_tensor_set(filename);
            local tap = frame_nodes.sink_taps(name, std.length(pnode.oports), sink);
            pg.shuntline(pnode, tap);
            

    local dump_inode_maybe(tier, inode, pnode) =
        local wrap="tensors-%(tier)s-%(tpcid)s-%(itype)s-%(iname)s.pkl";
        local nodump = std.length(std.findSubstr(tier, dump)) == 0;
        if nodump
        then pnode
        else
            local name = tier + '_' + inode.type + '_' + inode.name;
            local filename = wrap % {tier: tier, tpcid: std.toString(tpcid), itype:inode.type, iname:inode.name};
            local sink = tio.pickle_tensor_set(filename);
            local tap = frame_nodes.sink_taps(name, std.length(pnode.oports), sink);
            pg.shuntline(pnode, tap);

    local pnode_wrappers = {
        SPNGCrossViews: function(inode, pnode) dump_inode_maybe("crossviews", inode, pnode),
        SPNGCrossViewsExtract: function(inode, pnode) dump_inode_maybe("mps", inode, pnode),
        SPNGTensorForward: function(inode, pnode) dump_inode_maybe("dnnroi", inode, pnode),
    };
    local wpg = pg + { pnode:: wrapped_pnode(pnode_wrappers) };
            

    local tpc_nodes = tpc_mod(tpc, control, pg=wpg);
    local roi_nodes = roi_mod(control, pg=wpg);
    local fan_nodes = fans_mod(control);

    local source = io.frame_array_source(input);
    // 1->4
    local unpack = dump_views_maybe("unpack", tpc_nodes.frame_direct_unpack);
    // 4->3
    local decon = dump_views_maybe("decon", tpc_nodes.response_decon);

    // The targets_list implicitly reflects the view_crossed array, so FIXME:
    // parameterize its construction on view_crossed.
    local decon_fans = fan_nodes.fanout_select(tpc.name+"_decon_", N=3,
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

    local dnnroi_filter = dump_views_maybe("filterdnnroi",
                                           tpc_nodes.time_filter("dnnroi", views=crossed_view_indices));

    // 2->2
    local dnnroi_filter_stage = pg.shuntlines([
        decon_fans.targets.dnnroi,
        dnnroi_filter,
    ]);

    local gauss_filter = dump_views_maybe("filtergauss", tpc_nodes.time_filter("gauss", views=[0,1,2]));

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
    local signals=dump_views_maybe("signals", pg.crossline(view_signals));

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
    

    pg.main(graph, engine, plugins=["WireCellSpng"], uses=controls.uses)


