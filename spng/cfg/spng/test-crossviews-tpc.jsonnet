// This produces a final graph to input depos and output crossview tensors

local jobs = import "spng/jobs.jsonnet";
local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local io = import "spng/io.jsonnet";
local tio = import "spng/torchio.jsonnet";
local fans = import "spng/fans.jsonnet";
local frame = import "spng/frame.jsonnet";
local cv = import "spng/crossviews.jsonnet";
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
    // rest are 3->3
    local wiener_rebin = dump_views_maybe("rebinwiener", tpc_nodes.downsampler("wiener"));
    local wiener_filter = dump_views_maybe("filterwiener", tpc_nodes.time_filter("wiener"));
    local threshold = dump_views_maybe("threshold", tpc_nodes.threshold);

    // until here which is 3->1
    local repack = tpc_nodes.frame_set_repack;
    local sink = tio.pickle_tensor_set(output);


    local decon_fans = fans.fanout_select(tpc.name+"_decon_", N=3,
                                          targets_list=[
                                              ["wiener","gauss","dnnroi"],
                                              ["wiener","gauss","dnnroi"],
                                              ["wiener","gauss"]]);

    // 0->3 produce per-view decond
    local decon_stage = pg.shuntlines([
        source,
        unpack,
        decon,
        decon_fans.sink,
    ]);
    
    // 3->3 partial cross views.
    local crossviews = tpc_nodes.crossfan(view_crossed);

    // 3->3 decon->thresholded initial ROIs.
    local crossviews_stage = pg.shuntlines([
        decon_fans.targets.wiener,
        wiener_filter,
        wiener_rebin,
        threshold,
        crossviews,
    ]);
    // crossviews.oports[2] is W ROI.


    local dnnroi_rebin = dump_views_maybe("rebindnnroi", tpc_nodes.downsampler("ddnnroi", views=[0,1]));
    local dnnroi_filter = dump_views_maybe("filterdnnroi", tpc_nodes.time_filter("dnnroi", views=[0,1]));

    // 2->2
    local dnnroi_filter_stage = pg.shuntlines([
        decon_fans.targets.dnnroi,
        dnnroi_filter,
        dnnroi_rebin,
    ]);

    local dnnroi_subnet(view_layer) =
        local exstack = cv.extract_stack_one(tpc, view_layer);
        pg.intern(
            centernodes=[exstack],
            edges=[
                pg.edge(dnnroi_rebin, exstack, view_layer, 0),
                pg.edge(crossviews, exstack, view_layer, 1)
            ]);


    local dnnrois = [
        dnnroi_subnet(0), 
        dnnroi_subnet(1), 
    ];

    local gauss_filter = dump_views_maybe("filterdgauss", tpc_nodes.time_filter("gauss", views=[0,1,2]));

    local gauss_filter_stage = pg.shuntlines([
        decon_fans.targets.gauss,
        gauss_filter,
    ]);



    // local end_stage = pg.shuntlines([
    //     exconns,
    //     repack,
    //     sink
    // ]);
    // local graph = pg.components([decon_stage, crossviews_stage, dnnroi_filter_stage, exconns, end_stage]);
    local graph = pg.components([decon_stage, crossviews_stage, dnnroi_filter_stage, gauss_filter_stage]+dnnrois);
    // local graph = pg.components([dnnroi_filter_stage] + dnnrois);

    pg.main(graph, engine, plugins=["WireCellSpng"])

