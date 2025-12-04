// This produces a final graph to input depos and output crossview tensors

local jobs = import "spng/jobs.jsonnet";
local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local io = import "spng/io.jsonnet";
local tio = import "spng/torchio.jsonnet";
local frame = import "spng/frame.jsonnet";
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
    local rebin = dump_views_maybe("rebin", tpc_nodes.downsampler("decon"));
    local filter = dump_views_maybe("filter", tpc_nodes.time_filter("wiener"));
    local threshold = dump_views_maybe("threshold", tpc_nodes.threshold);
    //local crossviews = dump_views_maybe("crossviews", tpc_nodes.crossfan(view_crossed));
    // this is ultimate output right now, so don't offer it for dumping
    local crossviews = tpc_nodes.crossfan(view_crossed);
    // until here which is 3->1
    local repack = tpc_nodes.frame_set_repack;
    local sink = tio.pickle_tensor_set(output);
    local pipeline = [
        source,
        unpack,
        decon,
        filter,
        rebin,
        threshold,
        crossviews,
        repack,
        sink
    ];
    local graph = pg.shuntlines(pipeline);

    pg.main(graph, engine, plugins=["WireCellSpng"])
