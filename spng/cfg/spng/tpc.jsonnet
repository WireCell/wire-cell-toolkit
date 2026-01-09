// Flexible ways to build config objects for nodes or tightly interrelated
// subgraphs that depend on a TPC (anode) context object.  The subgraphs are
// Pnodes with one TPC context.  See view.jsonnet for per-view context.
//
local wc = import "wirecell.jsonnet";
local real_pg = import "pgraph.jsonnet";

local fans_mod = import "fans.jsonnet";
local frame_mod = import "frame.jsonnet";
local decon_mod = import "decon.jsonnet";
local cv_mod = import "crossviews.jsonnet";
// local dnnroi_mod = import "dnnroi.jsonnet";

local tlas = import "tlas.jsonnet";



/// Return an object with various TDM nodes and subgraphs for one tpc
///
function(tpc, control, pg) {

    // Pull these out to help debug dumping.
    tpc: tpc,
    anode: tpc.anode,

    local cv = cv_mod(control, pg=pg),
    local fans = fans_mod(control),
    local frame = frame_mod(control),
    local decon = decon_mod(tpc, control),
    // local dnnroi = dnnroi_mod(control),

    /// We use notation:
    ///
    /// Node: [nin]InType -> OutType[nout]
    ///
    /// A "set" means a port passees ITorchTensorSet.
    /// A "tensort" means a port passes ITorchTensor.

    // Node: [1]IFrame -> set[1]: This converts an ADC IFrame into corresponding
    // TDM tensor set.
    frame_to_tdm: frame.to_tdm(tpc),

    // Node: [1]set -> set[2]: frame tensor set fanout eg to send one to output
    // of graph for merging with SPNG results.
    frame_set_fanout: fans.fanout(tpc.name+"frame"),

    // Node: [1]set -> tensor[ngroup]: frame tensor set -> per-group tensors.
    // This breaks up the ADC frame into individual ADC tensors grouped by
    // channel groups.
    frame_set_unpack: frame.tensorset_unpacker(tpc, extra_name="spng"),

    // Node: [1]set -> tensor[ngroup]: directly connect without a fanout.
    frame_direct_unpack: pg.pipeline([$.frame_to_tdm, $.frame_set_unpack]),

    // Node: [ngroup]tensor -> tensor[nview]: apply per-group response decon and
    // merge any split groups to produce a subgraph with 3 per-view tensors of
    // decon signal.  Note, only FR*ER response channel filter is decon'ed.  Not
    // time filter nor "RC" response decon is included.
    response_decon: decon.group_decon_view,

    response_decon_simple: decon.group_decon_simple,

    // Node: [nview]tensor -> tensor [nview] downsample
    downsampler(name, downsample_factor=4, views=[0,1,2]):: pg.crossline([pg.pnode({
        type:'SPNGResampler',
        name: tpc.name + name + "_downsample_v" + std.toString(view),
        data: {
            ratio: 1.0/downsample_factor,
            // fixme: do we need integral norm?
        } + control
    }, nin=1, nout=1) for view in views]),

    // Node: [nview]tensor -> tensor [nview] rebin (interval domain downsample)
    rebinner(name, factor=4, norm="interpolation", views=[0,1,2]):: pg.crossline([pg.pnode({
        type:'SPNGRebinner',
        name: tpc.name + name + "_rebin_v" + std.toString(view),
        data: {
            factor: factor,
            norm: norm,
        } + control
    }, nin=1, nout=1) for view in views]),


    // Node: [nview]tensor -> tensor[nview].  Apply a time filter.
    //
    // They are interpreted as:
    //
    // - wiener :: provides input to thresholding prior to cross views
    //
    // - dnnroi :: provides the non-cross-views image for DNNROI style inference.
    //
    // - gauss :: provides final signals to which ROIs are applied.
    local filter_names = ["gauss", "wiener", "dnnroi"],
    time_filter(filter_name, views=[0,1,2])::
        decon.time_filter_views(filter_name, views=views),
        
    local time_filters = {
        [filter_name]: $.time_filter(filter_name)
        for filter_name in filter_names
    },

    // Node: [nview]tensor -> tensor[nview] apply threshold to produce boolean tensor
    threshold: cv.threshold_views(tpc),

    // Node: [nview]tensor->tensor[nview].  3 inputs are fanned to three
    // crossviews each of three inputs making a fully-connected net.
    crossfan(view_crossed=[1,1,0]): cv.crossfan(tpc, view_crossed=view_crossed),






    // // Node: 2x3 -> 3.  Connect the dnnroi filters and output of crossviews
    // dnnroi_prepare: dnnroi.prepare(tpc.name, time_filters.dnnroi, $.crossfan),
    
    // // Node: 3 -> 3.  Connect dnnroi pre with forward inference.
    // dnnroi_forward: dnnroi.forward(tpc.name, $.dnnroi_prepare),

    // // Node: 2x3 -> 3.  Apply rois.
    // dnnroi_apply: dnnroi.apply(tpc.name, time_filters.gauss, $.dnnroi_forward),

    // Node: nview->1 set.  Collect view tensors into a tensor set.
    frame_set_repack: frame.tensorset_view_repacker(tpc.name),
    frame_set_repack_1: frame.tensorset_view_repacker(tpc.name, multiplicity=1),


    // Node: 2->1.  Collect original tensor set and the repack with signals.
    frame_set_fanin: fans.fanin(tpc.name+"fanin", type='TensorSet'),

    // Node: 1->1
    frame_from_tdm: frame.from_tdm(tpc),


    // The "connected bypass" here is a bit weird given "bypass" itself is a bit
    // weird.  See frame.bypass for a very detailed comments to understand the
    // topology of the subgraph.
    //
    // The bypass's job is to bracket a graph with fanout on input and fanin on
    // output in order to allow data to go along one fan finger and bypass some
    // subgraph connecting the other finger.
    //
    // Here we connect the available ports of the basic bypass with nodes.
    //
    // The returned pnode is conceptually closed and caller should instead
    // connect to these nodes:
    //
    // Connect IFrame source to .frame_to_tdm.
    // Connect IFrame sink to .frame_from_tdm.
    //
    // Connect .frame_set_unpack to a 4-group sink.
    // Connect .frame_set_repack to a 3-view source.
    //
    // Note, the frame tensor set from the IFrame is sent along fan finger 0 so
    // that it is first in the final fanin as required by TDM.
    frame_bypass: frame.connect_bypass(frame.bypass(tpc.name+"bypass", 0, 0),
                                       sources=[$.frame_to_tdm, $.frame_set_repack],
                                       targets=[$.frame_set_unpack, $.frame_from_tdm]),



}
