// Flexible ways to build subgraphs with TDM nodes.
//
// Example cmd line test
//
// @code{sh}
// wcpy pgraph dotify -A stage=graph_crossview -A finalize=true --no-services --no-params test-tdm.pdf
// @endcode


local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

local fans = import "fans.jsonnet";
local frame = import "frame.jsonnet";
local decon = import "decon.jsonnet";
local cv = import "crossviews.jsonnet";
local dnnroi = import "dnnroi.jsonnet";

local tlas = import "tlas.jsonnet";

/// Return an object with various TDM subgraphs for one tpc
function(tpc, control) {

    // allow dumping of some details
    tpc: tpc,
    anode: tpc.anode,

    // Node: 1->1: IFrame -> frame tensor set.  This converts an ADC IFrame into
    // corresponding TDM tensor set.
    frame_to_tdm: frame.to_tdm(tpc),

    // Node: 1->2: frame tensor set fanout eg to send one to output of graph for
    // merging with SPNG results.
    frame_set_fanout: fans.fanout(tpc.name+"frame"),

    // Node: 1->ngroup: frame tensor set -> per-group tensors.  This breaks up the
    // ADC frame into individual ADC tensors grouped by channel groups. 
    frame_set_unpack: frame.tensorset_unpacker(tpc),

    // Node: 1->ngroup: directly connect without a fanout. 
    frame_direct_unpack: pg.pipeline([$.frame_to_tdm, $.frame_set_unpack]),

    // Node: ngroup->nview: apply per-group response decon and merge any split groups
    // to produce a subgraph with 3 per-view tensors of decon signal.  Note,
    // only FR*ER response channel filter is decon'ed.  Not time filter nor "RC"
    // response decon is included.
    response_decon: decon.group_decon_view(tpc),



    unpack_to_decon: pg.shuntline($.frame_set_unpack, $.response_decon),

    // FIXME: this is temporary.  We definitely do NOT want to leave this as-is.
    // It is mathematically correct to apply the downsample after the FR*ER
    // decon.  This means not having to do one downsample per filter.  In fact,
    // it is even faster if we apply the downsample while FR*ER decon is still
    // in Fourier space.  However, any downsampling done upstream of filtering
    // requires the filter kernel to have its sample period adjusted.  For right
    // now we do it after filtering.
    downsampler(name):: pg.crossline([pg.pnode({
        type:'Resampler',
        name: tpc.name + name + "_downsample_v" + std.toString(view),
        data: {
            ratio: 0.25,
            // fixme: do we need integral norm?
        },
    }, nin=1, nout=1) for view in [0,1,2]]),



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
        [filter]: pg.shuntline(decon.time_filter_views(tpc, filter),
                               $.downsampler(filter))
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
    
    crossviews_repack: pg.shuntline($.threshold_crossviews, $.frame_set_repack),

    // Node: 1->3. A subgraph consuming IFrame and producing 3 crossview
    // tensors.  It excludes gauss and dnnroi filtered processing.
    graph_crossview: pg.intern(
        innodes=[$.frame_direct_unpack],
        centernodes=[
            $.unpack_to_decon,
            $.crossview_filters,
            $.filter_threshold,
            $.threshold_crossviews,
        ],
        outnodes=[
            $.crossviews_repack
        ],
        // note, edges are handled already
    ),

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


    // all: pg.components([$.frame_filters, $.threshold_cv, $.dnnroi_prepare]),



}
