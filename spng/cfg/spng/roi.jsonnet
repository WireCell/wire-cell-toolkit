// Construct ROI related subgraphs including running dnnroi

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local util = import "util.jsonnet";

function(control={})
{
    /// 3->1 overall, but returns an object with various sinks and one source
    /// that bracket a DNNROI-like ROI finding sugraph.
    ///
    /// It produces an object with attributes giving pnodes that need to be
    /// connected to inputs/output by the caller:
    ///
    /// - crossviews :: a 1 port sink accepting a crossviews image for the view
    ///
    /// - gauss :: a 1 port sink accepting "gauss" filtered full image for the
    /// view.  This is assumed to be at original sampling and not "rebinned".
    ///
    /// - dense :: a 1 port sink accepting "dense" full image which is typically
    /// loose LF filtered decon.  It must match what the forward model was
    /// trained with.  This is assumed to be at original sampling and not
    /// "rebinned".
    ///
    /// - signals :: a 1 port source producing final signals for the view.  This
    /// is a full image at original sampling and not "rebinned".
    ///
    /// Caller should connect the "sinks" and "source".
    ///
    /// Arguments
    /// @param name A name unique to the view.
    /// @param forward The configuration object for the DNNROI ITensorForward.
    /// @param mps A list of MPn types.
    /// @param rebin The downsample factor to apply prior to forward inference.
    /// @param scale The multiplicative scale applied to all inputs to DNNROI forward.
    ///
    /// The subgraph is very DNNROI specific and bakes in pre/post processing
    /// around the "forward" node.  It is assumed the source to be connected to
    dnnroi_subgraph(name, forward, mps=["mp2", "mp3"], rebin=4, scale=1.0/4000.0)::
        local this_name = name;
        local nmps = std.length(mps);

        // 1->nmps.
        local ex = pg.pnode({
            type: 'SPNGCrossViewsExtract',
            name: this_name,
            data: {
                extraction: mps
            } + control
        }, nin=1, nout=nmps);

        /// 1->1 rebin the "dense" image
        local rebinner = pg.pnode({
            type: "SPNGRebinner",
            name: this_name+"_rebin",
            data: {
                norm: "interpolation",
                factor: rebin,
            } + control
        }, nin=1, nout=1);

        /// nmps+1 -> 1
        /// Stack dense + mps.  An iport will be exposed to connect to "dense".
        /// FIXME: batching may be a problem.  dim:-3 stacks along the dimension
        /// prior to (chan,tick) which would make (nbatch, 3, nchan, ntick).
        /// But, what does DNNROI ingest expect?
        local stack = pg.pnode({
            type:'SPNGReduce',
            name: this_name+"_stack",
            data: {
                operation: "stack",
                dim: -3,
                multiplicity: nmps+1,
            } + control
        }, nin=nmps+1, nout=1);

        /// 1 -> 1
        /// DNNROI scales all three images by 4000.0.
        local pre = pg.pnode({
            type: 'SPNGTransform',
            name: this_name+'_pre',
            data: {
                operations: [
                    /// DNNROI traditionally scales all three images down by
                    /// this seemingly arbitrary value.  The "normalize"
                    /// operation might be more portable?
                    { operation: "scale", scalar: scale },
                    /// DNNROI traditionally works in (ntick,nchan) for last two
                    /// dims while WCT works in (nchan,ntick) so we must bookend
                    /// the forward with a transpose.
                    { operation: "transpose", dims: [-2,-1] }
                ]
            } + control
        }, nin=1, nout=1);

        /// Do inference (eg, DNNROI) or other type of "forward"
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
                nominal: 0.5    // FIXME: a study is needed to best set this
            } + control
        }, nin=1, nout=1);
        // We use rebinner following threshold.  However, we may want to reverse
        // the order in which case we would use resample instead.  
        local unbinner = pg.pnode({
            type: "SPNGRebinner",
            name: this_name+"_unbin",
            data: {
                norm: "maximum",
                factor: -rebin,     // FIXME: must match upstream resampler.
            } + control
        }, nin=1, nout=1);
        // ApplyROI: gauss+roi in, signals out.  this provides the oport.
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
            data: control
        }, nin=1, nout=1);
        local pipe = pg.pipeline([rebinner, stack, pre, fwd, post, thresh, unbinner, mul, rbl]);
        {
            /// This pnode will hold all the nodes we made here, exposing only
            /// the final output.
            signals: pg.intern(outnodes=[pipe]),

            /// A "sink" pnode that exposes the iport of the extract node.  This
            /// connects to stack, avoiding iport[0] which will have the input
            /// rebin connected below.
            crossviews: pg.intern(centernodes=[ex, stack],
                                  iports=[ex.iports[0]],
                                  edges=[
                                      pg.edge(ex, stack, mpind, mpind+1)
                                      for mpind in wc.iota(nmps)
                                  ]),

            /// A "sink" pnode that exposes the iport of the apply ROI multiply that was left open.
            gauss: pg.intern(iports=[mul.iports[1]]),

            /// A "sink" pnode that exposes the  0th iport of the stack op.
            dense: pg.intern(innodes=[rebinner],
                             centernodes=[stack],
                             edges=[pg.edge(rebinner, stack)]),
        },
    
    /// A 1-port source of signal tensors.
    ///
    /// This constructs inputs to forward_subgraph().
    ///
    /// The first three objects are assumed to be 1-node sources.
    ///
    /// For views that get ROIs with simple thresholds see connect_threshold().
    ///
    /// The name must be unique to the view.
    ///
    /// Note, the input "sources" are not exposed in the returned pnode as
    /// innodes.  Typically, they do actually have their own inputs and it is
    /// assumed the caller connects them in some way.
    connect_dnnroi(name, forward,
                    crossviews, gauss, dense, forward,
                    mps=["mp2", "mp3"])::
        local guts = $.dnnroi_subgraph(name, forward, mps=mps);
        pg.intern(centernodes=[crossviews, dense, gauss,
                               guts.crossviews, guts.dense, guts.gauss],
                  outnodes=[guts.signals],
                  edges=[
                      pg.edge(crossviews, guts.crossviews),
                      pg.edge(gauss,      guts.gauss),
                      pg.edge(dense,      guts.dense),
                  ]),
        

    /// A 1-port source of signal tensors.
    ///
    /// This connects inputs to a Rebaseline.
    ///
    /// For connection using a DNNROI-like "forward", see connect_forward().
    ///
    /// Name must be unique to the view.
    connect_threshold(name, threshold, gauss)::
        // gauss+roi in, signals out.  this provides the oport.
        local mul = pg.pnode({
            type: 'SPNGReduce',
            name: name,
            data: {
                operation: 'mul',
            } + control,
        }, nin=2, nout=1);
        local rbl = pg.pnode({
            type: 'SPNGRebaseliner',
            name: name,
            data: { } + control
        }, nin=1, nout=1);
        pg.intern(centernodes=[threshold, gauss, mul],
                  outnodes=[rbl],
                  edges=[
                      pg.edge(gauss,     mul, 0, 0),
                      pg.edge(threshold, mul, 0, 1),
                      pg.edge(mul, rbl)
                  ]),
}
