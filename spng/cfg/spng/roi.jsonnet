// Construct ROI related subgraphs including running dnnroi

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local fans = import "fans.jsonnet";
local util = import "util.jsonnet";

{
    /// A 3->1 pnode configuring DNNROI-like ROI finding to produce signals.
    ///
    /// It produces an object with attributes giving pnodes that need to be
    /// connected to inputs/output by the caller:
    ///
    /// - crossviews :: a 1 port sink accepting a crossviews image for the view
    ///
    /// - gauss :: a 1 port sink accepting "gauss" filtered iamge for the view.
    ///
    /// - dense :: a 1 port sink accepting "dense" for the view data stacking
    ///   and input to dnnroi, usually loose LF filtered decond.
    ///
    /// - signals :: a 1 port source producing final signals for the view
    ///
    /// Caller should connect the "sinks" and "source".
    ///
    /// The name must be unique to the view.
    forward_subgraph(name, forward, mps=["mp2", "mp3"])::
        local this_name = name;
        local nmps = std.length(mps);
        // 1->nmps.
        local ex = pg.pnode({
            type: 'SPNGCrossViewsExtract',
            name: this_name,
            data: {
                extraction: mps
            },
        }, nin=1, nout=nmps);
        /// Stack dense + mps.   This will provide the iport for the "dense" pnode
        local op = util.reduce_one("stack", dim=-3, multiplicity=nmps+1, name=this_name);
        local fwd = pg.pnode({
            type: 'SPNGForward',
            name: this_name,
            data: {
                forward: wc.tn(forward),
            },
        }, nin=1, nout=1, uses=[forward]);
        // gauss+roi in, signals out.  this provides the oport.
        local mul = pg.pnode({
            type: 'SPNGReduce',
            name: this_name+"roi",
            data: {
                operation: 'mul',
            },
        }, nin=2, nout=1);
        local rbl = pg.pnode({
            type: 'SPNGRebaseliner',
            name: this_name,
            data: { }
        }, nin=1, nout=1);
        {
            /// This pnode will hold all the nodes we made here, exposing only
            /// the rbl's oport.
            signals: pg.intern(centernodes=[ex, op, fwd],
                               outnodes=[rbl],
                               edges=[
                                   // The extracted MPs into the stack op
                                   pg.edge(ex,op,mpind, mpind+1)
                                   for mpind in wc.iota(nmps)
                               ] + [
                                   pg.edge(op, fwd),
                                   pg.edge(fwd, mul, 0, 1),
                                   pg.edge(mul, rbl)
                               ]),

            /// A "sink" pnode that exposes the iport of the extract node
            crossviews: pg.intern(iports=[ex.iports[0]]),

            /// A "sink" pnode that exposes the 0th iport of the multiply that applies ROI
            gauss: pg.intern(iports=[mul.iports[0]]),

            /// A "sink" pnode that exposes the  0th iport of the stack op.
            dense: pg.intern(iports=[op.iports[0]]),
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
    connect_forward(name, forward,
                    crossviews, gauss, dense, forward,
                    mps=["mp2", "mp3"])::
        local guts = $.forward_subgraph(name, forward, mps=mps);
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
            },
        }, nin=2, nout=1);
        local rbl = pg.pnode({
            type: 'SPNGRebaseliner',
            name: name,
            data: { }
        }, nin=1, nout=1);
        pg.intern(centernodes=[threshold, gauss, mul],
                  outnodes=[rbl],
                  edges=[
                      pg.edge(gauss,     mul, 0, 0),
                      pg.edge(threshold, mul, 0, 1),
                      pg.edge(mul, rbl)
                  ]),
}
