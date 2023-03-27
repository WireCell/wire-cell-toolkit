// A function on an anode which gives factories for creating a
// configuration subgraph for imaging.

local util = import "util.jsonnet";
local pg = import "pgraph.jsonnet";
local wc = import "wirecell.jsonnet";

function (anode) {

    local ident = util.idents(anode),
    

    // The recommended slicer.  There is almost no reasonable default
    // to each detector variant may as well create the MaskSlices
    // directly.
    slicing :: function(tag, errtag, min_tbin=0, max_tbin=0, ext="",
                        span=4, active_planes=[0,1,2], masked_planes=[], dummy_planes=[]) 
        pg.pnode({
            type: "MaskSlices",
            name: ident+ext,
            data: {
                tag: tag,
                error_tag: errtag,
                tick_span: span,
                anode: wc.tn(anode),
                min_tbin: min_tbin,
                max_tbin: max_tbin,
                active_planes: active_planes,
                masked_planes: masked_planes,
                dummy_planes: dummy_planes,
            },
        }, nin=1, nout=1, uses=[anode]),


    // This slicing does not handle charge uncertainty nor bad channels.
    simple_slicing :: function(span=4, tag="") pg.pnode({
        type: "SumSlices",
        name: ident,
        data: {
            tag: tag,
            tick_span: span,
            anode: wc.tn(anode),
        },
    }, nin=1, nout=1, uses=[anode]),

    // A function that returns a tiling node for single faced anodes.
    // When used in the context of many tilings per anode face, a
    // unique name extension "ext" must be given.
    single_tiling :: function(face, ext="") pg.pnode({
        type: "GridTiling",
        name: "%s-%d%s"%[ident, face, ext],
        data: {
            anode: wc.tn(anode),
            face: face,
        }
    }, nin=1, nout=1, uses=[anode]),

    // A subgraph to do tiling for wrapped (two faced) anodes.  
    // When used in the context of many tilings per anode face, a
    // unique name extension "ext" must be given.
    wrapped_tiling :: function(ext="") {

        local slice_fanout = pg.pnode({
            type: "SliceFanout",
            name: ident+ext,
            data: { multiplicity: 2 },
        }, nin=1, nout=2),

        local tilings = [$.single_tiling(face, ext) for face in [0,1]],

        local blobsync = pg.pnode({
            type: "BlobSetSync",
            name: ident+ext,
            data: { multiplicity: 2 }
        }, nin=2, nout=1),

        ret: pg.intern(
            innodes=[slice_fanout],
            outnodes=[blobsync],
            centernodes=tilings,
            edges=
                [pg.edge(slice_fanout, tilings[n], n, 0) for n in [0,1]] +
                [pg.edge(tilings[n], blobsync, 0, n) for n in [0,1]],
            name=ident+ext),
    }.ret,

    // A tiling subgraph matching the anode faces.  The "uname" must
    // be unique across all application of this function to a subgraph
    // servicing a given anode.
    tiling :: function(uname="")
        if std.length(std.prune(anode.data.faces)) == 2 then
            $.wrapped_tiling(uname)
        else if std.type(anode.data.faces[0]) == "null" then
            $.single_tiling(1, uname)
        else
            $.single_tiling(0, uname),
    
    // Blobify: produce a subgraph which consumes a frame a, performs
    // slicing on the slices performs tiling finally producing a
    // blobset.
    // 
    // The subgraph may itself be composed of subgraphs each formed by
    // replicating a nominally linear pipeline segment into multiple
    // parallel pipelines with identical code but unique configuration
    // parameters.  There are three possible points of such pipeline
    // splitting.  Starting from downstream / inner-most points in the
    // graph and working upstream / outward:
    //
    // - face :: A tiling operates on a single anode face.  A single
    // time slice from a wrapped anode must be processed by two
    // separate tiling nodes each focused on one anode face.
    //
    // - slicing :: The slicer preceding one or a pair of tiling must
    // be configured with arrays of plane identifiers that determine
    // how channel activity in a slice is determined.  Typically only
    // one category of activity should be set for any given plane
    // though the code allows for a plane to reside in multiple
    // categories.  The categories, in their order of application,
    // are:
    // 
    //   - active :: assign the sum of all values from the channel
    //   across the current slice as the slice activity for that
    //   channel.
    // 
    //   - dummy :: assign a fixed, configured value and uncertainy as
    //   the activity for ALL channels in the plane.  Default is small
    //   value with large uncertainty.
    // 
    //   - masked :: assign a configured value and uncertainty as the
    //   slice activity for any channel with an entry in the frame's
    //   "bad" channel mask map that overlaps the slice.  Default
    //   small value with large uncertainty.
    // 
    // - strategy :: A logical grouping of a number of slicing
    // subgraphs.  For example, all slicing which span some
    // permutation of plane categorizations may comprise a strategy.
    // 
    // Some predefined strategies named:
    blobify_strategies : {
        // All 3view and 2view with active channels plus third view
        // with active or bad channels.
        live: [
            {
                active: [0,1,2],
            },
            {
                active: [0,1],
                masked:   [2],
            },
            {
                active: [1,2],
                masked:   [0],
            },
            {
                active: [0,2],
                masked:   [1],
            },
        ],
        // The pseudo-inverse of "live".  Only consider bad channels
        // in two views.  Third view poses no constraints on blob
        // forming by making all channels active with fixed values.
        dead: [
            {
                dummy:    [2],
                masked: [0,1],
            },
            {
                dummy:    [0],
                masked: [1,2],
            },
            {
                dummy:    [1],
                masked: [0,2],
            },
        ]
    },

    // This grouping may also be used to supply a second fanout layer
    // simply to work around fan size limits.
    //
    // Where any of these splittings result in more than one parallel
    // pipeline, the multiple pipelines are bounded by fanout/fanin
    // pairs.  Where there is no splitting (eg, tiling with single
    // faced anodes) the fans are omitted as a minor optimization.

    // Make blobify subgraph consuming frame and producing blobset.
    // It will have a singel slicer and one or two tiling depending on
    // number of anode faces.

    // - uname :: Unique name in context of applying function on a given anode
    // - tag, errtag :: The trace tags for input signal and signal error traces
    // - planes :: object giving the plane categorization arrays
    // - trange :: time slice tbin min/max range in ticks
    // - span :: the time slice span in ticks

    blobify_single :: function(uname, tag, errtag,
                               planes={}, trange=[0,0], span=4)
        local p = { active: [], masked: [], dummy: [] } + planes;
        pg.pipeline([$.slicing(tag, errtag, trange[0], trange[1], uname,
                               span, p.active, p.masked, p.dummy),
                     $.tiling(uname)], "blobify_single_"+uname),

    // Make a full blobify subgraph over a slicing strategy array of
    // slicing objects.
    blobify :: function(uname, tag, errtag,
                        strategy=[], trange=[0,0], span=4)
        local nstrats = std.length(strategy);
        local pipes = [
            $.blobify_single('%s-s%d'%[uname,sind], tag, errtag, strategy[sind], trange, span)
            for sind in std.range(0,nstrats-1)
        ];
        if nstrats == 1
        then pipes[0]
        else
        local fanout = pg.pnode({
            type:'FrameFanout',
            name: uname,
            data: { multiplicity: nstrats },
        }, nin=1, nout=nstrats);
        local fanin = pg.pnode({
            type:'BlobSetMerge',
            name: uname,
            data: { multiplicity: nstrats },
        }, nin=nstrats, nout=1);
        pg.fan.pipe('FrameFanout', pipes, 'BlobSetMerge', uname),

    // A node that clusters blobs between slices w/in the given number
    // of slice spans.  Makes blob-blob edges
    clustering :: function(spans=1.0) pg.pnode({
        type: "BlobClustering",
        name: ident,
        data:  { spans : spans }
    }, nin=1, nout=1),    
    
    // A node that groups wires+channels into "measurements".  Makes
    // blob-measure edges.
    grouping :: function() pg.pnode({
        type: "BlobGrouping",
        name: ident,
        data:  { }
    }, nin=1, nout=1),

    // The recommended charge solving
    charge_solving :: function(meas_val_thresh=0.0,
                               meas_err_thresh=1.0e9,
                               blob_val_thresh=0,
                               whiten=true)
        pg.pnode({
            type: "ChargeSolving",
            name: ident,
            data:  {
                "meas_value_threshold": meas_val_thresh,
                "meas_error_threshold": meas_err_thresh,
                "blob_value_threshold": blob_val_thresh,
                "whiten": whiten,
            }
        }, nin=1, nout=1),

    // This solver is simplistic, prefer charge solving.
    blob_solving :: function(threshold=0.0) pg.pnode({
        type: "BlobSolving",
        name: ident,
        data:  { threshold: threshold }
    }, nin=1, nout=1),    

    // A function that projects blobs back to frames.  
    reframing :: function(tag="") pg.pnode({
        type: "BlobReframer",
        name: ident,
        data: {
            frame_tag: tag,
        }
    }, nin=1, nout=1),
}

