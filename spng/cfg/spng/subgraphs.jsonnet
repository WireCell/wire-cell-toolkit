// This produces a large factory object giving a menus of many subgraph
// construction methods with minimal Jsonnet dependencies.

local wc = import "wirecell.jsonnet";
local real_pg = import "pgraph.jsonnet";
local frame_js = import "spng/frame.jsonnet";
local fans_js = import "spng/fans.jsonnet";


// Make a subgraphs factory object.
//
// @param tpc An detconf object describing the TPC context
// @param control A control object with verbosity, device, etc.
// @param pg A pgraph-like module.
// @param context_name A name to combine with the tpc name for all pnodes made.
//
// These general categories of configuration objects are returned:
//
// - a "service" component
// - an atomic pnode
// - a per-view pipeline
// - an all-view "crossline" with 3 inputs and/or 3 outputs
// - a structure of named "pseudo sinks and sources".
//
// A pseudo sink/source is a pnode that the caller should connect.  It is
// "pseudo" as it only appears to be a sink or a source to the caller while in
// fact the sinks and/or the sources hold edges connecting the sinks and
// sources.  This is done as a way to collect ports into separate, meaningful
// sets and facilitate the user making further connections.
function(tpc, control={}, pg=real_pg, context_name="") {

    local fans = fans_js(control),


    // A function to make a unique name in the context of a method.
    this_name(extra_name, meth_name="") :: tpc.name+context_name+meth_name+extra_name,

    //
    //  Decon related
    //


    /// A kernel providing a Fourier-space filter component (not a pnode).
    /// The axis config objects are, eg, made by tpc.filter_axis().
    filter_kernel(axis_config,  // list of per-axis config objects
                  extra_name="")::
        {
            type: "SPNGFilterKernel",
            name: $.this_name(extra_name),
            data: {
                axis: axis_config,
            } + control,
        },
    
    
    /// A kernel providing a response built from FR*ER (not a ponode)
    response_kernel(view_index, // the index of the view in the FR to use
                    extra_name="")::
        {
            type: "SPNGResponseKernel",
            name: $.this_name(extra_name, "v" + std.toString(view_index)),
            data: {
                field_response: wc.tn(tpc.fr),
                elec_response: wc.tn(tpc.er),
                // Take care that this should merely cover where ER is non-zero
                elec_duration: tpc.er.data.nticks*tpc.er.data.tick,
                period: tpc.adc.tick, // sample period of kernel time dimension.
                plane_index: view_index,
                // Scale to units of ADC
                scale: -1 / tpc.adc.lsb_voltage,
            } + control,
            uses: [tpc.fr, tpc.er],
        },

    /// A deconvolution kernel.
    decon_kernel(filt,          // a kernel providing the filter (numerator)
                 resp,          // a kernel providing the response (denominator)
                 extra_name="")::
        {
            type: "SPNGDeconKernel",
            name: $.this_name(extra_name),
            data: {
                filter: wc.tn(filt),
                response: wc.tn(resp),
            } + control,
            uses: [filt, resp]
        },

    /// A decon kernel with filters and responses determined from the tpc group.
    tpc_group_decon_kernel(group_index, extra_name="")::
        local group = tpc.view_groups[group_index];
        local name = $.this_name(extra_name, group.name);
        local plane_index = group.view_index;
        local fk = $.filter_kernel([tpc.filters[plane_index].channel.decon,{kind:"none"}],
                                   extra_name=extra_name);
        local rk = $.response_kernel(plane_index, extra_name=extra_name);
        $.decon_kernel(fk, rk, extra_name=extra_name),

    /// Return pnode convolve intput data with a given kernel
    kernel_convolve(kernel,     // the object providing the convolution kernel 
                    axes_config, // list of per-axis config objects
                    // See KernelConvolve for following options
                    tag="",
                    datapath_format="",
                    extra_name="")::
        pg.pnode({
            type: "SPNGKernelConvolve",
            name: $.this_name(extra_name),
            data: {
                kernel: wc.tn(kernel),
                axis: axes_config,
                tag: tag,
                datapath_format: datapath_format,
                faster: true,       // use "faster DFT size"
            } + control
        }, nin=1, nout=1, uses=[kernel]),

    /// Return kernel convolve pnode for one FR*ER decon node for one group in a
    /// TPC.  The node will be named uniquely for that context.  If more than
    /// one decon in the same tpc+group, pass extra_name.
    tpc_group_decon_frer(group_index, kernel, tag="", extra_name="")::
        local group = tpc.view_groups[group_index];
        local name = tpc.name + group.name + extra_name;
        //local co = $.convo_options(group.connection, group.view_index);
        local co_channel = [
            {cyclic: false, crop: -2},  // one plane
            {cyclic: false, crop: -2},  // two concatenated planes
            {cyclic: true,  crop:  0}, // two wrapped planes
        ][group.connection];
        local co_time = {
            cyclic: false, // time should always be linear
            crop: 0,       // don't crop yet, leave padding for filters
            baseline: true,
            roll: 0,
        };
        $.kernel_convolve(kernel, [co_channel, co_time],
                          datapath_format='/traces/group/' + std.toString(group_index) + "/KernelConvolve/" + name,
                          tag="decon", extra_name=group.name+extra_name),

    /// Return nodes[0] if length one else return a subgraph with a fanin that stacks the tensors.
    tpc_group_fanin(nodes, extra_name="")::
        local n = std.length(nodes);
        local fanin = pg.pnode({
            type: 'SPNGReduce',
            name: $.this_name(extra_name, std.join('-', [node.name for node in nodes])),
            data: {
                multiplicity: n,
                operation: "cat",
                dim: -2         // concatenate along channel dimension
            } + control
        }, nin=n, nout=1);
        
        if n == 1
        then nodes[0]
        else pg.intern(innodes=nodes, outnodes=[fanin],
                       edges=[pg.edge(it.value, fanin, 0, it.index) for it in wc.enumerate(nodes)]),
            

    /// Given array of nodes that span view groups, return array that span views
    /// fixme: move this out to a more generic jsonnet.
    local collect_groups_by_view(nodes)=
        [$.tpc_group_fanin([
            nodes[gi.index]
            for gi in wc.enumerate(tpc.view_groups)
            if gi.value.view_index == view_index])
         for view_index in [0,1,2]],

    /// Produce an ngroup->nview subgraph doing FR*ER decon.
    tpc_group_decon_views(extra_name=""):
        local groups = tpc.view_groups;
        local ngroups = std.length(groups);

        local decon_kernels = [$.tpc_group_decon_kernel(group_index, extra_name=extra_name)
                               for group_index in wc.iota(ngroups)];
        local decon_frers = [$.tpc_group_decon_frer(it.index,
                                                    decon_kernels[it.value.view_index],
                                                    extra_name=extra_name)
                             for it in wc.enumerate(groups)
                            ];
        local view_frers = collect_groups_by_view(decon_frers);
        pg.intern(innodes=decon_frers, outnodes=view_frers),




    /// Pnode to convert frame to TDM tensor set.
    frame_to_tdm(extra_name=""):: frame_js(control).to_tdm(tpc, extra_name=extra_name),


    /// An input tensor set to N-tensor outputs specifically for tpc groups.
    /// Map a tensors in a set by an iteration of a datapath pattern over a set
    /// of a count to a fan out of individual tensors.
    tpc_group_unpacker(
                       // How to locate each tensor, %d is interpolated as an index
                       datapath_pattern="/frames/\\d+/tags/null/rules/0/groups/%d/traces",
                       extra_name="")::
        local ngroups = std.length(tpc.view_groups);
        pg.pnode({
            type: 'SPNGTorchSetUnpacker',
            name: $.this_name(extra_name, "groups"),
            data: {
                selections: [{datapath:datapath_pattern % group} for group in wc.iota(ngroups) ],
            } + control
        }, nin=1, nout=ngroups),


    /// A 1->3 subgraph with frame input and 3 per-view decon output.
    frame_decon(extra_name=""):: pg.shuntlines([
        $.frame_to_tdm(extra_name=extra_name),
        $.tpc_group_unpacker(extra_name=extra_name),
        $.tpc_group_decon_views(extra_name=extra_name)]),


    /// An object providing pseudo sources for faning out from a 3-view pnode
    /// (eg decon) to multi-view nodes to use as sources for feeding initial_roi
    /// ("wiener" filter) and dnnroi_dense_Views ("dnnroi" filter).  This omits
    /// "gauss", see fanout_for_inference().  The dnnroi_views list which views
    /// provide dnnroi fodder.
    fanout_for_training(crossed_views=[1,1,0], extra_name=extra_name)::
        fans.fanout_select($.this_name(extra_name, "dnnroi_training"),
                           std.length(crossed_views), targets_list=[
                               if is_crossed == 1
                               then ["wiener", "dense"]
                               else ["wiener"]
                               for is_crossed in crossed_views]),


    // Like fanout_for_training but include gauss to connect with apply_roi_views.
    fanout_for_inference(crossed_views=[1,1,0], extra_name=extra_name)::
        fans.fanout_select($.this_name(extra_name, "dnnroi_training"),
                           std.length(crossed_views), targets_list=[
                               if is_crossed == 1
                               then ["wiener", "gauss", "dense"]
                               else ["wiener", "gauss"]
                               for is_crossed in crossed_views]),

    /// One filter along time dimension for one view
    time_filter_view(filter, view_index, extra_name="")::
        local vis = std.toString(view_index);
        local name = $.this_name(extra_name, 'v'+vis);
        local fk = $.filter_kernel([{kind:"none"}, tpc.filters[view_index].time[filter]], extra_name=extra_name);
        // local co = $.convo_options(0, view_index);
        local co_channel = {padding: "none", dft: false};
        local co_time = {padding: "none", dft: true} + tpc.filters[view_index].time.options;
        $.kernel_convolve(fk, [co_channel, co_time],
                          datapath_format='/traces/view/' + vis + "/KernelConvolve/" + name,
                          tag=filter, extra_name=extra_name),


    /// Apply a time filter across each of the views for a 3->3 subgraph.
    time_filter_views(filter="gauss", views=[0,1,2], extra_name="")::
        local nodes = [
            $.time_filter_view(filter, view_index, extra_name="v"+std.toString(view_index)+"_"+filter+extra_name)
            for view_index in views];
        pg.crossline(nodes),
        
    /// A 3->3 subgraph that downsamples its input along the time dimension
    downsample_views(downsample_factor=4, views=[0,1,2], extra_name=""):: pg.crossline([pg.pnode({
        local meth_name = "_downsample"+std.toString(downsample_factor)+"_v" + std.toString(view),
        type:'SPNGResampler',
        name: $.this_name(extra_name, meth_name),
        data: {
            ratio: 1.0/downsample_factor,
        } + control
    }, nin=1, nout=1) for view in views]),


    /// A 3->3 subgraph that applies a tpc cross view thresholds
    cross_threshold_views(extra_name=""):: pg.crossline([pg.pnode({
        local meth_name = "_cross_v" + std.toString(it.index),
        type: 'SPNGThreshold',
        name: $.this_name(extra_name, meth_name),
        data: it.value + control,
    }, nin=1, nout=1) for it in wc.enumerate(tpc.crossview_thresholds)]),


    /// A 3->3 subgraph that performs tight ROI finding on decon inputs and
    /// outputs giving the tight ROIs.  Fixme: support doing rebin on decon
    /// once which needs adjusting the filters to the different sampling.
    tight_roi(rebin=4, views=[0,1,2], extra_name="")::
        pg.shuntlines([
            $.time_filter_views("wiener", views=views, extra_name=extra_name),
            $.downsample_views(rebin, views=views, extra_name=extra_name),
            $.cross_threshold_views(extra_name=extra_name)]),


    scale_views(scale, views=[0,1,2], extra_name=""):: pg.crossline([pg.pnode({
        local meth_name = 'v'+std.toString(view)+"_scale",
        type:'SPNGTransform',
        name: $.this_name(extra_name, meth_name),
        data: {
            operations: [
                { operation: "scale", scalar: scale },
            ],                
        } + control,
    }, nin=1, nout=1) for view in views]),

    /// An N->N subgraph applying the "dnnroi" filter and rebinning.  N is the
    /// length of "views" and the "views" array gives the view indices for the
    /// ports.  When the "dense" data (usually looseLF) is upstream of dnnroi,
    /// we want views 0, 1.
    dnnroi_dense_views(scale=1.0/4000, views=[0,1], rebin=4, extra_name="")::
        pg.shuntlines([
            $.time_filter_views("dnnroi", views=views, extra_name=extra_name),
            $.downsample_views(rebin, views=views, extra_name=extra_name),
            $.scale_views(scale, views, extra_name="_dnnroi_dense"+extra_name),
        ]),


    /// A subgraph with one input and one output port passing ITorchTensorSet
    /// though CellViews.  
    cellviews_tensorset(out_views=[0,1],   // order and views to output
                        chunk_size=0,      // large detectors may want to set this
                        uvw_index=[0,1,2], // order of per-view Boolean tensors in set
                        extra_name="")::
        pg.pnode({
            type: 'SPNGCellViews',
            name: tpc.name + extra_name,
            data: {
                anode: wc.tn(tpc.anode),
                face_idents: tpc.faces,
                uvw_index: uvw_index,
                out_views: out_views,
                chunk_size: chunk_size,
            },
        }, nin=1, nout=1, uses=[tpc.anode]),
    
        

    /// A subgraph with three input ports and N output ports.  Input ports
    /// accept, in order, U, V and W ITorchTensor.  Output ports provide
    /// ITorchTensor for each target view in given order.  Large detectors may
    /// want to limit memory usage by specifying a chunk_size.
    cellviews_tensors(out_views=[0,1], chunk_size=0, extra_name="")::
        local this_name = $.this_name(extra_name, '_cellviews');
        local packer = pg.pnode({
            type: 'SPNGTensorPacker',
            name: this_name,
            data: {
                multiplicity: 3
            } + control
        }, nin=3, nout=1);
        local cellviews = $.cellviews_tensorset(tpc, out_views, chunk_size, extra_name=extra_name);
        local nout = std.length(out_views);
        local unpacker = pg.pnode({
            type: 'SPNGTensorUnpacker',
            name: this_name,
            data: {
                selections: [{index: ind} for ind in wc.iota(nout)],
            } + control
        }, nin=1, nout=nout);
        pg.pipeline([packer, cellviews, unpacker]),

    
    /// Make a transform crossline
    transform_views(operation, scalar=0, dims=[], views=[0,1,2], extra_name=""):: pg.crossline([pg.pnode({
        type: 'SPNGTransform',
        name: $.this_name(extra_name, 'v'+std.toString(view) + '_'+operation),
        data: {
            operation: operation,
            scalar: scalar,
            dims: dims,
        },
    }, nin=1, nout=1) for view in views]),
        
    /// Return list of pnode that apply a 2->1 reduction per view.
    reduce_views_list(operation, dim=0, views=[0,1,2], extra_name=""):: [pg.pnode({
        type: 'SPNGTransform',
        name: $.this_name(extra_name, 'v'+std.toString(view) + '_'+operation),
        data: {
            operation: operation,
            dim: dim,
        },
    }, nin=2, nout=1) for view in views],

    /// Return object with pseudo sinks/sources:
    /// - mp_sink :: an N-port sink of cellviews tensors (with feature dimension)
    /// - dense_sink :: an N-port sink of "dense' tensors (with no feature dimension)
    /// - source :: an N-port source of dense added to the feature dimension
    dnnroi_stack_features(views=[0,1], extra_name=""):: {
        local features = $.transform_views("unsqueeze", dims=[1], views=views, extra_name=extra_name),
        local stacks = $.reduce_views_list("cat", dim=-3, views=views, extra_name=extra_name),
        dense_sink: features,
        source: pg.intern(outnodes=stacks, edges=[
            pg.edge(features, stacks[vi.index], vi.index, 1)
            for vi in wc.enumerate(views)]),
        mp_sink: pg.intern(iports=[s.iports[0] for s in stacks]),
    },        

    // Connect source mp and dense to the dnnroi stack operation.  Return pnode
    // that acts as a source of dnnroi input fodder.  The inputs are carried by
    // this object but the user must still handle connecting any iports of the
    // inputs.
    connect_dnnroi_stack(mp, dense, views=[0,1], extra_name="")::
        local sf = $.dnnroi_stack_features(views, extra_name);
        pg.intern(centernodes=[pg.shuntline(mp, sf.mp_sink), pg.shuntline(dense, sf.dense_sink)],
                  outnodes=[sf.source]),


    /// An Nview->Nview applying the "gauss" filter and rebinning.
    gauss_dense_views(views=[0,1,2], rebin=4, extra_name="")::
        pg.shuntlines([
            $.time_filter_views("gauss", views=views, extra_name=extra_name),
            $.downsample_views(rebin, views=views, extra_name=extra_name)]),

    /// Return object with pseudo sinks/sources:
    /// - roi_sink :: An Nview port sink to connect to a source of ROIs
    /// - dense_sink :: An Nview port sink to connect to a source of dense (eg "gauss")
    /// - source :: An Nview port source of signals after ROIs are applied to dense.
    //apply_roi_views(views=[0,1,2], extra_name="")::

}

