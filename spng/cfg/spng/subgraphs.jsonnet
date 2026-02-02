// This produces a large factory object giving a menus of many subgraph
// construction methods with minimal Jsonnet dependencies.

local wc = import "wirecell.jsonnet";
local real_pg = import "pgraph.jsonnet";
local frame_js = import "spng/frame.jsonnet";
local fans_js = import "spng/fans.jsonnet";


/**
 * Make a subgraphs factory object.
 * 
 * These general categories of configuration objects are returned:
 *
 * - a "service" component
 * - an atomic pnode
 * - a per-view pipeline
 * - an all-view "crossline" with 3 inputs and/or 3 outputs
 * - a structure of named "pseudo sinks and sources".
 *
 * A pseudo sink/source is a pnode that the caller should connect.  It is
 * "pseudo" as it only appears to be a sink or a source to the caller while in
 * fact the sinks and/or the sources hold edges connecting the sinks and
 * sources.  This is done as a way to collect ports into separate, meaningful
 * sets and facilitate the user making further connections.
 * 
 * @name subgraphs
 * @param tpc An detconf object describing the TPC context
 * @param control A control object with verbosity, device, etc.
 * @param pg A pgraph-like module.
 * @param context_name A name to combine with the tpc name for all pnodes made.
 */
function(tpc, control={}, pg=real_pg, context_name="") {
  
    /// The subgraphs object is mostly independent from other Jsonnet's.  The
    /// exceptions are here and they are made available to the caller.
    local fans = fans_js(control),
    local frame = frame_js(control),
    wc: wc,
    pg: pg,
    fans: fans,
    frame: frame,

    
    /**
     * A function to make a unique name in the context of a method.
     *
     * Produce instance names in a regular way.
     *
     * @name this_name
     * @param extra_name The extra name passed to most functions to qualify a caller's context.
     * @param meth_name Qualify the context of a method in this object.
     * @return A name based on arguments and tpc.name.
     */
    this_name(extra_name, meth_name="") :: tpc.name+context_name+meth_name+extra_name,



    /// Deposet to IFrame of "true signal" aka "depo flux splat".
    splat(extra_name="")::
        local splat=pg.pnode({
            type: "DepoFluxSplat",
            name: $.this_name(extra_name),
            data: {
                anode: wc.tn(tpc.anode),
                field_response: wc.tn(tpc.fr),
                sparse: true,
                window_start: tpc.ductor.start_time,
                window_duration: tpc.ductor.readout_duration,
                tick: tpc.adc.tick,
                reference_time: 0.0,
                // Run wirecell-gen morse-* to find these numbers that match the extra
                // spread the sigproc induces.
            } + tpc.splat,
        }, nin=1, nout=1, uses=[tpc.anode, tpc.fr]);
        local reframer = pg.pnode({
            type: "Reframer",
            name: $.this_name(extra_name, '_splat'),
            data: {
                anode: wc.tn(tpc.anode),
                nticks: tpc.adc.readout_nticks,
            },
        }, nin=1, nout=1, uses=[tpc.anode]);
        pg.pipeline([splat, reframer]),

    
    // A crossline to resample and scale.  This is probably best applied across
    // "groups" (connected views).
    resample_crossline(ratio=1.0/4.0, scale=1.0, multiplicity=4, extra_name="")::
        pg.crossline([
            local this_name = $.this_name(extra_name, "_v"+std.toString(view));
            local res = pg.pnode({
                type:'SPNGResampler',
                name: this_name,
                data: {
                    ratio: ratio,
                } + control
            }, nin=1, nout=1);
            local scaler = pg.pnode({
                type:'SPNGTransform',
                name: this_name,
                data: {
                    operations: [
                        { operation: "scale", scalar: scale },
                    ],                
                } + control,
            }, nin=1, nout=1);
            pg.pipeline([res, scaler])
            for view in wc.iota(multiplicity)]),

    
    /// A subgraph to splat and apply resampling and output a frame.
    splat_frame(ratio=1.0/4.0, scale=1.0, extra_name="_splat")::
        local splat = $.splat(extra_name=extra_name);
        local unpack = $.tpc_group_unpacker(extra_name=extra_name);
        local ngroups = std.length(unpack.oports);
        local body = pg.shuntlines([
            unpack,
            $.resample_crossline(ratio=ratio, scale=scale, multiplicity=ngroups, extra_name=extra_name),
            $.tensor_packer(multiplicity=ngroups, extra_name=extra_name)
        ]);
        local resample = $.wrap_bypass(body, extra_name);
        pg.pipeline([
            splat,
            $.frame_to_tdm(extra_name=extra_name),
            resample,
            $.tdm_to_frame(extra_name=extra_name),
        ]),
        

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
        local name = $.this_name(extra_name, '_'+group.name);
        local plane_index = group.view_index;
        local fk = $.filter_kernel([tpc.filters[plane_index].channel.decon,{kind:"none"}],
                                   extra_name=extra_name);
        local rk = $.response_kernel(plane_index, extra_name=extra_name);
        $.decon_kernel(fk, rk, extra_name='_'+group.name+extra_name),

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
        local name = tpc.name + '_' + group.name + extra_name;
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
                          tag="decon", extra_name='_'+group.name+extra_name),

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




    /// Convert an IFrame to a TDM tensor set.
    frame_to_tdm(tag="", groups=null, extra_name=""):: pg.pnode({
        local the_groups = if std.type(groups) == "null"
                           then tpc.view_groups
                           else groups,
        type: 'SPNGFrameToTdm',
        name: tpc.name + extra_name,
        data: {
            anode: wc.tn(tpc.anode),
            rules: [{
                tag: tag,
                groups: [{
                    wpids: g.signed_wpids(tpc.ident),
                }, for g in the_groups]
            }],                 // just one rule
        } + control
    }, nin=1, nout=1, uses=[tpc.anode]),

    /// Convert a TDM tensor set to an IFrame.  This provides a simplified
    /// configuration where a single rule is used to match traces tensors.
    tdm_to_frame(extra_name="", traces_tag="signal", chid_tag="null"):: pg.pnode({
        type: 'SPNGTdmToFrame',
        name: tpc.name + extra_name,
        data: {
            // Locate the original frame object (just metadata)
            frame: {datapath: "/frames/\\d+/frame"},
            // Rules to locate tensor to include as tagged trace sets.
            tagged_traces: [ {
                // eg datapath of /frames/0/tags/gauss/groups/0/traces
                traces: { tag: traces_tag },
                // eg datapath of /frames/0/tags/null/rules/0/groups/0/chids
                chids: { tag: chid_tag },
            }],
            // chmasks: ...
            
        } + control
    }, nin=1, nout=1, uses=[]),


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
            name: $.this_name(extra_name, "_groups"),
            data: {
                selections: [{datapath:datapath_pattern % group} for group in wc.iota(ngroups) ],
            } + control
        }, nin=1, nout=ngroups),


    /// A 1->3 subgraph with TDM frame tensor set input and 3 per-view decon output.
    frame_decon(extra_name=""):: pg.shuntlines([
        $.tpc_group_unpacker(extra_name=extra_name),
        $.tpc_group_decon_views(extra_name=extra_name)]),


    /// An object providing pseudo sources for faning out from a 3-view pnode
    /// (eg decon) to multi-view nodes to use as sources for feeding initial_roi
    /// ("wiener" filter) and dnnroi_dense_Views ("dnnroi" filter).  This omits
    /// "gauss", see fanout_for_inference().  The dnnroi_views list which views
    /// provide dnnroi fodder.
    fanout_for_dnnroi_training(crossed_views=[1,1,0], extra_name=extra_name)::
        fans.fanout_select($.this_name(extra_name, "_dnnroi_training"),
                           std.length(crossed_views), targets_list=[
                               if is_crossed == 1
                               then ["wiener", "dense"]
                               else ["wiener"]
                               for is_crossed in crossed_views]),


    // Like fanout_for_training but include gauss to connect with apply_roi_views.
    fanout_for_dnnroi_inference(crossed_views=[1,1,0], extra_name=extra_name)::
        fans.fanout_select($.this_name(extra_name, "_dnnroi_inference"),
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
            $.downsample_views(rebin, views=views, extra_name='_tight'+extra_name),
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
            //$.downsample_views(rebin, views=views, extra_name='_dense'+extra_name),
            $.rebin_views(rebin, views=views, extra_name='_dense'+extra_name),
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
            } + control,
        }, nin=1, nout=1, uses=[tpc.anode]),
    
        
    tensor_packer(multiplicity=3, extra_name="")::
        pg.pnode({
            type: 'SPNGTensorPacker',
            name: $.this_name(extra_name, '_pack' + std.toString(multiplicity)),
            data: {
                multiplicity: multiplicity
            } + control
        }, nin=multiplicity, nout=1),


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
        local cellviews = $.cellviews_tensorset(out_views=out_views, chunk_size=chunk_size, extra_name=extra_name);
        local nout = std.length(out_views);
        local unpacker = pg.pnode({
            type: 'SPNGTorchSetUnpacker',
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
            operations: [{
                operation: operation,
                scalar: scalar,
                dims: dims,
            }]
        } + control,
    }, nin=1, nout=1) for view in views]),
        
    /// Return list of pnode that apply a 2->1 reduction per view.
    reduce_views_list(operation, dim=0, views=[0,1,2], extra_name=""):: [pg.pnode({
        type: 'SPNGReduce',
        name: $.this_name(extra_name, 'v'+std.toString(view) + '_'+operation),
        data: {
            operation: operation,
            dim: dim,
        } + control,
    }, nin=2, nout=1) for view in views],

    /// Return object with pseudo sinks/sources:
    /// - mp_sink :: an N-port sink of cellviews tensors (with feature dimension)
    /// - dense_sink :: an N-port sink of "dense' tensors (with no feature dimension)
    /// - source :: an N-port source of dense added to the feature dimension
    ///
    /// Note, this puts the dense as the first feature prior to the MPs.
    /// The order of the MPs depends on where they come from, eg CellViews.
    dnnroi_stack_features(views=[0,1], extra_name=""):: {
        local features = $.transform_views("unsqueeze", dims=[-3], views=views, extra_name=extra_name),
        local stacks = $.reduce_views_list("cat", dim=-3, views=views, extra_name=extra_name),
        dense_sink: features,
        source: pg.intern(outnodes=stacks, edges=[
            pg.edge(features, stacks[vi.index], vi.index, 0)
            for vi in wc.enumerate(views)]),
        mp_sink: pg.intern(iports=[s.iports[1] for s in stacks]),
    },        

    // Connect source mp and dense to the dnnroi stack operation.  Return pnode
    // that acts as a source of dnnroi input fodder.  The inputs are carried by
    // this object but the user must still handle connecting any iports of the
    // inputs.
    connect_dnnroi_stack(mp, dense, views=[0,1], extra_name="")::
        local sf = $.dnnroi_stack_features(views, extra_name);
        pg.intern(centernodes=[pg.shuntline(mp, sf.mp_sink), pg.shuntline(dense, sf.dense_sink)],
                  outnodes=[sf.source]),

    /// [1]tensor set -> tensor[2].
    /// Input is TDM frame tensor set.  Output is the 3-feature tensors for input to dnnroi training.
    ///
    dnnroi_training_preface(crossed_views = [1,1,0], rebin=4, extra_name="")::
        local sg1 = $.frame_decon(extra_name=extra_name);

        local decon_fan = $.fanout_for_dnnroi_training(crossed_views,extra_name=extra_name);
        local sg1_connection = pg.shuntline(sg1, decon_fan.sink);

        // [3]tensor -> tensor[3]
        local sg2 = $.tight_roi(rebin=rebin, extra_name=extra_name);

        local dnnroi_views = [vi.index for vi in wc.enumerate(crossed_views) if vi.value == 1];

        // [3]tensor -> tensor[2]
        local sg3 = $.cellviews_tensors(out_views=dnnroi_views, chunk_size=0, extra_name=extra_name);

        // [3]tensor -> tensor[2] combo of two above
        local sg23 = pg.shuntline(sg2, sg3);
        local sg23_connection = pg.shuntline(decon_fan.targets.wiener, sg23);

        // [2]tensor -> tensor[2]
        local sg4 = $.dnnroi_dense_views(views=dnnroi_views, rebin=rebin, extra_name=extra_name);
        local sg4_connection = pg.shuntline(decon_fan.targets.dense, sg4);

        // mp_sink:[3]tensor + dense_sink:[2]tensor(decon) -> source:tensor[2]
        local sg5 = $.connect_dnnroi_stack(sg3, sg4, views=dnnroi_views, extra_name=extra_name);

        pg.intern(innodes=[sg1], outnodes=[sg5],
                  centernodes=[sg1_connection, sg23_connection, sg4_connection]),

    /// Wrap each of source's oports with an "expand" operation on tensor
    /// dimension dim of size multiplicity.
    expand_dimension(source, dim=-3, multiplicity=3, operation="unbind", extra_name="")::
        local nports=std.length(source.oports);
        local unbinds = [ pg.pnode({
            type:'SPNGExpand',
            name: $.this_name(extra_name, 'v'+std.toString(ind)+'_'+operation),
            data: {
                operation: operation,
                dim: dim,
                multiplicity: multiplicity,
            } + control,
        }, nin=1, nout=multiplicity) for ind in wc.iota(nports)];
        pg.intern(innodes=[source],
                  outnodes=unbinds,
                  edges=[
                      pg.edge(source, unbinds[ind], ind, 0)
                      for ind in wc.iota(nports)
                  ]),

    /// An Nview->Nview applying the "gauss" filter and rebinning.
    ///
    /// Note, unlike wiener and dnnroi, we do not rebin by default.  With no
    /// rebin, ROIs must be unrebinned prior to being applied.
    gauss_dense_views(views=[0,1,2], rebin=1, extra_name="")::
        if rebin == 1 
        then $.time_filter_views("gauss", views=views, extra_name=extra_name)
        else pg.shuntlines([
            $.time_filter_views("gauss", views=views, extra_name=extra_name),
            $.downsample_views(rebin, views=views, extra_name=extra_name)]),

    // 
    dnnroi_inference_preface(crossed_views = [1,1,0], rebin=4, extra_name="")::
        local sg1 = $.frame_decon(extra_name=extra_name);

        // Give .sink and .targets.{gauss, wiener, dense}
        local decon_fan = $.fanout_for_dnnroi_inference(crossed_views,extra_name=extra_name);
        local sg1_cap = pg.shuntline(sg1, decon_fan.sink);

        // Must extract any rois NOT crossed so that we can expose them to caller.
        local roi_fan = fans.fanout_select($.this_name(extra_name, "_initrois"),
                                           std.length(crossed_views), targets_list=[
                                               if is_crossed == 1
                                               then ["shunt"]
                                               else ["shunt", "extract"]
                                               for is_crossed in crossed_views]);

        // Initial ROI block
        local sg2_block = pg.shuntlines([decon_fan.targets.wiener,
                                         $.tight_roi(extra_name=extra_name),
                                         roi_fan.sink]);


        local dnnroi_views = [vi.index for vi in wc.enumerate(crossed_views) if vi.value == 1];

        // [3]tensor -> tensor[2]
        local sg3 = $.cellviews_tensors(out_views=dnnroi_views, chunk_size=0, extra_name=extra_name);
        local sg3_feed = pg.shuntline(roi_fan.targets.shunt, sg3);

        // [2]tensor -> tensor[2]
        local sg4 = $.dnnroi_dense_views(views=dnnroi_views, rebin=rebin, extra_name=extra_name);
        local sg4_feed = pg.shuntline(decon_fan.targets.dense, sg4);

        // mp_sink:[3]tensor + dense_sink:[2]tensor(decon) -> source:tensor[2]
        local sg5 = $.connect_dnnroi_stack(sg3, sg4, views=dnnroi_views, extra_name=extra_name);

        local sg6 = pg.shuntline(decon_fan.targets.gauss,
                                 $.gauss_dense_views(views=wc.iota(std.length(crossed_views)),
                                                     rebin=1, extra_name=extra_name));
        {
            decon_sink: pg.intern(innodes=[sg1],
                                  centernodes=[
                                      sg1_cap,
                                      sg2_block,
                                      sg3_feed,
                                      sg4_feed,
                                      ]),
            fodder: sg5,
            gauss: sg6,
            rois: roi_fan.targets.extract
        },

    dnnroi_model(modelfile="unet-l23-cosmic500-e50.ts"):: {
        type:'SPNGTensorForwardTS',
        name:modelfile,     // explicitly do not use any context name to enable sharing.
        data: {
            ts_filename:modelfile,
        } + control
    },

    /// A single view dnnroi forward node.
    ///
    /// Note: OG DNNROI requires some seemingly arbitrary pre/post transforms.
    /// Future DNNs or if/when retraining DNNROI it would be nice to not require
    /// them.
    ///
    /// FIXME: There is likely a mismatch in the feature order between what is
    /// provided and what OG DNNROI expects.
    dnnroi_forward_view(model, view, extra_name="")::
        local prefix = "v"+std.toString(view)+"_dnnroi";
        local pre = pg.pnode({
            type: 'SPNGTransform',
            name: $.this_name(extra_name, prefix + "_pre"),
            data: {
                operations: [
                    { operation: "transpose", dims: [-2,-1] }
                ]
            } + control
        }, nin=1, nout=1);

        local fwd = pg.pnode({
            type: 'SPNGTensorForward',
            name: $.this_name(extra_name, prefix+"_fwd"),
            data: {
                forward: wc.tn(model),
                nbatch: 1,
            } + control
        }, nin=1, nout=1, uses=[model]);

        local post = pg.pnode({
            type: 'SPNGTransform',
            name: $.this_name(extra_name, prefix+"_post"),
            data: {
                operations: [
                    /// DNNROI traditionally works in (ntick,nchan) for last two
                    /// dims while WCT works in (nchan,ntick) so we must bookend
                    /// the forward with a transpose.
                    { operation: "transpose", dims: [-2,-1] },
                    // in principle could threshold here but see comments, next.

                    // We have to force nbatch=1 on TensorForward because DNNROI
                    // expects a batch dim.  But we don't actually want it on
                    // output so remove it now.  Also, DNNROI keeps the "feature
                    // dimension".
                    { operation: "squeeze", dims: [1, 0] },
                ],
            } + control
        }, nin=1, nout=1);

        // Convert DNNROI output to Boolean mask.
        local thresh_name = $.this_name(extra_name, prefix+'_roi');
        local thresh = pg.pnode({
            type: 'SPNGThreshold',
            name: thresh_name,
            data: {
                nominal: 0.5,   // FIXME: a study is needed to best set this
                tag: "roi",
                datapath_format: "/traces/Threshold/" + thresh_name,
            } + control
        }, nin=1, nout=1);

        pg.pipeline([pre, fwd, post, thresh]),

    /// Ncrossed->Ncrossed Make DNNROI-forward subgraph.  Input is initial
    /// dnnroi "fodder" tensors with three features.  Output are DNNROI ROI
    /// tensors.
    dnnroi_forward_views(modelfile="unet-l23-cosmic500-e50.ts", crossed_views=[1,1,0], extra_name="")::
        local model = $.dnnroi_model(modelfile);
        pg.crossline([$.dnnroi_forward_view(model, view=vi.index, extra_name=extra_name)
                      for vi in wc.enumerate(crossed_views)
                      if vi.value == 1]),

    /// 3->3 do interval space unrebin by given factor.
    unbin_views(rebin=4, views=[0,1,2], tag="", extra_name="")::
        pg.crossline([
            pg.pnode({
                local this_name = $.this_name(extra_name, 'v'+std.toString(view)+'_unbin'+std.toString(rebin)),
                type: "SPNGRebinner",
                name: this_name,
                data: {
                    norm: "maximum",
                    factor: -rebin,     // FIXME: must match upstream resampler.
                    tag: tag,
                    datapath_format: "/traces/Rebinner/" + this_name
                } + control
            }, nin=1, nout=1) for view in views]),
            
    rebin_views(rebin=4, views=[0,1,2], tag="", extra_name="")::
        pg.crossline([
            pg.pnode({
                local this_name = $.this_name(extra_name, 'v'+std.toString(view)+'_rebin'+std.toString(rebin)),
                type: "SPNGRebinner",
                name: this_name,
                data: {
                    norm: "maximum",
                    factor: rebin,     // FIXME: must match upstream resampler.
                    tag: tag,
                    datapath_format: "/traces/Rebinner/" + this_name
                } + control
            }, nin=1, nout=1) for view in views]),
            


    /// Return pnode for single view that applies ROIs.  It takes in a dense (eg
    /// gauss) on iport=0 and roi on iport=1 and outputs signal on its only
    /// oport=0.
    applyroi_view(view, extra_name="")::
        local name = $.this_name(extra_name, "v"+std.toString(view)+"_applyroi");
        local mul = pg.pnode({
            type: 'SPNGReduce',
            name: name,
            data: {
                operation: 'mul',
            } + control
        }, nin=2, nout=1);
        local rbl = pg.pnode({
            type: 'SPNGRebaseliner',
            name: name,
            data: {
                tag: "signal",
                datapath_format: "/traces/Rebaseliner/" + name
            } + control
        }, nin=1, nout=1);
        pg.pipeline([mul, rbl]),

    /// Return object with .roi_sink, .dense_sink and .signal_source with apply
    /// ROIs in the middle.  All inputs/outputs are size that of views.
    applyroi_views(views=[0,1,2], extra_name="")::
        local perview = [$.applyroi_view(view, extra_name) for view in views];
        {
            dense_sink: pg.intern(iports=[pv.iports[0] for pv in perview]),
            roi_sink: pg.intern(iports=[pv.iports[1] for pv in perview]),
            signal_source: pg.intern(outnodes=perview),
        },


    /// [3]tensor -> tensor[3]
    ///
    /// Input decon, output signals using cell basis initial rois and DNNROI
    /// inference for crossed views.  Intermediate waveforms are rebinned.
    dnnroi_inference(modelfile="unet-l23-cosmic500-e50", crossed_views=[1,1,0], rebin=4, extra_name="")::
        local all_views = wc.iota(std.length(crossed_views));
        
        // For inference, sg1_infer is just the beginning.  Must feed .fodder to
        // dnnroi forward and then that plus .rois plus .gauss into apply_roi.
        local sg1_infer = $.dnnroi_inference_preface(rebin=rebin, extra_name="_INFER");

        local dnnroifwd = $.dnnroi_forward_views(modelfile=modelfile, 
                                                  crossed_views=crossed_views, extra_name="_FORWARD");
        local sg1_dnnroi = pg.shuntline(sg1_infer.fodder, dnnroifwd);

        // Merge all sources of ROIs.
        ///
        // BIG FAT WARNING: this should be done with knowledge of crossed_views to
        // get the port ordering correct.  As is, it assumes dnnroi rois all come
        // before the views that have only initial rois.
        local fat_rois = pg.crossline([sg1_dnnroi, sg1_infer.rois]);
        local unbin = $.unbin_views(rebin=rebin, views=all_views, extra_name="_ROIS");
        local rois = pg.shuntline(fat_rois, unbin);

        /// Gives .roi_sink, .dense_sink and .signal_source
        local applyrois = $.applyroi_views(views=all_views, extra_name="_APPLYROIS");
        local rois_cap = pg.shuntline(rois, applyrois.roi_sink);
        local dense_cap = pg.shuntline(sg1_infer.gauss, applyrois.dense_sink);
        pg.intern(innodes=[sg1_infer.decon_sink],
                  outnodes=[applyrois.signal_source],
                  centernodes=[rois_cap, dense_cap]),


    // Return object that keeps source's iports and caps off its oports by
    // attaching a file sink saving in WCT "tensor file" format..
    attach_tensor_file_sink_views(source, filename, extra_name="", prefix="")::
        local mult = std.length(source.oports);
        local this_name = tpc.name + extra_name;
        local pack = pg.pnode({
            type: 'SPNGTensorPacker',
            name: this_name,
            data: {
                multiplicity: mult,
            } + control
        }, nin=mult, nout=1);
        local ttt = pg.pnode({
            type: 'TorchToTensor',
            name: this_name,
            data: {},
        }, nin=1, nout=1);
        local sink = pg.pnode({
            type: 'TensorFileSink',
            name: this_name,
            data: {
                outname: filename,
                prefix: prefix,
            },
        }, nin=1, nout=0);
        pg.intern(innodes=[source],
                  centernodes=[pack, ttt, sink],
                  edges=[
                      pg.edge(source, pack, port.index, port.index)
                      for port in wc.enumerate(source.oports)
                  ] + [
                      pg.edge(pack, ttt),
                      pg.edge(ttt, sink)
                  ]),

    // Like attach_file_sink_views() but source's oports are fanned out with one
    // fan going to a sink and the other exposing the oports.
    attach_tensor_file_tap_views(source, filename, extra_name="", prefix="")::
        local fos = fans.fanout_shuntline(std.length(source.oports), nout=2, extra_name='_tap'+extra_name);
        local sink = $.attach_tensor_file_sink_views(fos.sources[1], filename, extra_name, prefix);
        local head = pg.shuntline(source, fos.sink);
        pg.intern(innodes=[head], outnodes=[fos.sources[0]], centernodes=[source, sink, head]),


    /**
     * Wrap a graph in a bypass.
     *
     * The input frame tensor set is sent into body and bypasses the body graph.
     * The output tensor set from the body is merged into the bypassed tensor set and provides the final output.
     *
     * @param body A subgraph that is bookended with a frame tensor set sink and source.
     * @return A subgraph with same bookends.
     * 
     */
    wrap_bypass(body, extra_name="")::
        local name = $.this_name(extra_name, '_bypass');
        local bypass = frame.bypass(name);
        pg.intern(iports=[bypass.iports[0]],
                  oports=[bypass.oports[1]],
                  centernodes=[body, bypass],
                  edges=[
                      pg.edge(bypass, body),
                      pg.edge(body, bypass, 0, 1),
                  ]),
}


