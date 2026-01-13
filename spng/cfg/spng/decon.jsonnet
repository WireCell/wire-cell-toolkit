local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

function(tpc, control={})
{
    /// A kernel providing a Fourier-space filter.
    filter_kernel(name,         // unique name
                  axis_config,  // list of per-axis config objects
                  /// See DeconKernel for following options
                  debug_filename="")::
        {
            type: "SPNGFilterKernel",
            name: name,
            data: {
                axis: axis_config,
                debug_filename: debug_filename,
            } + control,
        },

    /// A kernel providing a response built from FR*ER
    response_kernel(plane_index, // the index of the plane in the FR to use
                    extra_name="",
                    debug_filename="")::
        {
            type: "SPNGResponseKernel",
            name: tpc.name + "v" + std.toString(plane_index) + extra_name,
            data: {
                field_response: wc.tn(tpc.fr),
                elec_response: wc.tn(tpc.er),
                // Take care that this should merely cover where ER is non-zero
                elec_duration: tpc.er.data.nticks*tpc.er.data.tick,
                period: tpc.adc.tick, // sample period of kernel time dimension.
                plane_index: plane_index,
                // Scale to units of ADC
                scale: -1 / tpc.adc.lsb_voltage,
                debug_filename: debug_filename
            } + control,
            uses: [tpc.fr, tpc.er],
        },


    /// A deconvolution kernel.
    decon_kernel(name,          // unique name
                 filt,          // a kernel providing the filter (numerator)
                 resp,          // a kernel providing the response (denominator)
                 /// See DeconKernel for following options:
                 debug_filename="")::
        {
            type: "SPNGDeconKernel",
            name: name,
            data: {
                filter: wc.tn(filt),
                response: wc.tn(resp),
                debug_filename: debug_filename,
            } + control,
            uses: [filt, resp]
        },
    

    /// Return a decon kernel for one tpc and group index.  If this decon kernel
    /// is to be used for this group in all TPCs, the extra name can remain
    /// empty.  If unique decon kernels are needed for this group index across
    /// many TPCs, give a unique extra name.
    tpc_group_decon_kernel(group_index, extra_name="")::
        local group = tpc.view_groups[group_index];
        local name = group.name + extra_name;
        local plane_index = group.view_index;
        local fk = $.filter_kernel(name, [tpc.filters[plane_index].channel.decon,{kind:"none"}]);
        local rk = $.response_kernel(plane_index, extra_name=name);
        $.decon_kernel(name, fk, rk),


    /// Configure node to convolve intput data with a given kernel
    kernel_convolve(name,       // unique name
                    kernel,     // the object providing the convolution kernel 
                    axes_config, // list of per-axis config objects
                    // See KernelConvolve for following options
                    tag="",
                    datapath_format="",
                    debug_filename="")::
        pg.pnode({
            type: "SPNGKernelConvolve",
            name: name,
            data: {
                kernel: wc.tn(kernel),
                axis: axes_config,
                tag: tag,
                datapath_format: datapath_format,
                faster: true,       // use "faster DFT size"
                debug_filename: debug_filename,
            } + control
        }, nin=1, nout=1, uses=[kernel]),


    /// Configuration for one FR*ER decon node for one group in a TPC.  The node
    /// will be named uniquely for that context.  If more than one decon in the
    /// same tpc+group, pass extra_name.
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
        $.kernel_convolve(name, kernel, [co_channel, co_time],
                          datapath_format='/traces/group/' + std.toString(group_index) + "/KernelConvolve/" + name,
                          tag="decon"),

    
    /// Return nodes[0] if length one else return a subgraph with a fanin that stacks the tensors.
    group_fanin(nodes)::
        local n = std.length(nodes);
        local fanin = pg.pnode({
            type: 'SPNGReduce',
            name: std.join('-', [node.name for node in nodes]),
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
    collect_groups_by_view(nodes)::
        [$.group_fanin([
            nodes[gi.index]
            for gi in wc.enumerate(tpc.view_groups)
            if gi.value.view_index == view_index])
         for view_index in [0,1,2]],

    /// Produce an ngroup->nview subgraph doing FR*ER decon.
    group_decon_view:
        local groups = tpc.view_groups;
        local ngroups = std.length(groups);

        local decon_kernels = [$.tpc_group_decon_kernel(group_index, extra_name="")
                               for group_index in wc.iota(ngroups)];
        local decon_frers = [$.tpc_group_decon_frer(it.index,
                                                    decon_kernels[it.value.view_index])
                             for it in wc.enumerate(groups)
                            ];
        local view_frers = $.collect_groups_by_view(decon_frers);
        pg.intern(innodes=decon_frers, outnodes=view_frers),


    /// Produce an ngroup->nview subgraph doing FR*ER decon.
    group_decon_simple:
        local groups = tpc.view_groups;
        local ngroups = std.length(groups);

        local decon_kernels = [$.tpc_group_decon_kernel(group_index, extra_name="")
                               for group_index in wc.iota(ngroups)];
        local decon_frers = [$.tpc_group_decon_frer(it.index,
                                                    decon_kernels[it.value.view_index])
                             for it in wc.enumerate(groups)
                            ];
        local view_frers = $.group_fanin(decon_frers);
        pg.intern(innodes=decon_frers, outnodes=[view_frers]),

    /// One filter along time dimension
    time_filter_one(filter, view_index, extra_name="")::
        local vis = std.toString(view_index);
        local name = tpc.name + 'v' + vis + extra_name;
        local fk = $.filter_kernel(name, [{kind:"none"}, tpc.filters[view_index].time[filter]]);
        // local co = $.convo_options(0, view_index);
        local co_channel = {padding: "none", dft: false};
        local co_time = {padding: "none", dft: true} + tpc.filters[view_index].time.options;
        $.kernel_convolve(name, fk, [co_channel, co_time],
                          datapath_format='/traces/view/' + vis + "/KernelConvolve/" + name,
                          tag=filter),


    /// Apply a time filter across each of the views for a 3->3 subgraph.
    time_filter_views(filter="gauss", views=[0,1,2])::
        local nodes = [
            $.time_filter_one(filter, view_index, extra_name=filter)
            for view_index in views];
        pg.crossline(nodes),
        
    


}
