local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

{
    /// A kernel providing a Fourier-space filter.
    filter_kernel(name,         // unique name
                  axis_config,  // list of per-axis config objects
                  /// See DeconKernel for following options
                  debug_filename="",
                  control={})::
        {
            type: "SPNGFilterKernel",
            name: name,
            data: {
                axis: axis_config,
                debug_filename: debug_filename,
            } + wc.object_with(control, ["verbosity", "device"]),
        },

    /// A kernel providing a response built from FR*ER
    response_kernel(tpc,
                    plane_index, // the index of the plane in the FR to use
                    extra_name="",
                    debug_filename="",
                    control={})::
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
            } + wc.object_with(control, ["verbosity", "device"]),
            uses: [tpc.fr, tpc.er],
        },


    /// A deconvolution kernel.
    decon_kernel(name,          // unique name
                 filt,          // a kernel providing the filter (numerator)
                 resp,          // a kernel providing the response (denominator)
                 /// See DeconKernel for following options:
                 debug_filename="",
                 control={})::
        {
            type: "SPNGDeconKernel",
            name: name,
            data: {
                filter: wc.tn(filt),
                response: wc.tn(resp),
                debug_filename: debug_filename,
            } + wc.object_with(control, ["verbosity", "device"]),
            uses: [filt, resp]
        },
    

    /// Return a decon kernel for one tpc and group index.  If this decon kernel
    /// is to be used for this group in all TPCs, the extra name can remain
    /// empty.  If unique decon kernels are needed for this group index across
    /// many TPCs, give a unique extra name.
    tpc_group_decon_kernel(tpc, group_index, extra_name="", control={})::
        local group = tpc.view_groups[group_index];
        local name = group.name + extra_name;
        local plane_index = group.view_index;
        local fk = $.filter_kernel(name, [tpc.filters[plane_index].channel.decon,{kind:"none"}], control=control);
        local rk = $.response_kernel(tpc, plane_index, extra_name=name, control=control);
        $.decon_kernel(name, fk, rk, control=control),


    /// Configure node to convolve intput data with a given kernel
    kernel_convolve(name,       // unique name
                    kernel,     // the object providing the convolution kernel 
                    axes_config, // list of per-axis config objects
                    // See KernelConvolve for following options
                    tag="",
                    datapath_format="",
                    debug_filename="",
                    control={})::
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
            } + wc.object_with(control, ["verbosity", "device"]),
        }, nin=1, nout=1, uses=[kernel]),

    /// Determine how to convolve each axis based on channel connection and if
    /// streaming mode is used and for each kind of convolution.
    convo_options(connection, streaming=false):: {
        local decon_channel = [
            {cyclic: false, crop: -2},  // one plane
            {cyclic: false, crop: -2},  // two concatenated planes
            {cyclic: true,  crop:  0}, // two wrapped planes
        ][connection],
        local decon_time =
            if streaming
            then {cyclic: false, crop: 0, baseline: true/*, roll_mode: "decon"*/}
            else {cyclic: false, crop: 0, baseline: true/*, roll_mode: "decon"*/},

        local nothing = {padding: "none", dft: false},
        local just_filter = {padding: "none", dft: true},

        // For FR*ER deconvolution, regardless of filters.
        decon_resp: [decon_channel, decon_time],
        // For follow-on of a time filter, eg "gauss".
        filter_time: [nothing, just_filter],
        // For future RC decon
        decon_time: [nothing, decon_time],

        nothing: nothing,
        just_filter: just_filter,
    },



    /// Configuration for one FR*ER decon node for one group in a TPC.  The node
    /// will be named uniquely for that context.  If more than one decon in the
    /// same tpc+group, pass extra_name.
    tpc_group_decon_frer(tpc, group_index, kernel, tag="", extra_name="", streaming=false, control={})::
        local group = tpc.view_groups[group_index];
        local name = tpc.name + group.name + extra_name;
        local co = $.convo_options(group.connection, streaming);
        $.kernel_convolve(name, kernel, co.decon_resp, tag="", control=control),

    
    /// Return nodes[0] if length one else return a subgraph with a fanin that stacks the tensors.
    group_fanin(nodes, control={})::
        local n = std.length(nodes);
        local fanin = pg.pnode({
            type: 'SPNGReduce',
            name: std.join('-', [node.name for node in nodes]),
            data: {
                multiplicity: n,
                operation: "cat",
                dim: -2         // concatenate along channel dimension
            } + wc.object_with(control, ["verbosity", "device"]),
        }, nin=n, nout=1);
        
        if n == 1
        then nodes[0]
        else pg.intern(innodes=nodes, outnodes=[fanin],
                       edges=[pg.edge(it.value, fanin, 0, it.index) for it in wc.enumerate(nodes)]),
            


    /// Given array of nodes that span view groups, return array that span views
    collect_groups_by_view(tpc, nodes, control={})::
        [$.group_fanin([
            nodes[gi.index]
            for gi in wc.enumerate(tpc.view_groups)
            if gi.value.view_index == view_index], control=control)
         for view_index in [0,1,2]],

    /// Produce an ngroup->nview subgraph doing FR*ER decon.
    group_decon_view(tpc, control={})::
        local groups = tpc.view_groups;
        local ngroups = std.length(groups);
        //local response_kernels = [$.response_kernel(tpc, view_index, control=control) for view_index in [0,1,2]];

        local decon_kernels = [$.tpc_group_decon_kernel(tpc, group_index, extra_name="", control=control)
                               for group_index in wc.iota(ngroups)];
        local decon_frers = [$.tpc_group_decon_frer(tpc, it.index,
                                                    decon_kernels[it.value.view_index], control=control)
                             for it in wc.enumerate(groups)
                            ];
        local view_frers = $.collect_groups_by_view(tpc, decon_frers, control=control);
        pg.intern(innodes=decon_frers, outnodes=view_frers),
                             
    /// One filter along time dimension
    time_filter_one(tpc, filter, view_index, extra_name="", streaming=false, control={})::
        local name = tpc.name + 'v' + std.toString(view_index) + extra_name;
        local fk = $.filter_kernel(name, [{kind:"none"}, tpc.filters[view_index].time[filter]], control=control);
        local co = $.convo_options(0, streaming);
        $.kernel_convolve(name, fk, co.filter_time, tag="", control=control),


    /// Apply a time filter across each of the views for a 3->3 subgraph.
    time_filter_views(tpc, filter="gauss", views=[0,1,2], streaming=false, control={})::
        local nodes = [
            $.time_filter_one(tpc, filter, view_index, extra_name=filter, streaming=streaming, control=control)
            for view_index in views];
        pg.crossline(nodes),
        
    


}
