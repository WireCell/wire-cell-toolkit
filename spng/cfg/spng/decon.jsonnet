local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

local converters = import "converters.jsonnet";

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
            },
        },

    /// A kernel providing a response built from FR*ER
    response_kernel(tpc,
                    plane_index, // the index of the plane in the FR to use
                    extra_name="",
                    debug_filename="")::
        {
            type: "SPNGResponseKernel",
            name: tpc.name + extra_name,
            data: {
                field_response: wc.tn(tpc.fr),
                elec_response: wc.tn(tpc.er),
                // Take care that this should merely cover where ER is non-zero
                elec_duration: tpc.er.data.nticks*tpc.er.data.tick,
                period: tpc.adc.tick, // sample period of kernel time dimension.
                plane_index: plane_index,
                // Scale to units of ADC and negate to give positive signals.
                scale: -1 / tpc.adc.vin_per_count,
                debug_filename: debug_filename
            },
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
            },
            uses: [filt, resp]
        },
    

    /// Return a decon kernel for one tpc and group index.  If this decon kernel
    /// is to be used for this group in all TPCs, the extra name can remain
    /// empty.  If unique decon kernels are needed for this group index across
    /// many TPCs, give a unique extra name.
    tpc_group_decon_kernel(tpc, group_index, extra_name="")::
        local name = "g" + std.toString(group_index) + extra_name;
        local plane_index = converters.tpc_group_planes(tpc)[group_index];
        local fk = $.filter_kernel(name, [tpc.filters[plane_index].channel.decon,{kind:"none"}]);
        local rk = $.response_kernel(tpc, plane_index, extra_name=name);
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
            },
        }, nin=1, nout=1, uses=[kernel]),




    

}
