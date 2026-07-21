/// A test job using the TDM retrofitted nodes.
///
/// This supports whatever detectors that omni supports.
///
/// TLAs:
///
/// input :: a depos file, optional, default to the muon-depos.npz file
/// output :: optional, if given must have a "%s" to name the tier.
/// anodeid :: optional, default is 3
/// device :: optional, default to cpu


local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local omnimap = import "omni/all.jsonnet";
local io = import "fileio.jsonnet";
local detectors = import "detectors.jsonnet"; // in wire-cell-data

/// Make a sink of frames through a body.  If output is given, tap a frame file sink
local frame_sink(body, output="", name="", tier="sig", digitize=false) = 
    local sink = if output == ""
                 then pg.pnode({ type: "DumpFrames", name: name }, nin=1, nout=0)
                 else io.frame_file_sink(output%tier, digitize=digitize);
    pg.pipeline([body, sink]);


local frame_tap(filename, digitize=true) =
    local frame_fout = pg.pnode({
        type: 'FrameFanout',
        name: filename,
        data: {
            multiplicity: 2
        }
    }, nin=1, nout=2);
    local sink = io.frame_file_sink(filename, digitize=digitize);
    pg.intern(innodes=[frame_fout],
              centernodes=[sink],
              outnodes=[frame_fout],
              edges=[pg.edge(frame_fout, sink, 1, 0)]);


/// Return a single sink that is actually a fanout to two sink
local frame_sink_two(one, two, name="") =
    local fout = pg.pnode({
        type: 'FrameFanout',
        name: name,
        data: {
            multiplicity: 2,
        },
    }, nin=1, nout=2);
    pg.intern(innodes=[fout], outnodes=[one, two], edges=[
        pg.edge(fout, one, 0, 0),
        pg.edge(fout, two, 1, 0)]);


// Define a "matrix" lookup table between a plane label and all labels.  Note,
// below we have 2 "groups" for W plane, but "plane" continues to mean just the
// three.
local plane_ids = {
    u: [0, "u", "U", "ind"],
    v: [1, "v", "V", "ind"],
    w: [2, "w", "W", "col"], 
    "0": self.u,
    "1": self.v,
    "2": self.w,
    ind: self.u + self.v,
    col: self.w,

    is(pid, key) :: std.member(self[std.toString(key)], pid),
};

/// Return string suitable for including in a configuration object .name to make
/// it unique to an anode/plane context.  In configurations where a distinction
/// across anodes is not required, do not pass anode.  Otherwise, pass anode
/// object or other anode identifier.
local plane_context_name(plane, anode=null) =
    if std.type(anode) == "null" then
        "p"+std.toString(plane)
    else if std.type(anode) == "object" then
        "a"+std.toString(anode.data.ident)+"p"+std.toString(plane)
    else
        "a"+std.toString(anode)+"p"+std.toString(plane);

// Make a generic 1D filter config.
local filter(scale, power=2, kind="lowpass", period=500*wc.ns, ignore_baseline=true) = {
    scale: scale,
    power: power,
    kind: kind,
    period: period,
    ignore_baseline: ignore_baseline
};

// This defines canonical 1D config parameter sets each with a hard-wired
// interpretation for signal processing deconvolution.
//
// This object is independent from the implementation of the components that
// perform decon.  It may be used to make OSP or SPNG config.
//
// In cases where the job does not distinquish different filters for different
// anodes, the anode will be null.  O.w. it may be an anode config object.
// 
// FIXME: this info needs to move to omni.  Here we gives a generic function
// implemented for PSP.
local decon_filters(plane_index, anode=null) = {

    plane_index: plane_index,
    anode: anode,

    /// This current implementation does not care about "anode".  When
    /// generalized to other detectors (eg PDHD with "bad" APA1) the "anode" may
    /// be used to distinquish different filter configs.

    // Channel filter applied as part of decon.
    //
    // Channel dimension from cfg/pgrapher/experiment/pdsp/sp-filters.jsonnet.
    // There it is called "wire filter" and a "wf()" function defines defaults
    // and separate filters for induction and collection are defined.  The
    // "sigma" there is the "scale" here.  Note, the filter baseline (DC zero
    // frequency) is preserved.  Also, note the "period" is a unitless 1.0 which
    // means one channel per sample.
    channel:
        if plane_ids.is(plane_index, "ind")
        then filter(scale=1.0 / wc.sqrtpi * 0.75, period=1.0, ignore_baseline=false)  // ind
        else filter(scale=1.0 / wc.sqrtpi * 10.0, period=1.0, ignore_baseline=false), // col
    // original pdsp has 0.75 and 10.0

    // Weiner (inspired) filter.
    //
    // This is one of two conventional time filters applied as part of decon.
    // It maximizes signal/noise.
    // 
    // For just decon, time dimension has two possibilities, the "Gauss" and the
    // "Wiener" (inspired) filters.  This also comes from sp-filters.jsonnet
    // where "tight" and "wide" variants of both are defined.  OSP hard-wires
    // the instance names "Wiener_tight_{U,V,W}".  Bad for user, good for
    // archaeology.
    wiener:
        if plane_ids.is(plane_index, "u") then
            filter(scale=0.1487880 * wc.megahertz, power=3.76194)
        else if plane_ids.is(plane_index, "v") then
            filter(scale=0.1596568 * wc.megahertz, power=4.36125)
        else                        // w
            filter(scale=0.1362300 * wc.megahertz, power=3.35324),

    // Gauss filter.
    //
    // This is the other conventional time filter applied as part of decon.  It
    // preserves signal.
    //
    // OSP lets the name for the Gauss filter be user defined with the default
    // "Gaus_wide" (sic, note, the typo!).  OSP weirdly devotes a dozen lines to
    // per-plane settings which are guaranteed identical so we only define a
    // single Gauss filter.
    gauss: filter(scale=0.12 * wc.megahertz),

    /// The identity filter config.
    flat_ac: { kind: "flat", ignore_baseline: true }, // AC coupled, DC is zero
    flat_dc: { kind: "flat", ignore_baseline: false }, // DC is 1.0

    /// Like identity but rely on unsqueeze to populate. 
    none: { kind: "none", },

    /// fixme: add more to cover ROI needs
};

local filter_kernel(name, axis) = {

    type: "SPNGFilterKernel",
    name: name,
    data: {
        axis: axis, 
        debug_filename: "filter-kernel-%s.pkl" % name
    },
    // no uses
};

local response_kernel(plane_index, fr, er, adc, name="") = {
    type: "SPNGResponseKernel",
    name: fr.name+er.name + name,
    data: {
        field_response: wc.tn(fr),
        elec_response: wc.tn(er),
        elec_duration: 20*wc.us, // just enough to cover non-zero value.
        period: adc.tick,        // sample period of kernel time dimension.
        plane_index: plane_index,
        // Negate to give positive signals.
        scale: -1 * adc.gain * (1 << adc.resolution) / (adc.fullscale[1] - adc.fullscale[0]),
        debug_filename: "response-kernel-%s.pkl" % name,
    },
    uses: [fr, er],         // NOT anode.
};


/// Bind an instance of the above generic filter configuration to
/// a specific implementation of the components.
///
local decon_kernel(filt, resp, name) = {
    type: "SPNGDeconKernel",
    name: name,
    data: {
        filter: wc.tn(filt),
        response: wc.tn(resp),
        debug_filename: "decon-kernel-%s.pkl" % name
    },
    uses: [filt, resp]
};

local convo_options(plane_index) = {
    decon_channel: [
        {cyclic: true,  crop:  0}, // u, cycle, wrapped
        {cyclic: true,  crop:  0}, // v, cyclic, wrapped
        {cyclic: false, crop: -2}  // w, linear
    ][plane_index],
    decon_time: {cyclic: false, crop: 0, baseline: true, roll_mode: "decon"},

    nothing: {padding: "none", dft: false},
    just_filter: {padding: "none", dft: true},

    axis: {
        // For FR*ER deconvolution, regardless of filters.
        decon_resp: [$.decon_channel, $.decon_time],
        // For follow-on of a time filter, eg "gauss".
        filter_time: [$.nothing, $.just_filter],
        // For future RC decon
        decon_time: [$.nothing, $.decon_time],
    },
};

/// Return config object for a KernelConvovle.  
local convo_node(kernel, axis, tag, name="") =

    // one decon per filter type
    pg.pnode({
        type: "SPNGKernelConvolve",
        name: name,
        data: {
            kernel: wc.tn(kernel),
            axis: axis,
            tag: tag,
            datapath_format: "/frames/{ident}/tags/{tag}/groups/{group}/traces",
            faster: true,       // use "faster DFT size"
            debug_filename: "kernel-convolve%s-{ident}.pkl"%name,
        },
    }, nin=1, nout=1, uses=[kernel]);


/// Return a compound node representing a source that simulates using input depos
local sim_source_node(depos_filename, anode, omni, nticks=6000) = 
              
    local depo_source = io.depo_file_source(depos_filename);

    # FIXME: remove omni and make this file stand-alone.
    local drift = omni.drift(anode);
    local signal = omni.signal(anode);
    local noise = omni.noise(anode);
    local digitize = omni.digitize(anode);

    // This is what sets the ADC waveform readout window.
    local reframer = pg.pnode({
        type: 'Reframer',
        name: anode.data.ident,
        data: {
            anode: wc.tn(anode),
            nticks: nticks,
        },
    }, nin=1, nout=1, uses=[anode]);

    pg.pipeline([depo_source, drift, signal, reframer, noise, digitize]);



/// Return a "plane group" object.  The "layer" MUST be only one.
local plane_group(layer, faces, anode, nchan) = {
    wpids: [wc.WirePlaneId(layer, face, anode.data.ident) for face in faces],
    plane: wc.wpid_layer_to_index(layer),
    nchan: nchan,
};


/// Return an array of "plane group" objects that determines how a frame is
/// split into parallel-processed tensors.
local plane_groups(detector, anode) =
    // fixme: some day, this should become dependent on "detector".  For now, we assume APA.
    local objs = [
        plane_group(wc.Ulayer, [0,1], anode, 800),
        plane_group(wc.Vlayer, [0,1], anode, 800),
        plane_group(wc.Wlayer, [0], anode, 480),
        plane_group(wc.Wlayer, [1], anode, 480),
    ];
    {
        objects: objs,
        wpids: [o.wpids for o in objs],
        planes: [o.plane for o in objs],
        nchans: [o.nchan for o in objs],
        count: std.length(objs),
    };


local adc_params(detector, anode=null) =
    // FIXME: make dependent on detector name and potentially anode.  For now,
    // it's hard-wired to return APA info.  And we use same object as in OG
    // params.adc
    {

        // The sample period of the ADC.  This assumes we live in a resampled
        // universe and not the 512ns reality of PDHD and part of PDVD.
        tick: 500*wc.ns,

        // The number of samples in an ADC readout.  Note, in general, this
        // number should NOT be used anywhere else.  It is NOT a size of any
        // convolution.  Specifically, it is NOT the size of the ER.
        nticks: 6000,

        // A relative, unitless gain.  Not ER gain.
        gain: 1.0,

        // The resolution (bits) of the ADC.  Use "modern" number (PDHD, FDHD).
        // FIXME: we probably still have 12bits in the sim part.
        resolution: 14,

        // The voltage range as [min,max] of the ADC, eg min voltage
        // counts 0 ADC, max counts 2^resolution-1.
        fullscale: [0.2*wc.volt, 1.6*wc.volt],
    };


// Return fr, er, rc responess for the detector.
local responses(detector, anode, adc_tick=500*wc.ns) = {
    local det = detectors[detector],

    fr : {
        type: 'FieldResponse',
        name: detector,
        data: {
            filename: det.fields[0], // "dune-garfield-1d565.json.bz2",
        },
    },

    // Fixme: un-hard-wire,  this is for APA
    er : {
        type: 'ColdElecResponse',
        name: detector,
        data: {
            // note, the default sampling (tick/nticks) is ignored so we don't
            // bother configuring it.
            shaping: 2.2*wc.us,
            gain : 14.0*wc.mV/wc.fC,
            postgain: 1.0,
        },
    },

    // Fixme: un-hard-wire,  this is for APA
    local rc = {
        type: 'RCResponse',
        name: detector,
        data: {
            width: 1.1*wc.ms,
        }
    },

    rc: [rc,rc],
};


/// A frame->(decon)->frame compound node that runs the original SPNG decon.
local orig_decon_node(detector, anode, adc, device="cpu") =
    local res = responses(detector, anode, adc.tick);

    local groups = plane_groups(detector, anode);

    local fanout = pg.pnode({
        type: 'FrameToTorchSetFanout',
        name: "fanout",
        data: {
            anode: wc.tn(anode),
            expected_nticks: adc.nticks,
            output_groups: groups.wpids,
            unsqueeze_output: true,
        },
    }, nin=1, nout=groups.count, uses=[anode]);

    // helper to make a per group (iplane) pieline
    local decon_pipeline(iplane) = 

        local torch_frer = {
            type: "TorchFRERSpectrum",
            name: "torch_frer%d_plane%d" % [anode.data.ident, iplane],
            uses: [res.fr, res.er],
            data: {
                field_response: wc.tn(res.fr),
                elec_response: wc.tn(res.er),
                fr_plane_id: if iplane > 2 then 2 else iplane,
                ADC_mV: 11702142857.142859,
                gain: 1.0,
                default_nchans : groups.nchans[iplane],
                default_nticks: 6000,
                readout_period: 500.0, #512.0,
                extra_scale: 1.0,
                debug_force_cpu: device == "cpu",
                
            }
        };
        local wire_filter_ind = if iplane < 2 then 0 else 1;
        local the_wire_filter = 
            local wf = {
                type: 'HfFilter',
                name: ['Wire_ind', 'Wire_col'][wire_filter_ind],
                data: {
                    max_freq: 1,  // warning: units
                    power: 2,
                    flag: false,
                    sigma: [ 1.0 / wc.sqrtpi * 0.75,  1.0 / wc.sqrtpi * 10.0][wire_filter_ind]
                }
            };
            {
                type: "Torch1DSpectrum",
                name: "orig-wire-filter%d" % wire_filter_ind,
                uses: [wf],
                data: {
                    spectra: [ wc.tn(wf) ],
                    device: device,
                },
            };

        local decon = pg.pnode({
            type: 'SPNGDecon',
            name: 'spng_decon_apa%d_plane%d' % [anode.data.ident, iplane],
            data: {
                frer_spectrum: wc.tn(torch_frer),
                wire_filter: wc.tn(the_wire_filter), #put in if statement
                coarse_time_offset: 1000,
                debug_no_frer: false,
                debug_no_wire_filter: false,
                debug_no_roll: false,
                debug_force_cpu: device == "cpu",
                pad_wire_domain: (iplane > 1), #Non-periodic planes get padded
                use_fft_best_length: true,
                unsqueeze_input: false,
            },
        }, nin=1, nout=1, uses=[torch_frer, the_wire_filter]);
        local gaus_filter = {
            type: 'HfFilter',
            name: 'Gaus_wide',
            data: {
                max_freq: 0.001,  // warning: units
                power: 2,
                flag: true,
                sigma: 0.00012,  // caller should provide
                use_negative_freqs: false,
                device: device,
            }
        };
        local torch_gaus_filter = {
            type: "Torch1DSpectrum",
            name: "torch_1dspec_gaus",
            uses: [gaus_filter],
            data: {
                spectra: [
                    wc.tn(gaus_filter),
                ],
                device: device,
            },
        };
        local gaus = pg.pnode({
            type: 'SPNGApply1DSpectrum',
            name: 'spng_gaus_apa%d_plane%d' % [anode.data.ident, iplane],
            data: {
                base_spectrum_name: wc.tn(torch_gaus_filter),
                dimension: 2,
                output_set_tag: "HfGausWide",
            },
        }, nin=1, nout=1, uses=[torch_gaus_filter]);

        pg.pipeline([decon, gaus]);

    local fanin =  pg.pnode({
        type: 'TorchTensorSetToFrameFanin',
        name: "",
        data: {
            anode: wc.tn(anode),
            input_groups: groups.wpids,
        },
    }, nin=groups.count, nout=1, uses=[anode]);


    local giota = wc.iota(groups.count);
    local pipes = [decon_pipeline(n) for n in giota];

    local body = pg.intern(innodes=[fanout], centernodes=pipes, outnodes=[fanin],
                           edges=[
                               pg.edge(fanout, pipes[n], n, 0) for n in giota
                           ] + [
                               pg.edge(pipes[n], fanin, 0, n) for n in giota
                           ]);
    body;

    

/// Return a frame->[decon]->frame compound node that does TDM decon.
local tdm_decon_node(detector, anode, adc, device="cpu") = 
    local verbosity = 2;          // loud debug logging
    local res = responses(detector, anode, adc.tick);

    local plane_filters = [decon_filters(plane) for plane in wc.iota(3)];

    local groups = plane_groups(detector, anode);

    // fan size.  one extra sends the input frame to the output
    local multiplicity = groups.count + 1;

    local frametotdm = pg.pnode({
        type: 'SPNGFrameToTdm',
        name: "totdm",
        data: {
            anode: wc.tn(anode),
            verbosity: verbosity,
            rules: [{
                tag: "",
                groups: [{
                    wpids: groups.wpids[g],
                }, for g in wc.iota(groups.count)]}],
        },
    }, nin=1, nout=1, uses=[anode]);
    
    // One goes into SP, the other joins the result at the end of the graph to
    // pass through any metadata and non-trace tensors.
    local fanpass = pg.pnode({
        type: 'SPNGFanoutTensorSets',
        name: "fanpass",
        data: {
            multiplicity: 2,
            verbosity: verbosity,
        },
    }, nin=1, nout=multiplicity);

    local unpack = pg.pnode({
        type: 'SPNGTorchSetUnpacker',
        name: "unpack",
        data: {
            selections: [ {datapath:"/frames/\\d+/tags/null/rules/0/groups/%d/traces" % group}
                          for group in wc.iota(groups.count) ],
        },
    }, nin=1, nout=groups.count);


    /// Loop over "groups" making nodes that operate on a single traces tensor.

    /// This applies the gauss filter as part of the decon filter.
    local decon_monolith = [
        local plane_index = groups.planes[group];
        local tag  = "gauss";
        local gname = 'g%d' % group;
        local name = plane_context_name(plane_index, anode) + gname;
        local df = plane_filters[plane_index];
        local filt = filter_kernel(name, [df.channel, df.gauss]);
        local resp = response_kernel(plane_index, res.fr, res.er, adc, name);
        local kern = decon_kernel(filt, resp, name);

        // local axis = convo_options(plane_index).axis.decon_resp;
        local co = convo_options(plane_index);

        // convo_node(kern, axis, tag, name)
        convo_node(kern, [co.decon_channel, co.decon_time], tag, name)

        for group in wc.iota(groups.count)
    ];
    // local pipes = decon_monolith;

    /// This applies gauss filtering as a follow-on.
    local decon_staged = [
        local plane_index = groups.planes[group];
        local gname = 'g%d' % group;
        local name = plane_context_name(plane_index, anode) + gname;

        local df = plane_filters[plane_index];
        // local coa = convo_options(plane_index).axis;
        local co = convo_options(plane_index);

        // The response decon
        local decon_tag = 'decon';
        local decon_name = decon_tag+'-'+name;
        local decon_filt = filter_kernel(decon_name, [df.channel, df.none]);
        local decon_resp = response_kernel(plane_index, res.fr, res.er, adc, decon_name);
        local decon_kern = decon_kernel(decon_filt, decon_resp, decon_name);
        local decon_node = convo_node(decon_kern, [co.decon_channel, co.decon_time], decon_tag, decon_name);

        // In this case we put gauss on its own to make room for future addition
        // of other filters.  Note, we could (should?) keep gauss in the decon
        // and then apply other filters "normalized" by gauss.
        local gauss_tag = 'gauss';
        local gauss_name = gauss_tag+'-'+name;
        local gauss_filt = filter_kernel(gauss_name, [df.none, df.gauss]);
        local gauss_node = convo_node(gauss_filt, [co.nothing, co.just_filter], gauss_tag, gauss_name);

        pg.pipeline([decon_node, gauss_node])

        for group in wc.iota(groups.count)
    ];
    local pipes = decon_staged;

    local repack = pg.pnode({
        type: 'SPNGTorchPacker',
        name: "repack",
        data: {
            multiplicity: groups.count,
        },
    }, nin=groups.count, nout=1);


   local fanin = pg.pipeline([
        pg.pnode({
            type: 'SPNGFaninTensorSets',
            name: "fanin",
            data: {
                multiplicity: 2,
                verbosity: verbosity,
            },
        }, nin=2, nout=1),

        pg.pnode({
            type: 'SPNGTdmToFrame',
            name: "fromtdm",
            data: {
                verbosity: verbosity,
                frame: {datapath: "/frames/\\d+/frame"},
                #frame: {datapath: "/frames/0/frame"},
                // chmasks: ...
                tagged_traces: [ {
                    // eg datapath of /frames/0/tags/gauss/groups/0/traces
                    traces: { tag: "gauss" },
                    // eg datapath of /frames/0/tags/null/rules/0/groups/0/chids
                    chids: { tag: "null" },
                }]
            },
        }, nin=1, nout=1, uses=[]),
    ]);

    // 1->2
    local frontend = pg.pipeline([frametotdm, fanpass]);

    // 1->4->1
    local groupfan = pg.intern(
        innodes=[unpack], centernodes=pipes, outnodes=[repack],
        edges=[
            pg.edge(unpack, pipes[n], n, 0)
            for n in wc.iota(groups.count)
        ] + [
            pg.edge(pipes[n], repack, 0, n)
            for n in wc.iota(groups.count)
        ]);

    // 1->2->{1, 4->1}->1
    local body = pg.intern(
        innodes=[frontend], centernodes=[groupfan], outnodes=[fanin],
        edges=[
            // slightly subtle ordering issue.  The set carrying the "frame"
            // parent tensor must be first.
            pg.edge(frontend, fanin, 0, 0),
            pg.edge(frontend, groupfan, 1, 0),
            pg.edge(groupfan, fanin, 0, 1)
        ]);

    body;



/// Note, different input may require selecting a different anodeid to get any activity.
function(input="test/data/muon-depos.npz", output="", anodeid="3", device="cpu", which="tdm")
    
    // One fine day, make this a CLI argument.  Until then, we try to make it
    // parameterize this config but it's mostly a show.
    local detector="pdhd";


    # FIXME: remove omni and make this file stand-alone.
    local omni = omnimap[detector];
    local aid = std.parseInt(anodeid);
    local anode = omni.anodes[aid];

    local adc = adc_params(detector, anode);

    local sim_source = sim_source_node(input, anode, omni, adc.nticks);

    local source = if output == ""
                   then sim_source
                   else pg.pipeline([sim_source, frame_tap(output%"adc")]);
    
    local tdm_decon = tdm_decon_node(detector, anode, adc, device);
    local tdm_sink = frame_sink(tdm_decon, output, "tdm", "tdm-sig");

    local orig_decon = orig_decon_node(detector, anode, adc, device);
    local orig_sink = frame_sink(orig_decon, output, "orig", "orig-sig");

    local both_sink = frame_sink_two(tdm_sink, orig_sink, "both");

    local sinks = { tdm: tdm_sink, orig: orig_sink, both: both_sink };
    local graph = pg.pipeline([source, sinks[which]]);
    pg.main(graph, 'Pgrapher', plugins=["WireCellSpng"])

