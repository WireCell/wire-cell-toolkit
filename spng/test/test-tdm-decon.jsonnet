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
        then filter(scale=1.0 / wc.sqrtpi * 0.75, period=1.0)  // ind
        else filter(scale=1.0 / wc.sqrtpi * 10.0, period=1.0), // col

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
        else                        // u
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

    /// fixme: add more to cover ROI needs
};

/// Bind an instance of the above generic filter configuration to
/// a specific implementation of the components.
///
local decon_kernel(filters, fr, er, which="gauss") = {
    /// A name that uniquely identifies the context.
    local name = plane_context_name(filters.plane_index, filters.anode),

    // two types of filter
    local filter = {
        type: "SPNGFilterKernel",
        name: name + which,
        data: {
            axis: [ filters.channel, filters[which] ],
            // axis: [ {kind: "identity"}, {kind: "identity"} ],
            // axis: [ {kind: "identity"}, filters[which] ],
            debug_filename: "filter-kernel-%s.pkl" % (name+which)
        },
        // no uses
    },

    local response = {
        type: "SPNGResponseKernel",
        name: fr.name+er.name + name,
        data: {
            field_response: wc.tn(fr),
            elec_response: wc.tn(er),
            plane_index: filters.plane_index,
            debug_filename: "response-kernel-%s.pkl" % (name+which)
        },
        uses: [fr, er],         // NOT anode.
    },

    type: "SPNGDeconKernel",
    name: which+"-decon-"+name,
    data: {
        filter: wc.tn(filter),
        response: wc.tn(response),
        debug_filename: "decon-kernel-%s.pkl" % (name+which)
    },
    uses: [response, filter]
};

/// Return config object for a KernelConvovle.  Use extra name to make otherwise
/// identical instances (eg for one plane split into groups)
local convo_node(kernel, plane_index, extra_name="", which="") =
    // For channel dimmension, wrapped wire planes are cyclic.
    local channel_options = [
        {cyclic: true,  crop:  0}, // u, cycle, wrapped
        {cyclic: true,  crop:  0}, // v, cyclic, wrapped
        {cyclic: false, crop: -2}  // w, linear
    ][plane_index];
    /// Would NOT crop in chunked-streaming mode.
    // local time_options = {cyclic: false, crop: -2, roll: 1000, roll_mode: "decon"};
    local time_options = {cyclic: false, crop: 0};
    /// Note, may need to change KC to roll then crop.

    // one decon per filter type
    pg.pnode({
        type: "SPNGKernelConvolve",
        name: kernel.name + extra_name,
        data: {
            kernel: wc.tn(kernel),
            axis: [
                channel_options,
                time_options,
            ],
            tag: which,
            datapath_format: "/frames/{ident}/tags/{tag}/groups/{group}/traces",
            //faster: true,       // use "faster DFT size"
        },
    }, nin=1, nout=1, uses=[kernel]);



/// Note, different input may require selecting a different anodeid to get any activity.
function(input="test/data/muon-depos.npz", output="", anodeid="3", device="cpu")

    local adc_tick = 500*wc.ns;
    local adc_nticks = 6000;

    local depo_source = io.depo_file_source(input);

    # FIXME: remove omni and make this file stand-alone.
    local detector="pdhd";
    local omni = omnimap[detector];
    local aid = std.parseInt(anodeid);

    local anode = omni.anodes[aid];
    local drift = omni.drift(anode);
    local signal = omni.signal(anode);
    local noise = omni.noise(anode);
    local digitize = omni.digitize(anode);

    // This is what sets the ADC waveform readout window.
    local reframer = pg.pnode({
        type: 'Reframer',
        name: anodeid,
        data: {
            anode: wc.tn(anode),
            nticks: adc_nticks,
        },
    }, nin=1, nout=1, uses=[anode]);

    local verbose = 2;

    local fr = {
        type: 'FieldResponse',
        name: detector,
        data: {
            filename: "dune-garfield-1d565.json.bz2",
        },
    };
    local er = {
        type: 'ColdElecResponse',
        name: detector,
        data: {
            tick: adc_tick,
            // enough to cover the ER, DO NOT PUT FULL READOUT SIZE HERE
            nticks: 40,
            shaping: 2.2*wc.us,
            gain : 14.0*wc.mV/wc.fC,
            postgain: 1.0,
        },
    };
    local rcr = {
        type: 'RCResponse',
        name: detector,
        data: {
            width: 1.1*wc.ms,
        }
    };
    local rc = [rcr, rcr];

    local plane_filters = [decon_filters(plane) for plane in wc.iota(3)];

    local groups = [
        [wc.WirePlaneId(wc.Ulayer, 0, anode.data.ident),
         wc.WirePlaneId(wc.Ulayer, 1, anode.data.ident)],
        
        [wc.WirePlaneId(wc.Vlayer, 0, anode.data.ident),
         wc.WirePlaneId(wc.Vlayer, 1, anode.data.ident)],
        
        [wc.WirePlaneId(wc.Wlayer, 0, anode.data.ident)],
        
        [wc.WirePlaneId(wc.Wlayer, 1, anode.data.ident)],
    ];
    local group_plane = [0, 1, 2, 2];

    local ngroups = std.length(groups);
    // fan size.  one extra sends the input frame to the output
    local multiplicity = ngroups + 1;

    local frametotdm = pg.pnode({
        type: 'SPNGFrameToTdm',
        name: "totdm",
        data: {
            anode: wc.tn(anode),
            verbose: verbose,
            rules: [{
                tag: "",
                groups: [{
                    wpids: groups[g],
                }, for g in wc.iota(ngroups)]}],
        },
    }, nin=1, nout=1, uses=[anode]);
    
    // One goes into SP, the other joins the result at the end of the graph to
    // pass through any metadata and non-trace tensors.
    local fanpass = pg.pnode({
        type: 'SPNGFanoutNode',
        name: "fanpass",
        data: {
            multiplicity: 2,
            verbose: verbose,
        },
    }, nin=1, nout=multiplicity);

    local unpack = pg.pnode({
        type: 'SPNGTorchSetUnpacker',
        name: "unpack",
        data: {
            selections: [ {datapath:"/frames/\\d+/tags/null/rules/0/groups/%d/traces" % group}
                          for group in wc.iota(ngroups) ],
        },
    }, nin=1, nout=ngroups);


    /// Loop over "groups" making nodes that operate on a single traces tensor.
    local pipes = [
        local plane = group_plane[group];
        local which  = "gauss";
        // to start, just decon, later turn this into a deeper pipeline
        convo_node(decon_kernel(decon_filters(plane), fr, er, which), plane, '-group%d'%group, which)

        for group in wc.iota(ngroups)
    ];
    
    local repack = pg.pnode({
        type: 'SPNGTorchPacker',
        name: "repack",
        data: {
            multiplicity: ngroups,
        },
    }, nin=ngroups, nout=1);

    local fanin = pg.pipeline([
        pg.pnode({
            type: 'SPNGFaninNode',
            name: "fanin",
            data: {
                multiplicity: 2,
                verbose: verbose,
            },
        }, nin=2, nout=1),

        pg.pnode({
            type: 'SPNGTdmToFrame',
            name: "fromtdm",
            data: {
                verbose: verbose,
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
            for n in wc.iota(ngroups)
        ] + [
            pg.edge(pipes[n], repack, 0, n)
            for n in wc.iota(ngroups)
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

    local sim_source = pg.pipeline([depo_source, drift, signal, reframer, noise, digitize]);
    local frame_fout = pg.pnode({
        type: 'FrameFanout',
        name: '',
        data: {
            multiplicity: 2
        }
    }, nin=1, nout=2);
    local adc_sink = io.frame_file_sink(output % "adc", digitize=true);
    local adc_tap = pg.intern(innodes=[frame_fout],
                              centernodes=[adc_sink],
                              outnodes=[frame_fout],
                              edges=[pg.edge(frame_fout, adc_sink, 1, 0)]);

    local source = if output == ""
                   then sim_source
                   else pg.pipeline([sim_source, adc_tap]);
    
    local sink = if output == ""
                 then pg.pnode({ type: "DumpFrames", name: "" }, nin=1, nout=0)
                 else io.frame_file_sink(output % "sig", digitize=false);
    // else io.frame_tensor_file_sink(output, digitize=false, mode="dense");

    local graph = pg.pipeline([source, body, sink]);
    pg.main(graph, 'Pgrapher', plugins=["WireCellSpng"])

