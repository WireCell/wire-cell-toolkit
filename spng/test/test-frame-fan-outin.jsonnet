///  A test job converting an IFrame to tensors and back.

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local omnimap = import "omni/all.jsonnet";
local io = import "fileio.jsonnet";

/// WARNING: is this detector agnostic? I don't think so.

local spng_filters = import "spng_filters.jsonnet";

/// Note, different input may require selecting a different anodeid to get any activity.
function(name="pdhd", input="test/data/muon-depos.npz", output="test-frame-fan-outin.npz", paradigm="orig", anodeid="3", device="cpu")
    local source = io.depo_file_source(input);
    local omni = omnimap[name];
    local aid = std.parseInt(anodeid);
    local anode = omni.anodes[aid];
    local drift = omni.drift(anode);
    local signal = omni.signal(anode);
    local noise = omni.noise(anode);
    local digitize = omni.digitize(anode);
    local verbose = 2;
    local filters = spng_filters(device);

    local groups = [
        [wc.WirePlaneId(wc.Ulayer, 0, anode.data.ident),
         wc.WirePlaneId(wc.Ulayer, 1, anode.data.ident)],
        
        [wc.WirePlaneId(wc.Vlayer, 0, anode.data.ident),
         wc.WirePlaneId(wc.Vlayer, 1, anode.data.ident)],
        
        [wc.WirePlaneId(wc.Wlayer, 0, anode.data.ident)],
        
        [wc.WirePlaneId(wc.Wlayer, 1, anode.data.ident)],
    ];

    local all_fans = {
        local ngroups = std.length(groups),

        // This option runs the initial mainline, non-TDM SPNG nodes
        orig: {
            multiplicity: ngroups,
            fanout : pg.pnode({
                type: 'FrameToTorchSetFanout',
                name: "fanout",
                data: {
                    anode: wc.tn(anode),
                    expected_nticks: 6000,
                    output_groups: groups,
                    unsqueeze_output: true,
                },
            }, nin=1, nout=ngroups, uses=[anode]),
            fanin : pg.pnode({
                type: 'TorchTensorSetToFrameFanin',
                name: "",
                data: {
                    anode: wc.tn(anode),
                    input_groups: groups,
                },
            }, nin=ngroups, nout=1, uses=[anode]),

        },

        // This option runs the TDM-compliant nodes.
        //
        // It fans out ngroups holding just the traces plus one extra shuttling
        // the rest forward.
        tdm: {

            local multiplicity = ngroups + 1, // fan size
            multiplicity: multiplicity,

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
            }, nin=1, nout=1, uses=[anode]),

            local fanout = pg.pnode({
                type: 'SPNGFanoutNode',
                name: "fanout",
                data: {
                    multiplicity: multiplicity,
                    verbose: verbose,
                },
            }, nin=1, nout=multiplicity),


            // A special selection for all-but-traces will go first in the fan.
            local header = pg.pnode({
                type: 'SPNGFunctionNode',
                name: "header",
                data: {
                    verbose: verbose,
                    tensor_selection:  [
                        { reject: "/frames/\\d+/tags/null/rules/0/groups/\\d+/traces" }
                    ],
                    keep_unselected: true,
                    select_parents: false,
                    combine_policy: "produced_only",
                },
            }, nin=1, nout=1, uses=[]),
        
            local res = omni.responses(anode, "sp"),
            local fr = res.fr[0],
            local er = res.er[0],  // unused???
            /// Things taking nchans should instead take an anode and figure this out in C++.
            local bad_bad_bad_nchans = [800, 800, 480, 480],

            /// Loop over groups making nodes that operate on a single traces tensor.
            local traces = [

                local frer = {
                    type: "TorchFRERSpectrum",
                    name: "frer-a%02d-p%d"%[aid,group],
                    data: {
                        field_response: wc.tn(fr),
                        elec_response: wc.tn(er),
                        fr_plane_id: if group > 2 then 2 else group,
                        ADC_mV: 1.0/omni.adc.lsb_voltage,
                        gain: 1.0,
                        default_nchans : bad_bad_bad_nchans[group],
                        default_nticks: 6000,
                        readout_period: 500*wc.ns,
                        extra_scale: 1.0,
                        device: device,
                    },
                    uses: [fr, er],
                };
                local wire_filter = 
                    if group < 2
                    then filters.torch_wire_filters[0]
                    else filters.torch_wire_filters[1];


                pg.pipeline([

                    // select
                    pg.pnode({
                        type: 'SPNGFunctionNode',
                        name: "group%02d"%group,
                        data: {
                            verbose: verbose,
                            tensor_selection:  [
                                { accept: "/frames/\\d+/frame" },
                                { accept: "/frames/\\d+/tags/null/rules/0/groups/%d/traces" % group }
                            ],
                            keep_unselected: false,
                            select_parents: false,
                            combine_policy: "produced_only",
                        },
                    }, nin=1, nout=1, uses=[]), 

                    // decon
                    pg.pnode({
                        type: 'SPNGDecon',
                        name: "decon%02d"%group,
                        data: {
                            tensor_index: 1, // "frame" is 0.
                            frer_spectrum: wc.tn(frer),
                            wire_filter: wc.tn(wire_filter),
                            coarse_time_offset: 1000,
                            pad_wire_domain: (group > 1), #Non-periodic planes get padded
                            use_fft_best_length: true,
                            unsqueeze_input: true,
                            device: device,
                        },
                    }, nin=1, nout=1, uses=[frer, wire_filter]),
                ]) for group in wc.iota(ngroups)],
        
            local pipes = [header] + traces,

            // At the end of this fan are tensor sets with just "frame" and one "traces" tensors.
            fanout: pg.pipeline([frametotdm, 
                                 pg.intern(innodes=[fanout], outnodes=pipes, edges=[
                                     pg.edge(fanout, pipes[n], n, 0) for n in wc.iota(multiplicity)])]),
            fanin: pg.pipeline([
                pg.pnode({
                    type: 'SPNGFaninNode',
                    name: "fanin",
                    data: {
                        multiplicity: multiplicity,
                        verbose: verbose,
                    },
                }, nin=multiplicity, nout=1),

                // big fat FIXME: need something to combine traces/chids that
                // were separated into "groups" (planes).  Make a new component
                // and do not this to TdmToFrame.

                pg.pnode({
                    type: 'SPNGTdmToFrame',
                    name: "fromtdm",
                    data: {
                        verbose: verbose,
                    },
                }, nin=1, nout=1, uses=[]),
            ]),
        },
    };
    local fans = all_fans[paradigm];
    local body = pg.intern(innodes=[fans.fanout],
                           outnodes=[fans.fanin],
                           edges=[pg.edge(fans.fanout, fans.fanin, n, n)
                                  for n in wc.iota(fans.multiplicity)]);
    //local sink = io.frame_file_sink(output);
    local sink = pg.pnode({ type: "DumpFrames" }, nin=1, nout=0);
    local graph = pg.pipeline([source, drift, signal, noise, digitize, body, sink]);
    pg.main(graph, 'Pgrapher', plugins=["WireCellSpng"])
