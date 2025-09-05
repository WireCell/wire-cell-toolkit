// TODO -- brief descrip

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';
local spng_filters = import 'spng_filters.jsonnet';

function(tools, debug_force_cpu=false, apply_gaus=true, do_roi_filters=false, do_collate_apa=false, do_mp_finding=false) {
    // make_spng :: function(tools, debug_force_cpu=false, apply_gaus=true, do_roi_filters=false, do_collate_apa=false) {

        local filter_settings = {
            debug_force_cpu: debug_force_cpu,
        },
        local filters = spng_filters(filter_settings),

        local make_fanout(anode, name=null) = {
            ret : g.pnode({
                type: 'FrameToTorchSetFanout',
                name:
                    if std.type(name) == 'null'
                    then anode.name + 'torchfanout%d' % anode.data.ident
                    else name,
                data: {
                    anode: wc.tn(anode),
                    expected_nticks: 6000,

                    output_groups: [
                        [wc.WirePlaneId(wc.Ulayer, 0, anode.data.ident),
                        wc.WirePlaneId(wc.Ulayer, 1, anode.data.ident)],
                        
                        [wc.WirePlaneId(wc.Vlayer, 0, anode.data.ident),
                        wc.WirePlaneId(wc.Vlayer, 1, anode.data.ident)],
                        
                        [wc.WirePlaneId(wc.Wlayer, 0, anode.data.ident)],
                        
                        [wc.WirePlaneId(wc.Wlayer, 1, anode.data.ident)],
                    ],
                    debug_force_cpu: debug_force_cpu,
                    unsqueeze_output: true,


                }
            }, nin=1, nout=4, uses=[anode]),
        }.ret,

        # TODO -- Abstract away
        local nchans = [800, 800, 480, 480],

        local make_replicator_post_tight(anode, iplane) = {
            ret : g.pnode({
                type: 'TorchTensorSetReplicator',
                name: 'post_tight_replicator_%d_%d' % [anode.data.ident, iplane],
                data: {
                    multiplicity: 5,
                }
            }, nin=1, nout=5, uses=[anode])
        }.ret,

        local make_replicator_post_gaus(anode, iplane) = {
            ret : g.pnode({
                type: 'TorchTensorSetReplicator',
                name: 'post_gaus_replicator_%d_%d' % [anode.data.ident, iplane],
                data: {multiplicity: 2,}

            }, nin=1, nout=2, uses=[anode])
        }.ret,
        
        local make_replicator_post_gaus_simple(anode, iplane) = {
            ret : g.pnode({
                type: 'TorchTensorSetReplicator',
                name: 'post_gaus_replicator_%d_%d' % [anode.data.ident, iplane],
                data: {multiplicity: (if iplane > 1 then 2 else 6)}
            }, nin=1, nout=(if iplane > 1 then 2 else 6), uses=[anode])
        }.ret,

        local make_pipeline(anode, iplane, apply_gaus=true, do_roi_filters=false) = {
            local the_field = if std.length(tools.fields) > 1 then tools.fields[anode.data.ident] else tools.fields[0],
            local torch_frer = {
                type: "TorchFRERSpectrum",
                name: "torch_frer%d_plane%d" % [anode.data.ident, iplane],
                uses: [
                    the_field,
                    tools.elec_resp
                ],
                data: {
                    field_response: wc.tn(the_field),#"FieldResponse:field%d"% anode.data.ident,
                    fr_plane_id: if iplane > 2 then 2 else iplane,
                    ADC_mV: 11702142857.142859,
                    gain: 1.0,
                    default_nchans : nchans[iplane],
                    default_nticks: 6000,
                    readout_period: 500.0, #512.0,
                    extra_scale: 1.0,
                    debug_force_cpu: debug_force_cpu,

                }
            },

            local the_wire_filter = if iplane > 1 then filters.torch_wire_filters[1] else filters.torch_wire_filters[0],

            local spng_decon = g.pnode({
                type: 'SPNGDecon',
                name: 'spng_decon_apa%d_plane%d' % [anode.data.ident, iplane],
                data: {
                    frer_spectrum: wc.tn(torch_frer),
                    wire_filter: wc.tn(the_wire_filter), #put in if statement
                    coarse_time_offset: 1000,
                    debug_no_frer: false,
                    debug_no_wire_filter: false,
                    debug_no_roll: false,
                    debug_force_cpu: debug_force_cpu,
                    pad_wire_domain: (iplane > 1), #Non-periodic planes get padded
                    use_fft_best_length: true,
                    unsqueeze_input: false,
                },
            },
            nin=1, nout=1,
            uses=[torch_frer, the_wire_filter]),

            local threshold_rois = g.pnode({
                type: 'SPNGThresholdROIs',
                name: 'spng_threshold_rois_apa%d_plane%d' % [anode.data.ident, iplane],
                data: {
                    unsqueeze_input: false,
                    threshold_rms_factor: 3.0,
                    debug_force_cpu: debug_force_cpu,
                }
            },
            nin=1, nout=1, uses=[]),

            local the_wiener_tight = if iplane < 3 then filters.wiener_tight_filters[iplane] else filters.wiener_tight_filters[2],
            local the_torch_wiener_tight = if iplane < 3 then filters.torch_wiener_tight_filters[iplane] else filters.torch_wiener_tight_filters[2],
            local the_torch_wiener_wide = if iplane < 3 then filters.torch_wiener_wide_filters[iplane] else filters.torch_wiener_wide_filters[2],

            local torch_roi_tight = {
                    type: "Torch1DSpectrum",
                    name: "torch_1dspec_roi_tight",
                    uses: [the_wiener_tight, filters.ROI_tight_lf],
                    data: {
                        spectra: [
                            wc.tn(the_wiener_tight),
                            wc.tn(filters.ROI_tight_lf),
                        ],
                        debug_force_cpu: debug_force_cpu,
                    },
            },
            local torch_roi_tighter = {
                    type: "Torch1DSpectrum",
                    name: "torch_1dspec_roi_tighter",
                    uses: [the_wiener_tight, filters.ROI_tighter_lf],
                    data: {
                        spectra: [
                            wc.tn(the_wiener_tight),
                            wc.tn(filters.ROI_tighter_lf),
                        ],
                        debug_force_cpu: debug_force_cpu,
                    },
            },

            
            local torch_roi_loose = {
                    type: "Torch1DSpectrum",
                    name: "torch_1dspec_roi_loose",
                    uses: [the_wiener_tight, filters.ROI_loose_lf, filters.ROI_tight_lf],
                    data: {
                        spectra: [
                            wc.tn(the_wiener_tight),
                            wc.tn(filters.ROI_loose_lf),
                            wc.tn(filters.ROI_tight_lf),
                        ],
                        debug_force_cpu: debug_force_cpu,
                    },
            },


            local apply_wiener_tight = g.pnode({
                type: 'SPNGApply1DSpectrum',
                name: 'spng_wiener_tight_apa%d_plane%d' % [anode.data.ident, iplane],
                data: {
                    base_spectrum_name: wc.tn(the_torch_wiener_tight),
                    dimension: 2,
                    // target_tensor: "HfGausWide",
                    output_set_tag: "WienerTight",
                },
            },
            nin=1, nout=1,
            uses=[the_torch_wiener_tight]),
            
            local apply_wiener_wide = g.pnode({
                type: 'SPNGApply1DSpectrum',
                name: 'spng_wiener_wide_apa%d_plane%d' % [anode.data.ident, iplane],
                data: {
                    base_spectrum_name: wc.tn(the_torch_wiener_wide),
                    dimension: 2,
                    // target_tensor: "HfGausWide",
                    output_set_tag: "WienerWide",
                },
            },
            nin=1, nout=1,
            uses=[the_torch_wiener_wide]),
            
            local apply_loose_roi = g.pnode({
                type: 'SPNGApply1DSpectrum',
                name: 'spng_loose_roi_apa%d_plane%d' % [anode.data.ident, iplane],
                data: {
                    base_spectrum_name: wc.tn(torch_roi_loose),
                    dimension: 2,
                    // target_tensor: 'HfGausWide',
                    output_set_tag: 'ROILoose',
                },
            },
            nin=1, nout=1,
            uses=[torch_roi_loose]),
            
            local apply_tight_roi = g.pnode({
                type: 'SPNGApply1DSpectrum',
                name: 'spng_tight_roi_apa%d_plane%d' % [anode.data.ident, iplane],
                data: {
                    base_spectrum_name: wc.tn(torch_roi_tight),
                    dimension: 2,
                    // target_tensor: 'HfGausWide',
                    output_set_tag: 'ROITight',
                },
            },
            nin=1, nout=1,
            uses=[torch_roi_tight]),

            local apply_tighter_roi = g.pnode({
                type: 'SPNGApply1DSpectrum',
                name: 'spng_tighter_roi_apa%d_plane%d' % [anode.data.ident, iplane],
                data: {
                    base_spectrum_name: wc.tn(torch_roi_tighter),
                    dimension: 2,
                    // target_tensor: 'HfGausWide',
                    output_set_tag: 'ROITighter',
                },
            },
            nin=1, nout=1,
            uses=[torch_roi_tighter]),

            local spng_gaus_app = g.pnode({
                type: 'SPNGApply1DSpectrum',
                name: 'spng_gaus_apa%d_plane%d' % [anode.data.ident, iplane],
                data: {
                    base_spectrum_name: wc.tn(filters.torch_gaus_filter), #put in if statement
                    dimension: 2,
                    // target_tensor: "Default",
                    output_set_tag: "HfGausWide",
                },
            },
            nin=1, nout=1,
            uses=[filters.torch_gaus_filter]),

            local torch_to_tensor = g.pnode({
                type: 'TorchToTensor',
                name: 'torchtotensor_%d_%d' % [anode.data.ident, iplane],
                data: {},
            }, nin=1, nout=1),

            local tagger = g.pnode({
                type: 'SPNGTorchTensorSetTagger',
                name: 'test_tagger_%d_%d' % [anode.data.ident, iplane],
                data: {
                    allow_retagging: false,
                    tag_list: {
                        a: '1',
                        b: '2',
                    },
                }
            }, nin=1, nout=1),

            local tensor_sink = g.pnode({
                type: 'TensorFileSink',
                name: 'tfsink_%d_%d' % [anode.data.ident, iplane],
                data: {
                    outname: 'testout_fan_apa%d_plane%d.tar' % [anode.data.ident, iplane],
                    prefix: ''
                },
            }, nin=1, nout=0),


            local decon_and_gaus = [spng_decon] + (if apply_gaus then [spng_gaus_app] else []),
            local convert_and_sink = [torch_to_tensor, tensor_sink],

            local post_gaus_replicator = make_replicator_post_gaus(anode, iplane),
            local post_tight_replicator = if iplane < 2 then make_replicator_post_tight(anode, iplane) else null,

            local post_gaus_replicator_simple = make_replicator_post_gaus_simple(anode, iplane),

            local collator = g.pnode({
                type: 'TorchTensorSetCollator',
                name: 'collate_%d_%d' % [anode.data.ident, iplane],
                data: {
                    output_set_tag: 'collated_%d' % iplane,
                    // multiplicity: (if iplane < 2 then 5 else 2),
                    multiplicity: (if iplane < 2 then 6 else 2),
                },
            // }, nin=if iplane < 2 then 5 else 2, nout=1),
            }, nin=if iplane < 2 then 6 else 2, nout=1),

            local post_gaus_filters = if iplane < 2 then g.intern(
                innodes=[post_gaus_replicator_simple],
                outnodes=[collator],
                centernodes=[apply_wiener_tight, apply_wiener_wide, apply_tight_roi, apply_tighter_roi, apply_loose_roi],
                edges = [
                    g.edge(post_gaus_replicator_simple, apply_wiener_tight, 0, 0),
                    g.edge(apply_wiener_tight, collator, 0, 0),

                    g.edge(post_gaus_replicator_simple, apply_wiener_wide, 1, 0),
                    g.edge(apply_wiener_wide, collator, 0, 1),

                    g.edge(post_gaus_replicator_simple, apply_tight_roi, 2, 0),
                    g.edge(apply_tight_roi, collator, 0, 2),
                    
                    g.edge(post_gaus_replicator_simple, apply_tighter_roi, 3, 0),
                    g.edge(apply_tighter_roi, collator, 0, 3),

                    g.edge(post_gaus_replicator_simple, apply_loose_roi, 4, 0),
                    g.edge(apply_loose_roi, collator, 0, 4),

                    g.edge(post_gaus_replicator_simple, collator, 5, 5),
                ],
            ) else g.intern(
                innodes=[post_gaus_replicator_simple],
                outnodes=[collator],
                centernodes=[apply_wiener_tight, apply_wiener_wide],
                edges = [
                    g.edge(post_gaus_replicator_simple, apply_wiener_tight, 0, 0),
                    g.edge(apply_wiener_tight, collator, 0, 0),

                    g.edge(post_gaus_replicator_simple, apply_wiener_wide, 1, 0),
                    g.edge(apply_wiener_wide, collator, 0, 1),

                    // g.edge(post_gaus_replicator_simple, collator, 2, 2)

                ],
            ),


            local full_pipeline = (
                decon_and_gaus + 
                (if do_roi_filters then [post_gaus_filters] else [])
                
            ),

            ret : g.pipeline(
                full_pipeline
            ),
            

        }.ret,

        // local tf_fans = [make_fanout(a) for a in tools.anodes],
        
        local spng_fanout(anode, apply_gaus=true, do_roi_filters=false, do_collate_apa=false) = {

            #FrameToTorchSetFanout
            local tf_fan = make_fanout(anode),

            #SPNGDecon + Output stuff -- per plane
            local pipelines = [
                make_pipeline(anode, iplane, apply_gaus, do_roi_filters)
                for iplane in std.range(0, 3)
            ],
            local torch_to_tensors_planes = [g.pnode({
                type: 'TorchToTensor',
                name: 'torchtotensor_%d_%d' % [anode.data.ident, iplane],
                data: {},
            }, nin=1, nout=1) for iplane in std.range(0,3)],

            local tensor_sinks_planes = [g.pnode({
                type: 'TensorFileSink',
                name: 'tfsink_%d_%d' % [anode.data.ident, iplane],
                data: {
                    outname: 'testout_fan_apa%d_plane%d.tar' % [anode.data.ident, iplane],
                    prefix: ''
                },
            }, nin=1, nout=0) for iplane in std.range(0,3)],

            local full_pipelines = [
                g.pipeline(
                    [pipelines[iplane], torch_to_tensors_planes[iplane], tensor_sinks_planes[iplane]]
                )
            for iplane in std.range(0,3)],

            local apa_collator = g.pnode({
                type: 'TorchTensorSetCollator',
                name: 'apa_collate_%d' % [anode.data.ident],
                data: {
                    output_set_tag: 'apa_collated',
                    multiplicity:4,
                },
            }, nin=4, nout=1),

            local torch_to_tensors_apa = g.pnode({
                type: 'TorchToTensor',
                name: 'torchtotensor_%d' % [anode.data.ident],
                data: {},
            }, nin=1, nout=1),
            local tensor_sinks_apa = g.pnode({
                type: 'TensorFileSink',
                name: 'tfsink_%d' % [anode.data.ident],
                data: {
                    outname: 'testout_fan_apa%d.tar' % [anode.data.ident],
                    prefix: ''
                },
            }, nin=1, nout=0),

            local apa_pipeline = g.pipeline(
                [apa_collator, torch_to_tensors_apa, tensor_sinks_apa]
            ),

            ret : (if !do_collate_apa then
                g.intern(
                    innodes=[tf_fan],
                    // centernodes=[]
                    outnodes=full_pipelines,
                    edges=[
                        g.edge(tf_fan, full_pipelines[0], 0),
                        g.edge(tf_fan, full_pipelines[1], 1),
                        g.edge(tf_fan, full_pipelines[2], 2),
                        g.edge(tf_fan, full_pipelines[3], 3),
                    ]
                ) else
                g.intern(
                    innodes=[tf_fan],
                    centernodes=pipelines,
                    outnodes=[apa_pipeline],
                    edges=[
                        g.edge(tf_fan, pipelines[0], 0),
                        g.edge(tf_fan, pipelines[1], 1),
                        g.edge(tf_fan, pipelines[2], 2),
                        g.edge(tf_fan, pipelines[3], 3),
                        g.edge(pipelines[0], apa_pipeline, 0, 0),
                        g.edge(pipelines[1], apa_pipeline, 0, 1),
                        g.edge(pipelines[2], apa_pipeline, 0, 2),
                        g.edge(pipelines[3], apa_pipeline, 0, 3),
                    ]
                )
            ),
        }.ret,
        
        spng_decons : [
            spng_fanout(anode, apply_gaus=apply_gaus, do_roi_filters=do_roi_filters, do_collate_apa=do_collate_apa) for anode in tools.anodes
        ],
   
        stacked_spng : {
            local tf_fans = [make_fanout(anode) for anode in tools.anodes],

            local u_stacker =  g.pnode({
                type: 'TorchTensorSetStacker',
                name: 'u_stacker',
                data: {
                    output_set_tag: 'u_stacked',
                    multiplicity:4,
                },
            }, nin=4, nout=1),
            local v_stacker =  g.pnode({
                type: 'TorchTensorSetStacker',
                name: 'v_stacker',
                data: {
                    output_set_tag: 'v_stacked',
                    multiplicity:4,
                },
            }, nin=4, nout=1),
            local w_stacker =  g.pnode({
                type: 'TorchTensorSetStacker',
                name: 'w_stacker',
                data: {
                    output_set_tag: 'w_stacked',
                    multiplicity:8,
                },
            }, nin=8, nout=1),


            local u_unstacker =  g.pnode({
                type: 'TorchTensorSetUnstacker',
                name: 'u_unstacker',
                data: {
                    output_set_tag: 'u_stacked',
                    multiplicity:4,
                },
            }, nin=1, nout=4),
            local v_unstacker =  g.pnode({
                type: 'TorchTensorSetUnstacker',
                name: 'v_unstacker',
                data: {
                    output_set_tag: 'v_stacked',
                    multiplicity:4,
                },
            }, nin=1, nout=4),
            
            local u_unstacker_face1 =  g.pnode({
                type: 'TorchTensorSetUnstacker',
                name: 'u_unstacker_face1',
                data: {
                    output_set_tag: 'u_stacked_face1',
                    multiplicity:4,
                },
            }, nin=1, nout=4),
            local v_unstacker_face1 =  g.pnode({
                type: 'TorchTensorSetUnstacker',
                name: 'v_unstacker_face1',
                data: {
                    output_set_tag: 'v_stacked_face1',
                    multiplicity:4,
                },
            }, nin=1, nout=4),

            local w_unstacker =  g.pnode({
                type: 'TorchTensorSetUnstacker',
                name: 'w_unstacker',
                data: {
                    output_set_tag: 'w_stacked',
                    multiplicity:8,
                },
            }, nin=1, nout=8),


            local torch_to_tensors = [
                g.pnode({
                    type: 'TorchToTensor',
                    name: 'torchtotensor_%s' % [i],
                    data: {},
                }, nin=1, nout=1) for i in ['u', 'v', 'w']
            ],

            local torch_to_tensors_unstacked = [
                ([
                    g.pnode({
                        type: 'TorchToTensor',
                        name: 'torchtotensor_%s_apa%d' % [i, anode.data.ident],
                        data: {},
                    }, nin=1, nout=1) for i in ['u', 'v']
                ] + [
                    g.pnode({
                        type: 'TorchToTensor',
                        name: 'torchtotensor_w%d_apa%d' % [i, anode.data.ident],
                        data: {},
                    }, nin=1, nout=1) for i in std.range(0,1)
                ]) for anode in tools.anodes
            ],
            local torch_to_tensors_mp_finding = [
                g.pnode({
                    type: 'TorchToTensor',
                    name: 'torchtotensor_mp_finding_%d' % [anode.data.ident],
                    data: {},
                }, nin=1, nout=1) for anode in tools.anodes
            ],

            local tensor_sinks_unstacked = [
                ([
                    g.pnode({
                        type: 'TensorFileSink',
                        name: 'tfsink_%s_apa%d' % [i, anode.data.ident],
                        data: {
                            outname: 'testout_fan_plane_%s_apa%d.tar' % [i, anode.data.ident],
                            prefix: ''
                        },
                    }, nin=1, nout=0) for i in ['u', 'v']
                ] + [
                    g.pnode({
                        type: 'TensorFileSink',
                        name: 'tfsink_w%d_apa%d' % [i, anode.data.ident],
                        data: {
                            outname: 'testout_fan_plane_w%d_apa%d.tar' % [i, anode.data.ident],
                            prefix: ''
                        },
                    }, nin=1, nout=0) for i in [0,1]
                ]) for anode in tools.anodes
            ],

            local tensor_sinks = [g.pnode({
                type: 'TensorFileSink',
                name: 'tfsink_%s' % i,
                data: {
                    outname: 'testout_fan_plane_%s.tar' % i,
                    prefix: ''
                },
            }, nin=1, nout=0) for i in ['u', 'v', 'w']],

            local torch_roi_loose_spectra = [{
                    type: "Torch1DSpectrum",
                    name: "torch_1dspec_roi_loose_%d" % iplane,
                    uses: [filters.wiener_tight_filters[iplane], filters.ROI_loose_lf],
                    data: {
                        spectra: [
                            wc.tn(filters.wiener_tight_filters[iplane]),
                            wc.tn(filters.ROI_loose_lf),
                            // wc.tn(filters.ROI_tight_lf),
                        ],
                        debug_force_cpu: debug_force_cpu,
                    },
            } for iplane in std.range(0,2)],

            local apply_loose_rois = [g.pnode({
                type: 'SPNGApply1DSpectrum',
                name: 'spng_loose_roi_plane%d' % iplane,
                data: {
                    base_spectrum_name: wc.tn(torch_roi_loose_spectra[iplane]),
                    dimension: 2,
                    // target_tensor: 'HfGausWide',
                    output_set_tag: 'ROILoose',
                },
            },
            nin=1, nout=1,
            uses=[torch_roi_loose_spectra[iplane]]) for iplane in std.range(0,2)],


            local torch_roi_tight_spectra = [{
                    type: "Torch1DSpectrum",
                    name: "torch_1dspec_roi_tight",
                    uses: [filters.wiener_tight_filters[iplane], filters.ROI_tight_lf],
                    data: {
                        spectra: [
                            wc.tn(filters.wiener_tight_filters[iplane]),
                            wc.tn(filters.ROI_tight_lf),
                        ],
                        debug_force_cpu: debug_force_cpu,
                    },
            } for iplane in std.range(0,2)],

            local apply_tight_rois = [g.pnode({
                type: 'SPNGApply1DSpectrum',
                name: 'spng_tight_roi_plane%d' % iplane,
                data: {
                    base_spectrum_name: wc.tn(torch_roi_tight_spectra[iplane]),
                    dimension: 2,
                    // target_tensor: 'HfGausWide',
                    output_set_tag: 'ROITight',
                },
            },
            nin=1, nout=1,
            uses=[torch_roi_tight_spectra[iplane]]) for iplane in std.range(0,2)],

            local fans_to_stacks = g.intern(
                innodes=tf_fans,
                centernodes=[],
                outnodes=[u_stacker, v_stacker, w_stacker],
                edges=[

                    g.edge(tf_fans[0], u_stacker, 0, 0),
                    g.edge(tf_fans[0], v_stacker, 1, 0),
                    g.edge(tf_fans[0], w_stacker, 2, 0),
                    g.edge(tf_fans[0], w_stacker, 3, 1),

                    g.edge(tf_fans[1], u_stacker, 0, 1),
                    g.edge(tf_fans[1], v_stacker, 1, 1),
                    g.edge(tf_fans[1], w_stacker, 2, 2),
                    g.edge(tf_fans[1], w_stacker, 3, 3),

                    g.edge(tf_fans[2], u_stacker, 0, 2),
                    g.edge(tf_fans[2], v_stacker, 1, 2),
                    g.edge(tf_fans[2], w_stacker, 2, 4),
                    g.edge(tf_fans[2], w_stacker, 3, 5),

                    g.edge(tf_fans[3], u_stacker, 0, 3),
                    g.edge(tf_fans[3], v_stacker, 1, 3),
                    g.edge(tf_fans[3], w_stacker, 2, 6),
                    g.edge(tf_fans[3], w_stacker, 3, 7),

                ]
            ),

            local the_field = tools.fields[0],

            local torch_frers = [{
                type: "TorchFRERSpectrum",
                name: "torch_frer_plane%d" % iplane,
                uses: [
                    the_field,
                    tools.elec_resp
                ],
                data: {
                    field_response: wc.tn(the_field),#"FieldResponse:field%d"% anode.data.ident,
                    fr_plane_id: if iplane > 2 then 2 else iplane,
                    ADC_mV: 11702142857.142859,
                    gain: 1.0,
                    default_nchans : nchans[iplane],
                    default_nticks: 6000,
                    readout_period: 500.0,
                    extra_scale: 1.0,
                    debug_force_cpu: debug_force_cpu,

                }
            } for iplane in std.range(0, 2)],
            // local the_wire_filters = [if iplane > 1 then filters.torch_wire_filters[1] else filters.torch_wire_filters[0],
        
            local spng_decons = [
                g.pnode({
                    type: 'SPNGDecon',
                    name: 'spng_decon__plane%d' % iplane,
                    data: {
                        frer_spectrum: wc.tn(torch_frers[iplane]),
                        wire_filter: wc.tn(if iplane < 2 then filters.torch_wire_filters[0] else filters.torch_wire_filters[1]),
                        coarse_time_offset: 1000,
                        debug_no_frer: false,
                        debug_no_wire_filter: false,
                        debug_no_roll: false,
                        debug_force_cpu: debug_force_cpu,
                        pad_wire_domain: (iplane > 1), #Non-periodic planes get padded
                        use_fft_best_length: true,
                        unsqueeze_input: false,
                    },
                },
                nin=1, nout=1,
                uses=[torch_frers[iplane], if iplane < 2 then filters.torch_wire_filters[0] else filters.torch_wire_filters[1]])
                for iplane in std.range(0, 2)
            ],

            local post_decon_replicators = [
                g.pnode({
                    type: 'TorchTensorSetReplicator',
                    name: 'post_decon_replicator_%d' % iplane,
                    data: {multiplicity: 2,}

                }, nin=1, nout=2)
                for iplane in std.range(0,2)
            ],

            local spng_gaus_apps = [g.pnode({
                type: 'SPNGApply1DSpectrum',
                name: 'spng_gaus_plane%d' % iplane,
                data: {
                    base_spectrum_name: wc.tn(filters.torch_gaus_filter), #put in if statement
                    dimension: 2,
                    // target_tensor: "Default",
                    output_set_tag: "HfGausWide",
                },
            },
            nin=1, nout=1,
            uses=[filters.torch_gaus_filter]) for iplane in std.range(0,2)],

            local threshold_rois = [g.pnode({
                type: 'SPNGThresholdROIs',
                name: 'spng_threshold_rois_plane%d' % iplane,
                data: {
                    unsqueeze_input: false,
                    // threshold_rms_factor: 4.5,
                    threshold_rms_factor: 3.0,
                    debug_force_cpu: debug_force_cpu,
                }
            },
            nin=1, nout=1, uses=[]) for iplane in std.range(0,2)],


            local mp_finding_replicators = [
                g.pnode({
                    type: 'TorchTensorSetReplicator',
                    name: '%s' % iplane,
                    data: {multiplicity: 2,}

                }, nin=1, nout=2)
                for iplane in ['u', 'v']
            ],


            local collators_for_mp_finding = [
                g.pnode({
                    type: 'TorchTensorSetCollator',
                    name: 'collate_mp_finding_%d' % anode.data.ident,
                    data: {
                        output_set_tag: 'collated_mp_finding_%d' %  anode.data.ident,
                        multiplicity: 3,
                    },
                }, nin=3, nout=1) for anode in tools.anodes
            ],

            local collators_for_mp_finding_face1 = [
                g.pnode({
                    type: 'TorchTensorSetCollator',
                    name: 'collate_mp_finding_%d_face1' % anode.data.ident,
                    data: {
                        output_set_tag: 'collated_mp_finding_%d_face1' %  anode.data.ident,
                        multiplicity: 3,
                    },
                }, nin=3, nout=1) for anode in tools.anodes
            ],

            // local collators_for_dnnroi_u = [
            //     g.pnode({
            //         type: 'TorchTensorSetCollator',
            //         name: 'collate_dnn_roi_u',
            //         data: {
            //             output_set_tag: 'collated_dnn_roi_u_%d' % anode.data.ident,
            //             multiplicity: 2,
            //         },
            //     }, nin=2, nout=1) // for anode in tools.anodes
            // ],

            local mp_finding = [
                g.pnode({
                    type: 'SPNGNoTileMPCoincidence',
                    name: 'mp_finding_%d' % anode.data.ident,
                    data: {
                        anode: wc.tn(anode),
                        face_index: 0,
                        target_plane_index: 0,
                        aux_plane_l_index: 2,
                        aux_plane_m_index: 1,
                        output_torch_name: "mp_finding_%d_tensors.pt" % anode.data.ident,
                    },
                }, nin=1, nout=1, uses=[anode]) for anode in tools.anodes
            ],

            local mp_finding_face1 = [
                g.pnode({
                    type: 'SPNGNoTileMPCoincidence',
                    name: 'mp_finding_face1_%d' % anode.data.ident,
                    data: {
                        anode: wc.tn(anode),
                        face_index: 1,
                        target_plane_index: 0,
                        aux_plane_l_index: 2,
                        aux_plane_m_index: 1,
                        output_torch_name: "mp_finding_face1_%d_tensors.pt" % anode.data.ident,
                    },
                }, nin=1, nout=1, uses=[anode]) for anode in tools.anodes
            ],

            local roi_application = [
                g.pnode({
                    type: 'SPNGApplyROI',
                    name: 'roi_app_%s' % plane,
                    data: {
                        ROI_tensor_index: 0,
                        value_tensor_index: 1,
                        output_set_tag: 'roi_applied_%s' %plane,
                    },
                }, nin=2, nout=1) for plane in ['u', 'v', 'w']
            ],

            local mp_finding_ors = [
                g.pnode({
                    type: 'SPNGSimpleOpOr',
                    name: 'mp_finding_or_%d' % anode.data.ident,
                    data: {
                        // ROI_tensor_index: 0,
                        // value_tensor_index: 1,
                        output_set_tag: 'mp_finding_or_%d' %anode.data.ident,
                    },
                }, nin=2, nout=1) for anode in tools.anodes
            ],

            local post_mp_stackers =  [g.pnode({
                type: 'TorchTensorSetStacker',
                name: 'mp_stacker_%d' % anode.data.ident,
                data: {
                    output_set_tag: 'mp_stacked_%d' % anode.data.ident,
                    multiplicity:2,
                },
            }, nin=2, nout=1) for anode in tools.anodes
            ],

            // local dnn_rois = [
            //     g.pnode({
            //         type: 'SPNGDNNROI',
            //         name: '_%d' % anode.data.ident,
            //         data: {
            //         },
            //     }, nin=1, nout=1) for anode in tools.anodes
            // ],

            local tensor_sinks_mp_finding = [g.pnode({
                type: 'TensorFileSink',
                name: 'tfsink_mp_finding_%d' % anode.data.ident,
                data: {
                    outname: 'testout_mp_finding_apa%d.tar' % anode.data.ident,
                    prefix: ''
                },
            }, nin=1, nout=0) for anode in tools.anodes],

            // local run_roi = false,
            local centers = torch_to_tensors + spng_decons + 
                    spng_gaus_apps + apply_loose_rois + threshold_rois + post_decon_replicators + roi_application,
            
            local edges = [
                g.edge(fans_to_stacks, spng_decons[0], 0),
                g.edge(fans_to_stacks, spng_decons[1], 1),
                g.edge(fans_to_stacks, spng_decons[2], 2),

                g.edge(spng_decons[0], post_decon_replicators[0]),
                g.edge(spng_decons[1], post_decon_replicators[1]),
                g.edge(spng_decons[2], post_decon_replicators[2]),

                g.edge(post_decon_replicators[0], apply_loose_rois[0], 0),
                g.edge(post_decon_replicators[1], apply_loose_rois[1], 0),
                g.edge(post_decon_replicators[2], apply_loose_rois[2], 0),

                g.edge(post_decon_replicators[0], spng_gaus_apps[0], 1),
                g.edge(post_decon_replicators[1], spng_gaus_apps[1], 1),
                g.edge(post_decon_replicators[2], spng_gaus_apps[2], 1),

                // g.edge(spng_decons[0], apply_loose_rois[0]),
                // g.edge(spng_decons[1], apply_loose_rois[1]),
                // g.edge(spng_decons[2], apply_loose_rois[2]),


                g.edge(apply_loose_rois[0], threshold_rois[0]),
                g.edge(apply_loose_rois[1], threshold_rois[1]),
                g.edge(apply_loose_rois[2], threshold_rois[2]),

                // g.edge(threshold_rois[0], torch_to_tensors[0]),
                // g.edge(threshold_rois[1], torch_to_tensors[1]),
                // g.edge(threshold_rois[2], torch_to_tensors[2]),

                // g.edge(apply_loose_rois[0], torch_to_tensors[0]),
                // g.edge(apply_loose_rois[1], torch_to_tensors[1]),
                // g.edge(apply_loose_rois[2], torch_to_tensors[2]),

                // g.edge(apply_loose_rois[0], threshold_rois[0]),
                // g.edge(apply_loose_rois[1], threshold_rois[1]),
                // g.edge(apply_loose_rois[2], threshold_rois[2]),

                g.edge(threshold_rois[0], roi_application[0], 0, 0),
                g.edge(threshold_rois[1], roi_application[1], 0, 0),
                g.edge(threshold_rois[2], roi_application[2], 0, 0),

                g.edge(spng_gaus_apps[0], roi_application[0], 0, 1),
                g.edge(spng_gaus_apps[1], roi_application[1], 0, 1),
                g.edge(spng_gaus_apps[2], roi_application[2], 0, 1),

                // g.edge(threshold_rois[0], torch_to_tensors[0]),
                // g.edge(threshold_rois[1], torch_to_tensors[1]),
                // g.edge(threshold_rois[2], torch_to_tensors[2]),
   
                g.edge(roi_application[0], torch_to_tensors[0]),
                g.edge(roi_application[1], torch_to_tensors[1]),
                g.edge(roi_application[2], torch_to_tensors[2]),

                g.edge(torch_to_tensors[0], tensor_sinks[0]),
                g.edge(torch_to_tensors[1], tensor_sinks[1]),
                g.edge(torch_to_tensors[2], tensor_sinks[2]),
            ],

            local mp_finding_centers = torch_to_tensors + spng_decons + 
                    spng_gaus_apps + apply_loose_rois + apply_tight_rois + threshold_rois +
                    [u_unstacker, v_unstacker, w_unstacker, u_unstacker_face1, v_unstacker_face1] +
                    std.flattenArrays(torch_to_tensors_unstacked) +
                    collators_for_mp_finding + collators_for_mp_finding_face1 + torch_to_tensors_mp_finding + mp_finding + mp_finding_replicators +
                    mp_finding_face1 + post_mp_stackers,
            local mp_finding_edges = [
                    g.edge(fans_to_stacks, spng_decons[0], 0),
                    g.edge(fans_to_stacks, spng_decons[1], 1),
                    g.edge(fans_to_stacks, spng_decons[2], 2),

                    // g.edge(spng_decons[0], post_decon_replicators[0]),
                    // g.edge(spng_decons[1], post_decon_replicators[1]),
                    // g.edge(spng_decons[2], post_decon_replicators[2]),


                    // g.edge(spng_decons[0], apply_loose_rois[0]),
                    // g.edge(spng_decons[1], apply_loose_rois[1]),
                    // g.edge(spng_decons[2], apply_loose_rois[2]),

                    // g.edge(post_decon_replicators[0], apply_tight_rois[0], 0),
                    // g.edge(post_decon_replicators[0], apply_loose_rois[0], 1),

                    // g.edge(post_decon_replicators[1], apply_tight_rois[1], 0),
                    // g.edge(post_decon_replicators[2], apply_tight_rois[2], 0),
                    g.edge(spng_decons[0], apply_tight_rois[0], 0),
                    g.edge(spng_decons[1], apply_tight_rois[1], 0),
                    g.edge(spng_decons[2], apply_tight_rois[2], 0),
                    
                    // g.edge(post_decon_replicators[0], spng_gaus_apps[0], 1),
                    // g.edge(post_decon_replicators[1], spng_gaus_apps[1], 1),
                    // g.edge(post_decon_replicators[2], spng_gaus_apps[2], 1),

                    // g.edge(post_decon_replicators[1], apply_loose_rois[1], 1),

                    g.edge(apply_tight_rois[0], threshold_rois[0]),
                    g.edge(apply_tight_rois[1], threshold_rois[1]),
                    g.edge(apply_tight_rois[2], threshold_rois[2]),

                    // g.edge(threshold_rois[0], u_unstacker),
                    // g.edge(threshold_rois[1], v_unstacker),
                    g.edge(threshold_rois[2], w_unstacker),

                    g.edge(threshold_rois[0], mp_finding_replicators[0]),//U
                    g.edge(threshold_rois[1], mp_finding_replicators[1]),//V

                    g.edge(mp_finding_replicators[0], u_unstacker, 0, 0),
                    g.edge(mp_finding_replicators[1], v_unstacker, 0, 0),
                    g.edge(mp_finding_replicators[0], u_unstacker_face1, 1, 0),
                    g.edge(mp_finding_replicators[1], v_unstacker_face1, 1, 0),

                    g.edge(u_unstacker, collators_for_mp_finding[0], 0, 0),
                    g.edge(v_unstacker, collators_for_mp_finding[0], 0, 1),
                    g.edge(w_unstacker, collators_for_mp_finding[0], 0, 2),
                    g.edge(u_unstacker_face1, collators_for_mp_finding_face1[0], 0, 0),
                    g.edge(v_unstacker_face1, collators_for_mp_finding_face1[0], 0, 1),
                    g.edge(w_unstacker, collators_for_mp_finding_face1[0], 1, 2),

                    g.edge(u_unstacker, collators_for_mp_finding[1], 1, 0),
                    g.edge(v_unstacker, collators_for_mp_finding[1], 1, 1),
                    g.edge(w_unstacker, collators_for_mp_finding[1], 2, 2),
                    g.edge(u_unstacker_face1, collators_for_mp_finding_face1[1], 1, 0),
                    g.edge(v_unstacker_face1, collators_for_mp_finding_face1[1], 1, 1),
                    g.edge(w_unstacker, collators_for_mp_finding_face1[1], 3, 2),

                    g.edge(u_unstacker, collators_for_mp_finding[2], 2, 0),
                    g.edge(v_unstacker, collators_for_mp_finding[2], 2, 1),
                    g.edge(w_unstacker, collators_for_mp_finding[2], 4, 2),
                    g.edge(u_unstacker_face1, collators_for_mp_finding_face1[2], 2, 0),
                    g.edge(v_unstacker_face1, collators_for_mp_finding_face1[2], 2, 1),
                    g.edge(w_unstacker, collators_for_mp_finding_face1[2], 5, 2),

                    g.edge(u_unstacker, collators_for_mp_finding[3], 3, 0),
                    g.edge(v_unstacker, collators_for_mp_finding[3], 3, 1),
                    g.edge(w_unstacker, collators_for_mp_finding[3], 6, 2),
                    g.edge(u_unstacker_face1, collators_for_mp_finding_face1[3], 3, 0),
                    g.edge(v_unstacker_face1, collators_for_mp_finding_face1[3], 3, 1),
                    g.edge(w_unstacker, collators_for_mp_finding_face1[3], 7, 2),

                    // g.edge(u_unstacker, collators_for_mp_finding[1], 1, 0),
                    // g.edge(v_unstacker, collators_for_mp_finding[1], 1, 1),
                    // g.edge(w_unstacker, collators_for_mp_finding[1], 2, 2),
                    // g.edge(w_unstacker, collators_for_mp_finding[1], 3, 3),

                    // g.edge(u_unstacker, collators_for_mp_finding[2], 2, 0),
                    // g.edge(v_unstacker, collators_for_mp_finding[2], 2, 1),
                    // g.edge(w_unstacker, collators_for_mp_finding[2], 4, 2),
                    // g.edge(w_unstacker, collators_for_mp_finding[2], 5, 3),

                    // g.edge(u_unstacker, collators_for_mp_finding[3], 3, 0),
                    // g.edge(v_unstacker, collators_for_mp_finding[3], 3, 1),
                    // g.edge(w_unstacker, collators_for_mp_finding[3], 6, 2),
                    // g.edge(w_unstacker, collators_for_mp_finding[3], 7, 3),

                    g.edge(collators_for_mp_finding[0], mp_finding[0]),
                    g.edge(collators_for_mp_finding[1], mp_finding[1]),
                    g.edge(collators_for_mp_finding[2], mp_finding[2]),
                    g.edge(collators_for_mp_finding[3], mp_finding[3]),

                    g.edge(collators_for_mp_finding_face1[0], mp_finding_face1[0]),
                    g.edge(collators_for_mp_finding_face1[1], mp_finding_face1[1]),
                    g.edge(collators_for_mp_finding_face1[2], mp_finding_face1[2]),
                    g.edge(collators_for_mp_finding_face1[3], mp_finding_face1[3]),

                    // g.edge(mp_finding[0], mp_finding_ors[0], 0, 0),
                    // g.edge(mp_finding[1], mp_finding_ors[1], 0, 0),
                    // g.edge(mp_finding[2], mp_finding_ors[2], 0, 0),
                    // g.edge(mp_finding[3], mp_finding_ors[3], 0, 0),

                    // g.edge(mp_finding_face1[0], mp_finding_ors[0], 0, 1),
                    // g.edge(mp_finding_face1[1], mp_finding_ors[1], 0, 1),
                    // g.edge(mp_finding_face1[2], mp_finding_ors[2], 0, 1),
                    // g.edge(mp_finding_face1[3], mp_finding_ors[3], 0, 1),

                    // g.edge(mp_finding_ors[0], torch_to_tensors_mp_finding[0]),
                    // g.edge(mp_finding_ors[1], torch_to_tensors_mp_finding[1]),
                    // g.edge(mp_finding_ors[2], torch_to_tensors_mp_finding[2]),
                    // g.edge(mp_finding_ors[3], torch_to_tensors_mp_finding[3]),

                    g.edge(mp_finding[0], post_mp_stackers[0], 0, 0),
                    g.edge(mp_finding[1], post_mp_stackers[1], 0, 0),
                    g.edge(mp_finding[2], post_mp_stackers[2], 0, 0),
                    g.edge(mp_finding[3], post_mp_stackers[3], 0, 0),

                    g.edge(mp_finding_face1[0], post_mp_stackers[0], 0, 1),
                    g.edge(mp_finding_face1[1], post_mp_stackers[1], 0, 1),
                    g.edge(mp_finding_face1[2], post_mp_stackers[2], 0, 1),
                    g.edge(mp_finding_face1[3], post_mp_stackers[3], 0, 1),

                    g.edge(post_mp_stackers[0], torch_to_tensors_mp_finding[0]),
                    g.edge(post_mp_stackers[1], torch_to_tensors_mp_finding[1]),
                    g.edge(post_mp_stackers[2], torch_to_tensors_mp_finding[2]),
                    g.edge(post_mp_stackers[3], torch_to_tensors_mp_finding[3]),

                    // g.edge(mp_finding[0], torch_to_tensors_mp_finding[0]),
                    // g.edge(mp_finding[1], torch_to_tensors_mp_finding[1]),
                    // g.edge(mp_finding[2], torch_to_tensors_mp_finding[2]),
                    // g.edge(mp_finding[3], torch_to_tensors_mp_finding[3]),

                    g.edge(torch_to_tensors_mp_finding[0], tensor_sinks_mp_finding[0]),
                    g.edge(torch_to_tensors_mp_finding[1], tensor_sinks_mp_finding[1]),
                    g.edge(torch_to_tensors_mp_finding[2], tensor_sinks_mp_finding[2]),
                    g.edge(torch_to_tensors_mp_finding[3], tensor_sinks_mp_finding[3]),

            ],

            ret : g.intern(
                innodes=[fans_to_stacks],
                centernodes=(if do_mp_finding then mp_finding_centers else centers),
                outnodes=(if do_mp_finding then tensor_sinks_mp_finding else tensor_sinks), #std.flattenArrays(tensor_sinks_unstacked),
                edges=(if do_mp_finding then mp_finding_edges else edges),
            ),
        }.ret,
}