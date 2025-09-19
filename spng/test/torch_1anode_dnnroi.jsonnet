// TODO -- brief descrip

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';
local spng_filters = import 'spng_filters.jsonnet';

function(
    tools, debug_force_cpu=false,
    ts_model_file="/nfs/data/1/abashyal/spng/spng_dev_050525/Pytorch-UNet/ts-model-2.3/unet-l23-cosmic500-e50.ts",
) {
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

        // local make_replicator_post_tight(anode, iplane) = {
        //     ret : g.pnode({
        //         type: 'TorchTensorSetReplicator',
        //         name: 'post_tight_replicator_%d_%d' % [anode.data.ident, iplane],
        //         data: {
        //             multiplicity: 5,
        //         }
        //     }, nin=1, nout=5, uses=[anode])
        // }.ret,

        // local make_replicator_post_gaus(anode, iplane) = {
        //     ret : g.pnode({
        //         type: 'TorchTensorSetReplicator',
        //         name: 'post_gaus_replicator_%d_%d' % [anode.data.ident, iplane],
        //         data: {multiplicity: 2,}

        //     }, nin=1, nout=2, uses=[anode])
        // }.ret,
        
        // local make_replicator_post_gaus_simple(anode, iplane) = {
        //     ret : g.pnode({
        //         type: 'TorchTensorSetReplicator',
        //         name: 'post_gaus_replicator_%d_%d' % [anode.data.ident, iplane],
        //         data: {multiplicity: (if iplane > 1 then 2 else 6)}
        //     }, nin=1, nout=(if iplane > 1 then 2 else 6), uses=[anode])
        // }.ret,

    
        
   
        stacked_spng : {
            //TODO: Model information should be coming from wct-framesrouce jsonnet
            local SPNGTorchService = {
                type: "SPNGTorchService",
                name: "dnnroi",
                data:{
                    model: ts_model_file,//"/nfs/data/1/abashyal/spng/spng_dev_050525/Pytorch-UNet/ts-model-2.3/unet-l23-cosmic500-e50.ts",
                    device: (if debug_force_cpu then 'cpu' else 'gpu'),
                }
            },
            local tf_fans = make_fanout(tools.anodes[0]),
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
            local w_unstacker =  g.pnode({
                type: 'TorchTensorSetUnstacker',
                name: 'w_unstacker',
                data: {
                    output_set_tag: 'w_stacked',
                    multiplicity:8,
                },
            }, nin=1, nout=8),


            // local torch_to_tensors = [
            //     g.pnode({
            //         type: 'TorchToTensor',
            //         name: 'torchtotensor_%s' % [i],
            //         data: {},
            //     }, nin=1, nout=1) for i in ['u', 'v', 'w']
            // ],

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
            local torch_to_tensors = [
                g.pnode({
                    type: 'TorchToTensor',
                    name: 'torchtotensor_%s' % plane,
                    data: {},
                }, nin=1, nout=1) for plane in ['u', 'v', 'w1', 'w2']
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

            // local tensor_sinks = [g.pnode({
            //     type: 'TensorFileSink',
            //     name: 'tfsink_%s' % i,
            //     data: {
            //         outname: 'testout_fan_plane_%s.tar' % i,
            //         prefix: ''
            //     },
            // }, nin=1, nout=0) for i in ['u', 'v', 'w']],

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
            } for iplane in std.range(0,1)],

            local do_loose_roi_filters = [g.pnode({
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
            uses=[torch_roi_loose_spectra[iplane]]) for iplane in std.range(0,1)],


            local torch_roi_tight_spectra = [{
                    type: "Torch1DSpectrum",
                    name: "torch_1dspec_roi_tight_%d" % iplane,
                    uses: [filters.wiener_tight_filters[iplane], filters.ROI_tight_lf],
                    data: {
                        spectra: [
                            wc.tn(filters.wiener_tight_filters[iplane]),
                            wc.tn(filters.ROI_tight_lf),
                        ],
                        debug_force_cpu: debug_force_cpu,
                    },
            } for iplane in std.range(0,2)],

            local do_tight_roi_filters = [g.pnode({
                type: 'SPNGApply1DSpectrum',
                name: 'spng_tight_roi_plane%d' % iplane,
                data: {
                    base_spectrum_name: wc.tn(torch_roi_tight_spectra[
                        if iplane < 3 then iplane else 2
                    ]),
                    dimension: 2,
                    // target_tensor: 'HfGausWide',
                    output_set_tag: 'ROITight',
                },
            },
            nin=1, nout=1,
            uses=[torch_roi_tight_spectra[if iplane < 3 then iplane else 2]]) for iplane in std.range(0,3)],


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
            } for iplane in std.range(0, 3)],
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
                for iplane in std.range(0, 3)
            ],

            local post_decon_replicators = [
                g.pnode({
                    type: 'TorchTensorSetReplicator',
                    name: 'post_decon_replicator_%d' % iplane,
                    data: {multiplicity: (if iplane < 2 then 3 else 2),}

                }, nin=1, nout=(if iplane < 2 then 3 else 2))
                for iplane in std.range(0, 3)
            ],


            local post_tight_replicators = [
                g.pnode({
                    type: 'TorchTensorSetReplicator',
                    name: 'post_tight_replicator_%d' % iplane,
                    data: {multiplicity: (if iplane > 1 then 3 else 2),}

                }, nin=1, nout=(if iplane > 1 then 3 else 2))
                for iplane in std.range(0,3)
            ],

            local do_gaus_filters = [g.pnode({
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
            uses=[filters.torch_gaus_filter]) for iplane in std.range(0,3)],

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
            nin=1, nout=1, uses=[]) for iplane in std.range(0,3)],

            local collators_for_mp_finding = [
                g.pnode({
                    type: 'TorchTensorSetCollator',
                    name: 'collate_mp_finding_%s' % plane,
                    data: {
                        output_set_tag: 'collated_mp_finding_%s' %  plane,
                        multiplicity: 4,
                    },
                }, nin=4, nout=1) for plane in ['u', 'v']],


            local collators_for_dnn_roi = [
                g.pnode({
                    type: 'TorchTensorSetCollator',
                    name: 'collate_dnn_roi_%s' % plane,
                    data: {
                        output_set_tag: 'collated_dnn_roi_%s' %  plane,
                        multiplicity: 2,
                    },
                }, nin=4, nout=1) for plane in ['u', 'v']],

            local mp_finding = [
                g.pnode({
                    type: 'SPNGNoTileMPCoincidence',
                    name: 'mp_finding_%s' % plane,
                    data: {
                        anode: wc.tn(tools.anodes[0]),
                        target_plane_index: 0,
                        aux_plane_l_index: 2,
                        aux_plane_m_index: 1,
                        debug_force_cpu: debug_force_cpu,
                        // output_torch_name: "mp_finding_%d_tensors.pt" % anode.data.ident,
                    },
                }, nin=1, nout=1, uses=[tools.anodes[0]]) for plane in ['u', 'v']
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
                }, nin=2, nout=1) for plane in ['u', 'v', 'w1', 'w2']
            ],

            local dnn_rois = [
                g.pnode({
                    type: 'SPNGDNNROI',
                    name: 'dnnroi_%s' % plane,
                    data: {
                        
                        plane: plane,
                        input_scale: 1.0/4000,
                        input_offset: 0.0,
                        mask_threshold: 0.5,
                        output_scale: 1.0,
                        output_offset: 0.0,
                        nchunks: 4,
                        forward: wc.tn(SPNGTorchService),

                    },
                }, nin=1, nout=1, uses=[SPNGTorchService]) for plane in ['u', 'v']
            ],

            local tensor_sinks = [g.pnode({
                type: 'TensorFileSink',
                name: 'tfsink_mp_finding_%s' % plane,
                data: {
                    outname: 'testout_mp_finding_%s.tar' % plane,
                    prefix: ''
                },
            }, nin=1, nout=0) for  plane in ['u', 'v', 'w1', 'w2']],


            local mp_finding_centers = torch_to_tensors + spng_decons +  post_decon_replicators +
                    do_gaus_filters + do_loose_roi_filters + do_tight_roi_filters +
                    threshold_rois + post_tight_replicators + roi_application +
                    collators_for_mp_finding + torch_to_tensors + mp_finding
                    + collators_for_dnn_roi + dnn_rois + tensor_sinks,
                    
            local mp_finding_edges = [
                    g.edge(tf_fans, spng_decons[0], 0),
                    g.edge(tf_fans, spng_decons[1], 1),
                    g.edge(tf_fans, spng_decons[2], 2),
                    g.edge(tf_fans, spng_decons[3], 3),


                    g.edge(spng_decons[0], post_decon_replicators[0]), //U
                    g.edge(spng_decons[1], post_decon_replicators[1]), //V
                    g.edge(spng_decons[2], post_decon_replicators[2]), //W1
                    g.edge(spng_decons[3], post_decon_replicators[3]), //W2

                    g.edge(post_decon_replicators[0], do_tight_roi_filters[0], 0),//U
                    g.edge(post_decon_replicators[0], do_loose_roi_filters[0], 1),//U
                    g.edge(post_decon_replicators[0], do_gaus_filters[0], 2),  //U

                    g.edge(post_decon_replicators[1], do_tight_roi_filters[1], 0),//V
                    g.edge(post_decon_replicators[1], do_loose_roi_filters[1], 1),//V
                    g.edge(post_decon_replicators[1], do_gaus_filters[1], 2),  //V

                    g.edge(post_decon_replicators[2], do_tight_roi_filters[2], 0),//W1
                    g.edge(post_decon_replicators[2], do_gaus_filters[2], 1),  //W1

                    g.edge(post_decon_replicators[3], do_tight_roi_filters[3], 0),//W2
                    g.edge(post_decon_replicators[3], do_gaus_filters[3], 1),  //W2

                    g.edge(do_tight_roi_filters[0], threshold_rois[0]),
                    g.edge(do_tight_roi_filters[1], threshold_rois[1]),
                    g.edge(do_tight_roi_filters[2], threshold_rois[2]),
                    g.edge(do_tight_roi_filters[3], threshold_rois[3]),

                    g.edge(threshold_rois[0], post_tight_replicators[0]),
                    g.edge(threshold_rois[1], post_tight_replicators[1]),
                    g.edge(threshold_rois[2], post_tight_replicators[2]),
                    g.edge(threshold_rois[3], post_tight_replicators[3]),

                    g.edge(post_tight_replicators[2], roi_application[2], 2, 0),//W1
                    g.edge(post_tight_replicators[3], roi_application[3], 2, 0),//W2
                    g.edge(do_gaus_filters[2], roi_application[2], 0, 1),//W1
                    g.edge(do_gaus_filters[3], roi_application[3], 0, 1),//W2

                    g.edge(post_tight_replicators[0], collators_for_mp_finding[0], 0, 0),//U
                    g.edge(post_tight_replicators[1], collators_for_mp_finding[0], 0, 1),//U
                    g.edge(post_tight_replicators[2], collators_for_mp_finding[0], 0, 2),//U
                    g.edge(post_tight_replicators[3], collators_for_mp_finding[0], 0, 3),//U

                    g.edge(post_tight_replicators[0], collators_for_mp_finding[1], 1, 0),//V
                    g.edge(post_tight_replicators[1], collators_for_mp_finding[1], 1, 1),//V
                    g.edge(post_tight_replicators[2], collators_for_mp_finding[1], 1, 2),//V
                    g.edge(post_tight_replicators[3], collators_for_mp_finding[1], 1, 3),//V

                    g.edge(collators_for_mp_finding[0], mp_finding[0]),
                    g.edge(collators_for_mp_finding[1], mp_finding[1]),

                    //Now, we have MP-found tensors for the U Plane for each APA
                    //Need something to collate MP-finding & loose lf
                    //then pass to DNNROI 
                    //For now let's save the output of DNNROI to a file

                    g.edge(mp_finding[0], collators_for_dnn_roi[0], 0, 0),//U
                    g.edge(do_loose_roi_filters[0], collators_for_dnn_roi[0], 0, 1),//U

                    g.edge(mp_finding[1], collators_for_dnn_roi[1], 0, 0),//V
                    g.edge(do_loose_roi_filters[1], collators_for_dnn_roi[1], 0, 1),//V

                    // g.edge(collators_for_dnn_roi[0], torch_to_tensors[0]),
                    g.edge(collators_for_dnn_roi[0], dnn_rois[0]),
                    g.edge(dnn_rois[0], roi_application[0], 0, 0),//U
                    g.edge(do_gaus_filters[0], roi_application[0], 0, 1),//U
                    g.edge(roi_application[0], torch_to_tensors[0]),

                    // g.edge(collators_for_dnn_roi[1], torch_to_tensors[1]),
                    g.edge(collators_for_dnn_roi[1], dnn_rois[1]),
                    g.edge(dnn_rois[1], roi_application[1], 0, 0),//V
                    g.edge(do_gaus_filters[1], roi_application[1], 0, 1),//V
                    g.edge(roi_application[1], torch_to_tensors[1]),


                    g.edge(roi_application[2], torch_to_tensors[2]),
                    g.edge(roi_application[3], torch_to_tensors[3]),

                    g.edge(torch_to_tensors[0], tensor_sinks[0]),
                    g.edge(torch_to_tensors[1], tensor_sinks[1]),
                    g.edge(torch_to_tensors[2], tensor_sinks[2]),
                    g.edge(torch_to_tensors[3], tensor_sinks[3]),

            ],

            ret : g.intern(
                innodes=[tf_fans],
                centernodes=(mp_finding_centers),
                outnodes=(tensor_sinks),
                edges=(mp_finding_edges),
            ),
        }.ret,
}