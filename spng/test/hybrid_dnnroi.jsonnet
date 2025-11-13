// This jsonnet will use the perfect_pdhd/sp.jsonnet for signal processing 
// And uses DNNROIProcess.cxx for DNN-ROI processing. 
//Jsonnet to validate the DNNROI Processing by taking off the Deconvolution step from SPNG

// function definition
function(
    input_file='tensor_frames.npz',
    output_path='spng%s.tar',     // MUST give %s for gpu and NOT for cpu.
    device='cpu',
    ts_model_file='/nfs/data/1/abashyal/spng/model_files/Pytorch-UNet/ts-model-2.3/unet-l23-cosmic500-e50.ts',
    plane='u',
    anode_index=0,
){
    local g = import 'pgraph.jsonnet',
    local f = import 'pgrapher/common/funcs.jsonnet',
    local wc = import 'wirecell.jsonnet',

    local io = import 'pgrapher/common/fileio.jsonnet',
    local fileio = import 'layers/high/fileio.jsonnet',
    local tools_maker = import 'pgrapher/common/tools.jsonnet',

    local base = import 'perfect_pdhd/simparams.jsonnet',
    local sim_maker = import 'perfect_pdhd/sim.jsonnet',
    local perfect = import 'perfect_pdhd/chndb-base.jsonnet',
    local nf_maker = import 'perfect_pdhd/nf.jsonnet',
    local sp_maker = import 'perfect_pdhd/sp.jsonnet',
    
    local params = base {
      lar: super.lar {
            // Longitudinal diffusion constant
            DL :  6.2 * wc.cm2/wc.s,
            // Transverse diffusion constant
            DT : 16.3 * wc.cm2/wc.s,
            lifetime : 50*wc.ms,
            drift_speed : 1.565*wc.mm/wc.us,
      },
    },
    local tools = tools_maker(params),
    local nanodes = std.length(tools.anodes),

    //we are not doing the simulation but instead using pre-saved framesrouce
    local frame_input = (fileio.frame_tensor_file_source(input_file)),

    //graph is:
    //channelselector --> OmnibusSigProc --> FrameToTorch --> DNNROI (U,V) --> FanIn --> Retagger --> FrameTensorSink

    local selectors = [
        g.pnode({
            type: 'ChannelSelector',
            name: 'channelselect-%s' % tools.anodes[n].name,
            data:{
                anode: wc.tn(tools.anodes[n]),
                tags: ['raw%d' % tools.anodes[n].data.ident],
                channels: std.range(2560*n, 2560*(n+1)-1)
            },
        },nin=1,nout=1) for n in std.range(0, nanodes-1)
    ],

    local sp_override = {
        sparse: false,
        use_roi_debug_mode: true,
        use_roi_refinement: true,
        save_negtive_charge: false,
        use_multi_plane_protection: true,
        process_planes: [0, 1, 2],
        debug_no_frer: false,
        debug_no_wire_filter: false,
        break_roi_loop1_tag: "",
        break_roi_loop2_tag: "",
        shrink_roi_tag: "",
        extend_roi_tag: "",
        // decon_charge_tag: "",
        do_not_mp_protect_traditional: true, // do_not_mp_protect_traditional to 
                                            // make a clear ref, defualt is false
        mp_tick_resolution: 4,
    },

    //create a local config function for that contains argument for SPNGDNNROI nodes
    local config = {
        nchans: 2560,
        nticks: 6000,
        tick_per_slice: 4,
        nchunks: 1,
        input_scale: 0.00025,
        input_offset: 0.0,
        output_scale: 1.0,
        output_offset: 0.0,
        all_preprocess: true,
    },

    local sp = sp_maker(params, tools, sp_override),
    local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes],

    //need to have frame to torch conversion node
    local frame_to_torch = 
        g.pnode({
            type: 'FrameToTorch',
            data: {},
        }, nin=1, nout=1),

    local intags = [
        'loose_lf%d' % tools.anodes[anode_index].data.ident,
        'mp2_roi%d' % tools.anodes[anode_index].data.ident,
        'mp3_roi%d' % tools.anodes[anode_index].data.ident
    ],

    local frame_to_tensorset = g.pnode({
        type: 'FrameToTorchSet',
        data: {
            intags: intags,
            anode: wc.tn(tools.anodes[anode_index]),
            plane: (if plane == 'u' then 0 else if plane == 'v' then 1 else 2),
            nticks: config.nticks,
        },
    }, nin=1, nout=1),

    local torch_packer = g.pnode({
        type: "SPNGTorchPacker",
        name: "torch_tensor_packer_%s" % plane,
        data: {
            'multiplicity':1,
        },
    }, nin=1, nout=1),

    //define the TorchService node helper
    local SPNGTorchService={
        type: "SPNGTorchService",
        name: "dnnroi-service",
        data:{
            model: ts_model_file,
            device:device,
        }
    },



    //create the toolkit pipe for dnnroi
    local dnnroi = import 'spng_helpers/dnn-roi.jsonnet',
    local toolkit_pipe = g.pipeline([
        sp_pipes[0],
        //frame_to_torch,
        frame_to_tensorset,
        //torch_packer,
        dnnroi(SPNGTorchService, plane, prefix="dnnroi", config=config)
    ]),

    //frame output node
    local frame_output = 
        g.pnode(
            {
                type: "TensorFileSink",
                name: "framesink-dnnroi-%s" % plane,
                data: {
                    outname: output_path % plane,
                    prefix: ''
                },
            }, nin=1, nout=0),

    //internal graph of input, processing and output
    local graph = g.intern(
        innodes = [frame_input],
        centernodes = selectors + [toolkit_pipe],
        outnodes = [frame_output],
        edges = [
            g.edge(frame_input,selectors[0]),
            g.edge(selectors[0],toolkit_pipe),
            g.edge(toolkit_pipe,frame_output),
        ]
    ),

    local app = {
        type: 'Pgrapher',
        data:{
            edges: g.edges(graph),
        },
    },

    local cmdline = {
        type: "wire-cell",
        data:{
            plugins: [
                "WireCellGen",
                "WireCellPgraph",
                "WireCellSio",
                "WireCellSigProc",
                "WireCellSpng"],
            apps: ["Pgrapher"]
        }
    },

    ret: [cmdline] + g.uses(graph) + [app]

}.ret
