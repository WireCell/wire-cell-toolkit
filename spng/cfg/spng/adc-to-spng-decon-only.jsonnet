// This produces a final graph to input ADC frame and output SPNG signal frame.

local wc = import "wirecell.jsonnet";
//local pg = import "pgraph.jsonnet";
local wrap = import "wrap.jsonnet";

local io = import "spng/io.jsonnet";
local tio = import "spng/torchio.jsonnet";
local detconf = import "spng/detconf.jsonnet";
local detector = import "spng/detector.jsonnet";
local control_js = import "spng/control.jsonnet";

local sg_js = import "spng/subgraphs.jsonnet";

// The TLAs:
//
// @param input The name of a file in WCT "depo file" format, usually .npz.
// @param output The output file name.
// @param detname The name of a supported detector, default "pdhd".
// @param tpcid The TPC ID number.
// @param engine The name of the graph execution engine, default Pgrapher or TbbFlow.
// @param device The name of the device for SPNG nodes, default "cpu" or "gpu", "gpu1", etc. 
// @param dump List of graph points to dump to files.
// @param verbosity The verbosity level for additional logging.
//
// The only required TLA is "input".
//
// Notes:
//
// - Input depos should be arranged to populate the given tpcid.

function(input,
         output="spng.npz",
         detname='pdhd',
         tpcid=0,
         engine='Pgrapher',
         device='cpu',
         dump="",
         verbosity=0)
    
    local controls = control_js(device=device, verbosity=wc.intify(verbosity));
    local det = detector.subset(detconf[detname], [wc.numberify(tpcid)]);

    // Make a Graphviz graph to help understand what SPNGType name and what
    // string to match against the instance name.
    local dump_finders = {
        SPNGThreshold: {
            "cross": "_cross_",
        },
        SPNGTransform: {
            "scale": "_scale_dnnroi",
            "unsqueeze": "_unsqueeze",
            "pre": "_dnnroi_pre",
        },
        SPNGKernelConvolve: {
            "looself": "_dnnroi",
            "wiener": "_wiener",
            "decon": "_group_",
            "gauss": "_gauss",
        },
        SPNGCellViews: {
            "cellviews": "tpc"
        },
        SPNGResampler: {
            "wresample": "_tight"
        },
        SPNGReduce: {
            "applyroi": "_applyroi",
        },
        SPNGRebaseliner: {
            "rebaseline": "_applyroi",
        },
    };
    local is_set(itype) = std.get({'SPNGCellViews': true}, itype, false);
    local make_wrapper = tio.pickle_wrapper("spng-%(itype)s-%(iname)s.pkl", is_set);
    local pg = wrap.replace_pnode(std.split(dump, ","),
                                  dump_finders,
                                  make_wrapper);

    local sg = sg_js(det.tpcs[0], controls.config, pg);


    // True if care about cross view info for the view.
    local rebin = 4;

    local source = io.frame_array_source(input);
    local sink = io.frame_array_any_sink(output);

    local head = sg.frame_to_tdm(extra_name="_TOTDM");
    local tail = sg.tdm_to_frame(extra_name="_FROMTDM", traces_tag='dnnroi');

    local infer = sg.simple_decon(rebin=rebin);
                                      
    local pack = sg.tensor_packer(extra_name="_signals");
    local guts = pg.shuntlines([infer, pack]);
    local body = sg.wrap_bypass(guts);


    local graph = pg.pipeline([source, head, body, tail, sink]);

    pg.main(graph, engine, plugins=["WireCellSpng"])

