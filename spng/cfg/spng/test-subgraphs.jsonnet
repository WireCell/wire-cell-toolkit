
// This produces a final graph to input depos and output crossview tensors

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

local io = import "spng/io.jsonnet";
local detconf = import "spng/detconf.jsonnet";
local detector = import "spng/detector.jsonnet";
local control_js = import "spng/control.jsonnet";

local sg_js = import "spng/subgraphs.jsonnet";

function(input, output="test-subgraphs.npz", detname='pdhd', tpcids=[], engine='Pgrapher', device='cpu', verbosity=0)
    
    local controls = control_js(device=device, verbosity=wc.intify(verbosity));
    local det = detector.subset(detconf[detname], tpcids);
    local sg = sg_js(det.tpcs[0], controls.config, pg);

    // True if care about cross view info for the view.
    local crossed_views = [1,1,0];
    local all_views = wc.iota(std.length(crossed_views));
    local rebin = 4;

    local dnnroi_model_file = "unet-l23-cosmic500-e50.ts";

    local source = io.frame_array_source(input);
    local sink = io.frame_array_any_sink(output);

    local head = sg.frame_to_tdm(extra_name="_TOTDM");
    local tail = sg.tdm_to_frame(extra_name="_FROMTDM");

    local infer = sg.dnnroi_inference(modelfile=dnnroi_model_file,
                                          rebin=rebin, crossed_views=crossed_views);
    local pack = sg.tensor_packer(extra_name="_signals");
    local guts = pg.shuntlines([infer, pack]);
    local body = sg.wrap_bypass(guts);


    local graph = pg.pipeline([source, head, body, tail, sink]);

    pg.main(graph, engine, plugins=["WireCellSpng"])


