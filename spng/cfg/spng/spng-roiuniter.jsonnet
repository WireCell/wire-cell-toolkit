// A per-TPC SPNG builder for the direct dense->ROI "roiuniter" DNN graph.
//
// This produces a function (tpc, control) -> a 1->1 node that consumes an ADC
// IFrame from one anode and produces a signal IFrame.  It is a drop-in for the
// spng_maker parameter of spng/det.jsonnet and mirrors the graph built by
// spng/adc-to-spng.jsonnet: it runs subgraphs.dnnroi_inference_cat() which
// deconvolves, concatenates the per-view gauss and dense images along the
// channel axis, runs ONE DNN forward over the concatenated dense image, unbins
// the resulting ROI and applies it to the concatenated gauss.  There are no
// crossviews, no MP2/MP3 cell views and no per-view DNN calls.
//
// @param model_file Path to the roiuniter TorchScript (.ts) model.  Required.
// @param do_transpose Whether the model expects a (ntick,nchan) transpose.  The
//        current roiuniter models do not, matching adc-to-spng.jsonnet.
// @param rebin Downsample factor of the "dense" branch fed to the DNN.
// @param scale Multiplicative scale applied to the dense DNN input.

local wc = import "wirecell.jsonnet";
local real_pg = import "pgraph.jsonnet";
local sg_js = import "spng/subgraphs.jsonnet";

function(model_file, do_transpose=false, rebin=4, scale=1.0/4000, pg=real_pg)
    function(tpc, control)
        local sg = sg_js(tpc, control, pg);
        local head = sg.frame_to_tdm(extra_name="_TOTDM");
        local tail = sg.tdm_to_frame(extra_name="_FROMTDM");  // traces_tag="signal"
        local infer = sg.dnnroi_inference_cat(modelfile=model_file,
                                              rebin=rebin, scale=scale,
                                              do_transpose=do_transpose);
        local pack = sg.tensor_packer(multiplicity=1, extra_name="_signals");
        local body = sg.wrap_bypass(pg.shuntlines([infer, pack]));
        pg.pipeline([head, body, tail])
