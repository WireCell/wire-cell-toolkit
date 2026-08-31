local wc = import "wirecell.jsonnet";
local pg = import 'pgraph.jsonnet';
local det_mod = import "spng/det.jsonnet";
local io = import "spng/io.jsonnet";

local control_mod = import "spng/control.jsonnet";

local detconf = import "spng/detconf.jsonnet";

local sg_js = import "spng/subgraphs.jsonnet";
local roiuniter = import "spng/spng-roiuniter.jsonnet";

// WCT job for the "spcmp" workflow.  Like spngdir.jsonnet it produces three
// signal outputs (splat, osp, spng) from one depo input.
//
// By default (outstage="roiuniter") the "spng" output is produced by the direct
// dense->ROI "roiuniter" DNN graph, exactly as spngdir.jsonnet and
// spng/adc-to-spng.jsonnet do.
//
// For validating the signal processing UPSTREAM of ROI finding, outstage may
// instead select a bare decon stage ("decon", "gauss", "dnnroi", "wiener",
// "tight").  Then the spng chain ends at that stage with no crossviews, ROI DNN,
// apply-ROI or rebaseliner (cf spng/adc-to-spng-decon-only.jsonnet).
//
// The TLAs:
//
// @param input The name of a file in WCT "depo file" format, usually .npz.
// @param model_file Path to the roiuniter TorchScript (.ts) model.  Required
//        when outstage="roiuniter"; ignored otherwise.
// @param outpat The output file name pattern with format variables.
// @param detname The name of a supported detector, default "pdhd".
// @param engine The name of the graph execution engine, default Pgrapher or TbbFlow.
// @param device The name of the device for SPNG nodes, default "cpu" or "gpu", "gpu1", etc.
// @param outstage "roiuniter" (default, the full new SP graph) or a decon stage
//        ("decon", "gauss", "dnnroi", "wiener", "tight") for upstream validation.
// @param rebin Downsample factor applied by the filtered decon stages (gauss/
//        dnnroi/wiener/tight).  Ignored by "decon" and "roiuniter".  Use 1 to
//        keep full time sampling.
//
// The required TLA is "input" (and "model_file" for outstage="roiuniter").
//
// The outpat must include these format variables:
// - %(tier)s will be filled with the output type: splat, osp, spng
//
// Note, this hard-wires use of TPC ID 0.  Input depos should be arranged to
// populate this TPC.
//
function(input,
         model_file="",
         outpat="test-det-%(tier)s.npz",
         detname='pdhd',
         engine='Pgrapher',
         device='cpu',
         outstage='roiuniter',
         rebin=1,
         verbosity=0)

    local tpcids = [0];

    // Per-job-unique prefix for OSP's optional 2D-spectra debug dumps.  Every
    // wire-cell job in the workflow shares a working directory and the OSP dump
    // filename is <prefix>_anode<N>_plane<P>.npz, so a fixed prefix makes
    // concurrent jobs race on the same files (corrupting them).  Derive the
    // prefix from outpat so each job writes its own set.
    local osp_dump_prefix = std.strReplace(outpat % {tier: "osp_dump"}, ".npz", "");

    local controls = control_mod(device=device, verbosity=wc.intify(verbosity));
    local det = detconf.get(detname, tpcids, sp_dump_prefix=osp_dump_prefix);

    // A per-tpc SPNG builder that ends at a bare decon stage (no crossviews/ROI/
    // DNN).  Mirrors spng/adc-to-spng-decon-only.jsonnet: frame_to_tdm ->
    // bypass(decon -> pack) -> tdm_to_frame, tagged by the chosen outstage.
    local decon_only_maker(outstage, rebin) =
        function(tpc, control)
            local sg = sg_js(tpc, control, pg);
            local head = sg.frame_to_tdm(extra_name="_TOTDM");
            local tail = sg.tdm_to_frame(extra_name="_FROMTDM", traces_tag=outstage);
            local infer = sg.simple_decon(rebin=wc.intify(rebin), output=outstage);
            local pack = sg.tensor_packer(extra_name="_signals");
            local body = sg.wrap_bypass(pg.shuntlines([infer, pack]));
            pg.pipeline([head, body, tail]);

    // Select the spng branch builder: the full roiuniter DNN graph (default) or
    // a bare decon stage for upstream validation.
    local spng_maker =
        if outstage == 'roiuniter'
        then roiuniter(model_file, do_transpose=false)
        else decon_only_maker(outstage, rebin);

    // make source node
    local source = io.depo_source(input);

    // Sink a node that has one oport per TPC across the det.
    local sink_det(src, tier) = pg.shuntline(src,  pg.crossline(
        [io.frame_array_sink(outpat % {tier: tier})
         for tpc in det.tpcs]));


    local guts = det_mod(det, controls.config,
                         spng_maker=spng_maker).kitchen_sink;

    local head = pg.pipeline([source, guts.depo_sink]);
    local splat = sink_det(guts.splat_source, "splat");
    local osp = sink_det(guts.osp_source, "osp");
    local spng = sink_det(guts.spng_source, "spng");

    local graph = pg.components([head, splat, osp, spng]);

    pg.main(graph, app=engine,
            plugins=["WireCellSpng", "WireCellSigProc", "WireCellGen", "WireCellPytorch"],
            uses=controls.uses)
