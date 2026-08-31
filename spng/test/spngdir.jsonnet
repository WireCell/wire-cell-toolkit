local wc = import "wirecell.jsonnet";
local pg = import 'pgraph.jsonnet';
local det_mod = import "spng/det.jsonnet";
local io = import "spng/io.jsonnet";

local control_mod = import "spng/control.jsonnet";

local detconf = import "spng/detconf.jsonnet";

local roiuniter = import "spng/spng-roiuniter.jsonnet";

// This is a variant of test-det.jsonnet that locks in using the "kitchen sink" graph and hooks up I/O.
//
// The "spng" output is produced by the direct dense->ROI "roiuniter" DNN graph
// (cf spng/adc-to-spng.jsonnet), not the older crossviews/MP2/MP3/DNNROI chain.
//
// The TLAs:
//
// @param input The name of a file in WCT "depo file" format, usually .npz.
// @param model_file Path to the roiuniter TorchScript (.ts) ROI DNN model.
// @param outpat The output file name pattern with format variables.
// @param detname The name of a supported detector, default "pdhd".
// @param engine The name of the graph execution engine, default Pgrapher or TbbFlow.
// @param device The name of the device for SPNG nodes, default "cpu" or "gpu", "gpu1", etc.
//
// The required TLAs are "input" and "model_file".
//
// The outpat must include these format variables:
// - %(tier)s will be filled with the output type: splat, osp, spng
//
// Note, this hard-wires use of TPC ID 0.  Input depos should be arranged to
// populate this TPC.
//
function(input,
         model_file,
         outpat="test-det-%(tier)s.npz",
         detname='pdhd',
         engine='Pgrapher',
         device='cpu',
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

    // make source node
    local source = io.depo_source(input);

    // Sink a node that has one oport per TPC across the det.
    local sink_det(src, tier) = pg.shuntline(src,  pg.crossline(
        [io.frame_array_sink(outpat % {tier: tier})
         for tpc in det.tpcs]));


    local spng_maker = roiuniter(model_file, do_transpose=false);
    local guts = det_mod(det, controls.config, spng_maker=spng_maker).kitchen_sink;

    local head = pg.pipeline([source, guts.depo_sink]);
    local splat = sink_det(guts.splat_source, "splat");
    local osp = sink_det(guts.osp_source, "osp");
    local spng = sink_det(guts.spng_source, "spng");

    local graph = pg.components([head, splat, osp, spng]);

    pg.main(graph, app=engine,
            plugins=["WireCellSpng", "WireCellSigProc", "WireCellGen", "WireCellPytorch"],
            uses=controls.uses)

