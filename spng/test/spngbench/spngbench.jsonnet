// Factored SPNG-vs-OSP compute benchmark configuration.
//
// This is a deliberately *factored* variant of spngdir.jsonnet.  Where
// spngdir.jsonnet builds the "kitchen sink" graph that produces splat + OSP +
// SPNG signals from one depo input in a single wire-cell job, this config
// builds exactly ONE signal-processing chain (OSP or SPNG) so the two can be
// timed in isolation.
//
// Both chains share the same simulation front-end (drift + detsim) which feeds
// ADC frames to the chosen SP chain.  The per-node Timer log lines let the
// shared sim cost be separated from the SP cost during analysis, so timing the
// full depos->signals chain here is fair for an OSP-vs-SPNG comparison.
//
// The det.jsonnet module already exposes the factored pipelines we need:
//   depos_to_osp  = drift -> detsim -> OSP
//   depos_to_spng = drift -> detsim -> SPNG
//
// TLAs:
//
// @param input      Depo file (WCT depo format, usually .npz).
// @param model_file Path to the roiuniter TorchScript (.ts) ROI DNN model.
// @param output     Output signal frame file (single TPC; PDHD uses TPC 0).
// @param which      Which SP chain to build: "osp" or "spng".
// @param detname    Supported detector name, default "pdhd".
// @param device     Device for GPU-capable nodes: "cpu", "gpu", "gpu0", "gpu1".
// @param engine     Graph execution engine: "Pgrapher" (serial) or "TbbFlow".
// @param wc_cores   Wire-Cell TBB parallelism.  Only used when engine=TbbFlow;
//                   sets TbbDataFlowGraph.max_threads (0 = unlimited).
// @param verbosity  SPNG verbosity.
//
// NOTE: torch intra-op threads are NOT set here.  They are controlled at run
// time via the OMP_NUM_THREADS environment variable (honored by both the OSP
// TorchService and libtorch used by SPNG).
//
// This hard-wires TPC ID 0 (PDHD).  Input depos should populate that TPC.

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local io = import "spng/io.jsonnet";
local control_mod = import "spng/control.jsonnet";
local detconf = import "spng/detconf.jsonnet";
local det_mod = import "spng/det.jsonnet";
local roiuniter = import "spng/spng-roiuniter.jsonnet";

function(input,
         model_file,
         output="spngbench-out.npz",
         which="spng",
         detname='pdhd',
         device='cpu',
         engine='Pgrapher',
         wc_cores=1,
         verbosity=0)

    assert which == "osp" || which == "spng"
           : "spngbench: 'which' must be 'osp' or 'spng', got " + which;

    local tpcids = [0];

    // Per-job-unique prefix for OSP's optional 2D-spectra debug dumps (see
    // spngdir.jsonnet).  Derive from output so concurrent jobs never collide.
    local osp_dump_prefix = std.strReplace(output, ".npz", "") + "_osp_dump";

    local controls = control_mod(device=device, verbosity=wc.intify(verbosity));
    local det = detconf.get(detname, tpcids, sp_dump_prefix=osp_dump_prefix);

    local source = io.depo_source(input);

    local spng_maker = roiuniter(model_file, do_transpose=false);
    local guts = det_mod(det, controls.config, spng_maker=spng_maker);

    // The factored [1]IDepoSet -> IFrame[ntpcs] chain for the chosen SP kind.
    local chain = if which == "osp" then guts.depos_to_osp else guts.depos_to_spng;

    // Terminate every SP output port with a frame file sink.  PDHD uses a single
    // TPC so there is one output port and one sink; the crossline pattern keeps
    // this correct if a multi-TPC detector is used later.
    local ntpcs = std.length(det.tpcs);
    local sinks = [
        io.frame_array_sink(if ntpcs == 1 then output
                            else std.strReplace(output, ".npz", "-tpc" + std.toString(t) + ".npz"))
        for t in wc.iota(ntpcs)
    ];

    local graph =
        if ntpcs == 1
        then pg.pipeline([source, chain, sinks[0]])
        else pg.intern(innodes=[pg.pipeline([source, chain])],
                       outnodes=sinks,
                       edges=[pg.edge(chain, sinks[t], t, 0) for t in wc.iota(ntpcs)]);

    // When running under TbbFlow, add a configured TbbDataFlowGraph component so
    // its max_threads bounds Wire-Cell's node-level parallelism.  pg.main does
    // not configure this component, so we inject it via the 'uses' list.
    // summary>=1 makes TbbDataFlowGraph log per-node "calls=.. time=.. core=.."
    // timing (at debug for 1, info for >=2).  spngbench parses these for the
    // per-node accounting under TbbFlow (Pgrapher emits its own "Timer:" lines).
    local tbb_dfp = {
        type: "TbbDataFlowGraph",
        name: "",
        data: { max_threads: wc.intify(wc_cores), summary: 1 },
    };
    local extra_uses = if engine == "TbbFlow" then [tbb_dfp] else [];

    pg.main(graph, app=engine,
            plugins=["WireCellSpng", "WireCellSigProc", "WireCellGen", "WireCellPytorch"],
            uses=controls.uses + extra_uses)
