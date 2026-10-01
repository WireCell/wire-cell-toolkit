// dune10kt-1x2x6/wct-pr-perevt.jsonnet -- per-event neutrino pattern-recognition job for DUNE FD-HD
// 1x2x6 (pr.jsonnet; wcp-porting-validation wcfm/docs/25).
//
// Input: the all-TPC pctree that wct-clustering.jsonnet writes (save_tensors, TensorFileSink prefix
// 'clustering_').  Output in output_dir: tracking-pr.root (T_tagger, T_kine, trajectory/dQ/dx trees),
// mabc-pr.zip (Bee: clustering-pr, track_fit, shower_track, vertices, the "mc" particle flow) and
// calib-pr-evt<N>.json (pr_display).
//
//   wcsonnet -A input=work/000001_1311/pctree-evt1311.tar.gz -S 'anode_indices=[0,2]' \
//            -A output_dir=work/000001_1311/pr -S run=1 -S subrun=1 -S event=1311 \
//            -o pr.json pgrapher/experiment/dune10kt-1x2x6/wct-pr-perevt.jsonnet
//   wire-cell -l stderr -c pr.json
//
// anode_indices must be the anodes the clustering ran with (the pctree's ctpc_a<A>f<F> sets).
// Light-less SIMULATION only: the first stage bundles the whole event (ClusteringBundleNoLight).
local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local P = import 'pgrapher/experiment/dune10kt-1x2x6/clus_params.jsonnet';
local tools_all = tools_maker(P.params);

// SBND's production order without the Q/L, cosmic-tagger and BDT stages (pr.jsonnet header):
// bundle_no_light before switch_scope; fiducialutils before any tagger; tagger_output after
// tracking_visitor.  protect_bundle (+ its steiner_refresh) is NOT in the default: on FD it splits
// a muon at the z gap between APA columns and the neutrino tagger then fits only the downstream
// piece -- vertex 128 cm (1311) and 115 cm (1308) from the true vertex, 0.6-0.7 cm without it
// (wcfm/docs/25).  Name 'protect_bundle', 'steiner_refresh' (in that order, after
// tagger_check_fc) to run SBND's full order.
local default_pipeline = [
    'bundle_no_light', 'switch_scope', 'steiner', 'fiducialutils', 'tagger_check_fc',
    'tagger_check_neutrino', 'tracking_visitor', 'tagger_output', 'pr_display',
];

function(
    input = 'pctree.tar.gz',
    anode_indices = std.range(0, std.length(tools_all.anodes) - 1),
    output_dir = '',
    run = 1,
    subrun = 1,
    event = 1,
    pipeline_names = default_pipeline,
    trackfitting_config = 'pgrapher/experiment/dune10kt-1x2x6/dune10kt_track_fitting.json',
    pr_bee = true,
    bee_detector = P.bee_detector,
    // '' keeps the inert dump_mode sink (trash-pr.tar.gz); a path writes the post-PR pctree.
    save_tensors = '',
)

local anodes = [tools_all.anodes[i] for i in anode_indices];

local source = g.pnode({
    type: 'TensorFileSource',
    name: 'clus_pctree',
    data: { inname: input, prefix: 'clustering_' },
}, nin=0, nout=1);

local pds = (import 'pgrapher/experiment/dune10kt-1x2x6/particle_dataset.jsonnet')();
local pr_maker = (import 'pgrapher/experiment/dune10kt-1x2x6/pr.jsonnet')(
    output_dir=output_dir, runNo=run, subRunNo=subrun, eventNo=event);
local pr = pr_maker.pr(anodes, pds, pipeline_names, dump=true, tensor_outname=save_tensors,
                       pr_bee=pr_bee, trackfitting_config_file=trackfitting_config,
                       bee_detector=bee_detector);

local graph = g.pipeline([source, pr]);

local app = { type: 'Pgrapher', data: { edges: g.edges(graph) } };
local cmdline = {
    type: 'wire-cell',
    data: {
        plugins: ['WireCellGen', 'WireCellPgraph', 'WireCellAux', 'WireCellSio',
                  'WireCellSigProc', 'WireCellImg', 'WireCellClus', 'WireCellRoot'],
        apps: ['Pgrapher'],
    },
};

[cmdline] + g.uses(graph) + [app]
