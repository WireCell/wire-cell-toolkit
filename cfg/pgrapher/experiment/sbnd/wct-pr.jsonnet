// STEP 2 of the SBND 2-step chain (ai-helper issue 33): pattern recognition,
// standalone (no art), re-runnable on the step-1 tar written by the LArSoft job
// wcls-img-clus-matching.jsonnet:
//
//   TensorFileSource(qlpctree.tar.gz) -> clus_pr MABC (sbnd-pr-stage.jsonnet) -> sink
//
// The PR MABC is built by the SAME definition the LArSoft 1-step chain uses
// (sbnd-pr-stage.jsonnet: pr() defaults = the SBND production operating point,
// the 15-visitor pipeline, the beam gate), so its compiled config equals the
// 1-step's clus_pr up to output file names.
//
//   wire-cell -l stdout -L debug -c pgrapher/experiment/sbnd/wct-pr.jsonnet \
//       --tla-str input=qlpctree.tar.gz --tla-str reality=data
//
// Multi-event: the tar holds every event of the step-1 lar job.  The per-event
// outputs (tracking-pr.root) go to output_dir/<evt_subdir>/ with evt_subdir a
// boost::format template of the set ident (the art event number), default
// 'pr_evt%1%' -- the writers create those directories themselves.
// evt_subdir='' writes ./tracking-pr.root, correct only for a one-event tar.
// RSE: every MABC takes run/subrun/event from the set metadata stamped in step 1
// (rse_from_metadata; rse_from_ident is the fallback), exactly as in the 1-step.
// Bee: all events into one zip (bee_outname), layers named as in the 1-step's
// mabc.zip (clustering-pr, track_fit, shower_track, vertices, mc, tagger_*).
// Truth: on MC the step-1 truth_nu / truth_pf tables ride in the tar; the PR MABC
// publishes them and tracking-pr.root gets T_truth_nu / T_truth_pf.
function(input='qlpctree.tar.gz',
         output_dir='.',
         reality='data',
         evt_subdir='pr_evt%1%',
         enable_tracking_root=true,
         bee_outname='mabc-pr.zip',
         // '' => the terminal sink is a dump_mode no-op; a path => write the
         // post-PR point-cloud tree there.
         pr_tensor_outname='',
         // ai-helper issue 35: record every DL-vertex network call (exact input,
         // payload, decision, MC truth vertex) as T_dlvtx_call / T_dlvtx_cloud.
         dl_vtx_dump=false)

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local params = import 'pgrapher/experiment/sbnd/simparams.jsonnet';
local tools_all = tools_maker(params);
local tools = tools_all {anodes: [tools_all.anodes[n] for n in [0, 1]]};

local clus = import 'pgrapher/experiment/sbnd/clus.jsonnet';
// event_from_ident is clus.jsonnet's required partner of evt_subdir (the per-event
// file names).  It never decides the RSE here: the step-1 metadata carries it, and
// MABC's precedence is metadata > ident > config.
local clus_maker = clus(output_dir=output_dir, rse_from_ident=true, rse_from_metadata=true,
                        event_from_ident=evt_subdir != '',
                        reality=reality, evt_subdir=evt_subdir);
local pr_stage = import 'pgrapher/experiment/sbnd/sbnd-pr-stage.jsonnet';

local source = g.pnode({
    type: 'TensorFileSource',
    name: 'ql_pctree',
    data: { inname: input, prefix: 'clustering_' },
}, nin=0, nout=1);

local bee = {
    type: 'BeeSink',
    name: 'mabc_pr',
    data: { outname: bee_outname, initial_index: 0 },
};

// tagger_bee: the tagger verdict Bee sets (tagger_stm/_tgm/_fc/_lm), which the
// 1-step's art-side labeler_tagger writes, come from TaggerBeeVisitor here.
local pr_node = pr_stage.node(clus_maker, tools.anodes, bee, enable_tracking_root, tagger_bee=true,
                              dl_vtx_dump=dl_vtx_dump);

local sink = g.pnode({
    type: 'TensorFileSink',
    name: 'pr',
    data: {
        outname: if pr_tensor_outname == '' then 'trash-pr.tar.gz' else pr_tensor_outname,
        prefix: 'clustering_',
        dump_mode: pr_tensor_outname == '',
    },
}, nin=1, nout=0);

local graph = g.intern(
    innodes=[source],
    centernodes=[pr_node],
    outnodes=[sink],
    edges=[g.edge(source, pr_node, 0, 0), g.edge(pr_node, sink, 0, 0)],
);

local app = {
    type: 'Pgrapher',
    data: { edges: g.edges(graph) },
};

// WireCellRoot: the PR output visitors, the BDT scorers and SCEFieldTH3 (in the
// DetectorVolumes uses of clus.jsonnet, see wct-pr-perevt.jsonnet).
local cmdline = {
    type: 'wire-cell',
    data: {
        plugins: ['WireCellGen', 'WireCellPgraph', 'WireCellAux', 'WireCellSio',
                  'WireCellSigProc', 'WireCellImg', 'WireCellClus', 'WireCellRoot'],
        apps: ['Pgrapher'],
    },
};

// The Bee sink is named by the PR MABC (and TaggerBeeVisitor) but pr() does not list
// it in its uses -- in the 1-step chain the per-APA MABCs bring the shared sink in --
// so it is configured explicitly here, ahead of its users.
[cmdline, bee] + g.uses(graph) + [app]
