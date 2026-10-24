// CANONICAL in-tree SBND clustering-only job — single source of truth.
// Promoted from sbnd_xin 2026-07-27 (sbnd_xin/docs/64_cfg-sync.md);
// sbnd_xin/wct-clustering.jsonnet is now a one-line re-export of this file.
// Runnable directly:
//   wire-cell -c pgrapher/experiment/sbnd/wct-clustering.jsonnet --tla-...
// (WIRECELL_PATH must contain toolkit/cfg + wire-cell-data).
//
// Run SBND per-APA and all-APA blob clustering — standalone (no LArSoft).
//
// Input:  icluster-apa<N>-active.npz, icluster-apa<N>-masked.npz  (in 'input' dir)
// Output: mabc-<anode>-face0.zip per APA + mabc-all-apa.zip  (in 'output_dir')
//
// Usage (called from run_clus_evt.sh):
//   wire-cell \
//     --tla-str  input=work/evt2 \
//     --tla-code anode_indices='[0,1]' \
//     --tla-str  output_dir=work/evt2 \
//     --tla-code run=0 --tla-code subrun=0 --tla-code event=2 \
//     --tla-str  reality=sim \
//     --tla-code DL=4.0 --tla-code DT=8.8 \
//     --tla-code lifetime=35 --tla-code driftSpeed=1.563 \
//     -c wct-clustering.jsonnet

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';

function(
    input         = '.',
    anode_indices = [0, 1],
    output_dir    = '.',
    run           = 0,
    subrun        = 0,
    event         = 0,
    reality       = 'sim',
    // DL/DT: sbndcode's production diffusion (wcsimsp_sbnd.fcl), restored
    // 2026-07-27, reverting the 6.5781/13.1349 retune of 2026-07-25 (docs/66).
    // These feed simparams.lar -> the sim Drifter; no reco component reads them.
    DL            = 4.0,
    DT            = 8.8,
    lifetime      = 35.0,
    driftSpeed    = 1.563,
    // pr_ml doc 19/20 T6: per-point truth labels (see clus.jsonnet clus_per_face)
    // and the pctree output.  false / '' / 0 => the clus.jsonnet defaults and the
    // all_apa default sink => byte-identical compiled config.
    truth_point_label = false,
    truth_depos = '',
    truth_time_offset_us = 0,
    truth_table = '',
    // doc 20 sec 6.1: fold larg4's negative EM-daughter ids to the ancestor (C++ default false => key omitted => byte-identical)
    truth_fold_negative_ids = false,
    // pr_ml doc 23 sec 10: TruthPointLabeler trajectories from the face's drift volume only / one per Geant4 particle; ON by default for truth-label jobs (owner 2026-10-22; C++ default false; false => key omitted => the pre-fix labels)
    truth_polyline_drift_volume_only = true,
    truth_polyline_per_particle = true,
    // pr_ml doc 24 sec 7: TruthPointLabeler / DepoBlobAssoc point_zero_width -- a depo drifted backward from between the response plane and the wires has zero extent and reached no blob (real charge labelled ghost); ON by default for truth-label jobs (owner 2026-10-23, decision 1A; C++ default false; false => key omitted => the pre-fix labels; unused unless truth_point_label)
    truth_point_zero_width = true,
    // pr_ml doc 25 sec 2.8 (owner 2026-10-24): TruthPointLabeler shared test per PF node (shared_by_pf_node) and
    // PFMerger em_fold (EM daughters of a delta / Michel node merge into it).  C++ default false for both;
    // false => key omitted => byte-identical compiled config; unused unless truth_point_label.
    truth_shared_by_pf_node = false,
    truth_em_fold = false,
    // pr_ml doc 25 sec 4 (owner 2026-10-24): the truth window adds the SP charge time filter width, in us, in quadrature
    // to each depo's longitudinal width (DepoBlobAssoc extra_smear_long; PDHD / PDVD / FD-HD Gaus_wide 0.12 MHz -> 1.326,
    // SBND 0.10 MHz -> 1.592).  C++ default 0; 0 => key omitted => byte-identical; unused unless truth_point_label.
    truth_extra_smear_long_us = 0,
    // pr_ml doc 25 sec 3.2 (owner 2026-10-24): the end flag from the point's own PF node -- 1 / 2 within
    // truth_end_radius_mm (3-D) of the node's start / end (track / delta / Michel: both; showers: start).
    // C++ default false / 10 mm; false => keys omitted => byte-identical; unused unless truth_point_label.
    truth_end_by_pf_node = false,
    truth_end_radius_mm = 10,
    // '' keeps the all-APA stage's default sink; a path writes the post-clustering
    // point-cloud tree there (all_apa tensor_outname, as the Q/L job's save_tensors).
    save_tensors = '',
)
    // Build params inside function so all physics values are TLAs.
    local base = import 'pgrapher/experiment/sbnd/simparams.jsonnet';
    local params = base {
        lar: super.lar {
            DL:          DL * wc.cm2 / wc.s,
            DT:          DT * wc.cm2 / wc.s,
            lifetime:    lifetime * wc.ms,
            drift_speed: driftSpeed * wc.mm / wc.us,
        },
    };

    local tools_all = tools_maker(params);
    local anodes  = [tools_all.anodes[i] for i in anode_indices];
    local nanodes = std.length(anodes);

    local cluster_source(fname) = g.pnode({
        type: 'ClusterFileSource',
        name: fname,
        data: {
            inname: fname,
            anodes: [wc.tn(a) for a in anodes],
        },
    }, nin=0, nout=1, uses=anodes);

    local active_files = ['%s/icluster-apa%d-active.npz' % [input, a.data.ident] for a in anodes];
    local masked_files = ['%s/icluster-apa%d-masked.npz' % [input, a.data.ident] for a in anodes];
    local active_clusters = [cluster_source(f) for f in active_files];
    local masked_clusters = [cluster_source(f) for f in masked_files];

    local clus_mod = import 'clus.jsonnet';
    local clus_maker = clus_mod(
        output_dir=output_dir,
        runNo=run,
        subRunNo=subrun,
        eventNo=event,
        truth_point_label=truth_point_label, truth_depos=truth_depos,
        truth_time_offset=truth_time_offset_us * wc.us, truth_table=truth_table,
        truth_fold_negative_ids=truth_fold_negative_ids,
        truth_polyline_drift_volume_only=truth_polyline_drift_volume_only, truth_polyline_per_particle=truth_polyline_per_particle,
                        truth_point_zero_width=truth_point_zero_width,
                        truth_shared_by_pf_node=truth_shared_by_pf_node, truth_em_fold=truth_em_fold,
                        truth_extra_smear_long=truth_extra_smear_long_us * wc.us,
                        truth_end_by_pf_node=truth_end_by_pf_node, truth_end_radius=truth_end_radius_mm * wc.mm);
    local clus_pipes = [clus_maker.per_apa(anodes[n], dump=false)
                        for n in std.range(0, nanodes - 1)];

    local img_clus_pipe = [g.intern(
        innodes=  [active_clusters[n], masked_clusters[n]],
        centernodes= [],
        outnodes= [clus_pipes[n]],
        edges=    [
            g.edge(active_clusters[n], clus_pipes[n], 0, 0),
            g.edge(masked_clusters[n], clus_pipes[n], 0, 1),
        ]
    ) for n in std.range(0, nanodes - 1)];

    local clus_all = clus_maker.all_apa(anodes, tensor_outname=save_tensors);

    local graph = g.intern(
        innodes=   img_clus_pipe,
        outnodes=  [clus_all],
        edges=     [g.edge(img_clus_pipe[i], clus_all, 0, i)
                    for i in std.range(0, nanodes - 1)]
    );

    local app = {
        type: 'Pgrapher',
        data: { edges: g.edges(graph) },
    };

    local cmdline = {
        type: 'wire-cell',
        data: {
            plugins: ['WireCellGen', 'WireCellPgraph', 'WireCellSio', 'WireCellSigProc',
                      'WireCellImg', 'WireCellRoot', 'WireCellTbb', 'WireCellClus']
                     + (if truth_point_label then ['WireCellTruth'] else []),  // pr_ml doc 19 T6; off => the same list
            apps: ['Pgrapher'],
        },
    };

    [cmdline] + g.uses(graph) + [app]
