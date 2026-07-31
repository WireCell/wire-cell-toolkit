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
        eventNo=event);
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

    local clus_all = clus_maker.all_apa(anodes);

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
                      'WireCellImg', 'WireCellRoot', 'WireCellTbb', 'WireCellClus'],
            apps: ['Pgrapher'],
        },
    };

    [cmdline] + g.uses(graph) + [app]
