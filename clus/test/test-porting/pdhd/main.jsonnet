local g = import "pgraph.jsonnet";
local f = import "pgrapher/common/funcs.jsonnet";
local wc = import "wirecell.jsonnet";

local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';

function(input="")

    // local input = std.extVar('input');

    local reality = 'data';
    local data_params = import 'pgrapher/experiment/pdhd/params.jsonnet';
    local simu_params = import 'pgrapher/experiment/pdhd/simparams.jsonnet';
    local params = if reality == 'data' then data_params else simu_params;
    local tools_maker = import 'pgrapher/common/tools.jsonnet';
    local tools = tools_maker(params);
    local anodes = tools.anodes;
    local nanodes = std.length(tools.anodes);


    local cluster_source(fname) = g.pnode({
        type: "ClusterFileSource",
        name: fname,
        data: {
            inname: fname,
            anodes: [wc.tn(a) for a in anodes],
        }
    }, nin=0, nout=1, uses=anodes);
    local active_files = [ "%s/clusters-apa-apa%d-ms-active.tar.gz"%[input, a.data.ident] for a in anodes];
    local masked_files = [ "%s/clusters-apa-apa%d-ms-masked.tar.gz"%[input, a.data.ident] for a in anodes];
    local active_clusters = [cluster_source(f) for f in active_files];
    local masked_clusters = [cluster_source(f) for f in masked_files];

    local clus = import 'clus.jsonnet';
    local clus_maker = clus();
    // local clus_pipes = [clus_maker.per_face(tools.anodes[n]) for n in std.range(0, std.length(tools.anodes) - 1)];
    local clus_pipes = [clus_maker.per_apa(tools.anodes[n], dump=false) for n in std.range(0, std.length(tools.anodes) - 1)];

    local img_clus_pipe = [g.intern(
        innodes = active_clusters + masked_clusters,
        centernodes = [],
        outnodes = [clus_pipes[n]],
        edges = [
            g.edge(active_clusters[n], clus_pipes[n], 0, 0),
            g.edge(masked_clusters[n], clus_pipes[n], 0, 1),
        ]
    )
                           for n in std.range(0, std.length(tools.anodes) - 1)];

    local clus_all_apa = clus_maker.all_apa(tools.anodes);

    local parallel_graph = 
        {
            local begin = img_clus_pipe,
            local end = clus_maker.all_apa(tools.anodes),
            ret :: g.intern(
                innodes=[begin],
                outnodes=[end],
                edges=[g.edge(begin, end, i, i) for i in std.range(0, nanodes-1)]
            ),
        }.ret;

    local graph = g.intern(
        innodes=img_clus_pipe,
        outnodes=[clus_all_apa],
        edges=[g.edge(img_clus_pipe[i], clus_all_apa, 0, i) for i in std.range(0, nanodes-1)]
    );

    local app = {
        type: 'Pgrapher', //Pgrapher, TbbFlow
        data: {
            edges: g.edges(graph),
        },
    };

    local cmdline = {
        type: "wire-cell",
        data: {
            plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio", "WireCellSigProc", "WireCellImg", "WireCellRoot", "WireCellTbb", "WireCellClus"],
            apps: ["Pgrapher"]
        }
    };

    [cmdline] + g.uses(graph) + [app]
// img_clus_pipe
// clus_maker.all_apa(tools.anodes)
