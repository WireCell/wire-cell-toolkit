// dune10kt-1x2x6/wct-clustering.jsonnet -- clustering job for DUNE FD-HD 1x2x6.
// Promoted VERBATIM from wcp-porting-validation wcfm/wct-clustering.jsonnet (wcfm/docs/02 sec 3,
// docs/25), itself forked BY DUPLICATION from pdhd/wct-clustering.jsonnet, minus the Q/L branch
// (no light here), with the two-face-volume topology of clus.jsonnet:
//
//   per anode: ClusterFileSource(active) -> ClusterFanout(2) -+-> per_face(anode, 0) -> \
//              ClusterFileSource(masked) -> ClusterFanout(2) -+-> per_face(anode, 1) ->  \
//   per face:  PointTreeMerging(all run anodes) -> per_group("groupf<F>")               -> all_tpc -> pctree
//
//   wire-cell --tla-str input=work/000001_1 --tla-code anode_indices=[8,10] \
//             --tla-str output_dir=work/000001_1 --tla-code event=1 \
//             --tla-str save_tensors=work/000001_1/pctree-evt1.tar.gz \
//             -c pgrapher/experiment/dune10kt-1x2x6/wct-clustering.jsonnet
// The pctree (save_tensors) is the input of wct-pr-perevt.jsonnet.

local g = import "pgraph.jsonnet";
local wc = import "wirecell.jsonnet";
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local P = import 'pgrapher/experiment/dune10kt-1x2x6/clus_params.jsonnet';
local params = P.params;
local tools_all = tools_maker(params);

function(
    input = ".",
    anode_indices = std.range(0, std.length(tools_all.anodes) - 1),
    output_dir = '',
    run = 1,
    subrun = 1,
    event = 1,
    time_offset = 0,
    // '' keeps the inert dump_mode sink (trash-all-apa.tar.gz); a path writes the pctree.
    save_tensors = '',
    wrapped_channel_charge = true,
    clus_save_assoc_id = false,
)

local anodes = [tools_all.anodes[i] for i in anode_indices];
local nanodes = std.length(anodes);
local iota = std.range(0, nanodes - 1);

local cluster_source(fname) = g.pnode({
    type: "ClusterFileSource",
    name: fname,
    data: {
        inname: fname,
        anodes: [wc.tn(a) for a in anodes],
    }
}, nin=0, nout=1, uses=anodes);

local clus = import 'pgrapher/experiment/dune10kt-1x2x6/clus.jsonnet';
local clus_maker = clus(output_dir=output_dir, runNo=run, subRunNo=subrun, eventNo=event,
                        time_offset=time_offset, trigger_offset=0);

local actives = [cluster_source("%s/clusters-apa-anode%d-ms-active.tar.gz"%[input, a.data.ident]) for a in anodes];
local maskeds = [cluster_source("%s/clusters-apa-anode%d-ms-masked.tar.gz"%[input, a.data.ident]) for a in anodes];
local fan_live = [g.pnode({ type: 'ClusterFanout', name: 'fanlive-anode%d' % a.data.ident, data: { multiplicity: 2 } }, nin=1, nout=2) for a in anodes];
local fan_dead = [g.pnode({ type: 'ClusterFanout', name: 'fandead-anode%d' % a.data.ident, data: { multiplicity: 2 } }, nin=1, nout=2) for a in anodes];
local face_pipes = [[clus_maker.per_face(anodes[i], face=face, dump=false, wrapped_channel_charge=wrapped_channel_charge) for i in iota] for face in [0, 1]];
local groups = P.face_groups(anodes);
local group_pipes = [clus_maker.per_group(anodes, gd.name, gd.face, dump=false, save_assoc_id=clus_save_assoc_id) for gd in groups];
local clus_all_tpc = clus_maker.all_tpc(anodes, ngroups=std.length(groups), tensor_outname=save_tensors, save_assoc_id=clus_save_assoc_id);

local graph = g.intern(
    innodes = actives + maskeds,
    centernodes = fan_live + fan_dead + face_pipes[0] + face_pipes[1] + group_pipes,
    outnodes = [clus_all_tpc],
    edges =
        [g.edge(actives[i], fan_live[i], 0, 0) for i in iota] +
        [g.edge(maskeds[i], fan_dead[i], 0, 0) for i in iota] +
        std.flattenArrays([[g.edge(fan_live[i], face_pipes[face][i], face, 0), g.edge(fan_dead[i], face_pipes[face][i], face, 1)]
                           for i in iota for face in [0, 1]]) +
        std.flattenArrays([[g.edge(face_pipes[gi][i], group_pipes[gi], 0, i) for i in iota] for gi in std.range(0, std.length(groups) - 1)]) +
        [g.edge(group_pipes[gi], clus_all_tpc, 0, gi) for gi in std.range(0, std.length(groups) - 1)]
);

local app = { type: 'Pgrapher', data: { edges: g.edges(graph) } };
local cmdline = {
    type: "wire-cell",
    data: {
        plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio", "WireCellSigProc", "WireCellImg", "WireCellTbb", "WireCellClus", "WireCellAux"],
        apps: ["Pgrapher"]
    }
};

[cmdline] + g.uses(graph) + [app]
