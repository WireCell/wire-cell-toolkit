local frame = import "frame.jsonnet";
local detconf = import "detconf.jsonnet";
local det = detconf.pdhd;


local tests = {
    local tpc = det.tpcs[0],
    tpc: tpc,
    connections: tpc.connections,
    view_groups: tpc.view_groups,

    frame_to_tdm_views: frame.to_tdm(tpc),

    frame_tensorset_unpacker: frame.tensorset_unpacker(tpc),

    group_view_indices: [g.view_index for g in tpc.view_groups],
    

};

function(which="")
    if which ==""
    then tests
    else tests[which]

