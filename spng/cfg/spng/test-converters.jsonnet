local converters = import "converters.jsonnet";
local detconf = import "detconf.jsonnet";
local det = detconf.pdhd;


local tests = {
    local tpc = det.tpcs[0],
    tpc: tpc,
    connections: tpc.connections,

    each: {
        [std.toString(i)]:converters.view_group_wpids(tpc, i)
        for i in [0,1,2]
    },

    groups: converters.tpc_group_wpids(tpc),

    frame_to_tdm_views: converters.frame_to_tdm_views(tpc),

    frame_tensorset_unpacker: converters.frame_tensorset_unpacker(tpc),

    plane_indices: [converters.tpc_group_planes(tpc)[group_index]
                    for group_index in [0,1,2,3]],
                    

};

function(which="")
    if which ==""
    then tests
    else tests[which]

