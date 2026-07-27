local decon = import "decon.jsonnet";
local detconf = import "detconf.jsonnet";
local det = detconf.pdhd;

local tests = {
    filter_kernel: decon.filter_kernel("name", []),
    local tpc = det.tpcs[0],
    tpc: tpc,
    local groups = tpc.view_groups,
    groups: groups,

    response_kernels: [decon.response_kernel(tpc, index) for index in [0,1,2]],

    tpc_group_decon_kernels: [decon.tpc_group_decon_kernel(tpc, group_index, extra_name="")
                              for group_index in [0,1,2,3]],


    tpc_group_decon_frers: [
        decon.tpc_group_decon_frer(tpc, group_index,
                                   $.response_kernels[groups[group_index].view_index])
        for group_index in [0, 1, 2, 3]],

    group_fanin0: decon.group_fanin($.tpc_group_decon_frers[:1]),
    group_fanin2: decon.group_fanin($.tpc_group_decon_frers[2:]),

    grouped_decon: decon.collect_groups_by_view(tpc, $.tpc_group_decon_frers),

    group_decon_view: decon.group_decon_view(tpc),

    time_filter_one: decon.time_filter_one(tpc, "gauss", 0),
    time_filter_views: decon.time_filter_views(tpc, "gauss"),
};

function(which="")
    if which ==""
    then tests
    else tests[which]

