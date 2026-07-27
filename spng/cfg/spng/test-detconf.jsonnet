local detconf = import "detconf.jsonnet";
local det = detconf.pdhd;


local tests = {
    local tpc = det.tpcs[0],
    tpc: tpc,
    anode: tpc.anode,
    fr: tpc.fr,
    er: tpc.er,
    rcs: tpc.rcs,
    view_groups: tpc.view_groups,
    filters: tpc.filters,
    one_channel_filter: tpc.filters[0].channel.decon,
};

function(which="")
    if which ==""
    then tests
    else tests[which]

