local decon = import "decon.jsonnet";
local detconf = import "detconf.jsonnet";
local det = detconf.pdhd;

local tests = {
    filter_kernel: decon.filter_kernel("name", []),
    local tpc = det.tpcs[0],
    tpc: tpc,

    response_kernels: [decon.response_kernel(tpc, index) for index in [0,1,2]],

    tpc_group_decon_kernels: [decon.tpc_group_decon_kernel(tpc, group_index, extra_name="")
                              for group_index in [0,1,2,3]],

};

function(which="")
    if which ==""
    then tests
    else tests[which]

