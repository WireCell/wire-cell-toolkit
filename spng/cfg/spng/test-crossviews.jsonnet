local cv = import "crossviews.jsonnet";
local decon = import "decon.jsonnet";
local detconf = import "detconf.jsonnet";
local det = detconf.pdhd;

local tests = {
    local tpc = det.tpcs[0],
    tpc: tpc,
    tv: cv.threshold_views(tpc),

    cv0: cv.crossview(tpc, 0),

    cf: cv.crossfan(tpc, $.tv),
};
function(which="")
    if which ==""
    then tests
    else tests[which]

