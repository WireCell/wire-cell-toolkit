// combine drift, and a parameterized model of detsim+SP
//
local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local drift = import "drift.jsonnet";
local sim = import "sim.jsonnet";


// 1->1 deposet to frame of signals.
function(det, control={})
    local drifter = drift(det, control).drifter;
    local the_sims = [sim(tpc, control) for tpc in det.tpcs];
    local pipes = [pg.pipeline([
        the_sim.splat,
        the_sim.reframer("_splat"),
    ]) for the_sim in the_sims];

    local dr = if std.length(pipes) == 1
               then pipes[0]
               else pg.fan.fanout('DepoSetFanout', pipes);


    pg.pipeline([drifter, dr])
