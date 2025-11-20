// combine drift, detector response and noise
//
// This can take a list of TPCs to use and will default to all in det if not given.
//
// It returns a 1 deposet -> N_tpc frame subgraph

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local drift = import "drift.jsonnet";
local sim = import "sim.jsonnet";

function(det, control={})
    local drifter = drift(det, control).drifter;
    local the_sims = [sim(tpc, control) for tpc in det.tpcs];
    local pipes = [pg.pipeline([
        the_sim.depo_transform,
        the_sim.reframer,
        the_sim.addnoise_empirical,
        the_sim.digitizer
    ]) for the_sim in the_sims];

    local dr = if std.length(pipes) == 1
               then pipes[0]
               else pg.fan.fanout('DepoSetFanout', pipes);


    pg.pipeline([drifter, dr])
