// combine drift, detector response and noise

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local drift = import "drift.jsonnet";
local sim = import "sim.jsonnet";

function(det, tpc, control) {

    drifter: drift(det, control).drifter, 

    local the_sim = sim(tpc, control),

    drifted_to_adc: pg.pipeline([
        the_sim.depo_transform,
        the_sim.reframer,
        the_sim.addnoise_empirical,
        the_sim.digitizer
    ]),

    depos_to_adc: pg.pipeline([$.drifter, $.drifted_to_adc]),

}
