// drifter-graph.jsonnet
local pg = import "pgraph.jsonnet";
local drift = import "drifter.jsonnet";

pg.main(drift.drifter())
