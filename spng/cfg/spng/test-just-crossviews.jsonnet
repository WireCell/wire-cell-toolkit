// debug
local pg = import "pgraph.jsonnet";
local crossviews = import "spng/crossviews.jsonnet";
local detconf = import "spng/detconf.jsonnet";
local cv = import "spng/crossviews.jsonnet";
local det = detconf["pdhd"];
local graph = cv.crossfan(det.tpcs[0]);
pg.main(graph)
