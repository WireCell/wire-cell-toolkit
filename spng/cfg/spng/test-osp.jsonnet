local detconf = import "detconf.jsonnet";
local det = detconf.pdhd;
local pgraph = import "pgraph.jsonnet";

pgraph.main(det.tpcs[0].osp_subgraph)

