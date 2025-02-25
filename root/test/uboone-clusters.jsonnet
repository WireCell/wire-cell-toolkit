// This will configure WCT to load in a Uboone WCP result file and exercise a
// chain that includes forming blobs, clusters, optical data and cluster-flash
// association.  It ends with a simple sink with logging.
//
// The tier of data sunk is a super-set of what is produced by WCT
// imaging+clustering.  This job produces an equivalent clustering pc-tree but
// adds the optical data and association.
//
// Input is a ROOT file (or files) providing Uboone TTrees Trun, TC, T_light, T_match.


local high = import "layers/high.jsonnet";
local wc = high.wc;
local pg = high.pg;
local ub = import "uboone.jsonnet";

local graph(infiles, datapath="pointtrees/%d/uboone") = pg.pipeline([
    ub.multiplex_blob_views(infiles, 'live', ["uvw","uv","vw","wu"]),
    ub.UbooneClusterSource(infiles, datapath=datapath),
    ub.ClusterFlashDump(datapath=datapath)
]);

local extra_plugins = ["WireCellAux", "WireCellRoot", "WireCellClus"];

function(infiles="uboone.root")
    local g = graph(wc.listify(infiles));
    ub.main(g, "Pgrapher", extra_plugins)
    
