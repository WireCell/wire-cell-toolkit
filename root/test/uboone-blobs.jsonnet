// Clustering testing with Uboone and blobs originating from WCP.
//
// This provides config for various "kinds" of subgraphs to be run separately
// and be connected together via files.:
// 
// - live :: load WCP "live" blobs from root file, produce npz
// - live :: load WCP "dead" blobs from root file, produce npz
// - clus :: read back live and dead npz files produce bee and tensor files
//
// Best to run full graph via uboone-blobs.smake

local ub = import "uboone.jsonnet";
local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";


local default_files = {
    root: "result_5384_130_6501.root",
    live: "live-clus.npz",
    dead: "dead-clus.npz",
    tens: "live-dead.npz",
    beez: "live-dead.bee"
};

local live(iname, oname) = pg.pipeline([
    ub.multiplex_blob_views(iname, "live", ["uvw","uv","vw","wu"]),
    ub.BlobClustering("live"),
    // BlobGrouping("0"),

    // "standard":
    // ProjectionDeghosting("1"),
    // BlobGrouping("1"), ChargeSolving("1a","uniform"), LocalGeomClustering("1"), ChargeSolving("1b","uboone"),
    // InSliceDeghosting("1",1),
    // ProjectionDeghosting("2"),
    // BlobGrouping("2"), ChargeSolving("2a","uniform"), LocalGeomClustering("2"), ChargeSolving("2b","uboone"),
    // InSliceDeghosting("2",2),
    // BlobGrouping("3"), ChargeSolving("3a","uniform"), LocalGeomClustering("3"), ChargeSolving("3b","uboone"),
    // InSliceDeghosting("3",3),
    ub.GlobalGeomClustering(""),
    ub.ClusterFileSink(oname),
]);

local dead(iname, oname) = pg.pipeline([
    ub.multiplex_blob_views(iname, "dead", ["uv","vw","wu"]),
    ub.BlobClustering("dead"),
    ub.GlobalGeomClustering("", "dead_clus"),
    ub.ClusterFileSink(oname),
]);

local clus = function(iname, oname)
    local live_dead = std.split(iname, ",");
    local beezip_tensor = std.split(oname, ",");
    pg.pipeline([
        ub.point_tree_source(live_dead[0], live_dead[1]),
        ub.MultiAlgBlobClustering(beezip_tensor[0]),
        ub.TensorFileSink(beezip_tensor[1])]);
    
local extra_plugins = ["WireCellRoot","WireCellClus"];

function(iname, oname, kind /*clus or live or dead*/)
    if kind == "live" then
        ub.main(live(iname, oname), "Pgrapher", extra_plugins)
    else if kind=="dead" then
        ub.main(dead(iname, oname), "Pgrapher", extra_plugins)
    else if kind=="clus" then
        ub.main(clus(iname, oname), "Pgrapher", extra_plugins)
    else
        error "uknown kind: "+kind

