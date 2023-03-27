// Generate disconnected blobify subgraphs for different detectors.
// For each detector, the subgraph is constructed over the predefiend
// blobify strategies.  Note, this combination of multiple detectors
// in the same graph is not something a "real" job would likely ever
// do.  It is done here to make a single dotify PDF showing
// everything.

local high = import "layers/high.jsonnet";
local low = import "layers/low.jsonnet";

// Make a composite of subgraphs over all defined strategies
local one_det = function(detname, variant="nominal") 
    local det = high.mid(detname, variant);
    local anode = det.anodes()[0];
    local img = low.img(anode);
    
    local graphs = [
        img.blobify("%s-%s-%s"%[detname, variant, strat],
                    "gauss", "gauss_err",
                    img.blobify_strategies[strat])
        for strat in std.objectFields(img.blobify_strategies)];
                    
    // graphs[0];
    // high.pg.intern(centernodes=graphs, name='%s-%s'%[detname,variant]);
    high.pg.fan.pipe('FrameFanout', graphs, 'BlobSetMerge', name='%s-%s'%[detname,variant]);


local many_det() =
    local graphs = [one_det(detname) for detname in ["uboone", "pdsp"]];
    high.pg.intern(centernodes=graphs, name='world');

high.main(many_det(), 'TbbFlow')
// high.main(one_det('pdsp'), 'TbbFlow')
// one_det('uboone')             
// local det = high.mid('uboone');
// local anode = det.anodes()[0];
// local img = low.img(anode);
// img.blobify_strategies["live"]
// img.blobify_single("blah", "gauss", "gauss_err", img.blobify_strategies["live"][0])
// img.blobify("blah", "gauss", "gauss_err", img.blobify_strategies["live"])

