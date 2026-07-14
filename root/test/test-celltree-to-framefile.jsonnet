// This converts a celltree file to a WCT frame file
local high = import "layers/high.jsonnet";

// This perhaps only works on uboone.
function(infile="celltree.root", outfile="frames.npz", eventid="6501")
         // detector="uboone", variant="nominal")

    local trace_tags = ["wiener", "gauss"];

    local source = high.fio.celltree_file_source(
        infile, eventid,
        branches = ["calibWiener", "calibGaussian"],
        frame_tags = ["gauss","wiener"],
        trace_tags = trace_tags,
        extra_params = {in_branch_thresholds: ["channelThreshold", "channelThreshold"]});

    local target = high.fio.frame_file_sink(outfile, trace_tags);
        
    local graph = high.pg.pipeline([source, target]);
    high.main(graph, extra_plugins=['WireCellRoot'])
