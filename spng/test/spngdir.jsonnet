// This configures jobs run the by spngdir Snakemake workflow.
//


local wc = import "wirecell.jsonnet";
local omnimap = import "omni/all.jsonnet";

// Return a task object for possible tasks in an anode context.
//
// Each task consists of one or more objects from omni method calls and
// conceptually represents a 1-to-1 functional node/subgraph.
//
local tasks_factory(anode, omni) = {
    tasks : {
        // depo to frame
        splat: omni.splat(anode),
        
        
    }

}


// Return a file source by interpreting the stage.
local source(omni, stage, input) = 
    if std.member(["drift","splat","sim"], stage)
    then omni.depo_file_source(input)
    else omni.frame_file_source(input);



function(name, input, output, tasks):
    local omni = omnimap[name];
    local stages = wc.listify(tasks);
    .... continue .....
map stage over anodes.
