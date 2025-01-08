// This defines and enforces the API that every detector must provide.
//
// It "works" merely by forwarding calls to the a detector implementation (imp)
// in order to enforce the arguments.  For methods that return a DFP subgraph
// the comments def expectation for input/output edges.
//
// This file should not be modified except to extend support to novel subgraph
// types.
//
// See mids/base/api.jsonnet for a default implementation that also gives more
// guidance on what detector implementations must provide.

local util = import "low/util.jsonnet";

function (imp)
    local call_it(meth)={
        [meth]: function(anode, name=null)
            imp[meth](anode, name=util.assure_name(name, anode))
    };
    local call_it_tags(meth)={
        [meth]: function(anode, name=null, tags=[])
            imp[meth](anode, name=util.assure_name(name, anode), tags=tags)
    };
{

    // Return list of anode configuration objects.
    anodes() : imp.anodes(),

    // Input DepoSet (original depos) output DepoSet (drifted depos). 
    drifter(name="") : imp.drifter(name),

  track_depos : imp.track_depos

}
+ call_it("signal")
+ call_it("noise")
+ call_it("digitizer")
+ call_it("splat")
+ call_it("nf")
+ call_it("sp")
+ call_it_tags("reframer")
+ call_it("img")
        
    
