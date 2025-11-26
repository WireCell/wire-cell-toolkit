// Build drifter.  See also sim.jsonnet and detsim.jsonnet.

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

function(det, control) {
        
    ntpcs: std.length(det.tpcs),

    // Collect unique xregions from the anodes.  It is because the drifter
    // operates at whole-detector level that these stages can not of a
    // single_tpc scope.
    xregions: 
        wc.unique_objects(std.prune(std.flattenArrays([tpc.anode.data.faces for tpc in det.tpcs ]))),


    /// Create a per-depo drifter.  Pretty much never use this bare.

    // A drifter of depo sets.
    drifter: 
        local depo_drifter = {
            type: 'Drifter',
            name: det.name,
            data: det.lar + {
                xregions: $.xregions,
                time_offset: 0.0,     // fixme: need to expose this to user
                fluctuate: true,      // fixme: this too
                scale_factor: 1.0,    // fixme: and this
                rng: wc.tn(control.rng),
            },
            uses: [control.rng],
        };
        pg.pnode({
            type: "DepoSetDrifter",
            name: det.name,
            data: {
                drifter: wc.tn(depo_drifter)
            },
        }, nin=1, nout=1, uses=[pg.pnode(depo_drifter, nin=1, nout=1)]),


}
