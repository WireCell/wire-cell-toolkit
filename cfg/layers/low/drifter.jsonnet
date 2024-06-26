/**
 A low level method to produce a drifter subgraph given
 detector-specific parameter values.

 Arguments:

 - rnd :: an IRandom configuration object
 
 - volumes :: an array of "volume" objects which define the volume
   relevant to each anode.  Each volume object must have a .faces
   attribute with two elements each holding an object with trio of
   plane X-locations or a null.  The trio object must have three
   attributes: "anode", "response" and "cathode" giving the location
   of these planes along the global X-coordinate (same in which the
   wire endpoints are expressed).

 - offset :: added to the central time of each depo.

 - fluctuate :: if true, apply statistical fluctuations related to
   electron absorption occuring during the drift.

*/
local pg = import "pgraph.jsonnet";
local wc = import "wirecell.jsonnet";
local u = import "util.jsonnet";
function(random, xregions, lar, offset=0, fluctuate=true, name="")
    pg.pnode({
        type: 'DepoSetDrifter',
        name: name,
        data: { drifter: "Drifter:"+name }
    }, nin=1, nout=1, uses = [
        pg.pnode({
            type: 'Drifter',
            name: name,
            data: lar {
                rng: wc.tn(random),
                xregions: xregions, // for help, see: low.util.driftsToXregions(drifts)
                time_offset: offset,
                fluctuate: fluctuate,
            },
        }, nin=1, nout=1, uses=[random])
    ])
        
