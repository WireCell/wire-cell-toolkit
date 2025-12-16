// FIXME: WARNING: this is a temporary construction and will go away.  End goal
// is to specify information for OSP + DNNROI which is independent of the "tpc"
// and "control" objects to allow varied ways of constructing the OSP and DNNR
// subgraphs based on high order information (eg, which wire planes to have
// DNNROI applied to).  For expediency to get a working OSP chain, we have this
// layer provide ready-to-connect subgraph.


local pg = import "pgraph.jsonnet";
local sp_mod = import "sp.jsonnet";
local dnnroi_mod = import "dnnroi.jsonnet";

// return OSP + DNNROI subgraph
function(tpc) {
    
    local sp = sp_mod(tpc),
    local dnnroi = dnnroi_mod(tpc),
    sp: sp, 
    dnnroi: dnnroi,
    osp: pg.pipeline([sp, dnnroi])
}
