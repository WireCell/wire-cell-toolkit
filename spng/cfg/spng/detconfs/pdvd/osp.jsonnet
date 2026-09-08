// PDVD OSP + DNNROI subgraph maker.  Same shape as spng/detconfs/pdhd/osp.jsonnet.
//
// detector.subset() eagerly calls tpc._osp_subgraphs(tpc, device), so every
// detector that flows through subset() (including the SPNG-native adc-to-spng
// job) must provide this even when only the SP node is used.

local pg = import 'pgraph.jsonnet';
local sp_mod = import 'sp.jsonnet';
local dnnroi_mod = import 'dnnroi.jsonnet';

// return OSP + DNNROI subgraph
function(tpc, device='cpu') {

    local sp = sp_mod(tpc),
    local dnnroi = dnnroi_mod(tpc, device=device),
    sp: sp,
    dnnroi: dnnroi,
    osp: pg.pipeline([sp, dnnroi]),
}
