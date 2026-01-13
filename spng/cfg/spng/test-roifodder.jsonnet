// This produces an intermediate subgraph that spans splat based roi production.

local pg = import "pgraph.jsonnet";

local detconf = import "spng/detconf.jsonnet";
local control_mod = import "spng/control.jsonnet";
local roifodder_js = import "spng/roifodder.jsonnet";

function(detname='pdhd', tpcid=0, engine='Pgrapher', device='cpu', verbosity=0)
    
    local controls = control_mod(device=device, verbosity=verbosity);
    local det = detconf.get(detname, [tpcid]);
    local tpc = det.tpcs[0];
    local graph = roifodder_js(tpc, controls.config);

    pg.main(graph, app=engine, plugins=["WireCellSpng"], uses=controls.uses)
    //graph
    //tpc.view_groups


    
