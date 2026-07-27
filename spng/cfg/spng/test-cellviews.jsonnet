
// This produces a final graph to input depos and output crossview tensors

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

local detconf = import "spng/detconf.jsonnet";
local detector = import "spng/detector.jsonnet";
local control_js = import "spng/control.jsonnet";

local cv_js = import "spng/cellviews.jsonnet";

function(detname='pdhd', tpcids=[], engine='Pgrapher', device='cpu', verbosity=0)
    
    local controls = control_js(device=device, verbosity=wc.intify(verbosity));
    local det = detector.subset(detconf[detname], tpcids);
    local cv = cv_js(controls.config, pg);
    local graph = cv.cellviews_tensors(det.tpcs[0]);
    pg.main(graph, engine, plugins=["WireCellSpng"])
    // graph
