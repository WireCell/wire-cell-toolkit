// return a mid-level API object

local nf = import "api/nf.jsonnet";
local sp = import "api/sp.jsonnet";
local sim = import "api/sim.jsonnet";
local img = import "api/img.jsonnet";

// Create a mid API object.  No options supported.
function(services, params, options={}) {
    nf : nf(services, params, options),
    sp : sp(services, params, options),
    img : img(services, params)
} + sim(services, params)
