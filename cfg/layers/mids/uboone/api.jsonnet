// return a mid-level API object

local nf = import "api/nf.jsonnet";
local sp = import "api/sp.jsonnet";
local sim = import "api/sim.jsonnet";
local img = import "api/img.jsonnet";

function(services, params, options) {
    nf : nf(services, params),
    sp : sp(services, params, options),
    img : img(services, params),
} + sim(services, params)

