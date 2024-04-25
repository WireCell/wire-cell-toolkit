// return a mid-level API object

local low = import "../../low.jsonnet";
local sim = import "api/sim.jsonnet";
local sigproc = import "api/sigproc.jsonnet";
local img = import "api/img.jsonnet";

// Create a mid API object.  No options supported.
function(services, params, options={}) {

    local pars = std.mergePatch(params, std.get(options, "params", {})),

    anodes :: function()
        low.anodes(pars.geometry.drifts, pars.geometry.wires_file),

    drifter :: function(name="")
        low.drifter(services.random,
                    low.util.driftsToXregions(pars.geometry.drifts),
                    pars.lar, name=name),

    // track_depos, signal, noise, digitizer
    sim :: sim(services, pars),

    // nf, sp, dnnroi
    sigproc :: sigproc(services, pars, options),

    img :: img(services, pars),
}

    
