// This file defines the base or default "mid" layer API implementation.
//
// People looking to configure specific detectors should read it but very likely
// should NOT modify it.  The comments give expectation for each subgraph.  In
// particular Frames passed between subgraphs HAVE NO TRACE TAGS.
//
// This is called through midapi.jsonnet to assure interface correctness.
//
// API methods with no implementation will ASSERT and thus must be provided by a
// detector implementation.

local low = import "../../low.jsonnet";
local tn = low.wc.tn;
local pg = low.pg;

function(services, params, options={})
{
    // Return a list of all anode config objects.
    anodes() :
        low.anodes(params.geometry.volumes, params.geometry.wires_file),


    // A node that returns a depo set drifter.
    drifter(name="") :
        low.drifter(services.random,
                    low.util.driftsToXregions(params.geometry.volumes),
                    params.lar, name=name),

    // DepoSet to voltage Frame.  NO TRACE TAGS.
    signal(anode, name) :
        std.assertEqual("signal" == "implemented"),
    
    // Input voltage Frame, output noisy voltage Frame.  NO TRACE TAGS.
    noise(anode, name) :
        std.assertEqual("signal" == "implemented"),
    
    // Voltage Frame to ADC Frame. NO TRACE TAGS
    digitizer(anode, name) :
        pg.pnode({
            type: 'Digitizer',
            name: name,
            data: params.digi {
                anode: tn(anode),
            }
        }, nin=1, nout=1, uses=[anode]),
    
    // Pull out the responses as a short hand variable.
    local sim_res = low.resps(params).sim,

    // DepoSet to signal Frame (an ersatz sim+SP).  NO TRACE TAGS.
    splat(anode, name) :
        pg.pnode({
            type: 'DepoFluxSplat',
            name: name,
            data: params.splat {
                anode: tn(anode),
                field_response: tn(sim_res.fr), // note, only uses (period, origin)
            },
        }, nin=1, nout=1, uses=[anode, sim_res.fr]),

    // noisy ADC frame to quiet ADC frame.  NO TRACE TAGS.
    nf(anode, name) :
        std.assertEqual("nf" == "implemented"),

    // ADC frame to signal frame.  NO TRACE TAGS.
    sp(anode, name) :
        std.assertEqual("sp" == "implemented"),

    // Signal frame to better signal Frame.  NO TRACE TAGS
    dnnroi(anode, name) :
        std.assertEqual("dnnroi" == "implemented"),

    // Frame to Frame
    reframer(anode, name, tags=[]) :
        low.reframer(params, anode, name=name, tags=tags),

    // frame to blob set
    img(anode, name) :
        std.assertEqual("img" == "implemented"),

}
