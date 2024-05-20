local low = import "../../../low.jsonnet";

// These are big and unique to PDSP so we have to write them out
// exhaustively.
local nf = import "nf.jsonnet";
local sp = import "sp.jsonnet";

function(services, params, options={}) {

    // API method
    nf :: nf(services, params, options),

    // API method
    sp :: sp(services, params, options),

    // API method
    dnnroi :: function(anode)
        low.dnnroi(anode, params.dnnroi.model_file,
                   services.platform, params.dnnroi.output_scale),

}
