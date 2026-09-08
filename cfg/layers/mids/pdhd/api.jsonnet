// The mid-level API for pdhd.
//
// This reuses the PDSP mid implementation: the api/ sub-tree (sim, sp, frs,
// sp-filters, nf, img, chndb-*, channel-*) is parameterized by the variant
// params, and only the detector-specific values live in variants/ and
// api/sp-filters.jsonnet.
//
// SCOPE: the drift/splat/sim/sp (OSP) path used by omnijob + test/scripts/spdir
// is supported and validated.  The noise-filter (nf) and imaging (img) chains
// are inherited verbatim from PDSP -- their channel-map / chndb data
// (channel-groups, chndb-*, channel-info, actual-bad-channels, rel-gain) are
// PDSP-specific and DORMANT here (spdir never calls nf/img).  Replace them with
// pdhd channel maps before using the nf/img API for this detector.


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
