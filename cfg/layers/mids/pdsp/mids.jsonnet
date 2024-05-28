// The "mid level" API factory for "pdsp" detector (ProtoDUNE-SP).

// The "real" code is elsewhere
local api = import "api.jsonnet";


// Import all variant parameter packs to allow for dict like lookup.
local variants = {
    nominal: import "variants/nominal.jsonnet",
} + import "variants/ssss.jsonnet";

function(services, variant="", options={})
    api(services, variants[variant], options)
