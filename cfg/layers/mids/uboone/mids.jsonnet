// The "mid level" API factory for "uboone" detector (MicroBooNE).

// The "real" code is elsewhere
local api = import "api.jsonnet";

// Import all variant parameter packs to allow for dict like lookup.
local variants = import "variants.jsonnet";

function(services, variant="", options={})
    api(services, variants[variant], options)

