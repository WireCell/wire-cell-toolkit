// The "mid level" API factory for "pdsp" detector (ProtoDUNE-SP).

// The "real" code is elsewhere
local api = import "api.jsonnet";

local variants = import "variants.jsonnet";

function(services, variant="", options={})
    api(services, variants[variant], options)
