local high = import "../high.jsonnet";

local un = high.params("uboone", "nominal");
local ua = high.params("uboone", "actual");

[
    std.assertEqual(un.detector_name, "uboone"),
    std.assertEqual(un.variant_name, "nominal"),
    std.assertEqual(un.detector_data.wires, "microboone-celltree-wires-v2.1.json.bz2"),

    std.assertEqual(ua.detector_name, "uboone"),
    std.assertEqual(ua.variant_name, "actual"),
    std.assertEqual(ua.detector_data.wires, "microboone-celltree-wires-v2.1.json.bz2"),
    
]    
    
