local basebase = import "../mids/base/variants.jsonnet";
local base = basebase.nominal {
    detector_name: "base",
    variant_name: "nominal",
    structure_name: "nominal"
};
    
std.assertEqual("base", base.detector_data.detname)

